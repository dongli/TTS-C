! **************************************************************************** !
! Module TTS - Trajectory Tracking Scheme                                      !
!                                                                              !
! Description:                                                                 !
!   Trajectory  tracking  scheme  (TTS)  is an  algorithm  for solving  linear !
!   advection  problem.  Its  object  is the linear  advection  equation under !
!   Lagrangian  framework.  The  main  ingradients are  characteristic   curve !
!   (trajectory),  some kind of interpolation of background fluid velocity and !
!   some  kind of  numerical  methods for solving ODE (trajectory equation and !
!   physical attribute equation).                                              !
!   The   highlights   of this  method  is that it can keep the  discontinuity !
!   in solution very well.                                                     !
!                                                                              !
! Authors:                                                                     !
!   DONG Li, dongli@lasg.iap.ac.cn                                             !
!                                                                              !
! References:                                                                  !
! **************************************************************************** !

module TTS

    use MsgManager
    use RunManager
    use TracerManager
    use MeshManager
    use FlowManager

    implicit none

    private

    public TTS_AdvectTracer

contains

    ! ************************************************************************ !
    ! TTS_AdvectTracer                                              !
    ! Purpose:                                                                 !
    !   Advect tracers for one time step.                                      !
    ! Details:                                                                 !
    !   Four-order Runge-Kutta numerical method is applied.                    !
    ! ************************************************************************ !

    subroutine TTS_AdvectTracer(tracerId, dt)
        integer, intent(in) :: tracerId
        real(8), intent(in) :: dt

        integer numSample
        real(8) v(numDim), div

        real(8) x(numDim), q
        type(Location) loc

        real(8) x_tmp(numDim), q_tmp, dt05
        type(Location) loc_tmp
        real(8) k1_x(numDim), k2_x(numDim), k3_x(numDim), k4_x(numDim)
        real(8) vt1(numDim), vt2(numDim), vt3(numDim), vt4(numDim)
        real(8) k1_q, k2_q, k3_q, k4_q

        integer i

        call MsgManager_RecordSpeaker("TTS_AdvectTracer")

        dt05 = dt*0.5d0

        call TracerManager_GetSampleNumber(tracerId, numSample)

        ! numerical integration
        if (TracerManager_IsConservative(tracerId)) then
            do i = 1, numSample
                call TracerManager_GetSampleCoordinate(tracerId, i, x)
                call TracerManager_GetSampleQuantity(tracerId, i, q)
                call TracerManager_GetSampleLocation(tracerId, i, loc)
                ! ************ four-order Runge-Kutta ODE solver ************* !
                if (MeshManager_IsNearPole(loc)) then
                    ! first step
                    call FlowManager_GetTransformedVelocity(x, loc, vt1, "old")
                    call FlowManager_GetDivergence(x, loc, div, "old")
                    k1_q = -q*div
                    ! second step
                    call MeshManager_MoveNearPole(x, x_tmp, vt1, dt05, loc)
                    call MeshManager_LocationCheck(x_tmp, loc_tmp)
                    q_tmp = q+dt05*k1_q
                    call FlowManager_GetTransformedVelocity(x_tmp, loc_tmp, vt2, "half")
                    call FlowManager_GetDivergence(x_tmp, loc_tmp, div, "half")
                    k2_q = -q_tmp*div
                    ! third step
                    call MeshManager_MoveNearPole(x, x_tmp, vt2, dt05, loc)
                    call MeshManager_LocationCheck(x_tmp, loc_tmp)
                    q_tmp = q+dt05*k2_q
                    call FlowManager_GetTransformedVelocity(x_tmp, loc_tmp, vt3, "half")
                    call FlowManager_GetDivergence(x_tmp, loc_tmp, div, "half")
                    k3_q = -q_tmp*div
                    ! fourth step
                    call MeshManager_MoveNearPole(x, x_tmp, vt3, dt, loc)
                    call MeshManager_LocationCheck(x_tmp, loc_tmp)
                    q_tmp = q+dt*k3_q
                    call FlowManager_GetTransformedVelocity(x_tmp, loc_tmp, vt4, "new")
                    call FlowManager_GetDivergence(x_tmp, loc_tmp, div, "new")
                    k4_q = -q_tmp*div
                    ! summrize
                    v = (vt1+2.0d0*vt2+2.0d0*vt3+vt4)/6.0d0
                    call MeshManager_MoveNearPole(x, x_tmp, v, dt, loc)
                else
                    ! first step
                    call FlowManager_GetVelocity(x, loc, k1_x, "old")
                    call FlowManager_GetDivergence(x, loc, div, "old")
                    k1_q = -q*div
                    ! second step
                    call MeshManager_Move(x, x_tmp, k1_x, dt05, loc) ! <--- loc
                    call MeshManager_LocationCheck(x_tmp, loc_tmp)
                    q_tmp = q+dt05*k1_q
                    call FlowManager_GetVelocity(x_tmp, loc_tmp, k2_x, "half")
                    call FlowManager_GetDivergence(x_tmp, loc_tmp, div, "half")
                    k2_q = -q_tmp*div
                    ! third step
                    call MeshManager_Move(x, x_tmp, k2_x, dt05, loc) ! <--- loc
                    call MeshManager_LocationCheck(x_tmp, loc_tmp)
                    q_tmp = q+dt05*k2_q
                    call FlowManager_GetVelocity(x_tmp, loc_tmp, k3_x, "half")
                    call FlowManager_GetDivergence(x_tmp, loc_tmp, div, "half")
                    k3_q = -q_tmp*div
                    ! fourth step
                    call MeshManager_Move(x, x_tmp, k3_x, dt, loc) ! <--- loc
                    call MeshManager_LocationCheck(x_tmp, loc_tmp)
                    q_tmp = q+dt*k3_q
                    call FlowManager_GetVelocity(x_tmp, loc_tmp, k4_x, "new")
                    call FlowManager_GetDivergence(x_tmp, loc_tmp, div, "new")
                    k4_q = -q_tmp*div
                    ! summerize
                    v = (k1_x+2.0d0*k2_x+2.0d0*k3_x+k4_x)/6.0d0
                    call MeshManager_Move(x, x_tmp, v, dt, loc)
                end if
#if (defined LAGRANGE_TO_EULER)
                call MeshManager_LocationCheck(x_tmp, loc_tmp, tracerId, i)
#else
                call MeshManager_LocationCheck(x_tmp, loc_tmp)
#endif
                q_tmp = q+dt*(k1_q+2.0d0*k2_q+2.0d0*k3_q+k4_q)/6.0d0
                ! end
                call TracerManager_SaveSampleCoordinate(tracerId, i, x_tmp)
                call TracerManager_SaveSampleQuantity(tracerId, i, q_tmp)
                call TracerManager_SaveSampleLocation(tracerId, i, loc_tmp)
            end do
        else if (TracerManager_IsAdvective(tracerId)) then
            print *, ":( Advective form TTS is not finished."
            call RunManager_EndRun
            do i = 1, numSample
                call TracerManager_GetSampleCoordinate(tracerId, i, x)
                call TracerManager_GetSampleLocation(tracerId, i, loc)
                ! ************ four-order Runge-Kutta ODE solver ************* !
                ! first step
                call FlowManager_GetVelocity(x, loc, v, "old")
                k1_x = v
                ! second step
                call MeshManager_Move(x, x_tmp, k1_x, 0.5*dt, loc)
                call FlowManager_GetVelocity(x_tmp, loc, v, "half")
                k2_x = v
                ! third step
                call MeshManager_Move(x, x_tmp, k2_x, 0.5*dt, loc)
                call FlowManager_GetVelocity(x_tmp, loc, v, "half")
                k3_x = v
                ! fourth step
                call MeshManager_Move(x, x_tmp, k3_x, dt, loc)
                call FlowManager_GetVelocity(x_tmp, loc, v, "new")
                k4_x = v
                ! summerize
                v = (k1_x+2*k2_x+2*k3_x+k4_x)/6.0d0
                call MeshManager_Move(x, x_tmp, v, dt, loc)
#if (defined LAGRANGE_TO_EULER)
                call MeshManager_LocationCheck(x_tmp, loc_tmp, tracerId, i)
#else
                call MeshManager_LocationCheck(x_tmp, loc_tmp)
#endif
                x = x_tmp
                ! end
                call TracerManager_SaveSampleCoordinate(tracerId, i, x)
                call TracerManager_SaveSampleLocation(tracerId, i, loc)
            end do
        end if

        call MsgManager_DeleteSpeaker

        return
    end subroutine TTS_AdvectTracer

end module TTS
