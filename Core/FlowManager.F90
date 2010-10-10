! *************************************************************************** !
! FlowManager module                                                          !
!                                                                             !
! Description:                                                                !
!                                                                             !
! Author:                                                                     !
!                                                                             !
! *************************************************************************** !

module FlowManager

    use MsgManager
    use RunManager
    use CommonTypes
    use SphereService
    use FloatingPoint
    use MeshManager
    use NFWrap

    implicit none

    private

    public FlowManager_Init
    public FlowManager_Final
    public FlowManager_Link
    public FlowManager_BeforeAdvance
    public FlowManager_AfterAdvance
    public FlowManager_GetVelocity
    public FlowManager_GetTransformedVelocity
    public FlowManager_GetDivergence
    public FlowManager_Output

    type(Quantity) u, v, div

    type, extends(PolarCapQuantity) :: PolarCapSpeed
        real(8), allocatable :: sinLon(:), cosLon(:)
        ! polar stereo transformed speed --->
        real(8) sinLat(2), sinLat2(2)
        type(TwoTimeLevel), allocatable :: transform(:,:)
        ! <---
    contains
        procedure :: initExtras => PolarCapSpeed_initExtras
    end type PolarCapSpeed

    type(PolarCapSpeed) u_pcap, v_pcap

    type(PolarCapQuantity) div_pcap

    type(FileCard) fcard

contains

    subroutine FlowManager_Init(filePath)
        character(*), intent(in) :: filePath

        procedure(), pointer :: p
        integer i, j

        call MsgManager_RecordSpeaker("FlowManager_Init")

        u%name = "Zonal wind"
        u%unit = "m s-1"
        call MeshManager_LinkMesh("Zonal half mesh", u)
        allocate(u%val(0:u%m%numLon+1,u%m%numLat))
        call u_pcap%stitch(u)
        call u_pcap%initExtras


        v%name = "Meridianal wind"
        v%unit = "m s-1"
        call MeshManager_LinkMesh("Meridianal half mesh", v)
        allocate(v%val(0:v%m%numLon+1,v%m%numLat))
        call v_pcap%stitch(v, "half offset")
        call v_pcap%initExtras


        div%name = "Divergence"
        div%unit = "s-1"
        call MeshManager_LinkMesh("Full mesh", div)
        allocate(div%val(0:div%m%numLon+1,div%m%numLat))
        call div_pcap%stitch(div)

        ! create output file --->
        call NFWrap_CreateSpherical2D(filePath, div%m%numLon, div%m%numLat, fcard)
        call NFWrap_Output1DVar(fcard, "lon", div%m%lon(1:div%m%numLon)*Rad2Deg)
        call NFWrap_Output1DVar(fcard, "lat", div%m%lat*Rad2Deg)
        call NFWrap_New2DVar(fcard, "u", "lon", "lat", u%name, u%unit, timeVariant=.true.)
        call NFWrap_New2DVar(fcard, "v", "lon", "lat", v%name, v%unit, timeVariant=.true.)
        call NFWrap_New2DVar(fcard, "div", "lon", "lat", div%name, div%unit, timeVariant=.true.)
        ! <---

#if (defined FC_GFORTRAN)
        p => FlowManager_Final
        call RunManager_RegisterOperation("EndRun", "FlowManager", "Final", p)
#elif (defined FC_IFORT)
        call RunManager_RegisterOperation("EndRun", "FlowManager", "Final", FlowManager_Final)
#endif

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine FlowManager_Init
    
    subroutine FlowManager_Final

        call MsgManager_RecordSpeaker("FlowManager_Final")

        call NFWrap_Close(fcard)

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker
    
    end subroutine FlowManager_Final
    
    subroutine FlowManager_Link(ui, vi)
        real(8), intent(in) :: ui(:,:), vi(:,:)

        integer i, j, k

        call MsgManager_RecordSpeaker("FlowManager_Link")

        ! normal area --->
        do j = 1, u%m%numLat
            do i = 1, u%m%numLon
                call u%val(i,j)%link(ui(i,j))
            end do
            call u%val(0,j)%mirror(u%val(u%m%numLon,j))
            call u%val(u%m%numLon+1,j)%mirror(u%val(1,j))
        end do

        do j = 1, v%m%numLat
            do i = 1, v%m%numLon
                call v%val(i,j)%link(vi(i,j))
            end do
            call v%val(0,j)%mirror(v%val(v%m%numLon,j))
            call v%val(v%m%numLon+1,j)%mirror(v%val(1,j))
        end do

        do j = 1, div%m%numLat
            do i = 1, div%m%numLon
                call div%val(i,j)%init
            end do
            call div%val(0,j)%mirror(div%val(div%m%numLon,j))
            call div%val(div%m%numLon+1,j)%mirror(div%val(1,j))
        end do
        call CalcDivergence
        ! <---

        ! polar cap --->
        do j = north, south
            k = merge(1, u%m%numLat, j == north)
            do i = 1, u_pcap%numLon
                call u_pcap%val(i,j)%mirror(u%val(i,k))
                call u_pcap%transform(i,j)%init
            end do

            do i = 1, v_pcap%numLon
                call v_pcap%val(i,j)%init
                call v_pcap%transform(i,j)%init
            end do

            k = merge(1, div%m%numLat, j == north)
            do i = 1, div_pcap%numLon
                call div_pcap%val(i,j)%mirror(div%val(i,k))
            end do
        end do
        
        call CalcPolarCapV
        call CalcPolarStereoTransformVelocity
        ! <--- 

        call FlowManager_BeforeAdvance

        call MsgManager_DeleteSpeaker

    end subroutine FlowManager_Link
    
    subroutine FlowManager_BeforeAdvance
        integer i, j

        call MsgManager_RecordSpeaker("FlowManager_BeforeAdvance")

        ! normal area --->
        do j = 1, u%m%numLat
            do i = 1, u%m%numLon
                call u%val(i,j)%save
            end do
        end do

        do j = 1, v%m%numLat
            do i = 1, v%m%numLon
                call v%val(i,j)%save
            end do
        end do

        do j = 1, div%m%numLat
            do i = 1, div%m%numLon
                call div%val(i,j)%save
            end do
        end do
        ! <---

        ! polar cap --->
        do j = north, south
            do i = 1, u_pcap%numLon
                call u_pcap%transform(i,j)%save
            end do
            do i = 1, v_pcap%numLon
                call v_pcap%val(i,j)%save
                call v_pcap%transform(i,j)%save
            end do
        end do
        ! <---
    
        call MsgManager_DeleteSpeaker

    end subroutine FlowManager_BeforeAdvance
    
    subroutine FlowManager_AfterAdvance

        call MsgManager_RecordSpeaker("FlowManager_AfterAdvance")

        call CalcDivergence
        call CalcPolarCapV
        call CalcPolarStereoTransformVelocity

        call MsgManager_DeleteSpeaker

    end subroutine FlowManager_AfterAdvance

    subroutine FlowManager_GetVelocity(x, loc, velocity, time)
        real(8), intent(in) :: x(2)
        type(Location), intent(inout) :: loc
        real(8), intent(out) :: velocity(2)
        character(*), intent(in) :: time

        integer i, j
        integer i1, i2, i3, i4, j1, j2, j3, j4
        real(8) x1(2), x2(2), x3(2), x4(2)
        real(8) u1, u2, u3, u4
        real(8) v1, v2, v3, v4
        real(8) distance, weight, weightSum, sign
        real(8) sinLon, cosLon, sinLat, sinLat2

        call MsgManager_RecordSpeaker("FlowManager_GetVelocity")
    
        ! polar cap --->
        if (loc%inPolarCap) then
            j = loc%polej
            u2 = 0.0d0
            v2 = 0.0d0
            weightSum = 0.0d0
            do i = 1, u_pcap%numLon ! v_pcap has the same numLon
                x1 = [u_pcap%lon(i),u_pcap%lat(j)]
                call MeshManager_CalcDistance(x1, x, distance)
                if (distance < eps) then
                    select case (time)
                    case ("old")
                        velocity(1) = u_pcap%val(i,j)%getOld()
                        velocity(2) = v_pcap%val(i,j)%getOld()
                    case ("half")
                        velocity(1) = (u_pcap%val(i,j)%getOld()+ &
                                       u_pcap%val(i,j)%getNew())/2.0d0
                        velocity(2) = (v_pcap%val(i,j)%getOld()+ &
                                       v_pcap%val(i,j)%getNew())/2.0d0
                    case ("new")
                        velocity(1) = u_pcap%val(i,j)%getNew()
                        velocity(2) = v_pcap%val(i,j)%getNew()
                    end select
                    call MsgManager_DeleteSpeaker
                    return
                else
                    ! inverse square distance weighted
                    weight = 1.0d0/distance
                    weightSum = weightSum+weight
                    select case (time)
                    case ("old")
                        u1 = u_pcap%transform(i,j)%getOld()
                        v1 = v_pcap%transform(i,j)%getOld()
                    case ("half")
                        u1 = (u_pcap%transform(i,j)%getOld()+ &
                              u_pcap%transform(i,j)%getNew())/2.0d0
                        v1 = (v_pcap%transform(i,j)%getOld()+ &
                              v_pcap%transform(i,j)%getNew())/2.0d0
                    case ("new")
                        u1 = u_pcap%transform(i,j)%getNew()
                        v1 = v_pcap%transform(i,j)%getNew()
                    end select
                    u2 = u2+weight*u1
                    v2 = v2+weight*v1
                end if
            end do
            u2 = u2/weightSum
            v2 = v2/weightSum
            ! transform back
            sign = merge(1.0d0, -1.0d0, j == north)
            sinLon = sin(x(1))
            cosLon = cos(x(1))
            sinLat = sin(x(2))
            sinLat2 = sinLat**2.0d0
            velocity(1) = sign*(-sinLon*u2+cosLon*v2)*sinLat
            velocity(2) = sign*(-cosLon*u2-sinLon*v2)*sinLat2
            call MsgManager_DeleteSpeaker
            return
        end if
        ! <---
 
        ! normal area --->
        i1 = loc%i(u%m%id); j1 = loc%j(u%m%id)
        i2 = i1+1; j2 = j1
        i3 = i1;   j3 = j1+1
        i4 = i1+1; j4 = j1+1
        x1 = [u%m%lon(i1),u%m%lat(j1)]
        x2 = [u%m%lon(i2),u%m%lat(j2)]
        x3 = [u%m%lon(i3),u%m%lat(j3)]
        x4 = [u%m%lon(i4),u%m%lat(j4)]
        select case (time)
        case("old")
            u1 = u%val(i1,j1)%getOld()
            u2 = u%val(i2,j2)%getOld()
            u3 = u%val(i3,j3)%getOld()
            u4 = u%val(i4,j4)%getOld()
        case("half")
            u1 = (u%val(i1,j1)%getOld()+u%val(i1,j1)%getNew())/2.0d0
            u2 = (u%val(i2,j2)%getOld()+u%val(i2,j2)%getNew())/2.0d0
            u3 = (u%val(i3,j3)%getOld()+u%val(i3,j3)%getNew())/2.0d0
            u4 = (u%val(i4,j4)%getOld()+u%val(i4,j4)%getNew())/2.0d0
        case("new")
            u1 = u%val(i1,j1)%getNew()
            u2 = u%val(i2,j2)%getNew()
            u3 = u%val(i3,j3)%getNew()
            u4 = u%val(i4,j4)%getNew()
        end select
        call BilinearInterp(x1, x2, x3, x4, u1, u2, u3, u4, x, velocity(1))
        i1 = loc%i(v%m%id); j1 = loc%j(v%m%id)
        i2 = i1+1; j2 = j1
        i3 = i1;   j3 = j1+1
        i4 = i1+1; j4 = j1+1
        x1 = [v%m%lon(i1),v%m%lat(j1)]
        x2 = [v%m%lon(i2),v%m%lat(j2)]
        x3 = [v%m%lon(i3),v%m%lat(j3)]
        x4 = [v%m%lon(i4),v%m%lat(j4)]
        select case (time)
        case("old")
            v1 = v%val(i1,j1)%getOld()
            v2 = v%val(i2,j2)%getOld()
            v3 = v%val(i3,j3)%getOld()
            v4 = v%val(i4,j4)%getOld()
        case("half")
            v1 = (v%val(i1,j1)%getOld()+v%val(i1,j1)%getNew())/2.0d0
            v2 = (v%val(i2,j2)%getOld()+v%val(i2,j2)%getNew())/2.0d0
            v3 = (v%val(i3,j3)%getOld()+v%val(i3,j3)%getNew())/2.0d0
            v4 = (v%val(i4,j4)%getOld()+v%val(i4,j4)%getNew())/2.0d0
        case("new")
            v1 = v%val(i1,j1)%getNew()
            v2 = v%val(i2,j2)%getNew()
            v3 = v%val(i3,j3)%getNew()
            v4 = v%val(i4,j4)%getNew()
        end select
        call BilinearInterp(x1, x2, x3, x4, v1, v2, v3, v4, x, velocity(2))
        ! <---

        call MsgManager_DeleteSpeaker
    
    end subroutine FlowManager_GetVelocity
    
    subroutine FlowManager_GetTransformedVelocity(x, loc, vt, time)
        real(8), intent(in) :: x(2)
        type(Location), intent(inout) :: loc
        real(8), intent(out) :: vt(2)
        character(*), intent(in) :: time

        integer i, j
        integer i1, i2, i3, i4, j1, j2, j3, j4
        real(8) x1(2), x2(2), x3(2), x4(2)
        real(8) u1, u2, u3, u4
        real(8) v1, v2, v3, v4
        real(8) velocity(2), distance, weight, weightSum, sign
        real(8) sinLon, cosLon, sinLat, sinLat2

        call MsgManager_RecordSpeaker("FlowManager_GetVelocity")
    
        ! polar cap --->
        if (loc%inPolarCap) then
            j = loc%polej
            u2 = 0.0d0
            v2 = 0.0d0
            weightSum = 0.0d0
            do i = 1, u_pcap%numLon ! v_pcap has the same numLon
                x1 = [u_pcap%lon(i),u_pcap%lat(j)]
                call MeshManager_CalcDistance(x1, x, distance)
                if (distance < eps) then
                    select case (time)
                    case ("old")
                        u1 = u_pcap%val(i,j)%getOld()
                        v1 = v_pcap%val(i,j)%getOld()
                    case ("half")
                        u1 = (u_pcap%val(i,j)%getOld()+ &
                              u_pcap%val(i,j)%getNew())/2.0d0
                        v1 = (v_pcap%val(i,j)%getOld()+ &
                              v_pcap%val(i,j)%getNew())/2.0d0
                    case ("new")
                        u1 = u_pcap%val(i,j)%getNew()
                        v1 = v_pcap%val(i,j)%getNew()
                    end select
                    sign = merge(1.0d0, -1.0d0, j == north)
                    sinLon = sin(x(1))
                    cosLon = cos(x(1))
                    sinLat = sin(x(2))
                    sinLat = sinLat**2.0d0
                    vt(1) = sign*(-sinLon/sinLat*u1-cosLon/sinLat2*v1)
                    vt(2) = sign*( cosLon/sinLat*u1-sinLon/sinLat2*v1)
                    call MsgManager_DeleteSpeaker
                    return
                else
                    ! inverse square distance weighted
                    weight = 1.0d0/distance
                    weightSum = weightSum+weight
                    select case (time)
                    case ("old")
                        u1 = u_pcap%transform(i,j)%getOld()
                        v1 = v_pcap%transform(i,j)%getOld()
                    case ("half")
                        u1 = (u_pcap%transform(i,j)%getOld()+ &
                              u_pcap%transform(i,j)%getNew())/2.0d0
                        v1 = (v_pcap%transform(i,j)%getOld()+ &
                              v_pcap%transform(i,j)%getNew())/2.0d0
                    case ("new")
                        u1 = u_pcap%transform(i,j)%getNew()
                        v1 = v_pcap%transform(i,j)%getNew()
                    end select
                    u2 = u2+weight*u1
                    v2 = v2+weight*v1
                end if
            end do
            u2 = u2/weightSum
            v2 = v2/weightSum
            vt = [u2,v2]
            call MsgManager_DeleteSpeaker
            return
        end if
        ! <---
 
        ! normal area --->
        i1 = loc%i(u%m%id); j1 = loc%j(u%m%id)
        i2 = i1+1; j2 = j1
        i3 = i1;   j3 = j1+1
        i4 = i1+1; j4 = j1+1
        x1 = [u%m%lon(i1),u%m%lat(j1)]
        x2 = [u%m%lon(i2),u%m%lat(j2)]
        x3 = [u%m%lon(i3),u%m%lat(j3)]
        x4 = [u%m%lon(i4),u%m%lat(j4)]
        select case (time)
        case("old")
            u1 = u%val(i1,j1)%getOld()
            u2 = u%val(i2,j2)%getOld()
            u3 = u%val(i3,j3)%getOld()
            u4 = u%val(i4,j4)%getOld()
        case("half")
            u1 = (u%val(i1,j1)%getOld()+u%val(i1,j1)%getNew())/2.0d0
            u2 = (u%val(i2,j2)%getOld()+u%val(i2,j2)%getNew())/2.0d0
            u3 = (u%val(i3,j3)%getOld()+u%val(i3,j3)%getNew())/2.0d0
            u4 = (u%val(i4,j4)%getOld()+u%val(i4,j4)%getNew())/2.0d0
        case("new")
            u1 = u%val(i1,j1)%getNew()
            u2 = u%val(i2,j2)%getNew()
            u3 = u%val(i3,j3)%getNew()
            u4 = u%val(i4,j4)%getNew()
        end select
        call BilinearInterp(x1, x2, x3, x4, u1, u2, u3, u4, x, velocity(1))
        i1 = loc%i(v%m%id); j1 = loc%j(v%m%id)
        i2 = i1+1; j2 = j1
        i3 = i1;   j3 = j1+1
        i4 = i1+1; j4 = j1+1
        x1 = [v%m%lon(i1),v%m%lat(j1)]
        x2 = [v%m%lon(i2),v%m%lat(j2)]
        x3 = [v%m%lon(i3),v%m%lat(j3)]
        x4 = [v%m%lon(i4),v%m%lat(j4)]
        select case (time)
        case("old")
            v1 = v%val(i1,j1)%getOld()
            v2 = v%val(i2,j2)%getOld()
            v3 = v%val(i3,j3)%getOld()
            v4 = v%val(i4,j4)%getOld()
        case("half")
            v1 = (v%val(i1,j1)%getOld()+v%val(i1,j1)%getNew())/2.0d0
            v2 = (v%val(i2,j2)%getOld()+v%val(i2,j2)%getNew())/2.0d0
            v3 = (v%val(i3,j3)%getOld()+v%val(i3,j3)%getNew())/2.0d0
            v4 = (v%val(i4,j4)%getOld()+v%val(i4,j4)%getNew())/2.0d0
        case("new")
            v1 = v%val(i1,j1)%getNew()
            v2 = v%val(i2,j2)%getNew()
            v3 = v%val(i3,j3)%getNew()
            v4 = v%val(i4,j4)%getNew()
        end select
        call BilinearInterp(x1, x2, x3, x4, v1, v2, v3, v4, x, velocity(2))
        ! <---

        sign = merge(1.0d0, -1.0d0, x(2) > 0.0d0)
        sinLon = sin(x(1))
        cosLon = cos(x(1))
        sinLat = sin(x(2))
        sinLat2 = sinLat**2.0d0
        vt(1) = sign*(-sinLon/sinLat*velocity(1)-cosLon/sinLat2*velocity(2))
        vt(2) = sign*( cosLon/sinLat*velocity(1)-sinLon/sinLat2*velocity(2))

        call MsgManager_DeleteSpeaker
    
    end subroutine FlowManager_GetTransformedVelocity

    subroutine FlowManager_GetDivergence(x, loc, odiv, time)
        real(8), intent(in) :: x(2)
        type(Location), intent(in) :: loc
        real(8), intent(out) :: odiv ! output divergence
        character(*), intent(in) :: time ! "old", "half" or "new"

        integer i, j
        integer i1, i2, i3, i4, j1, j2, j3, j4
        real(8) x1(2), x2(2), x3(2), x4(2)
        real(8) div1, div2, div3, div4
        real(8) distance, weight, weightSum
    
        call MsgManager_RecordSpeaker("FlowManager_GetDivergence")
    
        ! polar cap --->
        if (loc%inPolarCap) then
            j = loc%polej
            odiv = 0.0d0
            weightSum = 0.0d0
            do i = 1, div_pcap%numLon
                x1 = [div_pcap%lon(i),div_pcap%lat(j)]
                call MeshManager_CalcDistance(x1, x, distance)
                if (distance < eps) then
                    select case (time)
                    case ("old")
                        odiv = div_pcap%val(i,j)%getOld()
                    case ("half")
                        odiv = (div_pcap%val(i,j)%getOld() &
                               +div_pcap%val(i,j)%getNew())/2.0d0
                    case ("new")
                        odiv = div_pcap%val(i,j)%getNew()
                    end select
                    call MsgManager_DeleteSpeaker
                    return
                else
                    ! inverse square distance weighted
                    weight = 1.0d0/distance!**2.0d0
                    weightSum = weightSum+weight
                    select case (time)
                    case ("old")
                        odiv = odiv+weight*div_pcap%val(i,j)%getOld()
                    case ("half")
                        odiv = odiv+weight*(div_pcap%val(i,j)%getOld() &
                                           +div_pcap%val(i,j)%getNew())/2.0d0
                    case ("new")
                        odiv = odiv+weight*div_pcap%val(i,j)%getNew()
                    end select
                end if
            end do
            odiv = odiv/weightSum
            call MsgManager_DeleteSpeaker
            return
        end if
        ! <---

        ! normal area --->
        i1 = loc%i(div%m%id); j1 = loc%j(div%m%id)
        i2 = i1+1; j2 = j1
        i3 = i1;   j3 = j1+1
        i4 = i1+1; j4 = j1+1
        x1 = [div%m%lon(i1),div%m%lat(j1)]
        x2 = [div%m%lon(i2),div%m%lat(j2)]
        x3 = [div%m%lon(i3),div%m%lat(j3)]
        x4 = [div%m%lon(i4),div%m%lat(j4)]
        select case (time)
        case("old")
            div1 = div%val(i1,j1)%getNew()
            div2 = div%val(i2,j2)%getNew()
            div3 = div%val(i3,j3)%getNew()
            div4 = div%val(i4,j4)%getNew()
        case("half")
            div1 = (div%val(i1,j1)%getOld()+div%val(i1,j1)%getNew())/2.0d0
            div2 = (div%val(i2,j2)%getOld()+div%val(i2,j2)%getNew())/2.0d0
            div3 = (div%val(i3,j3)%getOld()+div%val(i3,j3)%getNew())/2.0d0
            div4 = (div%val(i4,j4)%getOld()+div%val(i4,j4)%getNew())/2.0d0
        case("new")
            div1 = div%val(i1,j1)%getOld()
            div2 = div%val(i2,j2)%getOld()
            div3 = div%val(i3,j3)%getOld()
            div4 = div%val(i4,j4)%getOld()
        end select
        call BilinearInterp(x1, x2, x3, x4, div1, div2, div3, div4, x, odiv)
        ! <---
    
        call MsgManager_DeleteSpeaker
    
    end subroutine FlowManager_GetDivergence
    
    subroutine FlowManager_Output(timeStep, time)
        integer, intent(in) :: timeStep
        real(8), intent(in) :: time

        real(8) fullU(div%m%numLon,div%m%numLat)
        real(8) fullV(div%m%numLon,div%m%numLat)
        real(8) fullDiv(div%m%numLon,div%m%numLat)
        integer i, j

        call MsgManager_RecordSpeaker("FlowManager_Output")
    
        do j = 1, div%m%numLat
            do i = 1, div%m%numLon
                fullU(i,j) = (u%val(i,j)%getNew()+u%val(i+1,j)%getNew())/2.0d0
                fullV(i,j) = (v%val(i,j)%getNew()+v%val(i,j+1)%getNew())/2.0d0
                fullDiv(i,j) = div%val(i,j)%getNew()
            end do
        end do

        call NFWrap_Advance(fcard, timeStep, time)
        call NFWrap_Output2DVar(fcard, "u", fullU)
        call NFWrap_Output2DVar(fcard, "v", fullV)
        call NFWrap_Output2DVar(fcard, "div", fullDiv)
    
        call MsgManager_DeleteSpeaker
    
    end subroutine FlowManager_Output
    
    subroutine CalcDivergence
        integer i, j
        real(8) u1, u2, vCosLat1, vCosLat2, dudlon, dvCosLatdlat

        do j = 1, div%m%numLat
            do i = 1, div%m%numLon
                u1 = u%val(i,j)%getNew()
                u2 = u%val(i+1,j)%getNew()
                dudlon = (u2-u1)/dimInfo%dlon
                vCosLat1 = v%val(i,j+1)%getNew()*v%m%cosLat(j+1)
                vCosLat2 = v%val(i,j)%getNew()*v%m%cosLat(j)
                dvCosLatdlat = (vCosLat2-vCosLat1)/dimInfo%dlat
                call div%val(i,j)%setNew( &
                    (dudlon+dvCosLatdlat)/Re/div%m%cosLat(j))
            end do
        end do

    end subroutine CalcDivergence

    subroutine CalcPolarCapV
        integer i, j, i1, i2
        real(8) v1, v2

        do i = 1, v%m%numLon
            i1 = i-1; i2 = i
            j = 1
            v1 = (v%val(i1,j)%getNew()+v%val(i2,j)%getNew())/2.0d0
            j = 2
            v2 = (v%val(i1,j)%getNew()+v%val(i2,j)%getNew())/2.0d0
            call v_pcap%val(i,north)%setNew((v1+v2)/2.0d0)
            j = v%m%numLat
            v1 = (v%val(i1,j)%getNew()+v%val(i2,j)%getNew())/2.0d0
            j = v%m%numLat-1
            v2 = (v%val(i1,j)%getNew()+v%val(i2,j)%getNew())/2.0d0
            call v_pcap%val(i,south)%setNew((v1+v2)/2.0d0)
        end do
    
    end subroutine CalcPolarCapV
   
    subroutine PolarCapSpeed_initExtras(pcap)
        class(PolarCapSpeed), intent(inout) :: pcap

        integer i, j

        allocate(pcap%sinLon(pcap%numLon))
        allocate(pcap%cosLon(pcap%numLon))
        allocate(pcap%transform(pcap%numLon,north:south))

        do i = 1, pcap%numLon
            pcap%sinLon(i) = sin(pcap%lon(i))
            pcap%cosLon(i) = cos(pcap%lon(i))
        end do

        do j = north, south
            pcap%sinLat(j) = sin(pcap%lat(j))
            pcap%sinLat2(j) = sin(pcap%lat(j))**2.0d0
        end do

    end subroutine PolarCapSpeed_initExtras
    
    subroutine CalcPolarStereoTransformVelocity
        integer i, j

        j = north
        do i = 1, u_pcap%numLon
            call u_pcap%transform(i,j)%setNew( &
                -u_pcap%sinLon(i)/u_pcap%sinLat (j)*u_pcap%val(i,j)%getNew() &
                -u_pcap%cosLon(i)/u_pcap%sinLat2(j)*v_pcap%val(i,j)%getNew())
            call v_pcap%transform(i,j)%setNew( &
                 v_pcap%cosLon(i)/v_pcap%sinLat (j)*u_pcap%val(i,j)%getNew() &
                -v_pcap%sinLon(i)/v_pcap%sinLat2(j)*v_pcap%val(i,j)%getNew())
        end do
        j = south
        do i = 1, u_pcap%numLon
            call u_pcap%transform(i,j)%setNew( &
                 u_pcap%sinLon(i)/u_pcap%sinLat (j)*u_pcap%val(i,j)%getNew() &
                +u_pcap%cosLon(i)/u_pcap%sinLat2(j)*v_pcap%val(i,j)%getNew())
            call v_pcap%transform(i,j)%setNew( &
                -v_pcap%cosLon(i)/v_pcap%sinLat (j)*u_pcap%val(i,j)%getNew() &
                +v_pcap%sinLon(i)/v_pcap%sinLat2(j)*v_pcap%val(i,j)%getNew())
        end do
    
    end subroutine CalcPolarStereoTransformVelocity
 
    subroutine PolarStereoTransformVelocity(x0, v0, v1, polej)
        real(8), intent(in) :: x0(2), v0(2)
        real(8), intent(out) :: v1(2)
        integer, intent(in) :: polej

        real(8) sign, sinLon, cosLon, sinLat, sinLat2

        sign = merge(1.0d0, -1.0d0, polej == north)
        sinLon = sin(x0(1))
        cosLon = cos(x0(1))
        sinLat = sin(x0(2))
        sinLat2 = sinLat**2.0d0
        v1(1) = sign*(-sinLon/sinLat*v0(1)-cosLon/sinLat2*v0(2))
        v1(2) = sign*( cosLon/sinLat*v0(1)-sinLon/sinLat2*v0(2))

    end subroutine PolarStereoTransformVelocity

    subroutine BilinearInterp(x1, x2, x3, x4, f1, f2, f3, f4, x, f)
        real(8), intent(in) :: x1(2), x2(2), x3(2), x4(2)
        real(8), intent(in) :: f1, f2, f3, f4
        real(8), intent(in) :: x(2)
        real(8), intent(out) :: f

        real(8) xb1, xb2, yb1, yb2, xx, yy
        real(8) a, b, c, d

        xb1 = x1(1)
        xb2 = x2(1)
        yb1 = x1(2)
        yb2 = x3(2)

        xx = (x(1)-xb1)/(xb2-xb1)
        yy = (x(2)-yb1)/(yb2-yb1)

        a = f1
        b = f2-f1
        c = f3-f1
        d = f1-f2-f3+f4

        f = a+b*xx+c*yy+d*xx*yy

    end subroutine BilinearInterp

end module FlowManager
