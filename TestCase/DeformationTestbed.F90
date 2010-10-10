! *************************************************************************** !
! DeformationTestbed module                                                   !
!                                                                             !
! Description:                                                                !
!   This is the new class of deformational flow test cases. Four test cases   !
!   have  been  proposed  in Nair and Lauritzen (2010). The  flow fields of   !
!   "case 1", "case 2" and "case 4" are non-divergent.                        !
!                                                                             !
! Reference:                                                                  !
!   - Nair, R. D. and Lauritzen, Peter H., 2010: A Class of Deformational     !
!     Flow Test Cases for Linear Transport Problems, Journal of Computational !
!     Physics, 229, 8868-8887.                                                !
!                                                                             !
! Author:                                                                     !
!   DONG Li, dongli@lasg.iap.ac.cn                                            !
!                                                                             !
! *************************************************************************** !

module DeformationTestbed

    use MsgManager
    use RunManager
    use SphereService
    use MeshManager

    implicit none

    private

    public DeformationTestbed_Init
    public DeformationTestbed_Advance
    public DeformationTestbed_Final
    public DeformationTestbed_CalcInitCond
    public u, v

    character(20) caseId, initId

    real(8) T    ! time scale of one cycle
    real(8) lon1 ! the initial location
    real(8) lat1 ! the initial location
    real(8) lon2 ! the initial location
    real(8) lat2 ! the initial location

    real(8) b, c, hmax, b0 ! parameters for initial conditions
    real(8) r
    real(8) k, k05

    integer numLonU, numLatU

    real(8), allocatable :: lonU(:), latU(:)
    real(8), allocatable :: sin2LonU05(:), sin2LonU(:)
    real(8), allocatable :: sinLatU2(:), cosLatU(:)
    real(8), allocatable :: u(:,:)

    integer numLonV, numLatV

    real(8), allocatable :: lonV(:), latV(:)
    real(8), allocatable :: sinLonV(:), sinLonV2(:)
    real(8), allocatable :: cosLatV(:)
    real(8), allocatable :: v(:,:)

    ! rotated spherical coordinate system for crossing-pole motion
    real(8), parameter :: lonP = 0.0d0, latP = 0.0d0
    real(8), parameter :: sinLatP = sin(latP), cosLatP = cos(latP)
    real(8), allocatable :: lonUR(:,:), latUR(:,:)
    real(8), allocatable :: sinLonMinusLonPU(:), cosLonMinusLonPU(:), sinLatU(:)
    real(8), allocatable :: sin2LonU05R(:,:), sin2LonUR(:,:), sinLonUR(:,:), cosLonUR(:,:)
    real(8), allocatable :: sinLatU2R(:,:), sinLatUR(:,:), cosLatUR(:,:), cos2LatUR(:,:), cos3LatUR(:,:)
    real(8), allocatable :: lonVR(:,:), latVR(:,:)
    real(8), allocatable :: sin2LonV05R(:,:), sin2LonVR(:,:), sinLonVR(:,:), cosLonVR(:,:)
    real(8), allocatable :: sinLatV2R(:,:), sinLatVR(:,:), cosLatVR(:,:), cos2LatVR(:,:), cos3LatVR(:,:)

contains

    subroutine DeformationTestbed_Init(caseId_in, initId_in)
        character(*), intent(in) :: caseId_in
        character(*), intent(in) :: initId_in

        type(Mesh), pointer :: m1, m2
        procedure(), pointer :: p
        integer i, j

        call MsgManager_RecordSpeaker("DeformationTestbed_Init")

        caseId = caseId_in
        initId = initId_in

        T = 5.0d0
        r = Re/2.0d0

        select case (caseId)
        case ("case 1")
            k = 2.4d0
            k05 = k*0.5d0
            lon1 = PI
            lat1 = PI/3.0d0
            lon2 = PI
            lat2 = -PI/3.0d0
        case ("case 2")
            k = 2.0d0
            lon1 = 5.0d0/6.0d0*PI
            lat1 = 0.0d0
            lon2 = 7.0d0/6.0d0*PI
            lat2 = 0.0d0
        case ("case 3")
            k = 1.0d0
            k05 = k*0.5d0
            lon1 = PI05
            lat1 = PI/3.0d0
            lon2 = PI*1.5d0
            lat2 = PI/3.0d0
        case ("case 4")
            k = 2.0d0
            lon1 = 5.0d0/6.0d0*PI
            lat1 = 0.0d0
            lon2 = 7.0d0/6.0d0*PI
            lat2 = 0.0d0
        case default
            call MsgManager_Speak(Error, &
                "Unknow test case "//trim(caseId)//".")
            call RunManager_EndRun
        end select

        select case (initId)
        case ("cosine bells")
            hmax = 1.0d0
            b = 0.1d0
            c = 0.9d0
        case ("Gaussian hills")
            hmax = 1.0d0 ! hill height
            b0 = 5.0d0   ! hill width
        case ("slotted cylinders")
            b = 0.1d0
            c = 1.0d0
        case default
            call MsgManager_Speak(Error, &
                "Unknow initial condition "//trim(caseId)//".")
            call RunManager_EndRun
        end select

        m1 => MeshSelector("Zonal half mesh")
        m2 => MeshSelector("Meridianal half mesh")

        numLonU = m1%numLon
        numLatU = m1%numLat
        numLonV = m2%numLon
        numLatV = m2%numLat

        allocate(lonU(numLonU))
        allocate(sin2LonU05(numLonU))
        allocate(sin2LonU(numLonU))
        allocate(latU(numLatU))
        allocate(sinLatU2(numLatU))
        allocate(cosLatU(numLatU))
        allocate(u(numLonU,numLatU))
        allocate(lonV(numLonV))
        allocate(sinLonV(numLonV))
        allocate(sinLonV2(numLonV))
        allocate(latV(numLatV))
        allocate(cosLatV(numLatV))
        allocate(v(numLonV,numLatV))

        do i = 1, numLonU
            lonU(i) = m1%lon(i)
            sin2LonU05(i) = sin(lonU(i)*0.5d0)**2.0d0
            sin2LonU(i) = sin(lonU(i))**2.0d0
        end do
        do j = 1, numLatU
            latU(j) = m1%lat(j)
            sinLatU2(j) = sin(latU(j)*2.0d0)
            cosLatU(j) = cos(latU(j))
        end do

        do i = 1, numLonV
            lonV(i) = m2%lon(i)
            sinLonV(i) = sin(lonV(i))
            sinLonV2(i) = sin(lonV(i)*2.0d0)
        end do
        do j = 1, numLatV
            latV(j) = m2%lat(j)
            cosLatV(j) = cos(latV(j))
        end do

        ! for crossing-pole motion
        allocate(lonUR(numLonU,numLatU))
        allocate(sinLonMinusLonPU(numLonU))
        allocate(cosLonMinusLonPU(numLonU))
        allocate(sinLatU(numLatU))
        allocate(sin2LonU05R(numLonU,numLatU))
        allocate(sin2LonUR(numLonU,numLatU))
        allocate(sinLonUR(numLonU,numLatU))
        allocate(cosLonUR(numLonU,numLatU))
        allocate(latUR(numLonU,numLatU))
        allocate(sinLatU2R(numLonU,numLatU))
        allocate(sinLatUR(numLonU,numLatU))
        allocate(cosLatUR(numLonU,numLatU))
        allocate(cos2LatUR(numLonU,numLatU))
        allocate(cos3LatUR(numLonU,numLatU))
        allocate(lonVR(numLonV,numLatV))
        allocate(sin2LonV05R(numLonV,numLatV))
        allocate(sin2LonVR(numLonV,numLatV))
        allocate(sinLonVR(numLonV,numLatV))
        allocate(cosLonVR(numLonV,numLatV))
        allocate(latVR(numLonV,numLatV))
        allocate(sinLatV2R(numLonV,numLatV))
        allocate(sinLatVR(numLonV,numLatV))
        allocate(cosLatVR(numLonV,numLatV))
        allocate(cos2LatVR(numLonV,numLatV))
        allocate(cos3LatVR(numLonV,numLatV))

        do i = 1, numLonU
            sinLonMinusLonPU(i) = sin(lonU(i)-lonP)
            cosLonMinusLonPU(i) = cos(lonU(i)-lonP)
        end do
        do j = 1, numLatU
            sinLatU(j) = sin(latU(j))
        end do
        do j = 1, numLatU
            do i = 1, numLonU
                call RotationTransform(lonP, latP, lonU(i), latU(j), lonUR(i,j), latUR(i,j))
                sin2LonU05R(i,j) = sin(lonUR(i,j)*0.5d0)**2.0d0
                sin2LonUR(i,j) = sin(lonUR(i,j))**2.0d0
                sinLonUR(i,j) = sin(lonUR(i,j))
                cosLonUR(i,j) = cos(lonUR(i,j))
                sinLatU2R(i,j) = sin(latUR(i,j)*2.0d0)
                sinLatUR(i,j) = sin(latUR(i,j))
                cosLatUR(i,j) = cos(latUR(i,j))
                cos2LatUR(i,j) = cos(latUR(i,j))**2.0d0
                cos3LatUR(i,j) = cos(latUR(i,j))**3.0d0
            end do
        end do

        do j = 1, numLatV
            do i = 1, numLonV
                call RotationTransform(lonP, latP, lonV(i), latV(j), lonVR(i,j), latVR(i,j))
                sin2LonV05R(i,j) = sin(lonVR(i,j)*0.5d0)**2.0d0
                sin2LonVR(i,j) = sin(lonVR(i,j))**2.0d0
                sinLonVR(i,j) = sin(lonVR(i,j))
                cosLonVR(i,j) = cos(lonVR(i,j))
                sinLatV2R(i,j) = sin(latVR(i,j)*2.0d0)
                sinLatVR(i,j) = sin(latVR(i,j))
                cosLatVR(i,j) = cos(latVR(i,j))
                cos2LatVR(i,j) = cos(latVR(i,j))**2.0d0
                cos3LatVR(i,j) = cos(latVR(i,j))**3.0d0
            end do
        end do

        call DeformationTestbed_Advance(0.0d0)

        p => DeformationTestbed_Final
        call RunManager_RegisterOperation("EndRun", &
            "DeformationTestbed", "Final", p)

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker
    
    end subroutine DeformationTestbed_Init
    
    subroutine DeformationTestbed_Final
    
        call MsgManager_RecordSpeaker("DeformationTestbed_Final")

        deallocate(lonU)
        deallocate(sin2LonU05)
        deallocate(sin2LonU)
        deallocate(latU)
        deallocate(sinLatU2)
        deallocate(cosLatU)
        deallocate(u)
        deallocate(lonV)
        deallocate(sinLonV)
        deallocate(sinLonV2)
        deallocate(latV)
        deallocate(cosLatV)
        deallocate(v)

        ! for crossing-pole motion
        deallocate(lonUR)
        deallocate(sinLonMinusLonPU)
        deallocate(cosLonMinusLonPU)
        deallocate(sinLatU)
        deallocate(sin2LonU05R)
        deallocate(sin2LonUR)
        deallocate(sinLonUR)
        deallocate(cosLonUR)
        deallocate(latUR)
        deallocate(sinLatU2R)
        deallocate(sinLatUR)
        deallocate(cosLatUR)
        deallocate(cos2LatUR)
        deallocate(cos3LatUR)
        deallocate(lonVR)
        deallocate(sin2LonV05R)
        deallocate(sin2LonVR)
        deallocate(sinLonVR)
        deallocate(cosLonVR)
        deallocate(latVR)
        deallocate(sinLatV2R)
        deallocate(sinLatVR)
        deallocate(cosLatVR)
        deallocate(cos2LatVR)
        deallocate(cos3LatVR)

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine DeformationTestbed_Final
    
    subroutine DeformationTestbed_Advance(elapseTime)
        real(8), intent(in) :: elapseTime

        real(8) cost, lon
        real(8) uR, vR ! velocity in rotated spherical coordinate system
        real(8) tmp1, tmp2, tmp3, tmp4, vtmp
        integer i, j

        cost = cos(PI*elapseTime/T)
        
        select case (caseId)
        case ("case 1")
            do j = 1, numLatU
                do i = 1, numLonU
                    u(i,j) = k*sin2LonU05(i)*sinLatU2(j)*cost
                end do
            end do
            do j = 1, numLatV
                do i = 1, numLonV
                    v(i,j) = k05*sinLonV(i)*cosLatV(j)*cost
                end do
            end do
        case ("case 2")
            do j = 1, numLatU
                do i = 1, numLonU
                    u(i,j) = k*sin2LonU(i)*sinLatU2(j)*cost
                end do
            end do
            do j = 1, numLatV
                do i = 1, numLonV
                    v(i,j) = k*sinLonV2(i)*cosLatV(j)*cost
                end do
            end do
        case ("case 3")
            do j = 1, numLatU
                do i = 1, numLonU
                    uR = -k*sin2LonU05R(i,j)*sinLatU2R(i,j)*cos2LatUR(i,j)*cost
                    vR = k05*sinLonUR(i,j)*cos3LatUR(i,j)*cost
                    tmp1 = sinLatP*cosLatUR(i,j)
                    tmp2 = cosLatP*sinLatUR(i,j)*cosLonUR(i,j)
                    tmp3 = tmp1+tmp2
                    tmp4 = cosLatP*sinLonUR(i,j)
                    vtmp = tmp3*vR+tmp4*uR
                    vtmp = vtmp/cosLatU(j)
                    tmp1 = sinLatP*cosLatU(j)
                    tmp2 = cosLatP*sinLatU(j)*cosLonMinusLonPU(i)
                    tmp3 = tmp1-tmp2
                    u(i,j) = tmp3*vtmp-cosLatUR(i,j)*vR
                    u(i,j) = u(i,j)/cosLatP/sinLonMinusLonPU(i)
                end do
            end do
            do j = 1, numLatV
                do i = 1, numLonV
                    uR = -k*sin2LonV05R(i,j)*sinLatV2R(i,j)*cos2LatVR(i,j)*cost
                    vR = k05*sinLonVR(i,j)*cos3LatVR(i,j)*cost
                    tmp1 = sinLatP*cosLatVR(i,j)
                    tmp2 = cosLatP*sinLatVR(i,j)*cosLonVR(i,j)
                    tmp3 = tmp1+tmp2
                    tmp4 = cosLatP*sinLonVR(i,j)
                    v(i,j) = tmp3*vR+tmp4*uR
                    v(i,j) = v(i,j)/cosLatV(j)
                end do
            end do
        case ("case 4")
            do j = 1, numLatU
                do i = 1, numLonU
                    lon = lonU(i)-PI2*elapseTime/T
                    u(i,j) = k*sin(lon)**2.0d0*sinLatU2(j)*cost+PI2*Re*cosLatU(j)/T
                end do
            end do
            do j = 1, numLatV
                do i = 1, numLonV
                    lon = lonV(i)-PI2*elapseTime/T
                    v(i,j) = k*sin(lon*2.0d0)*cosLatV(j)*cost
                end do
            end do
        end select

    end subroutine DeformationTestbed_Advance
    
    subroutine DeformationTestbed_CalcInitCond(n, lon, lat, q)
        integer, intent(in) :: n
        real(8), intent(in) :: lon(n), lat(n)
        real(8), intent(out) :: q(n)

        real(8) r1, r2
        real(8) x, y, z, x1, y1, z1, x2, y2, z2
        real(8) dlon1, dlat1, dlon2, dlat2, q1, q2
        integer i

        select case (initId)
        case ("cosine bells")
            do i = 1, n
                call MeshManager_CalcDistance([lon(i),lat(i)], [lon1,lat1], r1)
                call MeshManager_CalcDistance([lon(i),lat(i)], [lon2,lat2], r2)
                if (r1 < r) then
                    q(i) = b+c*hmax*0.5d0*(1.0d0+cos(PI*r1/r))
                else if (r2 < r) then
                    q(i) = b+c*hmax*0.5d0*(1.0d0+cos(PI*r2/r))
                else
                    q(i) = b
                end if
            end do
        case ("Gaussian hills")
            x1 = Re*cos(lat1)*cos(lon1)
            y1 = Re*cos(lat1)*sin(lon1)
            z1 = Re*sin(lat1)
            x2 = Re*cos(lat2)*cos(lon2)
            y2 = Re*cos(lat2)*sin(lon2)
            z2 = Re*sin(lat2)
            do i = 1, n
                x = Re*cos(lat(i))*cos(lon(i))
                y = Re*cos(lat(i))*sin(lon(i))
                z = Re*sin(lat(i))
                q1 = hmax*exp(-b0*((x-x1)**2.0d0+(y-y1)**2.0d0+(z-z1)**2.0d0))
                q2 = hmax*exp(-b0*((x-x2)**2.0d0+(y-y2)**2.0d0+(z-z2)**2.0d0))
                q(i) = q1+q2
            end do
        case ("slotted cylinders")
            ! see eqn. 17 in Nair et al. (2010)
            do i = 1, n
                call MeshManager_CalcDistance([lon(i),lat(i)], [lon1,lat1], r1)
                call MeshManager_CalcDistance([lon(i),lat(i)], [lon2,lat2], r2)
                dlon1 = abs(lon(i)-lon1)
                dlon2 = abs(lon(i)-lon2)
                dlat1 = lat(i)-lat1
                dlat2 = lat(i)-lat2
                q(i) = b
                if (r1 <= r .and. dlon1 >= r/6.0d0) then 
                    q(i) = c
                end if
                if (r1 <= r .and. dlon1 < r/6.0d0 .and. dlat1 < -5.0d0/12.0d0*r) then
                    q(i) = c
                end if
                if (r2 <= r .and. dlon2 >= r/6.0d0) then
                    q(i) = c
                end if
                if (r2 <= r .and. dlon2 < r/6.0d0 .and. dlat2 > -5.0d0/12.0d0*r) then
                    q(i) = c
                end if
            end do
        end select

    end subroutine DeformationTestbed_CalcInitCond

end module DeformationTestbed

