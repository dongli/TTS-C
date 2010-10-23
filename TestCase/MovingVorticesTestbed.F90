module MovingVorticesTestbed

    use MsgManager
    use FloatingPoint
    use SphereService
    use RunManager
    use MeshManager

    implicit none

    private

    public MovingVorticesTestbed_Init
    public MovingVorticesTestbed_Advance
    public MovingVorticesTestbed_Final
    public MovingVorticesTestbed_CalcExactSolution
    public u, v, lonVortex, latVortex

    real(8), parameter :: gamma = 5.0d0 ! control the "stiffness" of the field
    real(8), parameter :: V0 = PI2*Re/12.0d0/86400.0d0
    real(8), parameter :: V0_Re = V0/Re
    real(8), parameter :: alpha  = PI05  ! angle between rotation axis and pole
    real(8), parameter :: rho0 = 3.0d0
    real(8), parameter :: lonRotate0 = PI
    real(8), parameter :: latRotate0 = PI05-alpha
    real(8), parameter :: lonVortex0 = PI05
    real(8), parameter :: latVortex0 = Equator

    real(8) lonVR0, latVR0
    real(8) lonVortex, latVortex

    integer numLonU, numLatU

    real(8), allocatable :: lonU(:), latU(:)
    real(8), allocatable :: rotateU(:,:)
    real(8), allocatable :: u(:,:)

    integer numLonV, numLatV

    real(8), allocatable :: lonV(:), latV(:)
    real(8), allocatable :: rotateV(:,:)
    real(8), allocatable :: v(:,:)

    interface MovingVorticesTestbed_CalcExactSolution
        module procedure CalcExactSolution_1
        module procedure CalcExactSolution_2
    end interface

contains

    subroutine MovingVorticesTestbed_Init
        type(Mesh), pointer :: m1, m2
        integer i, j

        procedure(), pointer :: p

        call MsgManager_RecordSpeaker("MovingVorticesTestbed_Init")

        m1 => MeshSelector("Zonal half mesh")
        m2 => MeshSelector("Meridianal half mesh")

        numLonU = m1%numLon
        numLatU = m1%numLat
        numLonV = m2%numLon
        numLatV = m2%numLat

        allocate(lonU(numLonU))
        allocate(latU(numLatU))
        allocate(rotateU(numLonU,numLatU))
        allocate(u(numLonU,numLatU))
        allocate(lonV(numLonV))
        allocate(latV(numLatV))
        allocate(rotateV(numLonV,numLatV))
        allocate(v(numLonV,numLatV))

        do i = 1, numLonU
            lonU(i) = m1%lon(i)
        end do
        do j = 1, numLatU
            latU(j) = m1%lat(j)
        end do

        do i = 1, numLonV
            lonV(i) = m2%lon(i)
        end do
        do j = 1, numLatV
            latV(j) = m2%lat(j)
        end do

        ! set background rotating velocity
        do j = 1, numLatU
            do i = 1, numLonU
                rotateU(i,j) = V0*(cos(latU(j))*cos(alpha)+sin(latU(j))*cos(lonU(i))*sin(alpha))
            end do
        end do
        do j = 1, numLatV
            do i = 1, numLonV
                rotateV(i,j) = -V0*sin(lonV(i))*sin(alpha)
            end do
        end do

        ! set the initial coordinate of vortex (1/2) in the rotated coordinate system
        call RotationTransform( &
            lonRotate0, latRotate0, lonVortex0, latVortex0, lonVR0, latVR0)
        lonVortex = lonVortex0
        latVortex = latVortex0

        call MovingVorticesTestbed_Advance(0.0d0)

        p => MovingVorticesTestbed_Final
        call RunManager_RegisterOperation("EndRun", &
            "MovingVorticesTestbed", "MovingVorticesTestbed_Final", p)

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine MovingVorticesTestbed_Init

    subroutine MovingVorticesTestbed_Advance(elapseTime)
        real(8), intent(in) :: elapseTime

        real(8) lon, dlon, latRotate, vortexU, vortexV
        integer i, j

        ! calculate the coordinate of "north" vortex center
        ! * coordinate in the rotated coordinate system
        lon = lonVR0+V0_Re*elapseTime
        if (lon > PI2) lon = lon-PI2
        ! * coordinate in the original coordinate system
        call InverseRotationTransform( &
            lonRotate0, latRotate0, lonVortex, latVortex, lon, latVR0)

        ! calculate the vortex velocity
        do j = 1, numLatU
            do i = 1, numLonU
                dlon = lonU(i)-lonVortex
                call RotationTransform( &
                    lonVortex, latVortex, lonU(i), latU(j), latR=latRotate)
                vortexU = ReOmega(latRotate)*(sin(latVortex)*cos(latU(j))- &
                    cos(latVortex)*cos(dlon)*sin(latU(j)))
                u(i,j) = rotateU(i,j)+vortexU
            end do
        end do
        do j = 1, numLatV
            do i = 1, numLonV
                dlon = lonV(i)-lonVortex
                call RotationTransform( &
                    lonVortex, latVortex, lonV(i), latV(j), latR=latRotate)
                vortexV = ReOmega(latRotate)*cos(latVortex)*sin(dlon)
                v(i,j) = rotateV(i,j)+vortexV
            end do
        end do

    end subroutine MovingVorticesTestbed_Advance

    subroutine MovingVorticesTestbed_Final

        call MsgManager_RecordSpeaker("MovingVorticesTestbed_Final")

        deallocate(lonU)
        deallocate(latU)
        deallocate(rotateU)
        deallocate(u)
        deallocate(lonV)
        deallocate(latV)
        deallocate(rotateV)
        deallocate(v)

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine MovingVorticesTestbed_Final

    subroutine CalcExactSolution_1(n, lon, lat, t, q)
        integer, intent(in) :: n
        real(8), intent(in) :: lon(n), lat(n), t
        real(8), intent(out) :: q(n)

        logical, save :: firstCall = .true.
        logical, save :: reverseFlag1 = .false.
        logical, save :: reverseFlag2 = .true.
        real(8), save :: prevLonVortex
        real(8) lonR, latR
        integer i

        call MsgManager_RecordSpeaker("CalcExactSolution_1")

        if (firstCall) then
            firstCall = .false.
        else
            if (abs(prevLonVortex-lonVortex) > PI*0.99d0) then
                call MsgManager_Speak(Notice, "Doing reverse.")
                reverseFlag1 = .not. reverseFlag2
            end if
        end if

        do i = 1, n
            call RotationTransform( &
                lonVortex, latVortex, lon(i), lat(i), lonR, latR)
            if (reverseFlag1 .eqv. reverseFlag2) then
                lonR = lonR+PI
                if (lonR > PI2) lonR = lonR-PI2
            end if
            q(i) = 1.0d0-tanh(rho(latR)/gamma*sin(lonR-omega(latR)*t))
        end do

        prevLonVortex = lonVortex

        call MsgManager_DeleteSpeaker

    end subroutine CalcExactSolution_1

    subroutine CalcExactSolution_2(nlon, nlat, lon, lat, t, q)
        integer, intent(in) :: nlon, nlat
        real(8), intent(in) :: lon(nlon), lat(nlat), t
        real(8), intent(out) :: q(nlon,nlat)

        logical, save :: firstCall = .true.
        logical, save :: reverseFlag = .false.
        real(8), save :: prevLonVortex
        real(8) lonR, latR
        integer i, j

        call MsgManager_RecordSpeaker("CalcExactSolution_2")

        if (firstCall) then
            firstCall = .false.
        else
            if (abs(prevLonVortex-lonVortex) > PI*0.99d0) then
                call MsgManager_Speak(Notice, "Doing reverse.")
                reverseFlag = .not. reverseFlag
            end if
        end if

        do j = 1, nlat
            do i = 1, nlon
                call RotationTransform( &
                   lonVortex, latVortex, lon(i), lat(j), lonR, latR)
                if (reverseFlag) then
                    lonR = lonR+PI
                    if (lonR > PI2) lonR = lonR-PI2
                end if
                q(i,j) = 1.0d0-tanh(rho(latR)/gamma*sin(lonR-omega(latR)*t))
            end do
        end do

        prevLonVortex = lonVortex

        call MsgManager_DeleteSpeaker

    end subroutine CalcExactSolution_2

    real(8) function rho(lat)
        real(8), intent(in) :: lat

        rho = rho0*cos(lat)

    end function rho

    real(8) function omega(latRotate)
        real(8), intent(in) :: latRotate

        real(8), parameter :: fac = 1.5d0*sqrt(3.0d0)
        real(8) r ! radial distance of vortex
        real(8) V ! tangential velocity

        r = rho(latRotate)
        if (abs(r) < eps) then
            omega = 0.0d0
        else
            V = V0*fac*tanh(r)/cosh(r)**2.0d0
            omega = V/r/Re
        end if

    end function omega

    real(8) function ReOmega(latRotate)
        real(8), intent(in) :: latRotate ! rotated latitude

        real(8), parameter :: fac = 1.5d0*sqrt(3.0d0)
        real(8) r ! radial distance of vortex
        real(8) V ! tangential velocity

        r = rho(latRotate)
        if (abs(r) < eps) then
            ReOmega = 0.0d0
        else
            V = V0*fac*tanh(r)/cosh(r)**2.0d0
            ReOmega = V/r
        end if

    end function ReOmega

end module MovingVorticesTestbed
