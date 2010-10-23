module SolidRotationTestbed

    use MsgManager
    use RunManager
    use SphereService
    use MeshManager

    implicit none

    private

    public SolidRotationTestbed_Init
    public SolidRotationTestbed_Final
    public SolidRotationTestbed_CalcExactSolution
    public u, v

    ! parameters for rotation
    real(8), parameter :: omega = PI2/12.0d0/86400.0d0
    real(8), parameter :: V0 = omega*Re
    real(8), parameter :: alpha  = PI05  ! angle between rotation axis and pole
    real(8), parameter :: lonRotate0 = PI
    real(8), parameter :: latRotate0 = PI05-alpha

    ! parameters for cosine hill
    real(8), parameter :: C0(2) = [PI05,PI05] ! initial center of cosine hill
    real(8), parameter :: R = Re/3.0d0        ! radius of cosine hill
    real(8), parameter :: PHI0 = 1000.0d0     ! amplitude of cosine hill

    integer numLonU, numLatU

    real(8), allocatable :: lonU(:), latU(:)
    real(8), allocatable :: u(:,:)

    integer numLonV, numLatV

    real(8), allocatable :: lonV(:), latV(:)
    real(8), allocatable :: v(:,:)

    real(8) CR0(2) ! rotated C0

    interface SolidRotationTestbed_CalcExactSolution
        module procedure CalcExactSolution_1
        module procedure CalcExactSolution_2
    end interface

contains

    subroutine SolidRotationTestbed_Init
        type(Mesh), pointer :: m1, m2
        integer i, j

        procedure(), pointer :: p

        call MsgManager_RecordSpeaker("SolidRotationTestbed_Init")

        m1 => MeshSelector("Zonal half mesh")
        m2 => MeshSelector("Meridianal half mesh")

        numLonU = m1%numLon
        numLatU = m1%numLat
        numLonV = m2%numLon
        numLatV = m2%numLat

        allocate(lonU(numLonU))
        allocate(latU(numLatU))
        allocate(u(numLonU,numLatU))
        allocate(lonV(numLonV))
        allocate(latV(numLatV))
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
        do j = 1, numLatU
            do i = 1, numLonU
                u(i,j) = V0*(cos(latU(j))*cos(alpha)+sin(latU(j))*cos(lonU(i))*sin(alpha))
            end do
        end do
        do j = 1, numLatV
            do i = 1, numLonV
                v(i,j) = -V0*sin(lonV(i))*sin(alpha)
            end do
        end do

        call RotationTransform(lonRotate0, latRotate0, C0(1), C0(2), CR0(1), CR0(2))

        p => SolidRotationTestbed_Final
        call RunManager_RegisterOperation("EndRun", "SolidRotationTestbed", "Final", p)

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine SolidRotationTestbed_Init

    subroutine SolidRotationTestbed_Final

        call MsgManager_RecordSpeaker("SolidRotationTestbed_Final")

        deallocate(lonU)
        deallocate(latU)
        deallocate(u)
        deallocate(lonV)
        deallocate(latV)
        deallocate(v)

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine SolidRotationTestbed_Final

    subroutine CalcExactSolution_1(n, lon, lat, t, q)
        integer, intent(in) :: n
        real(8), intent(in) :: lon(n), lat(n), t
        real(8), intent(out) :: q(n)

        real(8) lonC, latC, lonCR
        integer i

        lonCR = CR0(1)+omega*t
        if (lonCR > PI2) lonCR = lonCR-PI2

        call InverseRotationTransform(lonRotate0, latRotate0, lonC, latC, lonCR, CR0(2))

        do i = 1, n
            q(i) = cosinehill(lonC, latC, lon(i), lat(i))
        end do

    end subroutine CalcExactSolution_1

    subroutine CalcExactSolution_2(nlon, nlat, lon, lat, t, q)
        integer, intent(in) :: nlon, nlat
        real(8), intent(in) :: lon(nlon), lat(nlat), t
        real(8), intent(out) :: q(nlon,nlat)

        real(8) lonC, latC, lonCR, lonR, latR
        integer i, j

        lonCR = CR0(1)+omega*t
        if (lonCR > PI2) lonCR = lonCR-PI2

        call InverseRotationTransform(lonRotate0, latRotate0, lonC, latC, lonCR, CR0(2))

        do j = 1, nlat
            do i = 1, nlon
                q(i,j) = cosinehill(lonC, latC, lon(i), lat(j))
            end do
        end do

    end subroutine CalcExactSolution_2

    real(8) function cosinehill(lon0, lat0, lon, lat) result(res)
        real(8), intent(in) :: lon0, lat0 ! coordinate of cosine hill center
        real(8), intent(in) :: lon, lat   ! 

        real(8) d

        call MeshManager_CalcDistance([lon0,lat0], [lon,lat], d)

        if (d < R) then
            res = PHI0*(1.0+cos(PI*d/R))/2.0d0
        else
            res = 0.0d0
        end if

    end function cosinehill

end module SolidRotationTestbed
