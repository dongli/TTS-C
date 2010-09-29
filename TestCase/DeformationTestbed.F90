! *************************************************************************** !
! DeformationTestbed module                                                   !
!                                                                             !
! Description:                                                                !
!                                                                             !
! Author:                                                                     !
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

    integer caseId, initId
    real(8) b
    real(8) c
    real(8) lon1
    real(8) lat1
    real(8) lon2
    real(8) lat2
    real(8) r
    real(8) k, k05
    real(8) T

    integer numLonU, numLatU

    real(8), allocatable :: lonU(:), latU(:)
    real(8), allocatable :: sin2LonU05(:), sinLonU2(:)
    real(8), allocatable :: sinLatU2(:), cosLatU(:), cosLatU2(:)
    real(8), allocatable :: u(:,:)

    integer numLonV, numLatV

    real(8), allocatable :: lonV(:), latV(:)
    real(8), allocatable :: sinLonV(:), sin2LonV(:)
    real(8), allocatable :: cosLatV(:), cosLatV3(:)
    real(8), allocatable :: v(:,:)

contains

    subroutine DeformationTestbed_Init(caseId_in, initId_in)
        integer, intent(in) :: caseId_in
        integer, intent(in) :: initId_in

        type(Mesh), pointer :: m1, m2
        procedure(), pointer :: p
        integer i, j

        call MsgManager_RecordSpeaker("DeformationTestbed_Init")

        caseId = caseId_in
        initId = initId_in

        T = 5.0d0
        r = Re/2.0d0

        select case (caseId)
        case (1)
            k = 2.4d0
            k05 = k*0.5d0
            lon1 = PI
            lat1 = PI/3.0d0
            lon2 = PI
            lat2 = -PI/3.0d0
        case (2)
            k = 2.0d0
            lon1 = 5.0d0/6.0d0*PI
            lat1 = 0.0d0
            lon2 = 7.0d0/6.0d0*PI
            lat2 = 0.0d0
        case (3)
            k = 1.0d0
            k05 = k*0.5d0
            lon1 = 3.0d0/4.0d0*PI
            lat1 = 0.0d0
            lon2 = 5.0d0/4.0d0*PI
            lat2 = 0.0d0
        case (4)
            k = PI2*Re/T
            lon1 = 5.0d0/6.0d0*PI
            lat1 = 0.0d0
            lon2 = 7.0d0/6.0d0*PI
            lat2 = 0.0d0
        end select

        select case (initId)
        case (1)
            b = 0.1d0
            c = 0.9d0
        case (2)
            b = 5.0d0
        case (3)
            b = 0.1d0
            c = 1.0d0
        end select

        m1 => MeshSelector("Zonal half mesh")
        m2 => MeshSelector("Meridianal half mesh")

        numLonU = m1%numLon
        numLatU = m1%numLat
        numLonV = m2%numLon
        numLatV = m2%numLat

        allocate(lonU(numLonU))
        allocate(sin2LonU05(numLonU))
        allocate(sinLonU2(numLonU))
        allocate(latU(numLatU))
        allocate(sinLatU2(numLatU))
        allocate(cosLatU(numLatU))
        allocate(cosLatU2(numLatU))
        allocate(u(numLonU,numLatU))
        allocate(lonV(numLonV))
        allocate(sinLonV(numLonV))
        allocate(sin2LonV(numLonV))
        allocate(latV(numLatV))
        allocate(cosLatV(numLatV))
        allocate(cosLatV3(numLatV))
        allocate(v(numLonV,numLatV))

        do i = 1, numLonU
            lonU(i) = m1%lon(i)
            sin2LonU05(i) = sin(lonU(i)*0.5d0)**2.0d0
            sinLonU2(i) = sin(lonU(i))**2.0d0
        end do
        do j = 1, numLatU
            latU(j) = m1%lat(j)
            sinLatU2(j) = sin(latU(j)*2.0d0)
            cosLatU(j) = cos(latU(j))
            cosLatU2(j) = cos(latU(j))**2.0d0
        end do

        do i = 1, numLonV
            lonV(i) = m2%lon(i)
            sinLonV(i) = sin(lonV(i))
            sin2LonV(i) = sin(lonV(i)*2.0d0)
        end do
        do j = 1, numLatV
            latV(j) = m2%lat(j)
            cosLatV(j) = cos(latV(j))
            cosLatV3(j) = cos(latV(j))**3.0d0
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
        deallocate(sinLonU2)
        deallocate(latU)
        deallocate(sinLatU2)
        deallocate(cosLatU)
        deallocate(cosLatU2)
        deallocate(u)
        deallocate(lonV)
        deallocate(sinLonV)
        deallocate(sin2LonV)
        deallocate(latV)
        deallocate(cosLatV)
        deallocate(cosLatV3)
        deallocate(v)

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine DeformationTestbed_Final
    
    subroutine DeformationTestbed_Advance(elapseTime)
        real(8), intent(in) :: elapseTime

        real(8) temp, lon
        integer i, j

        temp = cos(PI*elapseTime/T)
        
        select case (caseId)
        case (1)
            do j = 1, numLatU
                do i = 1, numLonU
                    u(i,j) = k*sin2LonU05(i)*sinLatU2(j)*temp
                end do
            end do
            do j = 1, numLatV
                do i = 1, numLonV
                    v(i,j) = k05*sinLonV(i)*cosLatV(j)*temp
                end do
            end do
        case (2)
            do j = 1, numLatU
                do i = 1, numLonU
                    u(i,j) = k*sinLonU2(i)*sinLatU2(j)*temp
                end do
            end do
            do j = 1, numLatV
                do i = 1, numLonV
                    v(i,j) = k*sin2LonV(i)*cosLatV(j)*temp
                end do
            end do
        case (3)
            do j = 1, numLatU
                do i = 1, numLonU
                    u(i,j) = -k*sin2LonU05(i)*sinLatU2(j)*cosLatU2(j)*temp
                end do
            end do
            do j = 1, numLatV
                do i = 1, numLonV
                    v(i,j) = k05*sinLonV(i)*cosLatV3(j)*temp
                end do
            end do
        case (4)
            do j = 1, numLatU
                do i = 1, numLonU
                    lon = lonU(i)-PI2*elapseTime/T
                    u(i,j) = k*sin(lon)**2.0d0*sinLatU2(j)*temp+PI2*Re*cosLatU(j)/T
                end do
            end do
            do j = 1, numLatV
                do i = 1, numLonV
                    lon = lonV(i)-PI2*elapseTime/T
                    v(i,j) = k*sin(lon*2.0d0)*cosLatV(j)*temp
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
        case (1)
            do i = 1, n
                call MeshManager_CalcDistance([lon(i),lat(i)], [lon1,lat1], r1)
                call MeshManager_CalcDistance([lon(i),lat(i)], [lon2,lat2], r2)
                if (r1 < r) then
                    q(i) = b+c*0.5d0*(1.0d0+cos(PI*r1/r))
                else if (r2 < r) then
                    q(i) = b+c*0.5d0*(1.0d0+cos(PI*r2/r))
                else
                    q(i) = b
                end if
            end do
        case (2)
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
                q1 = exp(-b*((x-x1)**2.0d0+(y-y1)**2.0d0+(z-z1)**2.0d0))
                q2 = exp(-b*((x-x2)**2.0d0+(y-y2)**2.0d0+(z-z2)**2.0d0))
                q(i) = q1+q2
            end do
        case (3)
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
                if (r2 <= r .and. dlon2 < r/6.0d0 .and. dlat2 < -5.0d0/12.0d0*r) then
                    q(i) = c
                end if
            end do
        end select

    end subroutine DeformationTestbed_CalcInitCond

end module DeformationTestbed

