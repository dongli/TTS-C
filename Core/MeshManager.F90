! *************************************************************************** !
! MeshManager module (regular latitude-longitude version)                     !
!                                                                             !
! Description:                                                                !
!                                                                             !
! Author:                                                                     !
!                                                                             !
! *************************************************************************** !

module MeshManager

    use CommonTypes
    use SphereService
    use MsgManager
    use RunManager

    implicit none

    integer, parameter :: numDim = 2
    real(8), parameter :: PoleR = 1.5d0/Rad2Deg
    character(10) :: coordinateShortName(2) = ["lon",       "lat"     ]
    character(50) :: coordinateLongName(2)  = ["longitude", "latitude "]
    character(10) :: coordinateUnit(2)      = ["degree_E",  "degree_N"]

    type DimensionInformation
        integer numLon, numLat
        real(8) dlon, dlat
        real(8) dlon05, dlat05
    end type DimensionInformation
    
    type(DimensionInformation) dimInfo

    type Point
        real(8) lon, lat
    end type Point

    type Mesh
        integer id
        character(30) name
        integer numLon, numLat
        real(8), allocatable :: lon(:), lat(:)
        real(8), allocatable :: cosLat(:)
        real(8), allocatable :: area(:,:)
    contains
        procedure :: dump => Mesh_dump
    end type Mesh

    integer, parameter :: numMeshType = 4

    type(Mesh), target :: meshes(numMeshType)

    type Quantity
        character(30) name
        character(30) unit
        type(Mesh), pointer :: m
        type(TwoTimeLevel), allocatable :: val(:,:)
    contains
        procedure :: dump => Quantity_dump
    end type Quantity

    integer, parameter :: north = 1
    integer, parameter :: south = 2

    type PolarCapQuantity
        integer numLon
        real(8) lat(north:south)
        real(8), allocatable :: lon(:)
        type(TwoTimeLevel), allocatable :: val(:,:)
    contains
        procedure :: stitch => PolarCapQuantity_stitch
    end type PolarCapQuantity

    type Location
        integer i(numMeshType)
        integer j(numMeshType)
        logical inPolarCap
        integer polej
        logical onPole
    contains
        procedure :: dump => Location_dump
    end type Location

    type(OperationList), target :: LocationCheck
    
contains

    subroutine MeshManager_Init(numLon, numLat)
        integer, intent(in) :: numLon, numLat
        
        integer i, j

        call MsgManager_RecordSpeaker("MeshManager_Init")

        call MsgManager_AddConfig("MeshManager", "PoleR", real2str(PoleR*Rad2Deg, 7, 5))

        dimInfo%numLon = numLon
        dimInfo%numLat = numLat
        dimInfo%dlon = PI2/numLon
        dimInfo%dlat = PI/numLat
        dimInfo%dlon05 = dimInfo%dlon*0.5d0
        dimInfo%dlat05 = dimInfo%dlat*0.5d0

        ! Arakawa C staggered meshes (3 types)
        do i = 1, numMeshType
            meshes(i)%id = i
        end do

        ! .................................................................... A
        meshes(1)%name = "Full mesh"
        meshes(1)%numLon = numLon
        meshes(1)%numLat = numLat-1
        !                        ^^
        allocate(meshes(1)%lon(0:meshes(1)%numLon+1))
        !                      ^^                ^^
        allocate(meshes(1)%lat(meshes(1)%numLat))
        allocate(meshes(1)%cosLat(meshes(1)%numLat))
        do i = 0, meshes(1)%numLon+1
            meshes(1)%lon(i) = i*dimInfo%dlon
        end do
        do j = 1, meshes(1)%numLat
            meshes(1)%lat(j) = PI05-j*dimInfo%dlat
            meshes(1)%cosLat(j) = cos(meshes(1)%lat(j))
        end do
        allocate(meshes(1)%area(meshes(1)%numLon,0:meshes(1)%numLat+1))
        do j = 1, meshes(1)%numLat
            do i = 1, meshes(1)%numLon
                meshes(1)%area(i,j) = Re2*dimInfo%dlon* &
                    (sin(meshes(1)%lat(j)+dimInfo%dlat05)- &
                     sin(meshes(1)%lat(j)-dimInfo%dlat05))
            end do
        end do
        j = 0
        do i = 1, meshes(1)%numLon
            meshes(1)%area(i,j) = Re2*dimInfo%dlon*(1.0d0-sin(PI05-dimInfo%dlat05))
        end do
        j = meshes(1)%numLat+1
        do i = 1, meshes(1)%numLon
            meshes(1)%area(i,j) = Re2*dimInfo%dlon*(1.0d0-sin(PI05-dimInfo%dlat05))
        end do

        ! .................................................................... B
        meshes(2)%name = "Zonal half mesh"
        meshes(2)%numLon = numLon
        meshes(2)%numLat = numLat-1
        !                        ^^
        allocate(meshes(2)%lon(0:meshes(2)%numLon+1))
        !                      ^^                ^^
        allocate(meshes(2)%lat(meshes(2)%numLat))
        allocate(meshes(2)%cosLat(meshes(2)%numLat))
        do i = 0, meshes(2)%numLon+1
            meshes(2)%lon(i) = i*dimInfo%dlon-dimInfo%dlon05
        end do
        do j = 1, meshes(2)%numLat
            meshes(2)%lat(j) = PI05-j*dimInfo%dlat
            meshes(2)%cosLat(j) = cos(meshes(2)%lat(j))
        end do

        ! .................................................................... C
        meshes(3)%name = "Meridianal half mesh"
        meshes(3)%numLon = numLon
        meshes(3)%numLat = numLat
        allocate(meshes(3)%lon(0:meshes(3)%numLon+1))
        !                      ^^                ^^
        allocate(meshes(3)%lat(meshes(3)%numLat))
        allocate(meshes(3)%cosLat(meshes(3)%numLat))
        do i = 0, meshes(3)%numLon+1
            meshes(3)%lon(i) = i*dimInfo%dlon
        end do
        do j = 1, meshes(3)%numLat
            meshes(3)%lat(j) = PI05-j*dimInfo%dlat+dimInfo%dlat05
            meshes(3)%cosLat(j) = cos(meshes(3)%lat(j))
        end do

        ! .................................................................... D
        meshes(4)%name = "Double half mesh"
        meshes(4)%numLon = numLon
        meshes(4)%numLat = numLat
        allocate(meshes(4)%lon(0:meshes(4)%numLon+1))
        allocate(meshes(4)%lat(meshes(4)%numLat))
        do i = 0, meshes(4)%numLon+1
            meshes(4)%lon(i) = i*dimInfo%dlon-dimInfo%dlon05
        end do
        do j = 1, meshes(4)%numLat
            meshes(4)%lat(j) = PI05-j*dimInfo%dlat+dimInfo%dlat05
        end do

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine MeshManager_Init

    subroutine MeshManager_RegisterOperation(hostName, modName, subName, handle)
        character(*), intent(in) :: hostName, modName, subName
        procedure(), intent(in), pointer :: handle

        call MsgManager_RecordSpeaker("MeshManager_RegisterOperation")
    
        select case (hostName)
        case ("LocationCheck")
            call LocationCheck%register(modName, subName, handle)
        case default
            call MsgManager_Speak(Error, &
                "No host subroutine """// &
                trim(hostName)//""" for registration.")
            call MsgManager_DeleteSpeaker
            call RunManager_EndRun
        end select
#if (defined VERBOSE)
        call MsgManager_Speak(Notice, &
            "Add operation """// &
            trim(modName)//"::"//trim(subName)// &
            """ to """//trim(hostName)//""".")
#endif
        call MsgManager_DeleteSpeaker
    
    end subroutine MeshManager_RegisterOperation
    
    function MeshSelector(meshName) result(res)
        character(*), intent(in) :: meshName

        type(Mesh), pointer :: res

        integer i

        call MsgManager_RecordSpeaker("MeshSelector")

        do i = 1, numMeshType
            if (meshName == meshes(i)%name) then
                res => meshes(i)
                call MsgManager_DeleteSpeaker
                return
            end if
        end do

        call MsgManager_Speak(Error, &
            "MeshSelector encounter unknown mesh "// &
            """"//trim(meshName)//""".")
        call RunManager_EndRun

    end function MeshSelector
    
    subroutine MeshManager_LinkMesh(meshName, q)
        character(*), intent(in) :: meshName
        type(Quantity), intent(inout) :: q
    
        q%m => MeshSelector(meshName)

    end subroutine MeshManager_LinkMesh
    
    subroutine MeshManager_LocationCheck(x, loc, tracerId, sampleId)
        real(8), intent(in) :: x(2)
        type(Location), intent(inout) :: loc
        integer, intent(in), optional :: tracerId, sampleId

        real(8) temp
        type(SubHandle), pointer :: sh
        integer i

        call MsgManager_RecordSpeaker("MeshManager_LocationCheck")

        temp = PI05-x(2)
        temp = merge(0.0d0, temp, abs(temp) < 1.0d-6)

        ! full mesh
        loc%i(1) = min(floor(x(1)/dimInfo%dlon-1)+1, meshes(1)%numLon)
        loc%j(1) = min(floor(temp/dimInfo%dlat-1)+1, meshes(1)%numLat)
        if (loc%j(1) == -1) then
            call MsgManager_Speak(Error, "Bad thing happened. loc%j(1) = -1!")
            call RunManager_EndRun
        end if

        ! zonal half mesh
        loc%i(2) = min(floor(x(1)/dimInfo%dlon-0.5d0)+1, meshes(2)%numLon)
        loc%j(2) = loc%j(1)

        ! meridianal half mesh
        loc%i(3) = loc%i(1)
        loc%j(3) = floor(temp/dimInfo%dlat-0.5d0)+1

        ! double half mesh
        loc%i(4) = loc%i(2)
        loc%j(4) = loc%j(3)

        if (loc%j(1) == 0) then
            loc%inPolarCap = .true.
            loc%polej = north
        else if (loc%j(1) == meshes(1)%numLat) then
            loc%inPolarCap = .true.
            loc%polej = south
        else
            loc%inPolarCap = .false.
            loc%polej = -1
        end if

        if (PI05-abs(x(2)) < PoleR) then
            loc%onPole = .true.
            if (x(2) > 0.0d0) then
                loc%polej = north
            else
                loc%polej = south
            end if
        else
            loc%onPole = .false.
        end if

        if (present(tracerId) .and. present(sampleId)) then
            sh => LocationCheck%head
            do i = 1, LocationCheck%num
                call sh%handle(tracerId, sampleId, loc)
                sh => sh%next
            end do
        end if

        call MsgManager_DeleteSpeaker

    end subroutine MeshManager_LocationCheck

    logical function MeshManager_IsNearPole(loc) result(res)
        type(Location), intent(in) :: loc

        res = .false.
        if (loc%onPole) then
            res = .true.
        end if

    end function MeshManager_IsNearPole
    
    subroutine MeshManager_CalcDistance(x1, x2, d)
        real(8), intent(in) :: x1(2), x2(2)
        real(8), intent(out) :: d

        real(8) temp1, temp2, temp3
        real(8) dlon

        dlon = x1(1)-x2(1)
        temp1 = sin(x1(2))*sin(x2(2))
        temp2 = cos(x1(2))*cos(x2(2))*cos(dlon)
        temp3 = temp1+temp2
        temp3 = min(temp3, 1.0d0)
        temp3 = max(temp3, -1.0d0)
        d = Re*acos(temp3)

    end subroutine MeshManager_CalcDistance

    subroutine MeshManager_Move(x0, x1, v, dt, loc, sampleId)
        real(8), intent(in) :: x0(2), v(2), dt
        real(8), intent(out) :: x1(2)
        type(Location), intent(in) :: loc
        integer, intent(in), optional :: sampleId

        real(8) dx(2)

        call MsgManager_RecordSpeaker("MeshManager_Move")

        dx(1) = v(1)*dt/Re/cos(x0(2))
        dx(2) = v(2)*dt/Re

        x1 = x0+dx

        ! polar boundary check
        if (x1(2) > PI05) then
            x1(1) = PI+x0(1)-dx(1)
            x1(2) = PI-x0(2)-dx(2)
        end if
        if (x1(2) < -PI05) then
            x1(1) = PI+x0(1)-dx(1)
            x1(2) = -PI-x0(2)-dx(2)
        end if
        ! zonal boundary check
        if (x1(1) < 0.0d0) then
#if (defined DEBUG)
            if (x1(1)+PI2 < 0.0d0) then
                call MsgManager_Speak(Warning, &
                    "Tracer moves more than one "// &
                    "zonal circle in one time step.")
            end if
#endif
            x1(1) = PI2+mod(x1(1), PI2)
        end if
        if (x1(1) > PI2) then
#if (defined DEBUG)
            if (x1(1)-PI2 > PI2) then
                call MsgManager_Speak(Warning, &
                    "Tracer moves more than one "// &
                    "zonal circle in one time step.")
            end if
#endif
            x1(1) = mod(x1(1), PI2)
        end if

        call MsgManager_DeleteSpeaker

    end subroutine MeshManager_Move

    subroutine MeshManager_MoveNearPole(x0, x1, vt, dt, loc)
        real(8), intent(in) :: x0(2), vt(2), dt
        real(8), intent(out) :: x1(2)
        type(Location), intent(in) :: loc

        real(8) xt0(2), xt1(2)

        ! transform into local transformed coordinate system
        call PolarStereoTransformCoordinate(x0, xt0, loc%polej)
        ! move in local transformed coordinate system
        xt1 = xt0+vt*dt
        ! transform into spherical coordinate system
        call PolarStereoTransformCoordinateBack(xt1, x1, loc%polej)
   
    end subroutine MeshManager_MoveNearPole
    
    subroutine PolarStereoTransformCoordinate(x0, x1, polej)
        real(8), intent(in) :: x0(2)
        real(8), intent(out) :: x1(2)
        integer, intent(in) :: polej

        if (polej == north) then
            x1(1) = Re*cos(x0(1))/tan(x0(2))
            x1(2) = Re*sin(x0(1))/tan(x0(2))
        else
            x1(1) = -Re*cos(x0(1))/tan(x0(2))
            x1(2) = -Re*sin(x0(1))/tan(x0(2))
        end if

    end subroutine PolarStereoTransformCoordinate

    subroutine PolarStereoTransformCoordinateBack(x0, x1, polej)
        real(8), intent(in) :: x0(2)
        real(8), intent(out) :: x1(2)
        integer, intent(in) :: polej

        x1(1) = atan2(x0(2), x0(1))
        if (x1(1) < 0.0d0) x1(1) = PI2+x1(1)

        if (polej == north) then
            x1(2) = atan(Re/sqrt(x0(1)**2.0d0+x0(2)**2.0d0))
        else
            x1(2) = -atan(Re/sqrt(x0(1)**2.0d0+x0(2)**2.0d0))
        end if

    end subroutine PolarStereoTransformCoordinateBack
    
    subroutine PolarCapQuantity_stitch(pcap, normal, offset)
        class(PolarCapQuantity), intent(out) :: pcap
        type(Quantity), intent(in) :: normal
        character(*), intent(in), optional :: offset

        call MsgManager_RecordSpeaker("PolarCapQuantity_stitch")

        pcap%numLon = normal%m%numLon
        allocate(pcap%lon(normal%m%numLon))
        if (.not. present(offset)) then
            pcap%lon = normal%m%lon(1:normal%m%numLon)
            pcap%lat(north) = normal%m%lat(1)
            pcap%lat(south) = normal%m%lat(normal%m%numLat)
        else
            select case (offset)
            case ("half offset")
                pcap%lon = normal%m%lon(1:normal%m%numLon) &
                           -dimInfo%dlon05
                pcap%lat(north) = normal%m%lat(1)-dimInfo%dlat05
                pcap%lat(south) = normal%m%lat(normal%m%numLat)+dimInfo%dlat05
            case default
                call MsgManager_Speak(Error, &
                    "Unknown offset type "// &
                    """"//trim(offset)//""".")
                call RunManager_EndRun
            end select
        end if
        
        allocate(pcap%val(pcap%numLon,north:south))
    
        call MsgManager_DeleteSpeaker

    end subroutine PolarCapQuantity_stitch
    
    subroutine Mesh_dump(m)
        class(Mesh), intent(in) :: m

        write(*, *)
        write(*, "('Mesh information of ', A)") trim(m%name)
        write(*, "('  numLon: ', I5)") m%numLon
        write(*, "('  numLat: ', I5)") m%numLat
        write(*, "('  lon: ')")
        write(*, *) m%lon*Rad2Deg
        write(*, "('  lat: ')")
        write(*, *) m%lat*Rad2Deg
        write(*, *)

    end subroutine Mesh_dump
    
    subroutine Quantity_dump(q)
        class(Quantity), intent(in) :: q

        write(*, *)
        write(*, "('Quantity information of ""', A, '""')") trim(q%name)
        write(*, "('  Mesh type: ', A)") trim(q%m%name)
        write(*, *)
    
    end subroutine Quantity_dump
    
    subroutine Location_dump(loc)
        class(Location), intent(in) :: loc

        integer k

        write(*, "('Location dump result:')")
        do k = 1, numMeshType
            write(*, "('  Index of mesh ', A, ': ')", advance="no") trim(meshes(k)%name)
            write(*, "('(', A, ',', A, ')')") trim(int2str(loc%i(k))), trim(int2str(loc%j(k)))
        end do
        write(*, "('  In the polar cap: ', L6)") loc%inPolarCap
    
    end subroutine Location_dump

end module MeshManager

