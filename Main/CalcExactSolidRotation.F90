program CalcExactSolidRotation

    use MsgManager
    use RunManager
    use SphereService
    use MeshManager
    use SolidRotationTestbed
    use NFWrap

    implicit none

    integer numLon, numLat

    integer numTimeStep
    real(8) dt

    character(50) namelistFile
    character(50) outputFile

    real(8), allocatable :: qExact(:,:)

    logical flag
    integer k

    type(Mesh), pointer :: m
    type(FileCard) fcard
    real(8), allocatable :: lon(:), lat(:)

    namelist /MeshParams/ numLon, numLat
    namelist /TimeParams/ numTimeStep, dt

    call MsgManager_RecordSpeaker("CalcExactSolidRotation")

    call get_command_argument(1, namelistFile)
    call get_command_argument(2, outputFile)

    if (namelistFile == "") then
        call MsgManager_Speak(Error, &
            "A namelist file should be given as command argument.")
        call RunManager_EndRun
    end if

    inquire(file=namelistFile, exist=flag)
    if (.not. flag) then
        call MsgManager_Speak(Error, &
            "Namelist file """//trim(namelistFile)//""" does not exist.")
        call RunManager_EndRun
    end if

    open(11, file=namelistFile)
    read(11, nml=MeshParams)
    read(11, nml=TimeParams)
    close(11)

    call MeshManager_Init(numLon, numLat)
    call SolidRotationTestbed_Init

    m => MeshSelector("Full mesh")
    allocate(lon(m%numLon))
    lon = m%lon(1:m%numLon)
    allocate(lat(0:m%numLat+1))
    lat(0) = PI05
    lat(1:m%numLat) = m%lat
    lat(m%numLat+1) = -PI05
    call NFWrap_CreateSpherical2D(outputFile, m%numLon, m%numLat+2, fcard)
    call NFWrap_Output1DVar(fcard, "lon", lon*Rad2Deg)
    call NFWrap_Output1DVar(fcard, "lat", lat*Rad2Deg)
    call NFWrap_New2DVar(fcard, "q", "lon", "lat", "moisture", "g m-3", timeVariant=.true.)
    allocate(qExact(m%numLon,0:m%numLat+1))

    do k = 0, numTimeStep
        call NFWrap_Advance(fcard, k, k*dt)
        call SolidRotationTestbed_CalcExactSolution( &
            m%numLon, m%numLat+2, lon, lat, k*dt, qExact)
        call NFWrap_Output2DVar(fcard, "q", qExact)
        print *, "---> time step", k
    end do

    call NFWrap_Close(fcard)

    call MsgManager_Speak(Notice, "Finished.")
    call MsgManager_DeleteSpeaker

    call RunManager_EndRun

end program CalcExactSolidRotation
