program DriveMovingVortices

    use MsgManager
    use RunManager
    use NFWrap
    use SphereService
    use MeshManager
    use FlowManager
    use TracerManager
    use TTS
    use MovingVorticesTestbed, uTest => u, vTest => v
#if (defined LAGRANGE_TO_EULER)
    use LagrangeToEuler
#endif

    ! ************************************************************************ !
    !                          RUNNING CONFIGURATION
    ! ************************************************************************ !
    integer numLon, numLat

    integer numTimeStep
    real(8) dt

    integer numTracerLon, numTracerLat
    integer numTracer

    integer tracerId
    real(8) dx(2)
    real(8), allocatable :: x(:,:), q(:)

    character(256) namelistFile
    character(256) searchlistFile
    character(256) flowFile
    character(256) gridFile
    character(256) tracerFile
    type(FileCard) fcard

    logical flag
    integer i, j, k

    namelist /MeshParams/ numLon, numLat
    namelist /TimeParams/ numTimeStep, dt
    namelist /TracerParams/ numTracerLon, numTracerLat
    namelist /FileIOParams/ searchlistFile, flowFile, gridFile, tracerFile

    call MsgManager_RecordSpeaker("DriveMovingVortices")

    ! ************************************************************************ !
    !                           READ IN NAMELIST
    ! ************************************************************************ !

    call get_command_argument(1, namelistFile)

    if (namelistFile == "") then
        call MsgManager_Speak(Error, &
            "A namelist file is needed.")
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
    read(11, nml=TracerParams)
    read(11, nml=FileIOParams)
    close(11)

    numTracer = numTracerLon*numTracerLat

    ! ************************************************************************ !
    !                           MODULE INITIALIZATION
    ! ************************************************************************ !
    call MeshManager_Init(numLon, numLat)
    call FlowManager_Init(flowFile)
    call TracerManager_Init(1, tracerFile)
    call MovingVorticesTestbed_Init
#if (defined LAGRANGE_TO_EULER)
    call LagrangeToEuler_Init(searchlistFile, gridFile)
#endif
    call MsgManager_ShowConfig

    ! ************************************************************************ !
    !                      ADVECTION VELOCITY INITIALIZATION
    ! ************************************************************************ !
    call FlowManager_Link(uTest, vTest)
    call FlowManager_Output(0, 0.0d0)

    ! ************************************************************************ !
    !                           TRACER INITIALIZATION
    ! ************************************************************************ !
    ! set tracer coordinates
    if (.false.) then
        allocate(x(numDim,numTracer))
        dx(1) = PI2/numTracerLon
        dx(2) = PI/numTracerLat
        do j = 1, numTracerLat
            do i = 1, numTracerLon
                k = i+(j-1)*numTracerLon
                x(1,k) = dx(1)/2.0d0+(i-1)*dx(1)
                x(2,k) = PI05-dx(2)/2.0d0-(j-1)*dx(2)
            end do
        end do
    else
        print *, "#############################################################################"
        print *, "Using external data C02562.global.nc for tracer sample's initial coordinates."
        print *, "#############################################################################"
        call NFWrap_OpenForRead("C02562.global.nc", fcard)
        call NFWrap_GetDimSize(fcard, "grid_size", numTracer)
        allocate(x(numDim,numTracer))
        call NFWrap_Input1DVar(fcard, "grid_center_lon", x(1,:))
        call NFWrap_Input1DVar(fcard, "grid_center_lat", x(2,:))
        call NFWrap_Close(fcard)
    end if
    ! set tracer quantities
    allocate(q(numTracer))
    call TracerManager_RegisterTracer("q", numTracer, Conservative, &
        "moisture density", "g m-3", tracerId)
    call MovingVorticesTestbed_CalcExactSolution( &
        numTracer, x(1,:), x(2,:), 0.0d0, q)
    call TracerManager_InitialCondition(tracerId, numTracer, x, q)
    call TracerManager_Output(0, 0.0d0)
#if (defined LAGRANGE_TO_EULER)
    call LagrangeToEuler_Warmup
    call LagrangeToEuler_Output(0, 0.0d0)
#endif

    ! *********************************************************************** !
    !                              TIME INTEGRATION
    ! *********************************************************************** !
    write(*, "('*** Starting time integration ***')")
    do k = 1, numTimeStep
        call TTS_AdvectTracer(tracerId, dt)
#if (defined LAGRANGE_TO_EULER)
        call LagrangeToEuler_Run
        call LagrangeToEuler_Output(k, k*dt)
#endif
        call TracerManager_Output(k, k*dt)
        call FlowManager_BeforeAdvance
        call MovingVorticesTestbed_Advance(k*dt)
        call FlowManager_AfterAdvance
        call FlowManager_Output(k, k*dt)
        print *, "---> time step ", k
    end do

    call RunManager_EndRun

end program DriveMovingVortices
