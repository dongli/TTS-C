program TestFlowManager

    use RunManager
    use MsgManager
    use MeshManager
    use FlowManager
#if (defined SOLID_ROTATION)
    use SolidRotationTestbed, uTest => u, vTest => v
#elif (defined MOVING_VORTICES)
    use MovingVorticesTestbed, uTest => u, vTest => v
#elif (defined DEFORMATION)
    use DeformationTestbed, uTest => u, vTest => v
#endif
    use NFWrap

    implicit none

    integer numLon, numLat

    integer numTimeStep
    real(8) dt

    integer numTracerLon, numTracerLat

    character(256) namelistFile

    real(8) dlon, dlat
    real(8), allocatable :: lon(:), lat(:)
    real(8), allocatable :: u(:,:), v(:,:), div(:,:)

    type(FileCard) fcard
    type(Location) loc
    real(8) x(2), velocity(2)
    integer i, j, k

    namelist /MeshParams/ numLon, numLat
    namelist /TimeParams/ numTimeStep, dt
    namelist /TracerParams/ numTracerLon, numTracerLat

    call MsgManager_RecordSpeaker("TestFlowManager")

    call get_command_argument(1, namelistFile)

    open(10, file=namelistFile)
    read(10, nml=MeshParams)
    read(10, nml=TimeParams)
    read(10, nml=TracerParams)
    close(10)

    call MeshManager_Init(numLon, numLat)
    
    call FlowManager_Init("flow.nc")

#if (defined SOLID_ROTATION)
    call SolidRotationTestbed_Init
#elif (defined MOVING_VORTICES)
    call MovingVorticesTestbed_Init
#elif (defined DEFORMATION)
    call DeformationTestbed_Init(1, 1)
#endif

    call FlowManager_Link(uTest, vTest)
    call FlowManager_Output(0, 0.0d0)

    allocate(lon(numTracerLon))
    dlon = PI2/numTracerLon
    do i = 1, numTracerLon
        lon(i) = dlon/2.0d0+(i-1)*dlon
    end do

    allocate(lat(numTracerLat))
    dlat = PI/numTracerLat
    do j = 1, numTracerLat
        lat(j) = PI05-dlat/2.0d0-(j-1)*dlat
    end do

    call NFWrap_CreateSpherical2D("flow_interp.nc", numTracerLon, numTracerLat, fcard)
    call NFWrap_Output1DVar(fcard, "lon", lon*Rad2Deg)
    call NFWrap_Output1DVar(fcard, "lat", lat*Rad2Deg)
    call NFWrap_New2DVar(fcard, "u", "lon", "lat", "Zonal wind", "m s-1", timeVariant=.true.)
    call NFWrap_New2DVar(fcard, "v", "lon", "lat", "Meridianal wind", "m s-1", timeVariant=.true.)
    call NFWrap_New2DVar(fcard, "div", "lon", "lat", "Divergence", "s-1", timeVariant=.true.)

    allocate(u(numTracerLon,numTracerLat))
    allocate(v(numTracerLon,numTracerLat))
    allocate(div(numTracerLon,numTracerLat))
    do j = 1, numTracerLat
        do i = 1, numTracerLon
            x = [lon(i),lat(j)]
            call MeshManager_LocationCheck(x, loc)
            call FlowManager_GetVelocity(x, loc, velocity, "new")
            u(i,j) = velocity(1)
            v(i,j) = velocity(2)
            call FlowManager_GetDivergence(x, loc, div(i,j), "new")
        end do
    end do
    call NFWrap_Advance(fcard, 0, 0.0d0)
    call NFWrap_Output2DVar(fcard, "u", u)
    call NFWrap_Output2DVar(fcard, "v", v)
    call NFWrap_Output2DVar(fcard, "div", div)

#if (defined MOVING_VORTICES || defined DEFORMATION)
    do k = 1, numTimeStep
        do j = 1, numTracerLat
            do i = 1, numTracerLon
                x = [lon(i),lat(j)]
                call MeshManager_LocationCheck(x, loc)
                call FlowManager_GetVelocity(x, loc, velocity, "new")
                u(i,j) = velocity(1)
                v(i,j) = velocity(2)
                call FlowManager_GetDivergence(x, loc, div(i,j), "new")
            end do
        end do
        call NFWrap_Advance(fcard, k, k*dt)
        call NFWrap_Output2DVar(fcard, "u", u)
        call NFWrap_Output2DVar(fcard, "v", v)
        call NFWrap_Output2DVar(fcard, "div", div)
        call FlowManager_BeforeAdvance
#if (defined MOVING_VORTICES)
        call MovingVorticesTestbed_Advance(k*dt)
#elif (defined DEFORMATION)
        call DeformationTestbed_Advance(k*dt)
#endif
        call FlowManager_AfterAdvance
        call FlowManager_Output(k, k*dt)
        write(*, "('---> time step ', I)") k
    end do
#endif

    call NFWrap_Close(fcard)

    call MsgManager_Speak(Notice, "Finished.")
    call MsgManager_DeleteSpeaker

    call RunManager_EndRun

end program TestFlowManager
