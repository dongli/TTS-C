! *************************************************************************** !
! LagrangeToEuler module                                                      !
!                                                                             !
! Description:                                                                !
!                                                                             !
! Author:                                                                     !
!                                                                             !
! *************************************************************************** !

module LagrangeToEuler

    use MsgManager
    use RunManager
    use MeshManager
    use FlowManager
    use TracerManager
    use SphereService
    use FloatingPoint
    use NFWrap

    implicit none

    private

    public LagrangeToEuler_Init
    public LagrangeToEuler_Final
    public LagrangeToEuler_Warmup
    public LagrangeToEuler_Run
    public LagrangeToEuler_Output

    public CellMesh

    integer, parameter :: maxNumSample = 2000
    integer, parameter :: numQuadNormal = 4
    real(8), parameter :: dlonNormal = PI2/numQuadNormal
    integer, parameter :: numQuadPole = 4
    real(8), parameter :: dlonPole = PI2/numQuadPole
    integer maxSearchLevel

    type SearchList
        integer :: numPoint = 0
        integer, allocatable :: i(:), j(:)
    end type SearchList

    type TracerCounter
        integer :: numSample = 0
        integer sampleId(maxNumSample)
    end type TracerCounter

    type Cell
        type(Mesh), pointer :: m
        type(SearchList), allocatable :: sl(:)
        type(TracerCounter), allocatable :: tc(:)
        real(8), allocatable :: q(:)
    end type Cell

    type CellMesh
        type(Mesh), pointer :: m ! center grid
        type(Cell), allocatable :: c(:,:)
        real(8), allocatable :: totalMass(:)
    end type CellMesh

    type(CellMesh) cm

    type(FileCard) fcard

contains

    subroutine LagrangeToEuler_Init(srchlistFile, outputFile)
        character(*), intent(in) :: srchlistFile, outputFile

        logical flag
        character option
        real(8), allocatable :: lon(:), lat(:)
        procedure(), pointer :: p
        integer i, j, k, n

        call MsgManager_RecordSpeaker("LagrangeToEuler_Init")

        call MsgManager_AddConfig("LagrangeToEuler", "maxNumSample", int2str(maxNumSample))
        call MsgManager_AddConfig("LagrangeToEuler", "numQuadNormal", int2str(numQuadNormal))
        call MsgManager_AddConfig("LagrangeToEuler", "numQuadPole", int2str(numQuadPole))
    
        ! ==================================================================== !
        ! initialize the cell
        cm%m => MeshSelector("Full mesh")
        allocate(cm%c(cm%m%numLon,0:cm%m%numLat+1))
        n = TracerManager_GetMaxTracerClassNumber()
        ! j = 0 and j = numLat+1 is the extra polar cap, 
        ! only used for recording tracer samples, and storing pole quantity
        allocate(cm%totalMass(n))
        do j = 0, cm%m%numLat+1
        do i = 1, cm%m%numLon
            allocate(cm%c(i,j)%tc(n))
            allocate(cm%c(i,j)%q(n))
        end do
        end do
        ! read in the search list for each grid
        inquire(file=srchlistFile, exist=flag)
        if (.not. flag) then
            call MsgManager_Speak(Error, &
                "File """//trim(srchlistFile)//""" does not exist.")
            call RunManager_EndRun
        end if
        open(10, file=srchlistFile, form="unformatted", status="old", action="read")
        read(10) maxSearchLevel
        do j = 1, cm%m%numLat
        do i = 1, cm%m%numLon
            allocate(cm%c(i,j)%sl(maxSearchLevel))
        end do
        end do
        do j = 1, cm%m%numLat
        do i = 1, cm%m%numLon
            do k = 1, maxSearchLevel
                read(10) cm%c(i,j)%sl(k)%numPoint
                allocate(cm%c(i,j)%sl(k)%i(cm%c(i,j)%sl(k)%numPoint))
                allocate(cm%c(i,j)%sl(k)%j(cm%c(i,j)%sl(k)%numPoint))
                read(10) cm%c(i,j)%sl(k)%i
                read(10) cm%c(i,j)%sl(k)%j
            end do
        end do
        end do
        close(10)

        ! ==================================================================== !
        ! initialize the output NETCDF file
        inquire(file=outputFile, exist=flag)
        if (flag) then
            call MsgManager_Speak(Warning, &
                "File """//trim(outputFile)//""" exists already, overwrite it? (y/n)")
            read(*, *) option
            if (option == "n") then
                call RunManager_EndRun
            end if
        end if
        call NFWrap_CreateSpherical2D(outputFile, cm%m%numLon, cm%m%numLat+2, fcard)
        allocate(lon(cm%m%numLon))
        lon = cm%m%lon(1:cm%m%numLon)*Rad2Deg
        allocate(lat(0:cm%m%numLat+1))
        lat(0) = 90.0d0
        lat(1:cm%m%numLat) = cm%m%lat*Rad2Deg
        lat(cm%m%numLat+1) = -90.0d0
        call NFWrap_Output1DVar(fcard, "lon", lon)
        call NFWrap_Output1DVar(fcard, "lat", lat)
        call NFWrap_New2DVar(fcard, "area", "lon", "lat", &
            "grid cell area", "m2", .false.)
        call NFWrap_Output2DVar(fcard, "area", cm%m%area)

        ! ==================================================================== !
        ! register operation

#if (defined FC_GFORTRAN)
        p => SampleCounter
        call MeshManager_RegisterOperation("LocationCheck", &
            "LagrangeToEuler", "SampleCounter", p)
        p => RegisterTracer
        call TracerManager_RegisterOperation("RegisterTracer", &
            "LagrangeToEuler", "RegisterTracer", p)
        p => LagrangeToEuler_Final
        call RunManager_RegisterOperation("EndRun", &
            "LagrangeToEuler", "Final", p)
#elif (defined FC_IFORT)
        call MeshManager_RegisterOperation("LocationCheck", &
            "LagrangeToEuler", "SampleCounter", SampleCounter)
        call TracerManager_RegisterOperation("RegisterTracer", &
            "LagrangeToEuler", "RegisterTracer", RegisterTracer)
        call RunManager_RegisterOperation("EndRun", &
            "LagrangeToEuler", "Final", LagrangeToEuler_Final)
#endif

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker
    
    end subroutine LagrangeToEuler_Init
    
    subroutine LagrangeToEuler_Final
    
        call MsgManager_RecordSpeaker("LagrangeToEuler_Final")

        call NFWrap_Close(fcard)

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker
    
    end subroutine LagrangeToEuler_Final
    
    subroutine RegisterTracer(varName, longName, unitName)
        character(*), intent(in) :: varName, longName, unitName

        call MsgManager_RecordSpeaker("LagrangeToEuler::RegisterTracer")

        call NFWrap_New2DVar(fcard, varName, "lon", "lat", &
            longName, unitName, timeVariant=.true.)

        call MsgManager_DeleteSpeaker

    end subroutine RegisterTracer

    subroutine SampleCounter(tracerId, sampleId, loc)
        integer, intent(in) :: tracerId, sampleId
        type(Location), intent(in) :: loc
    
        integer i, j

        call MsgManager_RecordSpeaker("SampleCounter")
    
        i = merge(loc%i(4), cm%m%numLon, loc%i(4) /= 0)
        j = loc%j(4)
        cm%c(i,j)%tc(tracerId)%numSample = cm%c(i,j)%tc(tracerId)%numSample+1
        if (cm%c(i,j)%tc(tracerId)%numSample > maxNumSample) then
            call MsgManager_Speak(Error, &
                "Exceed maximum sample number.")
            call RunManager_EndRun
        end if
        cm%c(i,j)%tc(tracerId)%sampleId(cm%c(i,j)%tc%numSample) = sampleId

        call MsgManager_DeleteSpeaker
    
    end subroutine SampleCounter

    subroutine LagrangeToEuler_Warmup
        integer i, j, k

        call MsgManager_RecordSpeaker("LagrangeToEuler_Warmup")
    
        call LagrangeToEuler_IDWInterp

        do k = 1, TracerManager_GetTracerClassNumber()
            if (TracerManager_IsConservative(k)) then
                cm%totalMass(k) = 0.0d0
                do j = 0, cm%m%numLat+1
                do i = 1, cm%m%numLon
                    cm%totalMass(k) = cm%totalMass(k)+cm%c(i,j)%q(k)*cm%m%area(i,j)
                end do
                end do
            end if
        end do

        call MsgManager_DeleteSpeaker
    
    end subroutine LagrangeToEuler_Warmup

    subroutine LagrangeToEuler_Run
        real(8) totalMass, weight
        integer i, j, k

        call MsgManager_RecordSpeaker("LagrangeToEuler_Run")
    
        call LagrangeToEuler_IDWInterp

#if (!defined NO_RENORMALIZATION)
        ! renormalization for mass conservation
        do k = 1, TracerManager_GetTracerClassNumber()
            if (.not. TracerManager_IsConservative(k)) cycle
            totalMass = 0.0d0
            do j = 0, cm%m%numLat+1
            do i = 1, cm%m%numLon
                totalMass = totalMass+cm%c(i,j)%q(k)*cm%m%area(i,j)
            end do
            end do
            weight = cm%totalMass(k)/totalMass
            do j = 0, cm%m%numLat+1
            do i = 1, cm%m%numLon
                cm%c(i,j)%q(k) = cm%c(i,j)%q(k)*weight
            end do
            end do
        end do
#endif

        call MsgManager_DeleteSpeaker
    
    end subroutine LagrangeToEuler_Run
    
    subroutine LagrangeToEuler_IDWInterp
        real(8) x1(2), x2(2), distance, lonR
        real(8) weight, weightSum, q
        integer supSmpIdNormal(numQuadNormal)
        real(8) supSmpDisNormal(numQuadNormal)
        integer supSmpIdPole(numQuadPole)
        real(8) supSmpDisPole(numQuadPole)
        integer i, j, k, l, n
        integer lev, pnt, smp, quad, ii, jj

        call MsgManager_RecordSpeaker("LagrangeToEuler_IDWInterp")

        n = TracerManager_GetTracerClassNumber()

        ! normal area
        do j = 1, cm%m%numLat
        do i = 1, cm%m%numLon
            x1 = [cm%m%lon(i),cm%m%lat(j)]
            do l = 1, n
                supSmpIdNormal = -1
                do lev = 1, maxSearchLevel
                    do pnt = 1, cm%c(i,j)%sl(lev)%numPoint
                        ii = cm%c(i,j)%sl(lev)%i(pnt)
                        jj = cm%c(i,j)%sl(lev)%j(pnt)
                        if (jj < 0 .or. jj > cm%m%numLat+1) cycle
                        do smp = 1, cm%c(ii,jj)%tc(l)%numSample
                            k = cm%c(ii,jj)%tc(l)%sampleId(smp)
                            call TracerManager_GetSampleCoordinate(l, k, x2)
                            call MeshManager_CalcDistance(x1, x2, distance)
                            if (distance < eps) then
                                call TracerManager_GetSampleQuantity(l, k, cm%c(i,j)%q(l))
                                goto 101
                            end if
                            call RotationTransform(x1(1), x1(2), x2(1), x2(2), lonR=lonR)
                            quad = min(floor(lonR/dlonNormal)+1, numQuadNormal)
                            if (supSmpIdNormal(quad) == -1) then
                                supSmpIdNormal(quad) = k
                                supSmpDisNormal(quad) = distance
                            else
                                if (supSmpDisNormal(quad) > distance) then
                                    supSmpIdNormal(quad) = k
                                    supSmpDisNormal(quad) = distance
                                end if
                            end if
                        end do
                    end do
                    if (all(supSmpIdNormal /= -1)) exit
                    !if (lev == maxSearchLevel) then
                    !    call MsgManager_Speak(Warning, &
                    !        "No sufficient support samples in grid ("// &
                    !        trim(int2str(i))//","//trim(int2str(j))//").")
                    !end if
                end do
                cm%c(i,j)%q(l) = 0.0d0
                weightSum = 0.0d0
                do quad = 1, numQuadNormal
                    k = supSmpIdNormal(quad)
                    if (k == -1) cycle
                    call TracerManager_GetSampleQuantity(l, k, q)
                    weight = 1.0d0/supSmpDisNormal(quad)
                    weightSum = weightSum+weight
                    cm%c(i,j)%q(l) = cm%c(i,j)%q(l)+weight*q
                end do
#if (defined DEBUG)
                if (abs(weightSum) < eps) then
                    call MsgManager_Speak(Error, "weightSum equals zero!")
                    write(*, "('  Grid:  ', I3, ',', I3)") i, j
                    write(*, "('  Level: ', I3)") lev
                    call RunManager_EndRun
                end if
#endif
                cm%c(i,j)%q(l) = cm%c(i,j)%q(l)/weightSum
101             continue
            end do
        end do
        end do

        i = 1
        ! north pole
        j = 0
        do l = 1, n
            supSmpIdPole = -1
            do lev = 1, maxSearchLevel
                do pnt = 1, cm%m%numLon
                    ii = pnt
                    jj = lev-1
                    CheckSamples: do smp = 1, cm%c(ii,jj)%tc(l)%numSample
                        k = cm%c(ii,jj)%tc(l)%sampleId(smp)
                        call TracerManager_GetSampleCoordinate(l, k, x2)
                        distance = Re*(PI05-x2(2))
                        if (distance < eps) then
                            call TracerManager_GetSampleQuantity(l, k, cm%c(i,j)%q(l))
                            goto 102
                        end if
                        call RotationTransform(x1(1), x1(2), x2(1), x2(2), lonR=lonR)
                        quad = min(floor(lonR/dlonPole)+1, numQuadPole)
                        if (supSmpIdPole(quad) == -1) then
                            supSmpIdPole(quad) = k
                            supSmpDisPole(quad) = distance
                        else
                            if (supSmpDisPole(quad) > distance) then
                                supSmpIdPole(quad) = k
                                supSmpDisPole(quad) = distance
                            end if
                        end if
                    end do CheckSamples
                end do
                if (all(supSmpIdPole /= -1)) exit
                !if (lev == maxSearchLevel) then
                !    call MsgManager_Speak(Warning, &
                !        "No sufficient support samples ("// &
                !        trim(int2str(count(supSmpIdPole == -1)))// &
                !        ") in north pole.")
                !end if
            end do
            cm%c(i,j)%q(l) = 0.0d0
            weightSum = 0.0d0
            do quad = 1, numQuadPole
                k = supSmpIdPole(quad)
                if (k == -1) cycle
                call TracerManager_GetSampleQuantity(l, k, q)
                weight = 1.0d0/supSmpDisPole(quad)
                weightSum = weightSum+weight
                cm%c(i,j)%q(l) = cm%c(i,j)%q(l)+weight*q
            end do
#if (defined DEBUG)
            if (abs(weightSum) < eps) then
                call MsgManager_Speak(Error, "weightSum equals zero!")
                print *, supSmpDisPole
                print *, supSmpIdPole
                print *, i, j
                print *, lev
                call RunManager_EndRun
            end if
#endif
            cm%c(i,j)%q(l) = cm%c(i,j)%q(l)/weightSum
            do ii = 1, cm%m%numLon
                cm%c(ii,j)%q(l) = cm%c(i,j)%q(l)
            end do
102         continue
        end do

        ! south pole
        j = cm%m%numLat+1
        do l = 1, n
            supSmpIdPole = -1
            do lev = 1, maxSearchLevel
                do pnt = 1, cm%m%numLon
                    ii = pnt
                    jj = cm%m%numLat+2-lev
                    do smp = 1, cm%c(ii,jj)%tc(l)%numSample
                        k = cm%c(ii,jj)%tc(l)%sampleId(smp)
                        call TracerManager_GetSampleCoordinate(l, k, x2)
                        distance = Re*(PI05+x2(2))
                        if (distance < eps) then
                            call TracerManager_GetSampleQuantity(l, k, cm%c(i,j)%q(l))
                            goto 103
                        end if
                        call RotationTransform(x1(1), x1(2), x2(1), x2(2), lonR=lonR)
                        quad = min(floor(lonR/dlonPole)+1, numQuadPole)
                        if (supSmpIdPole(quad) == -1) then
                            supSmpIdPole(quad) = k
                            supSmpDisPole(quad) = distance
                        else
                            if (supSmpDisPole(quad) > distance) then
                                supSmpIdPole(quad) = k
                                supSmpDisPole(quad) = distance
                            end if
                        end if
                    end do
                end do
                if (all(supSmpIdPole /= -1)) exit
                !if (lev == maxSearchLevel) then
                !    call MsgManager_Speak(Warning, &
                !        "No sufficient support samples ("// &
                !        trim(int2str(count(supSmpIdPole == -1)))// &
                !        ") in south pole.")
                !end if
            end do
            cm%c(i,j)%q(l) = 0.0d0
            weightSum = 0.0d0
            do quad = 1, numQuadPole
                k = supSmpIdPole(quad)
                if (k == -1) cycle
                call TracerManager_GetSampleQuantity(l, k, q)
                weight = 1.0d0/supSmpDisPole(quad)
                weightSum = weightSum+weight
                cm%c(i,j)%q(l) = cm%c(i,j)%q(l)+weight*q
            end do
#if (defined DEBUG)
            if (abs(weightSum) < eps) then
                call MsgManager_Speak(Error, "weightSum equals zero!")
                print *, supSmpDisPole
                print *, supSmpIdPole
                print *, i, j
                print *, lev
                call RunManager_EndRun
            end if
#endif
            cm%c(i,j)%q(l) = cm%c(i,j)%q(l)/weightSum
            do ii = 1, cm%m%numLon
                cm%c(ii,j)%q(l) = cm%c(i,j)%q(l)
            end do
103         continue
        end do

        ! reset
        do j = 1, cm%m%numLat
        do i = 1, cm%m%numLon
            cm%c(i,j)%tc(:)%numSample = 0
        end do
        end do

        call MsgManager_DeleteSpeaker
    
    end subroutine LagrangeToEuler_IDWInterp
 
    subroutine LagrangeToEuler_Output(timeStep, time)
        integer, intent(in) :: timeStep
        real(8), intent(in) :: time
        
        real(8) q(cm%m%numLon,0:cm%m%numLat+1)
        integer i, j, k, n
        character(20) name

        n = TracerManager_GetTracerClassNumber()

        call NFWrap_Advance(fcard, timeStep, time)
        do k = 1, n
            do j = 0, cm%m%numLat+1
                do i = 1, cm%m%numLon
                    q(i,j) = cm%c(i,j)%q(k)
                end do
            end do
            name = TracerManager_GetTracerName(k)
            call NFWrap_Output2DVar(fcard, name, q)
        end do

    end subroutine LagrangeToEuler_Output
   
end module LagrangeToEuler
