module TracerManager

    use CommonTypes
    use MsgManager
    use RunManager
    use MeshManager
    use NFWrap
    use netcdf

    implicit none

    private

    public TracerManager_Init
    public TracerManager_Final

    public TracerManager_RegisterOperation
    public TracerManager_RegisterTracer

    public TracerManager_InitialCondition

    public TracerManager_Output

    public TracerManager_GetMaxTracerClassNumber
    public TracerManager_GetTracerClassNumber
    public TracerManager_GetTracerName

    public TracerManager_GetSampleNumber
    public TracerManager_GetSampleCoordinate
    public TracerManager_GetSampleQuantity
    public TracerManager_GetSampleLocation
    public TracerManager_SaveSampleCoordinate
    public TracerManager_SaveSampleQuantity
    public TracerManager_SaveSampleLocation
    public TracerManager_SaveSampleVelocity

    public TracerManager_IsAdvective
    public TracerManager_IsConservative

    public Conservative, Advective

    integer, parameter :: Conservative = 1
    integer, parameter :: Advective    = 2

    ! ======================================================================== !
    !                          MAIN DATA STRUCTURE
    type TracerSample ! ...................................................... 1
        real(8) xOld(3)
        real(8) xNew(3)
        real(8) qOld
        real(8) qNew
        real(8) v(3), div
        type(Location) locOld
        type(Location) locNew
    end type TracerSample

    type TracerClass ! ....................................................... 2
        integer id
        character(20) name
        integer type
        integer number
        type(TracerSample), allocatable :: samples(:)
    end type TracerClass

    ! ======================================================================== !
    !                                DATA
    integer maxNumTracerClass
    integer numTracerClass
    type(TracerClass), allocatable, target :: tracers(:)

    ! ======================================================================== !
    !                           AUXILLARY DATA
    type(FileCard) fcard

    character(10), parameter :: velocityShortName(3) = ["u", "v", "w"]

    type(OperationList) RegisterTracer
    type(OperationList), target :: InitialCondition

contains

    subroutine TracerManager_Init(maxNum, filePath)
        integer, intent(in) :: maxNum
        character(*), intent(in) :: filePath

        procedure(), pointer :: p

        call MsgManager_RecordSpeaker("TracerManager_Init")

        ! ==================================================================== !
        !                         MEMORY ALLOCATION
        maxNumTracerClass = maxNum
        allocate(tracers(maxNumTracerClass))
        numTracerClass = 0

        ! ==================================================================== !
        !                           NETCDF PART
        call NFWrap_CreateIrregular(filePath, fcard)

        p => TracerManager_Final
        call RunManager_RegisterOperation("EndRun", "TracerManager", "Final", p)

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine TracerManager_Init

    subroutine TracerManager_Final
        integer i, j

        call MsgManager_RecordSpeaker("TracerManager_Final")

        do i = 1, numTracerClass
            deallocate(tracers(i)%samples)
        end do
        deallocate(tracers)

        call NFWrap_Close(fcard)

        call MsgManager_Speak(Notice, "Finished.")
        call MsgManager_DeleteSpeaker

    end subroutine TracerManager_Final

    subroutine TracerManager_RegisterOperation(hostName, modName, subName, handle)
        character(*), intent(in) :: hostName, modName, subName
        procedure(), intent(in), pointer :: handle

        call MsgManager_RecordSpeaker("TracerManager_RegisterOperation")

        select case (hostName)
        case ("InitialCondition")
            call InitialCondition%register(modName, subName, handle)
        case ("RegisterTracer")
            call RegisterTracer%register(modName, subName, handle)
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

    end subroutine TracerManager_RegisterOperation

    subroutine TracerManager_RegisterTracer(varName, number, type, &
        longName, unitName, tracerId)
        character(*), intent(in) :: varName
        integer, intent(in) :: number, type
        character(*), intent(in) :: longName, unitName
        integer, intent(out) :: tracerId

        type(SubHandle), pointer :: sh
        integer i

        call MsgManager_RecordSpeaker("TracerManager_RegisterTracer")

        ! validation of registration
        if (CheckTracerExistence(trim(varName))) then
            call MsgManager_Speak(Warning, &
                "Tracer "//trim(varName)//" has been registered already.")
            return
        end if
        numTracerClass = numTracerClass+1
        if (numTracerClass > maxNumTracerClass) then
            call MsgManager_Speak(Error, &
                "Exceed maximum number of tracer classes!")
            call RunManager_EndRun
        end if
        if (type /= Conservative .and. type /= Advective) then
            call MsgManager_Speak(Error, &
                "Invalid tracer type, choose Conservative or Advective.")
            call RunManager_EndRun
        end if

        ! fill the meta data
        tracerId = numTracerClass
        tracers(tracerId)%name = varName
        tracers(tracerId)%number = number
        tracers(tracerId)%type = type
        allocate(tracers(tracerId)%samples(number))

        ! ==================================================================== !
        call NFWrap_NewDim(fcard, trim(varName)//"_num", number)
        do i = 1, numDim
            call NFWrap_New1DVar(fcard, &
                trim(varName)//"_"//trim(coordinateShortName(i)), &
                trim(varName)//"_num", coordinateLongName(i), &
                coordinateUnit(i), timeVariant=.true.)
        end do
        call NFWrap_New1DVar(fcard, varName, trim(varName)//"_num", &
            longName, unitName, timeVariant=.true.)

        sh => RegisterTracer%head
        do i = 1, RegisterTracer%num
            call sh%handle(varName, longName, unitName)
            sh => sh%next
        end do

        call MsgManager_Speak(Notice, &
            "New tracer """//trim(varName)//""" is added.")
        call MsgManager_DeleteSpeaker

    end subroutine TracerManager_RegisterTracer

    subroutine TracerManager_InitialCondition(tracerId, n, x, q)
        integer, intent(in) :: tracerId, n
        real(8), intent(in) :: x(numDim,n)
        real(8), intent(in) :: q(n)

        type(SubHandle), pointer :: sh
        integer i

        call MsgManager_RecordSpeaker("TracerManager_InitialCondition")

        do i = 1, tracers(tracerId)%number
            tracers(tracerId)%samples(i)%xOld(1:numDim) = x(:,i)
            tracers(tracerId)%samples(i)%xNew(1:numDim) = x(:,i)
            tracers(tracerId)%samples(i)%qOld = q(i)
            tracers(tracerId)%samples(i)%qNew = q(i)
#if (defined LAGRANGE_TO_EULER)
            call MeshManager_LocationCheck(x(:,i), &
                tracers(tracerId)%samples(i)%locNew, tracerId, i)
#else
            call MeshManager_LocationCheck(x(:,i), &
                tracers(tracerId)%samples(i)%locNew)
#endif
            tracers(tracerId)%samples(i)%locOld = &
                tracers(tracerId)%samples(i)%locNew
        end do

        sh => InitialCondition%head
        do i = 1, InitialCondition%num
            call sh%handle(tracerId, n, x, q)
            sh => sh%next
        end do

        call MsgManager_Speak(Notice, &
            "Tracer """//trim(tracers(tracerId)%name)//""" is initialized.")
        call MsgManager_DeleteSpeaker

    end subroutine TracerManager_InitialCondition

    integer function TracerManager_GetMaxTracerClassNumber() result(res)
    
        res = maxNumTracerClass
    
    end function TracerManager_GetMaxTracerClassNumber
    
    integer function TracerManager_GetTracerClassNumber() result(res)
    
        res = numTracerClass
    
    end function TracerManager_GetTracerClassNumber
    
    character(20) function TracerManager_GetTracerName(tracerId) result(res)
        integer, intent(in) :: tracerId
        res = tracers(tracerId)%name
    end function TracerManager_GetTracerName

    subroutine TracerManager_GetSampleNumber(tracerId, numSample)
        integer, intent(in) :: tracerId
        integer, intent(out) :: numSample

        numSample = tracers(tracerId)%number

    end subroutine TracerManager_GetSampleNumber

    subroutine TracerManager_GetSampleCoordinate(tracerId, sampleId, x)
        integer, intent(in) :: tracerId, sampleId
        real(8), intent(out) :: x(:)

        x = tracers(tracerId)%samples(sampleId)%xNew(1:numDim)

    end subroutine TracerManager_GetSampleCoordinate

    subroutine TracerManager_GetSampleQuantity(tracerId, sampleId, q)
        integer, intent(in) :: tracerId, sampleId
        real(8), intent(out) :: q

        q = tracers(tracerId)%samples(sampleId)%qNew

    end subroutine TracerManager_GetSampleQuantity

    subroutine TracerManager_GetSampleLocation(tracerId, sampleId, loc)
        integer, intent(in) :: tracerId, sampleId
        type(Location), intent(out) :: loc

        loc = tracers(tracerId)%samples(sampleId)%locNew

    end subroutine TracerManager_GetSampleLocation

    subroutine TracerManager_Output(timeStep, curTime)
        integer, intent(in) :: timeStep
        real(8), intent(in) :: curTime

        character(20) name
        integer i, j

        call MsgManager_RecordSpeaker("TracerManager_Output")

        call NFWrap_Advance(fcard, timeStep, curTime)
        do i = 1, numTracerClass
            name = tracers(i)%name
            do j = 1, numDim
                call NFWrap_Output1DVar(fcard, &
                    trim(name)//"_"//trim(coordinateShortName(j)), &
                    tracers(i)%samples(:)%xNew(j))
            end do
            call NFWrap_Output1DVar(fcard, name, tracers(i)%samples(:)%qNew)
        end do

        call MsgManager_DeleteSpeaker

    end subroutine TracerManager_Output

    subroutine TracerManager_SaveSampleCoordinate(tracerId, sampleId, x)
        integer, intent(in) :: tracerId, sampleId
        real(8), intent(in) :: x(:)

        tracers(tracerId)%samples(sampleId)%xOld(1:numDim) = &
            tracers(tracerId)%samples(sampleId)%xNew(1:numDim)
        tracers(tracerId)%samples(sampleId)%xNew(1:numDim) = x(1:numDim)

    end subroutine TracerManager_SaveSampleCoordinate

    subroutine TracerManager_SaveSampleQuantity(tracerId, sampleId, q)
        integer, intent(in) :: tracerId, sampleId
        real(8), intent(in) :: q

        tracers(tracerId)%samples(sampleId)%qOld = &
            tracers(tracerId)%samples(sampleId)%qNew
        tracers(tracerId)%samples(sampleId)%qNew = q

    end subroutine TracerManager_SaveSampleQuantity

    subroutine TracerManager_SaveSampleLocation(tracerId, sampleId, loc)
        integer, intent(in) :: tracerId, sampleId
        type(Location), intent(in) :: loc

        tracers(tracerId)%samples(sampleId)%locOld = &
            tracers(tracerId)%samples(sampleId)%locNew
        tracers(tracerId)%samples(sampleId)%locNew = loc

    end subroutine TracerManager_SaveSampleLocation

    subroutine TracerManager_SaveSampleVelocity(tracerId, sampleId, v)
        integer, intent(in) :: tracerId, sampleId
        real(8), intent(in) :: v(numDim)

        tracers(tracerId)%samples(sampleId)%v(:numDim) = v

    end subroutine TracerManager_SaveSampleVelocity

    logical function TracerManager_IsConservative(tracerId) result(res)
        integer, intent(in) :: tracerId

        res = .false.
        if (tracers(tracerId)%type == Conservative) res = .true.

    end function TracerManager_IsConservative

    logical function TracerManager_IsAdvective(tracerId) result(res)
        integer, intent(in) :: tracerId

        res = .false.
        if (tracers(tracerId)%type == Advective) res = .true.

    end function TracerManager_IsAdvective

    ! ************************************************************************ !
    ! INTERNAL SUBROUTINES (LEAVE THEM ALONE)
    ! ************************************************************************ !

    logical function CheckTracerExistence(name) result(existence)
        character(*), intent(in) :: name
\
        integer i

        do i = 1, numTracerClass
            if (tracers(i)%name == name) then
                existence = .true.
                return
            end if
        end do
        existence = .false.

    end function CheckTracerExistence

end module TracerManager
