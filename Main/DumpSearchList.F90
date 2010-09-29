! *************************************************************************** !
! DumpSearchList program                                                      !
!                                                                             !
! Description:                                                                !
!                                                                             !
! Author:                                                                     !
!                                                                             !
! *************************************************************************** !

program DumpSearchList

    use MsgManager
    use RunManager
    use SphereService

    implicit none

    character(100) namelistFile, outputFile

    integer numLon, numLat

    integer maxSearchLevel
    integer numPoint
    integer, allocatable :: itmp(:), jtmp(:)
    logical flag
    integer i, j, k, ic, jc, ii, jj

    namelist /MeshParams/ numLon, numLat

    call MsgManager_RecordSpeaker("DumpSearchList")

    call get_command_argument(1, namelistFile)
    call get_command_argument(2, outputFile)

    inquire(file=namelistFile, exist=flag)
    if (.not. flag) then
        call MsgManager_Speak(Error, &
            "File """//trim(namelistFile)//""" does not exist.")
        call RunManager_EndRun
    end if
    open(10, file=namelistFile)
    read(10, nml=MeshParams)
    close(10)

    inquire(file=outputFile, exist=flag)
    if (.not. flag) then
        call MsgManager_Speak(Error, &
            "File """//trim(outputFile)//""" does not exist.")
        call RunManager_EndRun
    end if
    open(10, file=outputFile, form="unformatted", action="read")
    read(10) maxSearchLevel
    write(*, "('Maximum search level: ', I5)") maxSearchLevel

    write(*, "('Input a grid index for printing:')")
    read(*, *) ic, jc

    do j = 1, jc-1
        do i = 1, numLon
            do k = 1, maxSearchLevel
                read(10) numPoint
                allocate(itmp(numPoint))
                read(10) itmp
                deallocate(itmp)
                allocate(jtmp(numPoint))
                read(10) jtmp
                deallocate(jtmp)
            end do
        end do
    end do
    do i = 1, ic-1
        do k = 1, maxSearchLevel
            read(10) numPoint
            allocate(itmp(numPoint))
            read(10) itmp
            deallocate(itmp)
            allocate(jtmp(numPoint))
            read(10) jtmp
            deallocate(jtmp)
        end do
    end do

    do k = 1, maxSearchLevel
        write(*, "('Search level ', I5)") k
        read(10) numPoint
        allocate(itmp(numPoint))
        read(10) itmp
        allocate(jtmp(numPoint))
        read(10) jtmp
        write(*, "(' number of points:', I5)") numPoint
        do i = 1, numPoint
            write(*, "('  ---> ', 2I5)") itmp(i), jtmp(i)
        end do
        read(*, *)
        deallocate(itmp)
        deallocate(jtmp)
    end do

    close(10)

    call MsgManager_Speak(Notice, "Finished.")
    call MsgManager_DeleteSpeaker

    call RunManager_EndRun

end program DumpSearchList

