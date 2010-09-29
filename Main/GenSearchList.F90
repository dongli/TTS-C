! *************************************************************************** !
! GenSearchList program                                                       !
!                                                                             !
! Description:                                                                !
!                                                                             !
! Author:                                                                     !
!                                                                             !
! *************************************************************************** !

program GenSearchList

    use MsgManager
    use RunManager
    use MeshManager

    implicit none

    integer numLon, numLat

    integer maxSearchLevel
    character(100) nlFilePath, outputFilePath
    real(8) sr1, sr2, dr, x1(2), x2(2), r
    logical, allocatable :: checked(:,:)
    integer, parameter :: maxNumGrid = 1000
    integer itmp(maxNumGrid), jtmp(maxNumGrid)
    integer i, j, k, l, n, lev

    type(Mesh), pointer :: m

    namelist /MeshParams/ numLon, numLat

    call MsgManager_RecordSpeaker("GenSearchList")

    call get_command_argument(1, nlFilePath)
    call get_command_argument(2, outputFilePath)

    open(10, file=nlFilePath)
    read(10, nml=MeshParams)
    close(10)

    call MeshManager_Init(numLon, numLat)

    m => MeshSelector("Full mesh")

    maxSearchLevel = 5

    open(10, file=outputFilePath, form="unformatted")

    write(10) maxSearchLevel

    allocate(checked(m%numLon,0:m%numLat+1))
    do j = 1, m%numLat
    print *, "latitude ", j
    do i = 1, m%numLon
        sr1 = -0.01d0
        sr2 = Re*dimInfo%dlat
        dr = Re*dimInfo%dlat
        checked = .false.
        x1 = [m%lon(i),m%lat(j)]
        do lev = 1, maxSearchLevel
            n = 0
            do l = 0, m%numLat+1 ! include the polar caps
            do k = 1, m%numLon
                if (checked(k,l)) cycle
                if (l == 0) then
                    x2 = [m%lon(k),PI05-dimInfo%dlat05]
                else if (l == m%numLat+1) then
                    x2 = [m%lon(k),-PI05+dimInfo%dlat05]
                else
                    x2 = [m%lon(k),m%lat(l)]
                end if
                call MeshManager_CalcDistance(x1, x2, r)
                if (r > sr1 .and. r <= sr2) then
                    checked(k,l) = .true.
                    n = n+1
                    if (n > maxNumGrid) then
                        call MsgManager_Speak(Error, &
                            "Exceed maximum number of temporary "// &
                            "search grids, try to enlarge it.")
                        call RunManager_EndRun
                    end if
                    itmp(n) = k
                    jtmp(n) = l
                end if
            end do
            end do
            write(10) n
            write(10) itmp(1:n)
            write(10) jtmp(1:n)
            sr1 = sr2
            sr2 = sr2+dr
        end do
    end do
    end do

    close(10)

    call MsgManager_Speak(Notice, &
        "Binary file """//trim(outputFilePath)//""" has been generated.")

    call MsgManager_Speak(Notice, "Finished.")
    call MsgManager_DeleteSpeaker

end program GenSearchList

