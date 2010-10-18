program CalcError_2

    use MsgManager
    use FloatingPoint
    use NFWrap
    use netcdf

    implicit none

    integer ncid, ierr
    integer timeDimId, lonDimId, latDimId
    integer timeVarId, lonVarId, latVarId
    integer areaVarId, qVarId
    integer errVarId, L1ErrorVarId, L2ErrorVarId, LInfErrorVarId

    character(256) fileName
    integer numTimeStep, numLon, numLat
    real(8), allocatable :: area(:,:)
    real(8), allocatable :: q(:,:,:)
    real(8), allocatable :: err(:,:)
    real(8) L1Error, L2Error, LInfError
    real(8) dq, q1, q2, qInf
    integer i, j, k

    call MsgManager_RecordSpeaker("CalcError_2")

    call get_command_argument(1, fileName)

    ierr = nf90_open(fileName, nf90_write, ncid)
    call NFWrap_HandleError(ierr)
    ierr = nf90_inq_dimid(ncid, "time", timeDimId)
    call NFWrap_HandleError(ierr)
    ierr = nf90_inquire_dimension(ncid, timeDimId, len=numTimeStep)
    call NFWrap_HandleError(ierr)
    ierr = nf90_inq_dimid(ncid, "lon", lonDimId)
    call NFWrap_HandleError(ierr)
    ierr = nf90_inquire_dimension(ncid, lonDimId, len=numLon)
    call NFWrap_HandleError(ierr)
    ierr = nf90_inq_dimid(ncid, "lat", latDimId)
    call NFWrap_HandleError(ierr)
    ierr = nf90_inquire_dimension(ncid, latDimId, len=numLat)
    call NFWrap_HandleError(ierr)

    allocate(area(numLon,numLat))
    allocate(q(numLon,numLat,2))

    ierr = nf90_inq_varid(ncid, "q", qVarId)
    call NFWrap_HandleError(ierr)
    ierr = nf90_get_var(ncid, qVarId, q(:,:,1), [1,1,1], [numLon,numLat,1])
    call NFWrap_HandleError(ierr)
    ierr = nf90_get_var(ncid, qVarId, q(:,:,2), [1,1,numTimeStep], [numLon,numLat,1])
    call NFWrap_HandleError(ierr)
    ierr = nf90_inq_varid(ncid, "area", areaVarId)
    call NFWrap_HandleError(ierr)
    ierr = nf90_get_var(ncid, areaVarId, area)
    call NFWrap_HandleError(ierr)

    allocate(err(numLon,numLat))

    L1Error = 0.0d0
    L2Error = 0.0d0
    LInfError = 0.0d0
    q1 = 0.0d0
    q2 = 0.0d0
    qInf = 0.0d0
    do j = 1, numLat
        do i = 1, numLon
            dq = abs(q(i,j,2)-q(i,j,1))
            err(i,j) = q(i,j,2)-q(i,j,1)!dq/q(i,j,k,1)*100
            L1Error = L1Error+dq*area(i,j)
            L2Error = L2Error+dq**2.0d0*area(i,j)
            if (dq > LInfError) LInfError = dq
            q1 = q1+q(i,j,1)*area(i,j)
            q2 = q2+q(i,j,1)**2.0d0*area(i,j)
            if (q(i,j,1) > qInf) qInf = q(i,j,1)
        end do
    end do
    L1Error = L1Error/q1
    L2Error = sqrt(L2Error/q2)
    LInfError = LInfError/qInf

    write(*, "('L1Error   = ', F20.15)") L1Error
    write(*, "('L2Error   = ', F20.15)") L2Error
    write(*, "('LInfError = ', F20.15)") LInfError

    ierr = nf90_inq_varid(ncid, "error", errVarId)
    if (ierr /= nf90_noerr) then
        ierr = nf90_redef(ncid)
        call NFWrap_HandleError(ierr)
        ierr = nf90_def_var(ncid, "error", nf90_double, &
            [lonDimId,latDimId], errVarId)
        call NFWrap_HandleError(ierr)
        !ierr = nf90_put_att(ncid, errVarId, &
        !    "long_name", "relative error")
        !call NFWrap_HandleError(ierr)
        !ierr = nf90_put_att(ncid, errVarId, &
        !    "units", "percent")
        !call NFWrap_HandleError(ierr)
        ierr = nf90_enddef(ncid)
        call NFWrap_HandleError(ierr)
    else
        ierr = nf90_inq_varid(ncid, "error", errVarId)
        call NFWrap_HandleError(ierr)
    end if
    ierr = nf90_put_var(ncid, errVarId, err)
    call NFWrap_HandleError(ierr)

    ierr = nf90_close(ncid)

    call MsgManager_Speak(Notice, "Finished.")
    call MsgManager_DeleteSpeaker

end program CalcError_2
