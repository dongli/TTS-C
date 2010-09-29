program CalcError_1

    use MsgManager
    use FloatingPoint
    use NFWrap
    use netcdf

    implicit none

    integer ncid(2), ierr
    integer timeDimId(2), lonDimId(2), latDimId(2)
    integer timeVarId(2), lonVarId(2), latVarId(2)
    integer areaVarId, qVarId(2)
    integer errVarId, L1ErrorVarId, L2ErrorVarId, LInfErrorVarId
    integer :: ncmode(2) = [nf90_nowrite,nf90_write]

    character(256) :: fileName(2)

    integer numTimeStep(2), numLon(2), numLat(2)
    real(8), allocatable :: time(:,:), area(:,:)
    real(8), allocatable :: lon(:,:), lat(:,:), q(:,:,:,:)

    real(8), allocatable :: err(:,:,:)
    real(8), allocatable :: L1Error(:), L2Error(:), LInfError(:)

    real(8) dq, q1, q2, qInf
    integer i, j, k

    call MsgManager_RecordSpeaker("CalcError_1")

    call get_command_argument(1, fileName(1))
    call get_command_argument(2, fileName(2))

    do i = 1, 2
        ierr = nf90_open(fileName(i), ncmode(i), ncid(i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_inq_dimid(ncid(i), "time", timeDimId(i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_inquire_dimension(ncid(i), timeDimId(i), len=numTimeStep(i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_inq_dimid(ncid(i), "lon", lonDimId(i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_inquire_dimension(ncid(i), lonDimId(i), len=numLon(i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_inq_dimid(ncid(i), "lat", latDimId(i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_inquire_dimension(ncid(i), latDimId(i), len=numLat(i))
        call NFWrap_HandleError(ierr)
    end do

    if (numTimeStep(1) /= numTimeStep(2)) then
        write(*, *) "Time step numbers do not match."
        stop
    end if

    if (numLon(1) /= numLon(2)) then
        write(*, *) "Longitude grid numbers do not match."
        stop
    end if

    if (numLat(1) /= numLat(2)) then
        write(*, *) "Latitude grid numbers do not match."
        stop
    end if

    allocate(time(numTimeStep(1),2))
    allocate(area(numLon(1),numLat(1)))
    allocate(lon(numLon(1),2))
    allocate(lat(numLat(1),2))
    allocate(q(numLon(1),numLat(1),numTimeStep(1),2))

    do i = 1, 2
        ierr = nf90_inq_varid(ncid(i), "time", timeVarId(i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_get_var(ncid(i), timeVarId(i), time(:,i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_inq_varid(ncid(i), "lon", lonVarId(i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_get_var(ncid(i), lonVarId(i), lon(:,i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_inq_varid(ncid(i), "lat", latVarId(i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_get_var(ncid(i), latVarId(i), lat(:,i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_inq_varid(ncid(i), "q", qVarId(i))
        call NFWrap_HandleError(ierr)
        ierr = nf90_get_var(ncid(i), qVarId(i), q(:,:,:,i))
        call NFWrap_HandleError(ierr)
    end do

    ierr = nf90_inq_varid(ncid(2), "area", areaVarId)
    call NFWrap_HandleError(ierr)
    ierr = nf90_get_var(ncid(2), areaVarId, area)
    call NFWrap_HandleError(ierr)

    do k = 1, numTimeStep(1)
        if (abs(time(k,1)-time(k,2)) > eps) then
            write(*, *) "Times do not match."
            stop
        end if
    end do

    do i = 1, numLon(1)
        if (abs(lon(i,1)-lon(i,2)) > eps) then
            write(*, *) "Longitudes do not match."
            stop
        end if
    end do

    do j = 1, numLat(1)
        if (abs(lat(j,1)-lat(j,2)) > eps) then
            write(*, *) "Latitudes do not match."
            stop
        end if
    end do

    allocate(err(numLon(1),numLat(1),numTimeStep(1)))
    allocate(L1Error(numTimeStep(1)))
    allocate(L2Error(numTimeStep(1)))
    allocate(LInfError(numTimeStep(1)))

    do k = 1, numTimeStep(1)
        L1Error(k) = 0.0d0
        L2Error(k) = 0.0d0
        LInfError(k) = 0.0d0
        q1 = 0.0d0
        q2 = 0.0d0
        qInf = 0.0d0
        do j = 1, numLat(1)
            do i = 1, numLon(1)
                dq = abs(q(i,j,k,2)-q(i,j,k,1))
                err(i,j,k) = q(i,j,k,2)-q(i,j,k,1)!dq/q(i,j,k,1)*100
                L1Error(k) = L1Error(k)+dq*area(i,j)
                L2Error(k) = L2Error(k)+dq**2.0d0*area(i,j)
                if (dq > LInfError(k)) LInfError(k) = dq
                q1 = q1+q(i,j,k,1)*area(i,j)
                q2 = q2+q(i,j,k,1)**2.0d0*area(i,j)
                if (q(i,j,k,1) > qInf) qInf = q(i,j,k,1)
            end do
        end do
        L1Error(k) = L1Error(k)/q1
        L2Error(k) = sqrt(L2Error(k)/q2)
        LInfError(k) = LInfError(k)/qInf
    end do

    ierr = nf90_inq_varid(ncid(2), "L1Error", L1ErrorVarId)
    if (ierr /= nf90_noerr) then
        ierr = nf90_redef(ncid(2))
        call NFWrap_HandleError(ierr)
        ierr = nf90_def_var(ncid(2), "L1Error", nf90_double, &
            [timeDimId(2)], L1ErrorVarId)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(ncid(2), L1ErrorVarId, &
            "long_name", "L1 normalized standard error")
        call NFWrap_HandleError(ierr)
        ierr = nf90_def_var(ncid(2), "L2Error", nf90_double, &
            [timeDimId(2)], L2ErrorVarId)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(ncid(2), L2ErrorVarId, &
            "long_name", "L2 normalized standard error")
        call NFWrap_HandleError(ierr)
        ierr = nf90_def_var(ncid(2), "LInfError", nf90_double, &
            [timeDimId(2)], LInfErrorVarId)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(ncid(2), LInfErrorVarId, &
            "long_name", "L_infinite normalized standard error")
        call NFWrap_HandleError(ierr)
        ierr = nf90_enddef(ncid(2))
        call NFWrap_HandleError(ierr)
    else
        ierr = nf90_inq_varid(ncid(2), "L1Error", L1ErrorVarId)
        call NFWrap_HandleError(ierr)
        ierr = nf90_inq_varid(ncid(2), "L2Error", L2ErrorVarId)
        call NFWrap_HandleError(ierr)
        ierr = nf90_inq_varid(ncid(2), "LInfError", LInfErrorVarId)
        call NFWrap_HandleError(ierr)
    end if

    ierr = nf90_inq_varid(ncid(2), "error", errVarId)
    if (ierr /= nf90_noerr) then
        ierr = nf90_redef(ncid(2))
        call NFWrap_HandleError(ierr)
        ierr = nf90_def_var(ncid(2), "error", nf90_double, &
            [lonDimId(2),latDimId(2),timeDimId(2)], errVarId)
        call NFWrap_HandleError(ierr)
        !ierr = nf90_put_att(ncid(2), errVarId, &
        !    "long_name", "relative error")
        !call NFWrap_HandleError(ierr)
        !ierr = nf90_put_att(ncid(2), errVarId, &
        !    "units", "percent")
        !call NFWrap_HandleError(ierr)
        ierr = nf90_enddef(ncid(2))
        call NFWrap_HandleError(ierr)
    else
        ierr = nf90_inq_varid(ncid(2), "error", errVarId)
        call NFWrap_HandleError(ierr)
    end if
    ierr = nf90_put_var(ncid(2), errVarId, err)
    call NFWrap_HandleError(ierr)
    ierr = nf90_put_var(ncid(2), L1ErrorVarId, L1Error)
    call NFWrap_HandleError(ierr)
    ierr = nf90_put_var(ncid(2), L2ErrorVarId, L2Error)
    call NFWrap_HandleError(ierr)
    ierr = nf90_put_var(ncid(2), LInfErrorVarId, LInfError)
    call NFWrap_HandleError(ierr)

    do i = 1, 2
        ierr = nf90_close(ncid(i))
        call NFWrap_HandleError(ierr)
    end do

    call MsgManager_Speak(Notice, "Finished.")
    call MsgManager_DeleteSpeaker

end program CalcError_1
