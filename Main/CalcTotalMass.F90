program CalcTotalMass

    use MsgManager
    use RunManager
    use NFWrap
    use netcdf

    implicit none

    integer ncid, ierr
    integer timeDimId, lonDimId, latDimId
    integer timeVarId, lonVarId, latVarId
    integer areaVarId, qVarId, massVarId

    character(256) fileName

    integer numTimeStep, numLon, numLat
    real(8), allocatable :: area(:,:)
    real(8), allocatable :: q(:,:,:)

    real(8), allocatable :: totalMass(:)

    integer i, j, k

    call MsgManager_RecordSpeaker("CalcTotalMass")

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
    allocate(q(numLon,numLat,numTimeStep))

    ierr = nf90_inq_varid(ncid, "q", qVarId)
    call NFWrap_HandleError(ierr)
    ierr = nf90_get_var(ncid, qVarId, q)
    call NFWrap_HandleError(ierr)

    ierr = nf90_inq_varid(ncid, "area", areaVarId)
    call NFWrap_HandleError(ierr)
    ierr = nf90_get_var(ncid, areaVarId, area)
    call NFWrap_HandleError(ierr)

    allocate(totalMass(numTimeStep))
    totalMass = 0.0d0

    do k = 1, numTimeStep
        do j = 1, numLat
        do i = 1, numLon
            totalMass(k) = totalMass(k)+q(i,j,k)*area(i,j)
        end do
        end do
    end do

    ierr = nf90_inq_varid(ncid, "mass", massVarId)
    if (ierr /= nf90_noerr) then
        ierr = nf90_redef(ncid)
        call NFWrap_HandleError(ierr)
        ierr = nf90_def_var(ncid, "mass", nf90_double, &
            [timeDimId], massVarId)
        call NFWrap_HandleError(ierr)
        ierr = nf90_put_att(ncid, massVarId, &
            "long_name", "total mass")
        call NFWrap_HandleError(ierr)
        ierr = nf90_enddef(ncid)
        call NFWrap_HandleError(ierr)
    else
        ierr = nf90_inq_varid(ncid, "mass", massVarId)
        call NFWrap_HandleError(ierr)
    end if
    ierr = nf90_put_var(ncid, massVarId, totalMass)
    call NFWrap_HandleError(ierr)

    ierr = nf90_close(ncid)
    call NFWrap_HandleError(ierr)

    call MsgManager_Speak(Notice, "Finished.")
    call MsgManager_DeleteSpeaker

    call RunManager_EndRun

end program CalcTotalMass
