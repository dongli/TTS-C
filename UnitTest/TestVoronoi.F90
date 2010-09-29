program TestVoronoi

    use MsgManager
    use SphereService
    use NFWrap

    implicit none

    type(FileCard) fcard
    integer numTracer
    real(8), allocatable :: lon(:), lat(:)

    integer i, j

    ! Voronoi cell vertices
    integer numVtx ! for each lon and lat
    integer, allocatable :: idxVtx(:)
    real(4), allocatable :: lonVtx(:), latVtx(:)

    ! circumradii
    integer nca ! actual number of circumcenters in lonVtx and latVtx
    real(8), allocatable :: rc(:)

    ! Internal variables for CSVORO
    integer, allocatable :: iwork(:)
    real(8), allocatable :: rwork(:)
    integer nc, ierr

    call NFWrap_OpenForRead("tracer.nc", fcard)
    call NFWrap_GetDimSize(fcard, "q_num", numTracer)
    allocate(lon(numTracer), lat(numTracer))
    call NFWrap_Input1DVar(fcard, "q_lon", lon, 1)
    call NFWrap_Input1DVar(fcard, "q_lat", lat, 1)
    lon = lon*Rad2Deg
    lat = lat*Rad2Deg
    write(*, "('sample of lon: ', 5F10.5)") lon(1:5)
    write(*, "('sample of lat: ', 5F10.5)") lat(1:5)

    allocate(iwork(numTracer*27))
    allocate(rwork(numTracer*13))
    nc = 2*numTracer
    allocate(lonVtx(nc), latVtx(nc))
    allocate(rc(nc))
    allocate(idxVtx(numTracer))

    call CSVORO(numTracer, real(lat), real(lon), 1, 1, iwork, rwork, nc, &
        latVtx, lonVtx, rc, nca, numVtx, idxVtx, ierr)
    do i = 1, numTracer
        call CSVORO(numTracer, real(lat), real(lon), i, 0, iwork, rwork, nc, &
            latVtx, lonVtx, rc, nca, numVtx, idxVtx, ierr)
        print *, "numVtx = ", numVtx
        print *, "Vtx = "
        do j = 1, numVtx
            print *, idxVtx(j), lonVtx(idxVtx(j)), latVtx(idxVtx(j))
        end do
        read(*, *)
    end do

end program TestVoronoi

