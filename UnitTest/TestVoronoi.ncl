load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  f = addfile("tracer.nc", "r")

  dims = getfilevardimsizes(f, "q")

  numTime = dims(0)
  numTracer = dims(1)

  idx = ispan(0, numTracer-1, 1)
  numSparseTracer = dimsizes(idx)

  ; constants
  PI = 4.*atan(1.)
  fac = 180./PI

  ; create working variables
  lonOut = new(2*numSparseTracer, float)
  latOut = new(2*numSparseTracer, float)
  circumradius = new(2*numSparseTracer, float)
  numActual = new(1, integer)
  numVertex = new(1, integer)
  idxVertex = new(numSparseTracer, integer)
  lon = new(2, float)
  lat = new(2, float)

  ; configure graphics
  wks = gsn_open_wks("pdf", "tracer_voronoi")
  
  mapres = True
  mapres@gsnFrame = False
  mapres@mpGreatCircleLinesOn = True
  mapres@mpGridAndLimbOn = True
  mapres@mpGridLineColor = "blue"
  mapres@mpGridSpacingF = 10.0
  ;mapres@mpProjection = "Stereographic"
  ;mapres@gsnPolar = "NH"
  ;mapres@mpMinLatF = 85.0
  ;mapres@mpProjection = "Satellite"

  trres = True
  trres@gsMarkerIndex = 1
  trres@gsMarkerColor = "black"

  vorores = True
  vorores@gsLineColor = "Red"

  do k = 100, 110;numTime
    map = gsn_csm_map(wks, mapres)

    lonIn = doubletofloat(f->q_lon(k,idx)*fac)
    latIn = doubletofloat(f->q_lat(k,idx)*fac)

    gsn_polymarker(wks, map, lonIn, latIn, trres)

    csvoro(latIn, lonIn, 0, 1, latOut, lonOut, circumradius, numActual, numVertex, idxVertex)
    do i = 0, numSparseTracer-1
      csvoro(latIn, lonIn, i, 0, latOut, lonOut, circumradius, numActual, numVertex, idxVertex)
      do j = 1, numVertex-1
        lon(0) = lonOut(idxVertex(j-1))
        lat(0) = latOut(idxVertex(j-1))
        lon(1) = lonOut(idxVertex(j))
        lat(1) = latOut(idxVertex(j))
        gsn_polyline(wks, map, lon, lat, vorores)
      end do
    end do

    frame(wks)
    system("echo time step "+k)
  end do

end

