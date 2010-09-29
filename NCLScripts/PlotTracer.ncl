load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  if (.not. isdefined("fn")) then
    fn = "tracer.nc"
  end if

  if (.not. isdefined("figtype")) then
    figtype = "x11"
  end if

  if (.not. isdefined("outfile")) then
    outfile = "TracerTrajectory"
  end if

  if (.not. isdefined("plottype")) then
    plottype = "NH"
  end if

  PI = 4.*atan(1.)
  fac = 180./PI
  
  f = addfile(fn, "r")
  
  time = f->time
  
  numTimeStep = dimsizes(time)
  
  numTracer = dimsizes(f->q(0,:))

  idx1 = ispan(0, numTracer/2, 10)
  ;idx2 = ispan(numTracer/2+1, numTracer-1, 1)

  wks = gsn_open_wks(figtype, outfile)
  
                                ; resource for map
  mpres = True
  mpres@gsnFrame = False
  mpres@gsnMaximize = True
  mpres@mpGridLatSpacingF = 1.5

  if (plottype .eq. "SH") then
    mpres@gsnPolar = "SH"
    mpres@mpMaxLatF = -80.0
  end if
  if (plottype .eq. "NH") then
    mpres@gsnPolar = "NH"
    mpres@mpMinLatF = 80.0
  end if

                                ; resource for tracer
  trres1 = True
  trres1@gsMarkerIndex = 1
  trres1@gsMarkerSizeF = 0.01
  trres1@gsMarkerColor = "black"

  trres2 = True
  trres2@gsMarkerIndex = 1
  trres2@gsMarkerSizeF = 0.01
  trres2@gsMarkerColor = "red"

  do i = 0, numTimeStep-1
    if (plottype .eq. "NH" .or. plottype .eq. "SH") then
      map = gsn_csm_map_polar(wks, mpres)
    end if
    if (plottype .eq. "ce") then
      map = gsn_csm_map_ce(wks, mpres)
    end if
  
    gsn_polymarker(wks, map, f->q_lon(i,idx1)*fac, f->q_lat(i,idx1)*fac, trres1)
    ;gsn_polymarker(wks, map, f->q_lon(i,idx2)*fac, f->q_lat(i,idx2)*fac, trres2)
  
    frame(wks)

    system("echo time step "+i)
  end do
  
  ;do j = 0, 9
  ;  system("echo "+j)
  ;  map = gsn_csm_map_polar(wks, mpres)
  ;  gsn_polymarker(wks, map, f->q_lon(2,0:j)*fac, f->q_lat(2,0:j)*fac, trres1)
  ;  frame(wks)
  ;end do

end
