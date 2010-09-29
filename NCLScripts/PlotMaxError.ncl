load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  if (.not. isdefined("fn")) then
    fn = "./grid.nc"
  end if

  f = addfile(fn, "r")

  numTimeStep = dimsizes(f->time)
  errorDims = dimsizes(f->error(0,:,:))

  wks = gsn_open_wks("x11", "MaximumError")
  gsn_define_colormap(wks, "BlWhRe")

                                ; maximum error spatial plot
  res = True
  res@gsnFrame = False
  res@gsnMaximize = True
  res@gsnAddCyclic = True
  res@cnFillOn = True
  res@cnLinesOn = False
  res@gsnSpreadColors = True
  res@mpFillOn = False
  res@mpOutlineOn = True
  res@lbLabelAngleF = 45.0

                                ; resource for marker
  mkres = True
  mkres@gsMarkerIndex = 1
  mkres@gsMarkerColor = "green"
  mkres@gsMarkerSizeF = 0.05

  do i = 0, numTimeStep-1

    temp = ndtooned(abs(f->error(i,:,:)))

    indices = ind_resolve(maxind(temp), errorDims)

    if (f->lat(indices(0,0)) .gt. 0.0) then
      res@gsnPolar = "NH"
    else
      res@gsnPolar = "SH"
    end if

    plot = gsn_csm_contour_map_polar(wks, f->error(i,:,:), res)

    gsn_polymarker(wks, plot, f->lon(indices(:,1)), f->lat(indices(:,0)), mkres)

    frame(wks)

    system("echo time step "+i)

  end do

end
