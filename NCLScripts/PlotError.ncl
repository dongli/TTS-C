load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  if (.not. isdefined("fn")) then
    fn = "./grid.nc"
  end if

  if (.not. isdefined("title")) then
    title = ""
  end if

  if (.not. isdefined("figtype")) then
    figtype = "x11"
  end if

  f = addfile(fn, "r")
  
  numTimeStep = dimsizes(f->time)
  
  error = new((/3,numTimeStep/), double)
  error(0,:) = f->L1Error
  error(1,:) = f->L2Error
  error(2,:) = f->LInfError

  wks = gsn_open_wks(figtype, "norm_error")

                                ; error-time plot
  res = True
  res@gsnMaximize = True
  res@tiMainString = title
  res@tiXAxisString = "Time (Days)"
  res@tiYAxisString = "Normalized errors"
  res@trYMaxF = 0.1
  res@trYMinF = 0.0
  res@pmLegendDisplayMode = "Always"
  res@pmLegendSide = "Top"
  res@pmLegendParallelPosF = 0.2
  res@pmLegendOrthogonalPosF = -0.4
  res@pmLegendWidthF = 0.2
  res@pmLegendHeightF = 0.15
  res@lgPerimOn = False
  res@lgLabelFontHeightF = .02
  res@xyExplicitLegendLabels = (/"L1","L2","LInf"/)

  res@xyLineThicknesses = (/3.,3.,3./)

  plot = gsn_csm_xy(wks, f->time/86400., error, res)

end
