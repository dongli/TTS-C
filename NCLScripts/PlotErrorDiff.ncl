load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  if (.not. isdefined("figtype")) then
    figtype = "x11"
  end if

  if (.not. isdefined("title")) then
    title = ""
  end if

  f1 = addfile(fn1, "r")
  f2 = addfile(fn2, "r")

  numTimeStep = dimsizes(f1->time)
  
  diff = new((/3,numTimeStep/), double)

  diff(0,:) = f1->L1Error-f2->L1Error
  diff(1,:) = f1->L2Error-f2->L2Error
  diff(2,:) = f1->LInfError-f2->LInfError

  wks = gsn_open_wks(figtype, "ErrorDifference")

  res = True
  res@tiMainString = title
  res@tiXAxisString = "Time (Days)"
  res@tiYAxisString = "Normalized errors"
  res@trYMaxF = max(diff)+0.1
  res@trYMinF = min(diff)-0.001
  res@pmLegendDisplayMode = "Always"
  res@pmLegendSide = "Top"
  res@pmLegendParallelPosF = 0.7
  res@pmLegendOrthogonalPosF = -0.4
  res@pmLegendWidthF = 0.1
  res@pmLegendHeightF = 0.15
  res@lgPerimOn = False
  res@lgLabelFontHeightF = .02
  res@xyExplicitLegendLabels = (/"L1 diff.","L2 diff.","LInf diff."/)
  res@xyLineThicknesses = (/3.,3.,3./)

  plot = gsn_csm_xy(wks, f1->time/86400., diff, res)

end
