load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; Note:
;   fn1 - The exact solution data file
;   fn2 - The simulated data file

begin

  ; set default value of arguments
  if (.not. isdefined("fn1")) then
    fn1 = "exact.nc"
  end if

  if (.not. isdefined("fn2")) then
    fn2 = "grid.nc"
  end if

  if (.not. isdefined("figtype")) then
    figtype = "x11"
  end if

  ; read in data file
  f1 = addfile(fn1, "r")
  f2 = addfile(fn2, "r")

  time = f1->time
  
  timeDims = dimsizes(time)
  numTimeStep = timeDims(0)

  ; configure graphics
  wks = gsn_open_wks(figtype, "animate")

  gsn_define_colormap(wks, "WhBlReWh")

  res = True
  res@gsnMaximize	= True
  res@gsnAddCyclic	= True
  res@mpFillOn		= False
  res@mpOutlineOn	= True
  res@mpPerimOn		= False
  res@mpProjection	= "Orthographic"
  res@mpGridAndLimbOn	= True
  res@cnFillOn		= True
  res@cnLinesOn		= False
  res@gsnSpreadColors	= True
  res@cnLevelSpacingF	= 0.05
  res@lbLabelStride	= 4
  res@lbOrientation	= "Vertical"
  res@gsnStringFontHeightF = 0.015

  txres = True
  txres@txFontQuality 	= "High"
  txres@txFontHeightF 	= 0.015
  txres@txJust 		= "CenterLeft"

  do i = 0, numTimeStep-1, 1
    gsn_text_ndc(wks, "time step: 1800s", 0.5, 0.94, txres)
    gsn_text_ndc(wks, "grid number: 240x120", 0.5, 0.92, txres)
    gsn_text_ndc(wks, "tracer number: 115,200", 0.5, 0.90, txres)
    
    res@mpCenterLonF = f1->vtxlon(i)
    res@mpCenterLatF = f1->vtxlat(i)
    
    if (i .eq. 0) then
      temp1 = f1->vtxlat(0)
      temp2 = f1->vtxlat(1)
    else
      temp1 = f1->vtxlat(i-1)
      temp2 = f1->vtxlat(i)
    end if
    if (temp1 .lt. temp2 .and. f1->vtxlat(i) .lt. 0.) then
      res@mpCenterRotF = 180.
    end if
    if (temp1 .gt. temp2 .and. f1->vtxlat(i) .gt. 0. .and. i .ne. 0) then
      res@mpCenterRotF = 360.
    end if

    plot = gsn_csm_contour_map(wks, f2->q(i,:,:), res)
    
    system("echo time step "+i)
  end do
  
end
