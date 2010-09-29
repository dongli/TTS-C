load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  if (.not. isdefined("fn")) then
    fn = "../grid.nc"
  end if

  if (.not. isdefined("title")) then
    title = ""
  end if

  if (.not. isdefined("figtype")) then
    figtype = "x11"
  end if

  if (.not. isdefined("outfile")) then
    outfile = "TotalMass"
  end if

  f = addfile(fn, "r")

  wks = gsn_open_wks(figtype, outfile)

  res = True
  res@gsnMaximize = True
  res@tiMainString = title
  res@tiXAxisString = "Time (Days)"
  res@tiYAxisString = "TotalMass"
  
  plot = gsn_csm_xy(wks, f->time/86400., f->mass, res)

end
