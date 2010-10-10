load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

    f = addfile("flow.nc", "r")

    numTimeStep = dimsizes(f->time)

    wks = gsn_open_wks("x11", "Flow")

    res = True
    res@mpCenterLonF = 180.
    res@vcRefMagnitudeF = 1
    res@vcRefLengthF = 0.02
    res@vcMinDistanceF = 0.017

    do i = 0, numTimeStep-1
        plot = gsn_csm_vector_map_polar(wks, f->u(i,:,:), f->v(i,:,:), res)
    end do

end

