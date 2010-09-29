load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    f1 = addfile("flow.nc", "r")
    f2 = addfile("flow_interp.nc", "r")

    numTimeStep = dimsizes(f1->time)

    wks = gsn_open_wks("x11", "TestMeshInterp")

    gsn_define_colormap(wks, "wgne15")

    resCommon = True
    resCommon@gsnDraw = False
    resCommon@gsnFrame = False

    ; ************************************************************************ ;
    ;                         VELOCITY FIGURE RESOUCES
    res1 = resCommon
    res1@vcLineArrowColor = "red"
    res1@vcRefMagnitudeF = 40.
    res1@vcRefLengthF = 0.02
    res1@vcMinDistanceF = 0.017

    res2 = res1
    res2@vcLineArrowColor = "black"

    ; ************************************************************************ ;
    ;                       DIVERGENCE FIGURE RESOUCES
    res3 = resCommon
    res3@cnFillOn = True
    res3@cnLinesOn = False
    res3@gsnSpreadColors = True
    res3@lbLabelAngleF = 0
    res3@lbOrientation = "vertical"

    ; ************************************************************************ ;
    ;                           PANEL RESOUCES
    resPanel = True
    resPanel@gsnMaximize = True

    do i = 0, numTimeStep-1
        ;velocityOrigin = gsn_csm_vector_map_ce(wks, f1->u(i,:,:), f1->v(i,:,:), res1)
        res1@gsnPolar = "NH"
        res1@mpMinLatF = 85.
        velocityOrigin = gsn_csm_vector_map_polar(wks, f1->u(i,:,:), f1->v(i,:,:), res1)
        velocityInterp = gsn_csm_vector(wks, f2->u(i,:,:), f2->v(i,:,:), res2)
        overlay(velocityOrigin, velocityInterp)
        ;divergenceInterp = gsn_csm_contour_map_ce(wks, f2->div(i,:,:), res3)
        res3@gsnPolar = "NH"
        res3@mpMinLatF = 85.
        divergenceInterp = gsn_csm_contour_map_polar(wks, f2->div(i,:,:), res3)
        gsn_panel(wks, (/velocityOrigin,divergenceInterp/), (/2,1/), resPanel)
    end do

end
