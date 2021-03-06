module SphereService

    use MsgManager
    use RunManager
    use FloatingPoint

    implicit none

    real(8), parameter :: Equator = 0.0d0
    real(8), parameter :: PI    = 4.0d0*atan(1.0d0)
    real(8), parameter :: PI2   = 2.0d0*PI
    real(8), parameter :: PI05  = 0.5d0*PI
    real(8), parameter :: PI025 = 0.25d0*PI
    real(8), parameter :: PI15  = 1.5d0*PI
    real(8), parameter :: Re    = 6371.229d3
    real(8), parameter :: Re2   = Re**2.0d0
    real(8), parameter :: Rad2Deg = 180.0d0/PI

contains

    subroutine RotationTransform(lonP, latP, lonO, latO, lonR, latR)
        real(8), intent(in) :: lonP, latP ! rotated pole coordinate
        real(8), intent(in) :: lonO, latO ! original coordinate
        real(8), intent(out), optional :: lonR, latR ! rotated coordinate

        real(8) cosLatO, sinLatO, cosLatP, sinLatP, sinDlon, cosDlon
        real(8) temp1, temp2, temp3, dlon

        call MsgManager_RecordSpeaker("RotationTransform")

        dlon = lonO-lonP
        cosLatO = cos(latO)
        sinLatO = sin(latO)
        cosLatP = cos(latP)
        sinLatP = sin(latP)
        sinDlon = sin(dlon)
        cosDlon = cos(dlon)
        if (present(lonR)) then
            temp1 = cosLatO*sinDlon
            temp2 = cosLatO*sinLatP*cosDlon-cosLatP*sinLatO
            lonR = atan2(temp1, temp2)
            if (lonR < 0.0d0) lonR = PI2+lonR
        end if
        if (present(latR)) then
            temp1 = sinLatO*sinLatP
            temp2 = cosLatO*cosLatP*cosDlon
            temp3 = temp1+temp2
            temp3 = min(max(temp3, -1.0d0), 1.0d0)
            latR = asin(temp3)
        end if

        call MsgManager_DeleteSpeaker

    end subroutine RotationTransform

    subroutine InverseRotationTransform(lonP, latP, lonO, latO, lonR, latR)
        real(8), intent(in) :: lonP, latP  ! rotated pole coordinate
        real(8), intent(out) :: lonO, latO ! original coordinate
        real(8), intent(in) :: lonR, latR  ! rotated coordinate

        real(8) temp1, temp2, temp3

        call MsgManager_RecordSpeaker("InverseRotationTransform")

        temp1 = cos(latR)*sin(lonR)
        temp2 = sin(latR)*cos(latP)+cos(latR)*cos(lonR)*sin(latP)
        if (abs(temp2) < eps) temp2 = 0.0d0 ! DONG Li: This trick is due to the inaccuracy of trigonometry calculation.
        lonO = atan2(temp1, temp2)
        lonO = lonP+lonO
#if (defined DEBUG)
        if (.not. FloatingPoint_Check(lonO, "lonO")) then
            call RunManager_EndRun
        end if
#endif
        if (lonO > PI2) lonO = lonO-PI2
        temp1 = sin(latR)*sin(latP)
        temp2 = cos(latR)*cos(latP)*cos(lonR)
        temp3 = temp1-temp2
        temp3 = min(max(temp3, -1.0d0), 1.0d0)
        latO = asin(temp3)
#if (defined DEBUG)
        if (.not. FloatingPoint_Check(latO, "latO")) then
            call RunManager_EndRun
        end if
#endif

        call MsgManager_DeleteSpeaker

    end subroutine InverseRotationTransform

end module SphereService
