module FloatingPoint

    use MsgManager

    implicit none

    real(8), parameter :: eps = 1.0D-10

contains

    logical function FloatingPoint_Check(number, name) result(res)
        real(8), intent(in) :: number
        character(*), intent(in) :: name

#if (defined FC_IFORT)
        integer fp_status

        fp_status = fp_class(number)

        res = .true.
        if (fp_status == 0) then
            call MsgManager_Speak(Error, &
                "Encounter signaling NaN with "//trim(name)//".")
            res = .false.
        end if
        if (fp_status == 1) then
            call MsgManager_Speak(Error, &
                "Encounter quiet NaN with "//trim(name)//".")
            res = .false.
        end if
#else
        res = .true.
        if (isnan(number)) then
            call MsgManager_Speak(Error, &
                "Encounter NaN with "//trim(name)//".")
            res = .false.
        end if
#endif

    end function FloatingPoint_Check

end module FloatingPoint
