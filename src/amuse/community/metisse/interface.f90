module metisseInterface
    use iso_c_binding
    implicit none
    contains

    function initialize(error)
        implicit none
        integer :: error
        integer :: initialize
        initialize = 0
    end function

    function teststar(mass_in, time, mass_out, error)
        use track_support
        use z_support
        implicit none
        real(dp) :: mass_in, mass_out, time
        integer :: error
        integer :: teststar

        real(dp) :: zpars(20)

        call initialize_front_end('main')
        !call METISSE_zcnsts(initial_Z,zpars,'','',error)
        teststar = 0
    end function

end module

