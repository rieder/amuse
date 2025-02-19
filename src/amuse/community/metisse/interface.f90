module metisseInterface
    use iso_c_binding
    use store_stars, only: stars
    use track_support
    use z_support
    implicit none
    type(stars) :: star_system
    real(c_double), allocatable :: mass_array(:)

    contains

    ! standard AMUSE interface functions:
    ! initialize
    ! commit_parameters
    ! commit_particles
    ! cleanup_code

    function initialize(error)
        implicit none
        integer :: error
        integer :: initialize

        real(c_double) :: zpars(20)

        initialize = -1

        ! Need to define this front end for METISSE
        call initialize_front_end("COSMIC")
        initial_Z = -1.0_c_double

        call METISSE_zcnsts(initial_Z,zpars,'','', error)
        if (error/=0) return

        write(*,*) "Number of tracks: ", number_of_tracks

        initialize = 0
    end function

    function cleanup_code()
        implicit none
        integer :: cleanup_code
        cleanup_code=0
    end function
    
    function commit_parameters()
        implicit none
        integer :: commit_parameters
        commit_parameters=0
    end function
    
    function commit_particles()
        implicit none
        integer :: commit_particles
        integer :: number_of_particles
        integer :: error
    
        error = get_number_of_particles(number_of_particles)
        allocate(mass_array(number_of_particles))
        mass_array = 0.0
    
        commit_particles=0
    end function

    ! setters / getters for tracks
    ! metallicity_dir (string)
    ! metallicity_dir_he (string)
    ! z_accuracy_limit (real)
    ! mass_accuracy_limit (real)

    function set_metallicity_dir(metallicity_dir_in)
        implicit none
        character(len=256) :: metallicity_dir_in
        integer :: set_metallicity_dir
        METALLICITY_DIR = metallicity_dir_in
        set_metallicity_dir = 0
    end function

    function get_metallicity_dir(metallicity_dir_out)
        implicit none
        character(len=256) :: metallicity_dir_out
        integer :: get_metallicity_dir
        metallicity_dir_out = METALLICITY_DIR
        get_metallicity_dir = 0
    end function

    function set_metallicity_dir_he(metallicity_dir_he_in)
        implicit none
        character(len=256) :: metallicity_dir_he_in
        integer :: set_metallicity_dir_he
        METALLICITY_DIR_HE = metallicity_dir_he_in
        set_metallicity_dir_he = 0
    end function

    function get_metallicity_dir_he(metallicity_dir_he_out)
        implicit none
        character(len=256) :: metallicity_dir_he_out
        integer :: get_metallicity_dir_he
        metallicity_dir_he_out = METALLICITY_DIR_HE
        get_metallicity_dir_he = 0
    end function

    function set_z_accuracy_limit(z_accuracy_limit_in)
        implicit none
        real(c_double) :: z_accuracy_limit_in
        integer :: set_z_accuracy_limit
        z_accuracy_limit = z_accuracy_limit_in
        set_z_accuracy_limit = 0
    end function

    function get_z_accuracy_limit(z_accuracy_limit_out)
        implicit none
        real(c_double) :: z_accuracy_limit_out
        integer :: get_z_accuracy_limit
        z_accuracy_limit_out = z_accuracy_limit
        get_z_accuracy_limit = 0
    end function

    function set_mass_accuracy_limit(mass_accuracy_limit_in)
        implicit none
        real(c_double) :: mass_accuracy_limit_in
        integer :: set_mass_accuracy_limit
        mass_accuracy_limit = mass_accuracy_limit_in
        set_mass_accuracy_limit = 0
    end function

    function get_mass_accuracy_limit(mass_accuracy_limit_out)
        implicit none
        real(c_double) :: mass_accuracy_limit_out
        integer :: get_mass_accuracy_limit
        mass_accuracy_limit_out = mass_accuracy_limit
        get_mass_accuracy_limit = 0
    end function

    ! setters / getters for misc controls
    ! verbose (bool)
    ! construct_postagb_track (bool)

    function set_verbose(verbose_in)
        implicit none
        logical :: verbose_in
        integer :: set_verbose
        verbose = verbose_in
        set_verbose = 0
    end function

    function get_verbose(verbose_out)
        implicit none
        logical :: verbose_out
        integer :: get_verbose
        verbose_out = verbose
        get_verbose = 0
    end function

    function set_construct_postagb_track(construct_postagb_track_in)
        implicit none
        logical :: construct_postagb_track_in
        integer :: set_construct_postagb_track
        construct_postagb_track = construct_postagb_track_in
        set_construct_postagb_track = 0
    end function

    function get_construct_postagb_track(construct_postagb_track_out)
        implicit none
        logical :: construct_postagb_track_out
        integer :: get_construct_postagb_track
        construct_postagb_track_out = construct_postagb_track
        get_construct_postagb_track = 0
    end function

    ! setters/getters for parameters
    ! initial_metallicity(real)
    ! wd_mass_scheme (string, 256)
    ! use_initial_final_mass_relation(bool)
    ! bhns_mass_scheme (string, 256)
    ! max_ns_mass (real)
    ! allow_electron_capture (bool)

    function set_initial_metallicity(initial_metallicity_in)
        implicit none
        real(c_double) :: initial_metallicity_in
        integer :: set_initial_metallicity
        initial_Z = initial_metallicity_in
        set_initial_metallicity = 0
    end function

    function get_initial_metallicity(initial_metallicity_out)
        implicit none
        real(c_double) :: initial_metallicity_out
        integer :: get_initial_metallicity
        initial_metallicity_out = initial_Z
        get_initial_metallicity = 0
    end function

    function set_wd_mass_scheme(wd_mass_scheme_in)
        implicit none
        character(len=256) :: wd_mass_scheme_in
        integer :: set_wd_mass_scheme
        WD_mass_scheme = wd_mass_scheme_in
        set_wd_mass_scheme = 0
    end function

    function get_wd_mass_scheme(wd_mass_scheme_out)
        implicit none
        character(len=256) :: wd_mass_scheme_out
        integer :: get_wd_mass_scheme
        wd_mass_scheme_out = WD_mass_scheme
        get_wd_mass_scheme = 0
    end function

    function set_use_initial_final_mass_relation(use_initial_final_mass_relation_in)
        implicit none
        logical :: use_initial_final_mass_relation_in
        integer :: set_use_initial_final_mass_relation
        use_initial_final_mass_relation = use_initial_final_mass_relation_in
        set_use_initial_final_mass_relation = 0
    end function

    function get_use_initial_final_mass_relation(use_initial_final_mass_relation_out)
        implicit none
        logical :: use_initial_final_mass_relation_out
        integer :: get_use_initial_final_mass_relation
        use_initial_final_mass_relation_out = use_initial_final_mass_relation
        get_use_initial_final_mass_relation = 0
    end function

    function set_bhns_mass_scheme(bhns_mass_scheme_in)
        implicit none
        character(len=256) :: bhns_mass_scheme_in
        integer :: set_bhns_mass_scheme
        BHNS_mass_scheme = bhns_mass_scheme_in
        set_bhns_mass_scheme = 0
    end function

    function get_bhns_mass_scheme(bhns_mass_scheme_out)
        implicit none
        character(len=256) :: bhns_mass_scheme_out
        integer :: get_bhns_mass_scheme
        bhns_mass_scheme_out = BHNS_mass_scheme
        get_bhns_mass_scheme = 0
    end function

    function set_max_ns_mass(max_ns_mass_in)
        implicit none
        real(c_double) :: max_ns_mass_in
        integer :: set_max_ns_mass
        max_NS_mass = max_ns_mass_in
        set_max_ns_mass = 0
    end function

    function get_max_ns_mass(max_ns_mass_out)
        implicit none
        real(c_double) :: max_ns_mass_out
        integer :: get_max_ns_mass
        max_ns_mass_out = max_NS_mass
        get_max_ns_mass = 0
    end function

    function set_allow_electron_capture(allow_electron_capture_in)
        implicit none
        logical :: allow_electron_capture_in
        integer :: set_allow_electron_capture
        allow_electron_capture = allow_electron_capture_in
        set_allow_electron_capture = 0
    end function

    function get_allow_electron_capture(allow_electron_capture_out)
        implicit none
        logical :: allow_electron_capture_out
        integer :: get_allow_electron_capture
        allow_electron_capture_out = allow_electron_capture
        get_allow_electron_capture = 0
    end function

    ! setters/getters for timestep control
    ! pts_1 to pts_3 (real)

    function set_time_step_pts_1(pts_1_in)
        implicit none
        real(c_double) :: pts_1_in
        integer :: set_time_step_pts_1
        pts_1 = pts_1_in
        set_time_step_pts_1 = 0
    end function

    function get_time_step_pts_1(pts_1_out)
        implicit none
        real(c_double) :: pts_1_out
        integer :: get_time_step_pts_1
        pts_1_out = pts_1
        get_time_step_pts_1 = 0
    end function

    function set_time_step_pts_2(pts_2_in)
        implicit none
        real(c_double) :: pts_2_in
        integer :: set_time_step_pts_2
        pts_2 = pts_2_in
        set_time_step_pts_2 = 0
    end function

    function get_time_step_pts_2(pts_2_out)
        implicit none
        real(c_double) :: pts_2_out
        integer :: get_time_step_pts_2
        pts_2_out = pts_2
        get_time_step_pts_2 = 0
    end function

    function set_time_step_pts_3(pts_3_in)
        implicit none
        real(c_double) :: pts_3_in
        integer :: set_time_step_pts_3
        pts_3 = pts_3_in
        set_time_step_pts_3 = 0
    end function

    function get_time_step_pts_3(pts_3_out)
        implicit none
        real(c_double) :: pts_3_out
        integer :: get_time_step_pts_3
        pts_3_out = pts_3
        get_time_step_pts_3 = 0
    end function

    ! particle management:
    ! new_particle, delete_particle

    function new_particle(index_of_the_particle, mass)
        implicit none
        integer, intent(inout) :: index_of_the_particle
        real(c_double), intent(inout) :: mass
        integer :: new_particle
        index_of_the_particle = star_system%new_star(mass)
    end function

    function delete_star(index_of_the_star)
        implicit none
        integer :: index_of_the_star
        integer :: delete_star
        call star_system%remove_star(index_of_the_star)
        delete_star = 0
    end function

    ! evolving stars:
    ! evolve_for, evolve_one_step
    ! evolve_model, evolve_stars
  
    function evolve_for(index_of_the_star, delta_t)
        implicit none
        integer :: index_of_the_star
        real(c_double) :: delta_t
        integer :: evolve_for
        evolve_for = 0
    end function
    
    function evolve_one_step(index_of_the_star)
        implicit none
        integer :: index_of_the_star
        integer :: evolve_one_step
        integer :: error
        real(c_double) :: time_step
        real(c_double) :: mass
        real(c_double) :: age
   
        write(*,*) 'evolve_one_step', index_of_the_star
        call star_system%get_time_step(index_of_the_star, time_step, error)
        call star_system%get_initial_mass(index_of_the_star, mass, error)
        call star_system%get_age(index_of_the_star, age, error)
        write(*,*) 'age, mass, time_step', age, mass, time_step
        call allocate_track(1, mass)
        write(*,*) 'allocate_track done'
        call evolv_metisse(mass, age + time_step, error, 1)
        write(*,*) 'evolv_metisse done'
        call dealloc_track()
        write(*,*) 'dealloc_track done'
        if (error /= 0) then
            call star_system%set_mass(index_of_the_star, mass, error)
            call star_system%set_age(index_of_the_star, age + time_step, error)
        end if
        write(*,*) 'evolve_one_step done'
        evolve_one_step = 0
    end function

    function evolve_model(t_end)
        implicit none
        real(c_double) :: t_end
        integer :: evolve_model
        integer :: number_of_stars
        integer :: ierr
        integer :: i
        real(c_double) :: mass

        call star_system%get_number_of_stars(number_of_stars)

        do i = 1, number_of_stars
            call star_system%get_initial_mass(i, mass, ierr)
            call allocate_track(1, mass)
            call evolv_metisse(mass, t_end, ierr, 1)
            call dealloc_track()
            if (ierr /= 0) then
                call star_system%set_mass(i, mass, ierr)
                call star_system%set_age(i, t_end, ierr)
            end if
        end do
        evolve_model = 0
    end function

    function evolve_stars(indices_of_stars, delta_time, error, n)
        implicit none
        integer, intent(in) :: n
        integer :: error
        integer, intent(in), dimension(n) :: indices_of_stars
        real(c_double) :: delta_time
        integer :: evolve_stars
        integer :: i
        integer :: number_of_stars_to_evolve

        number_of_stars_to_evolve = size(indices_of_stars)

        do i = 1, number_of_stars_to_evolve
            error = evolve_for(indices_of_stars(i), delta_time)
            if (error /= 0) then
                evolve_stars = error
                return
            end if
        end do

        evolve_stars = 0
    end function
  
  function get_age(index_of_the_star, age)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: age
    integer :: get_age
    call star_system%get_age(index_of_the_star, age, get_age)
  end function
  
  function get_luminosity(index_of_the_star, luminosity)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: luminosity
    integer :: get_luminosity
    call star_system%get_luminosity(index_of_the_star, luminosity, get_luminosity)
  end function
  
  function get_mass(index_of_the_star, mass)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: mass
    integer :: get_mass
    call star_system%get_mass(index_of_the_star, mass, get_mass)
  end function

  function get_metallicity(metallicity)
    implicit none
    real(c_double) :: metallicity
    integer :: get_metallicity
    ! what to do here depends on whether metallicity can be set for individual stars or only globally
    get_metallicity=0
  end function

  function get_epoch(index_of_the_star, epoch)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: epoch
    integer :: get_epoch
    call star_system%get_epoch(index_of_the_star, epoch, get_epoch)
  end function

  function get_core_mass(index_of_the_star, core_mass)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: core_mass
    integer :: get_core_mass
    call star_system%get_core_mass(index_of_the_star, core_mass, get_core_mass)
  end function

  function get_core_radius(index_of_the_star, core_radius)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: core_radius
    integer :: get_core_radius
    call star_system%get_core_radius(index_of_the_star, core_radius, get_core_radius)
  end function

  function get_convective_envelope_mass(index_of_the_star, convective_envelope_mass)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: convective_envelope_mass
    integer :: get_convective_envelope_mass
    call star_system%get_convective_envelope_mass(index_of_the_star, convective_envelope_mass, get_convective_envelope_mass)
  end function

  function get_convective_envelope_radius(index_of_the_star, convective_envelope_radius)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: convective_envelope_radius
    integer :: get_convective_envelope_radius
    call star_system%get_convective_envelope_radius(index_of_the_star, convective_envelope_radius, get_convective_envelope_radius)
  end function

  function get_CO_core_mass(index_of_the_star, CO_core_mass)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: CO_core_mass
    integer :: get_CO_core_mass
    call star_system%get_CO_core_mass(index_of_the_star, CO_core_mass, get_CO_core_mass)
  end function

  function get_main_sequence_lifetime(index_of_the_star, main_sequence_lifetime)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: main_sequence_lifetime
    integer :: get_main_sequence_lifetime
    call star_system%get_main_sequence_lifetime(index_of_the_star, main_sequence_lifetime, get_main_sequence_lifetime)
  end function

  function get_number_of_particles(number_of_particles)
    implicit none
    integer :: number_of_particles
    integer :: get_number_of_particles
    call star_system%get_number_of_stars(number_of_particles)
    get_number_of_particles = 0
  end function
  
  function get_radius(index_of_the_star, radius)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: radius
    integer :: get_radius
    call star_system%get_radius(index_of_the_star, radius, get_radius)
  end function
  
  function get_stellar_type(index_of_the_star, stellar_type)
    implicit none
    integer :: index_of_the_star, stellar_type
    integer :: get_stellar_type
    call star_system%get_stellar_type(index_of_the_star, stellar_type, get_stellar_type)
  end function
  
  function get_temperature(index_of_the_star, temperature)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: temperature
    integer :: get_temperature
    call star_system%get_temperature(index_of_the_star, temperature, get_temperature)
  end function
  
  function get_time_step(index_of_the_star, time_step)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: time_step
    integer :: get_time_step
    call star_system%get_time_step(index_of_the_star, time_step, get_time_step)
  end function

  function get_initial_mass(index_of_the_star, mass)
    implicit none
    integer :: index_of_the_star
    double precision :: mass
    integer :: get_initial_mass
    call star_system%get_initial_mass(index_of_the_star, mass, get_initial_mass)
  end function
  
  function initialize_code()
    implicit none
    integer :: initialize_code
    initialize_code=0
  end function

  
  function recommit_parameters()
    implicit none
    integer :: recommit_parameters
    recommit_parameters=0
  end function
  
  function recommit_particles()
    implicit none
    integer :: recommit_particles
    recommit_particles=0
  end function
  
  function set_metallicity(metallicity)
    implicit none
    real(c_double) :: metallicity
    integer :: set_metallicity
    set_metallicity=0
  end function
  
  
end module

