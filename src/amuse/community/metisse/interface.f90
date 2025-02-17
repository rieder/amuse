module metisseInterface
    use iso_c_binding
    use store_stars, only: stars
    use track_support
    use z_support
    implicit none
    type(stars) :: star_system
    real(c_double), allocatable :: mass_array(:)

    contains

    function initialize(error)
        implicit none
        integer :: error
        integer :: initialize

        real(c_double) :: zpars(20)

        initialize = -1

        ! Need to define this front end for METISSE
        call initialize_front_end("amuse")
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
  
  function delete_star(index_of_the_star)
    implicit none
    integer :: index_of_the_star
    integer :: delete_star
    delete_star=0
  end function
  
  function evolve_for(index_of_the_star, delta_t)
    implicit none
    integer :: index_of_the_star
    real(c_double) :: delta_t
    integer :: evolve_for
    evolve_for=0
  end function
  
  function evolve_one_step(index_of_the_star)
    implicit none
    integer :: index_of_the_star
    integer :: evolve_one_step
    evolve_one_step=0
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
  
  function new_particle(index_of_the_star, mass)
    implicit none
    integer, intent(inout) :: index_of_the_star
    real(c_double), intent(inout) :: mass
    integer :: new_particle
    index_of_the_star = star_system%new_star(mass)
    new_particle=0
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

