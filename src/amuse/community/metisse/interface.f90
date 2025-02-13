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
    double precision :: delta_t
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
    double precision :: age
    integer :: get_age
    get_age=0
  end function
  
  function get_luminosity(index_of_the_star, luminosity)
    implicit none
    integer :: index_of_the_star
    double precision :: luminosity
    integer :: get_luminosity
    get_luminosity=0
  end function
  
  function get_mass(index_of_the_star, mass)
    implicit none
    integer :: index_of_the_star
    double precision :: mass
    integer :: get_mass
    get_mass=0
  end function
  
  function get_metallicity(metallicity)
    implicit none
    double precision :: metallicity
    integer :: get_metallicity
    get_metallicity=0
  end function
  
  function get_number_of_particles(number_of_particles)
    implicit none
    integer :: number_of_particles
    integer :: get_number_of_particles
    get_number_of_particles=0
  end function
  
  function get_radius(index_of_the_star, radius)
    implicit none
    integer :: index_of_the_star
    double precision :: radius
    integer :: get_radius
    get_radius=0
  end function
  
  function get_stellar_type(index_of_the_star, stellar_type)
    implicit none
    integer :: index_of_the_star, stellar_type
    integer :: get_stellar_type
    get_stellar_type=0
  end function
  
  function get_temperature(index_of_the_star, temperature)
    implicit none
    integer :: index_of_the_star
    double precision :: temperature
    integer :: get_temperature
    get_temperature=0
  end function
  
  function get_time_step(index_of_the_star, time_step)
    implicit none
    integer :: index_of_the_star
    double precision :: time_step
    integer :: get_time_step
    get_time_step=0
  end function
  
  function initialize_code()
    implicit none
    integer :: initialize_code
    initialize_code=0
  end function
  
  function new_particle(index_of_the_star, mass)
    implicit none
    integer :: index_of_the_star
    double precision :: mass
    integer :: new_particle
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
    double precision :: metallicity
    integer :: set_metallicity
    set_metallicity=0
  end function
  
  
end module

