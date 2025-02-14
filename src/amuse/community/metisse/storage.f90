module store_stars
  use iso_c_binding
  implicit none

  type, public :: star
    private
    integer :: id
    real(c_double) :: mass
    real(c_double) :: age
    real(c_double) :: initial_mass
    real(c_double) :: time_step
    real(c_double) :: luminosity
    real(c_double) :: temperature
    real(c_double) :: metallicity
    real(c_double) :: radius
    integer :: stellar_type
  end type star

  type, public :: stars
    private
    type(star), allocatable :: star_array(:)
    integer :: num_stars = 0 ! number of stars in the system
    integer :: next_star_id = 1 ! the id of the next star, should only ever increase
  contains
    procedure, public :: new_star
    procedure, public :: remove_star
    procedure, private :: resize
    procedure, private :: lookup_star_id
    procedure, private :: get_property_double
    procedure, private :: get_property_int
    procedure, public :: get_mass
    procedure, public :: get_radius
    procedure, public :: get_age
    procedure, public :: get_time_step
    procedure, public :: get_temperature
    procedure, public :: get_luminosity
    procedure, public :: get_stellar_type
    procedure, public :: get_metallicity
    procedure, public :: get_number_of_stars
    
  end type stars

contains

  subroutine get_number_of_stars(self, number_of_stars)
    class(stars), intent(in) :: self
    integer, intent(out) :: number_of_stars
    number_of_stars = self%num_stars
  end subroutine

  function new_star(self, initial_mass) result(new_id)
    class(stars), intent(inout) :: self
    real(c_double), intent(in) :: initial_mass
    integer :: new_id
    integer :: i

    self%num_stars = self%num_stars + 1
    call self%resize(self%num_stars)
    i = self%num_stars
    new_id = self%next_star_id

    self%star_array(i)%id = new_id
    self%star_array(i)%mass = initial_mass
    self%star_array(i)%age = 0.0_c_double
    self%star_array(i)%initial_mass = initial_mass
    self%star_array(i)%time_step = 0.0_c_double
    self%star_array(i)%luminosity = 0.0_c_double
    self%star_array(i)%temperature = 0.0_c_double
    self%star_array(i)%metallicity = 0.0_c_double
    self%star_array(i)%radius = 0.0_c_double
    self%star_array(i)%stellar_type = 0

    self%next_star_id = new_id + 1

    write(*,*) "new star with id: ", new_id
    write(*,*) "mass: ", initial_mass

  end function new_star

  subroutine remove_star(self, id)
    class(stars), intent(inout) :: self
    integer, intent(in) :: id

    integer :: i

    write(*,*) "removing star with id: ", id
    do i = 1, self%num_stars
      write(*,*) "i: ", i, " id: ", self%star_array(i)%id
      if (self%star_array(i)%id == id) then
        write(*,*) "removing star with id: ", id
        if (i /= self%num_stars) then
          self%star_array(i:self%num_stars-1) = self%star_array(i+1:self%num_stars)
        end if
        self%num_stars = self%num_stars - 1
        write(*,*) "num_stars: ", self%num_stars
        call self%resize(self%num_stars)
        exit
      end if
    end do

  end subroutine remove_star

  subroutine resize(self, new_size)
    class(stars), intent(inout) :: self
    integer, intent(in) :: new_size
    type(star), allocatable :: temp(:)
    integer :: current_size, new_capacity
  
    if (new_size <= 0) then
      self%num_stars = 0
      if (allocated(self%star_array)) deallocate(self%star_array)
      return
    end if
  
    current_size = self%num_stars
    new_capacity = max(100, int((1.1 ** ceiling(log10(real(new_size)))) * 10))
  
    if (.not. allocated(self%star_array) .or. new_capacity > size(self%star_array)) then
      if (allocated(self%star_array)) then
        allocate(temp(current_size))
        temp = self%star_array
        deallocate(self%star_array)
      end if
      allocate(self%star_array(new_capacity))
      if (allocated(temp)) then
        self%star_array(1:current_size) = temp
        deallocate(temp)
      end if
    end if
  
  end subroutine resize

  function lookup_star_id(self, id) result(index_of_the_star)
    class(stars), intent(in) :: self
    integer, intent(in) :: id
    integer :: index_of_the_star
    integer :: i

    do i = 1, self%num_stars
      if (self%star_array(i)%id == id) then
        index_of_the_star = i
        return
      end if
    end do

    index_of_the_star = 0
  end function

  ! Getters for all the stellar properties
  subroutine get_property_double(self, id, property_name, value, error)
    class(stars), intent(in) :: self
    integer, intent(in) :: id
    character(len=*), intent(in) :: property_name
    real(c_double), intent(out) :: value
    integer :: i, error
  
    i = lookup_star_id(self, id)
    if (i == 0) then
      value = 0.0_c_double
      error = -1  ! star not found
      return
    end if
  
    select case (trim(property_name))
      case ('mass')
        value = self%star_array(i)%mass
      case ('age')
        value = self%star_array(i)%age
      case ('luminosity')
        value = self%star_array(i)%luminosity
      case ('temperature')
        value = self%star_array(i)%temperature
      case ('time_step')
        value = self%star_array(i)%time_step
      case ('metallicity')
        value = self%star_array(i)%metallicity
      case ('radius')
        value = self%star_array(i)%radius
      case default
        value = 0.0_c_double
        error = -2  ! property not found
    end select
    error = 0
    write(*,*) "get_property_double: ", id, property_name, value
  end subroutine get_property_double

  subroutine get_property_int(self, id, property_name, value, error)
    class(stars), intent(in) :: self
    integer, intent(in) :: id
    character(len=*), intent(in) :: property_name
    integer, intent(out) :: value
    integer :: i, error
  
    i = lookup_star_id(self, id)
    if (i == 0) then
      value = 0.0_c_double
      error = -1  ! star not found
      return
    end if
  
    select case (trim(property_name))
      case ('stellar_type')
        value = self%star_array(i)%stellar_type
      case default
        value = 0
        error = -2  ! property not found
        return
    end select
    error = 0
    write(*,*) "get_property_int: ", id, property_name, value
  end subroutine get_property_int

  ! Setters for all the stellar properties
  subroutine set_property_double(self, id, property_name, value, error)
    class(stars), intent(inout) :: self
    integer, intent(in) :: id
    character(len=*), intent(in) :: property_name
    real(c_double), intent(in) :: value
    integer :: i, error
  
    i = lookup_star_id(self, id)
    if (i == 0) then
      error = -1  ! star not found
      return
    end if
  
    select case (trim(property_name))
      case ('mass')
        self%star_array(i)%mass = value
      case ('age')
        error = -3  ! not settable
        return
      case ('luminosity')
        error = -3  ! not settable
        return
      case ('temperature')
        error = -3  ! not settable
        return
      case ('time_step')
        self%star_array(i)%time_step = value
      case ('metallicity')
        if (self%star_array(i)%age > 0.0_c_double) then
          error = -4  ! not settable after having evolved
          return
        else
          self%star_array(i)%metallicity = value
        end if
      case ('radius')
        error = -3  ! not settable
        return
      case default
        error = -2  ! property not found
        return
    end select
    error = 0
    write(*,*) "set_property_double: ", id, property_name, value
  end subroutine

  subroutine set_property_int(self, id, property_name, value, error)
    class(stars), intent(inout) :: self
    integer, intent(in) :: id
    character(len=*), intent(in) :: property_name
    integer, intent(in) :: value
    integer :: i, error
  
    i = lookup_star_id(self, id)
    if (i == 0) then
      error = -1  ! star not found
      return
    end if
  
    select case (trim(property_name))
      case ('stellar_type')
        error = -3  ! not settable
        return
      case default
        error = -2  ! property not found
        return
    end select
    error = 0
    write(*,*) "set_property_int: ", id, property_name, value
  end subroutine


  ! getters for all the stellar properties
  subroutine get_mass(self, id, mass, error)
    class(stars), intent(in) :: self
    integer, intent(in) :: id
    real(c_double), intent(out) :: mass
    integer :: error
    call get_property_double(self, id, 'mass', mass, error)
  end subroutine

  subroutine get_radius(self, id, radius, error)
    class(stars), intent(in) :: self
    integer, intent(in) :: id
    real(c_double), intent(out) :: radius
    integer :: error
    call get_property_double(self, id, 'radius', radius, error)
  end subroutine

  subroutine get_age(self, id, age, error)
    class(stars), intent(in) :: self
    integer, intent(in) :: id
    real(c_double), intent(out) :: age
    integer :: error
    call get_property_double(self, id, 'age', age, error)
  end subroutine

  subroutine get_luminosity(self, id, luminosity, error)
    class(stars), intent(in) :: self
    integer, intent(in) :: id
    real(c_double), intent(out) :: luminosity
    integer :: error
    call get_property_double(self, id, 'luminosity', luminosity, error)
  end subroutine

  subroutine get_temperature(self, id, temperature, error)
    class(stars), intent(in) :: self
    integer, intent(in) :: id
    real(c_double), intent(out) :: temperature
    integer :: error
    call get_property_double(self, id, 'temperature', temperature, error)
  end subroutine

  subroutine get_time_step(self, id, time_step, error)
    class(stars), intent(in) :: self
    integer, intent(in) :: id
    real(c_double), intent(out) :: time_step
    integer :: error
    call get_property_double(self, id, 'time_step', time_step, error)
  end subroutine

  subroutine get_metallicity(self, id, metallicity, error)
    class(stars), intent(in) :: self
    integer, intent(in) :: id
    real(c_double), intent(out) :: metallicity
    integer :: error
    call get_property_double(self, id, 'metallicity', metallicity, error)
  end subroutine

  subroutine get_stellar_type(self, id, stellar_type, error)
    class(stars), intent(in) :: self
    integer, intent(in) :: id
    integer, intent(out) :: stellar_type
    integer :: error
    call get_property_int(self, id, 'stellar_type', stellar_type, error)
  end subroutine

  ! setters for all the stellar properties that are settable
  subroutine set_mass(self, id, mass, error)
    class(stars), intent(inout) :: self
    integer, intent(in) :: id
    real(c_double), intent(in) :: mass
    integer :: error
    call set_property_double(self, id, 'mass', mass, error)
  end subroutine

  subroutine set_metallicity(self, id, metallicity, error)
    class(stars), intent(inout) :: self
    integer, intent(in) :: id
    real(c_double), intent(in) :: metallicity
    integer :: error
    call set_property_double(self, id, 'metallicity', metallicity, error)
  end subroutine

  subroutine set_time_step(self, id, time_step, error)
    class(stars), intent(inout) :: self
    integer, intent(in) :: id
    real(c_double), intent(in) :: time_step
    integer :: error
    call set_property_double(self, id, 'time_step', time_step, error)
  end subroutine


end module store_stars
