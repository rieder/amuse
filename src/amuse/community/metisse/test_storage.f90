program test_store_stars
  use iso_c_binding
  use store_stars
  implicit none

  type(stars) :: star_system
  integer :: new_ids(2)
  integer :: ids_to_remove(2)
  integer :: i, number_of_stars

  ! Test new_star
  new_ids(1) = star_system%new_star(1.0_c_double)
  new_ids(2) = star_system%new_star(2.0_c_double)
  if (new_ids(1) /= 1 .or. new_ids(2) /= 2) then
    error stop "new_star failed"
  end if

  ! Test remove_star
  ids_to_remove = [1, 2]
  call star_system%remove_star(ids_to_remove(1))
  call star_system%remove_star(ids_to_remove(2))
  call star_system%get_number_of_stars(number_of_stars)
  if (number_of_stars /= 0) then
    write (*, *) "Number of stars: ", number_of_stars
    error stop "remove_star failed"
  end if

  ! Test edge cases
  new_ids = star_system%new_star(1.0_c_double)
  if (new_ids(1) /= 1) then
    error stop "new_star failed with single star"
  end if

  ids_to_remove(1) = 1
  call star_system%remove_star(ids_to_remove(1))
  call star_system%get_number_of_stars(number_of_stars)  
  if (number_of_stars /= 0) then
    error stop "remove_star failed with single star"
  end if

  write (*, *) "All tests passed"

end program test_store_stars
