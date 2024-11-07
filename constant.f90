module constants
  implicit none
  integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
  & wp = selected_real_kind(15)
  real(wp), parameter :: pi = 4*atan(1._wp)
end module constants
