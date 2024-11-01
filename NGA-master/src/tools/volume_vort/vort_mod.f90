module vort_mod
  use parallel
  use string
  use precision
  use geometry

  implicit none
  integer :: nx1, ny1, nz1
  real(WP), dimension(:,:,:), pointer :: data8, ddata, vort_x, vort_y, vort_z, vort
  real(WP), dimension(:), pointer :: x1, y1, z1
  real(WP) :: xp1, xp2, p1, p2, num, den

end module vort_mod
