module planeStats_metric
  use planeStats
  implicit none

  ! Grid information
  integer :: nx1, ny1, nz1
  integer :: xper, yper, zper, icyl
  integer, parameter :: stp=1
  real(WP), dimension(:), pointer :: x_grid,  y_grid,  z_grid
  real(WP), dimension(:), pointer :: xm_grid, ym_grid, zm_grid
  real(WP), dimension(:), pointer :: dxi
  real(WP), dimension(:), pointer :: dyi
  real(WP), dimension(:), pointer :: dxmi
  real(WP), dimension(:), pointer :: dymi
  real(WP) :: dzi  
  real(WP), dimension(:,:,:), pointer :: grad_x, grad_y, grad_z

end module planeStats_metric


! ========================================================== !
! Initialize the module                                      !
! ========================================================== !
subroutine planeStats_metric_init
  use planeStats_metric
  use parser
  implicit none

  character(len=str_medium) :: fconfig, config
  integer :: iunit, ierr
  integer :: i,j,k

  ! Very important - index conventions
  !  -- From data file:   y_grid(1:ny1)
  !  -- From plane file : y(1:ny)

  ! Read the source config file
  call parser_read('Configuration file',fconfig)
  fconfig = trim(data_dir)//'/'//trim(fconfig)
  call BINARY_FILE_OPEN(iunit,trim(fconfig),"r",ierr)
  call BINARY_FILE_READ(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_READ(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_READ(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_READ(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_READ(iunit,nx1,1,kind(nx1),ierr)
  call BINARY_FILE_READ(iunit,ny1,1,kind(ny1),ierr)
  call BINARY_FILE_READ(iunit,nz1,1,kind(nz1),ierr)
  allocate(x_grid(nx1+1),y_grid(ny1+1),z_grid(nz1+1))
  call BINARY_FILE_READ(iunit,x_grid,nx+1,kind(x_grid),ierr)
  call BINARY_FILE_READ(iunit,y_grid,ny+1,kind(y_grid),ierr)
  call BINARY_FILE_READ(iunit,z_grid,nz+1,kind(z_grid),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Create the midpoints - using grid x- and y-indices
  allocate(xm_grid(nx1),ym_grid(ny1),zm_grid(nz1))
  do i=1,nx1
     xm_grid(i) = 0.5_WP*(x_grid(i)+x_grid(i+1))
  end do
  do j=1,ny1
     ym_grid(j) = 0.5_WP*(y_grid(j)+y_grid(j+1))
  end do
  do k=1,nz1
     zm_grid(k) = 0.5_WP*(z_grid(k)+z_grid(k+1))
  end do

  ! Compute the inverses
  allocate(dxi(1:nx1))
  allocate(dyi(1:ny1))
  allocate(dxmi(1:nx1-1))
  allocate(dymi(1:ny1-1))
  do i=1,nx1
     dxi(i) = 1.0_WP/(x_grid(i+1)-x_grid(i))
  end do
  do i=1,nx1-1
     dxmi(i) = 1.0_WP/(xm_grid(i+1)-xm_grid(i))
  end do
  do j=1,ny1
     dyi(j) = 1.0_WP/(y_grid(j+1)-y_grid(j))
  end do
  do j=1,ny1-1
     dymi(j) = 1.0_WP/(ym_grid(j+1)-ym_grid(j))
  end do
  dzi = 1.0_WP/(z_grid(2)-z_grid(1))

  ! need to match x-location with corresponding plane - has been interpolated!
  ! also, j-indices in grid don't correspond to j-indices in planes (only ny/2 is at same location)

  ! Allocate the gradients - use plane x- and y-indices
!!$  allocate(grad_x(pnmin_:pnmax_,1:ny,-stp:stp))
!!$  allocate(grad_y(pnmin_:pnmax_,1:ny,-stp:stp))
!!$  allocate(grad_z(pnmin_:pnmax_,1:ny,-stp:stp))
!!$  print *, ny, ny1
!!$  do j=1,ny
!!$     print *, y(j), y_grid(j+(ny1-ny)/2)
!!$  end do

  return
end subroutine planeStats_metric_init
