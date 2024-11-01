module simplejet
  use precision
  use param
  implicit none
  
  ! Instructions

  ! Length and diameter of the combustor
  real(WP) :: length, diameter
  ! Diameter of the pipe
  real(WP) :: jet_diameter
  ! Inflow length
  real(WP) :: inflow_length
  ! Other variables
  real(WP) :: wall_thickness, dy_wall, dx_inflow 
  ! variables for stretched grids
  logical :: stretch_jet, stretch_coflow, stretch_axial
  real(WP) :: ratio_j, ratio_c, ratio_a, A_j, A_c, A_a, C_j, C_c, C_a

  ! Number of points in radial direction
  integer :: inflow_points, jet_points, wall_points
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: ZMIX

  contains 
    subroutine calc_stretch_params(x0,xn,dx,dx1orN,nx,A,b,C)

      ! Calculates stretch parameters A b and C that give stretching from x0 to xn over nx points
      ! with a base dx as specified. If dx1orN is true the base dx is the step from the first to 
      ! second points, if dx1orN is false the base dx is the step from the second last to last.
      ! A guess value for b must be provided because a Newton solver is used

      implicit none 
      logical,  intent(in)    :: dx1orN
      integer,  intent(in)    :: nx
      real(WP), intent(in)    :: x0,xn,dx
      real(WP), intent(inout) :: b
      real(WP), intent(out)   :: A, C

      real(WP) :: error, dx_test, dEdb, error_old, b_old
      real(WP), parameter :: tol = 1e-10_WP
      integer :: ii 
      
      ii = 0
      error = huge(1.0_WP)
      error_old = 0.0_WP
      b_old = b - 1.0e-8_WP
      dEdb = 1.0e6_WP
      do while (error .ge. tol .and. ii .le. 101)
         ii = ii+1
         A = (xn - x0)/(b**nx - 1.0_WP)
         C = A - x0
         if (dx1orN) then 
            dx_test = (A*b-C) - x0
         else 
            dx_test = xn - (A*b**real(nx-1,WP)-C)
         end if
         error = dx_test-dx
         dEdb = (error-error_old)/(b - b_old)
         error_old = error
         b_old = b
         b = b - error/dEdb
         error = abs(error)
      end do
      A = (xn - x0)/(b**nx - 1.0_WP)
      C = A - x0
      
      if (ii.ge.100) print *, 'Warning: Newton iteration did not converge'

    end subroutine calc_stretch_params

end module simplejet

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine simplejet_grid
  use simplejet
  use parser
  implicit none

  integer :: i,j,k
  real(WP), parameter :: twoPi = 2.0_WP*acos(-1.0_WP)

  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  call parser_read('Domain length'  ,length  )
  call parser_read('Domain diameter',diameter)

  call parser_read('Jet diameter' ,jet_diameter)

  call parser_read('Inflow length',inflow_length)

  call parser_read('Inflow points',inflow_points)
  call parser_read('Jet points'   ,jet_points   )
  call parser_read('Wall points'  ,wall_points  )
  call parser_read('Wall thickness'  ,wall_thickness  )
    
  call parser_read('Axial stretching' , stretch_axial ,.false.)
  call parser_read('Jet stretching'   , stretch_jet   ,.false.)
  call parser_read('Coflow stretching', stretch_coflow,.false.)

  ! Set the periodicity
  xper = 0
  yper = 0
  zper = 1

  ! Cylindrical
  icyl = 1

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))

  ! Create the grid

  ! Inflow
  dx_inflow = inflow_length/real(inflow_points,WP)
  do i=1,inflow_points+1
     x(i) = real(i-1,WP)*dx_inflow- inflow_length
  end do
  ! Remainder of domain
  if (stretch_axial) then 
     ratio_a = 1.01_WP
     call calc_stretch_params(0.0_WP,length,dx_inflow,.true.,nx-inflow_points,A_a,ratio_a,C_a)
     do i=inflow_points+1,nx+1
        x(i) = A_a * ratio_a**real(i-1-inflow_points) - C_a
     end do
     print *,'Stretched axial direction, ratio = ', ratio_a
     print *,'dx_inf = ', dx_inflow
     print *,'dx_min = ', x(inflow_points+2) - x(inflow_points+1)
     print *,'dx_max = ', x(nx+1) - x(nx)
  else
     do i=inflow_points+1,nx+1
        x(i) = real(i-1-inflow_points)*length/real(nx-inflow_points,WP) 
     end do
     print *,'Unstretched axial direction'
     print *,'dx_inflow = ', dx_inflow
     print *,'dx_domain = ', length/real(nx-inflow_points,WP)
  end if

  print *, ''

  ! Central jet
  dy_wall = wall_thickness/real(wall_points,WP)
  if (stretch_jet) then
     ratio_j = 1.01_WP
     call calc_stretch_params(0.0_WP,0.5_WP*jet_diameter,dy_wall,.false.,jet_points,A_j,ratio_j,C_j)
     do j=1,jet_points+1
        y(j) = A_j * ratio_j**real(j-1) - C_j
     end do
     print *, 'Stretched radial direction in jet, ratio = ', ratio_j
     print *, 'dy_wall = ', dy_wall
     print *, 'dy_min  = ', y(jet_points+1) - y(jet_points)
     print *, 'dy_max  = ', y(2) - y(1)
  else
     do j=1,jet_points+1
        y(j) = real(j-1,WP)*jet_diameter/(2.0_WP*real(jet_points,WP))
     end do
     print *,'dy_jet = ', jet_diameter/(2.0_WP*real(jet_points,WP)), '(unstretched)'
     print *,'dy_wall= ', dy_wall
  end if

  ! Wall
  do j=jet_points+1,jet_points+wall_points+1
     y(j) = real(j-jet_points-1,WP) * dy_wall +jet_diameter*0.5_WP
  end do

  print *, ''

  ! Coflow
  if (stretch_coflow) then
     ratio_c = 1.01_WP
     call calc_stretch_params(0.5_WP*jet_diameter + wall_thickness, 0.5_WP*diameter, dy_wall, &
          .true., ny-jet_points-wall_points, A_c, ratio_c, C_c)
     do j=jet_points+wall_points+1,ny+1
        y(j) = A_c * ratio_c ** (j - jet_points - wall_points-1) - C_c
     end do
     print *,'Stretched radial direction in coflow, ratio = ', ratio_c
     print *, 'dy_min = ', y(jet_points+wall_points+2) - y(jet_points+wall_points+1)
     print *, 'dy_max = ', y(ny+1) - y(ny)
  else
     do j=jet_points+wall_points+1,ny+1
        y(j) = real(j - jet_points - wall_points -1) &
             * (0.5_WP*diameter-jet_diameter*0.5_WP-wall_thickness)/real(ny-wall_points-jet_points) &
             + jet_diameter*0.5_WP + wall_thickness
     end do
     print *, 'dy_coflow', (0.5_WP*diameter-jet_diameter*0.5_WP-wall_thickness)/real(ny-wall_points-jet_points), &
          '(unstretched)'
  end if

  print *, ''

  ! Circumferential
  do k=1,nz+1
     z(k) = real(k-1,WP)*twoPi/real(nz,WP)
  end do

  ! Create the mid points
  do i=1,nx
     xm(i)= 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j)= 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k)= 0.5_WP*(z(k)+z(k+1))
  end do

  ! Create the masks
  mask = 0
  mask(1:inflow_points,jet_points+1:jet_points+wall_points) = 1

  return
end subroutine simplejet_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine simplejet_data
  use simplejet
  use parser
  implicit none

  ! Allocate the array data
  nvar = 16
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  data(:,:,:, 1) = 0.0_WP;      names( 1) = 'U'
  data(:,:,:, 2) = 0.0_WP;      names( 2) = 'V'
  data(:,:,:, 3) = 0.0_WP;      names( 3) = 'W'
  data(:,:,:, 4) = 0.0_WP;      names( 4) = 'P'
  data(:,:,:, 5) = 1.179640_WP; names( 5) = 'RHO'
  data(:,:,:, 6) = 0.0_WP;      names( 6) = 'dRHO'
  data(:,:,:, 7) = 0.0_WP;      names( 7) = 'VISC'
  data(:,:,:, 8) = 0.0_WP;      names( 8) = 'DIFF'
  data(:,:,:, 9) = 0.0_WP;      names( 9) = 'ZMIX'
  data(:,:,:,10) = 0.0_WP;      names(10) = 'ZMIX2'
  data(:,:,:,11) = 0.0_WP;      names(11) = 'LM_VEL' 
  data(:,:,:,12) = 1.0_WP;      names(12) = 'MM_VEL' 
  data(:,:,:,13) = 0.0_WP;      names(13) = 'LM_ZMIX2' 
  data(:,:,:,14) = 1.0_WP;      names(14) = 'MM_ZMIX2' 
  data(:,:,:,15) = 0.0_WP;      names(15) = 'LM_ZMIX' 
  data(:,:,:,16) = 1.0_WP;      names(16) = 'MM_ZMIX' 

  return

end subroutine simplejet_data

