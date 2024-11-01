module PARAT
  use precision
  use param
  implicit none
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  real(WP), dimension(:,:,:), pointer :: ZMIX
  real(WP), dimension(:,:,:), pointer :: VISC
  real(WP), dimension(:,:,:), pointer :: DIFF
  real(WP), dimension(:,:,:), pointer :: PROG
  real(WP), dimension(:,:,:), pointer :: ENTH
  real(WP), dimension(:,:,:), pointer :: NOX

end module PARAT

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine PARAT_grid
  use PARAT
  use parser
  use math
  implicit none

  integer :: i,j,k
  integer :: jet_points, pilot_points, inner_wall_points, outer_wall_points, inflow_pts
  real(WP) :: jet_diameter, pilot_diameter
  real(WP) :: inner_wall, outer_wall, inflow_length
  real(WP) :: length, height, width, xs, s, ss
  integer :: mk1,mk2,mk3,mk4

  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  call parser_read('Combustor length',length)
  call parser_read('Combustor full height',height)
  call parser_read('Combustor full width',width)

  ! Geometrical parameters
  call parser_read('Jet diameter',jet_diameter)
  call parser_read('Pilot diameter',pilot_diameter)
  call parser_read('Inner wall thickness',inner_wall)
  call parser_read('Outer wall thickness',outer_wall)
  call parser_read('Inflow length',inflow_length)

!  jet_diameter   = 0.0182118_WP
!  pilot_diameter = 0.0285496_WP
!  inner_wall     = 0.00635_WP
!  outer_wall     = 0.00635_WP

  ! Points per feature
  call parser_read('Jet points',jet_points)
  call parser_read('Pilot points',pilot_points)
  call parser_read('Inner wall points',inner_wall_points)
  call parser_read('Outer wall points',outer_wall_points)
  call parser_read('Inflow points',inflow_pts)

  ! Stretching
  call parser_read('Stretching',s)
  call parser_read('Wall stretching',ss)
  call parser_read('X Stretch',xs)

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

  ! Create the x-grid
  ! Inflow (1D)
  do i=1,inflow_pts
     x(i) = (i-1)*(inflow_length)/inflow_pts - inflow_length
     !print*, 'inflow',i,x(i)
  end do
  ! Jet exit
  x(inflow_pts+1) = 0
  ! Flame
  do i=inflow_pts+2,nx+1
     x(i) = length*sinh((xs)*real(i-(inflow_pts+1),WP)/real((nx-inflow_pts),WP))/sinh(xs)
     !print*, 'combustor',i,x(i)
  end do

  ! Create the y-grid
  ! Fuel Jet
  y(1) = 0.0_WP
  do j=2,jet_points
     y(j) = jet_diameter/2.0_WP*tanh(s*(real(j,WP)-1.0_WP)/(real(jet_points,WP)-1.0_WP))/tanh(s)
  end do

  mk1 = jet_points

  ! Inner Wall
  do j=jet_points+1,jet_points+inner_wall_points
     y(j) = inner_wall*sinh((ss)*real(j-jet_points,WP)/real(inner_wall_points,WP))/sinh(ss) + jet_diameter/2.0_WP
  end do

  mk2 = jet_points+inner_wall_points

  ! Pilot
  do j=jet_points+inner_wall_points+1,jet_points+inner_wall_points+pilot_points
     y(j) = y(j-1) + (pilot_diameter/2.0_WP-jet_diameter/2.0_WP-inner_wall)/real(pilot_points,WP)
     
     !y(j) = 0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP* &
          !tanh(2.0_WP*(2.0_WP*real(j-25,WP)/real(32,WP)-1.0_WP))/tanh(2.0_WP) + &
          !inner_wall+jet_diameter/2.0_WP + &
          !0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP
  end do

  mk3 = jet_points+inner_wall_points+pilot_points

  ! Outer Wall
  print*, jet_points+inner_wall_points+pilot_points, outer_wall_points, pilot_diameter/2.0_WP
  do j=jet_points+inner_wall_points+pilot_points+1,jet_points+inner_wall_points+pilot_points+outer_wall_points
     y(j) = outer_wall*real(j-(jet_points+inner_wall_points+pilot_points),WP)/real(outer_wall_points,WP) + pilot_diameter/2.0_WP
  end do

  mk4 = jet_points+inner_wall_points+pilot_points+outer_wall_points

  ! Coflow
  do j=jet_points+inner_wall_points+pilot_points+outer_wall_points+1,ny+1
     y(j) = 9.792449895477940e-07_WP*1.094038245134641_WP**j + 0.013809808609566_WP
  end do

  do i=1,nx+1
     print*, i, x(i)
  end do

  do j=1,ny+1
     print*, j, y(j)
  end do

  ! Create the z-grid
  do k=1,nz+1
     z(k) = real(k-1,WP)*twoPi/real(nz,WP)
  end do

  ! Create the mask
  mask = 0
  
  do i=1,nx
     do j=1,ny
        if (j.lt.mk2 .and. & 
            j.ge.mk1 .and. & 
            x(i).lt.0.0_WP) then
           mask(i,j) = 1
        end if
        if (j.lt.mk4 .and. & 
             j.ge.mk3 .and. & 
             x(i).lt.0.0_WP) then
           mask(i,j) = 1
        end if
     end do
  end do

  return
end subroutine PARAT_grid


! =============================== !
! Create the variable array: Data !
! =============================== !
subroutine PARAT_data
  use PARAT
  use parser
  implicit none

  real(WP) :: rho_init
  
  ! Allocate the array data
  nvar = 12
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Initialize the variables and names: Gas-phase
  U    => data(:,:,:,1); names(1) = 'U'
  V    => data(:,:,:,2); names(2) = 'V'
  W    => data(:,:,:,3); names(3) = 'W'
  P    => data(:,:,:,4); names(4) = 'P'
  RHO  => data(:,:,:,5); names(5) = 'RHO'
  dRHO => data(:,:,:,6); names(6) = 'dRHO'
  ZMIX => data(:,:,:,7); names(7) = 'ZMIX'
  VISC => data(:,:,:,8); names(8) = 'VISC'
  DIFF => data(:,:,:,9); names(9) = 'DIFF'
  PROG => data(:,:,:,10); names(10) = 'PROG'
  ENTH => data(:,:,:,11); names(11) = 'ENTH'
  NOX => data(:,:,:,12); names(12) = 'NOX'
  
  ! Create them
  U    = 0.0_WP
  V    = 0.0_WP
  W    = 0.0_WP
  P    = 0.0_WP
  call parser_read('Initial density',rho_init)
  RHO  = rho_init
  dRHO = 0.0_WP
  ZMIX = 0.0_WP
  VISC = 0.0_WP
  DIFF = 0.0_WP
  PROG = 0.0_WP
  ENTH = 0.0_WP
  NOX = 0.0_WP
  
  return
end subroutine PARAT_data

