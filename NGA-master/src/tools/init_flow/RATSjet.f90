module RATSjet
  use precision
  use param
  implicit none

  ! Length, height, and width of the combustor
  real(WP) :: length,height,width
  ! Enclosed or not
  logical :: enclosed
  ! Height of the channel
  real(WP) :: channel_height_1,channel_height_2
  ! Thickness of the wall/number of points in wall
  real(WP) :: wall_thickness_1,wall_thickness_2
  ! Inflow length
  real(WP) :: inflow_length
  ! Number of points in each part
  integer :: channel_points_1,channel_points_2
  integer :: wall_points_1,wall_points_2
  integer :: inflow_pts
  
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

  
end module RATSjet

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine RATSjet_grid
  use RATSjet
  use parser
  implicit none
  
  integer :: i,j,k
  real(WP), parameter :: piovertwo = 2.0_WP*atan(1.0_WP)
  real(WP) :: s,dytmp,xs
  integer :: mk1,mk2,mk3,mk4,mk5,mk6,mk7,mk8
  
  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  call parser_read('Combustor length',length)
  call parser_read('Combustor full height',height)
  call parser_read('Combustor full width',width)
  
  call parser_read('Channel 1 full height',channel_height_1)
  call parser_read('Channel 2 full height',channel_height_2)
  call parser_read('Wall 1 thickness',wall_thickness_1)
  call parser_read('Wall 2 thickness',wall_thickness_2)
  call parser_read('Inflow length',inflow_length)

  call parser_read('Channel 1 points',channel_points_1)
  call parser_read('Channel 2 points',channel_points_2)
  call parser_read('Stretching',s)
  call parser_read('X Stretch',xs)
  call parser_read('Wall 1 points',wall_points_1)
  call parser_read('Wall 2 points',wall_points_2)
  call parser_read('Inflow points',inflow_pts)
  
  ! Set the periodicity
  xper = 0
  yper = 0
  zper = 1

  ! Cylindrical
  icyl = 0

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))
  
  ! Create the grid
  ! Inflow
  do i=1,inflow_pts
     !x(i) = inflow_length*tanh(xs*real(i-1.0_WP)/real(inflow_pts,WP))/tanh(xs)-inflow_length
     x(i) = (i-1)*(inflow_length)/inflow_pts - inflow_length
     print*, 'inflow',i,x(i)
  end do
  ! Jet exit
  x(inflow_pts+1) = 0
  ! "Combustor"
  do i=inflow_pts+2,nx+1
     x(i) = length*sinh((xs)*real(i-(inflow_pts+1),WP)/real((nx-inflow_pts),WP))/sinh(xs)
     print*, 'combustor',i,x(i)
  end do
  ! Central jet
  do j=ny/2+2,ny/2+channel_points_1/2+1
     !Refined
     y(j) = channel_height_1/2.0_WP*tanh(s*(2.0_WP*real(j-(ny/2+1-channel_points_1/2),WP)/real(channel_points_1,WP)-1.0_WP))/tanh(s)
     print*, 'centraljet',j, y(j)
     mk1 = ny/2+channel_points_1/2+1
     mk5 = ny+2-mk1
  end do
  ! Wall around central jet
  do j=ny/2+channel_points_1/2+2,ny/2+channel_points_1/2+wall_points_1+1
     ! Unrefined
     y(j) = wall_thickness_1*sinh((s+0.4_WP)*(real(j-(ny/2+1+channel_points_1/2),WP))/real(wall_points_1,WP))/sinh(s+0.4_WP)+channel_height_1/2.0_WP
     ! Relax center stretching
     !y(j) = wall_thickness_1*sinh((s)*(real(j-(ny/2+1+channel_points_1/2),WP))/real(wall_points_1,WP))/sinh(s)+channel_height_1/2.0_WP
     !y(j) = y(j-1) + wall_thickness_1/real(wall_points_1,WP)
     print*, 'Wall 1', j,y(j)
     mk2 = ny/2+channel_points_1/2+wall_points_1+1
     mk6 = ny+2-mk2
  end do
  ! Auxiliary jet
  do j=ny/2+channel_points_1/2+wall_points_1+2,ny/2+channel_points_1/2+wall_points_1+channel_points_2+1
     ! No stretching in this region
     y(j) = y(j-1) + channel_height_2/real(channel_points_2,WP)
    ! y(j) = channel_height_2/2.0_WP*tanh(s*(2.0_WP*real(j-(ny/2+1+channel_points_1/2+wall_points_1+channel_points_2/2))/real(channel_points_2,WP)))/tanh(s)+channel_height_1/2.0_WP+wall_thickness_1+channel_height_2/2.0_WP
     print*, 'outerjet', j,y(j)
     mk3 = ny/2+channel_points_1/2+wall_points_1+channel_points_2+1
     mk7 = ny+2-mk3
  end do
  ! Wall around aux jet
  do j=ny/2+channel_points_1/2+wall_points_1+channel_points_2+2,ny/2+channel_points_1/2+wall_points_1+channel_points_2+wall_points_2+1
     y(j) = y(j-1) + wall_thickness_2/real(wall_points_2,WP)
     print*, 'Wall2', j,y(j)
     mk4 = ny/2+channel_points_1/2+wall_points_1+channel_points_2+wall_points_2+1
     mk8 = ny+2-mk4
  end do
  ! Coflow
  do j=ny/2+channel_points_1/2+wall_points_1+channel_points_2+wall_points_2+2,ny+1
     ! 256 pts, width 0.15 m
     !y(j) = 0.0000016759385446_WP*1.0422537622_WP**j + 0.00526669367_WP
     ! 192 pts, width 0.1 m 
     !y(j) = 0.000000001033274214589254_WP*1.095041626036468_WP**j + 0.007895658898731_WP
     ! 384 pts, width 0.1 m w/ proper y+
     !y(j) = 0.0000002718041068171145_WP*1.031593916782195_WP**j + 0.006834833721650_WP
     ! 384 pts, width 0.1 m
     !y(j) = 0.000000001307843850259802_WP*1.045930083996151_WP**j + 0.007822777767870_WP
     ! 256 pts, width 0.2 m
     !y(j) = 0.000000365893335_WP*1.0496637244799_WP**j + 0.0059729158_WP
     ! 256 pts, width 0.3 m
     y(j) = 0.00000005180513011717143_WP*1.059415773381891_WP**j + 0.006633890487051_WP
     print*, 'coflow', j,y(j)
  end do
  ! Set mid-point to zero
  y(ny/2+1) = 0.0_WP
  ! Mirror the y-grid
  do j=ny/2+2,ny+1
     y(ny+2-j) = -y(j)
  end do
  do k=1,nz+1
     z(k) = (k-1)*width/nz - 0.5_WP*width
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
  do i=1,nx
     do j=1,ny/2
        if (j.lt.mk5 .and. & !y(j).lt.-0.5_WP*channel_height_1 .and. &
            j.ge.mk6 .and. & !y(j+wall_points_1).ge.-0.5_WP*channel_height_1 .and. &
            x(i).lt.0.0_WP) then
           mask(i,j) = 1
        end if
        if (j.lt.mk7 .and. & !y(j).lt.(-0.5_WP*channel_height_1-wall_thickness_1-channel_height_2) .and. &
             j.ge.mk8 .and. & !y(j+wall_points_2).ge.(-0.5_WP*channel_height_1-wall_thickness_1-channel_height_2) .and. &
             x(i).lt.0.0_WP) then
           mask(i,j) = 1
        end if
     end do
     do j=ny/2+1,ny
        if (j.ge.mk1 .and. & !y(j).gt.0.5_WP*channel_height_1 .and. &
            j.lt.mk2 .and. & !y(j-wall_points_1).le.0.5_WP*channel_height_1 .and. &
            x(i).lt.0.0_WP) then
           mask(i,j) = 1
        end if
        if (j.ge.mk3 .and. & !y(j).gt.(0.5_WP*channel_height_1+wall_thickness_1+channel_height_2) .and. &
             j.lt.mk4 .and. & !y(j-wall_points_2).le.(0.5_WP*channel_height_1+wall_thickness_1+channel_height_2) .and. &
             x(i).lt.0.0_WP) then
           mask(i,j) = 1
        end if
     end do
  end do
  
  return
end subroutine RATSjet_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine RATSjet_data
  use RATSjet
  use parser
  implicit none

  character(str_medium) :: data_type
  real(WP) :: rho_init

  ! Data type
  call parser_read('Data type',data_type,'cold')

  select case(trim(adjustl(data_type)))
  case ('cold')
     ! Allocate the array data
     nvar = 4
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U => data(:,:,:,1); names(1) = 'U'
     V => data(:,:,:,2); names(2) = 'V'
     W => data(:,:,:,3); names(3) = 'W'
     P => data(:,:,:,4); names(4) = 'P'
  
     ! Create them
     U = 0.0_WP
     V = 0.0_WP
     W = 0.0_WP
     P = 0.0_WP
  case ('passive mixing')
     ! Allocate the array data
     nvar = 5
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     ZMIX => data(:,:,:,5); names(5) = 'ZMIX'
  
     ! Create them
     U    = 0.0_WP
     V    = 0.0_WP
     W    = 0.0_WP
     P    = 0.0_WP
     ZMIX = 0.0_WP
  case ('cold mixing')
     ! Allocate the array data
     nvar = 7
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
     ZMIX => data(:,:,:,7); names(7) = 'ZMIX'
  
     ! Create them
     U    = 0.0_WP
     V    = 0.0_WP
     W    = 0.0_WP
     P    = 0.0_WP
     call parser_read('Initial density',rho_init)
     RHO  = rho_init
     dRHO = 0.0_WP
     ZMIX = 0.0_WP

case ('hot combustion')
     ! Allocate the array data
     nvar = 12
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
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
  case default
     print*, "Data type not recognized..."
  end select

  return
end subroutine RATSjet_data

!!$  dytmp = y(ny/2+channel_points_1/2+wall_points_1+channel_points_2+wall_points_2+1)
!!$  s = log((wall_thickness_2/real(wall_points_2,WP))/(height/2.0_WP-channel_height_1/2.0_WP-wall_thickness_1-channel_height_2-wall_thickness_2))/log(1.0_WP/real(ny/2-channel_points_1/2-wall_points_1-channel_points_2-wall_points_2,WP))
!!$  do j=ny/2+channel_points_1/2+wall_points_1+channel_points_2+wall_points_2+2,ny+1
!!$     y(j) = dytmp+(height/2.0_WP-dytmp)*(real(j-ny/2-channel_points_1/2-wall_points_1-channel_points_2-wall_points_2-1,WP)/real(ny/2-channel_points_1/2-wall_points_1-channel_points_2-wall_points_2,WP))**s
