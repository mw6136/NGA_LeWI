module jetcart2
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
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  real(WP), dimension(:,:,:), pointer :: ZMIX
  
end module jetcart2

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine jetcart2_grid
  use jetcart2
  use parser
  implicit none
  
  integer :: i,j,k
  real(WP), parameter :: piovertwo = 2.0_WP*atan(1.0_WP)
  real(WP) :: s,dytmp,sx,dxmin_in,dxmin
  integer :: i0
  real(WP) :: dy_start, dy_end, dy, dy_prev, dym, hw
  real(WP), dimension(0:1) :: y_tmp
  integer :: iunit, ierr
  logical :: uniform
  
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
  call parser_read('Wall 1 points',wall_points_1)
  call parser_read('Wall 2 points',wall_points_2)

  call parser_read('Minimum dx',dxmin_in)
  call parser_read('y-Stretching',s)

  ! For debugging metrics
  call parser_read('Uniform',uniform,.false.)
  
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
  
  ! x-grid (JFM 2/24/2016)
  i0 = nint(inflow_length/dxmin_in)
  ! Rescale dxmin to get the correct inflow length with an integer number of points
  dxmin = inflow_length/i0
  ! Stretching factor
  sx = (log(length)-log(dxmin))/log(real(nx-i0,WP))
  ! Create the grid
  do i=0,nx
     if (i.le.i0) then
        x(i+1) = -(i0-i)*dxmin
     else
        x(i+1) = (i-i0)**sx * length/(nx-i0)**sx
     end if
  end do
  ! Display some info
  if (irank.eq.iroot) then
     print *, '========== x-grid =========='
     print *, 'Original dx_min: ', dxmin_in
     print *, 'Rescaled dx_min: ', dxmin
     print *, 'Number of cells in inflow:    ', i0
     print *, 'Number of cells in combustor: ', nx-i0
     if (i0.ge.nx) print *, 'WARNING: Number of points in inflow greater than nx!'
     print *, 'x-Stretching: ', sx
     if (sx.lt.1) print *, 'WARNING: Stretching factor less than 1!'
     print *, 'Inflow    -- dx_max: ', x(2)-x(1)
     print *, '             dx_min: ', x(i0+1)-x(i0)
     print *, 'Combustor -- dx_min: ', x(i0+2)-x(i0+1)
     print *, '             dx_max: ', x(nx+1)-x(nx)
  end if

!!$  i0 = nx*inflow_length**(1/s)/(length**(1/s)+inflow_length**(1/s))
!!$  do i=0,nx
!!$     if (i.le.i0) then
!!$        x(i+1) = -(i0-i)**s*inflow_length/i0**s
!!$     else
!!$        x(i+1) = (i-i0)**s*length/(nx-i0)**s
!!$     end if
!!$  do i=1,nx+1
!!$     x(i) = (i-1)*(length+inflow_length)/nx - inflow_length
!!$  end do

  ! Central jet
  do j=ny/2+2,ny/2+channel_points_1/2+1
     ! Hyperbolic tangent profile
     y(j) = channel_height_1/2.0_WP*tanh(s*(2.0_WP*real(j-(ny/2+1-channel_points_1/2),WP)/real(channel_points_1,WP)-1.0_WP))/tanh(s)
  end do
  if (irank.eq.iroot) print *, y(ny/2+channel_points_1/2+1)

  ! Last dy in central jet
  dy_start = y(ny/2+channel_points_1/2+1) - y(ny/2+channel_points_1/2)
  ! First dy in auxiliary jet
  do j=0,1
     y_tmp(j) = channel_height_2/2.0_WP * &
          tanh(s*(2.0_WP*real(j+1-(channel_points_2/2))/real(channel_points_2,WP))) / tanh(s) &
          + channel_height_1/2.0_WP+wall_thickness_1+channel_height_2/2.0_WP
  end do
  dy_end = y_tmp(0) - (channel_height_1/2.0_WP+wall_thickness_1)

  ! Compute initial stretching in primary wall
  dym = 0.5_WP*(dy_start+dy_end)
  wall_points_1 = ceiling(wall_thickness_1/dym)
  print *, 'dy_start=',dy_start
  print *, 'dy_end  =',dy_end
  print *, 'dym     =',dym
  print *, 'wall pts=',wall_points_1
  dy_prev = dy_start
  do j=ny/2+channel_points_1/2+2,ny/2+channel_points_1/2+wall_points_1+1
     dy = dy_start + (dy_end-dy_start)*real(j-(ny/2+channel_points_1/2+1),WP)/real(wall_points_1,WP)
     y(j) = y(j-1) + dy
     if (irank.eq.iroot) print '(i5,ES15.5,ES15.5,ES15.5)', j, y(j), y(j)-y(j-1), y(j)-y(j-1)-dy_prev
     dy_prev=y(j)-y(j-1)
  end do
  
  ! Rescale
  ! Mean grid error per point
  hw = ((channel_height_1/2.0_WP+wall_thickness_1)-y(ny/2+channel_points_1/2.0_WP+wall_points_1+1))/real(wall_points_1,WP)
  print *, 'rescale=',hw
  dy_prev = dy_start
  do j=ny/2+channel_points_1/2+2,ny/2+channel_points_1/2+wall_points_1+1
     dy = dy_start + (dy_end-dy_start)*real(j-(ny/2+channel_points_1/2+1),WP)/real(wall_points_1,WP) + hw
     y(j) = y(j-1) + dy
     if (irank.eq.iroot) print '(i5,ES15.5,ES15.5,ES15.5)', j, y(j), y(j)-y(j-1), y(j)-y(j-1)-dy_prev
     dy_prev=y(j)-y(j-1)
  end do


!!$  s_wall = (dy_end-dy_start)/(wall_thickness_1+dy_start+dy_end)*w1
!!$  ! Wall around central jet
!!$  dy_prev = dy_start
!!$  wall_points_1 = 0
!!$  hw = dy_start
!!$  j = ny/2+channel_points_1/2+1
!!$  if (irank.eq.iroot) print *, 'dy_start=', dy_start
!!$  do while (y(j).lt.(channel_height_1/2.0_WP+wall_thickness_1-w2*dy_end))
!!$     dy = s_wall*hw + dy_start
!!$     y(j+1) = y(j) + dy
!!$     wall_points_1 = wall_points_1 + 1
!!$     hw = hw + dy
!!$     if (irank.eq.iroot) print '(i5,ES15.5,ES15.5,ES15.5)', wall_points_1, y(j+1), y(j+1)-y(j), y(j+1)-y(j)-dy_prev
!!$     dy_prev = y(j+1)-y(j)
!!$     j = j + 1
!!$  end do

  if (irank.eq.iroot) print '(a8,ES15.5,ES15.5)', 'dy_end=', dy_end, dy_end-dy_prev

!!$
!!$  wall_points_1 = int(wall_thickness_1/dy_start)
!!$  
!!$  do j=ny/2+channel_points_1/2+2,ny/2+channel_points_1/2+wall_points_1+1
!!$     y(j) = y(j-1) + wall_thickness_1/real(wall_points_1,WP)
!!$  end do

  if (irank.eq.iroot) print *, y(ny/2+channel_points_1/2+wall_points_1+1)

  ! Auxiliary jet
  do j=ny/2+channel_points_1/2+wall_points_1+2,ny/2+channel_points_1/2+wall_points_1+channel_points_2+1
     y(j) = channel_height_2/2.0_WP*tanh(s*(2.0_WP*real(j-(ny/2+1+channel_points_1/2+wall_points_1+channel_points_2/2))/real(channel_points_2,WP)))/tanh(s)+channel_height_1/2.0_WP+wall_thickness_1+channel_height_2/2.0_WP
  end do

  ! Wall around aux jet
  do j=ny/2+channel_points_1/2+wall_points_1+channel_points_2+2,ny/2+channel_points_1/2+wall_points_1+channel_points_2+wall_points_2+1
     y(j) = y(j-1) + wall_thickness_2/real(wall_points_2,WP)
  end do

  ! Coflow
  dytmp = y(ny/2+channel_points_1/2+wall_points_1+channel_points_2+wall_points_2+1)
  s = log((wall_thickness_2/real(wall_points_2,WP))/(height/2.0_WP-channel_height_1/2.0_WP-wall_thickness_1-channel_height_2-wall_thickness_2))/log(1.0_WP/real(ny/2-channel_points_1/2-wall_points_1-channel_points_2-wall_points_2,WP))
  do j=ny/2+channel_points_1/2+wall_points_1+channel_points_2+wall_points_2+2,ny+1
     y(j) = dytmp+(height/2.0_WP-dytmp)*(real(j-ny/2-channel_points_1/2-wall_points_1-channel_points_2-wall_points_2-1,WP)/real(ny/2-channel_points_1/2-wall_points_1-channel_points_2-wall_points_2,WP))**s
  end do

  ! For testing - uniform grid
  if (uniform) then
     call parser_read('Wall 1 points',wall_points_1)
     call parser_read('Wall 2 points',wall_points_2)
     do i=1,nx+1
        x(i) = (i-1)*(length+inflow_length)/nx - inflow_length
     end do

!!$     do j=ny/2+2,ny/2+channel_points_1/2+1
!!$        y(j) = (j-ny/2+1-channel_points_1/2)*channel_height_1/channel_points_1
!!$     end do
!!$     do j=channel_points_1/2+2,ny+1
!!$        y(j) = y(j-1) + (height-channel_height_1)/(ny-channel_points_1)
!!$     end do
     do j=ny/2+1,ny+1
        y(j) = (j-1)*height/ny - 0.5_WP*height
     end do
  end if

  ! Set mid-point to zero
  y(ny/2+1) = 0.0_WP

  ! Mirror the y-grid
  do j=ny/2+2,ny+1
     y(ny+2-j) = -y(j)
  end do
  do k=1,nz+1
     z(k) = (k-1)*width/nz - 0.5_WP*width
  end do

  ! Print the grid to output
  if (irank.eq.iroot) then
     open(iunit,file='grid_output',form="formatted",iostat=ierr,status="REPLACE")
     do j=2,ny+1
        write(iunit, "(i5,ES15.5,ES15.5)") j, y(j), y(j)-y(j-1)
     end do
     close(iunit)
  end if

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
        if (y(j).lt.-0.5_WP*channel_height_1 .and. &
            y(j+wall_points_1).ge.-0.5_WP*channel_height_1 .and. &
            x(i+1).le.0.0_WP) then
           mask(i,j) = 1
        end if
        if (y(j).lt.(-0.5_WP*channel_height_1-wall_thickness_1-channel_height_2) .and. &
             y(j+wall_points_2).ge.(-0.5_WP*channel_height_1-wall_thickness_1-channel_height_2) .and. &
             x(i+1).le.0.0_WP) then
           mask(i,j) = 1
        end if
     end do
     do j=ny/2+1,ny
        if (y(j).gt.0.5_WP*channel_height_1 .and. &
            y(j-wall_points_1).le.0.5_WP*channel_height_1 .and. &
            x(i+1).le.0.0_WP) then
           mask(i,j-1) = 1
        end if
        if (y(j).gt.(0.5_WP*channel_height_1+wall_thickness_1+channel_height_2) .and. &
             y(j-wall_points_2).le.(0.5_WP*channel_height_1+wall_thickness_1+channel_height_2) .and. &
             x(i+1).le.0.0_WP) then
           mask(i,j-1) = 1
        end if
     end do
  end do
  
  return
end subroutine jetcart2_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine jetcart2_data
  use jetcart2
  use parser
  implicit none

  integer :: i,j,k
  character(str_medium) :: data_type
  real(WP) :: rho_init, T_init, U_init, V_init, W_init, P_init
  real(WP) :: U_jet, U_coflow, U_bulk
  real(WP), parameter :: R_cst   = 8314.34_WP   ! J/[kmol K]

  ! Data type
  call parser_read('Data type',data_type,'cold')

  select case(trim(adjustl(data_type)))
  case ('cold')
     ! Allocate the array data
     nvar = 4
     allocate(data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1); names(1) = 'U'
     V => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,2); names(2) = 'V'
     W => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,3); names(3) = 'W'
     P => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,4); names(4) = 'P'
  
     ! Create them
     U = 0.0_WP
     V = 0.0_WP
     W = 0.0_WP
     P = 0.0_WP

  case ('passive mixing')
     ! Allocate the array data
     nvar = 5
     allocate(data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1); names(1) = 'U'
     V    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,2); names(2) = 'V'
     W    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,3); names(3) = 'W'
     P    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,4); names(4) = 'P'
     ZMIX => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,5); names(5) = 'ZMIX'
  
     ! Create them
     U    = 0.0_WP
     V    = 0.0_WP
     W    = 0.0_WP
     P    = 0.0_WP
     ZMIX = 0.0_WP

  case ('cold mixing')
     ! Allocate the array data
     nvar = 7
     allocate(data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1); names(1) = 'U'
     V    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,2); names(2) = 'V'
     W    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,3); names(3) = 'W'
     P    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,4); names(4) = 'P'
     RHO  => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,5); names(5) = 'RHO'
     dRHO => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,6); names(6) = 'dRHO'
     ZMIX => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,7); names(7) = 'ZMIX'
  
     ! Create them
     U    = 0.0_WP
     V    = 0.0_WP
     W    = 0.0_WP
     P    = 0.0_WP
     call parser_read('Initial density',rho_init)
     RHO  = rho_init
     dRHO = 0.0_WP
     ZMIX = 0.0_WP

  case ('finite chem H2')
     ! Allocate the array data
     ! 7 + N species + T 
     nvar = 17
     allocate(data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1); names(1) = 'U'
     V    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,2); names(2) = 'V'
     W    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,3); names(3) = 'W'
     P    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,4); names(4) = 'P'
     RHO  => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,5); names(5) = 'RHO'
     dRHO => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,6); names(6) = 'dRHO'
     ZMIX => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,17); names(17) = 'ZMIX'
     names(7)  = 'N2'
     names(8)  = 'H'
     names(9)  = 'O2'
     names(10) = 'O'
     names(11) = 'OH'
     names(12) = 'H2'
     names(13) = 'H2O'
     names(14) = 'HO2'
     names(15) = 'H2O2'
     names(16) = 'T'
       
     ! Create them
     call parser_read('U',U_init)
     call parser_read('V',V_init)
     call parser_read('W',W_init)
     call parser_read('Pressure',P_init) ! Thermodynamic pressure to compute density
     U = U_init
     V = V_init
     W = W_init
     P = 0.0_WP !Hydrodynamic pressure
     data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,7)  = 0.768_WP
     ! Use 1.0e-30 for AF Jacobian initial step
     data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,8)  = 1.0e-30
     data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,9)  = 0.232_WP
     data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,10) = 1.0e-30
     data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,11) = 1.0e-30
     data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,12) = 1.0e-30
     data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,13) = 1.0e-30
     data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,14) = 1.0e-30
     data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,15) = 1.0e-30
    
     call parser_read('Initial temperature',T_init)
     data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,16) = T_init
     RHO  = P_init * 28.84_WP / (T_init * R_cst)
     dRHO = 0.0_WP
     ZMIX = 0.0_WP 


  case ('finite chem H2 prod')
     ! Allocate the array data
     ! 7 + N species + T 
     nvar = 17
     allocate(data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1); names(1) = 'U'
     V    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,2); names(2) = 'V'
     W    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,3); names(3) = 'W'
     P    => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,4); names(4) = 'P'
     RHO  => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,5); names(5) = 'RHO'
     dRHO => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,6); names(6) = 'dRHO'
     ZMIX => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,17); names(17) = 'ZMIX'
     names(7)  = 'N2'
     names(8)  = 'H'
     names(9)  = 'O2'
     names(10) = 'O'
     names(11) = 'OH'
     names(12) = 'H2'
     names(13) = 'H2O'
     names(14) = 'HO2'
     names(15) = 'H2O2'
     names(16) = 'T'
       
     ! Create them
     call parser_read('U jet',U_jet)
     call parser_read('U coflow',U_coflow)
     call parser_read('U bulk',U_bulk)
     call parser_read('V',V_init)
     call parser_read('W',W_init)
     V = V_init
     W = W_init
     P = 0.0_WP !Hydrodynamic pressure
     dRHO = 0.0_WP

     ! U velocity
     U = 0.0_WP
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           if (y(j).gt.-0.5_WP*channel_height_1.and.y(j).lt.0.5_WP*channel_height_1) then
              U(imin_:imax_,j,k) = U_jet
           else if ((y(j).gt.(-0.5_WP*channel_height_1-wall_thickness_1-channel_height_2).and.y(j).lt.(-0.5_WP*channel_height_1-wall_thickness_1)) &
              .or.  (y(j).lt.( 0.5_WP*channel_height_1+wall_thickness_1+channel_height_2).and.y(j).gt.( 0.5_WP*channel_height_1+wall_thickness_1))) then
              U(imin_:imax_,j,k) = U_coflow
           else if ((y(j).lt.(-0.5_WP*channel_height_1-wall_thickness_1-channel_height_2-wall_thickness_2)) &
              .or.  (y(j).gt.( 0.5_WP*channel_height_1+wall_thickness_1+channel_height_2+wall_thickness_2))) then
              U(imin_:imax_,j,k) = U_bulk
           end if
        end do
     end do

     ! Scalars
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           if (y(j).gt.-0.5_WP*channel_height_1.and.y(j).lt.0.5_WP*channel_height_1) then
              data(imin_:imax_,j,k,7)  = 8.095800e-01_WP
              data(imin_:imax_,j,k,8)  = 1.0e-60_WP
              data(imin_:imax_,j,k,9)  = 1.695440e-01_WP
              data(imin_:imax_,j,k,10) = 1.0e-60_WP
              data(imin_:imax_,j,k,11) = 1.0e-60_WP
              data(imin_:imax_,j,k,12) = 2.136260e-02_WP
              data(imin_:imax_,j,k,13) = 1.0e-60_WP
              data(imin_:imax_,j,k,14) = 1.0e-60_WP
              data(imin_:imax_,j,k,15) = 1.0e-60_WP
              data(imin_:imax_,j,k,16) = 300.0_WP
              RHO(imin_:imax_,j,k)  = 0.9073_WP
              ZMIX(imin_:imax_,j,k) = 1.0_WP
           else
              data(imin_:imax_,j,k,7)  = 8.095800e-01_WP
              data(imin_:imax_,j,k,8)  = 8.907410e-06_WP
              data(imin_:imax_,j,k,9)  = 4.112050e-03_WP
              data(imin_:imax_,j,k,10) = 5.892620e-05_WP
              data(imin_:imax_,j,k,11) = 1.434460e-03_WP
              data(imin_:imax_,j,k,12) = 2.662840e-04_WP
              data(imin_:imax_,j,k,13) = 1.845390e-01_WP
              data(imin_:imax_,j,k,14) = 3.167530e-07_WP
              data(imin_:imax_,j,k,15) = 7.001110e-08_WP
              data(imin_:imax_,j,k,16) = 2047.49_WP
              RHO(imin_:imax_,j,k)  = 0.150703_WP
              ZMIX(imin_:imax_,j,k) = 0.0_WP 
           end if
        end do
     end do

  case default
     print*, "Data type not recognized..."
  end select

  return
end subroutine jetcart2_data

