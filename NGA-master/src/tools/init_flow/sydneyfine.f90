module sydneyfine
  use precision
  use param
  implicit none
  
  ! Pointers to variable in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P

end module sydneyfine

! ==================== !
! Create the grid/mesh !
! ==================== !
subroutine sydneyfine_grid
  use sydney
  use parser
  use math
  implicit none

  integer :: i,j,k
  real :: jet_diameter, pilot_diameter
  real :: inner_wall, outer_wall

  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('ny',ny)
  call parser_read('nz',nz)

  ! Geometrical parameters
  jet_diameter   = 0.0075_WP
  pilot_diameter = 0.018_WP
  inner_wall     = 0.00025_WP
  outer_wall     = 0.00020_WP

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
  if (nx.eq.384) then
     ! Inflow (1D)
     do i=1,17
        x(i) = jet_diameter/16.0_WP*real(i-1,WP)-jet_diameter
     end do
     ! Flame
     do i=17,nx+1
        x(i) = 0.044953516_WP*1.006_WP**i-0.049765611_WP
     end do
  else if (nx.eq.400) then
     ! Inflow (1D)
     do i=1,31
        x(i) = jet_diameter/30.0_WP*real(i-1,WP)-jet_diameter
     end do
     ! Flame
     do i=31,nx+1
        x(i) = 0.030562376_WP*1.006659_WP**i-0.037544032_WP
     end do
  else if (nx.eq.416) then
     ! Inflow (1D)
     do i=1,49
        x(i) = jet_diameter/48.0_WP*real(i-1,WP)-jet_diameter
     end do
     ! Flame
     do i=49,nx+1
        x(i) = 0.020235226_WP*1.005813_WP**i-0.026881472_WP
     end do
  else if (nx.eq.424) then
     ! Inflow (1D)
     do i=1,49
        x(i) = jet_diameter/48.0_WP*real(i-1,WP)-jet_diameter
     end do
     ! Flame
     do i=49,nx+1
        x(i) = 0.012591563_WP*1.00828_WP**i-0.018860627_WP
     end do
  else  
     write(*,'(a)') 'nx must be 384 or 400 or 416 or 424'
     call kill
  end if

  do i=1,nx+1
     print*, i,x(i)
  end do

  ! Create the y-grid
  if (ny.eq.128) then
     ! Fuel Jet
     y(1) = 0.0_WP
     do j=2,17
        y(j) = -0.00479441_WP*0.80993658_WP**j+0.00388317_WP
        !y(j) = -0.00255106_WP*0.85800087_WP**j+0.00218883_WP
     end do
     ! Inner Wall
     do j=17,25
        y(j) = inner_wall/8.0_WP*real(j-17,WP) + jet_diameter/2.0_WP
     end do
     ! Pilot
     do j=25,57
        y(j) = 0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP* &
             tanh(2.0_WP*(2.0_WP*real(j-25,WP)/real(32,WP)-1.0_WP))/tanh(2.0_WP) + &
             inner_wall+jet_diameter/2.0_WP + &
             0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP
     end do
     ! Outer Wall
     do j=57,63
        y(j) = outer_wall/6.0_WP*real(j-57,WP) + pilot_diameter/2.0_WP
     end do
     ! Coflow
     do j=63,129
        y(j) = 0.000000588244_WP*1.10501068_WP**j+0.00888257_WP
     end do
  else if (ny.eq.192) then
     ! Fuel Jet
     y(1) = 0.0_WP
     do j=2,25
        y(j) = -0.0044690845_WP*0.8689_WP**j+0.0038831875_WP
     end do
     ! Inner Wall
     do j=25,37
        y(j) = inner_wall/12.0_WP*real(j-25,WP) + jet_diameter/2.0_WP
     end do
     ! Pilot
     do j=37,85
        y(j) = 0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP* &
             tanh(2.0_WP*(2.0_WP*real(j-37,WP)/real(48,WP)-1.0_WP))/tanh(2.0_WP) + &
             inner_wall+jet_diameter/2.0_WP + &
             0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP
     end do
     ! Outer Wall
     do j=85,95
        y(j) = outer_wall/10.0_WP*real(j-85,WP) + pilot_diameter/2.0_WP
     end do
     ! Coflow
     do j=95,193
        y(j) = 0.00000042663327_WP*1.0708_WP**j+0.0089166331_WP
     end do
  else if (ny.eq.208) then
     ! Fuel Jet
     y(1) = 0.0_WP
     do j=2,51
        y(j) = -0.0043848107_WP*0.956_WP**j+0.0041918790_WP
     end do
     ! Inner Wall
     do j=51,63
        y(j) = inner_wall/12.0_WP*real(j-51,WP) + jet_diameter/2.0_WP
     end do
     ! Pilot
     do j=63,117
        y(j) = 0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP* &
             tanh(1.779_WP*(2.0_WP*real(j-63,WP)/real(54,WP)-1.0_WP))/tanh(1.779_WP) + &
             inner_wall+jet_diameter/2.0_WP + &
             0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP
     end do
     ! Outer Wall
     do j=117,127
        y(j) = outer_wall/10.0_WP*real(j-117,WP) + pilot_diameter/2.0_WP
     end do
     ! Coflow
     do j=127,209
        y(j) = 4.9912942E-09_WP*1.08812_WP**j+0.0089729268_WP
     end do
  else if (ny.eq.256) then
     ! Fuel Jet
     y(1) = 0.0_WP
     do j=2,65
        y(j) = -0.0045497033_WP*0.9709_WP**j+0.0044173070_WP
     end do
     ! Inner Wall
     do j=65,77
        y(j) = inner_wall/12.0_WP*real(j-65,WP) + jet_diameter/2.0_WP
     end do
     ! Pilot
     do j=77,157
        y(j) = 0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP* &
             tanh(1.4752_WP*(2.0_WP*real(j-77,WP)/real(80,WP)-1.0_WP))/tanh(1.4752_WP) + &
             inner_wall+jet_diameter/2.0_WP + &
             0.5_WP*(pilot_diameter-jet_diameter-2.0_WP*inner_wall)/2.0_WP
     end do
     ! Outer Wall
     do j=157,167
        y(j) = outer_wall/10.0_WP*real(j-157,WP) + pilot_diameter/2.0_WP
     end do
     ! Coflow
     do j=167,257
        y(j) = 0.00000000082512130_WP*1.07862_WP**j+0.0089455795_WP
     end do
  else
     write(*,'(a)') 'ny must be 192 or 208 or 256'
     call kill
  end if

  do j=1,ny+1
     print*, j, y(j)
  end do

  ! Create the z-grid
  do k=1,nz+1
     z(k) = real(k-1,WP)*twoPi/real(nz,WP)
  end do

  ! Create the mask
  mask = 0
  if ((ny.eq.192).and.(nx.eq.384)) then
     mask(1:16,25:36) = 1
     mask(1:16,85:94) = 1
  else if ((ny.eq.208).and.(nx.eq.400)) then
     mask(1:30,51:62) = 1
     mask(1:30,117:128) = 1
  else if ((ny.eq.256).and.(nx.eq.416)) then
     mask(1:48,65:76) = 1
     mask(1:48,157:166) = 1
  else if ((ny.eq.128).and.(nx.eq.424)) then
     mask(1:48,17:24) = 1
     mask(1:48,57:63) = 1
  else
     write(*,'(a)') '(nx and ny) must be (192 and 384) or (208 and 400) or (256 and 416)'
     call kill
  end if
  return
end subroutine sydneyfine_grid


! =============================== !
! Create the variable array: Data !
! =============================== !
subroutine sydneyfine_data
  use sydney
  implicit none
  
  ! Allocate the array data
  nvar = 11
  allocate(data(nx,ny,nz,nvar))
  allocate(names(nvar))
  
  ! Initialize the variables and names: Gas-phase
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
  data(:,:,:,11) = 0.0_WP;      names(11) = 'PROG'
  
  return
end subroutine sydneyfine_data

