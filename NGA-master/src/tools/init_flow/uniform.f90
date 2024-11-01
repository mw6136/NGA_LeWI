module uniform
  use precision
  use param
  implicit none

  ! Length of the domain
  real(WP) :: Lx,Ly,Lz

  ! Mean velocities
  real(WP) :: Umean,Vmean,Wmean,ZMIXin

  ! Pointers to variables in data
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: T
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: dRHO
  real(WP), dimension(:,:,:), pointer :: ZMIX

end module uniform


subroutine uniform_grid
  use uniform
  use parser
  implicit none

  integer :: i,j,k

  ! Read the size of the domain
  call parser_read('nx',nx)
  call parser_read('Lx',Lx)
  call parser_read('ny',ny)
  call parser_read('Ly',Ly)
  call parser_read('nz',nz)
  call parser_read('Lz',Lz)
  
  ! Set the periodicity
!!$  xper = 1
!!$  yper = 1
!!$  zper = 1
  
  call parser_read('xper',xper)
  call parser_read('yper',yper)
  call parser_read('zper',zper)

  ! Cartesian
  icyl = 0

  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  allocate(mask(nx,ny))

  ! Create the grid
  do i=1,nx+1
     x(i) = real(i-1,WP)*Lx/real(nx,WP)
  end do
  do j=1,ny+1
     y(j) = real(j-1,WP)*Ly/real(ny,WP)
  end do
  do k=1,nz+1
     z(k) = real(k-1,WP)*Lz/real(nz,WP)
  end do

  ! Create the mid points
  do i=1,nx
     xm(i) = 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j) = 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k) = 0.5_WP*(z(k)+z(k+1))
  end do

  ! Create the masks
  mask = 0

  return
end subroutine uniform_grid


subroutine uniform_data
  use uniform
  use parser
  implicit none

  character(str_medium) :: data_type
  real(WP) :: rho_init, T_init, P_init, xg, c, dt, DIFF, s
  real(WP), parameter :: R_cst   = 8314.34_WP   ! J/[kmol K]
  integer :: i,j,k
  real(WP), parameter :: Pi    = 3.1415926535897932385_WP
  real(WP), dimension(:,:), allocatable :: tmp_data
  logical :: nonuniform, discontinuous
  real(WP) :: dx1, dx2
  
  ! Data type
  call parser_read('Data type',data_type,'cold')

  ! Read the means
  call parser_read('Mean U Velocity',Umean)
  call parser_read('Mean V Velocity',Vmean)
  call parser_read('Mean W Velocity',Wmean)
  call parser_read('Initial scalar value',ZMIXin)

  select case(trim(adjustl(data_type)))
  case ('cold')
     ! Allocate the array data
     nvar = 5
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))

     ! Link the pointers
     U => data(:,:,:,1); names(1) = 'U'
     V => data(:,:,:,2); names(2) = 'V'
     W => data(:,:,:,3); names(3) = 'W'
     T => data(:,:,:,4); names(4) = 'T'
     ZMIX => data(:,:,:,5); names(5) = 'ZMIX'

     ! Create them
     U = Umean
     V = Vmean
     W = Wmean
     T = 298.0_WP
     ZMIX = ZMIXin
!!$  do k=1,nz
!!$     do j=1,ny
!!$        do i=1,nx
!!$           ZMIX(i,j,k) = ZMIXin + sin(2*Pi*x(i)/Lx)
!!$        end do
!!$     end do
!!$  end do


  case ('finite chem onestep')
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
     ZMIX => data(:,:,:,12); names(12) = 'ZMIX'
     names(7) = 'N2'
     names(9) = 'O2'
     names(8) = 'H2'
     names(10) = 'H2O'
     names(11) = 'T'
     ! Create them
     call parser_read('Pressure',P_init) ! Thermodynamic pressure to compute density
     call parser_read('Timestep size',dt)
     call parser_read('Initial diffusivity',DIFF)
     U = Umean
     V = Vmean
     W = Wmean
     P = 0.0_WP !Hydrodynamic pressure
     call parser_read('Initial temperature',T_init)
     data(:,:,:,11) = T_init



     !! Some properties of the Gaussian distribution:
     xg = 0.0125_WP !0.00625_WP! - 0.5_WP*dt*Umean
     c  = 0.0015625_WP
!!$     print *, x
!!$
!!$     U(1,:,:) = 0.25_WP
!!$     do i=2,nx
!!$        U(i,:,:) = U(i-1,1,1) + abs(exp(x(i)*(xg-0.5_WP*x(i))/c**2)*(0.25_WP*2.71828_WP**(0.5_WP*(x(i)-xg)**2/c**2)*x(i) - &
!!$             0.25_WP*(0.405578_WP+2.71828_WP**(0.5_WP*(x(i)-xg)**2/c**2))*x(i) - &
!!$             0.25_WP*2.71828_WP**(0.5_WP*(x(i)-xg)**2/c**2)*xg + &
!!$             0.25_WP*(0.405578_WP+2.71828_WP**(0.5_WP*(x(i)-xg)**2/c**2))*xg) / &
!!$             (c**2*(0.405578_WP+exp(0.5_WP*xg**2/c**2)))) / (1500.0_WP*342.0_WP)
!!$     end do
!!$
!!$     print *, U

!!$     do i=1,nx
!!$        U(i,:,:) = 0.25*(0.405578+2.71828**(0.5*(x(i)-xg)**2/c**2))*exp(-0.5*(x(i)-xg)**2/c**2+0.5*xg**2/c**2)/(0.405578+2.71828**(0.5*xg**2/c**2))
!!$     end do

     !! Try this with diffusion correction - THIS WORKS BETTER
     !! N2
!!$     do i=1,nx
!!$        data(i,:,:,7) = (0.9717_WP-0.768_WP) & 
!!$             * ( 1 + 0.5_WP*dt*(DIFF*((xm(i)-xg)**2/c**4 - 1.0_WP/c**2) + Umean*(xm(i)-xg)/c**2)) &
!!$             * exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.768_WP
!!$     end do
!!$     !! H2
!!$     do i=1,nx
!!$        data(i,:,:,8) = 0.0283_WP & 
!!$             * ( 1 + 0.5_WP*dt*(DIFF*((xm(i)-xg)**2/c**4 - 1.0_WP/c**2) + Umean*(xm(i)-xg)/c**2)) &
!!$             * exp(-(xm(i)-xg)**2/(2.0_WP*c**2))
!!$     end do
!!$     !! O2
!!$     do i=1,nx
!!$        data(i,:,:,9) = (0.0_WP-0.232_WP) & 
!!$             * ( 1 + 0.5_WP*dt*(DIFF*((xm(i)-xg)**2/c**4 - 1.0_WP/c**2) + Umean*(xm(i)-xg)/c**2)) &
!!$             * exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.232_WP
!!$     end do


     !! --- First, compute scalar / density field at n-1/2
     !! ---    Use this for imposed density / velocity cases.
     do i=1,nx
        data(i,:,:,7)  = (0.9717_WP-0.768_WP)*exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.768_WP
     end do
     do i=1,nx
        data(i,:,:,9)  = (0.0_WP-0.232_WP)*exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.232_WP
     end do
     do i=1,nx
        data(i,:,:,8) = 0.0283_WP*exp(-(xm(i)-xg)**2/(2.0_WP*c**2))
     end do
     data(:,:,:,10) = 0.0_WP

     !! Sigmoidal flame
     ! N2
!!$     data(:,:,:,7) = 0.74515373
!!$     do i=1,nx
!!$        ! H2
!!$        data(i,:,:,8) = 0.0285173405_WP*0.5_WP*(tanh(2.0_WP*pi*(0.00625_WP-xm(i))/0.00625_WP)+1.0_WP)
!!$     end do
!!$     do i=1,nx
!!$        ! O2
!!$        data(i,:,:,9) = 0.22632883_WP*0.5_WP*(tanh(2.0_WP*pi*(0.00625_WP-xm(i))/0.00625_WP)+1.0_WP)
!!$     end do
!!$     do i=1,nx
!!$        ! H2O
!!$        data(i,:,:,10) = 0.25484627_WP*0.5_WP*(tanh(2.0_WP*pi*(xm(i)-0.00625_WP)/0.00625_WP)+1.0_WP)
!!$     end do
!!$     do i=1,nx
!!$        ! T
!!$        data(i,:,:,11) = (1000.0_WP-T_init)*0.5_WP*(tanh(2.0_WP*pi*(xm(i)-0.00625_WP)/0.00625_WP)+1.0_WP)+T_init
!!$     end do
        

     ! Density field
!!$     RHO = 1.0_WP
     do i=1,nx
        RHO(i,:,:) =   P_init / (data(i,:,:,7)/28.02_WP &
             + data(i,:,:,9)/32.0_WP &
             + data(i,:,:,8)/2.016_WP + data(i,:,:,10)/18.016_WP) / (data(i,:,:,11) * R_cst)
     end do
     dRHO = 0.0_WP





     !! Note: we're not using ZMIXin
     ZMIX = ZMIXin 



     !!!!!!!!!! NONPREMIXED CASES

  case ('finite chem H2')
     ! Allocate the array data
     nvar = 17
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
     ZMIX => data(:,:,:,17); names(17) = 'ZMIX'
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
     call parser_read('Pressure',P_init) ! Thermodynamic pressure to compute density
     call parser_read('Initial temperature',T_init)
     U = Umean
     V = Vmean
     W = Wmean
     P = 0.0_WP !Hydrodynamic pressure
     xg = 0.0125
     c  = 0.0015625


     !! Use the following with detailed H2 chemistry:
     do i=1,nx
        data(i,:,:,7)  = (0.9717_WP-0.768_WP)*exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.768_WP
     end do
!!$     data(:,:,:,7)  = 0.74515373_WP !N2
     data(:,:,:,8)  = 1.0e-30 !0.0_WP !H
     do i=1,nx
        data(i,:,:,9)  = (0.0_WP-0.232_WP)*exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.232_WP
     end do
!!$     data(:,:,:,9)  = 0.22632883_WP !O2
     data(:,:,:,10) = 1.0e-30 !0.0_WP !O
     data(:,:,:,11) = 1.0e-30 !0.0_WP !OH
     do i=1,nx
        data(i,:,:,12) = 0.0283_WP*exp(-(xm(i)-xg)**2/(2.0_WP*c**2))
     end do
!!$     data(:,:,:,12) = 0.02851743_WP !H2

!!$     !! Sigmoidal flame
!!$     ! N2
!!$     data(:,:,:,7) = 0.74515373
!!$     do i=1,nx
!!$        ! H2
!!$        data(i,:,:,12) = 0.02851743_WP*0.5_WP*(tanh(2.0_WP*pi*(-xm(i)+0.01875_WP)/0.00625_WP)+1.0_WP)
!!$     end do
!!$     do i=1,nx
!!$        ! O2
!!$        data(i,:,:,9) = 0.22632883_WP*0.5_WP*(tanh(2.0_WP*pi*(-xm(i)+0.01875_WP)/0.00625_WP)+1.0_WP)
!!$     end do
!!$     do i=1,nx
!!$        ! H2O
!!$        data(i,:,:,13) = 0.254846266_WP*0.5_WP*(tanh(2.0_WP*pi*(-0.01875_WP+xm(i))/0.00625_WP)+1.0_WP)
!!$     end do
     do i=1,nx
        ! T
        data(i,:,:,16) = (1000.0_WP-T_init)*0.5_WP*(tanh(2.0_WP*pi*(-0.0125_WP+xm(i))/0.00625_WP)+1.0_WP)+T_init
     end do

     data(:,:,:,13) = 1.0e-30 !0.0_WP
     data(:,:,:,14) = 1.0e-30 !0.0_WP
     data(:,:,:,15) = 1.0e-30 !0.0_WP

    
!!$     call parser_read('Initial temperature',T_init)
!!$     data(:,:,:,16) = T_init
!!$     RHO  = P_init * 28.84_WP / (T_init * R_cst)
     do i=1,nx
        RHO(i,:,:) = P_init / (data(i,:,:,7)/28.02_WP &
             + data(i,:,:,9)/32.0_WP &
             + data(i,:,:,12)/2.016_WP + data(i,:,:,13)/18.016_WP) / (data(i,:,:,16) * R_cst)
     end do
     

     dRHO = 0.0_WP
     ZMIX = ZMIXin 






  case ('finite chem NDD-MX')
     ! Allocate the array data
     nvar = 39
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
     names(7)  = 'N2'
     names(8)  = 'H'
     names(9)  = 'O2'
     names(10) = 'O'
     names(11) = 'OH'
     names(12) = 'H2'
     names(13) = 'H2O'
     names(14) = 'CO'
     names(15) = 'CO2'
     names(16) = 'HCO'
     names(17) = 'TXXCH2'
     names(18) = 'CH2O'
     names(19) = 'CH3'
     names(20) = 'C2H2'
     names(21) = 'SXXCH2'
     names(22) = 'CH4'
     names(23) = 'C2H4'
     names(24) = 'C2H5'
     names(25) = 'C2H3'
     names(26) = 'C3H6'
     names(27) = 'AXXC3H5'
     names(28) = 'C7H14'
     names(29) = 'PXXC4H9'
     names(30) = 'C5H6'
     names(31) = 'C5H5'
     names(32) = 'A1CH3XH8'!truncated names
     names(33) = 'A1XXC6H6'
     names(34) = 'A1CH3CH3'
     names(35) = 'C9H19'
     names(36) = 'C12H25'
     names(37) = 'NXXC12H2'
     names(38) = 'T'
     ZMIX => data(:,:,:,39); names(39) = 'ZMIX'

     ! Create them
     call parser_read('Pressure',P_init) ! Thermodynamic pressure to compute density
     call parser_read('Initial temperature',T_init)
     U = Umean
     V = Vmean
     W = Wmean
     P = 0.0_WP !Hydrodynamic pressure
     xg = 0.0125
     c  = 0.0015625

     !! Create the species profiles
     do i=1,nx !N2
        data(i,:,:,7)  = (0.075249708510528_WP-0.768_WP)*exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.768_WP
     end do
     data(:,:,:,8)  = 1.0e-30 !0.0_WP !H
     do i=1,nx !O2
        data(i,:,:,9)  = (0.0_WP-0.232_WP)*exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.232_WP
     end do
     data(:,:,:,10) = 1.0e-30 !0.0_WP !O
     data(:,:,:,11) = 1.0e-30 !0.0_WP !OH
     data(:,:,:,12) = 1.0e-30 !0.0_WP !H2
     data(:,:,:,13) = 1.0e-30 !0.0_WP !H2O
     data(:,:,:,14) = 1.0e-30 !0.0_WP !CO
     data(:,:,:,15) = 1.0e-30 !0.0_WP !CO2
     data(:,:,:,16) = 1.0e-30 !0.0_WP !HCO
     data(:,:,:,17) = 1.0e-30 !0.0_WP !TXXCH2
     data(:,:,:,18) = 1.0e-30 !0.0_WP !CH2O
     data(:,:,:,19) = 1.0e-30 !0.0_WP !CH3
     data(:,:,:,20) = 1.0e-30 !0.0_WP !C2H2
     data(:,:,:,21) = 1.0e-30 !0.0_WP !SXXCH2
     data(:,:,:,22) = 1.0e-30 !0.0_WP !CH4
     data(:,:,:,23) = 1.0e-30 !0.0_WP !C2H4
     data(:,:,:,24) = 1.0e-30 !0.0_WP !C2H5
     data(:,:,:,25) = 1.0e-30 !0.0_WP !C2H3
     data(:,:,:,26) = 1.0e-30 !0.0_WP !C3H6
     data(:,:,:,27) = 1.0e-30 !0.0_WP !AXXC3H5
     data(:,:,:,28) = 1.0e-30 !0.0_WP !C7H14
     data(:,:,:,29) = 1.0e-30 !0.0_WP !PXXC4H9
     data(:,:,:,30) = 1.0e-30 !0.0_WP !C5H6
     data(:,:,:,31) = 1.0e-30 !0.0_WP !C5H5
     data(:,:,:,32) = 1.0e-30 !0.0_WP !A1CH3XH8C7
     data(:,:,:,33) = 1.0e-30 !0.0_WP !A1XXC6H6
     do i=1,nx !A1CH3CH3XH10C8 m-xylene
        data(i,:,:,34) = 0.235968450241449_WP*exp(-(xm(i)-xg)**2/(2.0_WP*c**2))
     end do
     data(:,:,:,35) = 1.0e-30 !0.0_WP !C9H19
     data(:,:,:,36) = 1.0e-30 !0.0_WP !C12H25 
     do i=1,nx !NXXC12H26 n-dodecane
        data(i,:,:,37) = 0.688781841248023_WP*exp(-(xm(i)-xg)**2/(2.0_WP*c**2))
     end do
     do i=1,nx !T
        data(i,:,:,38) = (1400.0_WP-T_init)*0.5_WP*(tanh(2.0_WP*pi*(-0.0125_WP+xm(i))/0.00625_WP)+1.0_WP)+T_init
     end do



     do i=1,nx
        RHO(i,:,:) = P_init / &
             ( data(i,:,:,7)/28.02_WP &
             + data(i,:,:,9)/32.0_WP &
             + data(i,:,:,12)/2.016_WP &
             + data(i,:,:,34)/106.16_WP &
             + data(i,:,:,37)/170.34_WP ) / &
             ( data(i,:,:,38) * R_cst )
     end do
     dRHO = 0.0_WP
     ZMIX = ZMIXin 



  case ('finite chem GRI')
     ! Allocate the array data
     nvar = 42
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
     names(7)  = 'N2'
     names(8)  = 'O'
     names(9)  = 'O2'
     names(10) = 'H'
     names(11) = 'OH'
     names(12) = 'H2'
     names(13) = 'HO2'
     names(14) = 'H2O2'
     names(15) = 'CH'
     names(16) = 'CO'
     names(17) = 'CH2'
     names(18) = 'HCO'
     names(19) = 'CH2YXCH2'
     names(20) = 'CH3'
     names(21) = 'CH2O'
     names(22) = 'CH4'
     names(23) = 'CO2'
     names(24) = 'CH2OH'
     names(25) = 'CH3OH'
     names(26) = 'C2H'
     names(27) = 'C2H2'
     names(28) = 'HCCO'
     names(29) = 'C2H3'
     names(30) = 'CH2CO'
     names(31) = 'C2H4'
     names(32) = 'C2H5'
     names(33) = 'C2H6'
     names(34) = 'H2O'
     names(35) = 'AR'
     names(36) = 'C'
     names(37) = 'CH2CHO'
     names(38) = 'CH3CHO'
     names(39) = 'C3H8'
     names(40) = 'C3H7'
     names(41) = 'T'
     ZMIX => data(:,:,:,42); names(42) = 'ZMIX'

     ! Create them
     call parser_read('Pressure',P_init) ! Thermodynamic pressure to compute density
     call parser_read('Initial temperature',T_init)
     U = Umean
     V = Vmean
     W = Wmean
     P = 0.0_WP !Hydrodynamic pressure
     xg = 0.0125
     c  = 0.0015625

     !! Create the species profiles
     do i=1,nx !N2
        data(i,:,:,7)  = (0.428257015345112_WP-0.768_WP)*exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.768_WP
     end do
     data(:,:,:,8)  = 1.0e-30 !0.0_WP !O
     do i=1,nx !O2
        data(i,:,:,9)  = (0.0_WP-0.232_WP)*exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.232_WP
     end do
     data(:,:,:,10) = 1.0e-30 !0.0_WP !H
     data(:,:,:,11) = 1.0e-30 !0.0_WP !OH
     data(:,:,:,12) = 1.0e-30 !0.0_WP !H2
     data(:,:,:,13) = 1.0e-30 !0.0_WP !HO2
     data(:,:,:,14) = 1.0e-30 !0.0_WP !H2O2
     data(:,:,:,15) = 1.0e-30 !0.0_WP !CH
     data(:,:,:,16) = 1.0e-30 !0.0_WP !CO
     data(:,:,:,17) = 1.0e-30 !0.0_WP !CH2
     data(:,:,:,18) = 1.0e-30 !0.0_WP !HCO
     data(:,:,:,19) = 1.0e-30 !0.0_WP !CH2YXCH2
     data(:,:,:,20) = 1.0e-30 !0.0_WP !CH3
     data(:,:,:,21) = 1.0e-30 !0.0_WP !CH2O
     do i=1,nx !CH4
        data(i,:,:,22) = 0.571742984654888_WP*exp(-(xm(i)-xg)**2/(2.0_WP*c**2))
     end do
     data(:,:,:,23) = 1.0e-30 !0.0_WP !CO2
     data(:,:,:,24) = 1.0e-30 !0.0_WP !CH2OH
     data(:,:,:,25) = 1.0e-30 !0.0_WP !CH3OH
     data(:,:,:,26) = 1.0e-30 !0.0_WP !C2H
     data(:,:,:,27) = 1.0e-30 !0.0_WP !C2H2
     data(:,:,:,28) = 1.0e-30 !0.0_WP !HCCO
     data(:,:,:,29) = 1.0e-30 !0.0_WP !C2H3
     data(:,:,:,30) = 1.0e-30 !0.0_WP !CH2CO
     data(:,:,:,31) = 1.0e-30 !0.0_WP !C2H4
     data(:,:,:,32) = 1.0e-30 !0.0_WP !C2H5
     data(:,:,:,33) = 1.0e-30 !0.0_WP !C2H6
     data(:,:,:,34) = 1.0e-30 !0.0_WP !H2O
     data(:,:,:,35) = 1.0e-30 !0.0_WP !AR
     data(:,:,:,36) = 1.0e-30 !0.0_WP !C 
     data(:,:,:,37) = 1.0e-30 !0.0_WP !CH2CHO
     data(:,:,:,38) = 1.0e-30 !0.0_WP !CH3CHO
     data(:,:,:,39) = 1.0e-30 !0.0_WP !C3H8
     data(:,:,:,40) = 1.0e-30 !0.0_WP !C3H7 
     do i=1,nx !T
        data(i,:,:,41) = (1530.0_WP-T_init)*0.5_WP*(tanh(2.0_WP*pi*(-0.0125_WP+xm(i))/0.00625_WP)+1.0_WP)+T_init
     end do
     !1530


     do i=1,nx
        RHO(i,:,:) = P_init / &
             ( data(i,:,:,7)/28.02_WP &
             + data(i,:,:,9)/32.0_WP &
             + data(i,:,:,22)/16.032_WP ) / &
             ( data(i,:,:,41) * R_cst )
     end do
     dRHO = 0.0_WP
     ZMIX = ZMIXin 



     !!!!!!!!!!!!!! FLAMEMASTER CASES

  case ('finite chem H2 prem FM')
     ! Allocate the array data
     nvar = 6+11
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     allocate(tmp_data(nx+1,26))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
     ZMIX => data(:,:,:,17); names(17) = 'ZMIX';
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
     call parser_read('Pressure',P_init) ! Thermodynamic pressure to compute density

     ! Read in from FlameMaster data file
     open(unit=13,file='/home/jmacart/sim/test/scalar_convergence/one_step/temporal_n256/03_rho_test/flamemaster/prem/h2_flame/H2_p01_0phi1_0000tu0300.kg.256.R1.1')

     do i=1,nx+1
!!$        read(13,"(26E15.10E2)") (tmp_data(i,j), j=1,26)
        read(13,*) (tmp_data(i,j), j=1,26)
     end do
     close(13)

     ! Create the grid from FlameMaster data
     do i=1,nx+1
        x(i) = tmp_data(i,1)
     end do
     ! Create the mid points
     do i=1,nx
        xm(i) = 0.5_WP*(x(i)+x(i+1))
     end do
     
     ! Read in profiles
     do i=1,nx
!!$        data(i,:,:, 1) = 0.5_WP*(tmp_data(i, 2)/tmp_data(i,23)+tmp_data(i+1, 2)/tmp_data(i+1,23)) !U
        data(i,:,:, 1) = tmp_data(i, 2)/tmp_data(i,23) !U
        data(i,:,:, 5) = 0.5_WP*(tmp_data(i,23)+tmp_data(i+1,23)) !RHO
        data(i,:,:, 7) = 0.5_WP*(tmp_data(i, 4)+tmp_data(i+1, 4)) !N2
        data(i,:,:, 8) = 0.5_WP*(tmp_data(i, 5)+tmp_data(i+1, 5)) !H
        data(i,:,:, 9) = 0.5_WP*(tmp_data(i, 6)+tmp_data(i+1, 6)) !O2
        data(i,:,:,10) = 0.5_WP*(tmp_data(i, 7)+tmp_data(i+1, 7)) !O
        data(i,:,:,11) = 0.5_WP*(tmp_data(i, 8)+tmp_data(i+1, 8)) !OH
        data(i,:,:,12) = 0.5_WP*(tmp_data(i, 9)+tmp_data(i+1, 9)) !H2
        data(i,:,:,13) = 0.5_WP*(tmp_data(i,10)+tmp_data(i+1,10)) !H2O
        data(i,:,:,14) = 0.5_WP*(tmp_data(i,11)+tmp_data(i+1,11)) !HO2
        data(i,:,:,15) = 0.5_WP*(tmp_data(i,12)+tmp_data(i+1,12)) !H2O2
        data(i,:,:,16) = 0.5_WP*(tmp_data(i, 3)+tmp_data(i+1, 3)) !T
     end do

     V = Vmean
     W = Wmean
     !P = 0.0_WP !Hydrodynamic pressure
     dRHO = 0.0_WP
     ZMIX = ZMIXin 

     deallocate(tmp_data)



  case ('finite chem NDD-MX prem FM')
     ! Allocate the array data
     nvar = 39
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     allocate(tmp_data(nx+1,70))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
     names(7)  = 'N2'
     names(8)  = 'H'
     names(9)  = 'O2'
     names(10) = 'O'
     names(11) = 'OH'
     names(12) = 'H2'
     names(13) = 'H2O'
     names(14) = 'CO'
     names(15) = 'CO2'
     names(16) = 'HCO'
     names(17) = 'TXXCH2'
     names(18) = 'CH2O'
     names(19) = 'CH3'
     names(20) = 'C2H2'
     names(21) = 'SXXCH2'
     names(22) = 'CH4'
     names(23) = 'C2H4'
     names(24) = 'C2H5'
     names(25) = 'C2H3'
     names(26) = 'C3H6'
     names(27) = 'AXXC3H5'
     names(28) = 'C7H14'
     names(29) = 'PXXC4H9'
     names(30) = 'C5H6'
     names(31) = 'C5H5'
     names(32) = 'A1CH3XH8'!truncated names
     names(33) = 'A1XXC6H6'
     names(34) = 'A1CH3CH3'
     names(35) = 'C9H19'
     names(36) = 'C12H25'
     names(37) = 'NXXC12H2'
     names(38) = 'T'
     ZMIX => data(:,:,:,39); names(39) = 'ZMIX'

     ! Create them
     call parser_read('Pressure',P_init) ! Thermodynamic pressure to compute density
     V = Vmean
     W = Wmean
     P = 0.0_WP !Hydrodynamic pressure

     ! Read in from FlameMaster data file
     open(unit=13,file='/home/jmacart/flamelets/prem/unstretched/NDD-MX/Output_NDD-MX/NX-C12H2_p01_0phi1_0000tu0300-R1.1-n256.kg')

     do i=1,nx+1
        read(13,*) (tmp_data(i,j), j=1,70)
     end do
     close(13)

     ! Create the grid from FlameMaster data
     do i=1,nx+1
        x(i) = tmp_data(i,1)
     end do
     ! Create the mid points
     do i=1,nx
        xm(i) = 0.5_WP*(x(i)+x(i+1))
     end do

     !! Create the species profiles
     do i=1,nx !N2
        data(i,:,:,1)  = tmp_data(i, 2)/tmp_data(i,67) !U
        data(i,:,:,5)  = 0.5_WP*(tmp_data(i,67)+tmp_data(i+1,67)) !RHO
        do j=4,34
           data(i,:,:,j+3) = 0.5_WP*(tmp_data(i,j)+tmp_data(i+1,j)) !Species (see above)
        end do
        data(i,:,:,38) = 0.5_WP*(tmp_data(i, 3)+tmp_data(i+1, 3)) !T
     end do

     dRHO = 0.0_WP
     ZMIX = ZMIXin 





  case ('finite chem GRI prem FM')
     ! Allocate the array data
     nvar = 42
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
     allocate(tmp_data(nx+1,76))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
     names(7)  = 'N2'
     names(8)  = 'O'
     names(9)  = 'O2'
     names(10) = 'H'
     names(11) = 'OH'
     names(12) = 'H2'
     names(13) = 'HO2'
     names(14) = 'H2O2'
     names(15) = 'CH'
     names(16) = 'CO'
     names(17) = 'CH2'
     names(18) = 'HCO'
     names(19) = 'CH2YXCH2'
     names(20) = 'CH3'
     names(21) = 'CH2O'
     names(22) = 'CH4'
     names(23) = 'CO2'
     names(24) = 'CH2OH'
     names(25) = 'CH3OH'
     names(26) = 'C2H'
     names(27) = 'C2H2'
     names(28) = 'HCCO'
     names(29) = 'C2H3'
     names(30) = 'CH2CO'
     names(31) = 'C2H4'
     names(32) = 'C2H5'
     names(33) = 'C2H6'
     names(34) = 'H2O'
     names(35) = 'AR'
     names(36) = 'C'
     names(37) = 'CH2CHO'
     names(38) = 'CH3CHO'
     names(39) = 'C3H8'
     names(40) = 'C3H7'
     names(41) = 'T'
     ZMIX => data(:,:,:,42); names(42) = 'ZMIX'

     ! Create them
     call parser_read('Pressure',P_init) ! Thermodynamic pressure to compute density
     !U = Umean
     V = Vmean
     W = Wmean
     P = 0.0_WP !Hydrodynamic pressure

     ! Read in from FlameMaster data file
     open(unit=13,file='/home/jmacart/flamelets/prem/unstretched/CH4/Output_CH4/CH4_p01_0phi1_0000tu0300-R1.1-n256.kg')

     do i=1,nx+1
        read(13,*) (tmp_data(i,j), j=1,76)
     end do
     close(13)

     ! Create the grid from FlameMaster data
     do i=1,nx+1
        x(i) = tmp_data(i,1)
     end do
     ! Create the mid points
     do i=1,nx
        xm(i) = 0.5_WP*(x(i)+x(i+1))
     end do

     !! Create the species profiles
     do i=1,nx !N2
        data(i,:,:,1)  = tmp_data(i, 2)/tmp_data(i,73) !U
        data(i,:,:,5)  = 0.5_WP*(tmp_data(i,73)+tmp_data(i+1,73)) !RHO
        do j=4,37
           data(i,:,:,j+3) = 0.5_WP*(tmp_data(i,j)+tmp_data(i+1,j)) !Species (see above)
        end do
        data(i,:,:,41) = 0.5_WP*(tmp_data(i, 3)+tmp_data(i+1, 3)) !T
     end do

     dRHO = 0.0_WP
     ZMIX = ZMIXin 







     !! NONUNIFORM GRID
case ('cold nonuniform')
     ! Allocate the array data
     nvar = 4
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     ZMIX => data(:,:,:,4); names(4) = 'ZMIX'
     ! Create them
     U = Umean
     V = Vmean
     W = Wmean
     T = 298.0_WP
     xg = 0.0125
     c  = 0.0015625

     ! Stretch the grid
     call parser_read('Nonuniform',nonuniform,.false.)
     call parser_read('Nonsmooth',discontinuous,.false.)
     if (nonuniform) then
        call parser_read('Stretching',s)
        do i=1,nx+1
           x(i) = Lx*tanh(s*real(i-1,WP)/real(nx,WP))/tanh(s)
        end do
     elseif (discontinuous) then
        dx1 = 0.5_WP*Lx/real(nx,WP)
        dx2 = 1.5_WP*Lx/real(nx,WP)
        !x(1) = 0
        !x(2) = x(1) + dx2
        do i=1,nx/2
           x(2*i-1) = real(2*(i-1))*Lx/real(nx,WP)
           x(2*i)   = x(2*i-1) + dx1
        end do
        x(nx+1) = Lx
     end if
     if (nonuniform.or.discontinuous) then
        ! Create the mid points
        do i=1,nx
           xm(i) = 0.5_WP*(x(i)+x(i+1))
        end do
     end if
     ! Print out the grid
     do i=1,nx+1
        print *, i, x(i)
     end do

!!$     do i=1,nx
!!$        data(i,:,:,4)  = (1.0_WP-0.232_WP)*exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.232_WP
!!$     end do


     ! For MMS
     data(:,:,:,4) = 0.5_WP


     !! OLD STUFF......................


  case ('finite chem H2 ign')
     ! Allocate the array data
     nvar = 6+11
     allocate(data(nx,ny,nz,nvar))
     allocate(names(nvar))
  
     ! Link the pointers
     U    => data(:,:,:,1); names(1) = 'U'
     V    => data(:,:,:,2); names(2) = 'V'
     W    => data(:,:,:,3); names(3) = 'W'
     P    => data(:,:,:,4); names(4) = 'P'
     RHO  => data(:,:,:,5); names(5) = 'RHO'
     dRHO => data(:,:,:,6); names(6) = 'dRHO'
     ZMIX => data(:,:,:,17); names(17) = 'ZMIX';

     !! Some properties of the Gaussian distribution:
     xg = 0.00625_WP
     c  = 0.0015625_WP
     do i=1,nx !N2
        data(i,:,:,7)  = (0.9717_WP-0.768_WP)*exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.768_WP
     end do
     data(:,:,:,8) = 0.0_WP !H
     do i=1,nx !O2
        data(i,:,:,9)  = (0.0_WP-0.232_WP)*exp(-(xm(i)-xg)**2/(2.0_WP*c**2)) + 0.232_WP
     end do
     data(:,:,:,10) = 0.0_WP !O
     data(:,:,:,11) = 0.0_WP !OH
     do i=1,nx !H2
        data(i,:,:,12) = 0.0283_WP*exp(-(xm(i)-xg)**2/(2.0_WP*c**2))
     end do
     data(:,:,:,13) = 0.0_WP !H2O
     data(:,:,:,14) = 0.0_WP !HO2
     data(:,:,:,15) = 0.0_WP !H2O2

     call parser_read('Initial temperature',T_init)
     data(:,:,:,16) = T_init

     call parser_read('Pressure',P_init) ! Thermodynamic pressure to compute density
     do i=1,nx
        RHO(i,:,:) =   P_init / (data(i,:,:,7)/28.02_WP &
             + data(i,:,:,9)/32.0_WP &
             + data(i,:,:,12)/2.016_WP) / (T_init * R_cst)
     end do
     U = Umean
     V = Vmean
     W = Wmean
     P = 0.0_WP !Hydrodynamic pressure
     dRHO = 0.0_WP
     ZMIX = 1.0_WP


  case default
     print*, "Data type not recognized..."     
  end select

  return
end subroutine uniform_data
