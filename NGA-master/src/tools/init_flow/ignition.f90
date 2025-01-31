module ignition
  use precision
  use param
  use parser
  implicit none
  include 'fftw3.f'
  
! This module creates a square box with the required dimensions 
! such that all the boundaries are periodic.Scalars are initialized 
! using info from a mechanism file.
! Mixture fraction and temperature fields can be specified as 
! uniform or fluctuating
  
  ! Length of the domain
  real(WP) :: L
  real(WP) :: dx
  character(str_medium) :: geometry
  
  ! Pointers to variable in data
  ! Velocity. Pressure, Density, Temperature, Mixture Fraction at every gp
  real(WP), dimension(:,:,:), pointer :: U
  real(WP), dimension(:,:,:), pointer :: V
  real(WP), dimension(:,:,:), pointer :: W
  real(WP), dimension(:,:,:), pointer :: P
  real(WP), dimension(:,:,:), pointer :: RHO
  real(WP), dimension(:,:,:), pointer :: T
  real(WP), dimension(:,:,:), pointer :: Z_mix

  ! Average molecular weight of the mixture at every gp
  real(WP), dimension(:,:,:), pointer :: W_mixture
  ! Mixture fraction sum at every gp
  real(WP), dimension(:,:,:), pointer :: sum_Y
  ! Equivalence ratio of the mixture
  real(WP), dimension(:,:,:), pointer :: equiv
  ! Total Number of species involved in the reaction mechanism              
  integer :: Nspecies
  ! Number of non steady species involved in the reaction mechanism   
  integer :: N_nons
  ! Index for mixture fraction, temp, scalar
  integer :: isc_Z, isc_T, isc_sc
  ! Complex buffer
  complex(WP), dimension(:,:,:), pointer :: Cbuf
  ! Real buffer
  real(WP),    dimension(:,:,:), pointer :: Rbuf, Rbuf_temp
  ! Scalar buffer
  real(WP),    dimension(:,:,:), pointer :: tmp_sc
  real(WP)  :: tmp_mean, tmp_rms
  real(WP)  :: min_val, max_val, rms_val, mean_val
end module ignition

! ==================== !
! Create the grid/mesh !
! ====================! 
subroutine ignition_grid
  use ignition
  implicit none
  
  integer :: i,j,k
  logical :: is_p

  ! Read in the size of the domain
  call parser_read('nx',nx)
  call parser_read('Domain size',L)
  call parser_read('Geometry',geometry)  ! Check if 1D, 2D or 3D

  write (*,*) ' '
  write (*,*) 'Number of grid pts = ',  nx
  write (*,*) 'Domain size = ',  L
  write (*,*) 'Geometry = ', trim(geometry)

  ! Set all directions to periodic (make values zero for non periodic)
  call parser_is_defined('xper',is_p)
  if (is_p) then
     call parser_read('xper',xper)
     write (*,*) 'xper : ', xper
  else
     xper = 1
  end if
  yper = 1
  zper = 1  

  ! Choose grid based on the geometry
  geometry=trim(geometry)
  select case (geometry)
  case ('1D')  
     ny = 1
     nz = 1
  case ('2D')  ! Box with nx=ny, nz=1
     ny = nx
     nz = 1
  case ('3D')  ! Box with nx=ny=nz
     ny = nx
     nz = ny
  case default
     write (*,*) 'Enter valid geometry type'
  end select

  ! Set to Cartesian grid
  icyl = 0
  
  ! Allocate the arrays
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(xm(nx),ym(ny),zm(nz))
  ! Controls whether there is a wall
  allocate(mask(nx,ny))
  
  ! Create the grid
  dx = L/nx
  do i=1,nx+1
     x(i) = (i-1)*L/nx
  end do
  do j=1,ny+1
     if (geometry.eq.'1D') then
	! so Ly is thin in proportion to the ny/nx ratio (1D case)
        ! 10 chosen arbritarily
        y(j) = (j-1)*(L/10)
     else
        y(j) = (j-1)*L/ny
     end if
  end do
  do k=1,nz+1
     if (geometry.eq.'3D') then
        z(k) = (k-1)*L/nz			
     else
        if (geometry.eq.'1D') then
           ! 10 chosen arbritarily
           z(k) = (k-1)*(L/10) 
       else
          ! so Lz is thin in proportion to the nz/nx ratio (2D case)
          z(k) = (k-1)*(L/nx)/nz
       end if
    end if
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
  ! Zero needed for periodic case. Set them to one for walls in domain
  mask = 0
  
  return
end subroutine ignition_grid

! ========================= !
! Create the variable array !
! ========================= !
subroutine ignition_data
  use ignition
  implicit none

  character (20), dimension(:), pointer :: NameSpecies   ! Species Names
  ! Check if fuel is defined explicitly
  logical  :: isFueldefined
  ! Molar mass of the different species in kg/kmol
  real(WP), dimension(:), pointer :: W_species
  ! Counters
  integer  :: i, j, k, n
  ! Indexes to relate species number to specific fuels, oxidizer, scalar
  integer  :: isc_fuel, isc_ox_1, isc_ox_2
  ! Initial thermodynamic pressure in the domain
  real(WP) :: P_init
  real(WP), parameter :: R_cst   = 8314.34_WP   ! J/[kmol K]
  ! Get oxidizer and fuel name
  character (len=str_short) ::  fuel,oxidizer,tmp_val
  ! Mean W_mixture over the entire domain  
  real(WP) :: W_mixture_m
  real(WP) :: mean_RHO
  ! molar ratio of n2 to o2 in air
  real(WP) :: air_ratio  = 3.76_WP
  ! Stoichometric mass fraction
  real(WP) :: Z_st
  ! Number of moles of fuel per mole of oxygen in oxidizer
  real(WP) :: moles_fuel
  ! Check if file is present
  logical  :: isdefined, isSpecies
  ! Mass fraction of species at 0,1 Z
  real(WP), dimension(:), pointer :: Z0_species, Z1_species
  real(WP) :: temp_val
  ! Properties of the species
  real(WP), dimension(:), pointer :: Cp_sp, h_sp ! Cp [J/(mol.K)]
  
  ! Get the total number of species from mechanism.f
  call GETNSPECIES(Nspecies)
  ! Get the total number of non steady state species from mechanism.f
  call GETNSPECS(N_nons)
  if (Nspecies.eq.0) then
     write (*,*) 'Please check your mechanism file. No variables were read in'
  end if
  
  ! Allocate the array data
  
  ! Additional variables are 
  ! Vel(3), Temp, Density, dRHO, Pressure, Mixture fraction
  ! Only non steady state species are tranported
  nvar = N_nons+8
  
  allocate(NameSpecies(Nspecies))     ! This contains the species name
  allocate(W_species(Nspecies))       ! This contains species molecular weight
  allocate(data(nx,ny,nz,nvar))       ! This is the main array used by code
  allocate(names(nvar))               ! This contains the name of the scalars
  allocate(W_mixture(nx,ny,nz))       ! Average molecular weight at each gp
  allocate(sum_Y(nx,ny,nz))           ! Mixture fraction sum at each gp
  allocate(equiv(nx,ny,nz))           ! Equivalence ratio at each gp
  allocate(Z0_Species(N_nons)) ! This contains the Y at Z=0
  allocate(Z1_Species(N_nons)) ! This contains the Y at Z=1
  allocate(tmp_sc(nx,ny,nz))          ! Temp scalar
  allocate(Cp_sp(N_nons))             ! Cp of species
  allocate(h_sp(N_nons))              ! h of species
  ! Initialize to zero
  names = ' '
  data = 0.0_WP
  Z0_Species = 0.0_WP
  Z1_Species = 0.0_WP
  isc_fuel=0
  isc_ox_1=0 
  isc_ox_2=0
  sum_Y = 0.0_WP

  ! The seven additional variables are registered with pointers
  ! Names are stored for future reference
  U    => data(:,:,:,1);  names(1)  = 'U'
  V    => data(:,:,:,2);  names(2)  = 'V'
  W    => data(:,:,:,3);  names(3)  = 'W'
  P    => data(:,:,:,4);  names(4)  = 'P'
  RHO  => data(:,:,:,5);  names(5)  = 'RHO'
  names(6)  = 'dRHO'
  isc_sc = 7;
  isc_T=N_nons+isc_sc; T => data(:,:,:,isc_T); names(isc_T) = 'T'
  isc_Z=N_nons+isc_sc+1; Z_mix => data(:,:,:,isc_Z); names(isc_Z) = 'ZMIX'

  ! Get all the species names from the mechanism.f file.
  ! Later only the non steady state species will be important
  call GETSPECIESNAMES(NameSpecies)
  ! Get all the species molar mass from the mechanism.f file. 
  ! Later only the non steady state species will be important
  call GETMOLARMASS(W_species)
  write (*,*) '------------------------------------------'
  write (*,*) 'Using the Detailed Ignition Chemistry case'
  write (*,*) 'Total Species: ', Nspecies
  write (*,*) 'Total Non-Steady Species: ', N_nons
  write (*,*) '------------------------------------------'
  ! Initialize velocty field (either uniform or fluctuating)
  write (*,*) 'Initializing Velocity'
  call ignition_velocity

  !==============================================================!
  ! Initialization options
  ! a) Fuel, oxidizer (air) at present, equivalence ratio(uniform)
  ! b) Mixture fraction (uniform, fluctuating)
  ! c) Individual species (uniform, fluctuating)
  !==============================================================!

  !==============================================================!
  !a)
  ! Read in fuel and oxidizer if they are defined
  ! Check if fuel is defined explicitly
  call parser_is_defined('Fuel',isFueldefined)
  if (isFueldefined) then
     write (*,*) '------------------------------------------'
     write (*,*) 'Initializing with fuel, equivalence ratio:'
     ! Reads in species to be treated as fuel
     call parser_read('Fuel',fuel)
     ! Reads in the oxidizer kind. (allows only air presently)
     call parser_read('Oxidizer',oxidizer)
     ! Link species names to scalars
     do i = 1, N_nons
        ! The first five names are defined above already
        names(i+isc_sc-1)= NameSpecies(i)
        if (NameSpecies(i).eq.trim(fuel)) then
           ! Set index to i for fuel species
           isc_fuel=i
           write (*,*) 'Fuel : ', NameSpecies(i)
        end if
        if (oxidizer.eq.'air') then
           if (NameSpecies(i).eq.'O2') then
              isc_ox_1=i  ! Set index to i for O2
           end if
           if (NameSpecies(i).eq.'N2') then
              isc_ox_2=i  ! Set index to i for N2
           end if
        else
           write (*,*) 'Code presently only handles air as oxidizer. Please enter that'
        end if
     end do
     
     ! Check if the species were found
     do i = 1, N_nons
        if (isc_fuel.eq.0) then
           write (*,*) 'Fuel not found. Please check the fuel name'
           stop
        end if
        if (oxidizer.eq.'air') then
           if (isc_ox_1.eq.0) then
              write (*,*) 'Oxidizer species O2 not found'
           end if
           if (isc_ox_2.eq.0) then
              write (*,*) 'Oxidizer species N2 not found.'
           end if
        end if
     end do

     ! Initialize field
     call ignition_scalar('Equiv')
     equiv=tmp_sc
     write(*,*) 'Equiv mean  : ',tmp_mean
     write(*,*) 'Equiv rms : ',tmp_rms
     tmp_sc = 0.0_WP
     tmp_mean = 0.0_WP
     tmp_rms = 0.0_WP
     
     ! Read in stoichiometric mixture fraction
     call parser_read('St n_fuel per O2', moles_fuel)
     Z_st=(W_species(isc_fuel)*moles_fuel)/&
          (W_species(isc_ox_1)+air_ratio*W_species(isc_ox_2)&
          +(W_species(isc_fuel)*moles_fuel))

     ! Calculate mixture fraction from equivalence ratio
     Z_mix=1.0_WP/(((1.0_WP-Z_st)/(equiv*Z_st))+1.0_WP)

     if (oxidizer.eq.'air') then
        do i=1,nx
           do j=1,ny
              do k=1,nz
                 ! Set mass fraction of fuel
                 data(i,j,k,isc_fuel+isc_sc-1)=(W_species(isc_fuel)*equiv(i,j,k)*moles_fuel)/&
                      (W_species(isc_ox_1)+air_ratio*W_species(isc_ox_2)+&
                      (W_species(isc_fuel)*equiv(i,j,k)*moles_fuel))
                 ! Set mass fraction of O2
                 data(i,j,k,isc_ox_1+isc_sc-1)=W_species(isc_ox_1)/&
                      (W_species(isc_ox_1)+air_ratio*W_species(isc_ox_2)+&
                      (W_species(isc_fuel)*equiv(i,j,k)*moles_fuel))
                 ! Set mass fraction of N2
                 data(i,j,k,isc_ox_2+isc_sc-1)=air_ratio*W_species(isc_ox_2)/&
                      (W_species(isc_ox_1)+air_ratio*W_species(isc_ox_2)+&
                      (W_species(isc_fuel)*equiv(i,j,k)*moles_fuel))
              end do
           end do
        end do
     end if
  else
     call parser_is_defined('Z mean',isdefined)
     if (isdefined) then
        !===========================
        !(b)
        ! Initialize mixture fraction    
        write (*,*) '------------------------------------------'
        write (*,*) 'Initializing with Mixture fraction :'
        do i =1, N_nons
           call parser_is_defined('Z0_'//trim(NameSpecies(i)),isdefined)
           if (isdefined) then
              call parser_read('Z0_'//trim(NameSpecies(i)), Z0_species(i))
              write (*,*) 'Z0 : '//trim(NameSpecies(i)),' ', Z0_species(i)
           end if
           call parser_is_defined('Z1_'//trim(NameSpecies(i)), isdefined)
           if (isdefined) then
              call parser_read('Z1_'//trim(NameSpecies(i)), Z1_species(i))
              write (*,*) 'Z1 : '//trim(NameSpecies(i)),' ', Z1_species(i)
           end if
        end do
        
        if (sum(Z0_species).ne.1.0_WP) then
           write (*,*) 'Mass fraction sum at Z = 0 is not equal to one. Please check'
        end if
        if (sum(Z1_species).ne.1.0_WP) then
           write (*,*) 'Mass fraction sum at Z = 1 is not equal to one. Please check'
        end if
        
        ! Initialize field
        call ignition_scalar('Z')
        Z_mix=tmp_sc
        write(*,*) 'Z mean : ',tmp_mean
        write(*,*) 'Z rms  : ',tmp_rms
        write(*,*) 'Z min  : ',minval(Z_mix)
        write(*,*) 'Z max  : ',maxval(Z_mix)
        tmp_sc = 0.0_WP
        tmp_mean = 0.0_WP
        tmp_rms = 0.0_WP
        
        do i=1, N_nons
           names(i+isc_sc-1)= trim(NameSpecies(i))
           data(:,:,:,isc_sc+i-1)=z0_species(i)*(1.0_WP-Z_mix)+z1_species(i)*Z_mix
        end do
        data(:,:,:,isc_Z)=Z_mix
     else
        !=====================
        !(c) Inidividual species
        write (*,*) '------------------------------------------'
        write (*,*) 'Initializing with individual species :'
        isSpecies=.false.
        do i = 1, N_nons
           names(i+isc_sc-1)= trim(NameSpecies(i))
           ! The first five names are defined above already
           ! Check if this scalar is defined
           call parser_is_defined(trim(names(i+isc_sc-1)) // ' mean' ,isdefined)
           if (isdefined) then
              isSpecies=.true.
             call ignition_scalar(names(i+isc_sc-1))
              data(:,:,:,i+isc_sc-1)=tmp_sc
              write(*,*) trim(names(i+isc_sc-1)), ' mean  : ',tmp_mean
              write(*,*) trim(names(i+isc_sc-1)), ' rms   : ',tmp_rms
              write(*,*)
              tmp_sc = 0.0_WP
              tmp_mean = 0.0_WP
              tmp_rms = 0.0_WP
           end if
        end do
        if (.not.isSpecies) then
           write (*,*) '------------------------------------------'
           write(*,*) 'NO SPECIES WERE INITIALIZED. CHECK!!'
        end if
     end if
  end if

  ! Clip mass fractions between 0 and 1
  ! Calculate sum of mass fractions at a point
  data(:,:,:,isc_sc:isc_sc+N_nons-1) =  min(max(data&
       (:,:,:,isc_sc:isc_sc+N_nons-1),0.0_WP),1.0_WP)
  do n = 1,N_nons
     sum_Y = sum_Y + data(:,:,:,n+isc_sc-1)
  end do

  ! Rescale the mass fractions to 1.0_WP
  if ((maxval(sum_Y).ne.1.0_WP).or.(minval(sum_Y).ne.1.0_WP)) then
     write (*,*) 'Rescaling mass fractions to make sum 1.0'
     do n = 1,N_nons
        data(:,:,:,n+isc_sc-1)=data(:,:,:,n+isc_sc-1)/sum_Y
     end do
  end if
  sum_Y = 0.0_WP
  do n = 1,N_nons
     sum_Y = sum_Y + data(:,:,:,n+isc_sc-1)
  end do

  ! Initialize temperature field (see subroutine for options)
  call ignition_temp

  ! Obtain average molecular weight at every grid point
  W_mixture=0.0_WP
  do n = 1,N_nons
     W_mixture=real(data(:,:,:,n+isc_sc-1)/W_species(n),WP)+W_mixture
  end do
  ! average molecular weight of the mixture in kg/kmol at every grid point

  W_mixture = sum_Y/W_mixture
  W_mixture_m = 0.0_WP
  do i=1,nx
     do j=1,ny
        do k=1,nz
           W_mixture_m=W_mixture(i,j,k)+W_mixture_m
        end do
     end do
  end do
  W_mixture_m = W_mixture_m/real(nx*ny*nz,WP) ! Average molecular weight 
  write (*,*) '------------------------------------------'
  write (*,*) 'W mixture mean : ', W_mixture_m
  ! Read initial thermodynamic pressure
  call parser_read('Pressure',P_init)

  ! Initialize the density field (Ideal gas assumption)
  do i=1,nx
     do j=1,ny
        do k=1,nz
           RHO(i,j,k) = P_init*W_mixture(i,j,k)/(data(i,j,k,isc_T)*R_cst)
        end do
     end do
  end do
  ! Get mean RHO
  mean_RHO = sum(RHO)/(nx*ny*nz)

  ! Rescale velocity to ensure continuity
  do i=1,nx
     do j=1,ny
        do k=1,nz
           U(i,j,k)  = U(i,j,k) * mean_RHO / RHO(i,j,k)
           V(i,j,k)  = V(i,j,k) * mean_RHO / RHO(i,j,k)
           W(i,j,k)  = W(i,j,k) * mean_RHO / RHO(i,j,k)
        end do
     end do
  end do

  ! Output the stats from the initialization if vel fluctuations are there
  call parser_read('Fluctuations',temp_val)
  if (temp_val.gt.0.0_WP)  call ignition_stat  

  ! Get pdfs
  tmp_val = 'Z'
  call get_1dpdf(Z_mix,41,tmp_val)
  tmp_val = 'T'
  call get_1dpdf(data(:,:,:,isc_T),41,tmp_val)

  return
end subroutine ignition_data


! ====================================== !
! Initialize the velocity field          !
! ====================================== !
subroutine ignition_velocity
  use ignition
  implicit none
    
  real(WP) :: Ut                        ! Velocity fluctuation
  real(WP) :: Uval, Wval, Vval          ! Mean u, v, w values
  character(len=str_short) :: spectrum  ! Spectrum type
  real(WP) :: le,ld,epsilon             ! Length scale, epsilon for the spectrum

  call parser_read('U',Uval)
  call parser_read('V',Vval)
  call parser_read('W',Wval)  
  call parser_read('Fluctuations',Ut)
  write (*,*) 'Velocity Fluctuation  :', Ut
  if (Ut.eq.0) then
     ! Uniform
     U = Uval
     V = Vval
     W = Wval
  else
     ld = 0.0_WP
     epsilon = 0.0_WP
     ! Velocity initialization for 2-D and 3-D case
     call parser_read('Spectrum form',spectrum)
     call parser_read('Energetic scale',le)
     if (trim(spectrum).eq.'VKP') then
        call parser_read('Dissipative scale',ld)
        call parser_read('Dissipation',epsilon)
     end if
     
     if (geometry.eq.'1D') then
        U = Uval
        V = Vval
        W = Wval
     else           
        ! Get velocity field assuming mean = 0
        call ignition_spectrum(spectrum, Ut, le, ld, epsilon,'V','velocity')
        ! add background velocity field to the fluctuations
        U = U + Uval
        V = V + Vval
        W = W + Wval
     end if
  end if
end subroutine ignition_velocity

! ================================================================= !
! Initialize the temperature field (Including temperature)          !
! ================================================================= !
! Field can  be initiailized for                                    !
! a) Temperature of all species (uniform or fluctuating)            !
! b) Temperature of individual species (uniform or fluctuating)     !
!    Temperature at gp is then calculated based on mixture enthalpy !
! ================================================================= !
subroutine ignition_temp
  use ignition
  implicit none
  
  ! Scalar name
  character(len=str_medium) :: name_s, name_s_t
  ! Temperory values
  real(WP) :: tmp_hsp_test, h_slope, old_t,old_h, tmp1
  ! Counters
  integer  :: n, i, j, k, iterate
  ! Check vars
  logical  :: isdefined, check_sp, check_sp_t, check_it
  ! temp array for temp calculation
  real(WP), dimension(:,:,:), pointer     :: tmp_hsp, tmp_T_t
  real(WP), dimension(:,:,:,:), pointer   :: tmp_T
  ! Array for enthalpy of species if needed
  real(WP), dimension(N_nons)    :: h_sp_tmp, Cp_sp_tmp

  rms_val=0.0_WP
  iterate=0
  tmp1=0.0_WP
  allocate(tmp_hsp(nx,ny,nz))
  allocate(tmp_T_t(nx,ny,nz))
  allocate(tmp_T(nx,ny,nz,N_nons))
  tmp_hsp = 0.0_WP
  tmp_T = 0.0_WP
  tmp_T_t = 0.0_WP
  ! Read name of scalar
  name_s = trim(names(isc_T))
  if (name_s.ne.'T') then
     write (*,*) '------------------------------------------'
     write (*,*) 'T IN SIMULATION NOT INITIALIZED CORRECTLY!'
  else
     write (*,*) '------------------------------------------'
     write (*,*) 'Initializing Temperature Field'   
  end if

  ! check if overall T is defined
  call parser_is_defined('T mean',isdefined)
  if (isdefined) then
     write(*,*) 'Overall T field defined'
     call ignition_scalar('T')
     data(:,:,:,isc_T) = tmp_sc
     write(*,*) 'T mean  : ',tmp_mean
     write(*,*) 'T rms   : ',tmp_rms
     write(*,*) 'T min  : ', minval(data(:,:,:,isc_T))
     write(*,*) 'T max   : ',maxval(data(:,:,:,isc_T))
!!$     do i=1,nx
!!$        data(i,:,:,isc_T) = tmp_mean + 3.0e2_WP*sin(2.0_WP*acos(-1.0_WP)*xm(i)/L)
!!$     end do
     tmp_sc = 0.0_WP
     tmp_mean = 0.0_WP
     tmp_rms = 0.0_WP
     n = nx/16
     data(n:nx,:,:,isc_T) = 2000.0_WP
  else
     data(:,:,:,isc_T) = 300.0_WP
  end if

  return
end subroutine ignition_temp

! ===================================================== !
! Initialize the scalar fields                          !
! ===================================================== !
subroutine ignition_scalar(name_s)
  use ignition
  implicit none
  
  ! Spectrum type
  character(len=str_short) :: spectrum
  ! Scalar name
  character(len=*):: name_s
  ! Length scales, epsilon
  real(WP) :: le,ld, epsilon
  ! ks/ko, kc/ks ratios
  real(WP) :: ksk0ratio,kcksratio
  ! Counters
  integer  :: i, j, k
  ! Pi value
  real(WP) :: pi, first_p, first_g, last_g, last_p, tmp1
  ! Check if min and max bounds are defined
  logical  :: is_max, is_min, isfluctuating
  rms_val=0.0_WP
  call parser_read(trim(name_s) // ' mean', mean_val)
  call parser_is_defined(trim(name_s) // ' Fluctuations',isfluctuating)
  if (isfluctuating) then
     call parser_read(trim(name_s) // ' Fluctuations', rms_val)
  end if
  if (geometry.eq.'1D') then
     ! Scalar initialization for 1-D as Sin profile   
     pi = acos(-1.0_WP) ! Create pi
     if (rms_val.gt.0.0_WP) then
        call parser_read(trim(name_s) //' front length',last_p)
        first_p=((L/2.0_WP)-last_p)/2.0_WP
        last_p=min(first_p+last_p,L/2.0_WP)
        first_g=floor(first_p/dx)
        last_g=floor(last_p/dx)
     else
        first_p=0.0_WP
        first_g=1
        last_g = -1
     end if
     tmp_sc=mean_val-rms_val/2.0_WP
     do i=first_g,last_g
        ! Linear profile of scalar S is the default for 1D fluctuating scalar
        tmp_sc(i,:,:) = tmp_sc(i,:,:)+rms_val*(i-first_g)/(last_g-first_g)
     end do
     do i=last_g+1, last_g+2*first_g
        tmp_sc(i,:,:) = tmp_sc(i,:,:)+rms_val
     end do
     do i=last_g+1+2*first_g, last_g+1+2*first_g+last_g-first_g
        ! Linear profile of scalar S is the default for 1D fluctuating scalar
        tmp_sc(i,:,:) = tmp_sc(i,:,:)-rms_val*(i-(last_g+1+2*first_g+last_g-first_g))/(last_g-first_g)
     end do
     tmp_mean = mean_val
     tmp_rms = rms_val
  else
     ! Check for 1 gp
     if (nx.eq.1) rms_val = 0.0_WP
     if (rms_val.eq.0.0_WP) then
        tmp_sc = mean_val
        tmp_mean = mean_val
        tmp_rms = 0.0_WP
     else
        ! Find min max
        max_val=10000.0_WP                 ! Default max bound
        min_val=0.0_WP                     ! Default min bound
        ! Check if a maximum value for scalar is defined. Read the value if it is
        call parser_is_defined(trim(name_s) // ' max',is_max)
        if (is_max) then
           call parser_read(trim(name_s) // ' max' , max_val)
        end if
        ! Check if a minimum value for scalar is defined. Read the value if it is
        call parser_is_defined(trim(name_s) // ' min',is_min)
        if (is_min) then
           call parser_read(trim(name_s) // ' min' , min_val)
        end if
        ! Allowed values ::: PP, D-Delta, VKP
        call parser_is_defined(trim(name_s) // ' Spectrum form',isfluctuating)
        if (isfluctuating) then
           call parser_read(trim(name_s) // ' Spectrum form',spectrum)
        else
           spectrum = 'D-Delta'
        end if
        if (trim(spectrum).eq.'PP') then
           call parser_read(trim(name_s)// ' Energetic scale',le)
           call ignition_spectrum(spectrum, rms_val, le, 0.0_WP, 0.0_WP,'Scalar',trim(name_s))
        elseif (trim(spectrum).eq.'VKP') then
           call parser_read(trim(name_s)// ' Dissipative scale',ld)
           call parser_read(trim(name_s)// ' Dissipation',epsilon)
           call ignition_spectrum(spectrum, rms_val, 0.0_WP, ld, epsilon,'Scalar',trim(name_s))
        elseif (trim(spectrum).eq.'D-Delta') then
           call parser_read(trim(name_s)// ' ks/ko',ksk0ratio)
           call parser_read(trim(name_s)// ' kc/ks',kcksratio)
           call ignition_doubledelta(ksk0ratio,kcksratio,trim(name_s))
        end if
        if (trim(spectrum).ne.'D-Delta') then
           tmp_sc=tmp_sc+mean_val
        end if
        ! Call bound, mean checking routine
        call ignition_scalarcheck
     end if
  end if

  return
end subroutine ignition_scalar

! ===================================================== !
! Initialize a double delta field                       !
! ===================================================== !
subroutine ignition_doubledelta(ksk0ratio,kcksratio,out_name)
  use ignition
  use precision
  use random
  implicit none
  real(WP),intent(in) :: ksk0ratio,kcksratio
  real(WP) :: pi,dk,kc,ks,kx,ky,kz,kk,f_phi
  integer :: i,j,k,nk,kk2, iter
  integer(KIND=8) :: plan_r2c,plan_c2r
  real(WP) :: rand
  complex(WP) :: ii=(0.0_WP,1.0_WP)
  character (len=*) :: out_name  
  logical :: is_min, is_max, right_mean, is_rms
  real(WP) :: temp1, temp2, factor, old_factor
  call random_init ! Initialize the random number generator
  pi = acos(-1.0_WP) ! Create pi
  ! Spectrum computation
  dk = 2.0_WP*pi/L
  ! Initialize in similar manner to Eswaran and Pope 1988
  ks = ksk0ratio*dk
  kc = kcksratio*ks
  nk = nx/2+1
  ! Inverse Fourier transform
  allocate(Cbuf(nk,ny,nz))
  allocate(Rbuf(nx,ny,nz))
  allocate(Rbuf_temp(nx,ny,nz))
  if (geometry.eq.'3D') then
     call dfftw_plan_dft_c2r_3d(plan_c2r,nx,ny,nz,Cbuf,Rbuf,FFTW_ESTIMATE)
     call dfftw_plan_dft_r2c_3d(plan_r2c,nx,ny,nz,Rbuf,Cbuf,FFTW_ESTIMATE)
  else if (geometry.eq.'2D') then
     call dfftw_plan_dft_c2r_2d(plan_c2r,nx,ny,Cbuf,Rbuf,FFTW_ESTIMATE)
     call dfftw_plan_dft_r2c_2d(plan_r2c,nx,ny,Rbuf,Cbuf,FFTW_ESTIMATE)  
  end if

  ! Compute the Fourier coefficients
  do k=1,nz
     do j=1,ny
        do i=1,nk
           ! Wavenumbers
           kx=real(i-1,WP)*dk
           ky=real(j-1,WP)*dk
           if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
           kz=real(k-1,WP)*dk
           if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
           kk =sqrt(kx**2+ky**2+kz**2)
           kk2=sqrt(kx**2+ky**2)

           if ((ks-dk/2.0_WP.le.kk).and.(kk.le.ks+dk/2.0_WP)) then
              f_phi = 1.0_WP
           else
              f_phi = 0.0_WP
           end if
           call random_number(rand)
           if (kk.lt.1e-10) then
              Cbuf(i,j,k) = 0.0_WP
           else
              Cbuf(i,j,k) = sqrt(f_phi/(4.0_WP*pi*kk**2))*exp(ii*2.0_WP*pi*rand)
           end if
        end do
     end do
  end do

  ! Oddball
  if (geometry.eq.'3D') then
     do k=2,nz
        do j=nk+1,ny
           Cbuf(1,j,k)=conjg(Cbuf(1,ny+2-j,nz+2-k))
        end do
     end do
     do k=nk+1,nz
        Cbuf(1,1,k)=conjg(Cbuf(1,1,nz+2-k))
     end do
  else if (geometry .eq.'2D') then
     do j=nk+1,ny
        Cbuf(1,j,1)=conjg(Cbuf(1,ny+2-j,1))
     end do
  end if

  ! Inverse Fourier transform
  call dfftw_execute(plan_c2r)

  call parser_is_defined(trim(out_name) // ' Fluctuations', is_rms)
  ! Check if a maximum value for scalar is defined. Read the value if it is
  call parser_is_defined(trim(out_name) // ' max',is_max)
  if (.not.is_max) then
     if (is_rms) then
        max_val = mean_val + rms_val
     else
        max_val = 0.0_WP
     end if
  end if
  ! Check if a minimum value for scalar is defined. Read the value if it is
  call parser_is_defined(trim(out_name) // ' min',is_min)
  if (.not.is_min) then
     if (is_rms) then
        min_val = mean_val - rms_val
     else
        min_val = 0.0_WP
     end if
  end if

  ! Needed for normalization of Rbuf
  temp1 = sum(Rbuf)/real(nx*ny*nz,WP)
  temp2 = 0.0_WP
  do k=1,nz
     do j=1,ny
        do i=1,nx
           temp2 = temp2 + (Rbuf(i,j,k) - temp1)**2
        end do
     end do
  end do
  temp2 = sqrt(temp2/real(nx*ny*nz,WP))
  ! Normalize
  Rbuf = (Rbuf-temp1)/temp2
  factor = 1.0_WP/((max_val-mean_val)/(max_val-min_val))
  old_factor = factor
  if (max_val.le.min_val) max_val = min_val+2.0_WP*rms_val
  iter = 0
  right_mean = .false.
  ! Iterate to get right mean
  do while(iter.lt.1000 .and. .not.right_mean)
     iter = iter+1
     Rbuf_temp = 0.0_WP
     temp1 = maxval(Rbuf) - (maxval(Rbuf) - minval(Rbuf))/factor
     ! Force 'double-delta' pdf on scalar field
     do k=1,nz
        do j=1,ny
           do i=1,nx
              if (Rbuf(i,j,k).le.temp1) then
                 Rbuf_temp(i,j,k) = min_val
              else
                 Rbuf_temp(i,j,k) = max_val
              end if
           end do
        end do
     end do
     temp2 = sum(Rbuf_temp)/real(nx*ny*nz,WP)
     if (abs(temp2-mean_val)/abs(mean_val).lt.0.001) then
        right_mean = .true.
     elseif (temp2.gt.mean_val) then
        factor = factor*1.01
        factor = 0.5_WP * (factor+old_factor)
        old_factor = factor
     else
        factor = factor/1.01*(1.0+real(iter/10000,WP))
        factor = 0.5_WP * (factor+old_factor)
        old_factor = factor
     end if     
  end do
  if (.not.right_mean) write (*,*) 'Increase no of iterations. Mean wrong'
  Rbuf = Rbuf_temp
  Rbuf_temp = 0.0_WP

  ! Fourier Transform and filter to smooth
  call dfftw_execute(plan_r2c)

  do k=1,nz
     do j=1,ny
        do i=1,nk
           ! Wavenumbers
           kx=real(i-1,WP)*dk
           ky=real(j-1,WP)*dk
           if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
           kz=real(k-1,WP)*dk
           if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
           kk =sqrt(kx**2+ky**2+kz**2)
           kk2=sqrt(kx**2+ky**2)

           ! Filter to remove high wavenumber components
           if (kk.le.kc) then
              Cbuf(i,j,k) = Cbuf(i,j,k) * 1.0_WP
           else
              Cbuf(i,j,k) = Cbuf(i,j,k) * (kc/kk)**2
           end if
        end do
     end do
  end do

  ! Oddball
  if (geometry.eq.'3D') then
     do k=2,nz
        do j=nk+1,ny
           Cbuf(1,j,k)=conjg(Cbuf(1,ny+2-j,nz+2-k))
        end do
     end do
     do k=nk+1,nz
        Cbuf(1,1,k)=conjg(Cbuf(1,1,nz+2-k))
     end do
  else if (geometry .eq.'2D') then
     do j=nk+1,ny
        Cbuf(1,j,1)=conjg(Cbuf(1,ny+2-j,1))
     end do
  end if

  ! Fourier Transform back to real
  call dfftw_execute(plan_c2r)

  do k=1,nz
     do j=1,ny
        do i=1,nx
           tmp_sc(i,j,k) = Rbuf(i,j,k)/real(nx*ny*nz,WP)
        end do
     end do
  end do

  ! Destroy the plans
  call dfftw_destroy_plan(plan_c2r)
  call dfftw_destroy_plan(plan_r2c)

  ! Clean up
  deallocate(Cbuf)
  deallocate(Rbuf)
  return
end subroutine ignition_doubledelta

! ===================================================== !
! Initialize PP or VKP spectrum for velocity or scalars !
! ===================================================== !
subroutine ignition_spectrum(spectrum, Ut, le, ld, epsilon, type, out_name)
  use ignition
  use random
  use fileio
  implicit none
    
  ! Turbulent velocity
  real(WP) :: Ut
  
  ! Spectrum type
  character(len=*) :: spectrum
  real(WP) :: le,ld,epsilon
  
  ! Spectrum computation
  real(WP) :: psr,ps1,ps2,ke,kd,dk,kc,kk,kx,ky,kz,kk2
  real(WP) :: alpha,spec_amp,eps,amp_disc,e_total,energy_spec
  real(WP), dimension(:,:), pointer :: spect
  complex(WP), dimension(:,:,:), pointer :: ak,bk
  integer  :: nk ! Cutoff wave number 
        
  ! Other
  integer :: i,j,k,ik,iunit
  complex(WP) :: ii=(0.0_WP,1.0_WP)
  real(WP) :: rand,pi
  
  ! Fourier coefficients
  integer(KIND=8) :: plan_r2c,plan_c2r
  complex(WP), dimension(:,:,:), pointer :: Uk,Vk,Wk

  ! Type of quantity
  character (len=*) :: type
  character (len=*) :: out_name  

  call random_init ! Initialize the random number generator
  pi = acos(-1.0_WP) ! Create pi
  ! Spectrum computation
  ke = 2.0_WP*pi/le
  if (trim(spectrum).eq.'VKP') kd = 2.0_WP*pi/ld
  dk = 2.0_WP*pi/L
  kc = real(nx/2,WP)*dk
  
  eps=ke/1000000.0_WP
  if (trim(spectrum).eq.'PP') then
     if (geometry.eq.'3D') then
        spec_amp = 16.0_WP*sqrt(2.0_WP/pi)*Ut**2/ke
     else if (geometry.eq.'2D') then
        spec_amp = (32.0_WP/3.0_WP)*sqrt(2.0_WP/pi)*Ut**2/ke
     end if
  else if (trim(spectrum).eq.'VKP') then
     alpha=1.5_WP
     spec_amp = 1.5_WP*Ut**5/epsilon
  end if
  amp_disc = sqrt(dk)**3

  ! Compute spectrum
  nk = nx/2+1
  allocate(ak(nk,ny,nz),bk(nk,ny,nz))
  do k=1,nz
     do j=1,ny
        do i=1,nk
           ! Random numbers
           call random_number(rand)
           psr=2.0_WP*pi*(rand-0.5_WP)
           call random_number(rand)
           ps1=2.0_WP*pi*(rand-0.5_WP)
           call random_number(rand)
           ps2=2.0_WP*pi*(rand-0.5_WP)
           ! Wavenumbers
           kx=real(i-1,WP)*dk
           ky=real(j-1,WP)*dk
           if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
           kz=real(k-1,WP)*dk
           if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
           kk=sqrt(kx**2+ky**2+kz**2)
           ! Spectrums
           if (trim(spectrum).eq.'PP') then
              energy_spec=spec_amp*(kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)
           else if (trim(spectrum).eq.'VKP') then
              energy_spec=spec_amp*(kk/ke)**4/(1.0_WP+(kk/ke)**2)&
                   **(17.0_WP/6.0_WP)*exp(-1.5_WP*alpha*(kk/kd)**(4.0_WP/3.0_WP))
           end if
           ! Coeff
           ak(i,j,k)=0.0_WP
           bk(i,j,k)=0.0_WP
           if ((kk.gt.eps).and.(kk.le.kc)) then
              if (type.eq.'V') then
                 if (geometry.eq.'3D') then
                    ak(i,j,k)=amp_disc*sqrt(energy_spec/(2.0_WP*pi*kk**2))*&
                         exp(ii*ps1)*cos(psr)
                    bk(i,j,k)=amp_disc*sqrt(energy_spec/(2.0_WP*pi*kk**2))*&
                         exp(ii*ps2)*sin(psr)
                 else if (geometry.eq.'2D') then
                    ak(i,j,k)=dk*sqrt(energy_spec/(1.0_WP*pi*kk**1))*exp(ii*ps1)
                 end if
              else if (type .eq. 'Scalar') then
                 if (geometry.eq.'3D') then
                    ak(i,j,k)=sqrt(1.0_WP/3.0_WP)*amp_disc*sqrt(energy_spec&
                         /(2.0_WP*pi*kk**2))*(exp(ii*ps1)*cos(psr)+exp(ii*ps2)*sin(psr))
                 else if (geometry.eq.'2D') then
                    ak(i,j,k)=(0.5_WP)*dk*sqrt(energy_spec/(0.5_WP*pi*kk**1))*exp(ii*ps1)
                 end if
              end if
           end if
        end do
     end do
  end do
  ! Output spectrum (1D) for comparison
  allocate(spect(2,nx+1))  ! (wavenumber, magnitude)
  e_total=0.0_WP
  do ik=1,nx+1
     kk = dk*(ik-1)
     spect(1,ik)=kk
     if (trim(spectrum).eq.'PP') then
        energy_spec=spec_amp*(kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)
     else if (trim(spectrum).eq.'VKP') then
        energy_spec=spec_amp*(kk/ke)**4/(1.0_WP+(kk/ke)**2)**(17.0_WP/6.0_WP)*&
             exp(-1.5_WP*alpha*(kk/kd)**(4.0_WP/3.0_WP))
     end if
     if ((kk.gt.eps).and.(kk.le.kc)) then
        spect(2,ik) = energy_spec
        e_total=e_total+dk*energy_spec
     else
        spect(2,ik) = 0.0_WP
     end if
  end do

  ! Write spectrum to file
  iunit=iopen()
  open(iunit,file = trim(out_name) // '_spectrum.analytic',form='formatted') 
  do i=1,nx+1
     if ((spect(1,i).ne.0.0_WP).and.(spect(2,i).ne.0.0_WP)) then
        write(iunit,*) spect(1,i),'    ',spect(2,i)  
     end if
  end do
  close(iclose(iunit))

  ! Compute 3D field
  allocate(Uk(nk,ny,nz))
  allocate(Vk(nk,ny,nz))
  allocate(Wk(nk,ny,nz))
  Uk=(0.0_WP,0.0_WP)
  Vk=(0.0_WP,0.0_WP)
  Wk=(0.0_WP,0.0_WP)

  if (type.eq.'Scalar') then
     Uk=ak  
  else
     ! Compute the Fourier coefficients
     do k=1,nz
        do j=1,ny
           do i=1,nk
              ! Wavenumbers
              kx=real(i-1,WP)*dk
              ky=real(j-1,WP)*dk
              if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
              kz=real(k-1,WP)*dk
              if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
              kk =sqrt(kx**2+ky**2+kz**2)
              kk2=sqrt(kx**2+ky**2)
              
              if ((kk.gt.eps).and.(kk.le.kc)) then
                 if (geometry .eq. '3D') then
                    Wk(i,j,k)=-bk(i,j,k)*kk2/kk
                 end if
                 if (kk2.lt.eps) then
                    Uk(i,j,k)=(ak(i,j,k)+bk(i,j,k))/sqrt(2.0_WP)
                    Vk(i,j,k)=(bk(i,j,k)-ak(i,j,k))/sqrt(2.0_WP)
                 else
                    Uk(i,j,k)=(ak(i,j,k)*kk*ky+bk(i,j,k)*kx*kz)/(kk*kk2)
                    Vk(i,j,k)=(bk(i,j,k)*ky*kz-ak(i,j,k)*kk*kx)/(kk*kk2)
                 end if
              end if
           end do
        end do
     end do
  end if
  
  ! Oddball
  if (geometry.eq.'3D') then
     do k=2,nz
        do j=nk+1,ny
           Uk(1,j,k)=conjg(Uk(1,ny+2-j,nz+2-k))
           Vk(1,j,k)=conjg(Vk(1,ny+2-j,nz+2-k))
           Wk(1,j,k)=conjg(Wk(1,ny+2-j,nz+2-k))
        end do
     end do
     do k=nk+1,nz
        Uk(1,1,k)=conjg(Uk(1,1,nz+2-k))
        Vk(1,1,k)=conjg(Vk(1,1,nz+2-k))
        Wk(1,1,k)=conjg(Wk(1,1,nz+2-k))
     end do
  else if (geometry .eq.'2D') then
     do j=nk+1,ny
        Uk(1,j,1)=conjg(Uk(1,ny+2-j,1))
        Vk(1,j,1)=conjg(Vk(1,ny+2-j,1))
     end do
  end if
  
  ! Inverse Fourier transform
  allocate(Cbuf(nk,ny,nz))
  allocate(Rbuf(nx,ny,nz))
  if (geometry.eq.'3D') then
     call dfftw_plan_dft_c2r_3d(plan_c2r,nx,ny,nz,Cbuf,Rbuf,FFTW_ESTIMATE)
     call dfftw_plan_dft_r2c_3d(plan_r2c,nx,ny,nz,Rbuf,Cbuf,FFTW_ESTIMATE)
  else if (geometry.eq.'2D') then
     call dfftw_plan_dft_c2r_2d(plan_c2r,nx,ny,Cbuf,Rbuf,FFTW_ESTIMATE)
     call dfftw_plan_dft_r2c_2d(plan_r2c,nx,ny,Rbuf,Cbuf,FFTW_ESTIMATE)  
  end if

  ! Execute the plans
  if (type .eq. 'V') then
      Cbuf = Uk
     call dfftw_execute(plan_c2r)
     U = Rbuf
     Cbuf = Vk 
     call dfftw_execute(plan_c2r)
     V = Rbuf
     Cbuf = Wk
     call dfftw_execute(plan_c2r)
     W = Rbuf
     ! Check the results for velocity
     call ignition_check(spectrum,Ut,ke)

  else if (type .eq. 'Scalar') then
     Cbuf = Uk
     call dfftw_execute(plan_c2r)
     tmp_sc = Rbuf
  end if
        
  ! Clean up
  deallocate(Uk)
  deallocate(Vk)
  deallocate(Wk)
  deallocate(ak)
  deallocate(bk)
end subroutine ignition_spectrum

!---------------------------------------!
! Check the stats
!---------------------------------------!
subroutine ignition_check(spectrum,Ut,ke)
  use ignition
  use fileio
  use parser
  implicit none
  integer :: iunit,i
  real(WP) :: pi
  real(WP) :: Ut,ke,mu_cin
  character(len=*) :: spectrum
  real(WP) :: kcin,epsi,l11,lt,re_turb,tau_epsi,l_k,tau_k
  real(WP), dimension(:,:), pointer :: spect
  
  pi = acos(-1.0_WP)
  call parser_read('Viscosity',mu_cin)
  if (geometry.eq.'2D') then
     kcin = Ut*Ut
  else
     kcin = 3.0_WP*0.5_WP*Ut*Ut
  end if
  if (trim(spectrum).eq.'PP') then
     epsi  = 15.0_WP*Ut*Ut*ke*ke*mu_cin/4.0_WP
     l11   = sqrt(2.0_WP*pi)/ke
  else if (trim(spectrum).eq.'VKP') then
     call parser_read('Dissipation',epsi)
  end if
  lt       = Ut**3.0_WP/epsi
  re_turb  = Ut**4.0_WP/mu_cin/epsi
  tau_epsi = kcin/epsi
  l_k      = (mu_cin**0.75_WP)/(epsi**0.25_WP)
  tau_k    = sqrt(mu_cin/epsi)
  
  write(*,*)  
  write(*,*)' ================================================ '
  write(*,*)' Debugging turbulent values based on input data              '
  write(*,*)' ------------------------------------------------ '
  if (trim(spectrum).eq.'PP') then
     write(*,*)' -Spectrum type ----------> Passot-Pouquet'
  else if (trim(spectrum).eq.'VKP') then
     write(*,*)' -Spectrum type ----------> Von Karman-Pao'
  end if
  write(*,*)' -Turbulent Reynolds '
  write(*,*)'   Re_t --------> ', re_turb
  write(*,*)' -Turbulent kinetic energy '
  write(*,*)'   k -----------> ', kcin
  write(*,*)' -Turbulent dissipation rate '
  write(*,*)'   epsilon -----> ', epsi
  write(*,*)' -Size of the biggest eddy '
  write(*,*)'   l_t ---------> ', lt
  write(*,*)'   C1 ----------> ', L/lt
  write(*,*)' -Eddy turn over time '
  write(*,*)'   tau ---------> ', tau_epsi
  write(*,*)' -Kolmogorov scale '
  write(*,*)'   l_k ---------> ', l_k
  write(*,*)'   C2 ----------> ', l_k/dx
  write(*,*)' -Kolmogorov time scale '
  write(*,*)'   tau_k -------> ', tau_k
  write(*,*)' -Reynolds lambda '
  write(*,*)'   Re_lambda ---> ', sqrt(re_turb*15.0_WP)
  if (trim(spectrum).eq.'PP') then
     write(*,*)' -Integral length scale '
     write(*,*)'   Lii --------> ', l11
     write(*,*)' -Turbulent Reynolds based on Lii '
     write(*,*)'   Re_Lii -----> ', Ut*l11/mu_cin
  end if
  write(*,*)' ======================================== '
  write(*,*)
  write(*,*)
  
  ! Recompute the spectrum
  allocate(spect(2,nx+1))
  spect = 0.0_WP
  call get_spectrum_actual(nx,ny,nz,L,U,spect)
  call get_spectrum_actual(nx,ny,nz,L,V,spect)
  call get_spectrum_actual(nx,ny,nz,L,W,spect)
  
  ! Output it
  iunit=iopen()
  open(iunit,file='final_velocity_spectrum.numeric',form='formatted')
  do i=1,nx+1
     if ((spect(1,i).ne.0.0_WP).and.(spect(2,i).ne.0.0_WP)) then
        write(11,*) spect(1,i),',',spect(2,i)
     end if
  end do
  close(iclose(iunit))
  return
end subroutine ignition_check

! ===================================================== !
! This routine checks the mean and bounds of scalar     !
! ===================================================== !
subroutine ignition_scalarcheck
  use ignition
  implicit none
  ! Temperory values
  real(WP) :: tmp1
  ! Counters
  integer  :: i, j, k, iterate
  ! Check vars
  logical  :: right_mean, right_min, right_max
  right_mean=.false.
  right_min = .false.
  right_max = .false.
  iterate = 0
  ! Calculate actual mean and rms. Shift mean accordingly.
  ! Continue this process until the right mean is obtained with the bounded values
  if (minval(tmp_sc).ge.min_val) right_min = .true.
  if (maxval(tmp_sc).le.max_val) right_max = .true.
  tmp_mean = sum(tmp_sc)/real(nx*ny*nz,WP)
  if (abs(tmp_mean-mean_val)/abs(mean_val).lt.0.001) right_mean = .true.
  do while (.not.right_mean .or. .not.right_min .or. .not.right_max.and.iterate.lt.1000)
     iterate = iterate+1
     tmp1 = 0.0_WP    
     if (.not.right_min) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 if (tmp_sc(i,j,k).lt.min_val) then
                    tmp1 = tmp1 - tmp_sc(i,j,k) + min_val
                    tmp_sc(i,j,k) = min_val
                 end if
              end do
           end do
        end do
        right_min = .true.
     end if
     if (.not.right_max) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 if (tmp_sc(i,j,k).gt.max_val) then
                    tmp1 = tmp1 - tmp_sc(i,j,k) + max_val
                    tmp_sc(i,j,k) = max_val
                 end if
              end do
           end do
        end do
        right_max = .true.
     end if
     if (tmp1/real(nx*ny*nz,WP) .gt. 1E-8) then
        right_min = .false.; right_max = .false.
        tmp_sc = tmp_sc - tmp1/real(nx*ny*nz,WP)
     end if
     tmp_mean = sum(tmp_sc)/real(nx*ny*nz,WP)
     if (abs(tmp_mean-mean_val)/abs(mean_val).lt.0.001) right_mean = .true.

     if (iterate.gt.100) then
        write (*,*) iterate, ' iteration, error remaining is ', tmp1
        write (*,*) tmp_mean, mean_val
     end if
  end do

  ! Calculate the rms valu
  tmp1 = 0.0_WP
  do k=1,nz
     do j=1,ny
        do i=1,nx
           tmp1  = tmp1 + (tmp_sc(i,j,k)-tmp_mean)**2
        end do
     end do
  end do
  tmp_rms = sqrt(tmp1/real(nx*ny*nz,WP))

  return
end subroutine ignition_scalarcheck

!==============================================================!
! Check turbulence statistics                                  !
!==============================================================!
subroutine ignition_stat
  use ignition
  use parser
  use fileio
  implicit none
  integer :: i,j,k
  real(WP) :: upvp,upwp,vpwp,Ut_final,kcin_final,div,div_max,div_b
  real(WP), dimension(3) :: um,umax,umin,uminval
  
  ! Computing mean velocity
  um=0.0_WP
  umax=-1.0e+70_WP
  umin=99.0e+70_WP
  uminval=99.0e+70_WP
  do k=1,nz
     do j=1,ny
        do i=1,nx
           um(1)=um(1)+U(i,j,k)
           umax(1)=max(umax(1),(U(i,j,k)))
           umin(1)=min(umin(1),(U(i,j,k)))
           uminval(1)=min(uminval(1),abs(U(i,j,k)))
           um(2)=um(2)+V(i,j,k)
           umax(2)=max(umax(2),(V(i,j,k)))
           umin(2)=min(umin(2),(V(i,j,k)))
           uminval(2)=min(uminval(2),abs(V(i,j,k)))
           um(3)=um(3)+W(i,j,k)
           umax(3)=max(umax(3),(W(i,j,k)))
           umin(3)=min(umin(3),(W(i,j,k)))
           uminval(3)=min(uminval(3),abs(W(i,j,k)))
        end do
     end do
  end do
  um(1)=um(1)/real(nx*ny*nz,WP)
  um(2)=um(2)/real(nx*ny*nz,WP)
  um(3)=um(3)/real(nx*ny*nz,WP)
  
  ! Computing cross correlations to confirm isotropy
  upvp=0.0_WP
  upwp=0.0_WP
  vpwp=0.0_WP
  do k=1,nz 
     do j=1,ny
        do i=1,nx
           upvp=upvp+(U(i,j,k)-um(1))*(V(i,j,k)-um(2))
           upwp=upwp+(U(i,j,k)-um(1))*(W(i,j,k)-um(3))
           vpwp=vpwp+(V(i,j,k)-um(2))*(W(i,j,k)-um(3))
        end do
     end do
  end do
  upvp = upvp / (nx*ny*nz)
  upwp = upwp / (nx*ny*nz)
  vpwp = vpwp / (nx*ny*nz)
  
  ! Computing the final Ut
  Ut_final = 0.0_WP
  do k=1,nz
     do j=1,ny
        do i=1,nx
           Ut_final = Ut_final + (U(i,j,k)-um(1))**2+(V(i,j,k)-um(2))**2+(W(i,j,k)-um(3))**2
        end do
     end do
  end do
  if (geometry.eq.'3D') then
     Ut_final = sqrt(Ut_final/nx/ny/nz/3.0_WP)
     kcin_final = 3.0_WP*0.5_WP*Ut_final*Ut_final
  else if (geometry.eq.'2D') then
     Ut_final = sqrt(Ut_final/nx/ny/2.0_WP)
     kcin_final = Ut_final*Ut_final
  end if

  ! Computing divergence values
  div_b=0.0_WP
  div_max=0.0_WP
  if (nz .gt. 1) then
     do k=2,nz-1
        do j=2,ny-1
           do i=2,nx-1
              div=abs(RHO(i+1,j,k)*U(i+1,j,k)-RHO(i,j,k)*U(i,j,k)+&
                   RHO(i,j+1,k)*V(i,j+1,k)-RHO(i,j,k)*V(i,j,k)+RHO(i,j,k+1)*W(i,j,k+1)-RHO(i,j,k)*W(i,j,k))/dx
              div_max=max(div_max,div)
              div_b=div_b+div
           end do
        end do
     end do
     div_b=div_b/(nx-2)/(ny-2)/(nz-2)
  else
     do j=2,ny-1
        do i=2,nx-1
           div=abs(RHO(i+1,j,1)*U(i+1,j,1)-RHO(i,j,1)*U(i,j,1)+&
                   RHO(i,j+1,1)*V(i,j+1,1)-RHO(i,j,1)*V(i,j,1))/dx
           div_max=max(div_max,div)
           div_b=div_b+div
        end do
     end do
     div_b=div_b/(nx-2)/(ny-2)
  end if

  ! Output it ------------------------------------------
  write(*,*)  
  write(*,*)' ======================================== '
  write(*,*)' Debugging turbulent values given by      '
  write(*,*)' the generated field                      '
  write(*,*)' ---------------------------------------- '
  write(*,*)' -Turbulent velocity '
  write(*,*)'   up ---------> ', Ut_final
  write(*,*)' -Turbulent kinetic energy '
  write(*,*)'   k ----------> ', kcin_final
  write(*,*)' -Cross correlations'
  write(*,*)'   <upvp> -----> ', upvp
  write(*,*)'   <upwp> -----> ', upwp
  write(*,*)'   <vpwp> -----> ', vpwp
  write(*,*)' -Velocity min, max and mean'
  write(*,*)'   u_min ------> ', umin(1)
  write(*,*)'   u_mean -----> ', um(1)
  write(*,*)'   u_max ------> ', umax(1)
  write(*,*)'   |u|_min ----> ', uminval(1)
  write(*,*)'   v_min ------> ', umin(2)
  write(*,*)'   v_mean -----> ', um(2)
  write(*,*)'   v_max ------> ', umax(2)
  write(*,*)'   |v|_min ----> ', uminval(2)
  write(*,*)'   w_min ------> ', umin(3)
  write(*,*)'   w_mean -----> ', um(3)
  write(*,*)'   w_max ------> ', umax(3)
  write(*,*)'   |w|_min ----> ', uminval(3)
  write(*,*)' -Average divergence -1 order'
  write(*,*)'   Div_b ------> ', div_b
  write(*,*)' -Maximum divergence'
  write(*,*)'   Div_max ----> ', div_max
  write(*,*)' ======================================== '
  write(*,*)
  write(*,*)
  
  return
end subroutine ignition_stat

!==============================================================!
!  Get the actual spectrum for velocity                        !
!==============================================================!
subroutine get_spectrum_actual(nx,ny,nz,L,A,S)
  use precision
  implicit none
  include 'fftw3.f'  
  integer, intent(in) :: nx,ny,nz
  real(WP), intent(in) :: L
  real(WP), dimension(nx,ny,nz), intent(in) :: A
  real(WP), dimension(2,nx+1), intent(out) :: S
  complex(WP), dimension(:,:,:), pointer :: in,out
  real(WP), dimension(:,:,:), pointer :: B
  real(WP) :: pi,dk,kc,eps,kx,ky,kz,kk
  integer(KIND=8) :: plan
  integer :: i,j,k,ik
  
  ! Create the plan
  allocate(in(nx,ny,nz))
  allocate(out(nx,ny,nz))
  call dfftw_plan_dft_3d(plan,nx,ny,nz,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
  
  ! FFT
  allocate(B(nx,ny,nz))
  do k=1,nz
     do j=1,ny
        do i=1,nx
           in(i,j,k) = A(i,j,k)
        end do
     end do
  end do
  call dfftw_execute(plan)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           if (nz==1) then
              out(i,j,k) = out(i,j,k)/(nx**2)
           else
              out(i,j,k) = out(i,j,k)/(nx**3)
           end if
           B(i,j,k) = sqrt(real(out(i,j,k)*conjg(out(i,j,k))))
        end do
     end do
  end do
  
  ! Compute the spectrum
  pi=acos(-1.0_WP)
  dk=2.0_WP*pi/L
  kc=pi*nx/L
  eps=kc/1000000.0_WP
  do k=1,nz
     do j=1,ny
        do i=1,nx
           ! Wavenumbers
           kx=real(i-1,WP)*dk
           if (i.gt.(nx/2+1)) kx=-real(nx-i+1,WP)*dk
           ky=real(j-1,WP)*dk
           if (j.gt.(ny/2+1)) ky=-real(ny-j+1,WP)*dk
           kz=real(k-1,WP)*dk
           if (k.gt.(nz/2+1)) kz=-real(nz-k+1,WP)*dk
           kk=sqrt(kx**2+ky**2+kz**2)
           ! Spectrum
           ik=1+idint(kk/dk+0.5_WP)
           if ((kk.gt.eps).and.(kk.le.kc)) then
              S(2,ik)=S(2,ik)+0.5_WP*(B(i,j,k)**2)/dk
           end if
        end do
     end do
  end do
  do ik=1,nx+1
     S(1,ik)=dk*(ik-1)
  end do
  
  ! Clean up
  call dfftw_destroy_plan(plan)
  deallocate(B)
  deallocate(out)
  deallocate(in)
  
  return
end subroutine get_spectrum_actual

! ============================================ !
! Get varibale pdf distribution                !
! ============================================ !
subroutine get_1dpdf(S,nbin,what)
  use ignition
  use fileio
  implicit none
  
  integer :: iunit,i,j,k,p1
  character(len=str_medium) :: file
  character(len=str_short) :: what
  integer :: nbin, n
  real(WP), dimension(1:2,1:nbin) :: pdf
  real(WP) :: minS,maxS
  real(WP), dimension(1:nx,1:ny,1:nz) :: S
  real(WP) :: dpdf, diff
  minS=minval(S)
  maxS=maxval(S)
  pdf = 0.0_WP
  dpdf = (maxS-minS)/(nbin-1)
  diff = 0.5_WP*dpdf
  if (minS.ne.maxS) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              do n = 1, nbin
                 ! For non edge values
                 if (abs(S(i,j,k)-minS-((n-1)*dpdf)).lt.diff) then
                    pdf(2,n) = pdf(2,n) + 1.0_WP/real(nx*ny*nz,WP)
                 end if
                 ! For edge value
                 if (S(i,j,k)-minS-((n-1)*dpdf).eq.diff) then
                    pdf(2,n+1) = pdf(2,n+1) + 1.0_WP/real(nx*ny*nz,WP)
                 end if
              end do
           end do
        end do
     end do

     ! Dump results to a file
     file='./pdf'//'-'//trim(what)//'.'
     iunit=iopen()
     open(iunit,file=trim(file),form='formatted')
     do p1=1,nbin
        pdf(1,p1) = (real((p1-1),WP)*(maxS-minS)/real((nbin-1),WP))+minS
        write(iunit,'(ES14.7,xx,ES14.7)') pdf(1,p1),pdf(2,p1)
     end do
     close(iclose(iunit))
  else
     write (*,*) trim(what)// ' has constant values'
  end if
  
  return
end subroutine get_1dpdf
