module combustion
  use precision
  use partition
  use data
  implicit none
  
  ! Default values in the case constant properties
  real(WP), dimension(:), pointer :: rho_input
  real(WP), dimension(:), pointer :: visc_input
  real(WP), dimension(:), pointer :: diff_input
  real(WP), dimension(:), pointer :: T_input
  real(WP) :: P_input
  real(WP) :: MM_input, GAMMA_input
  
  ! Number of filtering steps of dRHO
  integer :: nfilter
  
  ! Sum of dRHO for boundary
  real(WP) :: sum_dRHO,drho_old,drho_new
  
  ! Thermodynamic pressure (background pressure) => P_thermo(t) only
  real(WP) :: P_thermo
  real(WP) :: P_thermo_old
  
  ! For constant volume
  logical  :: is_constant_volume
  real(WP) :: RHO_0
  
  ! Simple links to transported scalars
  ! Mixture Fraction
  integer :: isc_ZMIX
  ! Mixture Fraction Squared (for the variance)
  integer :: isc_ZMIX2
  ! Mixture Fraction Variance
  integer :: isc_ZVAR
  ! Progress Variable
  integer :: isc_PROG
  ! Progress Variable Squared (for the variance)
  integer :: isc_PROG2
  ! Progress Variable Variance
  integer :: isc_PROGVAR
  ! Enthalpy
  integer :: isc_ENTH
  ! Internal energy
  integer :: isc_E
  ! Second Mixture Fraction 
  integer :: isc_ZSTAR
  ! Fuel Premixing Fraction, FMIX = ZSTAR / ZMIX
  integer :: isc_FMIX

  ! Whether or not to use each variable
  logical :: use_ZVAR       ! Nonprermixed FPV
  logical :: use_PROGVAR    ! Premixed FPV
  logical :: use_FMIX       ! Two Mixture Fraction FPV

  ! Temperature
  real(WP), dimension(:,:,:), pointer :: T
  ! Gamma
  real(WP), dimension(:,:,:), pointer :: Ggas
  ! Specific gas constant
  real(WP), dimension(:,:,:), pointer :: Rgas

  ! Other Flamelet modeling variables
  real(WP), dimension(:,:,:), pointer :: CHI   ! Scalar Dissipation Rate
  real(WP), dimension(:,:,:), pointer :: ZVAR, PROGVAR ! Scalar variances
  real(WP), dimension(:,:,:), pointer :: FMIX  ! Fuel premixing fraction

  ! Delta RHO
  real(WP), dimension(:,:,:), pointer :: delta_RHO

  ! Values to monitor
  real(WP) :: min_T, max_T, min_rho, max_rho, min_P, max_P
  integer  :: clip_drho

  ! Constants
  real(WP), parameter :: R_uni = 8.31434_WP

  ! For Multiple Mixture Fraction FPV - subfilter modeling
  logical :: use_FMIXsub
  real(WP) :: tolZ ! min value of Z for which F source term is calculated
  
  ! Whether or not to clip ZMIX
  logical :: clip_zmix

contains

  ! Compute the variance of a general variable
  subroutine combustion_update_variance(SCnow,VAR,isc_var,isc_var2,isc_varvar)
    implicit none
    
    integer :: i,j,k
    ! indices of the scalar, scalar squared, and scalar variance 
    integer, intent(in) :: isc_var, isc_var2, isc_varvar 
    real(WP), intent(in), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar) :: SCnow
    real(WP), intent(out), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: VAR
    real(WP) :: tvar
    
    if (isc_varvar.ne.0) then
       ! If the variance is present take it
       !$OMP PARALLEL DO
       do k=kmino_,kmaxo_
          do j=jmino_,jmaxo_
             do i=imino_,imaxo_
                VAR(i,j,k) = SCnow(i,j,k,isc_varvar)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (isc_var2.ne.0) then
       ! If scalar^2 is present use it to get scalar variance
       !$OMP PARALLEL DO
       do k=kmino_,kmaxo_
          do j=jmino_,jmaxo_
             do i=imino_,imaxo_
                VAR(i,j,k) = SCnow(i,j,k,isc_var2) - SCnow(i,j,k,isc_var)**2
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       ! Otherwise compute with dynamic model - newest SC used ! CHANGE HERE
       call sgsmodel_ZVAR(VAR,isc_var)
    end if
    
    ! Clip the computed Variance
    !$OMP PARALLEL DO PRIVATE(tvar)
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             tvar = max(0.0_WP,min(1.0_WP,SCnow(i,j,k,isc_var)))
             VAR(i,j,k) = max(0.0_WP,min(tvar*(1.0_WP-tvar), VAR(i,j,k)))
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    return
  end subroutine combustion_update_variance
  

  ! Compute the fuel premixing fraction
  subroutine combustion_FMIX(SCnow)
    implicit none
    
    integer :: i,j,k
    real(WP), intent(in), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar) :: SCnow
    
    if (isc_FMIX.ne.0) then
       ! If FMIX is present take it
       !$OMP PARALLEL DO
       do k=kmino_,kmaxo_
          do j=jmino_,jmaxo_
             do i=imino_,imaxo_
                FMIX(i,j,k) = SCnow(i,j,k,isc_FMIX)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (isc_ZSTAR.ne.0) then
       ! If ZSTAR is present use it to calculate FMIX
       !$OMP PARALLEL DO
       do k=kmino_,kmaxo_
          do j=jmino_,jmaxo_
             do i=imino_,imaxo_
                FMIX(i,j,k) = max(SCnow(i,j,k,isc_ZSTAR),1.0E-8_WP) &
                     / max(SCnow(i,j,k,isc_ZMIX), 1.0E-8_WP)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if
    
    ! Clip the computed FMIX
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             FMIX(i,j,k) = min(max(0.0_WP,FMIX(i,j,k)), 1.0_WP)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    return
  end subroutine combustion_FMIX

  ! Compute the scalar dissipation rate (algebraic model for LES)
  subroutine combustion_CHI
    use memory
    implicit none

    integer :: i,j,k
    
    call gradient_squared(SC(:,:,:,isc_ZMIX),CHI)
    !$OMP PARALLEL DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             CHI(i,j,k) = 2.0_WP*CHI(i,j,k)*DIFF(i,j,k,isc_ZMIX)/RHO(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    ! With soot, additional source terms in the variance transport equation show up in the dissipation rate model
    if (use_sgs .and. use_soot .and. trim(chemistry).eq.'chemtable') then
       call chemtable_lookup('Z2RHODOT',tmp1)
       !$OMP PARALLEL DO
       do k=kmino_,kmaxo_
          do j=jmino_,jmaxo_
             do i=imino_,imaxo_
                CHI(i,j,k) = CHI(i,j,k) - min(tmp1(i,j,k), 0.0_WP)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       call chemtable_lookup('RHODOT',tmp1)
       !$OMP PARALLEL DO
       do k=kmino_,kmaxo_
          do j=jmino_,jmaxo_
             do i=imino_,imaxo_
                CHI(i,j,k) = CHI(i,j,k) + min(tmp1(i,j,k), 0.0_WP)*SC(i,j,k,isc_ZMIX)**2/RHO(i,j,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       call chemtable_lookup('ZSRC_ZMIX',tmp1)
       !$OMP PARALLEL DO
       do k=kmino_,kmaxo_
          do j=jmino_,jmaxo_
             do i=imino_,imaxo_
                CHI(i,j,k) = CHI(i,j,k) + 2.0_WP*min(tmp1(i,j,k), 0.0_WP)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       call chemtable_lookup('SRC_ZMIX',tmp1)
       !$OMP PARALLEL DO
       do k=kmino_,kmaxo_
          do j=jmino_,jmaxo_
             do i=imino_,imaxo_
                CHI(i,j,k) = CHI(i,j,k) - 2.0_WP*min(tmp1(i,j,k), 0.0_WP)*SC(i,j,k,isc_ZMIX)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

!!$    ! If we use subfilter variance, include subfilter contribution in CHI
!!$    if (isc_ZVAR.ne.0) then
!!$       !$OMP PARALLEL DO
!!$       do k=kmino_,kmaxo_
!!$          do j=jmino_,jmaxo_
!!$             do i=imino_,imaxo_
!!$                CHI(i,j,k) = CHI(i,j,k) + 20.0_WP*(VISC(i,j,k)-VISCmol(i,j,k))/delta_3D(i,j)**2*ZVAR(i,j,k)
!!$             end do
!!$          end do
!!$       end do
!!$       !$OMP END PARALLEL DO
!!$    end if
    
    return
  end subroutine combustion_CHI
  
  ! Compute change of mass inside the domain
  subroutine combustion_sum_drho
    use parallel
    use masks
    implicit none
    real(WP) :: tmp
    integer :: i,j,k
    
    tmp = 0.0_WP
    !$OMP PARALLEL DO REDUCTION(+:tmp)
    do j=jmin_,jmax_
          do i=imin_,imax_
          if (mask(i,j).ne.1) then
             do k=kmin_,kmax_
                tmp = tmp + dRHO(i,j,k)*vol(i,j)
             end do
          end if
       end do
    end do
    !$OMP END PARALLEL DO
    
    call parallel_sum(tmp,sum_dRHO)
    
    return
  end subroutine combustion_sum_drho
  
  ! Compute mean density
  subroutine combustion_mean_density(RHOmean)
    use parallel
    use time_info
    implicit none
    
    real(WP), intent(out) :: RHOmean
    real(WP) :: my_RHOmean
    integer  :: i,j,k
    
    ! Recompute mean density
    my_RHOmean = 0.0_WP
    !$OMP PARALLEL DO REDUCTION(+:my_RHOmean)
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             my_RHOmean = my_RHOmean + RHO(i,j,k)*vol(i,j)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    call parallel_sum(my_RHOmean,RHOmean)
    RHOmean = RHOmean/vol_total
    
    return
  end subroutine combustion_mean_density
  
  ! Rescale density to ensure continuity
  subroutine combustion_rescale_density
    implicit none
    
    real(WP) :: RHOmean
    integer  :: i,j,k,isc
    
    ! If not constant volume return
    if (.not.is_constant_volume) return
    
    ! Recompute mean density
    call combustion_mean_density(RHOmean)
    
    ! Update Pthermo
    P_thermo = P_thermo*RHO_0/RHOmean
    
    !$OMP PARALLEL
    
    ! Update density
    !$OMP DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             RHO(i,j,k) = RHO(i,j,k) * (RHO_0/RHOmean)
          end do
       end do
    end do
    !$OMP END DO
    
    ! Update scalars
    do isc=1,nscalar
       !$OMP DO
       do k=kmino_,kmaxo_
          do j=jmino_,jmaxo_
             do i=imino_,imaxo_
                SC(i,j,k,isc) = SC(i,j,k,isc) * (RHOmean/RHO_0)
             end do
          end do
       end do
       !$OMP END DO
    end do

    !$OMP END PARALLEL

    return
  end subroutine combustion_rescale_density
  
  
  ! Interpolate the dRHO for large values of dRHO
  subroutine combustion_adv_dRHO
    use parallel
    use masks
    implicit none
    
    integer :: i,j,k,mycount
    real(WP) :: tmp1,tmp2,sum1,sum2,sdt
    
    tmp1 = 0.0_WP
    tmp2 = 0.0_WP
    !$OMP PARALLEL DO REDUCTION(+:tmp1,tmp2)
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             dRHO(i,j,k) = dRHO(i,j,k) * (1.0_WP-real(mask(i,j),WP))
             tmp1 = tmp1 + dRHO(i,j,k)
             tmp2 = tmp2 + dRHO(i,j,k)**2
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    call parallel_sum(tmp1,sum1)
    call parallel_sum(tmp2,sum2)
    sum1 = sum1/real(nx*ny*nz,WP)
    sum2 = sum2/real(nx*ny*nz,WP)
    sdt = sqrt(sum2-sum1**2)

    mycount = 0
    !$OMP PARALLEL DO PRIVATE(tmp1) REDUCTION(+:mycount)
    do k=kmin_,kmax_
       do j=jmin_,jmax_
          do i=imin_,imax_
             if (abs(dRHO(i,j,k)-sum1).gt.10.0_WP*sdt) then
                mycount = mycount + 1
                !dRHO(i,j,k) = &
                !     ( sum(vol(i-1:i+1,j-1:j+1)*sum(dRHO(i-1:i+1,j-1:j+1,k-1:k+1),dim=3)) &
                !     - vol(i,j)*dRHO(i,j,k)) / (3.0_WP*sum(vol(i-1:i+1,j-1:j+1))-vol(i,j))
                call filter_local_3D(dRHO(i-1:i+1,j-1:j+1,k-1:k+1),tmp1,i,j,'n')
                dRHO(i,j,k) = tmp1 * (1.0_WP-real(mask(i,j),WP))
                !dRHO(i,j,k) = &
                !     sum(vol(i-1:i+1,j-1:j+1)*sum(dRHO(i-1:i+1,j-1:j+1,k-1:k+1),dim=3)) &
                !     / (3.0_WP*sum(vol(i-1:i+1,j-1:j+1)))
             end if
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    call parallel_sum(mycount,clip_drho)
    
    ! Update boundaries
    call boundary_update_border(dRHO,'+','ym')
    
    return
  end subroutine combustion_adv_dRHO
  
end module combustion


! ============================================================ !
! Initialize the combustion module                             !
!                                                              ! 
! -> read chemtable                                            !
! -> initialize combustion model                               !
! -> allocate the arrays                                       !
! -> update ghost cells                                        !
! -> apply boundary conditions                                 !
!                                                              !
! Before: RHO or RHOold or dRHO correct only inside the domain !
!         -> imin_:imax_,jmin_:jmax_,kmin_:kmax_               !
! After : RHO and RHOold and dRHO correct everywhere           !
!         -> imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_         !
! ============================================================ !
subroutine combustion_init
  use combustion
  use memory
  use parser
  use borders
  implicit none
  
  integer :: isc,n,i,j,k,nmon
  character(len=str_medium) :: filename
  
  ! Create & Start the timer
  call timing_create('combustion')
  call timing_start ('combustion')

  ! Allocate some properties arrays
  allocate(delta_RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_))

  allocate(T   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Ggas(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Rgas(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))

  ! Link the variables to the given scalars
  isc_ZMIX  = 0
  isc_ZMIX2 = 0
  isc_ZVAR  = 0
  isc_PROG  = 0
  isc_PROG2  = 0
  isc_PROGVAR= 0
  isc_T     = 0
  isc_ENTH  = 0
  isc_E     = 0
  
  do isc=1,nscalar
     select case(trim(SC_name(isc)))
     case ('ZMIX')
        isc_ZMIX = isc
     case ('ZMIX2')
        isc_ZMIX2 = isc
     case ('ZVAR')
        isc_ZVAR = isc
     case ('PROG')
        isc_PROG = isc
     case ('PROG2')
        isc_PROG2 = isc
     case ('PROGV')
        isc_PROGVAR = isc
     case ('T')
        isc_T = isc
     case ('ENTH')
        isc_ENTH = isc
     case ('E')
        isc_E = isc
     case ('ZSTAR')
        isc_ZSTAR = isc
     case ('FMIX')
        isc_FMIX = isc
     end select
  end do

  ! Check for energy in compressible formulation
  if (compressible .and. isc_E.eq.0 .and. isc_T.eq.0) then
     call die('combustion_init: compressible solver requires E or T')
  end if

  ! Check for both enthalpy/temperature or internal energy/temperature equations
  if (isc_ENTH.ne.0 .and. isc_T.ne.0) then
     call die('combustion_init: cannot have both ENTH and T')
  end if
  if (isc_E.ne.0 .and. isc_T.ne.0) then
     call die('combustion_init: cannot have both E and T')
  end if
  
  ! Set the number of filtering to zero
  nfilter = 0
  clip_drho = 0

  ! Use mixture fraction variance
  call parser_read('Mixture fraction variance',use_ZVAR,.false.)
  if (.not.use_ZVAR .and. isc_ZMIX.gt.0 .and. use_sgs) then
     if (irank.eq.iroot) print *, "combustion_init: Not using ZVAR but have ZMIX and SGS model. This is unusual; are you sure?"
  end if
  if (.not.use_ZVAR .and. (isc_ZVAR.gt.0 .or. isc_ZMIX2.gt.0)) then
     call die('combustion_init: Cannot solve for ZVAR or ZMIX2 and exclude mixture fraction variance')
  end if

  ! Use progress variable variance
  call parser_read('Progress variable variance',use_PROGVAR,.false.)
  if (.not.use_PROGVAR .and. (isc_PROGVAR.gt.0 .or. isc_PROG2.gt.0)) then
     call die('combustion_init: Cannot solve for PROGVAR or PROG2 and exclude prog var variance')
  end if

  ! Use FMIX
  if ((isc_ZSTAR.gt.0 .or. isc_FMIX.gt.0)) then 
     use_FMIX = .true.
  else 
     use_FMIX = .false.
  end if

  ! Subfilter modeling for FMIX
  call parser_read('Subfilter FMIX source', use_FMIXsub,.true.)
  if ( (isc_FMIX.ne.0) .and. use_sgs .and. use_FMIXsub .and. (isc_ZSTAR.eq.0) ) then
     call die('combustion_init: cannot model subfilter source term for FMIX unless ZSTAR is present')
  end if
  if (isc_FMIX.ne.0) call parser_read('Min Z for F src term', tolZ)

  ! Allocate more arrays
  allocate(CHI (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  if (use_ZVAR) allocate(ZVAR(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  if (use_PROGVAR) allocate(PROGVAR(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  if (use_FMIX) allocate(FMIX(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))

  ! Prepare the place to store properties
  n = 1
  if (xper.eq.1) n = ninlet
  
  ! Allocate the arrays
  allocate(rho_input(n))
  allocate(visc_input(n))
  allocate(diff_input(n))
  allocate(T_input(n))

  ! Read the chemistry model from the input file
  ! -> none : constant properties
  ! -> chemtable : read every quantity from a file
  ! -> finite chemt : solve species transport equations
  call parser_read('Chemistry model',chemistry)
  select case(trim(chemistry))
     
  case ('none') ! *** CONSTANT PROPERTIES ***
     combust = .false.
     
     ! Read default values
     call parser_read('Viscosity',visc_input)
     call parser_read('Diffusivity',diff_input)
     if (.not.compressible) then
        call parser_read('Density',rho_input)
        call parser_read('Temperature',T_input)
        call parser_read('Pressure',P_input,1.0132e5_WP)
     end if

     if (compressible) then
        call parser_read('Molar mass',MM_input,0.0289_WP)
        Rgas = R_uni / MM_input
        call parser_read('Specific heat ratio',GAMMA_input,1.4_WP)
        Ggas = GAMMA_input

        ! Look for density
        if (.not.rho_present) then
           call die('combustion_init: compressible solver requires rho')
        end if
     end if
     
     ! Set RHO, VISC and DIFF
     call combustion_density
     call combustion_viscosity
     call combustion_diffusivity

     ! Update ghost cells
     if (isc_E.gt.0) call boundary_update_border(SC(:,:,:,isc_E),'+','ym')
     if (isc_T.gt.0) call boundary_update_border(SC(:,:,:,isc_T),'+','ym')

  case ('chemtable') ! *** CHEMTABLE ***
     combust = .true.
     
     ! Read default values
     call parser_read('Pressure',P_input,1.0133e5_WP)

     ! Read in the chemtable and initialize it
     call parser_read('Chemtable file',filename)
     call chemtable_init(filename)

     ! VISC/DIFF must be in data file
     if (.not.visc_present .or. .not.diff_present) then
        call die('combustion_init: for chemtable, VISC and DIFF must be in data')
     end if
     
     ! Read number of filtering from input
     call parser_read('dRHO filtering',nfilter,0)

     ! Read whether or not to clip ZMIX
     call parser_read('Z clipping',clip_zmix,.false.)

  case ('finite chem') ! *** FINITE RATE CHEMISTRY ***
     combust = .true.
     
     ! Check we have enough to use the model
     if (isc_T .eq.0 .or. .not.rho_present) then
        call die('combustion_init: finite chem model requires T and RHO')
     end if

     ! Initialize the finite rate chemistry module
     call finitechem_init

     ! Set VISC and DIFF
     call combustion_viscosity
     call combustion_diffusivity

     ! Get the thermodynamic pressure
     call parser_read('Pressure',P_input)

  case default
     call die('combustion_init: unknown chemistry model')
  end select
  
  ! Set CHI or variances
  if (isc_ZMIX.gt.0) call combustion_CHI
  if (use_ZVAR)    call combustion_update_variance(SC,ZVAR   ,isc_ZMIX,isc_ZMIX2,isc_ZVAR)
  if (use_PROGVAR) call combustion_update_variance(SC,PROGVAR,isc_PROG,isc_PROG2,isc_PROGVAR)
  if (use_FMIX) call combustion_FMIX(SC)

  ! Detect if closed/constant volume
  if (xper.eq.1) then
     is_constant_volume = .true.
     call combustion_mean_density(RHO_0)
  else
     is_constant_volume = .false.
  end if
  
  ! Initialize thermodynamic pressure
  P_thermo = P_input
  P_thermo_old = P_thermo
  
  ! If the arrays were not present in the data file, update them
  ! -- If present, evaluate the density in the boundary points
  if (.not. rho_present)  then
     call combustion_density
  else
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              tmp1(i,j,k) = RHO(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call combustion_density
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              RHO(i,j,k) = tmp1(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  if (.not. drho_present) then
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              dRHO(i,j,k) = 0.0_WP
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  ! Update all the borders
  call boundary_update_border(RHO ,'+','ym')
  call boundary_update_border(dRHO,'+','ym')
  
  ! Recompute the old density
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           RHOold(i,j,k) = RHO(i,j,k) - dRHO(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Compute the initial temperature
  call combustion_temperature

  ! Compute the initial pressure
  if (compressible) call combustion_pressure

  if (combust .or. compressible) then
     ! Create a new file to monitor at each timesteps
     nmon = 6
     if (compressible) nmon = nmon+1
     call monitor_create_file_step('combustion',nmon)
     call monitor_set_header (1,'min_T','r')
     call monitor_set_header (2,'max_T','r')
     call monitor_set_header (3,'min_rho','r')
     call monitor_set_header (4,'max_rho','r')
     call monitor_set_header (5,'sum_dRHO','r')
     if (.not.compressible) then
        call monitor_set_header (6,'P_thermo','r')
     else
        call monitor_set_header (6,'min_P','r')
        call monitor_set_header (7,'max_P','r')
     end if

     ! Create a new file to monitor at each subiteration
     call monitor_create_file_iter('convergence_combustion',2)
     call monitor_set_header (1,'res_RHO','r')
     call monitor_set_header (2,'clip_dRHO','i')
  end if

  ! Stop a timer
  call timing_stop('combustion')
  
  return
end subroutine combustion_init


! ==================================================== !
! Precompute some quantities for the combustion module !
! ==================================================== !
subroutine combustion_prestep
  use combustion
  implicit none

  integer :: i,j,k
  
  ! Start a timer
  call timing_start('combustion')
  
  ! Compute Viscosity, Diffusivity, and Temperature
  ! -- VISC/DIFF need to be in subiterations but synced with SGS
  call combustion_viscosity
!!$  call combustion_diffusivity
  call combustion_temperature
  
  if (clip_zmix) call combustion_clip_z

  if (nscalar.ne.0) then
     ! Compute ZVAR and CHI
     if (use_ZVAR)    call combustion_update_variance(SC,ZVAR   ,isc_ZMIX,isc_ZMIX2,isc_ZVAR)
     if (use_PROGVAR) call combustion_update_variance(SC,PROGVAR,isc_PROG,isc_PROG2,isc_PROGVAR)
     if (use_FMIX) call combustion_FMIX(SC)
     if (isc_ZMIX.ne.0) call combustion_CHI
  end if
  
  ! Save the old density
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           RHOold(i,j,k) = RHO(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  ! Save the old Pthermo
  P_thermo_old = P_thermo

  ! Save the old pressure
  if (compressible) then
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              Pold(i,j,k) = P(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  ! Stop a timer
  call timing_stop('combustion')
  
  return
end subroutine combustion_prestep


! =========================================== !
! Compute the main quantities for the solvers !
! And some monitoring informations            !
! =========================================== !
subroutine combustion_step
  use combustion
  use masks
  use parallel
  use scalar
  implicit none
  integer  :: i,j,k,isc
  real(WP) :: max_dRHO
  
  ! Start a timer
  call timing_start('combustion')

  if (.not.compressible) then
     ! Save current density
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              delta_RHO(i,j,k) = - RHO(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     
     ! Compute the new density
     call combustion_invert_density

     ! Rescale the pressure to ensure continuity
     call combustion_rescale_density

     ! Compute change in RHO
     if (nscalar.ne.0) then
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imino_,imaxo_
                 !SC(i,j,k,:) = SC(i,j,k,:)*RHO(i,j,k)
                 dRHO(i,j,k) = RHO(i,j,k) - RHOold(i,j,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
        do i=1,nfilter
           call filter_global_3D(dRHO,dRHO,'+','n')
        end do
        if (nfilter.eq.-1) call combustion_adv_dRHO
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imino_,imaxo_
                 RHO(i,j,k) = RHOold(i,j,k) + dRHO(i,j,k)
                 !SC(i,j,k,:) = SC(i,j,k,:)/RHO(i,j,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
     
     !$OMP PARALLEL
     
     ! Compute change of RHO over iteration
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              delta_RHO(i,j,k) = delta_RHO(i,j,k) + RHO(i,j,k)
           end do
        end do
     end do
     !$OMP END DO
     
     ! Set zero in the walls
     !$OMP DO
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           if (mask(i,j).ne.0) dRHO(i,j,:) = 0.0_WP
        end do
     end do
     !$OMP END DO
     
     !$OMP END PARALLEL

     call combustion_sum_drho
  else
     call density_step
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              if (mask(i,j).ne.0) then
                 dRHO(i,j,k) = 0.0_WP
              else
                 dRHO(i,j,k) = RHO(i,j,k) - RHOold(i,j,k)
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO

     ! Save the change in density
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              delta_RHO(i,j,k) = dRHO(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  ! Compute temperature for output
  call combustion_temperature

!!$  ! Compute diffusivity for transport (mid-point)
!!$  !$OMP PARALLEL
!!$  do isc=1,nscalar
!!$     !$OMP DO
!!$     do k=kmino_,kmaxo_
!!$        do j=jmino_,jmaxo_
!!$           do i=imino_,imaxo_
!!$              SC(i,j,k,isc) = 0.5_WP * (SC(i,j,k,isc) + SCold(i,j,k,isc))
!!$           end do
!!$        end do
!!$     end do
!!$     !$OMP END DO
!!$  end do
!!$  !$OMP END PARALLEL
!!$  call combustion_diffusivity
!!$  !$OMP PARALLEL
!!$  do isc=1,nscalar
!!$     !$OMP DO
!!$     do k=kmino_,kmaxo_
!!$        do j=jmino_,jmaxo_
!!$           do i=imino_,imaxo_
!!$              SC(i,j,k,isc) = 2.0_WP * SC(i,j,k,isc) - SCold(i,j,k,isc)
!!$           end do
!!$        end do
!!$     end do
!!$     !$OMP END DO
!!$  end do
!!$  !$OMP END PARALLEL

  ! Compute pressure if compressible
  if (compressible) call combustion_pressure

  ! Compute the heat release rate
  if (trim(chemistry).eq.'finite chem') call finitechem_HR

  ! Set the monitor values
  if (combust .or. compressible) then
     call parallel_max(maxval(abs(delta_RHO)),max_dRHO)
     call monitor_select_file('convergence_combustion')
     call monitor_set_single_value(1,max_dRHO)
     call monitor_set_single_value(2,real(clip_dRHO,WP))
  end if

  ! Stop a timer
  call timing_stop('combustion')

  return
end subroutine combustion_step


! ======================================================= !
! Source terms for scalar equations evaluated at midpoint !
! ======================================================= !
subroutine combustion_source_scalar(srcSC)
  use combustion
  use filter
  use memory
  use metric_velocity_conv
  use strainrate
  use masks
  implicit none

  integer :: i,j,k
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(inout) :: srcSC

  ! Start a timer
  call timing_start('combustion')

  ! Pressure-divergence
  if (compressible .or. is_constant_volume) then
     if (isc_E.gt.0) then
        !$OMP PARALLEL DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 srcSC(i,j,k,isc_E) = srcSC(i,j,k,isc_E) - &
                      dt_*0.5_WP*(Pold(i,j,k)+P(i,j,k))*( &
                      sum(divc_u(i,j,:)*U(i-stc1:i+stc2,j,k)) + &
                      sum(divc_v(i,j,:)*V(i,j-stc1:j+stc2,k)) + &
                      sum(divc_w(i,j,:)*W(i,j,k-stc1:k+stc2)))
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
     if (isc_T.gt.0) then
        if (compressible) then
           !$OMP PARALLEL DO
           do k=kmin_,kmax_
              do j=jmin_,jmax_
                 do i=imin_,imax_
                    if (mask(i,j).ne.1) then
                       srcSC(i,j,k,isc_T) = srcSC(i,j,k,isc_T) - &
                            dt_*(Ggas(i,j,k)-1.0_WP)/Rgas(i,j,k)* &
                            0.5_WP*(Pold(i,j,k)+P(i,j,k))*( &
                            sum(divc_u(i,j,:)*U(i-stc1:i+stc2,j,k)) + &
                            sum(divc_v(i,j,:)*V(i,j-stc1:j+stc2,k)) + &
                            sum(divc_w(i,j,:)*W(i,j,k-stc1:k+stc2)))
                    end if
                 end do
              end do
           end do
           !$OMP END PARALLEL DO
        else
           call finitechem_pressuredivergence(srcSC,SC)
        end if
     end if
  end if

  ! Viscous dissipation
  if (compressible) then
     if (isc_E.gt.0) then
        call strainrate_compute(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
        !$OMP PARALLEL DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 srcSC(i,j,k,isc_E) = srcSC(i,j,k,isc_E) + dt_*2.0_WP*VISCmol(i,j,k)*S(i,j,k)**2
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
     if (isc_T.gt.0) then  
        call strainrate_compute(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
        !$OMP PARALLEL DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 if (mask(i,j).ne.1) then
                    srcSC(i,j,k,isc_T) = srcSC(i,j,k,isc_T) + &
                         dt_*(Ggas(i,j,k)-1.0_WP)/Rgas(i,j,k)*2.0_WP*VISCmol(i,j,k)*S(i,j,k)**2
                 end if
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  end if

  ! Mixture fraction squared: dissipation
  ! Subfilter dissipation model from Mueller PhD thesis section 4.2.1
  if (isc_ZMIX2.ne.0 .and. use_sgs) then
     call gradient_squared(SC(:,:,:,isc_ZMIX),tmp1)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              srcSC(i,j,k,isc_ZMIX2) = srcSC(i,j,k,isc_ZMIX2) &
                   -dt_*2.0_WP*DIFFmol(i,j,k,isc_ZMIX)*tmp1(i,j,k) &
                   -dt_*20.0_WP*(VISC(i,j,k)-VISCmol(i,j,k))/delta_3D(i,j)**2*ZVAR(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  ! Subfilter variance: production and dissipation
  ! Subfilter dissipation model from Mueller PhD thesis section 4.2.1
  if (isc_ZVAR.ne.0 .and. use_sgs) then
     call gradient_squared(SC(:,:,:,isc_ZMIX),tmp1)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              srcSC(i,j,k,isc_ZVAR) = srcSC(i,j,k,isc_ZVAR) &
                   +dt_*2.0_WP*(DIFF(i,j,k,isc_ZMIX)-DIFFmol(i,j,k,isc_ZMIX))*tmp1(i,j,k) &
                   -dt_*20.0_WP*(VISC(i,j,k)-VISCmol(i,j,k))/delta_3D(i,j)**2*ZVAR(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  ! Chemtable source terms
  if (trim(chemistry).eq.'chemtable') call chemtable_source(srcSC)
  
  ! Progress variable squared: dissipation and source term ~(PROG*SRC_PROG)
  ! Subfilter dissipation model from Mueller PhD thesis section 4.2.1
  if (isc_PROG2.ne.0 .and. use_sgs) then
     call gradient_squared(SC(:,:,:,isc_PROG),tmp1)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              srcSC(i,j,k,isc_PROG2) = srcSC(i,j,k,isc_PROG2) &
                   -dt_*2.0_WP*DIFFmol(i,j,k,isc_PROG)*tmp1(i,j,k) &
                   -dt_*20.0_WP*(VISC(i,j,k)-VISCmol(i,j,k))/delta_3D(i,j)**2 *PROGVAR(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  ! Subfilter variance: production and dissipation
  ! Subfilter dissipation model from Mueller PhD thesis section 4.2.1
  if (isc_PROGVAR.ne.0 .and. use_sgs) then
     call gradient_squared(SC(:,:,:,isc_PROG),tmp1)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              srcSC(i,j,k,isc_PROGVAR) = srcSC(i,j,k,isc_PROGVAR) &
                   +dt_*2.0_WP*(DIFF(i,j,k,isc_PROG) - DIFFmol(i,j,k,isc_PROG))*tmp1(i,j,k) &
                   -dt_*20.0_WP*(VISC(i,j,k)-VISCmol(i,j,k))/delta_3D(i,j)**2 *PROGVAR(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  ! FMIX: Cross-dissipation term in transport equation
  if (isc_FMIX.ne.0) then
     call gradient_dotproduct( SC(:,:,:,isc_ZMIX), SC(:,:,:,isc_FMIX), tmp1)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_

              ! Resolved Part 
              ! don't compute for Z near 0, would require dividing by 0
              if (SC(i,j,k,isc_ZMIX) .ge. tolZ) then
                 srcSC(i,j,k,isc_FMIX) = srcSC(i,j,k,isc_FMIX) &
                      + dt_ * (2.0_WP / SC(i,j,k,isc_ZMIX) * tmp1(i,j,k) &
                            * DIFFmol(i,j,k,isc_ZSTAR) )
              end if
              ! Modeled subgrid part
              ! don't compute for Z near 0, would require dividing by 0
              if (use_sgs .and. use_FMIXsub) then
                 if (SC(i,j,k,isc_ZMIX) .ge. tolZ) then
                    srcSC(i,j,k,isc_FMIX) = srcSC(i,j,k,isc_FMIX) &
                         + dt_ * 20.0_WP * (VISC(i,j,k)-VISCmol(i,j,k))/delta_3D(i,j)**2 &
                               * (SC(i,j,k,isc_ZSTAR) - SC(i,j,k,isc_FMIX)*SC(i,j,k,isc_ZMIX)) &
                               / SC(i,j,k,isc_ZMIX) 
                 end if
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  
  ! Finitechem transport sources
  if (trim(chemistry).eq.'finite chem') call finitechem_source_transport(srcSC) 

  ! Stop a timer
  call timing_stop('combustion')

  return
end subroutine combustion_source_scalar


! =========================================================== !
! Compute chemistry source terms for Strang split solver      !
! =========================================================== !
subroutine combustion_source_scalar_full(RHO_,SC_)
  use combustion
  implicit none

  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(inout) :: RHO_
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(inout) :: SC_
  integer :: i,j,k

  ! Start a timer
  call timing_start('combustion')

  ! Chemical source terms
  if (trim(chemistry).eq.'finite chem') then
     call finitechem_source_chemistry(RHO_,SC_) 
  end if

  ! Stop a timer
  call timing_stop('combustion')

  return
end subroutine combustion_source_scalar_full


! =========================================================== !
! Compute all relevant scalars by inverting : rhoZ (and rhoC) !
! =========================================================== !
subroutine combustion_invert_density
  use combustion
  use borders
  implicit none
  integer :: i,j,k,nflow

  select case(trim(chemistry))
  case ('none')
     ! Constant density from input file
     if (xper.eq.1) then
        do nflow=1,ninlet
           do j=max(inlet(nflow)%jmino,jmino_),min(inlet(nflow)%jmaxo,jmaxo_)
              !$OMP PARALLEL DO
              do k=kmino_,kmaxo_
                 do i=imino_,imaxo_
                    RHO(i,j,k) = rho_input(nflow)
                 end do
              end do
              !$OMP END PARALLEL DO
           end do
        end do
     else
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imino_,imaxo_
                 RHO(i,j,k) = rho_input(1)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  case ('chemtable')
     ! Invert rhoSC to get density
     call chemtable_invert_density
  case ('finite chem')
     ! Get density from the equation of state
     call finitechem_density
  end select
  
  
  return
end subroutine combustion_invert_density


! ================================ !
! Compute the Density everywhere   !
! Including the ghost cells        !
! --Nothing to do for compressible !
! ================================ !
subroutine combustion_density
  use combustion
  use borders
  implicit none
  integer :: i,j,k,nflow

  select case(trim(chemistry))
  case ('none')
     if (.not.compressible) then
        ! Constant density from input file
        if (xper.eq.1) then
           do nflow=1,ninlet
              do j=max(inlet(nflow)%jmino,jmino_),min(inlet(nflow)%jmaxo,jmaxo_)
                 !$OMP PARALLEL DO
                 do k=kmino_,kmaxo_
                    do i=imino_,imaxo_
                       RHO(i,j,k) = rho_input(nflow)
                    end do
                 end do
                 !$OMP END PARALLEL DO
              end do
           end do
        else
           !$OMP PARALLEL DO
           do k=kmino_,kmaxo_
              do j=jmino_,jmaxo_
                 do i=imino_,imaxo_
                    RHO(i,j,k) = rho_input(1)
                 end do
              end do
           end do
           !$OMP END PARALLEL DO
        end if
     end if
  case ('chemtable')
     ! Get density from chemtable
     call chemtable_lookup('RHO',RHO)
  case ('finite chem')
     ! Get density from the equation of state
     call finitechem_density  
  end select
  
  return
end subroutine combustion_density


! ================================ !
! Compute the Viscosity everywhere !
! Including the ghost cells        !
! ================================ !
subroutine combustion_viscosity
  use combustion
  use borders
  implicit none
  integer :: i,j,k,nflow
  
  select case(trim(chemistry))
  case ('none')
     ! Constant viscosity from input file
     if (xper.eq.1) then
        do nflow=1,ninlet
           do j=max(inlet(nflow)%jmino,jmino_),min(inlet(nflow)%jmaxo,jmaxo_)
              !$OMP PARALLEL DO
              do k=kmino_,kmaxo_
                 do i=imino_,imaxo_
                    VISC(i,j,k) = visc_input(nflow)
                 end do
              end do
              !$OMP END PARALLEL DO
           end do
        end do
     else
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imino_,imaxo_
                 VISC(i,j,k) = visc_input(1)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  case ('chemtable')
     ! Get viscosity from chemtable
     call chemtable_lookup('VISC',VISC)
  case ('finite chem')
     ! Get viscosity from species transport properties
     call finitechem_viscosity          
  end select
  
  ! Store molecular viscosity
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           VISCmol(i,j,k) = VISC(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine combustion_viscosity


! ================================== !
! Compute the Diffusivity everywhere !
! Including the ghost cells          !
! ================================== !
subroutine combustion_diffusivity
  use combustion
  use memory
  use borders
  implicit none
  integer :: isc,nflow
  integer :: i,j,k
  
  ! Return if no scalar
  if (nscalar.eq.0) return

  select case(trim(chemistry))
  case ('none')
     ! Constant diffusivity from input file
     if (xper.eq.1) then
        !$OMP PARALLEL
        do isc=1,nscalar
           do nflow=1,ninlet
              do j=max(inlet(nflow)%jmino,jmino_),min(inlet(nflow)%jmaxo,jmaxo_)
                 !$OMP DO
                 do k=kmino_,kmaxo_
                    do i=imino_,imaxo_
                       DIFF(i,j,k,isc) = diff_input(nflow)
                    end do
                 end do
                 !$OMP END DO
              end do
           end do
        end do
        !$OMP END PARALLEL
     else
        !$OMP PARALLEL
        do isc=1,nscalar
           !$OMP DO
           do k=kmino_,kmaxo_
              do j=jmino_,jmaxo_
                 do i=imino_,imaxo_
                    DIFF(i,j,k,isc) = diff_input(1)
                 end do
              end do
           end do
           !$OMP END DO
        end do
        !$OMP END PARALLEL
     end if
     ! Modified diffusivity for internal energy equation
     if (isc_E.gt.0) then
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imino_,imaxo_
                 if (mask(i,j).ne.1) DIFF(i,j,k,isc_E) = Ggas(i,j,k)*DIFF(i,j,k,isc_E)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
     if ((compressible .or. is_constant_volume) .and. isc_T.gt.0) then
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imino_,imaxo_
                 if (mask(i,j).ne.1) DIFF(i,j,k,isc_T) = Ggas(i,j,k)*DIFF(i,j,k,isc_T)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
  case ('chemtable')
     ! Get diffusivity from chemtable
     call chemtable_lookup('DIFF',tmp1)
     ! Same diffusivity for all scalars
     ! -- Overwritten for soot, etc.
     do isc=1,nscalar
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           do j=jmino_,jmaxo_
              do i=imino_,imaxo_
                 DIFF(i,j,k,isc) = tmp1(i,j,k)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  case ('finite chem') 
     ! Get diffusivity from species transport properties
     call finitechem_diffusivity   
  end select

  !$OMP PARALLEL

  ! Store molecular diffusivity
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              DIFFmol(i,j,k,isc)=DIFF(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do

  !$OMP END PARALLEL
  
  return
end subroutine combustion_diffusivity


! ===================================== !
! Compute the Temperature by everywhere !
! Including the ghost cells             !
! ===================================== !
subroutine combustion_temperature
  use combustion
  use borders
  implicit none
  integer :: i,j,k,nflow
  
  select case(trim(chemistry))
  case ('none')
     if (.not.compressible) then
        ! Constant temperature from input file
        if (xper.eq.1) then
           do nflow=1,ninlet
              do j=max(inlet(nflow)%jmino,jmino_),min(inlet(nflow)%jmaxo,jmaxo_)
                 !$OMP PARALLEL DO
                 do k=kmino_,kmaxo_
                    do i=imino_,imaxo_
                       T(i,j,k) = T_input(nflow)
                    end do
                 end do
                 !$OMP END PARALLEL DO
              end do
           end do
        else
           !$OMP PARALLEL DO
           do k=kmino_,kmaxo_
              do j=jmino_,jmaxo_
                 do i=imino_,imaxo_
                    T(i,j,k) = T_input(1)
                 end do
              end do
           end do
           !$OMP END PARALLEL DO
        end if
     else
        ! Get temperature from internal energy or temperature itself
        if (isc_E.gt.0) then
           !$OMP PARALLEL DO
           do k=kmino_,kmaxo_
              do j=jmino_,jmaxo_
                 do i=imino_,imaxo_
                    T(i,j,k) = (Ggas(i,j,k)-1.0_WP)/Rgas(i,j,k)*SC(i,j,k,isc_E)
                 end do
              end do
           end do
           !$OMP END PARALLEL DO
        else
           !$OMP PARALLEL DO
           do k=kmino_,kmaxo_
              do j=jmino_,jmaxo_
                 do i=imino_,imaxo_
                    T(i,j,k) = SC(i,j,k,isc_T)
                 end do
              end do
           end do
           !$OMP END PARALLEL DO
        end if
     end if
  case ('chemtable')
     ! Get temperature from chemtable table
     call chemtable_lookup('T',T)
  case ('finite chem')
     ! Get temperature from temperature scalar
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              T(i,j,k) = SC(i,j,k,isc_T)   
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end select
  
  return
end subroutine combustion_temperature


! =============================================================== !
! Compute the pressure with ideal gas law for compressible solver !
! =============================================================== !
subroutine combustion_pressure
  use combustion
  use borders
  implicit none
  integer :: i,j,k
  
  ! Pressure from ideal gas law
  !$OMP PARALLEL DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           P(i,j,k) = RHO(i,j,k)*Rgas(i,j,k)*T(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine combustion_pressure


! ================================ !
! Compute the local speed of sound !
! ================================ !
subroutine combustion_soundspeed(a)
  use combustion
  implicit none

  real(WP), intent(inout), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: a
  integer :: i,j,k

  select case(trim(chemistry))
  case ('none')
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              a(i,j,k) = sqrt( Ggas(i,j,k) * P(i,j,k) / RHO(i,j,k) )
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end select

  return
end subroutine combustion_soundspeed


! ================================ !
! Compute the local speed of sound !
! ================================ !
subroutine combustion_soundspeed_old(a)
  use combustion
  implicit none

  real(WP), intent(inout), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: a
  integer :: i,j,k

  select case(trim(chemistry))
  case ('none')
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              a(i,j,k) = sqrt( Ggas(i,j,k) * Pold(i,j,k) / RHOold(i,j,k) )
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end select

  return
end subroutine combustion_soundspeed_old


! ========================== !
! Ignite the flow everywhere !
! ========================== !
subroutine combustion_ignite
  use combustion
  implicit none
  
  integer :: i,j,k
!!$  real(WP) :: sc_min,sc_max
  
  ! Only for tabulated chemistry
!  if (trim(chemistry).ne.'chemtable') return

!!$  ! Get the most burning flamelet
!!$  if (trim(combModel).eq.'FPVA') then
!!$     do k=kmino_,kmaxo_
!!$        do j=jmino_,jmaxo_
!!$           do i=imino_,imaxo_
!!$              call chemtable_prog_minmax(SC(i,j,k,isc_ZMIX),ZVAR(i,j,k),sc_min,sc_max)
!!$              SC(i,j,k,isc_PROG) = sc_max
!!$           end do
!!$        end do
!!$     end do
!!$  end if

!!$  ! Get the most burning, zero heat loss flamelet
!!$  if (trim(combModel).eq.'RFPVA') then
!!$     do k=kmino_,kmaxo_
!!$        do j=jmino_,jmaxo_
!!$           do i=imino_,imaxo_
!!$              call unsteady_chemtable_enth_minmax(SC(i,j,k,isc_ZMIX), ZVAR(i,j,k),sc_min,sc_max)
!!$              SC(i,j,k,isc_ENTH) = sc_max
!!$              call unsteady_chemtable_prog_minmax(SC(i,j,k,isc_ZMIX),ZVAR(i,j,k),SC(i,j,k,isc_ENTH),sc_min,sc_max)
!!$              SC(i,j,k,isc_PROG) = sc_max
!!$           end do
!!$        end do
!!$     end do
!!$  end if
!!$
!!$  ! To go from ZVAR to ZMIX2 via renaming
!!$  do k=kmino_,kmaxo_
!!$     do j=jmino_,jmaxo_
!!$        do i=imino_,imaxo_
!!$           SC(i,j,k,isc_ZMIX2) = SC(i,j,k,isc_ZMIX)*SC(i,j,k,isc_ZMIX)
!!$        end do
!!$     end do
!!$  end do


!!$  ! Clips mixture fraction in case of spurious generation
!!$  do k=kmino_,kmaxo_
!!$     do j=jmino_,jmaxo_
!!$        do i=imino_,imaxo_
!!$           SC(i,j,k,isc_ZMIX) = min(SC(i,j,k,isc_ZMIX), 0.05515383_WP)
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  ! Initializes enthalpy properly after adding H to data file
!!$  ! COMMENT OUT WHEN NOT IN USE
!!$  if (isc_ENTH.ne.0) then
!!$     do k=kmino_,kmaxo_
!!$        do j=jmino_,jmaxo_
!!$           do i=imino_,imaxo_
!!$              SC(i,j,k,isc_ENTH) = -4.4332e6_WP*SC(i,j,k,isc_ZMIX) + 1.547436e5_WP
!!$           end do
!!$        end do
!!$     end do
!!$  end if

  ! Recompute rho
  call combustion_density
  
  return
end subroutine combustion_ignite


! ============================= !
! Monitor the combustion module !
! ============================= !
subroutine combustion_monitor
  use combustion
  use masks
  implicit none
  
  integer  :: i,j,k
  real(WP) :: min_rho_,max_rho_,min_T_,max_T_,min_P_,max_P_
  
  ! Nothing to do if constant density
  if (.not.combust .and. .not.compressible) return
  
  ! Start a timer
  call timing_start('combustion')
  
  ! Min and max of rho and T
  min_rho_ = +huge(1.0_WP)
  max_rho_ = -huge(1.0_WP)
  min_T_   = +huge(1.0_WP)
  max_T_   = -huge(1.0_WP)

  !$OMP PARALLEL DO REDUCTION(min:min_rho_,min_T_) REDUCTION(max:max_rho_,max_T_)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).ne.1) then
              if (RHO(i,j,k).lt.min_rho_) min_rho_ = RHO(i,j,k)
              if (RHO(i,j,k).gt.max_rho_) max_rho_ = RHO(i,j,k)
              if (T(i,j,k)  .lt.min_T_)   min_T_   = T(i,j,k)
              if (T(i,j,k)  .gt.max_T_)   max_T_   = T(i,j,k)
           end if
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call parallel_max(-min_rho_,min_rho)
  min_rho = -min_rho
  call parallel_max( max_rho_,max_rho)
  call parallel_max(-min_T_,min_T)
  min_T = - min_T
  call parallel_max( max_T_,max_T)

  ! Min and max of P
  if (compressible) then
     min_P_ = +huge(1.0_WP)
     max_P_ = -huge(1.0_WP)
     
     !$OMP PARALLEL DO REDUCTION(min:min_P_) REDUCTION(max:max_P_)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              if (mask(i,j).ne.1) then
                 if (P(i,j,k)  .lt.min_P_)   min_P_   = P(i,j,k)
                 if (P(i,j,k)  .gt.max_P_)   max_P_   = P(i,j,k)
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call parallel_max(-min_P_,min_P)
     min_P = -min_P
     call parallel_max( max_P_,max_P)
  end if
  
  ! Transfer values to monitor
  call monitor_select_file('combustion')
  call monitor_set_single_value(1,min_T)
  call monitor_set_single_value(2,max_T)
  call monitor_set_single_value(3,min_rho)
  call monitor_set_single_value(4,max_rho)
  call monitor_set_single_value(5,sum_dRHO)
  if (.not.compressible) then
     call monitor_set_single_value(6,P_thermo)
  else
     call monitor_set_single_value(6,min_P)
     call monitor_set_single_value(7,max_P)
  end if

  ! Stop a timer
  call timing_stop('combustion')
  
  return
end subroutine combustion_monitor

!==============================!
!         Clip Z mix           !
!==============================!

subroutine combustion_clip_z
  use combustion
  implicit none

  integer :: i,j,k
  
  ! CLIPPING
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           SC(i,j,k,isc_ZMIX) = max(0.0_WP,min(1.0_WP, SC(i,j,k,isc_ZMIX)))
           SC(i,j,k,isc_ZMIX2) = max(0.0_WP,min(1.0_WP, SC(i,j,k,isc_ZMIX2)))
        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine combustion_clip_z


!==============================!
! Recalculate Ctable Variables !
!==============================!

subroutine combustion_poststep
  use combustion
  implicit none

  ! Start a timer
  call timing_start('combustion')
  
  if (trim(chemistry).eq.'chemtable') then
     if (isc_ZMIX.gt.0) call combustion_CHI
     if (use_ZVAR)    call combustion_update_variance(SC,ZVAR   ,isc_ZMIX,isc_ZMIX2,isc_ZVAR)
     if (use_PROGVAR) call combustion_update_variance(SC,PROGVAR,isc_PROG,isc_PROG2,isc_PROGVAR)
     if (use_FMIX)      call combustion_FMIX(SC)
  end if

  ! Stop a timer
  call timing_stop('combustion')

end subroutine combustion_poststep
