module stat_2d
  use stat
  implicit none

  ! Sampled variables
  integer :: stat_nvar, nstat_scalar
  character(len=str_medium), dimension(:), pointer :: stat_name
  character(len=str_short),  dimension(:), pointer :: stat_scalar

  ! 2D Statistics
  real(WP), dimension(:,:,:), pointer :: stat_xy
  real(WP), dimension(:,:,:), pointer :: buf
  real(WP) :: Delta_t
  
  ! Fileview
  integer :: fileview

  logical :: use_chi_z
  logical :: use_chi_c
  logical :: give_ksgs
  
contains
  
  subroutine stat_2d_init_names
    use data
    use combustion
    use string
    implicit none
    
    integer :: isc,ns
    character(len=str_short) :: name
    
    ! ================================================================================================ !
    ! Count the variables
    stat_nvar = 0
    ! Density 
    if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+2
    ! Velocity
    stat_nvar = stat_nvar+7
    if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+7
    ! Scalars
    do isc=1,nscalar
       name = SC_name(isc)
       if (name(1:2).ne.'S_') then
          stat_nvar = stat_nvar+4
          if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+4
       end if
    end do
    ! Scalar Variance
    if (isc_ZVAR.eq.0 .and. isc_ZMIX2.eq.0 .and. use_ZVAR) then
       stat_nvar = stat_nvar+2
       if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+2
    end if
    ! Scalar Dissipation Rate
    if (isc_ZMIX.ne.0) then
       stat_nvar = stat_nvar+2
       if (trim(chemistry).ne.'none') stat_nvar = stat_nvar+2
    end if
    ! Chi
    if (use_chi_z) then
       stat_nvar = stat_nvar+2
    end if
    if (use_chi_c) then
       stat_nvar = stat_nvar+2
    end if
    ! Combustion
    select case (trim(chemistry))
    case ('chemtable')
       stat_nvar = stat_nvar+8*nstat_scalar  ! 2*(species + T)*4 (mean, variance, 2x velocity correlation)
       !stat_nvar = stat_nvar+32
    end select
    ! Source terms
    if (use_chi_z.or.use_chi_c) then
       stat_nvar = stat_nvar+2
    end if
    ! SGS K
    if (give_ksgs) then
       stat_nvar = stat_nvar+1
    end if
    ! Soot
    if (use_soot) then
       stat_nvar = stat_nvar + 18
       if (.not.use_pah) stat_nvar = stat_nvar + 4
    end if
    ! Allocate
    allocate(stat_name(stat_nvar))
    
    ns = 0
    ! Density statistics
    if (trim(chemistry).ne.'none') then
       stat_name(ns+1) = 'RHO'
       stat_name(ns+2) = 'RHO^2'
       ns = ns+2
    end if
    ! Velocity statistics
    stat_name(ns+1) = 'U'
    stat_name(ns+2) = 'V'
    stat_name(ns+3) = 'W'
    stat_name(ns+4) = 'U^2'
    stat_name(ns+5) = 'V^2'
    stat_name(ns+6) = 'W^2'
    stat_name(ns+7) = 'UV'
    ns = ns+7
    if (trim(chemistry).ne.'none') then
       stat_name(ns+1) = 'rhoU'
       stat_name(ns+2) = 'rhoV'
       stat_name(ns+3) = 'rhoW'
       stat_name(ns+4) = 'rhoU^2'
       stat_name(ns+5) = 'rhoV^2'
       stat_name(ns+6) = 'rhoW^2'
       stat_name(ns+7) = 'rhoUV'
       ns = ns+7
    end if
    
    ! Scalars
    do isc=1,nscalar
       name = SC_name(isc)
       if (name(1:2).ne.'S_') then
          stat_name(ns+1) = 'SC-'  // trim(adjustl(sc_name(isc)))
          stat_name(ns+2) = 'SC^2-'// trim(adjustl(sc_name(isc)))
          stat_name(ns+3) = 'SC*U-'// trim(adjustl(sc_name(isc)))
          stat_name(ns+4) = 'SC*V-'// trim(adjustl(sc_name(isc)))
          ns = ns+4
          if (trim(chemistry).ne.'none') then
             stat_name(ns+1) = 'rhoS-'  // trim(adjustl(sc_name(isc)))
             stat_name(ns+2) = 'rhoS^2-'// trim(adjustl(sc_name(isc)))
             stat_name(ns+3) = 'rhoS*U-'// trim(adjustl(sc_name(isc)))
             stat_name(ns+4) = 'rhoS*V-'// trim(adjustl(sc_name(isc)))
             ns = ns+4
          end if
       end if
    end do

    ! Scalar Variance
    if (isc_ZVAR.eq.0 .and. isc_ZMIX2.eq.0 .and. use_ZVAR) then
       stat_name(ns+1) = 'ZVAR'
       stat_name(ns+2) = 'ZVAR^2'
       ns = ns+2
       if (trim(chemistry).ne.'none') then
          stat_name(ns+1) = 'rhoZVAR'
          stat_name(ns+2) = 'rhoZVAR^2'
          ns = ns+2
       end if
    end if

    ! Scalar Dissipation Rate
    if (isc_ZMIX.ne.0) then
       stat_name(ns+1) = 'CHI'
       stat_name(ns+2) = 'CHI^2'
       ns = ns+2
       if (trim(chemistry).ne.'none') then
          stat_name(ns+1) = 'rhoCHI'
          stat_name(ns+2) = 'rhoCHI^2'
          ns = ns+2
       end if
    end if

    ! Chi
    if (use_chi_z) then
       stat_name(ns+1) = 'CHI_ZZ'
       stat_name(ns+2) = 'CHI_HH'
       ns = ns+2
    end if
    if (use_chi_c) then
       stat_name(ns+1) = 'CHI_CC'
       stat_name(ns+2) = 'CHI_HH'
       ns = ns+2
    end if
    
    ! Combustion
    select case (trim(chemistry))
    case ('chemtable')
       do isc=1,nstat_scalar
          name = stat_scalar(isc)
          stat_name(ns+1) = 'SC-'  // trim(adjustl(name))
          stat_name(ns+2) = 'SC^2-'// trim(adjustl(name))
          stat_name(ns+3) = 'SC*U-'// trim(adjustl(name))
          stat_name(ns+4) = 'SC*V-'// trim(adjustl(name))
          stat_name(ns+5) = 'rhoS-'  // trim(adjustl(name))
          stat_name(ns+6) = 'rhoS^2-'// trim(adjustl(name))
          stat_name(ns+7) = 'rhoS*U-'// trim(adjustl(name))
          stat_name(ns+8) = 'rhoS*V-'// trim(adjustl(name))
          ns = ns+8
       end do
    end select

    ! Source Terms
    if (use_chi_z.or.use_chi_c) then
       stat_name(ns+1) = 'SRC_H'
       stat_name(ns+2) = 'SRC_C'
       ns = ns+2
    end if
    ! K_sgs
    if (give_ksgs) then
       stat_name(ns+1) = 'K_SGS'
       ns = ns+1
    end if
    
    ! Soot
    if (use_soot) then
       stat_name(ns+ 1) = 'fV'
       stat_name(ns+ 2) = 'fV^2'
       stat_name(ns+ 3) = 'N'
       stat_name(ns+ 4) = 'N^2'
       stat_name(ns+ 5) = 'dp'
       stat_name(ns+ 6) = 'dp^2'
       stat_name(ns+ 7) = 'np'
       stat_name(ns+ 8) = 'np^2'
       stat_name(ns+ 9) = 'intermit'
       stat_name(ns+10) = 'intermit^2'
       stat_name(ns+11) = 'dNdt_nucl'
       stat_name(ns+12) = 'dNdt_coag'
       stat_name(ns+13) = 'dNdt_ox'
       stat_name(ns+14) = 'dNdt_frag'
       stat_name(ns+15) = 'dfvdt_nucl'
       stat_name(ns+16) = 'dfvdt_cond'
       stat_name(ns+17) = 'dfvdt_sg'
       stat_name(ns+18) = 'dfvdt_ox'
       ns = ns+18
       if (.not.use_pah) then
          stat_name(ns+1) = 'Y_PAH'
          stat_name(ns+2) = 'rhoY_PAH'
          stat_name(ns+3) = 'Y_PAH^2'
          stat_name(ns+4) = 'rhoY_PAH^2'
          ns = ns+4
       end if
    end if

    return
  end subroutine stat_2d_init_names
  
end module stat_2d


! ================================== !
! Initialize the 2D statistic module !
! ================================== !
subroutine stat_2d_init
  use stat_2d
  use parallel
  use parser
  use combustion
  implicit none
  
  integer, dimension(2) :: gsizes,lsizes,start
  integer :: ierr
  
  ! Test if we can gather stats
  if (xper.eq.1) call die('stat_2d_init: 2D statistics in x impossible (x is periodic)')
  if (yper.eq.1) call die('stat_2d_init: 2D statistics in y impossible (y is periodic)')


  call parser_read('Track chi for Z and H',use_chi_z,.false.)
  call parser_read('Track chi for C and H',use_chi_c,.false.)
  call parser_read('Stats for k_sgs',give_ksgs,.false.)

  ! If using a chemtable, get the names of scalars to output
  if (trim(chemistry).eq.'chemtable') then
     call parser_getsize('Statistics scalar output',nstat_scalar)
     allocate(stat_scalar(nstat_scalar))
     call parser_read('Statistics scalar output',stat_scalar)
  end if
  
  ! Get the number of variables and names
  call stat_2d_init_names
  
  ! Allocate the storage space
  allocate(stat_xy(imin_:imax_,jmin_:jmax_,1:stat_nvar))
  allocate(buf(imin_:imax_,jmin_:jmax_,1:stat_nvar))
  stat_xy = 0.0_WP
  
  ! Generate the fileview
  gsizes(1) = nx
  lsizes(1) = nx_
  start(1)  = imin_-imin
  gsizes(2) = ny
  lsizes(2) = ny_
  start(2)  = jmin_-jmin
  call MPI_TYPE_CREATE_SUBARRAY(2,gsizes,lsizes,start,MPI_ORDER_FORTRAN,MPI_REAL_WP,fileview,ierr)
  call MPI_TYPE_COMMIT(fileview,ierr)
  
  ! Read the stat file
  call stat_2d_read
  
  return
end subroutine stat_2d_init


! ===================== !
! Sample the statistics !
! ===================== !
subroutine stat_2d_sample
  use stat_2d
  use data
  use combustion
  use sgsmodel
  use strainrate
  use filter
  use time_info
  use memory
  use soot
  use pah
  implicit none

  integer :: i,j,k,isc,ns,n, ns_loop1
  character(len=str_short) :: name

  ! Prepare the chemical variable if necessary
  if (trim(chemistry).eq.'chemtable') then
!!$     call chemtable_lookup('T'    ,tmp1)
!!$     call chemtable_lookup('O2'   ,tmp2)
!!$     call chemtable_lookup('H2'   ,tmp3)
!!$     call chemtable_lookup('H2O'  ,tmp4)
!!$     call chemtable_lookup('OH'   ,tmp5)
!!$     call chemtable_lookup('HO2'  ,tmp6)
!!$     call chemtable_lookup('H2O2' ,tmp7)
!!$     call chemtable_lookup('Y_F'  ,tmp2)
!!$     call chemtable_lookup('Y_O2' ,tmp3)
!!$     call chemtable_lookup('Y_CO' ,tmp4)
!!$     call chemtable_lookup('Y_CO2',tmp5)
!!$     call chemtable_lookup('Y_H2' ,tmp6)
!!$     call chemtable_lookup('Y_H2O',tmp7)
!!$     call chemtable_lookup('Y_OH' ,tmp8)
     if (use_soot .and. .not.use_pah) call chemtable_lookup('Y_PAH',tmp9)
  end if

  ! Calculate necessary scalar dissipation rates and source terms
  if (use_chi_c) then
     call gradient_squared(SC(:,:,:,isc_PROG),tmp10)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmp10(i,j,k) = 2.0_WP*tmp10(i,j,k)*DIFF(i,j,k,isc_PROG)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  if (use_chi_z.or.use_chi_c) then
     call gradient_squared(SC(:,:,:,isc_ENTH),tmp11)
     !$OMP PARALLEL DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmp11(i,j,k) = 2.0_WP*tmp11(i,j,k)*DIFF(i,j,k,isc_ENTH)
           end do
        end do
     end do
     !$OMP END PARALLEL DO

     call chemtable_lookup('SRC_RAD',tmp12)
     call chemtable_lookup('SRC_PROG',tmp13)
  end if

  if (give_ksgs) then
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              tmp14(i,j,k) = S(i,j,k) * Cs_visc(i,j,k) * delta_3D(i,j)
              tmp14(i,j,k) = 0.5_WP * tmp14(i,j,k)**2
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
     
  ! Gather the stats
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_

           ns = 0

           ! Density
           if (trim(chemistry).ne.'none') then
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*RHO(i,j,k)
              stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*RHO(i,j,k)**2
              ns = ns+2
           end if
           
           ! Velocity
           stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*U(i,j,k)
           stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*V(i,j,k)
           stat_xy(i,j,ns+3) = stat_xy(i,j,ns+3) + dt*W(i,j,k)
           stat_xy(i,j,ns+4) = stat_xy(i,j,ns+4) + dt*U(i,j,k)**2
           stat_xy(i,j,ns+5) = stat_xy(i,j,ns+5) + dt*V(i,j,k)**2
           stat_xy(i,j,ns+6) = stat_xy(i,j,ns+6) + dt*W(i,j,k)**2
           stat_xy(i,j,ns+7) = stat_xy(i,j,ns+7) + dt*U(i,j,k)*V(i,j,k)
           ns = ns+7
           
           if (trim(chemistry).ne.'none') then
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*rhoU(i,j,k)
              stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*rhoV(i,j,k)
              stat_xy(i,j,ns+3) = stat_xy(i,j,ns+3) + dt*rhoW(i,j,k)
              stat_xy(i,j,ns+4) = stat_xy(i,j,ns+4) + dt*rhoU(i,j,k)*U(i,j,k)
              stat_xy(i,j,ns+5) = stat_xy(i,j,ns+5) + dt*rhoV(i,j,k)*V(i,j,k)
              stat_xy(i,j,ns+6) = stat_xy(i,j,ns+6) + dt*rhoW(i,j,k)*W(i,j,k)
              stat_xy(i,j,ns+7) = stat_xy(i,j,ns+7) + dt*rhoU(i,j,k)*V(i,j,k)
              ns = ns+7
           end if
           
           ! Scalars
           do isc=1,nscalar
              name = SC_name(isc)
              if (name(1:2).ne.'S_') then
                 stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*SC(i,j,k,isc)
                 stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*SC(i,j,k,isc)**2
                 stat_xy(i,j,ns+3) = stat_xy(i,j,ns+3) + dt*SC(i,j,k,isc)*U(i,j,k) ! this might be incorrect (staggered grid) -- below, also
                 stat_xy(i,j,ns+4) = stat_xy(i,j,ns+4) + dt*SC(i,j,k,isc)*V(i,j,k)
                 ns = ns+4
                 if (trim(chemistry).ne.'none') then
                    stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*RHO(i,j,k)*SC(i,j,k,isc)
                    stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*RHO(i,j,k)*SC(i,j,k,isc)**2
                    stat_xy(i,j,ns+3) = stat_xy(i,j,ns+3) + dt*RHO(i,j,k)*SC(i,j,k,isc)*U(i,j,k)
                    stat_xy(i,j,ns+4) = stat_xy(i,j,ns+4) + dt*RHO(i,j,k)*SC(i,j,k,isc)*V(i,j,k)
                    ns = ns+4
                 end if
              end if
           end do

           ! Scalar Variance
           if (isc_ZVAR.eq.0 .and. isc_ZMIX2.eq.0 .and. use_ZVAR) then
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*ZVAR(i,j,k)
              stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*ZVAR(i,j,k)**2
              ns = ns+2
              if (trim(chemistry).ne.'none') then
                 stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*RHO(i,j,k)*ZVAR(i,j,k)
                 stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*RHO(i,j,k)*ZVAR(i,j,k)**2
                 ns = ns+2
              end if
           end if

           ! Scalar Dissipation Rate
           if (isc_ZMIX.ne.0) then
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*CHI(i,j,k)
              stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*CHI(i,j,k)**2
              ns = ns+2
              if (trim(chemistry).ne.'none') then
                 stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*RHO(i,j,k)*CHI(i,j,k)
                 stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*RHO(i,j,k)*CHI(i,j,k)**2
                 ns = ns+2
              end if
           end if

           ! Chi
           if (use_chi_z) then
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*CHI(i,j,k)
              stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*tmp11(i,j,k)
              ns = ns + 2
           end if
           if (use_chi_c) then
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*tmp10(i,j,k)
              stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*tmp11(i,j,k)
              ns = ns + 2
           end if

           ! need count for next loop
           ns_loop1 = ns

        end do
     end do
  end do
           
  ! Combustion
  select case (trim(chemistry))
  case ('chemtable')
     do isc=1,nstat_scalar
        call chemtable_lookup(trim(adjustl(stat_scalar(isc))) ,tmp1)

        ns = ns_loop1

        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_                 
                 stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*tmp1(i,j,k)
                 stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*tmp1(i,j,k)**2
                 stat_xy(i,j,ns+3) = stat_xy(i,j,ns+3) + dt*tmp1(i,j,k)*U(i,j,k)
                 stat_xy(i,j,ns+4) = stat_xy(i,j,ns+4) + dt*tmp1(i,j,k)*V(i,j,k)
                 stat_xy(i,j,ns+5) = stat_xy(i,j,ns+5) + dt*RHO(i,j,k)*tmp1(i,j,k)
                 stat_xy(i,j,ns+6) = stat_xy(i,j,ns+6) + dt*RHO(i,j,k)*tmp1(i,j,k)**2
                 stat_xy(i,j,ns+7) = stat_xy(i,j,ns+7) + dt*RHO(i,j,k)*tmp1(i,j,k)*U(i,j,k)
                 stat_xy(i,j,ns+8) = stat_xy(i,j,ns+8) + dt*RHO(i,j,k)*tmp1(i,j,k)*V(i,j,k)
              end do
           end do
        end do

        ns_loop1 = ns_loop1 + 8
     end do
  end select



  ! Src Terms
  if (use_chi_z .or. use_chi_c) then
     ns = ns_loop1
     
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*tmp12(i,j,k)
              stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*tmp13(i,j,k)
           end do
        end do
     end do

     ns_loop1 = ns_loop1 + 2
  end if

  ! KSGS
  if (give_ksgs) then
     ns = ns_loop1
     
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*tmp14(i,j,k)
           end do
        end do
     end do
     
     ns_loop1 = ns_loop1 + 1
  end if

           
  ! Soot
  if (use_soot) then
     
     ns = ns_loop1

     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              stat_xy(i,j,ns+ 1) = stat_xy(i,j,ns+ 1) + dt*volfrac(i,j,k)
              stat_xy(i,j,ns+ 2) = stat_xy(i,j,ns+ 2) + dt*volfrac(i,j,k)**2
              stat_xy(i,j,ns+ 3) = stat_xy(i,j,ns+ 3) + dt*numdens(i,j,k)
              stat_xy(i,j,ns+ 4) = stat_xy(i,j,ns+ 4) + dt*numdens(i,j,k)**2
              stat_xy(i,j,ns+ 5) = stat_xy(i,j,ns+ 5) + dt*partdiam(i,j,k)
              stat_xy(i,j,ns+ 6) = stat_xy(i,j,ns+ 6) + dt*partdiam(i,j,k)**2
              stat_xy(i,j,ns+ 7) = stat_xy(i,j,ns+ 7) + dt*partaggr(i,j,k)
              stat_xy(i,j,ns+ 8) = stat_xy(i,j,ns+ 8) + dt*partaggr(i,j,k)**2
              stat_xy(i,j,ns+ 9) = stat_xy(i,j,ns+ 9) + dt*intermit(i,j,k)
              stat_xy(i,j,ns+10) = stat_xy(i,j,ns+10) + dt*intermit(i,j,k)**2
              stat_xy(i,j,ns+11) = stat_xy(i,j,ns+11) + dt*Nsrc_nucl(i,j,k)
              stat_xy(i,j,ns+12) = stat_xy(i,j,ns+12) + dt*Nsrc_coag(i,j,k)
              stat_xy(i,j,ns+13) = stat_xy(i,j,ns+13) + dt*Nsrc_ox  (i,j,k)
              stat_xy(i,j,ns+14) = stat_xy(i,j,ns+14) + dt*Nsrc_frag(i,j,k)
              stat_xy(i,j,ns+15) = stat_xy(i,j,ns+15) + dt*FVsrc_nucl(i,j,k)
              stat_xy(i,j,ns+16) = stat_xy(i,j,ns+16) + dt*FVsrc_cond(i,j,k)
              stat_xy(i,j,ns+17) = stat_xy(i,j,ns+17) + dt*FVsrc_sg  (i,j,k)
              stat_xy(i,j,ns+18) = stat_xy(i,j,ns+18) + dt*FVsrc_ox  (i,j,k)
              ns = ns+18
              if (.not.use_pah) then
                 stat_xy(i,j,ns+1) = stat_xy(i,j,ns+1) + dt*tmp9(i,j,k)
                 stat_xy(i,j,ns+2) = stat_xy(i,j,ns+2) + dt*RHO(i,j,k)*tmp9(i,j,k)
                 stat_xy(i,j,ns+3) = stat_xy(i,j,ns+3) + dt*tmp9(i,j,k)**2
                 stat_xy(i,j,ns+4) = stat_xy(i,j,ns+4) + dt*RHO(i,j,k)*tmp9(i,j,k)**2
                 ns = ns+4
              end if
              
           end do
        end do
     end do
  end if
  
  Delta_t = Delta_t+dt
  
  return
end subroutine stat_2d_sample


! ================================= !
! Read the statistics from the disk !
! ================================= !
subroutine stat_2d_read
  use stat_2d
  use parallel
  implicit none
  
  real(WP) :: time
  integer, dimension(4) :: dims
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: disp
  character(len=str_medium)  :: name
  character(len=str_medium) :: filename
  integer :: var,ierr,ifile,data_size
  logical :: file_is_there
  
  ! -- GAS PHASE ---------------------------------------------------------------------
  filename = "stat/stat-2D"
  if (use_mpiiofs) filename = trim(mpiiofs) // ":" // trim(filename)
  inquire(file=filename, exist=file_is_there)

  if (file_is_there) then
     
     ! Open the file to write
     call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,ifile,ierr)
     
     ! Read dimensions from header
     call MPI_FILE_READ_ALL(ifile,dims,4,MPI_INTEGER,status,ierr)
     if ((dims(1).ne.nx) .or. (dims(2).ne.ny) .or. (dims(3).ne.1)) then
        print*, 'expected = ',nx,ny,1
        print*, 'stat = ',dims(1),dims(2),dims(3)
        call die('stat_2d_read: The size of the stat file is incorrect')
     end if
     if (dims(4).ne.stat_nvar) call die('stat_2d_read: Wrong number of variables in stat file')
     
     ! Read some headers
     call MPI_FILE_READ_ALL(ifile,Delta_t,1,MPI_REAL_WP,status,ierr)
     call MPI_FILE_READ_ALL(ifile,time,1,MPI_REAL_WP,status,ierr)
     
     ! Read variable names
     do var=1,stat_nvar
        call MPI_FILE_READ_ALL(ifile,name,str_medium,MPI_CHARACTER,status,ierr)
        if (name.ne.stat_name(var)) then
           call die('stat_2d_read: Variables names in stat are incorrect')
        end if
     end do
          
     ! Read each variables
     data_size = nx_*ny_
     do var=1,stat_nvar
        disp = 4*4 + str_medium*stat_nvar + 2*WP + real(var-1,WP)*nx*ny*WP
        call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,fileview,"native",MPI_INFO_NULL,ierr)
        call MPI_FILE_READ_ALL(ifile,buf(:,:,var),data_size,MPI_REAL_WP,status,ierr)
     end do
     
     ! Close the file
     call MPI_FILE_CLOSE(ifile,ierr)
     
     ! Recompute the stats
     stat_xy = buf*Delta_t*real(nz_)
     
  else
     
     ! Start from scratch
     Delta_t = 0.0_WP
     
  end if
  
  return
end subroutine stat_2d_read


! ================================ !
! Write the statistics to the disk !
! ================================ !
subroutine stat_2d_write
  use stat_2d
  use parallel
  use time_info
  implicit none
  
  integer, dimension(4) :: dims
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: disp
  character(len=str_medium) :: filename
  integer :: var,ierr,ifile,data_size,i,j
  logical :: file_is_there
  
  ! -- GAS PHASE ---------------------------------------------------------------------
  
  ! Gather the data
  call parallel_sum_dir(stat_xy,buf,'z')
  buf = buf / (Delta_t*real(nz))
  
  ! Open the file to write
  filename = "stat/stat-2D"
  if (use_mpiiofs) filename = trim(mpiiofs) // ":" // trim(filename)
  inquire(file=filename, exist=file_is_there)
  if (file_is_there.and.irank.eq.iroot) call MPI_FILE_DELETE(filename,MPI_INFO_NULL,ierr)
  call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,ifile,ierr)
  
  ! Write the headers
  if (irank.eq.iroot) then
     ! Write dimensions
     dims(1) = nx
     dims(2) = ny
     dims(3) = 1
     dims(4) = stat_nvar
     call MPI_FILE_WRITE(ifile,dims,4,MPI_INTEGER,status,ierr)
     call MPI_FILE_WRITE(ifile,Delta_t,1,MPI_REAL_WP,status,ierr)
     call MPI_FILE_WRITE(ifile,time,1,MPI_REAL_WP,status,ierr)
     ! Write variable names
     do var=1,stat_nvar
        call MPI_FILE_WRITE(ifile,stat_name(var),str_medium,MPI_CHARACTER,status,ierr)
     end do
  end if
  
  ! Write each variables
  data_size = nx_*ny_
  if (kproc.ne.1) data_size = 0
  do var=1,stat_nvar
     disp = 4*4 + str_medium*stat_nvar + 2*WP + real(var-1,WP)*nx*ny*WP
     call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,fileview,"native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(ifile,buf(:,:,var),data_size,MPI_REAL_WP,status,ierr)
  end do
  
  ! Close the file
  call MPI_FILE_CLOSE(ifile,ierr)
  
  return
end subroutine stat_2d_write
