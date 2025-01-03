module scalar
  use precision
  use partition
  use geometry
  use time_info
  use velocity
  implicit none
  
  ! Solution at old time step : n
  real(WP), dimension(:,:,:,:), pointer :: SCold

  ! Solution at mid-point of split system : n+1
  real(WP), dimension(:,:,:,:), pointer :: SC_
  real(WP), dimension(:,:,:,:), pointer :: SC_s
  real(WP), dimension(:,:,:),   pointer :: RHO_
  real(WP), dimension(:,:,:),   pointer :: RHO_s
  
  ! Source term between time n and n+1
  real(WP), dimension(:,:,:,:), pointer :: srcSCmid
  real(WP), dimension(:,:,:,:), pointer :: srcSCfull
  
  ! Residual between time n and n+1
  real(WP), dimension(:,:,:,:), pointer :: ResSC
  
  ! Values to monitor
  real(WP), dimension(:), pointer :: max_resSC
  real(WP), dimension(:), pointer :: ext_SC

  ! Formulation: Strang splitting or un-split
  logical :: strang_splitting

  ! N2 clipping for mass conservation
  logical :: clip_N2
  real(WP) :: clip_N2_val
  
end module scalar


! ==================================================== !
! Initialize the scalar module                         !
!                                                      !
! -> allocate the arrays                               !
! -> update ghost cells                                !
! -> apply boundary conditions                         !
! -> run specific init                                 !
!                                                      !
! Before: SC correct only inside the domain            !
!         -> imin_:imax_,jmin_:jmax_,kmin_:kmax_       !
! After : SC correct everywhere                        !
!         -> imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_ !
! ==================================================== !
subroutine scalar_init
  use scalar
  use parser
  use data
  use implicit
  implicit none
  
  integer :: isc,i,j,k, mpi_id,mpi_err
  
  ! If no scalar => exit
  if (nscalar.eq.0) return
  
  ! Create & Start the timer
  call timing_create('scalar')
  call timing_start ('scalar')
  
  ! Allocate arrays for old solution
  allocate(SCold    (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(SC_      (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(SC_s     (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(srcSCmid (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(srcSCfull(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(ResSC    (imin_ :imax_ ,jmin_ :jmax_ ,kmin_ :kmax_ ,nscalar))
  allocate(RHO_     (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(RHO_s    (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  
  ! The scalars fields read from the file were without ghost cells
  ! Update the ghost cells for periodicity and domain decomposition
  do isc=1,nscalar
     call boundary_update_border(SC(:,:,:,isc),'+','ym')
  end do

  ! Apply boundary conditions
  call boundary_scalar_dirichlet
  call boundary_scalar_outflow
  call boundary_scalar_neumann

  ! Allocate arrays for the total advecting momentum (thermophoresis, etc.)
  allocate(Ut(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(Vt(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(Wt(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))

  allocate(rhoUt(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(rhoVt(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))
  allocate(rhoWt(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nscalar))

  ! Determine if using the approximately factorized Jacobian
  call parser_read('AF Jacobian',AF_jac,.false.)
  call parser_read('Diagonal chem Jacobian',AF_diag_jac,.false.)
  if (AF_diag_jac) AF_jac = .true.
  
  ! Initialize the given scheme
  select case (trim(scalar_scheme))
  case ('quick')
     call scalar_quick_init
  case ('weno3')
     call scalar_weno3_init
  case ('weno5')
     call scalar_weno5_init
  case ('houc')
     call scalar_houc_init
  case default
     call die('Unknown scalar scheme specified')
  end select

  ! Determine if Strang splitting or not
  call parser_read('Strang splitting',strang_splitting,.false.)

  if (strang_splitting.and.AF_jac) then
     call die('scalar_init: Cannot request both Strang splitting and approximately factorized Jacobian.')
  end if
  
  ! Create new file to monitor at each iterations
  call monitor_create_file_step('scalar',2*nscalar)
  allocate(ext_SC(2*nscalar))
  do isc=1,nscalar
     call monitor_set_header(2*isc-1,'min_'//SC_name(isc),'r')
     call monitor_set_header(2*isc+0,'max_'//SC_name(isc),'r')
  end do
  
  ! Create new file to monitor at each subiterations
  call monitor_create_file_iter('convergence_scalar',nscalar)
  allocate(max_resSC(nscalar))
  do isc=1,nscalar
     call monitor_set_header(isc,'res_'//SC_name(isc),'r')
  end do

  ! Determine if doing mms
  call parser_read('MMS',mms,.false.)

  ! Determine if clipping N2
  call parser_is_defined('Max N2 value',clip_N2)
  if (clip_N2) call parser_read('Max N2 value',clip_N2_val)
  
  ! Stop the timer
  call timing_stop('scalar')
  
  return
end subroutine scalar_init


! ================================================== !
! PRE-TIMESTEP Routine                               !
!                                                    !
! -> Set up the iterative process                    !
! -> Compute the source term for the scalar equation !
! ================================================== !
subroutine scalar_prestep
  use velocity
  use scalar
  use data
  implicit none
  integer :: i,j,k,isc
  
  ! If no scalar => exit
  if (nscalar.eq.0) return
  
  ! Start the timer
  call timing_start('scalar')

  !$OMP PARALLEL

  ! Save the old scalar
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SCold(i,j,k,isc) = SC(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do
     
  ! Zero the source term
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              ! Mid source if part of C-N: CHEMTABLE
              ! IS THIS ACCURATE ENOUGH FOR SOOT?  NOT BASED ON VIDA RESULTS!
              srcSCmid (i,j,k,isc) = 0.0_WP
              ! Full source is split and computed prestep: FINITE RATE
              srcSCfull(i,j,k,isc) = 0.0_WP
           end do
        end do
     end do
     !$OMP END DO
  end do

  !$OMP END PARALLEL

  ! Source term from SGS -- LOOK AT THIS LATER
  call sgsmodel_src_sc(srcSCfull)
  
  ! Stop the timer
  call timing_stop('scalar')
  
  return
end subroutine scalar_prestep


! ========================================================== !
! ADVANCE the solution                                       !
!   -> second order in time                                  !
!   -> variable accuracy in space                            !
!   -> explicit prediction                                   !
!                                                            !
! Z(n+3/2,k+1) = Z(n+1/2) + dt*F(0.5*(Z(n+3/2,k)+Z(n+1/2))   !
!                   + 0.5*dt*dF/dZ*(Z(n+3/2,k+1)-Z(n+3/2,k)) !
! n : time step                                              !
! k : inner loop iteration                                   !
! Velocity field used : best approximation for U(n+1)        !
! ========================================================== !
subroutine scalar_step
  use scalar
  use data
  implicit none
  integer :: i,j,k,isc,n
  
  ! If no scalar => exit
  if (nscalar.eq.0) return
  
  ! Start the timer
  call timing_start('scalar')

  !$OMP PARALLEL

  ! Clip N2
  if (clip_N2) then
     do isc=1,nscalar
        if (trim(SC_name(isc)).eq.'N2') then
           !$OMP DO
           do k=kmino_,kmaxo_
              do j=jmino_,jmaxo_
                 do i=imino_,imaxo_
                    SC(i,j,k,isc) = min(SC(i,j,k,isc), clip_N2_val)
                 end do
              end do
           end do
           !$OMP END DO
        end if
     end do
  end if

  ! JFM 7/9/14
  ! Save the old 'n+3/2' scalar
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SC_s(i,j,k,isc) = SC(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do

  !! NEED TO EVALUATE DIFFUSIVITY AT N+1
  !! ALSO NEED VISCOSITY AT N+1/2
  !! -- NEED TO PUT SGSMODEL IN SUBITERATIONS AND SYNC THOSE IN TIME
  ! Molecular diffusivity
  !!$  ! Compute diffusivity for transport (mid-point)
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SC(i,j,k,isc) = 0.5_WP * (SC(i,j,k,isc) + SCold(i,j,k,isc))
           end do
        end do
     end do
     !$OMP END DO
  end do
  !$OMP END PARALLEL

  call combustion_diffusivity
  call pollutants_diffusivity

  ! Subfilter diffusivity
  ! moved here because also needs to be avaluated with scalars at N+1 (midpoint of SCALAR step)
  call sgsmodel_eddyDIFF
  
  !$OMP PARALLEL
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SC(i,j,k,isc) = 2.0_WP * SC(i,j,k,isc) - SCold(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do
  !$OMP END PARALLEL
  
  !$OMP PARALLEL
  ! Compute the combustion source terms
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              srcSCmid(i,j,k,isc) = 0.0_WP
           end do
        end do
     end do
     !$OMP END DO
  end do

  ! Compute the total advecting momemuntum for each scalar
  ! - Thermophoresis, inertial effects, etc.
  ! - rhoUt, rhoVt, rhoWt used in scalar solver and soot. Be careful when using
  !   Strang-split solver with soot, since these arrays are reused later.
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              rhoUt(i,j,k,isc) = rhoU(i,j,k)
              rhoVt(i,j,k,isc) = rhoV(i,j,k)
              rhoWt(i,j,k,isc) = rhoW(i,j,k)
           end do
        end do
     end do
     !$OMP END DO
  end do
  !$OMP END PARALLEL

  !call combustion_total_momentum
  !call pollutants_total_momentum
  !call spray_total_momentum

  ! Compute the total advecting velocity
  call scalar_rho_divide

  !$OMP PARALLEL
  ! Store (*)old values to (*)_ for input to scalar solver
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SC_(i,j,k,isc) = SCold(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do

  !$OMP DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           RHO_(i,j,k) = RHOold(i,j,k)
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Switch for Strang splitting
  if (strang_splitting) then

      ! 1. Set dt=dt/2 and compute transport in first half step
      dt_ = dt*0.5_WP

      !$OMP PARALLEL

      !! Save the current subiteration's n+3/2 density to RHO_s
      !$OMP DO
      do k=kmino_,kmaxo_
         do j=jmino_,jmaxo_
            do i=imino_,imaxo_
               RHO_s(i,j,k) = RHO(i,j,k)
            end do
         end do
      end do
      !$OMP END DO

      ! JFM 6/13/14
      ! Interpolate velocity and density for first half of scalar transport
      ! -- Compute mid point (n+3/4)
      ! MEM -- This needs to be changed; rhoUt is the total momentum for soot/sprays
      ! JFM -- Maybe set it back to the correct value after we're done w splitting
      do isc=1,nscalar
         !$OMP DO
         do k=kmino_,kmaxo_
            do j=jmino_,jmaxo_
               do i=imino_,imaxo_
                  rhoUt(i,j,k,isc) = 0.75_WP*rhoU(i,j,k) + 0.25_WP*rhoUold(i,j,k)
                  rhoVt(i,j,k,isc) = 0.75_WP*rhoV(i,j,k) + 0.25_WP*rhoVold(i,j,k)
                  rhoWt(i,j,k,isc) = 0.75_WP*rhoW(i,j,k) + 0.25_WP*rhoWold(i,j,k)
               end do
            end do
         end do
         !$OMP END DO
      end do


      ! -- Interpolate density at n+1
      !$OMP DO
      do k=kmino_,kmaxo_
         do j=jmino_,jmaxo_
            do i=imino_,imaxo_
               RHO(i,j,k) = 0.5_WP*(RHO(i,j,k)+RHOold(i,j,k))
            end do
         end do
      end do
      !$OMP END DO

      ! JFM 11/12/14
      ! Compute mid point for first half of scalar transport (n+3/4).
      ! Store it in the 'n+3/2' scalar
      do isc=1,nscalar
         !$OMP DO
         do k=kmino_,kmaxo_
            do j=jmino_,jmaxo_
               do i=imino_,imaxo_
                  SC(i,j,k,isc) = 0.25_WP*SC_s(i,j,k,isc) + 0.75_WP*SCold(i,j,k,isc)
               end do
            end do
         end do
         !$OMP END DO
      end do
      !$OMP END PARALLEL

      ! JFM 8/22/14
      ! Recompute diffusivity at n+3/4 for rhs evaluation (finite rate chemistry only)
      if (trim(chemistry).eq.'finite chem') call combustion_diffusivity

      ! Compute the source terms
      call combustion_source_scalar(srcSCmid)
      call pollutants_source_scalar(SC,srcSCmid,srcP)
      call spray_source_scalar(SC,srcSCmid,srcP)

      select case (trim(scalar_scheme))
      case ('quick')
         call scalar_quick_residual
         call scalar_quick_inverse
      case ('weno3')
         call scalar_weno3_residual
         call scalar_weno3_inverse
      case ('weno5')
         call scalar_weno5_residual
         call scalar_weno5_inverse
      case ('houc')
         call scalar_houc_residual
         call scalar_houc_inverse
      end select

      !$OMP PARALLEL

      ! Update the scalars after first half step, save to SC_ ("SC prime").
      ! SC^n+1_k+1 = SC^n+1_k + (delta SC)
      do isc=1,nscalar
         !$OMP DO
         do k=kmin_,kmax_
            do j=jmin_,jmax_
               do i=imin_,imax_
                  SC_(i,j,k,isc) = 2.0_WP*SC(i,j,k,isc)-SC_(i,j,k,isc) + ResSC(i,j,k,isc)
               end do
            end do
         end do
         !$OMP END DO
      end do
      !$OMP END PARALLEL

      ! Update the overlapped cells
      do isc=1,nscalar
         call boundary_update_border(SC_(:,:,:,isc),'+','ym')
      end do

      ! Update the physical boundaries on scalar with previous subiteration
      ! -- THIS NEEDS TO BE BETTER...
      do k=kmino_,kmaxo_
         if (k.lt.kmin .or. k.gt.kmax) then
            do j=jmino_,jmaxo_
               if (j.lt.jmin .or. j.gt.jmax) then
                  do i=imino_,imaxo_
                     if (i.lt.imin .or. i.gt.imax) then
                        SC_(i,j,k,:) = SC_s(i,j,k,:)
                     end if
                  end do
               end if
            end do
         end if
      end do


      ! 2. Compute the chemistry (CVODE in finitechem_source). Yields SC_ ("SC hat").
      call combustion_source_scalar_full(RHO,SC_)


      ! Update the physical boundaries on scalar with previous subiteration
      ! -- THIS NEEDS TO BE BETTER...
      do k=kmino_,kmaxo_
         if (k.lt.kmin .or. k.gt.kmax) then
            do j=jmino_,jmaxo_
               if (j.lt.jmin .or. j.gt.jmax) then
                  do i=imino_,imaxo_
                     if (i.lt.imin .or. i.gt.imax) then
                        SC_(i,j,k,:) = SC_s(i,j,k,:)
                     end if
                  end do
               end if
            end do
         end if
      end do

      ! Update the overlapped cells
      do isc=1,nscalar
         call boundary_update_border(SC_(:,:,:,isc),'+','ym')
      end do
      call boundary_update_border(RHO(:,:,:),'+','ym')


      ! 3. Compute transport in second half step
      ! JFM 6/13/14
      ! Interpolate density and extrapolate velocity for second half of scalar transport
      ! -- Extrapolate velocity to (n+5/4); store it to rhoUt
      !$OMP PARALLEL
      do isc=1,nscalar
         !$OMP DO
         do k=kmino_,kmaxo_
            do j=jmino_,jmaxo_
               do i=imino_,imaxo_
                  rhoUt(i,j,k,isc) = 1.25_WP*rhoU(i,j,k) - 0.25_WP*rhoUold(i,j,k)
                  rhoVt(i,j,k,isc) = 1.25_WP*rhoV(i,j,k) - 0.25_WP*rhoVold(i,j,k)
                  rhoWt(i,j,k,isc) = 1.25_WP*rhoW(i,j,k) - 0.25_WP*rhoWold(i,j,k)
               end do
            end do
         end do
         !$OMP END DO
      end do

      ! -- Recall saved density at n+3/2
      !$OMP DO
      do k=kmino_,kmaxo_
         do j=jmino_,jmaxo_
            do i=imino_,imaxo_
               ! RHO_: density at n+1
               RHO_(i,j,k) = RHO(i,j,k)
            end do
         end do
      end do
      !$OMP END DO
      !$OMP DO
      do k=kmino_,kmaxo_
         do j=jmino_,jmaxo_
            do i=imino_,imaxo_
               ! RHO: density at n+3/2
               RHO(i,j,k)  = RHO_s(i,j,k)
            end do
         end do
      end do
      !$OMP END DO

      ! JFM 3/11/15
      ! Compute mid point for second half of scalar transport (n+5/4)
      ! using saved values at 'n+3/2' (SC_s) and 'n+1/2' (SCold).
      ! Inconsistent if using the most current value of SC_{k+1}^{n+1} (SC_):
      ! Cannot use SC with chem to update transport in the operator-split framework.
      ! Store it in the 'n+3/2' scalar.
      do isc=1,nscalar
         !$OMP DO
         do k=kmino_,kmaxo_
            do j=jmino_,jmaxo_
               do i=imino_,imaxo_
                  !! This is key:
                  SC(i,j,k,isc) = 0.75_WP*SC_s(i,j,k,isc) + 0.25_WP*SCold(i,j,k,isc)
               end do
            end do
         end do
         !$OMP END DO
      end do
      !$OMP END PARALLEL

      ! JFM 8/22/14
      ! Recompute diffusivity at n+5/4 for rhs evaluation (finite rate chemistry only)
      if (trim(chemistry).eq.'finite chem') call combustion_diffusivity

      !$OMP PARALLEL

      ! JFM 8/26/14
      ! Recompute source terms at n+5/4 for rhs evaluation
      do isc=1,nscalar
         !$OMP DO
         do k=kmino_,kmaxo_
            do j=jmino_,jmaxo_
               do i=imino_,imaxo_
                  srcSCmid(i,j,k,isc) = 0.0_WP
               end do
            end do
         end do
         !$OMP END DO
      end do
      !$OMP END PARALLEL

      ! -think about compressible vs. incompressible!
      !    

  else
      ! Non split solver: use normal dt
      dt_ = dt

      !$OMP PARALLEL

      ! Compute mid point (n+1)
      ! Store it in the 'n+3/2' scalar
      do isc=1,nscalar
         !$OMP DO
         do k=kmino_,kmaxo_
            do j=jmino_,jmaxo_
               do i=imino_,imaxo_
                  SC(i,j,k,isc) = 0.5_WP*(SC(i,j,k,isc)+SCold(i,j,k,isc))
               end do
            end do
         end do
         !$OMP END DO
      end do
      !$OMP END PARALLEL

      ! JFM 9/2/14
      ! Recompute diffusivity at n+1 for rhs evaluation (finite rate chemistry only)
      if (trim(chemistry).eq.'finite chem') call combustion_diffusivity

  end if

  ! Compute the source terms using SC from previous subiteration
  call combustion_source_scalar(srcSCmid)
  call pollutants_source_scalar(SC,srcSCmid,srcP)
  call spray_source_scalar(SC,srcSCmid,srcP)

  ! Compute transport: non-split or second half step of split scheme
  select case (trim(scalar_scheme))
  case ('quick')
     call scalar_quick_residual
     call scalar_quick_inverse
  case ('weno3')
     call scalar_weno3_residual
     call scalar_weno3_inverse
  case ('weno5')
     call scalar_weno5_residual
     call scalar_weno5_inverse
  case ('houc')
     call scalar_houc_residual
     call scalar_houc_inverse
  end select

  !$OMP PARALLEL
  
  do isc=1,nscalar
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              SC(i,j,k,isc) = 2.0_WP*SC(i,j,k,isc)-SC_(i,j,k,isc) + ResSC(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do

  !$OMP END PARALLEL

  !MEM
  if (.not.(strang_splitting.or.AF_jac)) then
      ! Compute reaction source term
     call combustion_source_scalar_full(RHO,SC)
  end if

  ! Update the physical boundaries
  call boundary_scalar_dirichlet
  call boundary_scalar_neumann
  call boundary_scalar_outflow

  ! Update the overlapped cells
  do isc=1,nscalar
     call boundary_update_border(SC(:,:,:,isc),'+','ym')
  end do

  ! Compute max of residuals
  do isc=1,nscalar
     call parallel_max(maxval(abs(resSC(:,:,:,isc))),max_resSC(isc))
  end do
  
  ! Transfer values to monitor
  call monitor_select_file('convergence_scalar')
  call monitor_set_array_values(max_resSC)
  
  ! Stop the timer
  call timing_stop('scalar')
  
  return
end subroutine scalar_step


! =================================================== !
! Compute the total velocity components from momentum !
! =================================================== !
subroutine scalar_rho_divide
  use scalar
  use velocity
  use data
  use metric_generic
  implicit none

  integer :: isc,i,j,k

  !$OMP PARALLEL
  do isc=1,nscalar

     ! rhoU
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              Ut(i,j,k,isc) = rhoUt(i,j,k,isc) / sum(interp_sc_x(i,j,:)*RHOmid(i-st2:i+st1,j,k))
           end do
        end do
     end do
     !$OMP END DO
     
     ! rhoV
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              Vt(i,j,k,isc) = rhoVt(i,j,k,isc) / sum(interp_sc_y(i,j,:)*RHOmid(i,j-st2:j+st1,k))
           end do
        end do
     end do
     !$OMP END DO
  
     ! rhoW
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              Wt(i,j,k,isc) = rhoWt(i,j,k,isc) / sum(interp_sc_z(i,j,:)*RHOmid(i,j,k-st2:k+st1))
           end do
        end do
     end do
     !$OMP END DO

  end do
  !$OMP END PARALLEL

  ! CPU borders and periodicity
  do isc=1,nscalar
     call boundary_update_border(Ut(:,:,:,isc),'+','ym')
     call boundary_update_border(Vt(:,:,:,isc),'-','y')
     call boundary_update_border(Wt(:,:,:,isc),'-','ym')
  end do

  return
end subroutine scalar_rho_divide


! =================== !
! Monitor the scalars !
! =================== !
subroutine scalar_monitor
  use scalar
  use data
  use masks
  implicit none
  
  integer :: i, j, k, isc
  real(WP) :: min_SC_, max_SC_
  
  if (nscalar.eq.0) return
  
  ! Start the timer
  call timing_start('scalar')

  ! Get min/max
  do isc=1,nscalar
     min_SC_ = +huge(1.0_WP)
     max_SC_ = -huge(1.0_WP)

     !$OMP PARALLEL DO REDUCTION(min:min_SC_) REDUCTION(max:max_SC_)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              if (mask(i,j).ne.1) then
                 if (SC(i,j,k,isc).lt.min_SC_) min_SC_ = SC(i,j,k,isc)
                 if (SC(i,j,k,isc).gt.max_SC_) max_SC_ = SC(i,j,k,isc)
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call parallel_max(-min_SC_,ext_SC(2*isc-1))
     ext_SC(2*isc-1) = -ext_SC(2*isc-1)
     call parallel_max( max_SC_,ext_SC(2*isc+0))
  end do
  
  ! Transfer values to monitor
  call monitor_select_file('scalar')
  call monitor_set_array_values(ext_SC)
  
  ! Stop the timer
  call timing_stop('scalar')
  
  return
end subroutine scalar_monitor
