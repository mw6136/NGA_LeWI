module planeStats_budgets
  use planeStats_average

  real(WP), dimension(:,:,:), pointer :: tmpv_save_TKE_phys
  real(WP), dimension(:,:,:), pointer :: tmpv_save_TKE_cond
  real(WP), dimension(:,:,:), pointer :: tmpv_save_SCV_phys
  real(WP), dimension(:,:,:), pointer :: tmpv_save_SCV_cond
  integer :: isc_sc

  ! Output quantities
  ! TKE budgets - physical space
  real(WP), dimension(:,:), pointer :: E_conv_phys, E_div1_phys, E_div2_phys, E_div3_phys
  real(WP), dimension(:,:), pointer :: E_pres_phys, E_prod_phys, E_diss_phys, E_unst_phys
  real(WP), dimension(:,:), pointer :: E_TKE_phys, eta_phys
  real(WP), dimension(:,:), pointer :: E_pwork_phys, E_ptrans_phys, E_pdil_phys
  ! TKE budgets - conditional space
  real(WP), dimension(:,:), pointer :: E_conv_cond, E_div1_cond, E_div2_cond, E_div3_cond
  real(WP), dimension(:,:), pointer :: E_pres_cond, E_prod_cond, E_diss_cond, E_unst_cond
  real(WP), dimension(:,:), pointer :: E_TKE_cond, eta_cond

  ! Reynolds stress budget
  real(WP), dimension(:,:,:,:), pointer :: R_stress_budg
  real(WP), dimension(:,:,:,:), pointer :: vel_triple_corr

  ! Scalar variance budgets - physical space
  real(WP), dimension(:,:), pointer :: SCV_conv_phys, SCV_div1_phys, SCV_div2_phys
  real(WP), dimension(:,:), pointer :: SCV_prod_phys, SCV_diss_phys, SCV_diss2_phys, SCV_unst_phys
  real(WP), dimension(:,:), pointer :: SCV_src_phys,  SCV_phys
  ! Scalar variance budgets - conditional space
  real(WP), dimension(:,:), pointer :: SCV_conv_cond, SCV_div1_cond, SCV_div2_cond
  real(WP), dimension(:,:), pointer :: SCV_prod_cond, SCV_diss_cond, SCV_diss2_cond, SCV_unst_cond
  real(WP), dimension(:,:), pointer :: SCV_src_cond,  SCV_cond

end module planeStats_budgets

! ========================================================== !
! Compute the budgets in progress variable                   !
!    or mixture fraction space                               !
! ========================================================== !
subroutine planeStats_budgets_init
  use planeStats_budgets
  use parser
  implicit none

  ! Variables in module planeStats
  allocate(dTAUdxp  (1:nz,1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(dTAUdxp_c(1:nz,1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(dPdxp    (1:nz,1:ny,pnmin_:pnmax_,1:3))
  allocate(dPdxp_c  (1:nz,1:ny,pnmin_:pnmax_,1:3))

  ! Allocate computed variables
  ! Physical space
  allocate(Pf    (1:nz,1:ny,pnmin_:pnmax_))
  allocate(Ui    (1:nz,1:ny,pnmin_:pnmax_,1:3))
  allocate(dUdxi (1:nz,1:ny,pnmin_:pnmax_,1:3,1:3))

  ! Conditional space
  allocate(U_c      (1:nz,1:ny,pnmin_:pnmax_,1:3))
  allocate(dUdx_c   (1:nz,1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(SCi      (1:nz,1:ny,pnmin_:pnmax_,1:nscalar))
  allocate(SC_c     (1:nz,1:ny,pnmin_:pnmax_,1:nscalar))
  allocate(dSCdxi   (1:nz,1:ny,pnmin_:pnmax_,1:nscalar,1:3))
  allocate(dSCdx_c  (1:nz,1:ny,pnmin_:pnmax_,1:nscalar,1:3))

  ! TKE budgets - physical space
  allocate(E_unst_phys(1:ny,pnmin_:pnmax_))
  allocate(E_conv_phys(1:ny,pnmin_:pnmax_))
  allocate(E_div1_phys(1:ny,pnmin_:pnmax_))
  allocate(E_div2_phys(1:ny,pnmin_:pnmax_))
  allocate(E_div3_phys(1:ny,pnmin_:pnmax_))
  allocate(E_pres_phys(1:ny,pnmin_:pnmax_))
  allocate(E_prod_phys(1:ny,pnmin_:pnmax_))
  allocate(E_diss_phys(1:ny,pnmin_:pnmax_))
  allocate(E_TKE_phys (1:ny,pnmin_:pnmax_))
  allocate(eta_phys   (1:ny,pnmin_:pnmax_))
  allocate(E_pwork_phys (1:ny,pnmin_:pnmax_))
  allocate(E_ptrans_phys(1:ny,pnmin_:pnmax_))
  allocate(E_pdil_phys  (1:ny,pnmin_:pnmax_))
  E_unst_phys = 0.0_WP
  E_conv_phys = 0.0_WP
  E_div1_phys = 0.0_WP
  E_div2_phys = 0.0_WP
  E_div3_phys = 0.0_WP
  E_pres_phys = 0.0_WP
  E_prod_phys = 0.0_WP
  E_diss_phys = 0.0_WP
  E_TKE_phys  = 0.0_WP
  E_pwork_phys  = 0.0_WP
  E_ptrans_phys = 0.0_WP
  E_pdil_phys   = 0.0_WP
  eta_phys    = 0.0_WP
  ! TKE budgets - conditional space
  allocate(E_unst_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(E_conv_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(E_div1_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(E_div2_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(E_div3_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(E_pres_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(E_prod_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(E_diss_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(E_TKE_cond (1:nbins_cond,pnmin_:pnmax_))
  allocate(eta_cond   (1:nbins_cond,pnmin_:pnmax_))
  E_unst_cond = 0.0_WP
  E_conv_cond = 0.0_WP
  E_div1_cond = 0.0_WP
  E_div2_cond = 0.0_WP
  E_div3_cond = 0.0_WP
  E_pres_cond = 0.0_WP
  E_prod_cond = 0.0_WP
  E_diss_cond = 0.0_WP
  E_TKE_cond  = 0.0_WP
  eta_cond    = 0.0_WP

  ! Reynolds stress budget
  ! Six independent components
  ! Ten terms in budget
  allocate(R_stress_budg(1:ny,pnmin_:pnmax_,1:6,1:10))
  allocate(vel_triple_corr(1:ny,pnmin_:pnmax_,1:3,1:3))

  ! Scalar variance budgets - physical space
  allocate(SCV_unst_phys(1:ny,pnmin_:pnmax_))
  allocate(SCV_conv_phys(1:ny,pnmin_:pnmax_))
  allocate(SCV_div1_phys(1:ny,pnmin_:pnmax_))
  allocate(SCV_div2_phys(1:ny,pnmin_:pnmax_))
  allocate(SCV_prod_phys(1:ny,pnmin_:pnmax_))
  allocate(SCV_diss_phys(1:ny,pnmin_:pnmax_))
  allocate(SCV_diss2_phys(1:ny,pnmin_:pnmax_))
  allocate(SCV_src_phys (1:ny,pnmin_:pnmax_))
  allocate(SCV_phys     (1:ny,pnmin_:pnmax_))
  SCV_unst_phys = 0.0_WP
  SCV_conv_phys = 0.0_WP
  SCV_div1_phys = 0.0_WP
  SCV_div2_phys = 0.0_WP
  SCV_prod_phys = 0.0_WP
  SCV_diss_phys = 0.0_WP
  SCV_diss2_phys = 0.0_WP
  SCV_src_phys  = 0.0_WP
  SCV_phys      = 0.0_WP
  ! Scalar variance budgets - conditional space
  allocate(SCV_unst_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(SCV_conv_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(SCV_div1_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(SCV_div2_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(SCV_prod_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(SCV_diss_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(SCV_diss2_cond(1:nbins_cond,pnmin_:pnmax_))
  allocate(SCV_src_cond (1:nbins_cond,pnmin_:pnmax_))
  allocate(SCV_cond     (1:nbins_cond,pnmin_:pnmax_))
  SCV_unst_cond = 0.0_WP
  SCV_conv_cond = 0.0_WP
  SCV_div1_cond = 0.0_WP
  SCV_div2_cond = 0.0_WP
  SCV_prod_cond = 0.0_WP
  SCV_diss_cond = 0.0_WP
  SCV_diss2_cond = 0.0_WP
  SCV_src_cond  = 0.0_WP
  SCV_cond      = 0.0_WP

  ! Temporary arrays
  allocate(tmpv_save_TKE_phys(1:nz,1:ny,pnmin_:pnmax_))
  allocate(tmpv_save_TKE_cond(1:nz,1:ny,pnmin_:pnmax_))
  tmpv_save_TKE_phys = 0.0_WP
  tmpv_save_TKE_cond = 0.0_WP
  allocate(tmpv_save_SCV_phys(1:nz,1:ny,pnmin_:pnmax_))
  allocate(tmpv_save_SCV_cond(1:nz,1:ny,pnmin_:pnmax_))
  tmpv_save_SCV_phys = 0.0_WP
  tmpv_save_SCV_cond = 0.0_WP

  return
end subroutine planeStats_budgets_init


subroutine planeStats_budgets_compute
  use planeStats_budgets
  implicit none
  
  integer :: iunit, ierr
  integer :: iplane, bin, i, j, k, isc, ii, jj
  real(WP) :: tmp1, tmp2
  real(WP), dimension(1:ny) :: tmpv1, tmpv2, tmpv3, tmpv4
  real(WP), dimension(1:nbins_cond) :: tmpc1, tmpc2, tmpc3, tmpc4
  real(WP), dimension(1:nz,1:ny) :: tmpxy

  real(WP), dimension(:,:), pointer :: o_Pm

  allocate(o_Pm(1:nplanes,1:ny))

  ! Compute the conditional mapping
  call cond_map

  ! Evaluate diffusivity, viscosity, and chemical source from finitechem
  !call planeStats_finitechem_diffusivity
  !call planeStats_finitechem_viscosity
  do iplane=pnmin_,pnmax_
     do j=1,ny
        do k=1,nz
           call planeStats_finitechem_source(SC(k,j,iplane,1:N_tot+1),src_SC_fc(k,j,iplane,1:N_tot+1))
        end do
     end do
  end do

  ! Mean pressure for this timestep
  do iplane=pnmin_,pnmax_
     call zmean(P(:,:,iplane), Pm(:,iplane))
  end do

  ! Compute the system mean pressure
  do j=1,ny
     call MPI_GATHER(Pm(j,:), nplanes_, MPI_REAL_WP, o_Pm(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  if (irank.eq.iroot) then
     P_mean = 0.0_WP
     do j=1,ny
        P_mean = P_mean + sum(o_Pm(:,j))/real(nplanes,WP)
     end do
     P_mean = P_mean/real(ny,WP)
     print *, 'P_mean = ', P_mean
  end if
  call MPI_BCAST(P_mean, 1, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)

  ! FLUCTUATING QUANTITIES
  ! P
  do iplane=pnmin_,pnmax_
     do j=1,ny
        ! Physical space fluctuations using global (spatial+temporal) mean
        Pf(:,j,iplane) = P(:,j,iplane) - P_mean
     end do
  end do

  ! U
  do ii=1,3
     do iplane=pnmin_,pnmax_
        do j=1,ny
           do k=1,nz
              ! Save instantaneous
              Ui(k,j,iplane,ii) = U(k,j,iplane,ii)
           end do
        end do
     end do
  end do
  do ii=1,3
     do iplane=pnmin_,pnmax_
        do j=1,ny
           do k=1,nz
              ! Conditional fluctuations
              U_c(k,j,iplane,ii) = U(k,j,iplane,ii) - Um_c(bin_index(k,j,iplane),iplane,ii)
           end do
        end do
     end do
  end do
  do ii=1,3
     do iplane=pnmin_,pnmax_
        do j=1,ny
           ! Spatial fluctuations
           U(:,j,iplane,ii) = U(:,j,iplane,ii) - Um(j,iplane,ii)
        end do
     end do
  end do

  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           tmpxy = RHO(:,:,iplane)*U(:,:,iplane,ii)*U(:,:,iplane,ii)*U(:,:,iplane,jj)
           call zmean(tmpxy,tmpv1)
           vel_triple_corr(:,iplane,ii,jj) = tmpv1/RHOm(:,iplane)
        end do
     end do
  end do
  

  ! Velocity gradient fluctuations
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           dUdxi(:,:,iplane,ii,jj) = dUdx(:,:,iplane,ii,jj)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           do j=1,ny
              do k=1,nz
                 dUdx_c(k,j,iplane,ii,jj) = dUdx(k,j,iplane,ii,jj) - dUdxm_c(bin_index(k,j,iplane),iplane,ii,jj)
              end do
           end do
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           do j=1,ny
              dUdx(:,j,iplane,ii,jj) = dUdx(:,j,iplane,ii,jj) - dUdxm(j,iplane,ii,jj)
           end do
        end do
     end do
  end do

!!$  ! Update conditional mean convection term in conditional mean velocity equation
!!$  do ii=1,3
!!$     do iplane=pnmin_,pnmax_
!!$        tmpxy = 0.0_WP
!!$        do jj=1,3
!!$           tmpxy = tmpxy &
!!$                + RHO(:,:,iplane)*Ui(:,:,iplane,jj)*dUdxi(:,:,iplane,ii,jj) &
!!$                - RHO(:,:,iplane)*U (:,:,iplane,jj)*dUdx (:,:,iplane,ii,jj)
!!$        end do
!!$        call cmean(iplane, tmpxy, cons_mom_m_cond(:,iplane,1,ii), icount_all)
!!$     end do
!!$  end do
!!$
!!$  ! Update fluctuating term in mean velocity equations
!!$  do ii=1,3
!!$     do iplane=pnmin_,pnmax_
!!$        tmpxy = 0.0_WP
!!$        do jj=1,3
!!$           tmpxy = tmpxy &
!!$                + U(:,:,iplane,ii)*U(:,:,iplane,jj)*dRHOdx(:,:,iplane,jj) &
!!$                + RHO(:,:,iplane)*U(:,:,iplane,jj)*dUdx(:,:,iplane,ii,jj) &
!!$                + RHO(:,:,iplane)*U(:,:,iplane,ii)*dUdx(:,:,iplane,jj,jj)
!!$        end do
!!$        call condmean(iplane, tmpxy, cons_mom_m_phys(:,iplane,2,ii), cons_mom_m_cond(:,iplane,2,ii), icount_all)
!!$     end do
!!$  end do

!!$  do ii=1,3
!!$     do iplane=pnmin_,pnmax_
!!$        tmpxy = 0.0_WP
!!$        do jj=1,3
!!$           tmpxy = tmpxy &
!!$                + U_c(:,:,iplane,ii)*U_c(:,:,iplane,jj)*dRHOdx(:,:,iplane,jj) &
!!$                + RHO(:,:,iplane)*U_c(:,:,iplane,jj)*dUdx_c(:,:,iplane,ii,jj) &
!!$                + RHO(:,:,iplane)*U_c(:,:,iplane,ii)*dUdx_c(:,:,iplane,jj,jj)
!!$        end do
!!$        call cmean(iplane, tmpxy, tmpc1)
!!$        cons_mom_m_cond(:,iplane,2,ii) = (real(ntime_curr-1,WP)*cons_mom_m_cond(:,iplane,2,ii) + tmpc1)/real(ntime_curr,WP)
!!$     end do
!!$  end do

     
  ! Fluctuating pressure gradient
  ! dPdxp : fluctuating
  ! dPdx  : stays instantaneous
  do jj=1,3
     do iplane=pnmin_,pnmax_
        do j=1,ny
           dPdxp(:,j,iplane,jj) = dPdx(:,j,iplane,jj) - dPdxm(j,iplane,jj)
        end do
     end do
  end do
  do jj=1,3
     do iplane=pnmin_,pnmax_
        do j=1,ny
           do k=1,nz
              dPdxp_c(k,j,iplane,jj) = dPdx(k,j,iplane,jj) - dPdxm_c(bin_index(k,j,iplane),iplane,jj)
           end do
        end do
     end do
  end do

  ! Fluctuating stress tensor
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           do j=1,ny
              dTAUdxp(:,j,iplane,ii,jj) = dTAUdx(:,j,iplane,ii,jj) - dTAUdxm(j,iplane,ii,jj)
           end do
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           do j=1,ny
              do k=1,nz
                 dTAUdxp_c(k,j,iplane,ii,jj) = dTAUdx(k,j,iplane,ii,jj) - dTAUdxm_c(bin_index(k,j,iplane),iplane,ii,jj)
              end do
           end do
        end do
     end do
  end do

  
  ! TKE BUDGETS
  call planeStats_budgets_velocity_compute


  ! SCALARS
  ! Scalar fluctuations
  do isc=1,nscalar
     do iplane=pnmin_,pnmax_
        do j=1,ny
           do k=1,nz
              ! Save instantaneous
              SCi(k,j,iplane,isc) = SC(k,j,iplane,isc)
           end do
        end do
     end do
  end do
  do isc=1,nscalar
     do iplane=pnmin_,pnmax_
        do j=1,ny
           do k=1,nz
              ! Conditional fluctuations
              SC_c(k,j,iplane,isc) = SC(k,j,iplane,isc) - SCm_c(bin_index(k,j,iplane),iplane,isc)
           end do
        end do
     end do
  end do
  do isc=1,nscalar
     do iplane=pnmin_,pnmax_
        do j=1,ny
           ! Spatial fluctuations
           SC(:,j,iplane,isc) = SC(:,j,iplane,isc) - SCm(j,iplane,isc)
        end do
     end do
  end do

  ! Fluctuating scalar gradients
  if (use_dSC) then
     do jj=1,3
        do isc=1,nscalar
           do iplane=pnmin_,pnmax_
              do j=1,ny
                 do k=1,nz
                    ! Save instantaneous
                    dSCdxi(k,j,iplane,isc,jj) = dSCdx(k,j,iplane,isc,jj)
                 end do
              end do
           end do
        end do
     end do
     do jj=1,3
        do isc=1,nscalar
           do iplane=pnmin_,pnmax_
              do j=1,ny
                 do k=1,nz
                    ! Conditional fluctuations
                    dSCdx_c(k,j,iplane,isc,jj) = dSCdx(k,j,iplane,isc,jj) - dSCdxm_c(bin_index(k,j,iplane),iplane,isc,jj)
                 end do
              end do
           end do
        end do
     end do
     do jj=1,3
        do isc=1,nscalar
           do iplane=pnmin_,pnmax_
              do j=1,ny
                 ! Spatial fluctuations
                 dSCdx(:,j,iplane,isc,jj) = dSCdx(:,j,iplane,isc,jj) - dSCdxm(j,iplane,isc,jj)
              end do
           end do
        end do
     end do
  end if

  
!!$  ! Update fluctuating and mean terms in mean scalar equations
!!$  if (use_qSC) then
!!$     do isc=1,nscalar-2
!!$        do iplane=pnmin_,pnmax_
!!$           ! Mean convection term (1)
!!$           ! Conditional
!!$           tmpxy = 0.0_WP
!!$           do jj=1,3
!!$              tmpxy = tmpxy &
!!$                   + RHO(:,:,iplane)*Ui(:,:,iplane,jj)*dSCdxi(:,:,iplane,isc,jj) &
!!$                   - RHO(:,:,iplane)*U (:,:,iplane,jj)*dSCdx (:,:,iplane,isc,jj)
!!$           end do
!!$           call cmean(iplane, tmpxy, cons_SC_m_cond(:,iplane,1,isc), icount_qSC)
!!$           
!!$           ! Fluctuating term (2)
!!$           ! Physical space
!!$           tmpxy = 0.0_WP
!!$           do jj=1,3
!!$              tmpxy = tmpxy &
!!$                   + U(:,:,iplane,jj)*SC(:,:,iplane,isc)*dRHOdx(:,:,iplane,jj) &
!!$                   + RHO(:,:,iplane)*U(:,:,iplane,jj)*dSCdx(:,:,iplane,isc,jj) &
!!$                   + RHO(:,:,iplane)*SC(:,:,iplane,isc)*dUdx(:,:,iplane,jj,jj)
!!$           end do
!!$           call condmean(iplane, tmpxy, cons_SC_m_phys(:,iplane,2,isc), cons_SC_m_cond(:,iplane,2,isc), icount_qSC)
!!$           
!!$           ! Diffusion term (3)
!!$           tmpxy = diff_flux_div(:,:,iplane,isc)
!!$           call condmean(iplane, tmpxy, cons_SC_m_phys(:,iplane,3,isc), cons_SC_m_cond(:,iplane,3,isc), icount_qSC)
!!$
!!$           ! Source term (4)
!!$           tmpxy = src_SC(:,:,iplane,isc)
!!$           call condmean(iplane, tmpxy, cons_SC_m_phys(:,iplane,4,isc), cons_SC_m_cond(:,iplane,4,isc), icount_qSC)
!!$        end do
!!$     end do
!!$  end if
  

  ! SCALAR VARIANCE BUDGETS
  call planeStats_budgets_scalar_compute


  ! REYNOLDS STRESS DISTRIBUTIONS
  ! Spatial means
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           tmpxy = RHO(:,:,iplane)*U(:,:,iplane,ii)*U(:,:,iplane,jj)
           call condmean(iplane, tmpxy, R_stress(:,iplane,ii,jj), R_stress_c(:,iplane,ii,jj), icount_all)
        end do
     end do
  end do
  ! TRUE conditional mean - using conditional fluctuations
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           tmpxy = RHO(:,:,iplane)*U_c(:,:,iplane,ii)*U_c(:,:,iplane,jj)
           call cmean(iplane, tmpxy, R_stress_tc(:,iplane,ii,jj), icount_all)
        end do
     end do
  end do


  ! SCALAR FLUX DISTRIBUTIONS
  ! Spatial means
  do jj=1,3
     do isc=1,nscalar
        do iplane=pnmin_,pnmax_
           tmpxy = RHO(:,:,iplane)*SC(:,:,iplane,isc)*U(:,:,iplane,jj)
           call condmean(iplane, tmpxy, SC_flux(:,iplane,isc,jj), SC_flux_c(:,iplane,isc,jj), icount_all)
        end do
     end do
  end do
  ! TRUE Conditional mean - using conditional fluctuations
  do jj=1,3
     do isc=1,nscalar
        do iplane=pnmin_,pnmax_
           tmpxy = RHO(:,:,iplane)*SC_c(:,:,iplane,isc)*U_c(:,:,iplane,jj)
           call cmean(iplane, tmpxy, SC_flux_tc(:,iplane,isc,jj), icount_all)
        end do
     end do
  end do


  ! SCALAR VARIANCE DISTRIBUTIONS
  ! Spatial means
  do isc=1,nscalar
     do iplane=pnmin_,pnmax_
        tmpxy = RHO(:,:,iplane)*SC(:,:,iplane,isc)**2
        call condmean(iplane, tmpxy, SC_var(:,iplane,isc), SC_var_c(:,iplane,isc), icount_all)
     end do
  end do
  ! TRUE Conditional mean - using conditional fluctuations
  do isc=1,nscalar
     do iplane=pnmin_,pnmax_
        tmpxy = RHO(:,:,iplane)*SC_c(:,:,iplane,isc)**2
        call cmean(iplane, tmpxy, SC_var_tc(:,iplane,isc), icount_all)
     end do
  end do

  return
end subroutine planeStats_budgets_compute


subroutine planeStats_budgets_output
  use planeStats_budgets
  implicit none

  integer :: iunit, iplane, bin, i, j, k, isc, ii, jj, ierr
  real(WP) :: w1, w2, tmp1, tmp2, tmp3, U_half!, U_c, U_s, alpha, y_alpha1, y_alpha9, U_prev, U_next
  real(WP) :: tau_shear, y_f, y_f10, y_f90, y_OH, y_T
  real(WP) :: w1_10, w2_10, w1_90, w2_90
  integer :: ii_10, ii_90

  ! Arrays for output
  real(WP), dimension(:,:),   pointer :: o_y_condm
  real(WP), dimension(:),   pointer :: r_half, U_cl
  real(WP), dimension(:,:), pointer :: o_nu_out, o_nu_c_out
  real(WP), dimension(:,:,:), pointer :: o_Um_out, o_Um_c_out, o_dRHOdxm, o_dRHOdxm_c
  ! Turbulence parameters
  real(WP), dimension(:,:), pointer :: k_TKE, u_rms, Re_L, L11, tau_L, eta, tau_eta, t_L
  ! Velocity and stress tensor gradients
  real(WP), dimension(:,:,:), pointer :: o_dUdxm_out_phys
  real(WP), dimension(:,:,:), pointer :: o_dUdxm_out_cond
  real(WP), dimension(:,:,:), pointer :: o_dTAUdxm_out
  ! TKE budgets - physical space
  real(WP), dimension(:,:), pointer :: o_E_conv_phys, o_E_div1_phys, o_E_div2_phys, o_E_div3_phys
  real(WP), dimension(:,:), pointer :: o_E_pres_phys, o_E_prod_phys, o_E_diss_phys, o_E_unst_phys
  real(WP), dimension(:,:), pointer :: o_E_TKE_phys, o_eta_phys, o_RHO_phys
  real(WP), dimension(:,:), pointer :: o_E_pwork_phys, o_E_ptrans_phys, o_E_pdil_phys
  ! TKE budgets - conditional space
  real(WP), dimension(:,:), pointer :: o_E_conv_cond, o_E_div1_cond, o_E_div2_cond, o_E_div3_cond
  real(WP), dimension(:,:), pointer :: o_E_pres_cond, o_E_prod_cond, o_E_diss_cond, o_E_unst_cond
  real(WP), dimension(:,:), pointer :: o_E_TKE_cond, o_eta_cond, o_RHO_cond
  ! Reynolds stress budget
  real(WP), dimension(:,:,:,:), pointer :: o_R_stress_budg
  ! Scalars
  real(WP), dimension(:,:,:), pointer :: o_SC_phys, o_SC_cond
  real(WP), dimension(:,:,:,:), pointer :: o_dSCdx_out_phys
  ! Scalar variance budgets - physical space
  real(WP), dimension(:,:), pointer :: o_SCV_conv_phys, o_SCV_div1_phys, o_SCV_div2_phys
  real(WP), dimension(:,:), pointer :: o_SCV_prod_phys, o_SCV_diss_phys, o_SCV_diss2_phys, o_SCV_unst_phys
  real(WP), dimension(:,:), pointer :: o_SCV_src_phys,  o_SCV_phys
  ! Scalar variance budgets - conditional space
  real(WP), dimension(:,:), pointer :: o_SCV_conv_cond, o_SCV_div1_cond, o_SCV_div2_cond
  real(WP), dimension(:,:), pointer :: o_SCV_prod_cond, o_SCV_diss_cond, o_SCV_diss2_cond, o_SCV_unst_cond
  real(WP), dimension(:,:), pointer :: o_SCV_src_cond,  o_SCV_cond

  ! Reynolds stresses
  real(WP), dimension(:,:,:,:), pointer :: o_R_stress, o_R_stress_c, o_R_stress_tc
  ! Scalar fluxes
  real(WP), dimension(:,:,:,:), pointer :: o_SC_flux, o_SC_flux_c, o_SC_flux_tc
  ! Scalar variance
  real(WP), dimension(:,:,:), pointer :: o_SC_var, o_SC_var_c, o_SC_var_tc

  ! Allocate output arrays on the root process
  !if (irank.eq.iroot) then
  allocate(o_y_condm      (1:nplanes,1:ny))
     allocate(r_half(1:nplanes))
     allocate(U_cl  (1:nplanes))
     allocate(o_Um_out  (1:3,1:nplanes,1:ny))
     allocate(o_Um_c_out(1:3,1:nplanes,1:nbins_cond))
     allocate(o_dRHOdxm  (1:3,1:nplanes,1:ny))
     allocate(o_dRHOdxm_c(1:3,1:nplanes,1:nbins_cond))
     allocate(o_nu_out  (1:nplanes,1:ny))
     allocate(o_nu_c_out  (1:nplanes,1:nbins_cond))
     ! Turbulence parameters
     allocate(k_TKE(1:nplanes,1:ny))
     allocate(Re_L(1:nplanes,1:ny))
     allocate(L11 (1:nplanes,1:ny))
     allocate(tau_L(1:nplanes,1:ny))
     allocate(t_L(1:nplanes,1:ny))
     allocate(eta(1:nplanes,1:ny))
     allocate(tau_eta(1:nplanes,1:ny))
     allocate(u_rms(1:nplanes,1:ny))
     ! Velocity and stress tensor gradients
     allocate(o_dUdxm_out_phys(1:9,1:nplanes,1:ny))
     allocate(o_dUdxm_out_cond(1:9,1:nplanes,1:nbins_cond))
     allocate(o_dTAUdxm_out(1:9,1:nplanes,1:ny))
     ! TKE budgets - physical space
     allocate(o_E_unst_phys(1:nplanes,1:ny))
     allocate(o_E_conv_phys(1:nplanes,1:ny))
     allocate(o_E_div1_phys(1:nplanes,1:ny))
     allocate(o_E_div2_phys(1:nplanes,1:ny))
     allocate(o_E_div3_phys(1:nplanes,1:ny))
     allocate(o_E_pres_phys(1:nplanes,1:ny))
     allocate(o_E_prod_phys(1:nplanes,1:ny))
     allocate(o_E_diss_phys(1:nplanes,1:ny))
     allocate(o_E_TKE_phys (1:nplanes,1:ny))
     allocate(o_E_pwork_phys (1:nplanes,1:ny))
     allocate(o_E_ptrans_phys(1:nplanes,1:ny))
     allocate(o_E_pdil_phys  (1:nplanes,1:ny))
     allocate(o_eta_phys   (1:nplanes,1:ny))
     allocate(o_RHO_phys   (1:nplanes,1:ny))
     ! TKE budgets - conditional space
     allocate(o_E_unst_cond(1:nplanes,1:nbins_cond))
     allocate(o_E_conv_cond(1:nplanes,1:nbins_cond))
     allocate(o_E_div1_cond(1:nplanes,1:nbins_cond))
     allocate(o_E_div2_cond(1:nplanes,1:nbins_cond))
     allocate(o_E_div3_cond(1:nplanes,1:nbins_cond))
     allocate(o_E_pres_cond(1:nplanes,1:nbins_cond))
     allocate(o_E_prod_cond(1:nplanes,1:nbins_cond))
     allocate(o_E_diss_cond(1:nplanes,1:nbins_cond))
     allocate(o_E_TKE_cond (1:nplanes,1:nbins_cond))
     allocate(o_eta_cond   (1:nplanes,1:nbins_cond))
     allocate(o_RHO_cond   (1:nplanes,1:nbins_cond))
     ! Reynolds stress budget
     allocate(o_R_stress_budg(1:nplanes,1:ny,1:6,1:10))
     ! Scalars
     allocate(o_SC_phys(1:nplanes,1:ny,1:nscalar))
     allocate(o_SC_cond(1:nplanes,1:nbins_cond,1:nscalar))
     allocate(o_dSCdx_out_phys (1:nplanes,1:ny,1:nscalar,1:3))
     ! Scalar variance budgets - physical space
     allocate(o_SCV_unst_phys(1:nplanes,1:ny))
     allocate(o_SCV_conv_phys(1:nplanes,1:ny))
     allocate(o_SCV_div1_phys(1:nplanes,1:ny))
     allocate(o_SCV_div2_phys(1:nplanes,1:ny))
     allocate(o_SCV_prod_phys(1:nplanes,1:ny))
     allocate(o_SCV_diss_phys(1:nplanes,1:ny))
     allocate(o_SCV_diss2_phys(1:nplanes,1:ny))
     allocate(o_SCV_src_phys (1:nplanes,1:ny))
     allocate(o_SCV_phys     (1:nplanes,1:ny))
     ! Scalar variance budgets - conditional space
     allocate(o_SCV_unst_cond(1:nplanes,1:nbins_cond))
     allocate(o_SCV_conv_cond(1:nplanes,1:nbins_cond))
     allocate(o_SCV_div1_cond(1:nplanes,1:nbins_cond))
     allocate(o_SCV_div2_cond(1:nplanes,1:nbins_cond))
     allocate(o_SCV_prod_cond(1:nplanes,1:nbins_cond))
     allocate(o_SCV_diss_cond(1:nplanes,1:nbins_cond))
     allocate(o_SCV_diss2_cond(1:nplanes,1:nbins_cond))
     allocate(o_SCV_src_cond (1:nplanes,1:nbins_cond))
     allocate(o_SCV_cond     (1:nplanes,1:nbins_cond))
     
     ! Reynolds stresses
     allocate(o_R_stress   (1:nplanes,1:ny,1:3,1:3))
     allocate(o_R_stress_c (1:nplanes,1:nbins_cond,1:3,1:3))
     allocate(o_R_stress_tc(1:nplanes,1:nbins_cond,1:3,1:3))
     ! Scalar fluxes
     allocate(o_SC_flux   (1:nplanes,1:ny,1:nscalar,1:3))
     allocate(o_SC_flux_c (1:nplanes,1:nbins_cond,1:nscalar,1:3))
     allocate(o_SC_flux_tc(1:nplanes,1:nbins_cond,1:nscalar,1:3))
     ! Scalar variance
     allocate(o_SC_var   (1:nplanes,1:ny,1:nscalar))
     allocate(o_SC_var_c (1:nplanes,1:nbins_cond,1:nscalar))
     allocate(o_SC_var_tc(1:nplanes,1:nbins_cond,1:nscalar))

  ! Gather the data to the root process
  do j=1,ny
     call MPI_GATHER(y_condm(j,:), nplanes_, MPI_REAL_WP, o_y_condm(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do ii=1,3
     do j=1,ny
        call MPI_GATHER(Um(j,:,ii), nplanes_, MPI_REAL_WP, o_Um_out(ii,:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,nbins_cond
        call MPI_GATHER(Um_c(j,:,ii), nplanes_, MPI_REAL_WP, o_Um_c_out(ii,:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do
  do ii=1,3
     do j=1,ny
        call MPI_GATHER(dRHOdxm(j,:,ii), nplanes_, MPI_REAL_WP, o_dRHOdxm(ii,:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,nbins_cond
        call MPI_GATHER(dRHOdxm_c(j,:,ii), nplanes_, MPI_REAL_WP, o_dRHOdxm_c(ii,:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do
  do j=1,ny
     call MPI_GATHER(NU(j,:), nplanes_, MPI_REAL_WP, o_nu_out(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(NU_c(j,:), nplanes_, MPI_REAL_WP, o_nu_c_out(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do

  k = 1
  do ii=1,3
     do jj=1,3
        do j=1,ny
           call MPI_GATHER(dUdxm(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_dUdxm_out_phys(k,:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
        k = k + 1
     end do
  end do
  k = 1
  do ii=1,3
     do jj=1,3
        do j=1,nbins_cond
           call MPI_GATHER(dUdxm_c(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_dUdxm_out_cond(k,:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
        k = k + 1
     end do
  end do

  k = 1
  do ii=1,3
     do jj=1,3
        do j=1,ny
           call MPI_GATHER(dTAUdxm(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_dTAUdxm_out(k,:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
        k = k + 1
     end do
  end do
  
  do j=1,ny
     call MPI_GATHER(E_unst_phys(j,:), nplanes_, MPI_REAL_WP, o_E_unst_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(E_conv_phys(j,:), nplanes_, MPI_REAL_WP, o_E_conv_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(E_div1_phys(j,:), nplanes_, MPI_REAL_WP, o_E_div1_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(E_div2_phys(j,:), nplanes_, MPI_REAL_WP, o_E_div2_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(E_div3_phys(j,:), nplanes_, MPI_REAL_WP, o_E_div3_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(E_pres_phys(j,:), nplanes_, MPI_REAL_WP, o_E_pres_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(E_prod_phys(j,:), nplanes_, MPI_REAL_WP, o_E_prod_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(E_diss_phys(j,:), nplanes_, MPI_REAL_WP, o_E_diss_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(E_TKE_phys(j,:), nplanes_, MPI_REAL_WP, o_E_TKE_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(E_pwork_phys(j,:), nplanes_, MPI_REAL_WP, o_E_pwork_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(E_ptrans_phys(j,:), nplanes_, MPI_REAL_WP, o_E_ptrans_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(E_pdil_phys(j,:), nplanes_, MPI_REAL_WP, o_E_pdil_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(RHOm(j,:), nplanes_, MPI_REAL_WP, o_RHO_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(eta_phys(j,:), nplanes_, MPI_REAL_WP, o_eta_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(E_unst_cond(j,:), nplanes_, MPI_REAL_WP, o_E_unst_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(E_conv_cond(j,:), nplanes_, MPI_REAL_WP, o_E_conv_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(E_div1_cond(j,:), nplanes_, MPI_REAL_WP, o_E_div1_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(E_div2_cond(j,:), nplanes_, MPI_REAL_WP, o_E_div2_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(E_div3_cond(j,:), nplanes_, MPI_REAL_WP, o_E_div3_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(E_pres_cond(j,:), nplanes_, MPI_REAL_WP, o_E_pres_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(E_prod_cond(j,:), nplanes_, MPI_REAL_WP, o_E_prod_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(E_diss_cond(j,:), nplanes_, MPI_REAL_WP, o_E_diss_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(E_TKE_cond(j,:), nplanes_, MPI_REAL_WP, o_E_TKE_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(eta_cond(j,:), nplanes_, MPI_REAL_WP, o_eta_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do jj=1,10
     do ii=1,6
        do j=1,ny
           call MPI_GATHER(R_stress_budg(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_R_stress_budg(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  
  do j=1,nbins_cond
     call MPI_GATHER(RHOm_c(j,:), nplanes_, MPI_REAL_WP, o_RHO_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do isc=1,nscalar
     do j=1,ny
        call MPI_GATHER(SCm(j,:,isc), nplanes_, MPI_REAL_WP, o_SC_phys(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do
  do isc=1,nscalar
     do j=1,nbins_cond
        call MPI_GATHER(SCm_c(j,:,isc), nplanes_, MPI_REAL_WP, o_SC_cond(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do

  do ii=1,3
     do isc=1,nscalar
        do j=1,ny
           call MPI_GATHER(dSCdxm(j,:,isc,ii), nplanes_, MPI_REAL_WP, o_dSCdx_out_phys(:,j,isc,ii), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do

  if (use_qSC) then
     do j=1,ny
        call MPI_GATHER(SCV_unst_phys(j,:), nplanes_, MPI_REAL_WP, o_SCV_unst_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,ny
        call MPI_GATHER(SCV_conv_phys(j,:), nplanes_, MPI_REAL_WP, o_SCV_conv_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,ny
        call MPI_GATHER(SCV_div1_phys(j,:), nplanes_, MPI_REAL_WP, o_SCV_div1_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,ny
        call MPI_GATHER(SCV_div2_phys(j,:), nplanes_, MPI_REAL_WP, o_SCV_div2_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,ny
        call MPI_GATHER(SCV_prod_phys(j,:), nplanes_, MPI_REAL_WP, o_SCV_prod_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,ny
        call MPI_GATHER(SCV_diss_phys(j,:), nplanes_, MPI_REAL_WP, o_SCV_diss_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,ny
        call MPI_GATHER(SCV_diss2_phys(j,:), nplanes_, MPI_REAL_WP, o_SCV_diss2_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,ny
        call MPI_GATHER(SCV_src_phys(j,:), nplanes_, MPI_REAL_WP, o_SCV_src_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,ny
        call MPI_GATHER(SCV_phys(j,:), nplanes_, MPI_REAL_WP, o_SCV_phys(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do

     do j=1,nbins_cond
        call MPI_GATHER(SCV_unst_cond(j,:), nplanes_, MPI_REAL_WP, o_SCV_unst_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,nbins_cond
        call MPI_GATHER(SCV_conv_cond(j,:), nplanes_, MPI_REAL_WP, o_SCV_conv_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,nbins_cond
        call MPI_GATHER(SCV_div1_cond(j,:), nplanes_, MPI_REAL_WP, o_SCV_div1_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,nbins_cond
        call MPI_GATHER(SCV_div2_cond(j,:), nplanes_, MPI_REAL_WP, o_SCV_div2_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,nbins_cond
        call MPI_GATHER(SCV_prod_cond(j,:), nplanes_, MPI_REAL_WP, o_SCV_prod_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,nbins_cond
        call MPI_GATHER(SCV_diss_cond(j,:), nplanes_, MPI_REAL_WP, o_SCV_diss_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,nbins_cond
        call MPI_GATHER(SCV_diss2_cond(j,:), nplanes_, MPI_REAL_WP, o_SCV_diss2_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,nbins_cond
        call MPI_GATHER(SCV_src_cond(j,:), nplanes_, MPI_REAL_WP, o_SCV_src_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
     do j=1,nbins_cond
        call MPI_GATHER(SCV_cond(j,:), nplanes_, MPI_REAL_WP, o_SCV_cond(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end if

  ! Reynolds stresses
  do jj=1,3
     do ii=1,3
        do j=1,ny
           call MPI_GATHER(R_stress(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_R_stress(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=1,nbins_cond
           call MPI_GATHER(R_stress_c(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_R_stress_c(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=1,nbins_cond
           call MPI_GATHER(R_stress_tc(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_R_stress_tc(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do

  ! Scalar fluxes
  do jj=1,3
     do isc=1,nscalar
        do j=1,ny
           call MPI_GATHER(SC_flux(j,:,isc,jj), nplanes_, MPI_REAL_WP, o_SC_flux(:,j,isc,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do isc=1,nscalar
        do j=1,nbins_cond
           call MPI_GATHER(SC_flux_c(j,:,isc,jj), nplanes_, MPI_REAL_WP, o_SC_flux_c(:,j,isc,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do isc=1,nscalar
        do j=1,nbins_cond
           call MPI_GATHER(SC_flux_tc(j,:,isc,jj), nplanes_, MPI_REAL_WP, o_SC_flux_tc(:,j,isc,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do

  ! Scalar variance
  do isc=1,nscalar
     do j=1,ny
        call MPI_GATHER(SC_var(j,:,isc), nplanes_, MPI_REAL_WP, o_SC_var(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do
  do isc=1,nscalar
     do j=1,nbins_cond
        call MPI_GATHER(SC_var_c(j,:,isc), nplanes_, MPI_REAL_WP, o_SC_var_c(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do
  do isc=1,nscalar
     do j=1,nbins_cond
        call MPI_GATHER(SC_var_tc(j,:,isc), nplanes_, MPI_REAL_WP, o_SC_var_tc(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do

  ! Only the root process writes output
  if (irank.ne.iroot) return

  ! Mirror physical space stats about centerline
  ! DON'T MIRROR VECTORS & TENSORS
!!$  do iplane=1,nplanes
!!$     do j=ny/2+1,ny-1
!!$        o_RHO_phys(iplane,j) = 0.5_WP*(o_RHO_phys(iplane,j) + o_RHO_phys(iplane,ny-j))
!!$        o_Um_out(1,iplane,j) = 0.5_WP*(o_Um_out(1,iplane,j) + o_Um_out(1,iplane,ny-j))
!!$        o_nu_out(iplane,j) = 0.5_WP*(o_nu_out(iplane,j) + o_nu_out(iplane,ny-j))
!!$        o_eta_phys(iplane,j) = 0.5_WP*(o_eta_phys(iplane,j) + o_eta_phys(iplane,ny-j))
!!$        o_E_TKE_phys (iplane,j) = 0.5_WP*(o_E_TKE_phys (iplane,j) + o_E_TKE_phys (iplane,ny-j))
!!$        o_E_conv_phys(iplane,j) = 0.5_WP*(o_E_conv_phys(iplane,j) + o_E_conv_phys(iplane,ny-j))
!!$        o_E_div2_phys(iplane,j) = 0.5_WP*(o_E_div2_phys(iplane,j) + o_E_div2_phys(iplane,ny-j))
!!$        o_E_div3_phys(iplane,j) = 0.5_WP*(o_E_div3_phys(iplane,j) + o_E_div3_phys(iplane,ny-j))
!!$        o_E_pres_phys(iplane,j) = 0.5_WP*(o_E_pres_phys(iplane,j) + o_E_pres_phys(iplane,ny-j))
!!$        o_E_prod_phys(iplane,j) = 0.5_WP*(o_E_prod_phys(iplane,j) + o_E_prod_phys(iplane,ny-j))
!!$        o_E_diss_phys(iplane,j) = 0.5_WP*(o_E_diss_phys(iplane,j) + o_E_diss_phys(iplane,ny-j))
!!$        o_E_unst_phys(iplane,j) = 0.5_WP*(o_E_unst_phys(iplane,j) + o_E_unst_phys(iplane,ny-j))
!!$        do isc=1,nscalar
!!$           o_SC_phys(iplane,j,isc) = 0.5_WP*(o_SC_phys(iplane,j,isc) + o_SC_phys(iplane,ny-j,isc))
!!$        end do
!!$        do ii=1,1
!!$           do jj=1,4
!!$              o_cons_mom_phys(iplane,j,jj,ii) = 0.5_WP*(o_cons_mom_phys(iplane,j,jj,ii) + o_cons_mom_phys(iplane,ny-j,jj,ii))
!!$           end do
!!$        end do
!!$        if (use_qSC) then
!!$           o_SCV_phys (iplane,j) = 0.5_WP*(o_SCV_phys (iplane,j) + o_SCV_phys (iplane,ny-j))
!!$           o_SCV_unst_phys(iplane,j) = 0.5_WP*(o_SCV_unst_phys(iplane,j) + o_SCV_unst_phys(iplane,ny-j))
!!$           o_SCV_conv_phys(iplane,j) = 0.5_WP*(o_SCV_conv_phys(iplane,j) + o_SCV_conv_phys(iplane,ny-j))
!!$           o_SCV_div2_phys(iplane,j) = 0.5_WP*(o_SCV_div2_phys(iplane,j) + o_SCV_div2_phys(iplane,ny-j))
!!$           o_SCV_prod_phys(iplane,j) = 0.5_WP*(o_SCV_prod_phys(iplane,j) + o_SCV_prod_phys(iplane,ny-j))
!!$           o_SCV_diss_phys(iplane,j) = 0.5_WP*(o_SCV_diss_phys(iplane,j) + o_SCV_diss_phys(iplane,ny-j))
!!$           o_SCV_diss2_phys(iplane,j) = 0.5_WP*(o_SCV_diss2_phys(iplane,j) + o_SCV_diss2_phys(iplane,ny-j))
!!$           o_SCV_src_phys(iplane,j) = 0.5_WP*(o_SCV_src_phys(iplane,j) + o_SCV_src_phys(iplane,ny-j))
!!$           do ii=1,4
!!$              o_cons_SC_phys(iplane,j,ii) = 0.5_WP*(o_cons_SC_phys(iplane,j,ii) + o_cons_SC_phys(iplane,ny-j,ii))
!!$           end do
!!$        end if
!!$     end do
!!$  end do

  ! Compute centerline axial velocity, jet half width, and turbulence stats

  print '(a12,ES15.5,a2)', 'Total time: ', time(ntime_curr)-time(1), ' s'
  open(unit=iunit,  file=trim(output_name)//'_centerline_stats', action='write')
  print       '(17a15)' ,'x/H', 'r1/2 (m)', 'U_cl (m/s)', 'u_rms (m/s)', 'TKE (m^2/s^2)', 'L_cl (mm)', 'tau_L (1/s)', 'Re_L', 'L times', 'eta_cl (um)', 'tau_eta (1/s)', 'tau_s12 (1/s)', 'y_f/(r1/2)', 'y_f10/(r1/2)', 'y_f90/(r1/2)', 'y_OH/(r1/2)', 'y_T/(r1/2)'
  write(iunit,'(17a15)'),'x/H', 'r1/2 (m)', 'U_cl (m/s)', 'u_rms (m/s)', 'TKE (m^2/s^2)', 'L_cl (mm)', 'tau_L (1/s)', 'Re_L', 'L times', 'eta_cl (um)', 'tau_eta (1/s)', 'tau_s12 (1/2)', 'y_f/(r1/2)', 'y_f10/(r1/2)', 'y_f90/(r1/2)', 'y_OH/(r1/2)', 'y_T/(r1/2)'

!!$  ! 1/9/17: Normalize by: (moot point if we're not running high-velocity coflows?)
!!$  !    Distance : Velocity deficit half-width
!!$  !    Velocity : Velocity deficit
!!$  ! Characteristic convection velocity and velocity deficit
!!$  U_c = 0.5_WP*(U_bulk + U_coflow)
!!$  U_s = abs(U_bulk - U_coflow)

  r_half = -1.0_WP
  do iplane=1,nplanes
!!$     ! Characteristic shear layer width
!!$     alpha = 0.1
!!$     do j=ny/2+1,ny-1
!!$        U_prev = o_Um_out(iplane,j)
!!$        U_next = o_Um_out(iplane,j+1)
!!$        if ((U_prev.lt.(U_coflow+alpha*(U_bulk-U_prev))) .and. (U_next.ge.(U_coflow+alpha*(U_bulk-U_next)))) then
!!$           
!!$    
!!$     
     U_cl(iplane) = o_Um_out(1,iplane,ny/2)
     U_half = 0.5_WP*U_cl(iplane)
     ! March upward to find the half-width
     do j=ny/2+1,ny-1
        if ((o_Um_out(1,iplane,j).gt.U_half).and.(o_Um_out(1,iplane,j+1).le.U_half)) then
           w1 = (U_half-o_Um_out(1,iplane,j+1))/(o_Um_out(1,iplane,j)-o_Um_out(1,iplane,j+1))
           w2 = 1.0_WP-w1
           r_half(iplane) = w1*y(j) + w2*y(j+1)
        end if
     end do
     ! If it wasn't found, jet velocity is increasing outward. Use 2*U_cl criterion instead.
     if (r_half(iplane).lt.0.0_WP) then
        print *, 'iplane=',iplane,' using 2*U_cl normalization'
        do j=ny/2+1,ny-1
           if ((o_Um_out(1,iplane,j).lt.4.0_WP*U_half).and.(o_Um_out(1,iplane,j+1).ge.U_half)) then
              w1 = (4.0_WP*U_half-o_Um_out(1,iplane,j+1))/(o_Um_out(1,iplane,j)-o_Um_out(1,iplane,j+1))
              w2 = 1.0_WP-w1
              r_half(iplane) = w1*y(j) + w2*y(j+1)
           end if
        end do
     end if

     ! Compute velocity and turbulence stats at centerline
     !k_TKE = 0.5_WP*o_E_TKE_phys(iplane,ny/2)/o_RHO_phys(iplane,ny/2)
     !u_rms = sqrt(2.0_WP/3.0_WP*k_TKE)
     ! Integral scale stats
     !L_cl  = sqrt(k_TKE**3)/abs(o_E_diss_phys(iplane,ny/2)/o_RHO_phys(iplane,ny/2))
     !tau_L = u_rms/L_cl
     do j=1,ny
        ! Turbulence and fluctuations
        k_TKE(iplane,j) = 0.5_WP*o_E_TKE_phys(iplane,j)/o_RHO_phys(iplane,j)
        u_rms(iplane,j) = sqrt(2.0_WP/3.0_WP*k_TKE(iplane,j))
        ! Integral scale
        L11(iplane,j) = sqrt(k_TKE(iplane,j)**3)/abs(o_E_diss_phys(iplane,j)/o_RHO_phys(iplane,j))
        tau_L(iplane,j) = u_rms(iplane,j)/L11(iplane,j)
        Re_L(iplane,j)  = sqrt(k_TKE(iplane,j))*L11(iplane,j)/o_nu_out(iplane,j)
        t_L(iplane,j)   = (time(ntime_curr)-time(1))*u_rms(iplane,j)/L11(iplane,j)
        ! Kolmogorov scale stats
        eta(iplane,j)  = o_eta_phys(iplane,j)
        tau_eta(iplane,j) = sqrt(o_nu_out(iplane,j)/abs(o_E_diss_phys(iplane,j)/o_RHO_phys(iplane,j)))
     end do
     ! Look for mean flame location
     if (trim(cond_var).eq.'ZMIX') then
        do j=ny/2,ny-1
!!$           ! Conserved scalar mixture fraction
!!$           tmp1 = o_SC_phys(iplane,j  ,isc_ZMIX)
!!$           tmp2 = o_SC_phys(iplane,j+1,isc_ZMIX)
           ! Bilger (1988) mixture fraction (Peters, Turbulent Combustion, p. 175)
           tmp3 = 0.077942_WP/1.008_WP + 2.0_WP*0.232_WP/16.0_WP
           tmp1 = (o_SC_phys(iplane,j  ,isc_H2)/1.008_WP + 2.0_WP*(0.232_WP-o_SC_phys(iplane,j  ,isc_O2))/16.0_WP)/tmp3
           tmp2 = (o_SC_phys(iplane,j+1,isc_H2)/1.008_WP + 2.0_WP*(0.232_WP-o_SC_phys(iplane,j+1,isc_O2))/16.0_WP)/tmp3
           if ((tmp1.gt.flameloc_cond).and.(tmp2.le.flameloc_cond)) then
              ii = j
              w2 = (flameloc_cond-tmp1)/(tmp2-tmp1)
              w1 = 1.0_WP-w2
           end if
           ! for now..
           ii_10 = ny/2; ii_90=ny/2
           w2_10 = 0.0_WP; w1_10 = 0.0_WP
           w2_90 = 0.0_WP; w2_10 = 0.0_WP
        end do
!!$        print *, w1, w2, w1*o_SC_phys(iplane,ii,isc_ZMIX) + w2*o_SC_phys(iplane,ii+1,isc_ZMIX)
     elseif (trim(cond_var).eq.'PROG') then
        !print *, prog_upper, prog_lower
        do j=ny/2,ny
           tmp1 = (prog_upper - o_SC_phys(iplane,j  ,isc_PROG))/(prog_upper-prog_lower)
           tmp2 = (prog_upper - o_SC_phys(iplane,j+1,isc_PROG))/(prog_upper-prog_lower)
           if ((tmp1.lt.flameloc_cond).and.(tmp2.ge.flameloc_cond)) then
              ii = j
              w2 = (flameloc_cond-tmp1)/(tmp2-tmp1)
              w1 = 1.0_WP-w2
           end if
           if ((tmp1.lt.0.1_WP).and.(tmp2.ge.0.1_WP)) then
              ii_10 = j
              w2_10 = (0.1_WP-tmp1)/(tmp2-tmp1)
              w1_10 = 1.0_WP-w2_10
           end if
           if ((tmp1.lt.0.9_WP).and.(tmp2.ge.0.9_WP)) then
              ii_90 = j
              w2_90 = (0.9_WP-tmp1)/(tmp2-tmp1)
              w1_90 = 1.0_WP-w2_90
           end if
        end do
!!$        tmp1 = (prog_upper - o_SC_phys(iplane,ii  ,isc_PROG))/(prog_upper-prog_lower)
!!$        tmp2 = (prog_upper - o_SC_phys(iplane,ii+1,isc_PROG))/(prog_upper-prog_lower)
!!$        print *, w1, w2, w1*tmp1 + w2*tmp2
!!$        tmp1 = (prog_upper - o_SC_phys(iplane,ii_10  ,isc_PROG))/(prog_upper-prog_lower)
!!$        tmp2 = (prog_upper - o_SC_phys(iplane,ii_10+1,isc_PROG))/(prog_upper-prog_lower)
!!$        print *, w1_10, w2_10, w1_10*tmp1 + w2_10*tmp2
!!$        tmp1 = (prog_upper - o_SC_phys(iplane,ii_90  ,isc_PROG))/(prog_upper-prog_lower)
!!$        tmp2 = (prog_upper - o_SC_phys(iplane,ii_90+1,isc_PROG))/(prog_upper-prog_lower)
!!$        print *, w1_90, w2_90, w1_90*tmp1 + w2_90*tmp2
     end if
     ! Shear time scale at mean flame location
     tau_shear = w1*abs(dUdxm(iplane,ii,1,2)) + w2*abs(dUdxm(iplane,ii+1,1,2))

     ! Mean flame location -- normalized by jet half width
     !print *, ny, ii, ii_10, ii_90
     y_f   = (w1*y(ii) + w2*y(ii+1))!/r_half(iplane)
     y_f10 = (w1_10*y(ii_10) + w2_10*y(ii_10+1))!/r_half(iplane)
     y_f90 = (w1_90*y(ii_90) + w2_90*y(ii_90+1))!/r_half(iplane)

     ! Max OH location - y_OH
     y_OH = y(maxloc(o_SC_phys(iplane,ny/2:ny,isc_OH),1)+ny/2-1)/r_half(iplane)
     ! Max temp location - y_T
     y_T  = y(maxloc(o_SC_phys(iplane,ny/2:ny,isc_T ),1)+ny/2-1)/r_half(iplane)
     ! Print
     print       '(17ES15.6)' , x(iplane)/H_jet, r_half(iplane), U_cl(iplane), u_rms(iplane,ny/2), k_TKE(iplane,ny/2), L11(iplane,ny/2)*1.0e3, tau_L(iplane,ny/2), Re_L(iplane,ny/2), t_L(iplane,ny/2), eta(iplane,ny/2)*1.0e6, tau_eta(iplane,ny/2), tau_shear, y_f, y_f10, y_f90, y_OH, y_T
     write(iunit,'(17ES15.6)'), x(iplane)/H_jet, r_half(iplane), U_cl(iplane), u_rms(iplane,ny/2), k_TKE(iplane,ny/2), L11(iplane,ny/2)*1.0e3, tau_L(iplane,ny/2), Re_L(iplane,ny/2), t_L(iplane,ny/2), eta(iplane,ny/2)*1.0e6, tau_eta(iplane,ny/2), tau_shear, y_f, y_f10, y_f90, y_OH, y_T
  end do
  close(iunit)


  ! Write output files
  ! Mean axial velocity -- raw and normalized
  open(unit=iunit, file=trim(output_name)//'_U',   action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)/H_jet
        end if
        if (iplane.eq.nplanes) then
           write(iunit,'(ES22.13,ES22.13,ES22.13)') &
                y(j)/r_half(iplane), o_Um_out(1,iplane,j)/U_bulk, o_Um_out(1,iplane,j)/U_cl(iplane)
        else
           write(iunit,'(ES22.13,ES22.13,ES22.13)',advance='no') &
                y(j)/r_half(iplane), o_Um_out(1,iplane,j)/U_bulk, o_Um_out(1,iplane,j)/U_cl(iplane)
        end if
     end do
  end do
  close(iunit)
  ! Gradients - instantaneous and fluctuating - include pressure gradient term
  open(unit=iunit, file=trim(output_name)//'_dUdx_phys',   action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        write(iunit,'(ES22.13, ES22.13)',advance='no') r_half(iplane), U_cl(iplane)
        do ii=1,9
           if (iplane.eq.nplanes .and. ii.eq.9) then
              write(iunit,'(ES22.13)') &
                   o_dUdxm_out_phys(ii,iplane,j)
           else
              write(iunit,'(ES22.13)',advance='no') &
                   o_dUdxm_out_phys(ii,iplane,j)
           end if
        end do
     end do
  end do
  close(iunit)
  open(unit=iunit, file=trim(output_name)//'_dUdx_cond',   action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        write(iunit,'(ES22.13, ES22.13)',advance='no') r_half(iplane), U_cl(iplane)
        do ii=1,9
           if (iplane.eq.nplanes .and. ii.eq.9) then
              write(iunit,'(ES22.13)') &
                   o_dUdxm_out_cond(ii,iplane,j)
           else
              write(iunit,'(ES22.13)',advance='no') &
                   o_dUdxm_out_cond(ii,iplane,j)
           end if
        end do
     end do
  end do
  close(iunit)
  open(unit=iunit, file=trim(output_name)//'_dTAUdx', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        write(iunit,'(ES22.13, ES22.13)',advance='no') r_half(iplane), U_cl(iplane)
        do ii=1,9
           if (iplane.eq.nplanes .and. ii.eq.9) then
              write(iunit,'(ES22.13)') o_dTAUdxm_out(ii,iplane,j)
           else
              write(iunit,'(ES22.13)',advance='no') o_dTAUdxm_out(ii,iplane,j)
           end if
        end do
     end do
  end do
  close(iunit)

  ! Density in physical space
  open(unit=iunit, file=trim(output_name)//'_RHO_phys',   action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        if (iplane.eq.nplanes) then
           write(iunit,'(ES22.13)') o_RHO_phys(iplane,j)
        else
           write(iunit,'(ES22.13)',advance='no') o_RHO_phys(iplane,j)
        end if
     end do
  end do
  close(iunit)
  ! Density in conditional space
  open(unit=iunit, file=trim(output_name)//'_RHO_cond',   action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        if (iplane.eq.nplanes) then
           write(iunit,'(ES22.13)') o_RHO_cond(iplane,j)
        else
           write(iunit,'(ES22.13)',advance='no') o_RHO_cond(iplane,j)
        end if
     end do
  end do
  close(iunit)

  ! Kolmogorov scales --- WORK THIS INTO TKE??
  open(unit=iunit, file=trim(output_name)//'_eta_phys',   action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        if (iplane.eq.nplanes) then
           write(iunit,'(8ES22.13)') &
                ! r_1/2, RHO, NU, epsilon, eta, L11, Re_L, Zvar
                r_half(iplane), o_RHO_phys(iplane,j), o_nu_out(iplane,j), abs(o_E_diss_phys(iplane,j)/o_RHO_phys(iplane,j)), o_eta_phys(iplane,j), L11(iplane,j), Re_L(iplane,j), o_y_condm(iplane,j)
        else
           write(iunit,'(8ES22.13)',advance='no') &
                r_half(iplane), o_RHO_phys(iplane,j), o_nu_out(iplane,j), abs(o_E_diss_phys(iplane,j)/o_RHO_phys(iplane,j)), o_eta_phys(iplane,j), L11(iplane,j), Re_L(iplane,j), o_y_condm(iplane,j)
        end if
     end do
  end do
  close(iunit)

  open(unit=iunit, file=trim(output_name)//'_eta_cond',   action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        if (iplane.eq.nplanes) then
           write(iunit,'(5ES22.13)') &
                ! r_1/2, RHO, NU, epsilon, eta
                r_half(iplane), o_RHO_cond(iplane,j), o_nu_c_out(iplane,j), abs(o_E_diss_cond(iplane,j)/o_RHO_cond(iplane,j)), o_eta_cond(iplane,j)
        else
           write(iunit,'(5ES22.13)',advance='no') &
                r_half(iplane), o_RHO_cond(iplane,j), o_nu_c_out(iplane,j), abs(o_E_diss_cond(iplane,j)/o_RHO_cond(iplane,j)), o_eta_cond(iplane,j)
        end if
     end do
  end do
  close(iunit)

  ! TKE budgets in physical space
  open(unit=iunit, file=trim(output_name)//'_budgets_phys', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        if (iplane.lt.nplanes) then
           write(iunit,'(12ES22.13)',advance='no') &
                r_half(iplane), U_cl(iplane), o_RHO_phys(iplane,j), &
                o_E_TKE_phys(iplane,j),  o_E_unst_phys(iplane,j), o_E_conv_phys(iplane,j), o_E_div2_phys(iplane,j), &
                o_E_div3_phys(iplane,j), o_E_pres_phys(iplane,j), o_E_prod_phys(iplane,j), o_E_diss_phys(iplane,j), o_y_condm(iplane,j)
        else
           write(iunit,'(12ES22.13)') &
                r_half(iplane), U_cl(iplane), o_RHO_phys(iplane,j), &
                o_E_TKE_phys(iplane,j),  o_E_unst_phys(iplane,j), o_E_conv_phys(iplane,j), o_E_div2_phys(iplane,j), &
                o_E_div3_phys(iplane,j), o_E_pres_phys(iplane,j), o_E_prod_phys(iplane,j), o_E_diss_phys(iplane,j), o_y_condm(iplane,j)
        end if
     end do
  end do
  close(iunit)
  ! TKE budgets in physical space - with pressure decomposition
  open(unit=iunit, file=trim(output_name)//'_budgets_phys_pres', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        if (iplane.lt.nplanes) then
           write(iunit,'(15ES22.13)',advance='no') &
                r_half(iplane), U_cl(iplane), o_RHO_phys(iplane,j), &
                o_E_TKE_phys(iplane,j),   o_E_unst_phys(iplane,j),   o_E_conv_phys(iplane,j), o_E_div2_phys(iplane,j), &
                o_E_div3_phys(iplane,j),  o_E_pres_phys(iplane,j),   o_E_prod_phys(iplane,j), o_E_diss_phys(iplane,j), &
                o_E_pwork_phys(iplane,j), o_E_ptrans_phys(iplane,j), o_E_pdil_phys(iplane,j), o_y_condm(iplane,j)
        else
           write(iunit,'(15ES22.13)') &
                r_half(iplane), U_cl(iplane), o_RHO_phys(iplane,j), &
                o_E_TKE_phys(iplane,j),   o_E_unst_phys(iplane,j),   o_E_conv_phys(iplane,j), o_E_div2_phys(iplane,j), &
                o_E_div3_phys(iplane,j),  o_E_pres_phys(iplane,j),   o_E_prod_phys(iplane,j), o_E_diss_phys(iplane,j), &
                o_E_pwork_phys(iplane,j), o_E_ptrans_phys(iplane,j), o_E_pdil_phys(iplane,j), o_y_condm(iplane,j)
        end if
     end do
  end do
  close(iunit)
  ! TKE budgets in conditional space
  open(unit=iunit, file=trim(output_name)//'_budgets_cond', action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        if (iplane.eq.nplanes) then
           write(iunit,'(ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13)') &
                r_half(iplane), U_cl(iplane), o_RHO_cond(iplane,j), &
                o_E_TKE_cond(iplane,j),  o_E_unst_cond(iplane,j), o_E_conv_cond(iplane,j), o_E_div2_cond(iplane,j), &
                o_E_div3_cond(iplane,j), o_E_pres_cond(iplane,j), o_E_prod_cond(iplane,j), o_E_diss_cond(iplane,j)
        else
           write(iunit,'(ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13)',advance='no') &
                r_half(iplane), U_cl(iplane), o_RHO_cond(iplane,j), &
                o_E_TKE_cond(iplane,j),  o_E_unst_cond(iplane,j), o_E_conv_cond(iplane,j), o_E_div2_cond(iplane,j), &
                o_E_div3_cond(iplane,j), o_E_pres_cond(iplane,j), o_E_prod_cond(iplane,j), o_E_diss_cond(iplane,j)
        end if
     end do
  end do
  close(iunit)

  ! Scalar gradients
  ! Need a switch starting here
  open(unit=iunit, file=trim(output_name)//'_dSCdx_phys', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        write (iunit,'(ES22.13)',advance='no') r_half(iplane)
        do ii=1,3
           do isc=1,nscalar
              if (iplane.eq.nplanes .and. isc.eq.nscalar .and. ii.eq.3) then
                 write(iunit,'(ES22.13)') &
                      o_dSCdx_out_phys(iplane,j,isc,ii)!, o_dSCdxp_out_phys(iplane,j,isc,ii), o_dSCdxm_out_phys(iplane,j,isc,ii)
              else
                 write(iunit,'(ES22.13)',advance='no') &
                      o_dSCdx_out_phys(iplane,j,isc,ii)!, o_dSCdxp_out_phys(iplane,j,isc,ii), o_dSCdxm_out_phys(iplane,j,isc,ii)
              end if
           end do
        end do
     end do
  end do
  close(iunit)
!!$  open(unit=iunit, file=trim(output_name)//'_dSCdx_cond', action='write')
!!$  do j=1,nbins_cond
!!$     do iplane=1,nplanes
!!$        if (iplane.eq.1) then
!!$           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
!!$        end if
!!$        write (iunit,'(ES22.13)',advance='no') r_half(iplane)
!!$        do ii=1,3
!!$           do isc=1,nscalar
!!$              if (iplane.eq.nplanes .and. isc.eq.nscalar .and. ii.eq.3) then
!!$                 write(iunit,'(3ES22.13)') &
!!$                      o_dSCdx_out_cond(iplane,j,isc,ii), o_dSCdxp_out_cond(iplane,j,isc,ii), o_dSCdxm_out_cond(iplane,j,isc,ii)
!!$              else
!!$                 write(iunit,'(3ES22.13)',advance='no') &
!!$                      o_dSCdx_out_cond(iplane,j,isc,ii), o_dSCdxp_out_cond(iplane,j,isc,ii), o_dSCdxm_out_cond(iplane,j,isc,ii)
!!$              end if
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$  close(iunit)


  ! Scalar variance budgets in physical space
  if (use_qSC) then
     open(unit=iunit, file=trim(output_name)//'_SCVbudgets_phys', action='write')
     do j=1,ny
        do iplane=1,nplanes
           if (iplane.eq.1) then
              write(iunit,'(ES22.13)',advance='no') y(j)
           end if
           if (iplane.lt.nplanes) then
              write(iunit,'(13ES22.13)',advance='no') &
                   r_half(iplane), U_cl(iplane), o_RHO_phys(iplane,j), &
                   o_SCV_phys(iplane,j),      o_SCV_unst_phys(iplane,j), o_SCV_conv_phys(iplane,j), o_SCV_div1_phys(iplane,j), o_SCV_div2_phys(iplane,j), &
                   o_SCV_prod_phys(iplane,j), o_SCV_diss_phys(iplane,j), o_SCV_diss2_phys(iplane,j), o_SCV_src_phys(iplane,j), o_y_condm(iplane,j)
           else
              write(iunit,'(13ES22.13)') &
                   r_half(iplane), U_cl(iplane), o_RHO_phys(iplane,j), &
                   o_SCV_phys(iplane,j),      o_SCV_unst_phys(iplane,j), o_SCV_conv_phys(iplane,j), o_SCV_div1_phys(iplane,j), o_SCV_div2_phys(iplane,j), &
                   o_SCV_prod_phys(iplane,j), o_SCV_diss_phys(iplane,j), o_SCV_diss2_phys(iplane,j), o_SCV_src_phys(iplane,j), o_y_condm(iplane,j)
           end if
        end do
     end do
     close(iunit)
     ! Scalar variance budgets in conditional space
     open(unit=iunit, file=trim(output_name)//'_SCVbudgets_cond', action='write')
     do j=1,nbins_cond
        do iplane=1,nplanes
           if (iplane.eq.1) then
              write(iunit,'(ES22.13)',advance='no') bins_cond(j)
           end if
           if (iplane.lt.nplanes) then
              write(iunit,'(12ES22.13)',advance='no') &
                   r_half(iplane), U_cl(iplane), o_RHO_cond(iplane,j), &
                   o_SCV_cond(iplane,j),      o_SCV_unst_cond(iplane,j), o_SCV_conv_cond(iplane,j), o_SCV_div1_cond(iplane,j), o_SCV_div2_cond(iplane,j), &
                   o_SCV_prod_cond(iplane,j), o_SCV_diss_cond(iplane,j), o_SCV_diss2_cond(iplane,j), o_SCV_src_cond(iplane,j)
           else
              write(iunit,'(12ES22.13)') &
                   r_half(iplane), U_cl(iplane), o_RHO_cond(iplane,j), &
                   o_SCV_cond(iplane,j),      o_SCV_unst_cond(iplane,j), o_SCV_conv_cond(iplane,j), o_SCV_div1_cond(iplane,j), o_SCV_div2_cond(iplane,j), &
                   o_SCV_prod_cond(iplane,j), o_SCV_diss_cond(iplane,j), o_SCV_diss2_cond(iplane,j), o_SCV_src_cond(iplane,j)
           end if
        end do
     end do
     close(iunit)
  end if


  ! Reynolds stresses
  open(unit=iunit, file=trim(output_name)//'_R_stress_phys', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        write (iunit,'(5ES22.13)',advance='no') o_y_condm(iplane,j), r_half(iplane), U_cl(iplane), o_RHO_phys(iplane,j), o_E_TKE_phys(iplane,j)
        do ii=1,3
           do jj=1,3
              if (iplane.eq.nplanes .and. ii.eq.3 .and. jj.eq.3) then
                 write(iunit,'(ES22.13)') o_R_stress(iplane,j,ii,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') o_R_stress(iplane,j,ii,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)
  open(unit=iunit, file=trim(output_name)//'_R_stress_cond', action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        write (iunit,'(4ES22.13)',advance='no') r_half(iplane), U_cl(iplane), o_RHO_cond(iplane,j), o_E_TKE_cond(iplane,j)
        do ii=1,3
           do jj=1,3
              if (iplane.eq.nplanes .and. ii.eq.3 .and. jj.eq.3) then
                 write(iunit,'(2ES22.13)') &
                      o_R_stress_c(iplane,j,ii,jj), o_R_stress_tc(iplane,j,ii,jj)
              else
                 write(iunit,'(2ES22.13)',advance='no') &
                      o_R_stress_c(iplane,j,ii,jj), o_R_stress_tc(iplane,j,ii,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)

  ! Reynolds stress budget
  open(unit=iunit, file=trim(output_name)//'_R_stress_budg', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        write (iunit,'(5ES22.13)',advance='no') o_y_condm(iplane,j), r_half(iplane), U_cl(iplane), o_RHO_phys(iplane,j), o_E_TKE_phys(iplane,j)
        do ii=1,6
           do jj=1,10
              if (iplane.eq.nplanes .and. ii.eq.3 .and. jj.eq.3) then
                 write(iunit,'(ES22.13)') o_R_stress_budg(iplane,j,ii,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') o_R_stress_budg(iplane,j,ii,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)

  
  ! Scalar fluxes
  open(unit=iunit, file=trim(output_name)//'_SC_flux_phys', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        write (iunit,'(5ES22.13)',advance='no') o_y_condm(iplane,j), r_half(iplane), U_cl(iplane), o_RHO_phys(iplane,j), o_E_TKE_phys(iplane,j)
        do isc=1,nscalar
           write(iunit,'(ES22.13)',advance='no') o_SC_phys(iplane,j,isc)
           do jj=1,3
              if (iplane.eq.nplanes .and. isc.eq.nscalar .and. jj.eq.3) then
                 write(iunit,'(ES22.13)') o_SC_flux(iplane,j,isc,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') o_SC_flux(iplane,j,isc,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)
  open(unit=iunit, file=trim(output_name)//'_SC_flux_cond', action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        write (iunit,'(4ES22.13)',advance='no') r_half(iplane), U_cl(iplane), o_RHO_cond(iplane,j), o_E_TKE_cond(iplane,j)
        do isc=1,nscalar
           write(iunit,'(ES22.13)',advance='no') o_SC_cond(iplane,j,isc)
           do jj=1,3
              if (iplane.eq.nplanes .and. isc.eq.nscalar .and. jj.eq.3) then
                 write(iunit,'(2ES22.13)') &
                      o_SC_flux_c(iplane,j,isc,jj), o_SC_flux_tc(iplane,j,isc,jj)
              else
                 write(iunit,'(2ES22.13)',advance='no') &
                      o_SC_flux_c(iplane,j,isc,jj), o_SC_flux_tc(iplane,j,isc,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)
  
  ! Scalar variance
  open(unit=iunit, file=trim(output_name)//'_SC_var_phys', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        write (iunit,'(5ES22.13)',advance='no') o_y_condm(iplane,j), r_half(iplane), U_cl(iplane), o_RHO_phys(iplane,j), o_E_TKE_phys(iplane,j)
        do isc=1,nscalar
           write(iunit,'(ES22.13)',advance='no') o_SC_phys(iplane,j,isc)
           if (iplane.eq.nplanes .and. isc.eq.nscalar) then
              write(iunit,'(ES22.13)') o_SC_var(iplane,j,isc)
           else
              write(iunit,'(ES22.13)',advance='no') o_SC_var(iplane,j,isc)
           end if
        end do
     end do
  end do
  close(iunit)
  open(unit=iunit, file=trim(output_name)//'_SC_var_cond', action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        write (iunit,'(4ES22.13)',advance='no') r_half(iplane), U_cl(iplane), o_RHO_cond(iplane,j), o_E_TKE_cond(iplane,j)
        do isc=1,nscalar
           write(iunit,'(ES22.13)',advance='no') o_SC_cond(iplane,j,isc)
           if (iplane.eq.nplanes .and. isc.eq.nscalar) then
              write(iunit,'(2ES22.13)') &
                   o_SC_var_c(iplane,j,isc), o_SC_var_tc(iplane,j,isc)
           else
              write(iunit,'(2ES22.13)',advance='no') &
                   o_SC_var_c(iplane,j,isc), o_SC_var_tc(iplane,j,isc)
           end if
        end do
     end do
  end do
  close(iunit)

  ! For testing
  open(unit=iunit, file=trim(output_name)//'_vel_triple_corr', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        do ii=1,3
           do jj=1,3
              if (iplane.eq.nplanes .and. ii.eq.3 .and. jj.eq.3) then
                 write(iunit,'(ES22.13)') vel_triple_corr(iplane,j,ii,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') vel_triple_corr(j,iplane,ii,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)

  print *, 'Output done.'
  
  return
end subroutine planeStats_budgets_output
