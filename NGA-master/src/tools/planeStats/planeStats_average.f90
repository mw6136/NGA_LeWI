module planeStats_average
  use planeStats

  real(WP), dimension(:,:),     pointer :: y_condm

  real(WP), dimension(:,:,:,:), pointer :: dUdx_tmpv, dSCdx_tmpv
  real(WP), dimension(:,:,:,:), pointer :: dUdx_tmpc, dSCdx_tmpc

  ! For scalar flux and Reynolds stress tensor
  real(WP), dimension(:,:,:,:), pointer :: rhoUiUjm, rhoUiUjm_c, dUiUjm
  real(WP), dimension(:,:,:,:), pointer :: rhoUiSCm, rhoUiSCm_c, dUjSCm

!!$  ! Alignment
!!$  real(WP), dimension(:,:,:), pointer :: align_SCF_SCG
!!$  real(WP), dimension(:,:,:), pointer :: align_RST_SRT
!!$
!!$  ! Eigenvalues
!!$  real(WP), dimension(:,:,:), pointer :: eval_S,   eval_R
!!$  real(WP), dimension(:,:,:), pointer :: eval_S_d, eval_R_d

  ! Trace
  real(WP), dimension(:,:), pointer :: trace_S, trace_R, diss_rate
  real(WP), dimension(:,:), pointer :: trace_S_c, trace_R_c

  ! For dissipation rate
  real(WP), dimension(:,:,:,:), pointer :: taum
  real(WP), dimension(:,:), pointer :: tau_dUdxm

  ! RMS
  real(WP), dimension(:,:), pointer :: rms_S_local, rms_R_local
  
  ! For model evaluation
  real(WP), dimension(:,:,:), pointer :: DIFFm, src_SCm, diff_flux_divm
  real(WP), dimension(:,:,:), pointer :: DIFFm_fc, src_SCm_fc, src_SCm_test
  real(WP), dimension(:,:),   pointer :: VISCm_fc

  ! For verification
  real(WP), dimension(:,:,:,:), pointer :: tmpv_save_mom
  real(WP), dimension(:,:,:),   pointer :: tmpv_save_SC
  real(WP), dimension(:,:,:,:), pointer :: cons_mom_m_phys!, cons_mom_m_cond
  real(WP), dimension(:,:,:,:), pointer :: cons_SC_m_phys!,  cons_SC_m_cond

  ! Full & anisotropic Reynolds stress & strain rate tensors - compute in Round 2
  real(WP), dimension(:,:,:,:), pointer :: tau_ij, tau_ij_d, S_ij, S_ij_d
  real(WP), dimension(:,:,:,:), pointer :: tau_ij_c, tau_ij_dc, S_ij_c, S_ij_dc
  real(WP), dimension(:), pointer :: rms_tau, rms_Sij
  real(WP) :: min_tau, max_tau, min_Sij, max_Sij

  ! Reynolds stresses
  real(WP), dimension(:,:,:,:), pointer :: R_stress, R_stress_c, R_stress_tc
  ! Scalar fluxes
  real(WP), dimension(:,:,:,:), pointer :: SC_flux, SC_flux_c, SC_flux_tc
  ! Scalar variance
  real(WP), dimension(:,:,:), pointer :: SC_var, SC_var_c, SC_var_tc

end module planeStats_average

! ========================================================== !
! Initialize the module                                      !
! ========================================================== !
subroutine planeStats_average_init
  use planeStats_average
  implicit none

  call planeStats_finitechem_init

  ! Located in module planeStats

  allocate(y_condm(1:ny,pnmin_:pnmax_))
  allocate(dTAUdxm  (1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(dTAUdxm_c(1:nbins_cond,pnmin_:pnmax_,1:3,1:3))
  dTAUdxm   = 0.0_WP
  dTAUdxm_c = 0.0_WP
  allocate(dPdxm    (1:ny,pnmin_:pnmax_,1:3))
  allocate(dPdxm_c  (1:nbins_cond,pnmin_:pnmax_,1:3))
  dPdxm   = 0.0_WP
  dPdxm_c = 0.0_WP
  allocate(SCm  (1:ny,pnmin_:pnmax_,1:nscalar))
  allocate(SCm_c(1:nbins_cond,pnmin_:pnmax_,1:nscalar))
  SCm   = 0.0_WP
  SCm_c = 0.0_WP
  allocate(NU    (1:ny,pnmin_:pnmax_))
  allocate(Um    (1:ny,pnmin_:pnmax_,1:3))
  allocate(RHOm  (1:ny,pnmin_:pnmax_))
  allocate(Pm    (1:ny,pnmin_:pnmax_))
  allocate(dUdxm (1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(dSCdxm(1:ny,pnmin_:pnmax_,1:nscalar,1:3))
  allocate(dRHOdxm  (1:ny,pnmin_:pnmax_,1:3))

  allocate(dRHOdxm_c(1:nbins_cond,pnmin_:pnmax_,1:3))
  allocate(NU_c     (1:nbins_cond,pnmin_:pnmax_))
  allocate(RHOm_c   (1:nbins_cond,pnmin_:pnmax_))
  allocate(Pm_c     (1:nbins_cond,pnmin_:pnmax_))
  allocate(Um_c     (1:nbins_cond,pnmin_:pnmax_,1:3))
  allocate(dUdxm_c  (1:nbins_cond,pnmin_:pnmax_,1:3,1:3))
  allocate(dSCdxm_c (1:nbins_cond,pnmin_:pnmax_,1:nscalar,1:3))

  allocate(dUdx_tmpv (1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(dSCdx_tmpv(1:ny,pnmin_:pnmax_,1:nscalar,1:3))
  allocate(dUdx_tmpc (1:nbins_cond,pnmin_:pnmax_,1:3,1:3))
  allocate(dSCdx_tmpc(1:nbins_cond,pnmin_:pnmax_,1:nscalar,1:3))

  allocate(dUiUjm  (1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(dUjSCm  (1:ny,pnmin_:pnmax_,1:3,1:nscalar))

  allocate(rhoUiUjm  (1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(rhoUiUjm_c(1:nbins_cond,pnmin_:pnmax_,1:3,1:3))
  allocate(rhoUiSCm  (1:ny,pnmin_:pnmax_,1:nscalar,1:3))
  allocate(rhoUiSCm_c(1:nbins_cond,pnmin_:pnmax_,1:nscalar,1:3))

!!$  allocate(align_SCF_SCG(1:ny,pnmin_:pnmax_,1:nscalar))
!!$  allocate(align_RST_SRT(1:ny,pnmin_:pnmax_,1:3))
!!$
!!$  allocate(eval_S  (1:ny,pnmin_:pnmax_,1:3))
!!$  allocate(eval_R  (1:ny,pnmin_:pnmax_,1:3))
!!$  allocate(eval_S_d(1:ny,pnmin_:pnmax_,1:3))
!!$  allocate(eval_R_d(1:ny,pnmin_:pnmax_,1:3))

  allocate(trace_S(1:ny,pnmin_:pnmax_))
  allocate(trace_R(1:ny,pnmin_:pnmax_))
  allocate(diss_rate(1:ny,pnmin_:pnmax_))
  allocate(trace_S_c(1:nbins_cond,pnmin_:pnmax_))
  allocate(trace_R_c(1:nbins_cond,pnmin_:pnmax_))

  allocate(tau_dUdxm(1:ny,pnmin_:pnmax_))
  allocate(taum(1:ny,pnmin_:pnmax_,1:3,1:3))
  tau_dUdxm = 0.0_WP
  taum = 0.0_WP

  allocate(rms_S_local(1:ny,pnmin_:pnmax_))
  allocate(rms_R_local(1:ny,pnmin_:pnmax_))

  allocate(DIFFm     (1:ny,pnmin_:pnmax_,1:nscalar))
  allocate(DIFFm_fc  (1:ny,pnmin_:pnmax_,1:nscalar))
  allocate(src_SCm   (1:ny,pnmin_:pnmax_,1:nscalar-1))
  allocate(src_SCm_fc(1:ny,pnmin_:pnmax_,1:nscalar-1))
  allocate(src_SCm_test(1:ny,pnmin_:pnmax_,1:nscalar-1))
  allocate(VISCm_fc  (1:ny,pnmin_:pnmax_))
  allocate(diff_flux_divm(1:ny,pnmin_:pnmax_,1:nscalar))

  ! Reynolds stresses
  allocate(R_stress   (1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(R_stress_c (1:nbins_cond,pnmin_:pnmax_,1:3,1:3))
  allocate(R_stress_tc(1:nbins_cond,pnmin_:pnmax_,1:3,1:3))
  ! Scalar fluxes
  allocate(SC_flux   (1:ny,pnmin_:pnmax_,1:nscalar,1:3))
  allocate(SC_flux_c (1:nbins_cond,pnmin_:pnmax_,1:nscalar,1:3))
  allocate(SC_flux_tc(1:nbins_cond,pnmin_:pnmax_,1:nscalar,1:3))
  ! Scalar variance
  allocate(SC_var   (1:ny,pnmin_:pnmax_,1:nscalar))
  allocate(SC_var_c (1:nbins_cond,pnmin_:pnmax_,1:nscalar))
  allocate(SC_var_tc(1:nbins_cond,pnmin_:pnmax_,1:nscalar))

  ! Deviatoric Reynolds stress tensor & strain rate
  allocate(tau_ij   (1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(tau_ij_d (1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(tau_ij_c (1:nbins_cond,pnmin_:pnmax_,1:3,1:3))
  allocate(tau_ij_dc(1:nbins_cond,pnmin_:pnmax_,1:3,1:3))
  allocate(S_ij     (1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(S_ij_d   (1:ny,pnmin_:pnmax_,1:3,1:3))
  allocate(S_ij_c   (1:nbins_cond,pnmin_:pnmax_,1:3,1:3))
  allocate(S_ij_dc  (1:nbins_cond,pnmin_:pnmax_,1:3,1:3))
  tau_ij    = 0.0_WP
  tau_ij_d  = 0.0_WP
  tau_ij_dc = 0.0_WP
  tau_ij_c  = 0.0_WP
  S_ij      = 0.0_WP
  S_ij_d    = 0.0_WP
  S_ij_dc   = 0.0_WP
  S_ij_c    = 0.0_WP
  allocate(rms_tau(pnmin_:pnmax_))
  allocate(rms_Sij(pnmin_:pnmax_))

  ! For verification
  allocate(tmpv_save_mom(1:nz,1:ny,pnmin_:pnmax_,1:3))
  allocate(tmpv_save_SC (1:nz,1:ny,pnmin_:pnmax_))

  allocate(cons_mom_m_phys(1:ny,        pnmin_:pnmax_,1:4,1:3))
!!$  allocate(cons_mom_m_cond(1:nbins_cond,pnmin_:pnmax_,1:4,1:3))
  cons_mom_m_phys = 0.0_WP
!!$  cons_mom_m_cond = 0.0_WP
  allocate(cons_SC_m_phys (1:ny,        pnmin_:pnmax_,1:4,1:nscalar))
!!$  allocate(cons_SC_m_cond (1:nbins_cond,pnmin_:pnmax_,1:4,1:nscalar))
  cons_SC_m_phys  = 0.0_WP
!!$  cons_SC_m_cond  = 0.0_WP

  return
end subroutine planeStats_average_init


! ========================================================== !
! Compute the averages of all variables in space (physical   !
!    and conditional) and time                               !
! ========================================================== !
subroutine planeStats_average_all
  use planeStats_average
  implicit none

  integer :: iplane, ii, jj, isc, j, k
  real(WP), dimension(1:nz,1:ny) :: tmpxy, tau
  real(WP), dimension(1:ny) :: tmpv1, tmpv2
  real(WP), dimension(1:nbins_cond) :: tmpc1, tmpc2

  ! Compute the conditional mapping
  call cond_map

  ! Evaluate diffusivity, viscosity, and chemical source from finitechem
  call planeStats_finitechem_diffusivity
  call planeStats_finitechem_viscosity
  do iplane=pnmin_,pnmax_
     do j=1,ny
        do k=1,nz
           call planeStats_finitechem_source(SC(k,j,iplane,1:N_tot+1),src_SC_fc(k,j,iplane,1:N_tot+1))
        end do
     end do
  end do
  
  ! Compute the mean progress variable
  do iplane=pnmin_,pnmax_
     call zmean(RHO(:,:,iplane)*y_cond(:,:,iplane), tmpv1)
     y_condm(:,iplane) = (real(ntime_curr-1,WP)*y_condm(:,iplane) + tmpv1)/real(ntime_curr,WP)
  end do

  ! Compute the average density
  do iplane=pnmin_,pnmax_
     call condmean(iplane, RHO(:,:,iplane), RHOm(:,iplane), RHOm_c(:,iplane), icount_all)
  end do

  ! Averaged density gradients
  do ii=1,3
     do iplane=pnmin_,pnmax_
        call condmean(iplane, dRHOdx(:,:,iplane,ii), dRHOdxm(:,iplane,ii), dRHOdxm_c(:,iplane,ii), icount_all)
     end do
  end do

  ! VISC
  do iplane=pnmin_,pnmax_
     call condmean(iplane, VISC(:,:,iplane), NU(:,iplane), NU_c(:,iplane), icount_all)
  end do
  do iplane=pnmin_,pnmax_
     call zmean(VISC_fc(:,:,iplane), tmpv1)
     VISCm_fc(:,iplane) = (real(ntime_curr-1,WP)*VISCm_fc(:,iplane) + tmpv1)/real(ntime_curr,WP)
  end do
  
  ! U
  do ii=1,3
     do iplane=pnmin_,pnmax_
        tmpxy = RHO(:,:,iplane)*U(:,:,iplane,ii)
        call condmean(iplane, tmpxy, Um(:,iplane,ii), Um_c(:,iplane,ii), icount_all)
     end do
  end do

  ! Compute the average pressure on the planes and ny
  !   In average_finalize, compute the mean pressure across the system (Pmean)
  do iplane=pnmin_,pnmax_
     call condmean(iplane, P(:,:,iplane), Pm(:,iplane), Pm_c(:,iplane), icount_all)
  end do
  
  ! Average pressure gradients
  do ii=1,3
     do iplane=pnmin_,pnmax_
        call condmean(iplane, dPdx(:,:,iplane,ii), dPdxm(:,iplane,ii), dPdxm_c(:,iplane,ii), icount_all)
     end do
  end do
  
  ! Stress tensor gradients
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           call condmean(iplane, dTAUdx(:,:,iplane,ii,jj), dTAUdxm(:,iplane,ii,jj), dTAUdxm_c(:,iplane,ii,jj), icount_all)
        end do
     end do
  end do

  ! Velocity gradients
  ! All constituent quantities are instantaneous
  ! Compute components now; update average gradients in finalize
  ! RHO*U in Um
  ! d(RHO*U)
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           tmpxy = RHO(:,:,iplane)*dUdx(:,:,iplane,ii,jj) + U(:,:,iplane,ii)*dRHOdx(:,:,iplane,jj)
           call condmean(iplane, tmpxy, dUdx_tmpv(:,iplane,ii,jj), dUdx_tmpc(:,iplane,ii,jj), icount_all)
        end do
     end do
  end do

  ! Scalars
  do isc=1,nscalar
     do iplane=pnmin_,pnmax_
        tmpxy = RHO(:,:,iplane)*SC(:,:,iplane,isc)
        call condmean(iplane, tmpxy, SCm(:,iplane,isc), SCm_c(:,iplane,isc), icount_all)
     end do
  end do

  ! Scalar gradients
  ! RHO*SC in SCm
  if (use_dSC) then
     ! d(RHO*SC)
     do jj=1,3
        do isc=1,nscalar
           do iplane=pnmin_,pnmax_
              tmpxy = RHO(:,:,iplane)*dSCdx(:,:,iplane,isc,jj) + SC(:,:,iplane,isc)*dRHOdx(:,:,iplane,jj)
              call condmean(iplane, tmpxy, dSCdx_tmpv(:,iplane,isc,jj), dSCdx_tmpc(:,iplane,isc,jj), icount_dSC)
           end do
        end do
     end do
  end if

  ! For scalar flux
  ! RHO*U_i*SC
  do isc=1,nscalar
     do ii=1,3
        do iplane=pnmin_,pnmax_
           tmpxy = RHO(:,:,iplane)*U(:,:,iplane,ii)*SC(:,:,iplane,isc)
           call condmean(iplane, tmpxy, rhoUiSCm(:,iplane,isc,ii), rhoUiSCm_c(:,iplane,isc,ii), icount_all)
        end do
     end do
  end do

  ! For divergence of scalar flux
  do isc=1,nscalar
     do ii=1,3
        do iplane=pnmin_,pnmax_
           tmpxy = RHO(:,:,iplane) * &
                ( U(:,:,iplane,ii)*dSCdx(:,:,iplane,isc,ii) &
                + SC(:,:,iplane,isc)*dUdx(:,:,iplane,ii,ii) )
           call zmean(tmpxy, tmpv1)
           dUjSCm(:,iplane,ii,isc) = (real(ntime_curr-1,WP)*dUjSCm(:,iplane,ii,isc) + tmpv1)/real(ntime_curr,WP)
        end do
     end do
  end do

  ! For Reynolds stress tensor
  ! RHO*U_i*U_j
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           tmpxy = RHO(:,:,iplane)*U(:,:,iplane,ii)*U(:,:,iplane,jj)
           call condmean(iplane, tmpxy, rhoUiUjm(:,iplane,ii,jj), rhoUiUjm_c(:,iplane,ii,jj), icount_all)
        end do
     end do
  end do

  ! For divergence of Reynolds stress
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           tmpxy = RHO(:,:,iplane) * & 
                ( U(:,:,iplane,ii)*dUdx(:,:,iplane,jj,jj) &
                + U(:,:,iplane,jj)*dUdx(:,:,iplane,ii,jj) )
           call zmean(tmpxy, tmpv1)
           dUiUjm(:,iplane,ii,jj) = (real(ntime_curr-1,WP)*dUiUjm(:,iplane,ii,jj) + tmpv1)/real(ntime_curr,WP)
        end do
     end do
  end do

  ! For dissipation rate
  do iplane=pnmin_,pnmax_
     tmpxy = 0.0_WP
     do jj=1,3
        do ii=1,3
           if (ii.eq.jj) then
              tau = VISC(:,:,iplane)*( dUdx(:,:,iplane,ii,jj) + dUdx(:,:,iplane,jj,ii) &
                   - 2.0_WP/3.0_WP*(dUdx(:,:,iplane,1,1) + dUdx(:,:,iplane,2,2) + dUdx(:,:,iplane,3,3)) )
           else
              tau = VISC(:,:,iplane)*( dUdx(:,:,iplane,ii,jj) + dUdx(:,:,iplane,jj,ii) )
           end if
           ! Mean stress tensor
           call zmean(tau, tmpv1)
           taum(:,iplane,ii,jj) = (real(ntime_curr-1,WP)*taum(:,iplane,ii,jj) + tmpv1)/real(ntime_curr,WP)
           ! Mean product of stress tensor and velocity gradient
           tmpxy = tmpxy + tau*dUdx(:,:,iplane,ii,jj)
        end do
     end do
     call zmean(tmpxy, tmpv1)
     tau_dUdxm(:,iplane) = (real(ntime_curr-1,WP)*tau_dUdxm(:,iplane) + tmpv1)/real(ntime_curr,WP)
  end do

  ! For model evaluation - mean source terms and diffusivities
  if (use_dSC) then
     ! Source terms
     do isc=1,nscalar-1
        do iplane=pnmin_,pnmax_
           call zmean(src_SC(:,:,iplane,isc), tmpv1)
           src_SCm(:,iplane,isc) = (real(ntime_dSC-1,WP)*src_SCm(:,iplane,isc) + tmpv1)/real(ntime_dSC,WP)
        end do
     end do
  end if
  if (use_qSC) then
     ! Diffusivities and diffusive fluxes
     do isc=1,nscalar
        do iplane=pnmin_,pnmax_
           call zmean(DIFF(:,:,iplane,isc), tmpv1)
           DIFFm(:,iplane,isc) = (real(ntime_qSC-1,WP)*DIFFm(:,iplane,isc) + tmpv1)/real(ntime_qSC,WP)
        end do
     end do
     do isc=1,nscalar
        do iplane=pnmin_,pnmax_
           call zmean(diff_flux_div(:,:,iplane,isc), tmpv1)
           diff_flux_divm(:,iplane,isc) = (real(ntime_qSC-1,WP)*diff_flux_divm(:,iplane,isc) + tmpv1)/real(ntime_qSC,WP)
        end do
     end do
  end if

  ! Source terms from finitechem
  do isc=1,nscalar-1
     do iplane=pnmin_,pnmax_
        call zmean(src_SC_fc(:,:,iplane,isc), tmpv1)
        src_SCm_fc(:,iplane,isc) = (real(ntime_curr-1,WP)*src_SCm_fc(:,iplane,isc) + tmpv1)/real(ntime_curr,WP)
     end do
  end do
  ! Diffusivities from finitechem
  do isc=1,nscalar
     do iplane=pnmin_,pnmax_
        call zmean(DIFF_fc(:,:,iplane,isc), tmpv1)
        DIFFm_fc(:,iplane,isc) = (real(ntime_curr-1,WP)*DIFFm_fc(:,iplane,isc) + tmpv1)/real(ntime_curr,WP)
     end do
  end do

  return
end subroutine planeStats_average_all


subroutine planeStats_average_finalize
  use planeStats_average
  implicit none

  integer :: iplane, ii, jj, isc, j, ierr
  real(WP) :: tmp1, tmp2
  real(WP), dimension(1:ny) :: tmpv1
  real(WP), dimension(1:3,1:3) :: evec_S, evec_R
  real(WP), dimension(1:3)     :: eval
  real(WP), dimension(:,:), pointer :: o_Pm

  allocate(o_Pm(1:nplanes,1:ny))

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

  ! Complete the density weighting
  do iplane=pnmin_,pnmax_
     do j=1,ny
        if (RHOm(j,iplane).gt.0.0_WP) then
           y_condm(j,iplane) = y_condm(j,iplane)/RHOm(j,iplane)
        end if
     end do
  end do
  do iplane=pnmin_,pnmax_
     do j=1,ny
        if (RHOm(j,iplane).gt.0.0_WP) then
           NU  (j,iplane) = NU  (j,iplane)/RHOm  (j,iplane)
        end if
     end do
  end do
  do iplane=pnmin_,pnmax_
     do j=1,nbins_cond
        if (RHOm_c(j,iplane).gt.0.0_WP) then
           NU_c(j,iplane) = NU_c(j,iplane)/RHOm_c(j,iplane)
        end if
     end do
  end do
  do iplane=pnmin_,pnmax_
     do j=1,ny
        if (RHOm(j,iplane).gt.0.0_WP) then
           VISCm_fc(j,iplane) = VISCm_fc(j,iplane)/RHOm(j,iplane)
        end if
     end do
  end do

  ! Compute averaged velocity gradients
  !  -- Um is still mean(RHO*U)
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           do j=1,ny
              if (RHOm(j,iplane).gt.0.0_WP) then
                 dUdxm(j,iplane,ii,jj) = dUdx_tmpv(j,iplane,ii,jj)/RHOm(j,iplane) &
                      - Um(j,iplane,ii)*dRHOdxm(j,iplane,jj)/RHOm(j,iplane)**2
              end if
           end do
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           do j=1,nbins_cond
              if (RHOm_c(j,iplane).gt.0.0_WP) then
                 dUdxm_c(j,iplane,ii,jj) = dUdx_tmpc(j,iplane,ii,jj)/RHOm_c(j,iplane) &
                      - Um_c(j,iplane,ii)*dRHOdxm_c(j,iplane,jj)/RHOm_c(j,iplane)**2
              end if
           end do
        end do
     end do
  end do

  ! Compute averaged scalar gradients
  !  -- SCm is still mean(RHO*SC)
  if (use_dSC) then
     do jj=1,3
        do isc=1,nscalar
           do iplane=pnmin_,pnmax_
              do j=1,ny
                 if (RHOm(j,iplane).gt.0.0_WP) then
                    dSCdxm(j,iplane,isc,jj) = dSCdx_tmpv(j,iplane,isc,jj)/RHOm(j,iplane) &
                         - SCm(j,iplane,isc)*dRHOdxm(j,iplane,jj)/RHOm(j,iplane)**2
                 end if
              end do
           end do
        end do
     end do
     do jj=1,3
        do isc=1,nscalar
           do iplane=pnmin_,pnmax_
              do j=1,nbins_cond
                 if (RHOm_c(j,iplane).gt.0.0_WP) then
                    dSCdxm_c(j,iplane,isc,jj) = dSCdx_tmpc(j,iplane,isc,jj)/RHOm_c(j,iplane) &
                         - SCm_c(j,iplane,isc)*dRHOdxm_c(j,iplane,jj)/RHOm_c(j,iplane)**2
                 end if
              end do
           end do
        end do
     end do
  end if

  ! Mean velocity and scalars
  do ii=1,3
     do iplane=pnmin_,pnmax_
        do j=1,ny
           if (RHOm(j,iplane).gt.0.0_WP) then
              Um  (j,iplane,ii) = Um  (j,iplane,ii)/RHOm  (j,iplane)
           end if
        end do
     end do
     do iplane=pnmin_,pnmax_
        do j=1,nbins_cond
           if (RHOm_c(j,iplane).gt.0.0_WP) then
              Um_c(j,iplane,ii) = Um_c(j,iplane,ii)/RHOm_c(j,iplane)
           end if
        end do
     end do
  end do
  do isc=1,nscalar
     do iplane=pnmin_,pnmax_
        do j=1,ny
           if (RHOm(j,iplane).gt.0.0_WP) then
              SCm  (j,iplane,isc) = SCm  (j,iplane,isc)/RHOm  (j,iplane)
           end if
        end do
     end do
  end do
  do isc=1,nscalar
     do iplane=pnmin_,pnmax_
        do j=1,nbins_cond
           if (RHOm_c(j,iplane).gt.0.0_WP) then  
              SCm_c(j,iplane,isc) = SCm_c(j,iplane,isc)/RHOm_c(j,iplane)
           end if
        end do
     end do
  end do

  ! Diffusivity from NGA is rho*D
  if (use_qSC) then
     ! Diffusivities
     do isc=1,nscalar
        do iplane=pnmin_,pnmax_
           do j=1,ny
              if (RHOm(j,iplane).gt.0.0_WP) then
                 DIFFm(j,iplane,isc) = DIFFm(j,iplane,isc)/RHOm(j,iplane)
              end if
           end do
        end do
     end do
  end if
  do isc=1,nscalar
     do iplane=pnmin_,pnmax_
        do j=1,ny
           if (RHOm(j,iplane).gt.0.0_WP) then
              DIFFm_fc(j,iplane,isc) = DIFFm_fc(j,iplane,isc)/RHOm(j,iplane)
           end if
        end do
     end do
  end do

  ! Source terms??
  ! test for evaluation from mean scalar
  do iplane=pnmin_,pnmax_
     do j=1,ny
        call planeStats_finitechem_source(SCm(j,iplane,1:N_tot+1),src_SCm_test(j,iplane,1:N_tot+1))
     end do
  end do


  ! ALIGNMENT OF SCALAR FLUXES & GRADIENTS
  ! Compute scalar fluxes
  do jj=1,3
     do isc=1,nscalar
        do iplane=pnmin_,pnmax_
           SC_flux  (:,iplane,isc,jj) = rhoUiSCm  (:,iplane,isc,jj)/RHOm  (:,iplane) - Um  (:,iplane,jj)*SCm  (:,iplane,isc)
           SC_flux_c(:,iplane,isc,jj) = rhoUiSCm_c(:,iplane,isc,jj)/RHOm_c(:,iplane) - Um_c(:,iplane,jj)*SCm_c(:,iplane,isc)
        end do
     end do
  end do

!!$  ! Compute angle of alignment with scalar gradients (align_SCF_SCG(j,iplane,isc))
!!$  align_SCF_SCG = 0.0_WP
!!$  if (use_dSC) then
!!$     do isc=1,nscalar
!!$        do iplane=pnmin_,pnmax_
!!$           do j=1,ny
!!$              tmp1 = 0.0_WP
!!$              tmp2 = 0.0_WP
!!$              do jj=1,3
!!$                 align_SCF_SCG(j,iplane,isc) = align_SCF_SCG(j,iplane,isc) &
!!$                      + SC_flux(j,iplane,isc,jj)*dSCdxm(j,iplane,isc,jj)
!!$                 tmp1 = tmp1 + SC_flux(j,iplane,isc,jj)**2
!!$                 tmp2 = tmp2 + dSCdxm (j,iplane,isc,jj)**2
!!$              end do
!!$              tmp1 = sqrt(tmp1)
!!$              tmp2 = sqrt(tmp2)
!!$              align_SCF_SCG(j,iplane,isc) = align_SCF_SCG(j,iplane,isc)/(tmp1*tmp2)
!!$           end do
!!$        end do
!!$     end do
!!$  end if


  ! EIGENVALUES OF REYNOLDS STRESS & STRAIN-RATE TENSORS
  ! Compute full & anisotropic Reynolds stress tensors
  trace_R = 0.0_WP
  trace_R_c = 0.0_WP
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           tau_ij  (:,iplane,ii,jj) = rhoUiUjm  (:,iplane,ii,jj)/RHOm  (:,iplane) - Um  (:,iplane,ii)*Um  (:,iplane,jj)
           tau_ij_c(:,iplane,ii,jj) = rhoUiUjm_c(:,iplane,ii,jj)/RHOm_c(:,iplane) - Um_c(:,iplane,ii)*Um_c(:,iplane,jj)
        end do
     end do
     do iplane=pnmin_,pnmax_
        ! Trace of Reynolds stress tensor is 2*TKE
        trace_R  (:,iplane) = trace_R  (:,iplane) + tau_ij  (:,iplane,jj,jj)
        trace_R_c(:,iplane) = trace_R_c(:,iplane) + tau_ij_c(:,iplane,jj,jj)
     end do
  end do
  ! Compute anisotropic part
  tau_ij_d  = tau_ij
  tau_ij_dc = tau_ij_c
  do jj=1,3
     do iplane=pnmin_,pnmax_
        tau_ij_d (:,iplane,jj,jj) = tau_ij  (:,iplane,jj,jj) - 1.0_WP/3.0_WP*trace_R  (:,iplane)
        tau_ij_dc(:,iplane,jj,jj) = tau_ij_c(:,iplane,jj,jj) - 1.0_WP/3.0_WP*trace_R_c(:,iplane)
     end do
  end do
  ! Compute rms of anisotropic part
  do iplane=pnmin_,pnmax_
     tmp2 = 0.0_WP
     do j=1,ny
        tmp1 = 0.0_WP
        do jj=1,3
           do ii=1,3
              tmp1 = tmp1 + tau_ij_d(j,iplane,ii,jj)**2/9.0_WP
           end do
        end do
        tmp2 = tmp2 + tmp1
        rms_R_local(j,iplane) = sqrt(tmp1)
     end do
     tmp2 = sqrt(tmp2/real(ny,WP))
     rms_tau(iplane) = tmp2
  end do

  ! Compute full & anisotropic strain-rate tensors
  trace_S = 0.0_WP
  trace_S_c = 0.0_WP
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           S_ij  (:,iplane,ii,jj) = 0.5_WP*(dUdxm  (:,iplane,ii,jj) + dUdxm  (:,iplane,jj,ii))
           S_ij_c(:,iplane,ii,jj) = 0.5_WP*(dUdxm_c(:,iplane,ii,jj) + dUdxm_c(:,iplane,jj,ii))
        end do
     end do
     do iplane=pnmin_,pnmax_
        trace_S  (:,iplane) = trace_S  (:,iplane) + S_ij  (:,iplane,jj,jj)
        trace_S_c(:,iplane) = trace_S_c(:,iplane) + S_ij_c(:,iplane,jj,jj)
     end do
  end do
  ! Compute anisotropic part
  S_ij_d = S_ij
  do jj=1,3
     do iplane=pnmin_,pnmax_
        S_ij_d (:,iplane,jj,jj) = S_ij  (:,iplane,jj,jj) - 1.0_WP/3.0_WP*trace_S  (:,iplane)
        S_ij_dc(:,iplane,jj,jj) = S_ij_c(:,iplane,jj,jj) - 1.0_WP/3.0_WP*trace_S_c(:,iplane)
     end do
  end do
  ! Compute rms of anisotropic strain-rate tensor
  do iplane=pnmin_,pnmax_
     tmp2 = 0.0_WP
     do j=1,ny
        tmp1 = 0.0_WP
        do jj=1,3
           do ii=1,3
              tmp1 = tmp1 + S_ij_d(j,iplane,ii,jj)**2/9.0_WP
           end do
        end do
        tmp2 = tmp2 + tmp1
        rms_S_local(j,iplane) = sqrt(tmp1)
     end do
     tmp2 = sqrt(tmp2/real(ny,WP))
     rms_Sij(iplane) = tmp2
  end do

  ! Compute dissipation rate
  diss_rate = 0.0_WP
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           diss_rate(:,iplane) = diss_rate(:,iplane) + taum(:,iplane,ii,jj)*dUdxm(:,iplane,ii,jj)
        end do
     end do
  end do
  do iplane=pnmin_,pnmax_
     do j=1,ny
        if (RHOm(j,iplane).gt.0.0_WP) then
           diss_rate(j,iplane) = ( tau_dUdxm(j,iplane) - diss_rate(j,iplane) )/RHOm(j,iplane)
        end if
     end do
  end do

!!$  ! Compute eigenvalues and angles of alignment between corresponding eigenvectors 
!!$  !                           (align_RST_SRT(j,iplane,ii), ii=eigenvalue index)
!!$  align_RST_SRT = 0.0_WP
!!$  do iplane=pnmin_,pnmax_
!!$     do j=1,ny
!!$        ! Eigenvalues of full strain rate tensor
!!$        ! evec : columns are eigenvectors corresponding to eigenvalues in eval
!!$        do jj=1,3
!!$           do ii=1,3
!!$              evec_S(ii,jj) = S_ij(j,iplane,ii,jj)
!!$           end do
!!$        end do
!!$        call DSYEV('V', 'U', 3, evec_S, 3, eval, tmpv1, ny, ierr)
!!$        if (ierr.eq.0) then
!!$           do ii=1,3
!!$              eval_S(j,iplane,ii) = eval(ii)
!!$           end do
!!$        else
!!$           print *, 'Something went wrong with DSYEV (S_ij) on ', iplane, ierr
!!$        end if
!!$        
!!$        ! Eigenvalues of full Reynolds stress tensor
!!$        do jj=1,3
!!$           do ii=1,3
!!$              evec_R(ii,jj) = tau_ij(j,iplane,ii,jj)
!!$           end do
!!$        end do
!!$        call DSYEV('V', 'U', 3, evec_R, 3, eval, tmpv1, ny, ierr)
!!$        if (ierr.eq.0) then
!!$           do ii=1,3
!!$              eval_R(j,iplane,ii) = eval(ii)
!!$           end do
!!$        else
!!$           print *, 'Something went wrong with DSYEV (tau_ij) on ', iplane, ierr
!!$        end if
!!$
!!$        ! Eigenvalues of anisotropic strain rate tensor
!!$        do jj=1,3
!!$           do ii=1,3
!!$              evec_S(ii,jj) = S_ij_d(j,iplane,ii,jj)
!!$           end do
!!$        end do
!!$        call DSYEV('V', 'U', 3, evec_S, 3, eval, tmpv1, ny, ierr)
!!$        if (ierr.eq.0) then
!!$           do ii=1,3
!!$              eval_S_d(j,iplane,ii) = eval(ii)!*trace_R(j,iplane)/diss_rate(j,iplane)  !/trace_S(j,iplane) !rms_Sij(iplane)
!!$              ! Can't normalize by trace (dilatation) -- zero except within the flame
!!$           end do
!!$        else
!!$           print *, 'Something went wrong with DSYEV (S_ij_d) on ', iplane, ierr
!!$        end if
!!$        
!!$        ! Eigenvalues of anisotropic Reynolds stress tensor
!!$        do jj=1,3
!!$           do ii=1,3
!!$              evec_R(ii,jj) = tau_ij_d(j,iplane,ii,jj)
!!$           end do
!!$        end do
!!$        call DSYEV('V', 'U', 3, evec_R, 3, eval, tmpv1, ny, ierr)
!!$        if (ierr.eq.0) then
!!$           do ii=1,3
!!$              eval_R_d(j,iplane,ii) = eval(ii)!/ trace_R(j,iplane) !rms_tau(iplane)
!!$           end do
!!$        else
!!$           print *, 'Something went wrong with DSYEV (tau_ij_d) on ', iplane, ierr
!!$        end if
!!$        !if (iplane.eq.1.and.j.eq.1) print *, eval_R_d(j,iplane,:)
!!$
!!$
!!$
!!$        ! Alignment angle (ii=eigenvalue index)
!!$        ! THIS PROBABLY DOESN'T MAKE MUCH SENSE....
!!$        do ii=1,3
!!$           tmp1 = 0.0_WP
!!$           tmp2 = 0.0_WP
!!$           do jj=1,3                              !! NOTE NEGATIVE SIGN due to definition of tau_ij
!!$              align_RST_SRT(j,iplane,ii) = align_RST_SRT(j,iplane,ii) - evec_S(jj,ii)*evec_R(jj,ii)
!!$              tmp1 = tmp1 + evec_S(jj,ii)**2
!!$              tmp2 = tmp2 + evec_R(jj,ii)**2
!!$           end do
!!$           tmp1 = sqrt(tmp1)
!!$           tmp2 = sqrt(tmp2)
!!$           align_RST_SRT(j,iplane,ii) = align_RST_SRT(j,iplane,ii)/(tmp1*tmp2)
!!$        end do
!!$     end do
!!$  end do

  ! Compute density weighting before output
  do jj=1,3
     do isc=1,nscalar
        do iplane=pnmin_,pnmax_
           rhoUiSCm  (:,iplane,isc,jj) = rhoUiSCm  (:,iplane,isc,jj)/RHOm  (:,iplane)
           rhoUiSCm_c(:,iplane,isc,jj) = rhoUiSCm_c(:,iplane,isc,jj)/RHOm_c(:,iplane)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           rhoUiUjm  (:,iplane,ii,jj) = rhoUiUjm  (:,iplane,ii,jj)/RHOm  (:,iplane)
           rhoUiUjm_c(:,iplane,ii,jj) = rhoUiUjm_c(:,iplane,ii,jj)/RHOm_c(:,iplane)
        end do
     end do
  end do
  do isc=1,nscalar
     do ii=1,3
        do iplane=pnmin_,pnmax_
           dUjSCm(:,iplane,ii,isc) = dUjSCm(:,iplane,ii,isc)/RHOm(:,iplane)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do iplane=pnmin_,pnmax_
           dUiUjm(:,iplane,ii,jj) = dUiUjm(:,iplane,ii,jj)/RHOm(:,iplane)
        end do
     end do
  end do

  ! Save terms in mean equations - for verification
  ! Velocity
  do ii=1,3
     do iplane=pnmin_,pnmax_
        ! Unsteady term (0)

        ! Convection term (1)
        do jj=1,3
           cons_mom_m_phys(:,iplane,1,ii) = cons_mom_m_phys(:,iplane,1,ii) &
                - Um(:,iplane,ii) * Um(:,iplane,jj) * dRHOdxm(:,iplane,jj) &
                - Um(:,iplane,ii) * RHOm(:,iplane) * dUdxm(:,iplane,jj,jj) &
                - Um(:,iplane,jj) * RHOm(:,iplane) * dUdxm(:,iplane,ii,jj)
        end do

        ! Divergence of Reynolds stress (2)
        !   rhoUiUjm is now tilde(UiUj)
        do jj=1,3
           cons_mom_m_phys(:,iplane,2,ii) = cons_mom_m_phys(:,iplane,2,ii) - &
                ( dRHOdxm(:,iplane,jj) * ( rhoUiUjm(:,iplane,ii,jj) - Um(:,iplane,ii)*Um(:,iplane,jj) ) &
                + RHOm(:,iplane) * dUiUjm(:,iplane,ii,jj) &
                - RHOm(:,iplane) * Um(:,iplane,ii) * dUdxm(:,iplane,jj,jj) &
                - RHOm(:,iplane) * Um(:,iplane,jj) * dUdxm(:,iplane,ii,jj) )
        end do

        ! Pressure term (3)
        cons_mom_m_phys(:,iplane,3,ii) = -dPdxm(:,iplane,ii)

        ! Stress tensor gradient (4)
        do jj=1,3
           cons_mom_m_phys(:,iplane,4,ii) = cons_mom_m_phys(:,iplane,4,ii) + dTAUdxm(:,iplane,ii,jj)
        end do
     end do
  end do

  ! Scalars
  print *, ntime_curr, use_dSC, use_qSC
  if (use_qSC) then
     do isc=1,nscalar-1
        do iplane=pnmin_,pnmax_
           ! Convection term (1)
           do jj=1,3
              cons_SC_m_phys(:,iplane,1,isc) = cons_SC_m_phys(:,iplane,1,isc) - &
                   ( Um(:,iplane,jj) * SCm(:,iplane,isc) * dRHOdxm(:,iplane,jj) &
                   + Um(:,iplane,jj) * RHOm(:,iplane) * dSCdxm(:,iplane,isc,jj) &
                   + SCm(:,iplane,isc) * RHOm(:,iplane) * dUdxm(:,iplane,jj,jj) )
           end do

           ! Fluctuating term (2)
           do jj=1,3
              cons_SC_m_phys(:,iplane,2,isc) = cons_SC_m_phys(:,iplane,2,isc) - &
                   ( dRHOdxm(:,iplane,jj) * ( rhoUiSCm(:,iplane,isc,jj) - Um(:,iplane,jj)*SCm(:,iplane,isc) ) &
                   + RHOm(:,iplane) * dUjSCm(:,iplane,jj,isc) &
                   - RHOm(:,iplane) * Um(:,iplane,jj) * dSCdxm(:,iplane,isc,jj) &
                   - RHOm(:,iplane) * SCm(:,iplane,isc) * dUdxm(:,iplane,jj,jj) )
           end do

           ! Diffusion term (3)
           cons_SC_m_phys(:,iplane,3,isc) = diff_flux_divm(:,iplane,isc)

           ! Source term (4)
           cons_SC_m_phys(:,iplane,4,isc) = src_SCm_fc(:,iplane,isc)
        end do
     end do
  end if
  
  return
end subroutine planeStats_average_finalize


subroutine planeStats_average_output
  use planeStats_average
  implicit none

  integer :: iunit, iplane, ii, jj, isc, j, ierr
  
  ! Arrays to collect
  real(WP), dimension(:,:),   pointer :: o_y_condm
!!$  real(WP), dimension(:,:,:), pointer :: o_align_SCF_SCG
!!$  real(WP), dimension(:,:,:), pointer :: o_align_RST_SRT
!!$  real(WP), dimension(:,:,:), pointer :: o_eval_S_d, o_eval_R_d
!!$  real(WP), dimension(:,:,:), pointer :: o_eval_S,   o_eval_R
  real(WP), dimension(:,:),   pointer :: o_trace_S, o_trace_R, o_diss_rate, o_rms_S_local, o_rms_R_local
  real(WP), dimension(:,:),   pointer :: o_trace_S_c, o_trace_R_c

  real(WP), dimension(:,:,:), pointer   :: o_U_phys, o_U_cond
  real(WP), dimension(:,:,:), pointer   :: o_SC_phys, o_SC_cond
  real(WP), dimension(:,:,:,:), pointer :: o_SC_flux, o_dSCdxm
  real(WP), dimension(:,:,:,:), pointer :: o_SC_flux_c, o_dSCdxm_c
  real(WP), dimension(:,:,:,:), pointer :: o_tau_ij,   o_S_ij
  real(WP), dimension(:,:,:,:), pointer :: o_tau_ij_c, o_S_ij_c
  real(WP), dimension(:,:,:,:), pointer :: o_tau_ij_d, o_S_ij_d
  real(WP), dimension(:,:,:,:), pointer :: o_tau_ij_dc, o_S_ij_dc

  real(WP), dimension(:,:,:,:), pointer :: o_rhoUiSCm, o_rhoUiSCm_c
  real(WP), dimension(:,:,:,:), pointer :: o_rhoUiUjm, o_rhoUiUjm_c

  real(WP), dimension(:,:,:), pointer :: o_DIFFm, o_src_SCm,o_src_SCm_test

  real(WP), dimension(:,:,:,:), pointer :: o_cons_mom_m_phys, o_cons_SC_m_phys

  ! Allocate them
  allocate(o_y_condm      (1:nplanes,1:ny))
!!$  allocate(o_align_SCF_SCG(1:nplanes,1:ny,1:nscalar))
!!$  allocate(o_align_RST_SRT(1:nplanes,1:ny,1:3))
!!$  allocate(o_eval_S       (1:nplanes,1:ny,1:3))
!!$  allocate(o_eval_R       (1:nplanes,1:ny,1:3))
!!$  allocate(o_eval_S_d     (1:nplanes,1:ny,1:3))
!!$  allocate(o_eval_R_d     (1:nplanes,1:ny,1:3))
  allocate(o_trace_S  (1:nplanes,1:ny))
  allocate(o_trace_R  (1:nplanes,1:ny))
  allocate(o_trace_R_c(1:nplanes,1:nbins_cond))
  allocate(o_trace_S_c(1:nplanes,1:nbins_cond))
  allocate(o_diss_rate(1:nplanes,1:ny))
  allocate(o_rms_S_local(1:nplanes,1:ny))
  allocate(o_rms_R_local(1:nplanes,1:ny))

  allocate(o_SC_flux (1:nplanes,1:ny,1:nscalar,1:3))
  allocate(o_dSCdxm  (1:nplanes,1:ny,1:nscalar,1:3))
  allocate(o_rhoUiSCm(1:nplanes,1:ny,1:nscalar,1:3))

  allocate(o_tau_ij  (1:nplanes,1:ny,1:3,1:3))
  allocate(o_tau_ij_d(1:nplanes,1:ny,1:3,1:3))
  allocate(o_rhoUiUjm(1:nplanes,1:ny,1:3,1:3))
  allocate(o_tau_ij_c  (1:nplanes,1:nbins_cond,1:3,1:3))
  allocate(o_tau_ij_dc (1:nplanes,1:nbins_cond,1:3,1:3))
  allocate(o_rhoUiUjm_c(1:nplanes,1:nbins_cond,1:3,1:3))
  allocate(o_S_ij    (1:nplanes,1:ny,1:3,1:3))
  allocate(o_S_ij_d  (1:nplanes,1:ny,1:3,1:3))
  allocate(o_S_ij_c  (1:nplanes,1:nbins_cond,1:3,1:3))
  allocate(o_S_ij_dc (1:nplanes,1:nbins_cond,1:3,1:3))

  allocate(o_SC_flux_c (1:nplanes,1:nbins_cond,1:nscalar,1:3))
  allocate(o_dSCdxm_c  (1:nplanes,1:nbins_cond,1:nscalar,1:3))
  allocate(o_rhoUiSCm_c(1:nplanes,1:nbins_cond,1:nscalar,1:3))

  ! Velocity
  allocate(o_U_phys(1:nplanes,1:ny,1:3))
  allocate(o_U_cond(1:nplanes,1:nbins_cond,1:3))

  ! Scalars
  allocate(o_SC_phys(1:nplanes,1:ny,1:nscalar))
  allocate(o_SC_cond(1:nplanes,1:nbins_cond,1:nscalar))
  allocate(o_DIFFm  (1:nplanes,1:ny,1:nscalar-2))
  allocate(o_src_SCm(1:nplanes,1:ny,1:nscalar-1))
  allocate(o_src_SCm_test(1:nplanes,1:ny,1:nscalar-1))

  ! For verification
  allocate(o_cons_mom_m_phys(1:nplanes,1:ny,1:4,1:3))
  allocate(o_cons_SC_m_phys (1:nplanes,1:ny,1:4,1:nscalar-1))

  ! Collect from the processes
  do j=1,ny
     call MPI_GATHER(y_condm(j,:), nplanes_, MPI_REAL_WP, o_y_condm(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
!!$  do isc=1,nscalar
!!$     do j=1,ny
!!$        call MPI_GATHER(align_SCF_SCG(j,:,isc), nplanes_, MPI_REAL_WP, o_align_SCF_SCG(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
!!$     end do
!!$  end do
!!$  do ii=1,3
!!$     do j=1,ny
!!$        call MPI_GATHER(align_RST_SRT(j,:,ii), nplanes_, MPI_REAL_WP, o_align_RST_SRT(:,j,ii), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
!!$     end do
!!$  end do
!!$  do ii=1,3
!!$     do j=1,ny
!!$        call MPI_GATHER(eval_S(j,:,ii), nplanes_, MPI_REAL_WP, o_eval_S(:,j,ii), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
!!$     end do
!!$  end do
!!$  do ii=1,3
!!$     do j=1,ny
!!$        call MPI_GATHER(eval_R(j,:,ii), nplanes_, MPI_REAL_WP, o_eval_R(:,j,ii), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
!!$     end do
!!$  end do
!!$  do ii=1,3
!!$     do j=1,ny
!!$        call MPI_GATHER(eval_S_d(j,:,ii), nplanes_, MPI_REAL_WP, o_eval_S_d(:,j,ii), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
!!$     end do
!!$  end do
!!$  do ii=1,3
!!$     do j=1,ny
!!$        call MPI_GATHER(eval_R_d(j,:,ii), nplanes_, MPI_REAL_WP, o_eval_R_d(:,j,ii), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
!!$     end do
!!$  end do
  do j=1,ny
     call MPI_GATHER(trace_S(j,:), nplanes_, MPI_REAL_WP, o_trace_S(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(trace_S_c(j,:), nplanes_, MPI_REAL_WP, o_trace_S_c(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(trace_R(j,:), nplanes_, MPI_REAL_WP, o_trace_R(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,nbins_cond
     call MPI_GATHER(trace_R_c(j,:), nplanes_, MPI_REAL_WP, o_trace_R_c(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(diss_rate(j,:), nplanes_, MPI_REAL_WP, o_diss_rate(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(rms_S_local(j,:), nplanes_, MPI_REAL_WP, o_rms_S_local(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do j=1,ny
     call MPI_GATHER(rms_R_local(j,:), nplanes_, MPI_REAL_WP, o_rms_R_local(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do

  do jj=1,3
     do isc=1,nscalar
        do j=1,ny
           call MPI_GATHER(SC_flux(j,:,isc,jj), nplanes_, MPI_REAL_WP, o_SC_flux(:,j,isc,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do isc=1,nscalar
        do j=1,ny
           call MPI_GATHER(dSCdxm(j,:,isc,jj), nplanes_, MPI_REAL_WP, o_dSCdxm(:,j,isc,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=1,ny
           call MPI_GATHER(tau_ij(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_tau_ij(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=1,ny
           call MPI_GATHER(tau_ij_d(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_tau_ij_d(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=1,nbins_cond
           call MPI_GATHER(tau_ij_c(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_tau_ij_c(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=1,nbins_cond
           call MPI_GATHER(tau_ij_dc(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_tau_ij_dc(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=1,ny
           call MPI_GATHER(S_ij(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_S_ij(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=1,ny
           call MPI_GATHER(S_ij_d(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_S_ij_d(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=1,nbins_cond
           call MPI_GATHER(S_ij_c(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_S_ij_c(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=1,nbins_cond
           call MPI_GATHER(S_ij_dc(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_S_ij_dc(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
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
           call MPI_GATHER(dSCdxm_c(j,:,isc,jj), nplanes_, MPI_REAL_WP, o_dSCdxm_c(:,j,isc,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do

  do jj=1,3
     do isc=1,nscalar
        do j=1,ny
           call MPI_GATHER(rhoUiSCm(j,:,isc,jj), nplanes_, MPI_REAL_WP, o_rhoUiSCm(:,j,isc,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do isc=1,nscalar
        do j=1,nbins_cond
           call MPI_GATHER(rhoUiSCm_c(j,:,isc,jj), nplanes_, MPI_REAL_WP, o_rhoUiSCm_c(:,j,isc,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do

  do jj=1,3
     do ii=1,3
        do j=1,ny
           call MPI_GATHER(rhoUiUjm(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_rhoUiUjm(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=1,nbins_cond
           call MPI_GATHER(rhoUiUjm_c(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_rhoUiUjm_c(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  ! Velocity
  do jj=1,3
     do j=1,ny
        call MPI_GATHER(Um(j,:,jj), nplanes_, MPI_REAL_WP, o_U_phys(:,j,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do
  do jj=1,3
     do j=1,nbins_cond
        call MPI_GATHER(Um_c(j,:,jj), nplanes_, MPI_REAL_WP, o_U_cond(:,j,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do
  ! Scalars
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
  do isc=1,nscalar-2
     do j=1,ny
        call MPI_GATHER(DIFFm(j,:,isc), nplanes_, MPI_REAL_WP, o_DIFFm(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do
  do isc=1,nscalar-1
     do j=1,ny
        call MPI_GATHER(src_SCm_fc(j,:,isc), nplanes_, MPI_REAL_WP, o_src_SCm(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do
  do isc=1,nscalar-1
     do j=1,ny
        call MPI_GATHER(src_SCm_test(j,:,isc), nplanes_, MPI_REAL_WP, o_src_SCm_test(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
     end do
  end do



  ! For verification
  do ii=1,3
     do jj=1,4
        do j=1,ny
           call MPI_GATHER(cons_mom_m_phys(j,:,jj,ii), nplanes_, MPI_REAL_WP, o_cons_mom_m_phys(:,j,jj,ii), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  if (use_qSC) then
     do isc=1,nscalar-1
        do jj=1,4
           do j=1,ny
              call MPI_GATHER(cons_SC_m_phys(j,:,jj,isc), nplanes_, MPI_REAL_WP, o_cons_SC_m_phys(:,j,jj,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
           end do
        end do
     end do
  end if

  ! Only the root process writes output
  if (irank.ne.iroot) return
  
  ! Output the files

  ! Scalar fluxes & gradients
  open(unit=iunit, file=trim(output_name)//'_SCF_SCG', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j) ! y-coordinate
        end if
        write(iunit,'(ES22.13)',advance='no') o_y_condm(iplane,j) ! Progress variable
        do isc=1,nscalar
           do jj=1,3
              ! Scalar flux
              write(iunit,'(ES22.13)',advance='no') o_SC_flux(iplane,j,isc,jj)
           end do
           do jj=1,3
              ! Mean scalar gradient
              write(iunit,'(ES22.13)',advance='no') o_dSCdxm(iplane,j,isc,jj)
           end do
           do jj=1,3
              ! Raw second moment
              if (iplane.eq.nplanes .and. isc.eq.nscalar .and. jj.eq.3) then
                 write(iunit,'(ES22.13)') o_rhoUiSCm(iplane,j,isc,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') o_rhoUiSCm(iplane,j,isc,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)

  ! Scalar fluxes & gradients - Conditional space
  open(unit=iunit, file=trim(output_name)//'_SCF_SCG_cond', action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        do isc=1,nscalar
           do jj=1,3
              ! Scalar flux
              write(iunit,'(ES22.13)',advance='no') o_SC_flux_c(iplane,j,isc,jj)
           end do
           do jj=1,3
              ! Mean scalar gradient
              write(iunit,'(ES22.13)',advance='no') o_dSCdxm_c(iplane,j,isc,jj)
           end do
           do jj=1,3
              ! Raw second moment
              if (iplane.eq.nplanes .and. isc.eq.nscalar .and. jj.eq.3) then
                 write(iunit,'(ES22.13)') o_rhoUiSCm_c(iplane,j,isc,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') o_rhoUiSCm_c(iplane,j,isc,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)

!!$  ! Alignment of scalar fluxes & gradients
!!$  open(unit=iunit, file=trim(output_name)//'_align_SCF_SCG', action='write')
!!$  do j=1,ny
!!$     do iplane=1,nplanes
!!$        if (iplane.eq.1) then
!!$           write(iunit,'(ES22.13)',advance='no') y(j) ! y-coordinate
!!$        end if
!!$        write(iunit,'(ES22.13)',advance='no') o_y_condm(iplane,j) ! Progress variable
!!$        do isc=1,nscalar
!!$           if (iplane.eq.nplanes .and. isc.eq.nscalar) then
!!$              write(iunit,'(ES22.13)') o_align_SCF_SCG(iplane,j,isc)
!!$           else
!!$              write(iunit,'(ES22.13)',advance='no') o_align_SCF_SCG(iplane,j,isc)
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!!$  close(iunit)

!!$  ! Alignment of eigenvectors of Reynolds stress & strain-rate tensors
!!$  ! USE FOR ALIGNMENT ONLY
!!$  ! EIGENVALUES ARE FOR DEVIATORIC TENSORS
!!$  open(unit=iunit, file=trim(output_name)//'_align_RST_SRT', action='write')
!!$  do j=1,ny
!!$     do iplane=1,nplanes
!!$        if (iplane.eq.1) then
!!$           write(iunit,'(ES22.13)',advance='no') y(j) ! y-coordinate
!!$        end if
!!$        write(iunit,'(ES22.13)',advance='no') o_y_condm(iplane,j) ! Progress variable
!!$        do ii=1,3
!!$           if (iplane.eq.nplanes .and. ii.eq.3) then
!!$              write(iunit,'(3ES22.13)') &
!!$                   o_align_RST_SRT(iplane,j,ii), o_eval_S_d(iplane,j,ii), o_eval_R_d(iplane,j,ii)
!!$           else
!!$              write(iunit,'(3ES22.13)',advance='no') &
!!$                   o_align_RST_SRT(iplane,j,ii), o_eval_S_d(iplane,j,ii), o_eval_R_d(iplane,j,ii)
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!!$  close(iunit)
!!$
!!$  ! Eigenvalues, plus alignment of eigenvectors, of Reynolds stress & strain-rate tensors
!!$  ! NEW NORMALIZATION
!!$  ! EIGENVALUES ARE FOR DEVIATORIC TENSORS
!!$  open(unit=iunit, file=trim(output_name)//'_align_RST_SRT_norm', action='write')
!!$  do j=1,ny
!!$     do iplane=1,nplanes
!!$        if (iplane.eq.1) then
!!$           write(iunit,'(ES22.13)',advance='no') y(j) ! y-coordinate
!!$        end if
!!$        write(iunit,'(6ES22.13)',advance='no') &
!!$             o_y_condm(iplane,j),   o_trace_S(iplane,j),     o_trace_R(iplane,j), &
!!$             o_diss_rate(iplane,j), o_rms_S_local(iplane,j), o_rms_R_local(iplane,j)
!!$        do ii=1,3
!!$           if (iplane.eq.nplanes .and. ii.eq.3) then
!!$              write(iunit,'(3ES22.13)') &
!!$                   o_align_RST_SRT(iplane,j,ii), o_eval_S_d(iplane,j,ii), o_eval_R_d(iplane,j,ii)
!!$           else
!!$              write(iunit,'(3ES22.13)',advance='no') &
!!$                   o_align_RST_SRT(iplane,j,ii), o_eval_S_d(iplane,j,ii), o_eval_R_d(iplane,j,ii)
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!!$  close(iunit)

  ! All components of DEVIATORIC strain-rate and Reynolds stress tensors
  open(unit=iunit, file=trim(output_name)//'_RST_SRT', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j) ! y-coordinate
        end if
        write(iunit,'(6ES22.13)',advance='no') &
             o_y_condm(iplane,j),   o_trace_S(iplane,j),     o_trace_R(iplane,j), &
             o_diss_rate(iplane,j), o_rms_S_local(iplane,j), o_rms_R_local(iplane,j)
        do ii=1,3
           do jj=1,3
              write(iunit,'(ES22.13)',advance='no') o_tau_ij_d(iplane,j,ii,jj)
           end do
        end do
        do ii=1,3
           do jj=1,3
              if (iplane.eq.nplanes .and. ii.eq.3 .and. jj.eq.3) then
                 write(iunit,'(ES22.13)') o_S_ij_d(iplane,j,ii,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') o_S_ij_d(iplane,j,ii,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)



  ! All components of DEVIATORIC strain-rate and Reynolds stress tensors
  ! TRUE CONDITIONAL MEAN
  open(unit=iunit, file=trim(output_name)//'_RST_SRT_cond', action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j) ! y-coordinate
        end if
        write(iunit,'(2ES22.13)',advance='no') &
             o_trace_S_c(iplane,j),     o_trace_R_c(iplane,j)
        do ii=1,3
           do jj=1,3
              write(iunit,'(ES22.13)',advance='no') o_tau_ij_dc(iplane,j,ii,jj)
           end do
        end do
        do ii=1,3
           do jj=1,3
              if (iplane.eq.nplanes .and. ii.eq.3 .and. jj.eq.3) then
                 write(iunit,'(ES22.13)') o_S_ij_dc(iplane,j,ii,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') o_S_ij_dc(iplane,j,ii,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)



!!$  ! Eigenvalues of FULL Reynolds stress & strain-rate tensors
!!$  open(unit=iunit, file=trim(output_name)//'_eval_RST_SRT_full', action='write')
!!$  do j=1,ny
!!$     do iplane=1,nplanes
!!$        if (iplane.eq.1) then
!!$           write(iunit,'(ES22.13)',advance='no') y(j) ! y-coordinate
!!$        end if
!!$        write(iunit,'(ES22.13)',advance='no') o_y_condm(iplane,j)
!!$        do ii=1,3
!!$           if (iplane.eq.nplanes .and. ii.eq.3) then
!!$              write(iunit,'(2ES22.13)') &
!!$                   o_eval_S(iplane,j,ii), o_eval_R(iplane,j,ii)
!!$           else
!!$              write(iunit,'(2ES22.13)',advance='no') &
!!$                   o_eval_S(iplane,j,ii), o_eval_R(iplane,j,ii)
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!!$  close(iunit)



  ! All components of FULL strain-rate and Reynolds stress tensors
  open(unit=iunit, file=trim(output_name)//'_RST_SRT_full', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j) 
        end if
        write(iunit,'(ES22.13)',advance='no') o_y_condm(iplane,j)
        do ii=1,3
           do jj=1,3
              ! Reynolds stress
              write(iunit,'(ES22.13)',advance='no') o_tau_ij(iplane,j,ii,jj)
           end do
        end do
        do ii=1,3
           do jj=1,3
              ! Strain-rate tensor
              write(iunit,'(ES22.13)',advance='no') o_S_ij(iplane,j,ii,jj)
           end do
        end do
        do ii=1,3
           do jj=1,3
              ! Raw second moment
              if (iplane.eq.nplanes .and. ii.eq.3 .and. jj.eq.3) then
                 write(iunit,'(ES22.13)') o_rhoUiUjm(iplane,j,ii,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') o_rhoUiUjm(iplane,j,ii,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)


  ! TRUE CONDITIONAL MEAN
  open(unit=iunit, file=trim(output_name)//'_RST_SRT_full_cond', action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j) 
        end if
        do ii=1,3
           do jj=1,3
              ! Reynolds stress
              write(iunit,'(ES22.13)',advance='no') o_tau_ij_c(iplane,j,ii,jj)
            end do
        end do
        do ii=1,3
           do jj=1,3
              ! Strain-rate tensor
              write(iunit,'(ES22.13)',advance='no') o_S_ij_c(iplane,j,ii,jj)
           end do
        end do
        do ii=1,3
           do jj=1,3
              ! Raw second moment
              if (iplane.eq.nplanes .and. ii.eq.3 .and. jj.eq.3) then
                 write(iunit,'(ES22.13)') o_rhoUiUjm_c(iplane,j,ii,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') o_rhoUiUjm_c(iplane,j,ii,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)


  ! Velocity
  open(unit=iunit, file=trim(output_name)//'_U_phys', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        write(iunit,'(ES22.13)',advance='no') o_y_condm(iplane,j)
        do jj=1,3
           if (iplane.eq.nplanes.and.jj.eq.3) then
              write(iunit,'(ES22.13)') &
                   o_U_phys(iplane,j,jj)
           else
              write(iunit,'(ES22.13)',advance='no') &
                   o_U_phys(iplane,j,jj)
           end if
        end do
     end do
  end do
  close(iunit)

  open(unit=iunit, file=trim(output_name)//'_U_cond', action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        do jj=1,3
           if (iplane.eq.nplanes .and. jj.eq.3) then
              write(iunit,'(ES22.13)') &
                   o_U_cond(iplane,j,jj)
           else
              write(iunit,'(ES22.13)',advance='no') &
                   o_U_cond(iplane,j,jj)
           end if
        end do
     end do
  end do
  close(iunit)

  ! Scalars
  open(unit=iunit, file=trim(output_name)//'_SC_phys', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        write(iunit,'(ES22.13)',advance='no') o_y_condm(iplane,j)
        do isc=1,nscalar-2
           write(iunit,'(3ES22.13)',advance='no') &
                o_SC_phys(iplane,j,isc), o_src_SCm(iplane,j,isc), o_DIFFm(iplane,j,isc)
        end do
        isc = nscalar-1
        write(iunit,'(3ES22.13)',advance='no') &
             o_SC_phys(iplane,j,isc), o_src_SCm(iplane,j,isc), 0.0_WP
        isc = nscalar
        if (iplane.eq.nplanes) then
           write(iunit,'(3ES22.13)') &
                o_SC_phys(iplane,j,isc), 0.0_WP, 0.0_WP
        else
           write(iunit,'(3ES22.13)',advance='no') &
                o_SC_phys(iplane,j,isc), 0.0_WP, 0.0_WP
        end if
     end do
  end do
  close(iunit)

  open(unit=iunit, file=trim(output_name)//'_SC_cond', action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        do isc=1,nscalar
           if (iplane.eq.nplanes .and. isc.eq.nscalar) then
              write(iunit,'(ES22.13)') &
                   o_SC_cond(iplane,j,isc)!, o_src_SC_out_cond(iplane,j,isc)
           else
              write(iunit,'(ES22.13)',advance='no') &
                   o_SC_cond(iplane,j,isc)!, o_src_SC_out_cond(iplane,j,isc)
           end if
        end do
     end do
  end do
  close(iunit)

  ! Mean velocity equations
  open(unit=iunit, file=trim(output_name)//'_cons_mom_m_phys',   action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        write (iunit,'(ES22.13)',advance='no') o_y_condm(iplane,j)
        do ii=1,3
           do jj=1,4
              if (iplane.eq.nplanes.and.jj.eq.4.and.ii.eq.3) then
                 write(iunit,'(ES22.13)') o_cons_mom_m_phys(iplane,j,jj,ii)
              else
                 write(iunit,'(ES22.13)',advance='no') o_cons_mom_m_phys(iplane,j,jj,ii)
              end if
           end do
        end do
     end do
  end do
  close(iunit)
  
  if (use_qSC) then
     ! Mean scalar equations
     open(unit=iunit, file=trim(output_name)//'_cons_SC_m_phys',   action='write')
     do j=1,ny
        do iplane=1,nplanes
           if (iplane.eq.1) then
              write(iunit,'(ES22.13)',advance='no') y(j)
           end if
           write (iunit,'(ES22.13)',advance='no') o_y_condm(iplane,j)
           do isc=1,nscalar-1
              do jj=1,4
                 if (iplane.eq.nplanes.and.jj.eq.4.and.isc.eq.(nscalar-1)) then
                    write(iunit,'(ES22.13)') o_cons_SC_m_phys(iplane,j,jj,isc)
                 else
                    write(iunit,'(ES22.13)',advance='no') o_cons_SC_m_phys(iplane,j,jj,isc)
                 end if
              end do
           end do
        end do
     end do
     close(iunit)
  end if

  ! Scalars
  open(unit=iunit, file=trim(output_name)//'_src_SC_phys', action='write')
  do j=1,ny
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') y(j)
        end if
        write(iunit,'(ES22.13)',advance='no') o_y_condm(iplane,j)
        do isc=1,nscalar-1
           if (iplane.eq.nplanes.and.isc.eq.(nscalar-1)) then
              write(iunit,'(3ES22.13)') &
                   o_SC_phys(iplane,j,isc), o_src_SCm(iplane,j,isc), o_src_SCm_test(iplane,j,isc)
           else
              write(iunit,'(3ES22.13)',advance='no') &
                   o_SC_phys(iplane,j,isc), o_src_SCm(iplane,j,isc), o_src_SCm_test(iplane,j,isc)
           end if
        end do
     end do
  end do
  close(iunit)

  return
end subroutine planeStats_average_output
