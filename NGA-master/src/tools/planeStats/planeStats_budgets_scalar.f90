subroutine planeStats_budgets_scalar_compute
  use planeStats_budgets
  implicit none

  integer :: iplane, isc
  real(WP), dimension(1:ny) :: tmpv1
  real(WP), dimension(1:nbins_cond) :: tmpc1
  real(WP), dimension(1:nz,1:ny) :: tmpxy

  if (.not.use_qSC) return

  ! Just look at H2O for now
  isc_sc = isc_H2O   !isc_H2O

  do isc=1,nscalar-2
     if (isc.eq.isc_sc) then
        ! Scalar variance and unsteady term (RHO*SC''^2)
        ! Physical space
        do iplane=pnmin_,pnmax_
           ! Scalar variance
           tmpxy = RHO(:,:,iplane)*SC(:,:,iplane,isc)**2
           call condmean(iplane, tmpxy, SCV_phys(:,iplane), SCV_cond(:,iplane), icount_qSC)

           ! Save unsteady term
           if (ntime_qSC.gt.1) then
              tmpxy = (tmpxy - tmpv_save_SCV_phys(:,:,iplane))/(time(ntime_curr) - time(ntime_curr-1))
              call condmean(iplane, tmpxy, SCV_unst_phys(:,iplane), SCV_unst_cond(:,iplane), icount_qSC)
           end if

           ! Save for next time step
           tmpv_save_SCV_phys(:,:,iplane) = RHO(:,:,iplane)*SC(:,:,iplane,isc)**2
        end do

!!$        ! Conditional space
!!$        do iplane=pnmin_,pnmax_
!!$           ! Scalar variance
!!$           tmpxy = RHO(:,:,iplane)*SC_c(:,:,iplane,isc)**2
!!$           call cmean(iplane, tmpxy, tmpc1)
!!$           SCV_cond(:,iplane) = (real(ntime_qSC-1,WP)*SCV_cond(:,iplane) + tmpc1)/real(ntime_qSC,WP)
!!$
!!$           ! Save unsteady term
!!$           if (ntime_qSC.gt.1) then
!!$              tmpxy = (tmpxy - tmpv_save_SCV_cond(:,:,iplane))/(time(ntime_curr) - time(ntime_curr-1))
!!$              call cmean(iplane, tmpxy, tmpc1)
!!$              SCV_unst_cond(:,iplane) = (real(ntime_qSC-2,WP)*SCV_unst_cond(:,iplane) + tmpc1)/real(ntime_qSC-1,WP)
!!$           end if
!!$
!!$           ! Save for next time step
!!$           tmpv_save_SCV_cond(:,:,iplane) = RHO(:,:,iplane)*SC_c(:,:,iplane,isc)**2
!!$        end do

        ! Mean convection term
        do iplane=pnmin_,pnmax_
           call planeStats_budgets_scalar_conv(iplane,tmpv1)
           ! Physical space
           SCV_conv_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_conv_phys(:,iplane) + tmpv1)/real(ntime_qSC,WP)
!!$           ! Conditional space
!!$           SCV_conv_cond(:,iplane) = (real(ntime_qSC-1,WP)*SCV_conv_cond(:,iplane) + tmpc1)/real(ntime_qSC,WP)
        end do

        ! Divergence term 1 (triple product)
        do iplane=pnmin_,pnmax_
           call planeStats_budgets_scalar_div1(iplane,tmpv1)
           ! Physical space
           SCV_div1_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_div1_phys(:,iplane) + tmpv1)/real(ntime_qSC,WP)
!!$           ! Conditional space
!!$           SCV_div1_cond(:,iplane) = (real(ntime_qSC-1,WP)*SCV_div1_cond(:,iplane) + tmpc1)/real(ntime_qSC,WP)
        end do

        ! Divergence term 2 (diffusive flux)
        do iplane=pnmin_,pnmax_
           call planeStats_budgets_scalar_div2(iplane,tmpv1)
           ! Physical space
           SCV_div2_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_div2_phys(:,iplane) + tmpv1)/real(ntime_qSC,WP)
!!$           ! Conditional space
!!$           SCV_div2_cond(:,iplane) = (real(ntime_qSC-1,WP)*SCV_div2_cond(:,iplane) + tmpc1)/real(ntime_qSC,WP)
        end do

        ! Production term
        do iplane=pnmin_,pnmax_
           call planeStats_budgets_scalar_prod(iplane,tmpv1)
           ! Physical space
           SCV_prod_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_prod_phys(:,iplane) + tmpv1)/real(ntime_qSC,WP)
!!$           ! Conditional space
!!$           SCV_prod_cond(:,iplane) = (real(ntime_qSC-1,WP)*SCV_prod_cond(:,iplane) + tmpc1)/real(ntime_qSC,WP)
        end do

        ! Dissipation term -- USING DIFFUSIVE FLUX
        do iplane=pnmin_,pnmax_
           call planeStats_budgets_scalar_diss(iplane,tmpv1)
           ! Physical space
           SCV_diss_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_diss_phys(:,iplane) + tmpv1)/real(ntime_qSC,WP)
!!$           ! Conditional space
!!$           SCV_diss_cond(:,iplane) = (real(ntime_qSC-1,WP)*SCV_diss_cond(:,iplane) + tmpc1)/real(ntime_qSC,WP)
        end do

        ! Dissipation term -- USING DIFFUSIVITIES
        do iplane=pnmin_,pnmax_
           call planeStats_budgets_scalar_diss2(iplane,tmpv1)
           ! Physical space
           SCV_diss2_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_diss2_phys(:,iplane) + tmpv1)/real(ntime_qSC,WP)
!!$           ! Conditional space
!!$           SCV_diss2_cond(:,iplane) = (real(ntime_qSC-1,WP)*SCV_diss2_cond(:,iplane) + tmpc1)/real(ntime_qSC,WP)
        end do

        ! Chemical source term
        do iplane=pnmin_,pnmax_
           call planeStats_budgets_scalar_src(iplane,tmpv1)
           ! Physical space
           SCV_src_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_src_phys(:,iplane) + tmpv1)/real(ntime_qSC,WP)
!!$           ! Conditional space
!!$           SCV_src_cond(:,iplane) = (real(ntime_qSC-1,WP)*SCV_src_cond(:,iplane) + tmpc1)/real(ntime_qSC,WP)
        end do
     end if
  end do

  return
end subroutine planeStats_budgets_scalar_compute



! ========================================================== !
! Mean convection term                                       !
! ========================================================== !
subroutine planeStats_budgets_scalar_conv(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: j,k,jj
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c
  real(WP) :: SCV

!!$  ! Conditional space
!!$  do j=1,ny
!!$     do k=1,nz
!!$        ib  = bin_index(k,j,iplane)
!!$        SCV = SC_c(k,j,iplane,isc_sc)**2
!!$        ! Divergence of the Favre-averaged velocity field
!!$        out(k,j) = -SCV*(dUdxm_c(ib,iplane,1,1) + dUdxm_c(ib,iplane,2,2) + dUdxm_c(ib,iplane,3,3))
!!$        ! Advection of SCV with the Favre-averaged velocities
!!$        do jj=1,3
!!$           out(k,j) = out(k,j) - 2.0_WP*Um_c(ib,iplane,jj)*SC_c(k,j,iplane,isc_sc)*dSCdx_c(k,j,iplane,isc_sc,jj)
!!$        end do
!!$        out(k,j) = out(k,j)*RHO(k,j,iplane)
!!$        ! Density gradient term
!!$        do jj=1,3
!!$           out(k,j) = out(k,j) - SCV*Um_c(ib,iplane,jj)*dRHOdx(k,j,iplane,jj)
!!$        end do
!!$     end do
!!$  end do
!!$  call cmean(iplane, out, out_c)

  ! Physical space
  do j=1,ny
     do k=1,nz
        SCV = SC(k,j,iplane,isc_sc)**2
        ! Divergence of the Favre-averaged velocity field
        out(k,j) = -SCV*(dUdxm(j,iplane,1,1) + dUdxm(j,iplane,2,2) + dUdxm(j,iplane,3,3))
        ! Advection of SCV with the Favre-averaged velocities
        do jj=1,3
           out(k,j) = out(k,j) - 2.0_WP*Um(j,iplane,jj)*SC(k,j,iplane,isc_sc)*dSCdx(k,j,iplane,isc_sc,jj)
        end do
        out(k,j) = out(k,j)*RHO(k,j,iplane)
        ! Density gradient term
        do jj=1,3
           out(k,j) = out(k,j) - SCV*Um(j,iplane,jj)*dRHOdx(k,j,iplane,jj)
        end do
     end do
  end do
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, SCV_conv_cond(:,iplane), icount_qSC)
!!$  SCV_conv_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_conv_phys(:,iplane) + out_z)/real(ntime_qSC,WP)

  return
end subroutine planeStats_budgets_scalar_conv


! ========================================================== !
! Divergence term 1 (triple product)                         !
! ========================================================== !
subroutine planeStats_budgets_scalar_div1(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: j,k,jj
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c
  real(WP) :: SCV

!!$  ! Conditional space
!!$  do j=1,ny
!!$     do k=1,nz
!!$        SCV = SC_c(k,j,iplane,isc_sc)**2
!!$        out(k,j) = 0.0_WP
!!$        do jj=1,3
!!$           out(k,j) = out(k,j) - 2.0_WP*U_c(k,j,iplane,jj)*SC_c(k,j,iplane,isc_sc)*dSCdx_c(k,j,iplane,isc_sc,jj)
!!$        end do
!!$        out(k,j) = RHO(k,j,iplane)*(out(k,j) - SCV*(dUdx_c(k,j,iplane,1,1) + dUdx_c(k,j,iplane,2,2) + dUdx_c(k,j,iplane,3,3)))
!!$        do jj=1,3
!!$           out(k,j) = out(k,j) - SCV*U_c(k,j,iplane,jj)*dRHOdx(k,j,iplane,jj)
!!$        end do
!!$     end do
!!$  end do
!!$  call cmean(iplane, out, out_c)

  ! Physical space
  do j=1,ny
     do k=1,nz
        SCV = SC(k,j,iplane,isc_sc)**2
        out(k,j) = 0.0_WP
        do jj=1,3
           out(k,j) = out(k,j) - 2.0_WP*U(k,j,iplane,jj)*SC(k,j,iplane,isc_sc)*dSCdx(k,j,iplane,isc_sc,jj)
        end do
        out(k,j) = RHO(k,j,iplane)*(out(k,j) - SCV*(dUdx(k,j,iplane,1,1) + dUdx(k,j,iplane,2,2) + dUdx(k,j,iplane,3,3)))
        do jj=1,3
           out(k,j) = out(k,j) - SCV*U(k,j,iplane,jj)*dRHOdx(k,j,iplane,jj)
        end do
     end do
  end do
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, SCV_div1_cond(:,iplane), icount_qSC)
!!$  SCV_div1_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_div1_phys(:,iplane) + out_z)/real(ntime_qSC,WP)

  return
end subroutine planeStats_budgets_scalar_div1


! ========================================================== !
! Divergence term 2 (diffusive flux)                         !
! ========================================================== !
subroutine planeStats_budgets_scalar_div2(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: j,k,jj
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c

!!$  ! Conditional space
!!$  do j=1,ny
!!$     do k=1,nz
!!$        out(k,j) = 0.0_WP
!!$        do jj=1,3
!!$           out(k,j) = out(k,j) + 2.0_WP*diff_flux(k,j,iplane,isc_sc,jj)*dSCdx_c(k,j,iplane,isc_sc,jj)
!!$        end do
!!$        out(k,j) = out(k,j) + 2.0_WP*SC_c(k,j,iplane,isc_sc)*diff_flux_div(k,j,iplane,isc_sc)
!!$     end do
!!$  end do
!!$  call cmean(iplane, out, out_c)

  ! Physical space
  do j=1,ny
     do k=1,nz
        out(k,j) = 0.0_WP
        do jj=1,3
           out(k,j) = out(k,j) + 2.0_WP*diff_flux(k,j,iplane,isc_sc,jj)*dSCdx(k,j,iplane,isc_sc,jj)
        end do
        out(k,j) = out(k,j) + 2.0_WP*SC(k,j,iplane,isc_sc)*diff_flux_div(k,j,iplane,isc_sc)
     end do
  end do
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, SCV_div2_cond(:,iplane), icount_qSC)
!!$  SCV_div2_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_div2_phys(:,iplane) + out_z)/real(ntime_qSC,WP)

  return
end subroutine planeStats_budgets_scalar_div2


! ========================================================== !
! Production term                                            !
! ========================================================== !
subroutine planeStats_budgets_scalar_prod(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: j,k,jj
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c

!!$  ! Conditional space
!!$  do j=1,ny
!!$     do k=1,nz
!!$        ib = bin_index(k,j,iplane)
!!$        out(k,j) = 0.0_WP
!!$        do jj=1,3
!!$           out(k,j) = out(k,j) + U_c(k,j,iplane,jj)*SC_c(k,j,iplane,isc_sc)*dSCdxm_c(ib,iplane,isc_sc,jj)
!!$        end do
!!$        out(k,j) = -2.0_WP*RHO(k,j,iplane)*out(k,j)
!!$     end do
!!$  end do
!!$  call cmean(iplane, out, out_c)

  ! Physical space
  do j=1,ny
     do k=1,nz
        out(k,j) = 0.0_WP
        do jj=1,3
           out(k,j) = out(k,j) + U(k,j,iplane,jj)*SC(k,j,iplane,isc_sc)*dSCdxm(j,iplane,isc_sc,jj)
        end do
        out(k,j) = -2.0_WP*RHO(k,j,iplane)*out(k,j)
     end do
  end do
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, SCV_prod_cond(:,iplane), icount_qSC)
!!$  SCV_prod_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_prod_phys(:,iplane) + out_z)/real(ntime_qSC,WP)

  return
end subroutine planeStats_budgets_scalar_prod


! ========================================================== !
! Dissipation term -- USING DIFFUSIVE FLUX                   !
! ========================================================== !
subroutine planeStats_budgets_scalar_diss(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: j,k,jj
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c

!!$  ! Conditional space
!!$  do j=1,ny
!!$     do k=1,nz
!!$        out(k,j) = 0.0_WP
!!$        do jj=1,3
!!$           out(k,j) = out(k,j) - 2.0_WP*diff_flux(k,j,iplane,isc_sc,jj)*dSCdx_c(k,j,iplane,isc_sc,jj)
!!$        end do
!!$     end do
!!$  end do
!!$  call cmean(iplane, out, out_c)

  ! Physical space
  do j=1,ny
     do k=1,nz
        out(k,j) = 0.0_WP
        do jj=1,3
           out(k,j) = out(k,j) - 2.0_WP*diff_flux(k,j,iplane,isc_sc,jj)*dSCdx(k,j,iplane,isc_sc,jj)
        end do
     end do
  end do
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, SCV_diss_cond(:,iplane), icount_qSC)
!!$  SCV_diss_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_diss_phys(:,iplane) + out_z)/real(ntime_qSC,WP)

  return
end subroutine planeStats_budgets_scalar_diss


! ========================================================== !
! Dissipation term -- USING DIFFUSIVITIES                    !
! ========================================================== !
subroutine planeStats_budgets_scalar_diss2(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: j,k,jj
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c

!!$  ! Conditional space
!!$  do j=1,ny
!!$     do k=1,nz
!!$        out(k,j) = 0.0_WP
!!$        do jj=1,3
!!$           out(k,j) = out(k,j) &
!!$                - 2.0_WP*RHO(k,j,iplane)*DIFF(k,j,iplane,isc_sc)*dSCdx_c(k,j,iplane,isc_sc,jj)**2
!!$        end do
!!$     end do
!!$  end do
!!$  call cmean(iplane, out, out_c)

  ! Physical space
  do j=1,ny
     do k=1,nz
        out(k,j) = 0.0_WP
        do jj=1,3
           out(k,j) = out(k,j) &
                - 2.0_WP*RHO(k,j,iplane)*DIFF(k,j,iplane,isc_sc)*dSCdx(k,j,iplane,isc_sc,jj)**2
        end do
     end do
  end do
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, SCV_diss2_cond(:,iplane), icount_qSC)
!!$  SCV_diss2_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_diss2_phys(:,iplane) + out_z)/real(ntime_qSC,WP)

  return
end subroutine planeStats_budgets_scalar_diss2


! ========================================================== !
! Chemical source term                                       !
! ========================================================== !
subroutine planeStats_budgets_scalar_src(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: j,k
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c

!!$  ! Conditional space
!!$  do j=1,ny
!!$     do k=1,nz
!!$        out(k,j) = 2.0_WP*SC_c(k,j,iplane,isc_sc)*src_SC(k,j,iplane,isc_sc)
!!$     end do
!!$  end do
!!$  call cmean(iplane, out, out_c)

  ! Physical space
  do j=1,ny
     do k=1,nz
        out(k,j) = 2.0_WP*SC(k,j,iplane,isc_sc)*src_SC_fc(k,j,iplane,isc_sc)
        !out(k,j) = 2.0_WP*SC(k,j,iplane,isc_sc)*src_SC(k,j,iplane,isc_sc)
     end do
  end do
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, SCV_src_cond(:,iplane), icount_qSC)
!!$  SCV_src_phys(:,iplane) = (real(ntime_qSC-1,WP)*SCV_src_phys(:,iplane) + out_z)/real(ntime_qSC,WP)

  return
end subroutine planeStats_budgets_scalar_src
