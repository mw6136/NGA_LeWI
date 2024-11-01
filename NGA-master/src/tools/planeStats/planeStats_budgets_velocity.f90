subroutine planeStats_budgets_velocity_compute
  use planeStats_budgets
  implicit none

  integer :: iplane, ii, jj, rndx
  integer, dimension(1:6) :: indx, jndx
  real(WP), dimension(1:ny) :: tmpv1, tmpv2
  real(WP), dimension(1:nbins_cond) :: tmpc1
  real(WP), dimension(1:nz,1:ny) :: tmpxy

  ! Indices for six independent components of Reynolds stress
  indx = (/ 1, 2, 3, 1, 1, 2 /)
  jndx = (/ 1, 2, 3, 2, 3, 3 /)

  ! TKE and Unsteady Term (RHO*TKE)
  ! Physical space
  do iplane=pnmin_,pnmax_
     ! TKE
     tmpxy = RHO(:,:,iplane)*( U(:,:,iplane,1)**2 + U(:,:,iplane,2)**2 + U(:,:,iplane,3)**2 )
     call condmean(iplane, tmpxy, E_TKE_phys(:,iplane), E_TKE_cond(:,iplane), icount_all)

     ! Save unsteady term
     if (ntime_curr.gt.1) then
        tmpxy = (tmpxy - tmpv_save_TKE_phys(:,:,iplane))/(time(ntime_curr) - time(ntime_curr-1))
        call condmean(iplane, tmpxy, E_unst_phys(:,iplane), E_unst_cond(:,iplane), icount_all)
     end if

     ! Save for next time step
     tmpv_save_TKE_phys(:,:,iplane) = RHO(:,:,iplane)*( U(:,:,iplane,1)**2 + U(:,:,iplane,2)**2 + U(:,:,iplane,3)**2 )
  end do

!!$  ! Conditional space
!!$  do iplane=pnmin_,pnmax_
!!$     ! TKE
!!$     tmpxy = RHO(:,:,iplane)*( U_c(:,:,iplane,1)**2 + U_c(:,:,iplane,2)**2 + U_c(:,:,iplane,3)**2 )
!!$     call cmean(iplane, tmpxy, tmpc1)
!!$     E_TKE_cond(:,iplane) = (real(ntime_curr-1,WP)*E_TKE_cond(:,iplane) + tmpc1)/real(ntime_curr,WP)
!!$
!!$     ! Save unsteady term
!!$     if (ntime_curr.gt.1) then
!!$        tmpxy = (tmpxy - tmpv_save_TKE_cond(:,:,iplane))/(time(ntime_curr) - time(ntime_curr-1))
!!$        call cmean(iplane, tmpxy, tmpc1)
!!$        E_unst_cond(:,iplane) = (real(ntime_curr-2,WP)*E_unst_cond(:,iplane) + tmpc1)/real(ntime_curr-1,WP)
!!$     end if
!!$
!!$     ! Save for next time step
!!$     tmpv_save_TKE_cond(:,:,iplane) = RHO(:,:,iplane)*( U_c(:,:,iplane,1)**2 + U_c(:,:,iplane,2)**2 + U_c(:,:,iplane,3)**2 )
!!$  end do

  ! Mean convection term
  do iplane=pnmin_,pnmax_
     tmpv2 = 0.0_WP
     do rndx=1,6
        call planeStats_budgets_velocity_conv(iplane,indx(rndx),jndx(rndx),tmpv1)
        ! Reynolds stress
        R_stress_budg(:,iplane,rndx,1) = (real(ntime_curr-1,WP)*R_stress_budg(:,iplane,rndx,1) + tmpv1)/real(ntime_curr,WP)
        if (rndx.le.3) then
           tmpv2 = tmpv2 + tmpv1
        end if
     end do
     ! TKE
     ! Physical space
     E_conv_phys(:,iplane) = (real(ntime_curr-1,WP)*E_conv_phys(:,iplane) + tmpv2)/real(ntime_curr,WP)
!!$     ! Conditional space
!!$     E_conv_cond(:,iplane) = (real(ntime_curr-1,WP)*E_conv_cond(:,iplane) + tmpc1)/real(ntime_curr,WP)
  end do

  ! Divergence term 2 (fluctuating transport)
  do iplane=pnmin_,pnmax_
     call planeStats_budgets_velocity_div2(iplane,tmpv1)
     ! Physical space
     E_div2_phys(:,iplane) = (real(ntime_curr-1,WP)*E_div2_phys(:,iplane) + tmpv1)/real(ntime_curr,WP)
!!$     ! Conditional space
!!$     E_div2_cond(:,iplane) = (real(ntime_curr-1,WP)*E_div2_cond(:,iplane) + tmpc1)/real(ntime_curr,WP)
  end do

  ! Divergence term 3 (viscous stresses)
  do iplane=pnmin_,pnmax_
     call planeStats_budgets_velocity_div3(iplane,tmpv1)
     ! Physical space
     E_div3_phys(:,iplane) = (real(ntime_curr-1,WP)*E_div3_phys(:,iplane) + tmpv1)/real(ntime_curr,WP)
!!$     ! Conditional space
!!$     E_div3_cond(:,iplane) = (real(ntime_curr-1,WP)*E_div3_cond(:,iplane) + tmpc1)/real(ntime_curr,WP)
  end do

  ! Velocity-pressure gradient term (lumped)
  do iplane=pnmin_,pnmax_
     call planeStats_budgets_velocity_pres(iplane,tmpv1)
     ! Physical space
     E_pres_phys(:,iplane) = (real(ntime_curr-1,WP)*E_pres_phys(:,iplane) + tmpv1)/real(ntime_curr,WP)
!!$     ! Conditional space
!!$     E_pres_cond(:,iplane) = (real(ntime_curr-1,WP)*E_pres_cond(:,iplane) + tmpc1)/real(ntime_curr,WP)
  end do

  ! Pressure-work term
  do iplane=pnmin_,pnmax_
     call planeStats_budgets_velocity_pwork(iplane,tmpv1)
     E_pwork_phys(:,iplane) = (real(ntime_curr-1,WP)*E_pwork_phys(:,iplane) + tmpv1)/real(ntime_curr,WP)
  end do

  ! Pressure-transport term
  do iplane=pnmin_,pnmax_
     call planeStats_budgets_velocity_ptrans(iplane,tmpv1)
     E_ptrans_phys(:,iplane) = (real(ntime_curr-1,WP)*E_ptrans_phys(:,iplane) + tmpv1)/real(ntime_curr,WP)
  end do

  ! Pressure-dilatation term
  do iplane=pnmin_,pnmax_
     call planeStats_budgets_velocity_pdil(iplane,tmpv1)
     E_pdil_phys(:,iplane) = (real(ntime_curr-1,WP)*E_pdil_phys(:,iplane) + tmpv1)/real(ntime_curr,WP)
  end do

  ! TKE production term
  do iplane=pnmin_,pnmax_
     call planeStats_budgets_velocity_prod(iplane,tmpv1)
     ! Physical space
     E_prod_phys(:,iplane) = (real(ntime_curr-1,WP)*E_prod_phys(:,iplane) + tmpv1)/real(ntime_curr,WP)
!!$     ! Conditional space
!!$     E_prod_cond(:,iplane) = (real(ntime_curr-1,WP)*E_prod_cond(:,iplane) + tmpc1)/real(ntime_curr,WP)
  end do

  ! TKE dissipation term
  do iplane=pnmin_,pnmax_
     call planeStats_budgets_velocity_diss(iplane,tmpv1)
     ! Physical space
     E_diss_phys(:,iplane) = (real(ntime_curr-1,WP)*E_diss_phys(:,iplane) + tmpv1)/real(ntime_curr,WP)
!!$     ! Conditional space
!!$     E_diss_cond(:,iplane) = (real(ntime_curr-1,WP)*E_diss_cond(:,iplane) + tmpc1)/real(ntime_curr,WP)
     ! Kolmogorov scale
     eta_phys(:,iplane) = (real(ntime_curr-1,WP)*eta_phys(:,iplane) + sqrt(sqrt(NU  (:,iplane)**3/abs(tmpv1/RHOm  (:,iplane)))) )/real(ntime_curr,WP)
!!$     eta_cond(:,iplane) = (real(ntime_curr-1,WP)*eta_cond(:,iplane) + sqrt(sqrt(NU_c(:,iplane)**3/abs(tmpc1/RHOm_c(:,iplane)))) )/real(ntime_curr,WP)
  end do

  return
end subroutine planeStats_budgets_velocity_compute



! ========================================================== !
! Mean convection term                                       !
! ========================================================== !
subroutine planeStats_budgets_velocity_conv(iplane,ii,jj,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: j,k,kk,ii,jj
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c
  real(WP) :: tke

!!$  ! Conditional space
!!$  do j=1,ny
!!$     do k=1,nz
!!$        ib  = bin_index(k,j,iplane)
!!$        tke = U_c(k,j,iplane,1)**2 + U_c(k,j,iplane,2)**2 + U_c(k,j,iplane,3)**2
!!$        ! Divergence of Favre-averaged velocity field
!!$        out(k,j) = -tke*(dUdxm_c(ib,iplane,1,1) + dUdxm_c(ib,iplane,2,2) + dUdxm_c(ib,iplane,3,3))
!!$        ! Advection of TKE with the Favre-averaged velocities
!!$        do jj=1,3
!!$           do ii=1,3
!!$              out(k,j) = out(k,j) - 2.0_WP*Um_c(ib,iplane,jj)*U_c(k,j,iplane,ii)*dUdx_c(k,j,iplane,ii,jj)
!!$           end do
!!$        end do
!!$        out(k,j) = out(k,j)*RHO(k,j,iplane)
!!$        ! Density gradient term
!!$        do jj=1,3
!!$           out(k,j) = out(k,j) - tke*Um_c(ib,iplane,jj)*dRHOdx(k,j,iplane,jj)
!!$        end do
!!$     end do
!!$  end do
!!$  call cmean(iplane, out, out_c)
  
!!$  ! Physical space
!!$  do j=1,ny
!!$     do k=1,nz
!!$        tke = U(k,j,iplane,1)**2 + U(k,j,iplane,2)**2 + U(k,j,iplane,3)**2
!!$        ! Divergence of Favre-averaged velocity field
!!$        out(k,j) = -tke*(dUdxm(j,iplane,1,1) + dUdxm(j,iplane,2,2) + dUdxm(j,iplane,3,3))
!!$        ! Advection of TKE with the Favre-averaged velocities
!!$        do jj=1,3
!!$           do ii=1,3
!!$              out(k,j) = out(k,j) - 2.0_WP*Um(j,iplane,jj)*U(k,j,iplane,ii)*dUdx(k,j,iplane,ii,jj)
!!$           end do
!!$        end do
!!$        out(k,j) = out(k,j)*RHO(k,j,iplane)
!!$        ! Density gradient term
!!$        do jj=1,3
!!$           out(k,j) = out(k,j) - tke*Um(j,iplane,jj)*dRHOdx(k,j,iplane,jj)
!!$        end do
!!$     end do
!!$  end do
!!$  call condmean(iplane, out, out_z, out_c)

  ! Reynolds stress
  do j=1,ny
     do k=1,nz
        out(k,j) = 0.0_WP
        do kk=1,3
           out(k,j) = out(k,j) + dUdxm(j,iplane,kk,kk)
        end do
        out(k,j) = -U(k,j,iplane,ii)*U(k,j,iplane,jj)*out(k,j)
        do kk=1,3
           out(k,j) = out(k,j) - Um(j,iplane,kk)* &
                               ( U(k,j,iplane,ii)*dUdx(k,j,iplane,jj,kk) &
                               + U(k,j,iplane,jj)*dUdx(k,j,iplane,ii,kk) )
        end do
        out(k,j) = out(k,j)*RHO(k,j,iplane)
        do kk=1,3
           out(k,j) = out(k,j) - Um(j,iplane,kk)*U(k,j,iplane,ii)*U(k,j,iplane,jj)*dRHOdx(k,j,iplane,kk)
        end do
     end do
  end do
  call zmean(out, out_z)
!!$  out_c = 0.0_WP

  return
end subroutine planeStats_budgets_velocity_conv


!!$  ! Divergence term 1 (pressure-velocity)
!!$  do i1=pnmin_,pnmax_
!!$     do i2=1,nbins_cond
!!$        trans(1,:) = 2.0_WP*( U(i1,bin,1:nrec)*dPdx1 (i1,bin,1:nrec) &
!!$                            + P(i1,bin,1:nrec)*dUdx11(i1,bin,1:nrec) )
!!$        trans(2,:) = 2.0_WP*( V(i1,bin,1:nrec)*dPdx2 (i1,bin,1:nrec) &
!!$                            + P(i1,bin,1:nrec)*dUdx22(i1,bin,1:nrec) )
!!$        trans(3,:) = 2.0_WP*( W(i1,bin,1:nrec)*dPdx3 (i1,bin,1:nrec) &
!!$                            + P(i1,bin,1:nrec)*dUdx33(i1,bin,1:nrec) )
!!$
!!$        tmpv1 = trans(1,1:nrec) + trans(2,1:nrec) + trans(3,1:nrec)
!!$        call recmean(tmpv1, tmp1)
!!$        E_div1(i1,bin) = tmp1/RHOm(i1,bin)
!!$     end do
!!$  end do


! ========================================================== !
! Divergence term 2 (fluctuating transport)                  !
! ========================================================== !
subroutine planeStats_budgets_velocity_div2(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: j,k,ii,jj
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), dimension(1:ny), intent(out) :: out_z
!!$  real(WP), dimension(1:nbins_cond), intent(out) :: out_c
  real(WP) :: tke

!!$  ! Conditional space
!!$  do j=1,ny
!!$     do k=1,nz
!!$        tke = U_c(k,j,iplane,1)**2 + U_c(k,j,iplane,2)**2 + U_c(k,j,iplane,3)**2
!!$        out(k,j) = 0.0_WP
!!$        do jj=1,3
!!$           do ii=1,3
!!$              out(k,j) = out(k,j) - 2.0_WP*U_c(k,j,iplane,jj)*U_c(k,j,iplane,ii)*dUdx_c(k,j,iplane,ii,jj)
!!$           end do
!!$        end do
!!$        out(k,j) = RHO(k,j,iplane)*(out(k,j) - tke*(dUdx_c(k,j,iplane,1,1) + dUdx_c(k,j,iplane,2,2) + dUdx_c(k,j,iplane,3,3)))
!!$        do jj=1,3
!!$           out(k,j) = out(k,j) - tke*U_c(k,j,iplane,jj)*dRHOdx(k,j,iplane,jj)
!!$        end do
!!$     end do
!!$  end do
!!$  call cmean(iplane, out, out_c)

  ! Physical space
  do j=1,ny
     do k=1,nz
        tke = U(k,j,iplane,1)**2 + U(k,j,iplane,2)**2 + U(k,j,iplane,3)**2
        out(k,j) = 0.0_WP
        do jj=1,3
           do ii=1,3
              out(k,j) = out(k,j) - 2.0_WP*U(k,j,iplane,jj)*U(k,j,iplane,ii)*dUdx(k,j,iplane,ii,jj)
           end do
        end do
        out(k,j) = RHO(k,j,iplane)*(out(k,j) - tke*(dUdx(k,j,iplane,1,1) + dUdx(k,j,iplane,2,2) + dUdx(k,j,iplane,3,3)))
        do jj=1,3
           out(k,j) = out(k,j) - tke*U(k,j,iplane,jj)*dRHOdx(k,j,iplane,jj)
        end do
     end do
  end do
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, out_c)

  return
 end subroutine planeStats_budgets_velocity_div2


! ========================================================== !
! Divergence term 3 (viscous stresses)                       !
! ========================================================== !
subroutine planeStats_budgets_velocity_div3(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: ii,jj
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c
  real(WP), dimension(1:nz,1:ny) :: divg

  ! Conditional space
  divg = 2.0_WP/3.0_WP*( dUdxi(:,:,iplane,1,1) + dUdxi(:,:,iplane,2,2) + dUdxi(:,:,iplane,3,3) )
!!$  out = 0.0_WP
!!$  do jj=1,3
!!$     do ii=1,3
!!$        if (ii.eq.jj) then
!!$           out = out + dUdx_c(:,:,iplane,ii,jj)*(dUdxi(:,:,iplane,ii,jj) + dUdxi(:,:,iplane,jj,ii) - divg)
!!$        else
!!$           out = out + dUdx_c(:,:,iplane,ii,jj)*(dUdxi(:,:,iplane,ii,jj) + dUdxi(:,:,iplane,jj,ii))
!!$        end if
!!$     end do
!!$  end do
!!$  out = 2.0_WP*VISC(:,:,iplane)*out
!!$  do ii=1,3
!!$     out = out + 2.0_WP*U_c(:,:,iplane,ii)*( dTAUdx(:,:,iplane,ii,1) + dTAUdx(:,:,iplane,ii,2) + dTAUdx(:,:,iplane,ii,3) )
!!$  end do
!!$  call cmean(iplane, out, out_c)

  ! Physical space
  out = 0.0_WP
  do jj=1,3
     do ii=1,3
        if (ii.eq.jj) then
           out = out + dUdx(:,:,iplane,ii,jj)*(dUdxi(:,:,iplane,ii,jj) + dUdxi(:,:,iplane,jj,ii) - divg)
        else
           out = out + dUdx(:,:,iplane,ii,jj)*(dUdxi(:,:,iplane,ii,jj) + dUdxi(:,:,iplane,jj,ii))
        end if
     end do
  end do
  out = 2.0_WP*VISC(:,:,iplane)*out
  do ii=1,3
     out = out + 2.0_WP*U(:,:,iplane,ii)*( dTAUdx(:,:,iplane,ii,1) + dTAUdx(:,:,iplane,ii,2) + dTAUdx(:,:,iplane,ii,3) )
  end do
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, out_c)

  return
end subroutine planeStats_budgets_velocity_div3


!!$  ! Pressure divergence term
!!$  do i1=pnmin_,pnmax_
!!$     do bin=1,nbins_cond
!!$        trans(1,:) = 2.0_WP*P(i1,bin,1:n3)*dUdx11(i1,bin,1:n3)
!!$        trans(2,:) = 2.0_WP*P(i1,bin,1:n3)*dUdx22(i1,bin,1:n3)
!!$        trans(3,:) = 2.0_WP*P(i1,bin,1:n3)*dUdx33(i1,bin,1:n3)
!!$
!!$        tmpv1 = trans(1,1:n3) + trans(2,1:n3) + trans(3,1:n3)
!!$        call recmean(tmpv1, tmp1)
!!$        E_pres(i1,bin) = tmp1/RHOm(i1,bin)
!!$     end do
!!$  end do


! ========================================================== !
! Pressure terms                                             !
! ========================================================== !
! Velocity-pressure gradient: -2 bar(u_i''* dp/dx_i)
subroutine planeStats_budgets_velocity_pres(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: ii,j,k
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c

!!$  ! Conditional space
!!$  out = 0.0_WP
!!$  do ii=1,3
!!$     do j=1,ny
!!$        do k=1,nz
!!$           out(k,j) = out(k,j) - 2.0_WP*U_c(k,j,iplane,ii)*dPdx(k,j,iplane,ii)
!!$        end do
!!$     end do
!!$  end do
!!$  call cmean(iplane, out, out_c)

  ! Physical space
  out = 0.0_WP
  do ii=1,3
     do j=1,ny
        do k=1,nz
           out(k,j) = out(k,j) - 2.0_WP*U(k,j,iplane,ii)*dPdx(k,j,iplane,ii)
        end do
     end do
  end do
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, out_c)

  return
end subroutine planeStats_budgets_velocity_pres


! Pressure-work: -2 bar(u_i'') d bar(p)/dx_i
subroutine planeStats_budgets_velocity_pwork(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: ii,j,k
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z

  out = 0.0_WP
  do ii=1,3
     do j=1,ny
        do k=1,nz
           out(k,j) = out(k,j) - 2.0_WP*U(k,j,iplane,ii)*dPdxm(j,iplane,ii)
        end do
     end do
  end do
  call zmean(out, out_z)

  return
end subroutine planeStats_budgets_velocity_pwork

! Pressure-transport
subroutine planeStats_budgets_velocity_ptrans(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: ii,j,k
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z

  out = 0.0_WP
  do ii=1,3
     do j=1,ny
        do k=1,nz
           out(k,j) = out(k,j) - 2.0_WP*( U(k,j,iplane,ii)*dPdxp(k,j,iplane,ii) &
                                        + dUdx(k,j,iplane,ii,ii)*Pf(k,j,iplane) )
        end do
     end do
  end do
  call zmean(out, out_z)

  return
end subroutine planeStats_budgets_velocity_ptrans

! Pressure-dilatation
subroutine planeStats_budgets_velocity_pdil(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: ii,j,k
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z

  out = 0.0_WP
  do ii=1,3
     do j=1,ny
        do k=1,nz
           out(k,j) = out(k,j) + 2.0_WP*Pf(k,j,iplane)*dUdx(k,j,iplane,ii,ii)
        end do
     end do
  end do
  call zmean(out, out_z)

  return
end subroutine planeStats_budgets_velocity_pdil



! ========================================================== !
! TKE production term                                        !
! ========================================================== !
subroutine planeStats_budgets_velocity_prod(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: j,k,ii,jj
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c

!!$  ! Conditional space
!!$  do j=1,ny
!!$     do k=1,nz
!!$        ib = bin_index(k,j,iplane)
!!$        out(k,j) = 0.0_WP
!!$        do jj=1,3
!!$           do ii=1,3
!!$              out(k,j) = out(k,j) + U_c(k,j,iplane,ii)*U_c(k,j,iplane,jj)*dUdxm_c(ib,iplane,ii,jj)
!!$           end do
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
           do ii=1,3
              out(k,j) = out(k,j) + U(k,j,iplane,ii)*U(k,j,iplane,jj)*dUdxm(j,iplane,ii,jj)
           end do
        end do
        out(k,j) = -2.0_WP*RHO(k,j,iplane)*out(k,j)
     end do
  end do
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, out_c)

  return
end subroutine planeStats_budgets_velocity_prod


! ========================================================== !
! TKE dissipation term                                       !
! ========================================================== !
subroutine planeStats_budgets_velocity_diss(iplane,out_z)
  use planeStats_budgets
  implicit none

  integer, intent(in) :: iplane
  integer :: ii,jj
  real(WP), dimension(1:nz,1:ny) :: out
  real(WP), intent(out), dimension(1:ny) :: out_z
!!$  real(WP), intent(out), dimension(1:nbins_cond) :: out_c
  real(WP), dimension(1:nz,1:ny) :: divg

  ! Conditional space
  ! tau_{ik}*du_i''/dx_k (stress tensor is instantaneous)
  divg = 2.0_WP/3.0_WP*( dUdxi(:,:,iplane,1,1) + dUdxi(:,:,iplane,2,2) + dUdxi(:,:,iplane,3,3) )

!!$  out = 0.0_WP
!!$  do jj=1,3
!!$     do ii=1,3
!!$        if (ii.eq.jj) then
!!$           out = out + dUdx_c(:,:,iplane,ii,jj)*(dUdxi(:,:,iplane,ii,jj) + dUdxi(:,:,iplane,jj,ii) - divg)
!!$        else
!!$           out = out + dUdx_c(:,:,iplane,ii,jj)*(dUdxi(:,:,iplane,ii,jj) + dUdxi(:,:,iplane,jj,ii))
!!$        end if
!!$     end do
!!$  end do
!!$  out = -2.0_WP*VISC(:,:,iplane)*out
!!$  call cmean(iplane, out, out_c)

  ! Physical space
!!$  divg = 2.0_WP/3.0_WP*( dUdxi(:,:,iplane,1,1) + dUdxi(:,:,iplane,2,2) + dUdxi(:,:,iplane,3,3) )
  out = 0.0_WP
  do jj=1,3
     do ii=1,3
        if (ii.eq.jj) then
           out = out + dUdx(:,:,iplane,ii,jj)*(dUdxi(:,:,iplane,ii,jj) + dUdxi(:,:,iplane,jj,ii) - divg)
        else
           out = out + dUdx(:,:,iplane,ii,jj)*(dUdxi(:,:,iplane,ii,jj) + dUdxi(:,:,iplane,jj,ii))
        end if
     end do
  end do
  out = -2.0_WP*VISC(:,:,iplane)*out
  call zmean(out, out_z)
!!$  call condmean(iplane, out, out_z, out_c)

  return
end subroutine planeStats_budgets_velocity_diss
