! ========================================================== !
! Compute the spectral budgets using inhomogeneous data      !
! ========================================================== !
subroutine planeStats_spectrum
  use planeStats
  implicit none

  integer :: iunit, iprobe, irec
  integer :: i, j, k, n
  real(WP) :: dk, pi, tmp1, tmp2, tmp3, tmp4
  real(WP), dimension(:), pointer :: Sij_trace, tmpv1, tmpv2, tmpv3
  real(WP), parameter :: delta_L = 0.000443875_WP

  ! Spectrum workspaces
  real(WP), dimension(:,:),   pointer :: E_u1, E_u2, E_u3
  real(WP), dimension(:,:),   pointer :: eps
  real(WP), dimension(:,:),   pointer :: eta
  real(WP), dimension(:,:),   pointer :: t_eta

  ! Output quantities
  real(WP), dimension(:,:),   pointer :: E_mean, E_mean_s
  real(WP), dimension(:),     pointer :: eps_mean
  real(WP), dimension(:),     pointer :: eta_mean
  real(WP), dimension(:,:),   pointer :: kz, kz_s
  real(WP), dimension(:,:),   pointer :: trans
  real(WP), dimension(:,:,:), pointer :: E_span
  !real(WP), dimension(:),     pointer :: t_eta_mean

  ! Allocate remaining arrays
  allocate(eps     (ipmin_:ipmax_,1:nrec))
  allocate(eta     (ipmin_:ipmax_,1:nrec))
  allocate(t_eta   (ipmin_:ipmax_,1:nrec))
  allocate(E_u1    (1:nrec,1:nz))
  allocate(E_u2    (1:nrec,1:nz))
  allocate(E_u3    (1:nrec,1:nz))
  allocate(E_mean  (ipmin_:ipmax_,1:nz))
  allocate(E_mean_s(ipmin_:ipmax_,1:nz))
  allocate(kz      (ipmin_:ipmax_,1:nz))
  allocate(kz_s    (ipmin_:ipmax_,1:nz))
  allocate(eps_mean(ipmin_:ipmax_))
  allocate(eta_mean(ipmin_:ipmax_))
  allocate(trans   (3,1:nz))
  allocate(E_span (3,1:nrec,1:nz))
  allocate(E_conv  (ipmin_:ipmax_,1:nz))
  allocate(E_div1  (ipmin_:ipmax_,1:nz))
  allocate(E_div2  (ipmin_:ipmax_,1:nz))
  allocate(E_div3  (ipmin_:ipmax_,1:nz))
  allocate(E_pres  (ipmin_:ipmax_,1:nz))
  allocate(E_prod  (ipmin_:ipmax_,1:nz))
  allocate(E_diss  (ipmin_:ipmax_,1:nz))
  allocate(Sij_trace(1:nz))
  allocate(tmpv1    (1:nz))
  allocate(tmpv2    (1:nz))
  allocate(tmpv3    (1:nz))

  ! Compute mean and fluctuating quantities
  if (irank.eq.iroot) then
     print *, 'Computing mean and fluctuating quantities'
  end if
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(RHO(iprobe,irec,1:nz), RHOm(iprobe,irec))
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(VISC(iprobe,irec,1:nz), NU(iprobe,irec))
        NU(iprobe,irec) = NU(iprobe,irec)/RHOm(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(RHO(iprobe,irec,1:nz)*U(iprobe,irec,1:nz), Um(iprobe,irec))
        Um(iprobe,irec) = Um(iprobe,irec)/RHOm(iprobe,irec)
        U(iprobe,irec,1:nz) = U(iprobe,irec,1:nz) - Um(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(RHO(iprobe,irec,1:nz)*V(iprobe,irec,1:nz), Vm(iprobe,irec))
        Vm(iprobe,irec) = Vm(iprobe,irec)/RHOm(iprobe,irec)
        V(iprobe,irec,1:nz) = V(iprobe,irec,1:nz) - Vm(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(RHO(iprobe,irec,1:nz)*W(iprobe,irec,1:nz), Wm(iprobe,irec))
        Wm(iprobe,irec) = Wm(iprobe,irec)/RHOm(iprobe,irec)
        W(iprobe,irec,1:nz) = W(iprobe,irec,1:nz) - Wm(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(            RHO(iprobe,irec,1:nz)* dUdx11(iprobe,irec,1:nz), tmp1)
        call zmean(              U(iprobe,irec,1:nz)*dRHOdx1(iprobe,irec,1:nz), tmp2)
        call zmean(1.0_WP/(RHO(iprobe,irec,1:nz)**2)*dRHOdx1(iprobe,irec,1:nz), tmp3)
        dUdx11m(iprobe,irec) = 1.0_WP/RHOm(iprobe,irec)*(tmp1+tmp2) - Um(iprobe,irec)*RHOm(iprobe,irec)*tmp3
        dUdx11i(iprobe,irec,1:nz) = dUdx11(iprobe,irec,1:nz)
        dUdx11(iprobe,irec,1:nz) = dUdx11(iprobe,irec,1:nz) - dUdx11m(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(            RHO(iprobe,irec,1:nz)* dUdx12(iprobe,irec,1:nz), tmp1)
        call zmean(              U(iprobe,irec,1:nz)*dRHOdx2(iprobe,irec,1:nz), tmp2)
        call zmean(1.0_WP/(RHO(iprobe,irec,1:nz)**2)*dRHOdx2(iprobe,irec,1:nz), tmp3)
        dUdx12m(iprobe,irec) = 1.0_WP/RHOm(iprobe,irec)*(tmp1+tmp2) - Um(iprobe,irec)*RHOm(iprobe,irec)*tmp3
        dUdx12i(iprobe,irec,1:nz) = dUdx12(iprobe,irec,1:nz)
        dUdx12(iprobe,irec,1:nz) = dUdx12(iprobe,irec,1:nz) - dUdx12m(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(            RHO(iprobe,irec,1:nz)* dUdx13(iprobe,irec,1:nz), tmp1)
        call zmean(              U(iprobe,irec,1:nz)*dRHOdx3(iprobe,irec,1:nz), tmp2)
        call zmean(1.0_WP/(RHO(iprobe,irec,1:nz)**2)*dRHOdx3(iprobe,irec,1:nz), tmp3)
        dUdx13m(iprobe,irec) = 1.0_WP/RHOm(iprobe,irec)*(tmp1+tmp2) - Um(iprobe,irec)*RHOm(iprobe,irec)*tmp3
        dUdx13i(iprobe,irec,1:nz) = dUdx13(iprobe,irec,1:nz)
        dUdx13(iprobe,irec,1:nz) = dUdx13(iprobe,irec,1:nz) - dUdx13m(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(            RHO(iprobe,irec,1:nz)* dUdx21(iprobe,irec,1:nz), tmp1)
        call zmean(              V(iprobe,irec,1:nz)*dRHOdx1(iprobe,irec,1:nz), tmp2)
        call zmean(1.0_WP/(RHO(iprobe,irec,1:nz)**2)*dRHOdx1(iprobe,irec,1:nz), tmp3)
        dUdx21m(iprobe,irec) = 1.0_WP/RHOm(iprobe,irec)*(tmp1+tmp2) - Vm(iprobe,irec)*RHOm(iprobe,irec)*tmp3
        dUdx21i(iprobe,irec,1:nz) = dUdx21(iprobe,irec,1:nz)
        dUdx21(iprobe,irec,1:nz) = dUdx21(iprobe,irec,1:nz) - dUdx21m(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(            RHO(iprobe,irec,1:nz)* dUdx22(iprobe,irec,1:nz), tmp1)
        call zmean(              V(iprobe,irec,1:nz)*dRHOdx2(iprobe,irec,1:nz), tmp2)
        call zmean(1.0_WP/(RHO(iprobe,irec,1:nz)**2)*dRHOdx2(iprobe,irec,1:nz), tmp3)
        dUdx22m(iprobe,irec) = 1.0_WP/RHOm(iprobe,irec)*(tmp1+tmp2) - Vm(iprobe,irec)*RHOm(iprobe,irec)*tmp3
        dUdx22i(iprobe,irec,1:nz) = dUdx22(iprobe,irec,1:nz)
        dUdx22(iprobe,irec,1:nz) = dUdx22(iprobe,irec,1:nz) - dUdx22m(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(            RHO(iprobe,irec,1:nz)* dUdx23(iprobe,irec,1:nz), tmp1)
        call zmean(              V(iprobe,irec,1:nz)*dRHOdx3(iprobe,irec,1:nz), tmp2)
        call zmean(1.0_WP/(RHO(iprobe,irec,1:nz)**2)*dRHOdx3(iprobe,irec,1:nz), tmp3)
        dUdx23m(iprobe,irec) = 1.0_WP/RHOm(iprobe,irec)*(tmp1+tmp2) - Vm(iprobe,irec)*RHOm(iprobe,irec)*tmp3
        dUdx23i(iprobe,irec,1:nz) = dUdx23(iprobe,irec,1:nz)
        dUdx23(iprobe,irec,1:nz) = dUdx23(iprobe,irec,1:nz) - dUdx23m(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(            RHO(iprobe,irec,1:nz)* dUdx31(iprobe,irec,1:nz), tmp1)
        call zmean(              W(iprobe,irec,1:nz)*dRHOdx1(iprobe,irec,1:nz), tmp2)
        call zmean(1.0_WP/(RHO(iprobe,irec,1:nz)**2)*dRHOdx1(iprobe,irec,1:nz), tmp3)
        dUdx31m(iprobe,irec) = 1.0_WP/RHOm(iprobe,irec)*(tmp1+tmp2) - Wm(iprobe,irec)*RHOm(iprobe,irec)*tmp3
        dUdx31i(iprobe,irec,1:nz) = dUdx31(iprobe,irec,1:nz)
        dUdx31(iprobe,irec,1:nz) = dUdx31(iprobe,irec,1:nz) - dUdx31m(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(            RHO(iprobe,irec,1:nz)* dUdx32(iprobe,irec,1:nz), tmp1)
        call zmean(              W(iprobe,irec,1:nz)*dRHOdx2(iprobe,irec,1:nz), tmp2)
        call zmean(1.0_WP/(RHO(iprobe,irec,1:nz)**2)*dRHOdx2(iprobe,irec,1:nz), tmp3)
        dUdx32m(iprobe,irec) = 1.0_WP/RHOm(iprobe,irec)*(tmp1+tmp2) - Wm(iprobe,irec)*RHOm(iprobe,irec)*tmp3
        dUdx32i(iprobe,irec,1:nz) = dUdx32(iprobe,irec,1:nz)
        dUdx32(iprobe,irec,1:nz) = dUdx32(iprobe,irec,1:nz) - dUdx32m(iprobe,irec)
     end do
  end do
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call zmean(            RHO(iprobe,irec,1:nz)* dUdx33(iprobe,irec,1:nz), tmp1)
        call zmean(              W(iprobe,irec,1:nz)*dRHOdx3(iprobe,irec,1:nz), tmp2)
        call zmean(1.0_WP/(RHO(iprobe,irec,1:nz)**2)*dRHOdx3(iprobe,irec,1:nz), tmp3)
        dUdx33m(iprobe,irec) = 1.0_WP/RHOm(iprobe,irec)*(tmp1+tmp2) - Wm(iprobe,irec)*RHOm(iprobe,irec)*tmp3
        dUdx33i(iprobe,irec,1:nz) = dUdx33(iprobe,irec,1:nz)
        dUdx33(iprobe,irec,1:nz) = dUdx33(iprobe,irec,1:nz) - dUdx33m(iprobe,irec)
     end do
  end do
  
  ! Compute the dissipation rate from the strain rate tensor
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        Sij_trace = (dUdx11(iprobe,irec,:) + dUdx22(iprobe,irec,:) + dUdx33(iprobe,irec,:))/3.0_WP
        tmpv1     = (dUdx11(iprobe,irec,:) - Sij_trace)*dUdx11(iprobe,irec,:) &
                  + (dUdx22(iprobe,irec,:) - Sij_trace)*dUdx22(iprobe,irec,:) &
                  + (dUdx33(iprobe,irec,:) - Sij_trace)*dUdx33(iprobe,irec,:) &
                  + 0.5_WP*(dUdx12(iprobe,irec,:) + dUdx21(iprobe,irec,:))**2 &
                  + 0.5_WP*(dUdx13(iprobe,irec,:) + dUdx31(iprobe,irec,:))**2 &
                  + 0.5_WP*(dUdx23(iprobe,irec,:) + dUdx32(iprobe,irec,:))**2
        call zmean(VISC(iprobe,irec,1:nz)*tmpv1, eps(iprobe,irec))
        eps(iprobe,irec) = 2.0_WP*eps(iprobe,irec)/RHOm(iprobe,irec)
     end do
  end do

  ! Compute Kolmogorov scales
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        eta(iprobe,irec)   = (NU(iprobe,irec)**3/eps(iprobe,irec))**(0.25)
        t_eta(iprobe,irec) = (NU(iprobe,irec)/eps(iprobe,irec))**(0.5)
     end do
  end do

  ! Compute the energy spectrum
  if (irank.eq.iroot) then
     print *, 'Computing spectral quantities'
  end if
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        call get_spectrum_span(iprobe, irec, U(iprobe,irec,1:nz), E_u1(irec,:))
        call get_spectrum_span(iprobe, irec, V(iprobe,irec,1:nz), E_u2(irec,:))
        call get_spectrum_span(iprobe, irec, W(iprobe,irec,1:nz), E_u3(irec,:))
     end do
     call recmean(RHOm(iprobe,:), tmp4)
     do k=1,nz
        call recmean(E_u1(:,k), tmp1)
        call recmean(E_u2(:,k), tmp2)
        call recmean(E_u3(:,k), tmp3)
        E_mean(iprobe,k) = (tmp1 + tmp2 + tmp3)/(2.0_WP*tmp4)
     end do
  end do

  ! Compute transport terms in spectral space
  ! Mean convection term
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        tmpv1 = RHO(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)**2
        tmpv2 = RHO(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)**2
        tmpv3 = RHO(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)**2
        trans(1,:) = tmpv1*(dUdx11m(iprobe,irec) + dUdx22m(iprobe,irec) + dUdx33m(iprobe,irec))
        trans(2,:) = tmpv2*(dUdx11m(iprobe,irec) + dUdx22m(iprobe,irec) + dUdx33m(iprobe,irec))
        trans(3,:) = tmpv3*(dUdx11m(iprobe,irec) + dUdx22m(iprobe,irec) + dUdx33m(iprobe,irec))
        tmpv1 = RHO(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)*dUdx11(iprobe,irec,1:nz)
        tmpv2 = RHO(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)*dUdx12(iprobe,irec,1:nz)
        tmpv3 = RHO(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)*dUdx13(iprobe,irec,1:nz)
        trans(1,:) = trans(1,:) + 2.0_WP*(Um(iprobe,irec)*tmpv1 + Vm(iprobe,irec)*tmpv2 + Wm(iprobe,irec)*tmpv3)
        tmpv1 = RHO(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)*dUdx21(iprobe,irec,1:nz)
        tmpv2 = RHO(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)*dUdx22(iprobe,irec,1:nz)
        tmpv3 = RHO(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)*dUdx23(iprobe,irec,1:nz)
        trans(2,:) = trans(2,:) + 2.0_WP*(Um(iprobe,irec)*tmpv1 + Vm(iprobe,irec)*tmpv2 + Wm(iprobe,irec)*tmpv3)
        tmpv1 = RHO(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)*dUdx31(iprobe,irec,1:nz)
        tmpv2 = RHO(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)*dUdx32(iprobe,irec,1:nz)
        tmpv3 = RHO(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)*dUdx33(iprobe,irec,1:nz)
        trans(3,:) = trans(3,:) + 2.0_WP*(Um(iprobe,irec)*tmpv1 + Vm(iprobe,irec)*tmpv2 + Wm(iprobe,irec)*tmpv3)
        tmpv1 = U(iprobe,irec,1:nz)**2*dRHOdx1(iprobe,irec,1:nz)
        tmpv2 = U(iprobe,irec,1:nz)**2*dRHOdx2(iprobe,irec,1:nz)
        tmpv3 = U(iprobe,irec,1:nz)**2*dRHOdx3(iprobe,irec,1:nz)
        trans(1,:) = trans(1,:) + Um(iprobe,irec)*tmpv1 + Vm(iprobe,irec)*tmpv2 + Wm(iprobe,irec)*tmpv3
        tmpv1 = V(iprobe,irec,1:nz)**2*dRHOdx1(iprobe,irec,1:nz)
        tmpv2 = V(iprobe,irec,1:nz)**2*dRHOdx2(iprobe,irec,1:nz)
        tmpv3 = V(iprobe,irec,1:nz)**2*dRHOdx3(iprobe,irec,1:nz)
        trans(2,:) = trans(2,:) + Um(iprobe,irec)*tmpv1 + Vm(iprobe,irec)*tmpv2 + Wm(iprobe,irec)*tmpv3
        tmpv1 = W(iprobe,irec,1:nz)**2*dRHOdx1(iprobe,irec,1:nz)
        tmpv2 = W(iprobe,irec,1:nz)**2*dRHOdx2(iprobe,irec,1:nz)
        tmpv3 = W(iprobe,irec,1:nz)**2*dRHOdx3(iprobe,irec,1:nz)
        trans(3,:) = trans(3,:) + Um(iprobe,irec)*tmpv1 + Vm(iprobe,irec)*tmpv2 + Wm(iprobe,irec)*tmpv3
        ! Compute the spectrum
!!$        call get_spectrum_span(iprobe, irec, trans(1,1:nz), E_span(1,irec,:))
!!$        call get_spectrum_span(iprobe, irec, trans(2,1:nz), E_span(2,irec,:))
!!$        call get_spectrum_span(iprobe, irec, trans(3,1:nz), E_span(3,irec,:))
        tmpv1 = trans(1,1:nz) + trans(2,1:nz) + trans(3,1:nz)
        call get_spectrum_span(iprobe, irec, tmpv1, E_span(1,irec,:))
     end do
     call recmean(RHOm(iprobe,:), tmp4)
     do k=1,nz
        call recmean(E_span(1,:,k), tmp1)
!!$        call recmean(E_span(2,:,k), tmp2)
!!$        call recmean(E_span(3,:,k), tmp3)
!!$        E_conv(iprobe,k) = (tmp1 + tmp2 + tmp3)/(2.0_WP*tmp4)
        E_conv(iprobe,k) = tmp1/tmp4
     end do
  end do
  ! Divergence term 1 (pressure-velocity)
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        trans(1,:) = 2.0_WP*( U(iprobe,irec,1:nz)*dPdx1 (iprobe,irec,1:nz) &
                            + P(iprobe,irec,1:nz)*dUdx11(iprobe,irec,1:nz) )
        trans(2,:) = 2.0_WP*( V(iprobe,irec,1:nz)*dPdx2 (iprobe,irec,1:nz) &
                            + P(iprobe,irec,1:nz)*dUdx22(iprobe,irec,1:nz) )
        trans(3,:) = 2.0_WP*( W(iprobe,irec,1:nz)*dPdx3 (iprobe,irec,1:nz) &
                            + P(iprobe,irec,1:nz)*dUdx33(iprobe,irec,1:nz) )
        call get_spectrum_span(iprobe, irec, trans(1,1:nz)+trans(2,1:nz)+trans(3,1:nz), E_span(1,irec,:))
     end do
     call recmean(RHOm(iprobe,:), tmp4)
     do k=1,nz
        call recmean(E_span(1,:,k), tmp1)
        E_div1(iprobe,k) = tmp1/tmp4
     end do
  end do
  ! Divergence term 2 (fluctuating velocity)
  ! CHECK SUMMATIONS HERE
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        trans(1,:) = 2.0_WP*U(iprobe,irec,1:nz)**2*( U(iprobe,irec,1:nz)*dRHOdx1(iprobe,irec,1:nz) &
                                                   + V(iprobe,irec,1:nz)*dRHOdx2(iprobe,irec,1:nz) &
                                                   + W(iprobe,irec,1:nz)*dRHOdx3(iprobe,irec,1:nz) ) &
                   + 4.0_WP*RHO(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)*( U(iprobe,irec,1:nz)*dUdx11(iprobe,irec,1:nz) &
                                                                      + V(iprobe,irec,1:nz)*dUdx12(iprobe,irec,1:nz) &
                                                                      + W(iprobe,irec,1:nz)*dUdx13(iprobe,irec,1:nz) ) &
                   + 2.0_WP*RHO(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)**2*( dUdx11(iprobe,irec,1:nz) &
                                                                         + dUdx22(iprobe,irec,1:nz) &
                                                                         + dUdx33(iprobe,irec,1:nz) )
        trans(2,:) = 2.0_WP*V(iprobe,irec,1:nz)**2*( U(iprobe,irec,1:nz)*dRHOdx1(iprobe,irec,1:nz) &
                                                   + V(iprobe,irec,1:nz)*dRHOdx2(iprobe,irec,1:nz) &
                                                   + W(iprobe,irec,1:nz)*dRHOdx3(iprobe,irec,1:nz) ) &
                   + 4.0_WP*RHO(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)*( U(iprobe,irec,1:nz)*dUdx21(iprobe,irec,1:nz) &
                                                                      + V(iprobe,irec,1:nz)*dUdx22(iprobe,irec,1:nz) &
                                                                      + W(iprobe,irec,1:nz)*dUdx23(iprobe,irec,1:nz) ) &
                   + 2.0_WP*RHO(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)**2*( dUdx11(iprobe,irec,1:nz) &
                                                                         + dUdx22(iprobe,irec,1:nz) &
                                                                         + dUdx33(iprobe,irec,1:nz) )
        trans(3,:) = 2.0_WP*W(iprobe,irec,1:nz)**2*( U(iprobe,irec,1:nz)*dRHOdx1(iprobe,irec,1:nz) &
                                                   + V(iprobe,irec,1:nz)*dRHOdx2(iprobe,irec,1:nz) &
                                                   + W(iprobe,irec,1:nz)*dRHOdx3(iprobe,irec,1:nz) ) &
                   + 4.0_WP*RHO(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)*( U(iprobe,irec,1:nz)*dUdx31(iprobe,irec,1:nz) &
                                                                      + V(iprobe,irec,1:nz)*dUdx32(iprobe,irec,1:nz) &
                                                                      + W(iprobe,irec,1:nz)*dUdx33(iprobe,irec,1:nz) ) &
                   + 2.0_WP*RHO(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)**2*( dUdx11(iprobe,irec,1:nz) &
                                                                         + dUdx22(iprobe,irec,1:nz) &
                                                                         + dUdx33(iprobe,irec,1:nz) )
        call get_spectrum_span(iprobe, irec, trans(1,1:nz)+trans(2,1:nz)+trans(3,1:nz), E_span(1,irec,:))
     end do
     call recmean(RHOm(iprobe,:), tmp4)
     do k=1,nz
        call recmean(E_span(1,:,k), tmp1)
        E_div2(iprobe,k) = tmp1/tmp4
     end do
  end do
  ! Divergence term 3 (viscous stresses)
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        ! tmpv becomes the dissipation rate (same as used later)
        tmpv3 = 2.0_WP/3.0_WP*( dUdx11i(iprobe,irec,1:nz) + dUdx22i(iprobe,irec,1:nz) + dUdx33i(iprobe,irec,1:nz) )
        tmpv1 = 2.0_WP*VISC(iprobe,irec,1:nz)*( dUdx11(iprobe,irec,1:nz)*(dUdx11i(iprobe,irec,1:nz) + dUdx11i(iprobe,irec,1:nz) - tmpv3) &
                                              + dUdx12(iprobe,irec,1:nz)*(dUdx12i(iprobe,irec,1:nz) + dUdx21i(iprobe,irec,1:nz)) &
                                              + dUdx13(iprobe,irec,1:nz)*(dUdx13i(iprobe,irec,1:nz) + dUdx31i(iprobe,irec,1:nz)) )
        trans(1,:) = tmpv1 + 2.0_WP*U(iprobe,irec,1:nz)*( dTAUdx11(iprobe,irec,1:nz) &
                                                        + dTAUdx12(iprobe,irec,1:nz) &
                                                        + dTAUdx13(iprobe,irec,1:nz) )
        tmpv2 = 2.0_WP*VISC(iprobe,irec,1:nz)*( dUdx21(iprobe,irec,1:nz)*(dUdx21i(iprobe,irec,1:nz) + dUdx12i(iprobe,irec,1:nz)) &
                                              + dUdx22(iprobe,irec,1:nz)*(dUdx22i(iprobe,irec,1:nz) + dUdx22i(iprobe,irec,1:nz) - tmpv3) &
                                              + dUdx23(iprobe,irec,1:nz)*(dUdx23i(iprobe,irec,1:nz) + dUdx32i(iprobe,irec,1:nz)) )
        trans(2,:) = tmpv2 + 2.0_WP*V(iprobe,irec,1:nz)*( dTAUdx21(iprobe,irec,1:nz) &
                                                        + dTAUdx22(iprobe,irec,1:nz) &
                                                        + dTAUdx23(iprobe,irec,1:nz) )
        tmpv3 = 2.0_WP*VISC(iprobe,irec,1:nz)*( dUdx31(iprobe,irec,1:nz)*(dUdx31i(iprobe,irec,1:nz) + dUdx13i(iprobe,irec,1:nz)) &
                                              + dUdx32(iprobe,irec,1:nz)*(dUdx32i(iprobe,irec,1:nz) + dUdx23i(iprobe,irec,1:nz)) &
                                              + dUdx33(iprobe,irec,1:nz)*(dUdx33i(iprobe,irec,1:nz) + dUdx33i(iprobe,irec,1:nz) - tmpv3))
        trans(3,:) = tmpv3 + 2.0_WP*W(iprobe,irec,1:nz)*( dTAUdx31(iprobe,irec,1:nz) &
                                                        + dTAUdx32(iprobe,irec,1:nz) &
                                                        + dTAUdx33(iprobe,irec,1:nz) )
        call get_spectrum_span(iprobe, irec, trans(1,1:nz)+trans(2,1:nz)+trans(3,1:nz), E_span(1,irec,:))
     end do
     call recmean(RHOm(iprobe,:), tmp4)
     do k=1,nz
        call recmean(E_span(1,:,k), tmp1)
        E_div3(iprobe,k) = tmp1/tmp4
     end do
  end do
  ! Pressure divergence term
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        trans(1,:) = 2.0_WP*P(iprobe,irec,1:nz)*dUdx11(iprobe,irec,1:nz)
        trans(2,:) = 2.0_WP*P(iprobe,irec,1:nz)*dUdx22(iprobe,irec,1:nz)
        trans(3,:) = 2.0_WP*P(iprobe,irec,1:nz)*dUdx33(iprobe,irec,1:nz)
!!$        call get_spectrum_span(iprobe, irec, trans(1,1:nz), E_span(1,irec,:))
!!$        call get_spectrum_span(iprobe, irec, trans(2,1:nz), E_span(2,irec,:))
!!$        call get_spectrum_span(iprobe, irec, trans(3,1:nz), E_span(3,irec,:))
        call get_spectrum_span(iprobe, irec, trans(1,1:nz)+trans(2,1:nz)+trans(3,1:nz), E_span(1,irec,:))
     end do
     call recmean(RHOm(iprobe,:), tmp4)
     do k=1,nz
        call recmean(E_span(1,:,k), tmp1)
!!$        call recmean(E_span(2,:,k), tmp2)
!!$        call recmean(E_span(3,:,k), tmp3)
!!$        E_pres(iprobe,k) = (tmp1 + tmp2 + tmp3)/(2.0_WP*tmp4)
        E_pres(iprobe,k) = tmp1/tmp4
     end do
  end do
  ! TKE production term
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        tmpv1 = RHO(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)
        tmpv2 = RHO(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)
        tmpv3 = RHO(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)
        trans(1,:) = 2.0_WP*(tmpv1*dUdx11m(iprobe,irec) + tmpv2*dUdx12m(iprobe,irec) + tmpv3*dUdx13m(iprobe,irec))
        tmpv1 = RHO(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)
        tmpv2 = RHO(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)
        tmpv3 = RHO(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)
        trans(2,:) = 2.0_WP*(tmpv1*dUdx21m(iprobe,irec) + tmpv2*dUdx22m(iprobe,irec) + tmpv3*dUdx23m(iprobe,irec))
        tmpv1 = RHO(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)*U(iprobe,irec,1:nz)
        tmpv2 = RHO(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)*V(iprobe,irec,1:nz)
        tmpv3 = RHO(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)*W(iprobe,irec,1:nz)
        trans(3,:) = 2.0_WP*(tmpv1*dUdx31m(iprobe,irec) + tmpv2*dUdx32m(iprobe,irec) + tmpv3*dUdx33m(iprobe,irec))
        ! Compute the spectrum
!!$        call get_spectrum_span(iprobe, irec, trans(1,1:nz), E_span(1,irec,:))
!!$        call get_spectrum_span(iprobe, irec, trans(2,1:nz), E_span(2,irec,:))
!!$        call get_spectrum_span(iprobe, irec, trans(3,1:nz), E_span(3,irec,:))
        tmpv1 = trans(1,1:nz) + trans(2,1:nz) + trans(3,1:nz)
        call get_spectrum_span(iprobe, irec, tmpv1, E_span(1,irec,:))
     end do
     call recmean(RHOm(iprobe,:), tmp4)
     do k=1,nz
        call recmean(E_span(1,:,k), tmp1)
!!$        call recmean(E_span(2,:,k), tmp2)
!!$        call recmean(E_span(3,:,k), tmp3)
!!$        E_prod(iprobe,k) = (tmp1 + tmp2 + tmp3)/(2.0_WP*tmp4)
        E_prod(iprobe,k) = tmp1/tmp4
     end do
  end do
  ! TKE dissipation term
  do iprobe=ipmin_,ipmax_
     do irec=1,nrec
        tmpv3 = 2.0_WP/3.0_WP*( dUdx11i(iprobe,irec,1:nz) + dUdx22i(iprobe,irec,1:nz) + dUdx33i(iprobe,irec,1:nz) )
        tmpv1 = 2.0_WP*VISC(iprobe,irec,1:nz)*( dUdx11(iprobe,irec,1:nz)*(dUdx11i(iprobe,irec,1:nz) + dUdx11i(iprobe,irec,1:nz) - tmpv3) &
                                              + dUdx12(iprobe,irec,1:nz)*(dUdx12i(iprobe,irec,1:nz) + dUdx21i(iprobe,irec,1:nz)) &
                                              + dUdx13(iprobe,irec,1:nz)*(dUdx13i(iprobe,irec,1:nz) + dUdx31i(iprobe,irec,1:nz)) )
        tmpv2 = 2.0_WP*VISC(iprobe,irec,1:nz)*( dUdx21(iprobe,irec,1:nz)*(dUdx21i(iprobe,irec,1:nz) + dUdx12i(iprobe,irec,1:nz)) &
                                              + dUdx22(iprobe,irec,1:nz)*(dUdx22i(iprobe,irec,1:nz) + dUdx22i(iprobe,irec,1:nz) - tmpv3) &
                                              + dUdx23(iprobe,irec,1:nz)*(dUdx23i(iprobe,irec,1:nz) + dUdx32i(iprobe,irec,1:nz)) )
        tmpv3 = 2.0_WP*VISC(iprobe,irec,1:nz)*( dUdx31(iprobe,irec,1:nz)*(dUdx31i(iprobe,irec,1:nz) + dUdx13i(iprobe,irec,1:nz)) &
                                              + dUdx32(iprobe,irec,1:nz)*(dUdx32i(iprobe,irec,1:nz) + dUdx23i(iprobe,irec,1:nz)) &
                                              + dUdx33(iprobe,irec,1:nz)*(dUdx33i(iprobe,irec,1:nz) + dUdx33i(iprobe,irec,1:nz) - tmpv3))
!!$        call get_spectrum_span(iprobe, irec, tmpv1, E_span(1,irec,:))
!!$        call get_spectrum_span(iprobe, irec, tmpv2, E_span(2,irec,:))
!!$        call get_spectrum_span(iprobe, irec, tmpv3, E_span(3,irec,:))
        tmpv1 = tmpv1 + tmpv2 + tmpv3
        call get_spectrum_span(iprobe, irec, tmpv1, E_span(1,irec,:))
     end do
     call recmean(RHOm(iprobe,:), tmp4)
     do k=1,nz
        call recmean(E_span(1,:,k), tmp1)
!!$        call recmean(E_span(2,:,k), tmp2)
!!$        call recmean(E_span(3,:,k), tmp3)
!!$        E_diss(iprobe,k) = (tmp1 + tmp2 + tmp3)/(2.0_WP*tmp4)
        E_diss(iprobe,k) = tmp1/tmp4
     end do
  end do

  ! Average over the symmetry planes
  nPROG = nPROG/2
  iprobe = 1
  do i=1,nplanes
     do j=1,nPROG
        eps(iprobe,:) = 0.5_WP*(eps(iprobe,:) + eps(iprobe+nPROG,:))
        eta(iprobe,:) = 0.5_WP*(eta(iprobe,:) + eta(iprobe+nPROG,:))
        E_mean(iprobe,:) = 0.5_WP*(E_mean(iprobe,:) + E_mean(iprobe+nPROG,:))
        E_conv(iprobe,:) = 0.5_WP*(E_conv(iprobe,:) + E_conv(iprobe+nPROG,:))
        E_div1(iprobe,:) = 0.5_WP*(E_div1(iprobe,:) + E_div1(iprobe+nPROG,:))
        E_div2(iprobe,:) = 0.5_WP*(E_div2(iprobe,:) + E_div2(iprobe+nPROG,:))
        E_div3(iprobe,:) = 0.5_WP*(E_div3(iprobe,:) + E_div3(iprobe+nPROG,:))
        E_pres(iprobe,:) = 0.5_WP*(E_pres(iprobe,:) + E_pres(iprobe+nPROG,:))
        E_prod(iprobe,:) = 0.5_WP*(E_prod(iprobe,:) + E_prod(iprobe+nPROG,:))
        E_diss(iprobe,:) = 0.5_WP*(E_diss(iprobe,:) + E_diss(iprobe+nPROG,:))
        iprobe = iprobe+1
     end do
  end do
  nprobes = nprobes/2

  ! Compute wavenumbers
  pi = acos(-1.0_WP)
  dk = 2.0_WP*pi/Lz
  do k=1,nz
     kz(:,k) = dk*(k-1)
  end do

  ! Average over the records
  do iprobe=1,nprobes
     call recmean(eps(iprobe,1:nrec), eps_mean(iprobe))
     !print *, eps_mean(iprobe)
  end do
  do iprobe=1,nprobes
     call recmean(eta(iprobe,1:nrec), eta_mean(iprobe))
     !print *, eta_mean(iprobe)
  end do

  ! Scale everything
  do iprobe=1,nprobes
     kz_s(iprobe,:)     = kz(iprobe,:)*eta_mean(iprobe)
     E_mean_s(iprobe,:) = E_mean(iprobe,:)*eps_mean(iprobe)**(-2/3)*eta_mean(iprobe)**(-5/3)
  end do

  ! Write some data to stdout

  ! Save the spectrum
  if (irank.eq.iroot) then
     ! Gather the info from all processes
     print *, 'Writing the output files'
     open(unit=iunit, file=trim(output_name), action='write')
     do k=1,nz/2
        do iprobe=1,nprobes
           if (iprobe.lt.nprobes) then
              write(iunit,'(ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13)',advance='no') kz_s(iprobe,k), E_mean_s(iprobe,k), E_conv(iprobe,k), E_div1(iprobe,k), E_div2(iprobe,k), E_div3(iprobe,k), E_pres(iprobe,k), E_prod(iprobe,k), E_diss(iprobe,k)
           else
              write(iunit,'(ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13)') kz_s(iprobe,k), E_mean_s(iprobe,k), E_conv(iprobe,k), E_div1(iprobe,k), E_div2(iprobe,k), E_div3(iprobe,k), E_pres(iprobe,k), E_prod(iprobe,k), E_diss(iprobe,k)
           end if
        end do
     end do
     close(iunit)
  end if

  ! Flame thickness scaling
  do iprobe=1,nprobes
     kz_s(iprobe,:)     = kz(iprobe,:)*delta_L
     E_mean_s(iprobe,:) = E_mean(iprobe,:)*eps_mean(iprobe)**(-2/3)*delta_L**(-5/3)
  end do
  
  if (irank.eq.iroot) then
     ! FIX THIS FOR MPI!!
     open(unit=iunit, file=trim(output_name)//'-flame', action='write')
     do k=1,nz/2
        do iprobe=1,nprobes
           if (iprobe.lt.nprobes) then
              write(iunit,'(ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13)',advance='no') kz_s(iprobe,k), E_mean_s(iprobe,k), E_conv(iprobe,k), E_div1(iprobe,k), E_div2(iprobe,k), E_div3(iprobe,k), E_pres(iprobe,k), E_prod(iprobe,k), E_diss(iprobe,k)
           else
              write(iunit,'(ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13, ES22.13)') kz_s(iprobe,k), E_mean_s(iprobe,k), E_conv(iprobe,k), E_div1(iprobe,k), E_div2(iprobe,k), E_div3(iprobe,k), E_pres(iprobe,k), E_prod(iprobe,k), E_diss(iprobe,k)
           end if
        end do
     end do
     close(iunit)

     print *, 'Done'
  end if

end subroutine planeStats_spectrum


! ========================================================== !
! Compute the density-weighted spectrum of 1-D data          !
! ========================================================== !
subroutine get_spectrum_span(iprobe,irec,A,S)
  use planeStats
  implicit none
  include 'fftw3.f'

  real(WP), dimension(nz), intent(in) :: A
  real(WP), dimension(nz), intent(out) :: S
  complex(WP), dimension(:), pointer :: in, out1, out2
  integer(KIND=8) :: plan
  integer :: iprobe, irec, k
  complex(WP), parameter :: ii = (0.0_WP,1.0_WP)

  ! Create the plan
  allocate(in  (nz))
  allocate(out1(nz))
  allocate(out2(nz))
  call dfftw_plan_dft_1d(plan, nz, in, out2, FFTW_FORWARD, FFTW_ESTIMATE)

  ! FFT
  in = A*RHO(iprobe,irec,1:nz)
  call dfftw_execute(plan)
  out1 = out2/nz
  in = A
  call dfftw_execute(plan)
  out2 = out2/nz
  do k=1,nz
     !S(k) = (0.5_WP*( real(out1(k)*conjg(out2(k))) &
     S(k) = sqrt(abs(0.5_WP*( real(out1(k)*conjg(out2(k))) &
                            + real(out2(k)*conjg(out1(k))))))
  end do

  ! Clean up
  call dfftw_destroy_plan(plan)
  deallocate(in)
  deallocate(out1)
  deallocate(out2)

  return
end subroutine get_spectrum_span
