module volumeStats_velocity
  use volumeStats_average

  ! For velocity residuals
  real(WP), dimension(:,:,:),     pointer :: rhoU, rhoV, rhoW
  real(WP), dimension(:,:,:),     pointer :: rhoUi, rhoVi, rhoWi
  real(WP), dimension(:,:,:),     pointer :: FX, FY, FZ

end module volumeStats_velocity

subroutine volumeStats_velocity_init
  use volumeStats_velocity
  implicit none

  ! For velocity residuals
  allocate(FX     (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(FY     (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(FZ     (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(rhoU   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(rhoV   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(rhoW   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(rhoUi  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(rhoVi  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(rhoWi  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))

  return
end subroutine volumeStats_velocity_init

! ========================================================== !
! Compute the velocity budgets using NGA metrics             !
! ========================================================== !
subroutine volumeStats_velocity_compute
  use volumeStats_velocity
  use volumeStats_metric
  implicit none

  integer :: i, j, k, ii, jj, kk, isc, n, st
  real(WP) :: RHOi, rhs, tmp2, tmp3, tmp_nzi
  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_) :: tmpxyz,tmpxyz1,tmpxyz2,tmpxyz3,tmpxyz4,tmpxyz5,tmpxyz6,tmpxyz7,tmpxyz8
  real(WP), dimension(imin_:imax_,jmin_:jmax_) :: tmpxy1,tmpxy2,tmpxy3,tmpxy4,tmpxy5
  
  ! rho multiply
  !$OMP PARALLEL PRIVATE(j,i)
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           rhoU(i,j,k) = U(i,j,k,1) * sum(interp_sc_x(i,j,:)*RHO(i-st2:i+st1,j,k))
        end do
     end do
  end do
  !$OMP END DO
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           rhoV(i,j,k) = U(i,j,k,2) * sum(interp_sc_y(i,j,:)*RHO(i,j-st2:j+st1,k))
        end do
     end do
  end do
  !$OMP END DO
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           rhoW(i,j,k) = U(i,j,k,3) * sum(interp_sc_z(i,j,:)*RHO(i,j,k-st2:k+st1))
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Do the temporal update at the midpoint of the current and previous timesteps
  if (ntime_curr.gt.1) then
     st = 2
     tmp_nzi = 1.0_WP/real(kmax_-st-(kmin_+st)+1,WP)
     tmp3 = 1.0_WP/dt
     
     !$OMP PARALLEL PRIVATE(j,i)
     ! Density weighted for unconditional budgets
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmpxyz1(i,j,k) = -(rhoU(i,j,k) - U_save_temp(i,j,k,1))*tmp3
           end do
        end do
     end do
     !$OMP END DO
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmpxyz2(i,j,k) = -(rhoV(i,j,k) - U_save_temp(i,j,k,2))*tmp3
           end do
        end do
     end do
     !$OMP END DO
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmpxyz8(i,j,k) = -(rhoW(i,j,k) - U_save_temp(i,j,k,3))*tmp3
           end do
        end do
     end do
     !$OMP END DO
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmpxyz3(i,j,k) = -(RHO(i,j,k)*(Ui(i,j,k,1)**2 + Ui(i,j,k,2)**2 + Ui(i,j,k,3)**2) - E_save_temp(i,j,k))*tmp3
           end do
        end do
     end do
     !$OMP END DO
     ! Density-unweighted for conditional budgets
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmpxyz4(i,j,k) = -(U(i,j,k,1) - Uc_save_temp(i,j,k,1))*tmp3
           end do
        end do
     end do
     !$OMP END DO
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmpxyz5(i,j,k) = -(U(i,j,k,2) - Uc_save_temp(i,j,k,2))*tmp3
           end do
        end do
     end do
     !$OMP END DO
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              ! Nonconservative form for conditional velocity budgets
              tmpxyz6(i,j,k) = -(Ui(i,j,k,1)**2 + Ui(i,j,k,2)**2 + Ui(i,j,k,3)**2 - Ec_save_temp(i,j,k))*tmp3
           end do
        end do
     end do
     !$OMP END DO
     ! Just density for dRHO/dt
     !$OMP DO
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmpxyz7(i,j,k) = (RHO(i,j,k) - RHO_save_temp(i,j,k))*tmp3
           end do
        end do
     end do
     !$OMP END DO
     !$OMP END PARALLEL

     ! Physical average
     tmpxy1 = U_budget_phys_cf(:,:,1,0)*real(ntime_curr-2,WP)
     tmpxy2 = U_budget_phys_cf(:,:,2,0)*real(ntime_curr-2,WP)
     tmpxy3 = E_temp_m*real(ntime_curr-2,WP)
     tmpxy4 = U_temp_m(:,:,1)*real(ntime_curr-2,WP)
     tmpxy5 = U_temp_m(:,:,2)*real(ntime_curr-2,WP)
     tmp2 = 1.0_WP/real(ntime_curr-1,WP)
     do j=jmin_,jmax_
        do i=imin_,imax_
           U_budget_phys_cf(i,j,1,0) = (tmpxy1(i,j) + sum(tmpxyz1(i,j,kmin_+st:kmax_-st))*tmp_nzi)*tmp2
           U_budget_phys_cf(i,j,2,0) = (tmpxy2(i,j) + sum(tmpxyz2(i,j,kmin_+st:kmax_-st))*tmp_nzi)*tmp2 ! fixed the k-range here 1/8/18
           E_temp_m(i,j)             = (tmpxy3(i,j) + sum(tmpxyz3(i,j,kmin_+st:kmax_-st))*tmp_nzi)*tmp2
           U_temp_m(i,j,1)           = (tmpxy4(i,j) + sum(tmpxyz4(i,j,kmin_+st:kmax_-st))*tmp_nzi)*tmp2
           U_temp_m(i,j,2)           = (tmpxy5(i,j) + sum(tmpxyz5(i,j,kmin_+st:kmax_-st))*tmp_nzi)*tmp2
        end do
     end do

     ! Conditional average
     call cond_avg_z(tmpxyz1, U_cond_temp(:,:,:,1), 2)
     call cond_avg_z(tmpxyz2, U_cond_temp(:,:,:,2), 2)
     call cond_avg_z(tmpxyz8, U_cond_temp(:,:,:,3), 2)
     call cond_avg_z(tmpxyz3, E_cond_temp(:,:,:),   2)
     call cond_avg_z(tmpxyz7, RHO_cond_temp(:,:,:), 2)
  end if

  ! Save current velocity and KE for next time step
  !$OMP PARALLEL
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           U_save_temp(i,j,k,1) = rhoU(i,j,k)
           U_save_temp(i,j,k,2) = rhoV(i,j,k)
           U_save_temp(i,j,k,3) = rhoW(i,j,k)
           E_save_temp(i,j,k  ) = RHO(i,j,k)*(Ui(i,j,k,1)**2 + Ui(i,j,k,2)**2 + Ui(i,j,k,3)**2)
           Uc_save_temp(i,j,k,1) = U(i,j,k,1)
           Uc_save_temp(i,j,k,2) = U(i,j,k,2)
           Ec_save_temp(i,j,k  ) = (Ui(i,j,k,1)**2 + Ui(i,j,k,2)**2 + Ui(i,j,k,3)**2)
           RHO_save_temp(i,j,k)  = RHO(i,j,k)
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Enable below to test single-snapshot budgets
!!$  if (ntime_curr.eq.1) then
!!$     print *, 'first dt'
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           U_budget_phys_cf(i,j,1,2) = rhoU(i,j,20)
!!$           !U_budget_phys   (i,j,1,2) = RHO(i,j,20)*Ui(i,j,20,1)
!!$        end do
!!$     end do
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           U_budget_phys_cf(i,j,2,2) = rhoV(i,j,20)
!!$           !U_budget_phys   (i,j,2,2) = RHO(i,j,20)*Ui(i,j,20,2)
!!$        end do
!!$     end do
!!$  else if (ntime_curr.eq.3) then
!!$     print *, 'third dt'
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           U_budget_phys_cf(i,j,1,2) = -(rhoU(i,j,20) - U_budget_phys_cf(i,j,1,2))/(2.0_WP*dt)
!!$           !U_budget_phys   (i,j,1,2) = -(RHO(i,j,20)*Ui(i,j,20,1) - U_budget_phys(i,j,1,2))/(2.0_WP*dt)
!!$        end do
!!$     end do
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           U_budget_phys_cf(i,j,2,2) = -(rhoV(i,j,20) - U_budget_phys_cf(i,j,2,2))/(2.0_WP*dt)
!!$           !U_budget_phys   (i,j,2,2) = -(RHO(i,j,20)*Ui(i,j,20,2) - U_budget_phys(i,j,2,2))/(2.0_WP*dt)
!!$        end do
!!$     end do
!!$  end if


  ! x-Momentum ------------------------------------------------------------------------
  ! Convective part of velocity residual
  !$OMP PARALLEL PRIVATE(i,j,k,n,RHOi,st,rhs)
  !$OMP DO
  do kk=kmin_-stc1,kmax_+stc2
     do jj=jmin_-stc1,jmax_+stc2
        do ii=imin_-stc1,imax_+stc2

           i = ii-1; j = jj-1; k = kk-1;

           rhoUi(i,j,k) = sum(interp_u_xm(i,j,:)*rhoU(i-stc1:i+stc2,j,k))
           
           i = ii; j = jj; k = kk;
           
           rhoVi(i,j,k) = sum(interp_Jv_x(i,j,:)*rhoV(i-stc2:i+stc1,j,k))
           rhoWi(i,j,k) = sum(interp_Jw_x(i,j,:)*rhoW(i-stc2:i+stc1,j,k))
           
        end do
     end do
  end do
  !$OMP END DO
  
  ! Viscous part
  !$OMP DO
  do kk=kmin_-stv1,kmax_+stv2
     do jj=jmin_-stv1,jmax_+stv2
        do ii=imin_-stv1,imax_+stv2
           
           i = ii-1; j = jj-1; k = kk-1;
           
           FX(i,j,k) = &
                + 2.0_WP*VISC(i,j,k)*( &
                   + sum(grad_u_x(i,j,:)*U(i-stv1:i+stv2,j,k,1)) &
                   - 1.0_WP/3.0_WP*( sum(divv_u(i,j,:)*U(i-stv1:i+stv2,j,k,1)) &
                                   + sum(divv_v(i,j,:)*U(i,j-stv1:j+stv2,k,2)) &
                                   + sum(divv_w(i,j,:)*U(i,j,k-stv1:k+stv2,3))))
           
           i = ii; j = jj; k = kk;
           
           FY(i,j,k) = &
                + sum(interp_sc_xy(i,j,:,:)*VISC(i-st2:i+st1,j-st2:j+st1,k)) * &
                ( sum(grad_u_y(i,j,:)*U(i,j-stv2:j+stv1,k,1)) &
                + sum(grad_v_x(i,j,:)*U(i-stv2:i+stv1,j,k,2)) )
           
           FZ(i,j,k) = &
                + sum(interp_sc_xz(i,j,:,:)*VISC(i-st2:i+st1,j,k-st2:k+st1)) * &
                ( sum(grad_u_z(i,j,:)*U(i,j,k-stv2:k+stv1,1)) &
                + sum(grad_w_x(i,j,:)*U(i-stv2:i+stv1,j,k,3)) )
        end do
     end do
  end do
  !$OMP END DO

  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           rhs = 0.0_WP
           RHOi = sum(interp_sc_x(i,j,:)*RHO(i-st2:i+st1,j,k))

           ! Convective term
           do st=-stc2,stc1
              n = interp_xx(i,j,st)
              rhs = rhs - divc_xx(i,j,st) * rhoUi(i+st,j,k) * &
                   0.5_WP*(U(i+st+n+1,j,k,1)+U(i+st-n,j,k,1))
           end do
           do st=-stc1,stc2              
              n = interp_xy(i,j,st)
              rhs = rhs - divc_xy(i,j,st) * rhoVi(i,j+st,k) * &
                   0.5_WP*(U(i,j+st+n-1,k,1)+U(i,j+st-n,k,1))
              n = interp_xz(i,j,st)
              rhs = rhs - divc_xz(i,j,st) * rhoWi(i,j,k+st) * &
                   0.5_WP*(U(i,j,k+st+n-1,1)+U(i,j,k+st-n,1))
           end do

           tmpxyz(i,j,k) = rhs
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  call phys_avg_z(tmpxyz, U_budget_phys_cf(:,:,1,1))

  ! Enable commented blocks below to test unsteady term
!!$  if (ntime_curr.eq.2) then
!!$     print *, 'second dt'
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           U_budget_phys_cf(i,j,1,1) = tmpxyz(i,j,20)
!!$           !U_budget_phys   (i,j,1,1) = - RHO(i,j,20)*Ui(i,j,20,1)*(dUdx(i,j,20,1,1)+dUdx(i,j,20,2,2)+dUdx(i,j,20,3,3)) &
!!$           !     - RHO(i,j,20)*(Ui(i,j,20,1)*dUdx(i,j,20,1,1) + Ui(i,j,20,2)*dUdx(i,j,20,1,2) + Ui(i,j,20,3)*dUdx(i,j,20,1,3) ) &
!!$           !     - Ui(i,j,20,1)*( Ui(i,j,20,1)*dRHOdx(i,j,20,1) + Ui(i,j,20,2)*dRHOdx(i,j,20,2) + Ui(i,j,20,3)*dRHOdx(i,j,20,3) )
!!$        end do
!!$     end do
!!$  end if

  ! Pressure part of velocity residual
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           tmpxyz(i,j,k) = -sum(grad_x(i,j,:)*P(i-stc2:i+stc1,j,k))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call phys_avg_z(tmpxyz, U_budget_phys_cf(:,:,1,3))

  ! Enable commented blocks below to test unsteady term
!!$  if (ntime_curr.eq.2) then
!!$     print *, 'second dt'
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           U_budget_phys_cf(i,j,1,3) = tmpxyz(i,j,20)
!!$           !U_budget_phys   (i,j,1,3) = -dPdx(i,j,20,1)
!!$        end do
!!$     end do
!!$  end if

  ! Viscous part of velocity residual
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           tmpxyz(i,j,k) =  sum(divv_xx(i,j,:)*FX(i-stv2:i+stv1,j,k)) &
                + sum(divv_xy(i,j,:)*FY(i,j-stv1:j+stv2,k)) &
                + sum(divv_xz(i,j,:)*FZ(i,j,k-stv1:k+stv2)) 
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call phys_avg_z(tmpxyz, U_budget_phys_cf(:,:,1,4))

  ! Enable commented blocks below to test unsteady term
!!$  if (ntime_curr.eq.2) then
!!$     print *, 'second dt'
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           U_budget_phys_cf(i,j,1,4) = tmpxyz(i,j,20)
!!$           !U_budget_phys   (i,j,1,4) = dTAUdx(i,j,20,1,1) + dTAUdx(i,j,20,1,2) + dTAUdx(i,j,20,1,3)
!!$        end do
!!$     end do
!!$  end if

           




  ! y-Momentum ------------------------------------------------------------------------
  ! Convective part of velocity residual
  !$OMP PARALLEL PRIVATE(i,j,k,n,RHOi,st,rhs)
  !$OMP DO
  do kk=kmin_-stc1,kmax_+stc2
     do jj=jmin_-stc1,jmax_+stc2
        do ii=imin_-stc1,imax_+stc2
           
           i = ii-1; j = jj-1; k = kk-1;
           
           rhoVi(i,j,k) = sum(interp_v_ym(i,j,:)*rhoV(i,j-stc1:j+stc2,k))
           
           i = ii; j = jj; k = kk;
           
           rhoUi(i,j,k) = sum(interp_Ju_y(i,j,:)*rhoU(i,j-stc2:j+stc1,k))
           rhoWi(i,j,k) = sum(interp_Jw_y(i,j,:)*rhoW(i,j-stc2:j+stc1,k))
           
        end do
     end do
  end do
  !$OMP END DO
  
  ! Viscous part
  !$OMP DO
  do kk=kmin_-stv1,kmax_+stv2
     do jj=jmin_-stv1,jmax_+stv2
        do ii=imin_-stv1,imax_+stv2
           
           i = ii-1; j = jj-1; k = kk-1;
           
           FY(i,j,k) = &
                + 2.0_WP*VISC(i,j,k)*( &
                   + sum(grad_v_y(i,j,:)*U(i,j-stv1:j+stv2,k,2)) &
                   - 1.0_WP/3.0_WP*( sum(divv_u(i,j,:)*U(i-stv1:i+stv2,j,k,1)) &
                                   + sum(divv_v(i,j,:)*U(i,j-stv1:j+stv2,k,2)) &
                                   + sum(divv_w(i,j,:)*U(i,j,k-stv1:k+stv2,3))))
           
           i = ii; j = jj; k = kk;
           
           FX(i,j,k) = &
                + sum(interp_sc_xy(i,j,:,:)*VISC(i-st2:i+st1,j-st2:j+st1,k)) * &
                ( sum(grad_u_y(i,j,:)*U(i,j-stv2:j+stv1,k,1)) &
                + sum(grad_v_x(i,j,:)*U(i-stv2:i+stv1,j,k,2)) )
           
           FZ(i,j,k) = &
                + sum(interp_sc_yz(i,j,:,:)*VISC(i,j-st2:j+st1,k-st2:k+st1)) * &
                ( sum(grad_v_z(i,j,:)*U(i,j,k-stv2:k+stv1,2)) &
                + sum(grad_w_y(i,j,:)*U(i,j-stv2:j+stv1,k,3)) )
        end do
     end do
  end do
  !$OMP END DO
  
  ! Residual
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           rhs = 0.0_WP
           RHOi = sum(interp_sc_y(i,j,:)*RHO(i,j-st2:j+st1,k))
           
           ! Convective term
           do st=-stc2,stc1
              n = interp_yy(i,j,st)
              rhs = rhs - divc_yy(i,j,st) * rhoVi(i,j+st,k) * &
                   0.5_WP*(U(i,j+st+n+1,k,2)+U(i,j+st-n,k,2))
           end do
           do st=-stc1,stc2
              n = interp_yx(i,j,st)
              rhs = rhs - divc_yx(i,j,st) * rhoUi(i+st,j,k) * &
                   0.5_WP*(U(i+st+n-1,j,k,2)+U(i+st-n,j,k,2))
              n = interp_yz(i,j,st)
              rhs = rhs - divc_yz(i,j,st) * rhoWi(i,j,k+st) * &
                   0.5_WP*(U(i,j,k+st+n-1,2)+U(i,j,k+st-n,2))
           end do

           tmpxyz(i,j,k) = rhs
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  call phys_avg_z(tmpxyz, U_budget_phys_cf(:,:,2,1))

!!$  if (ntime_curr.eq.2) then
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           U_budget_phys_cf(i,j,2,1) = tmpxyz(i,j,20)
!!$           !U_budget_phys   (i,j,2,1) = - RHO(i,j,20)*Ui(i,j,20,2)*(dUdx(i,j,20,1,1)+dUdx(i,j,20,2,2)+dUdx(i,j,20,3,3)) &
!!$           !     - RHO(i,j,20)*(Ui(i,j,20,1)*dUdx(i,j,20,2,1) + Ui(i,j,20,2)*dUdx(i,j,20,2,2) + Ui(i,j,20,3)*dUdx(i,j,20,2,3) ) &
!!$           !     - Ui(i,j,20,2)*( Ui(i,j,20,1)*dRHOdx(i,j,20,1) + Ui(i,j,20,2)*dRHOdx(i,j,20,2) + Ui(i,j,20,3)*dRHOdx(i,j,20,3) )
!!$        end do
!!$     end do
!!$  end if

  ! Pressure part of velocity residual
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           tmpxyz(i,j,k) = -sum(grad_y(i,j,:)*P(i,j-stc2:j+stc1,k))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call phys_avg_z(tmpxyz, U_budget_phys_cf(:,:,2,3))

!!$  if (ntime_curr.eq.2) then
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           U_budget_phys_cf(i,j,2,3) = tmpxyz(i,j,20)
!!$           !U_budget_phys   (i,j,2,3) = -dPdx(i,j,20,2)
!!$        end do
!!$     end do
!!$  end if

  ! Viscous part of velocity residual
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           tmpxyz(i,j,k) = sum(divv_yx(i,j,:)*FX(i-stv1:i+stv2,j,k)) &
                     + sum(divv_yy(i,j,:)*FY(i,j-stv2:j+stv1,k)) &
                     + sum(divv_yz(i,j,:)*FZ(i,j,k-stv1:k+stv2)) 
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call phys_avg_z(tmpxyz, U_budget_phys_cf(:,:,2,4))

!!$  if (ntime_curr.eq.2) then
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           U_budget_phys_cf(i,j,2,4) = tmpxyz(i,j,20)
!!$           !U_budget_phys   (i,j,2,4) = dTAUdx(i,j,20,2,1) + dTAUdx(i,j,20,2,2) + dTAUdx(i,j,20,2,3)
!!$        end do
!!$     end do
!!$  end if


  return
end subroutine volumeStats_velocity_compute
