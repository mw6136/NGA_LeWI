module volumeStats_budgets
  use volumeStats_average

  ! Momentum budgets - physical space -- at cell centers
  !integer, parameter :: nterms_Up=4
  real(WP), dimension(:,:,:,:), pointer :: U_budget_phys

  ! TKE budgets - physical space
  integer, parameter :: nterms_Ep=6
  real(WP), dimension(:,:,:), pointer :: E_budget_phys
  real(WP), dimension(:,:),   pointer :: TKE_phys, eta_phys

  ! Scalar budgets - conditional
  integer, parameter :: nterms_Sc=7
  integer :: ns_out_min, ns_out_max
  real(WP), dimension(:,:,:,:,:), pointer :: SC_budget_cond

  ! Momentum budgets - conditional
  integer, parameter :: nterms_Uc=21
  real(WP), dimension(:,:,:,:,:), pointer :: U_budget_cond

  ! TKE budgets - conditional
  integer, parameter :: nterms_Ec=13
  real(WP), dimension(:,:,:,:),   pointer :: E_budget_cond

end module volumeStats_budgets


subroutine volumeStats_budgets_init
  use volumeStats_budgets
  implicit none

  ! Momentum budgets - physical space
  allocate(U_budget_phys(imin_:imax_,jmin_:jmax_,1:2,0:nterms_Up))
  U_budget_phys = 0.0_WP

  ! TKE budgets - physical space
  allocate(TKE_phys(imin_:imax_,jmin_:jmax_))
  TKE_phys = 0.0_WP
  allocate(E_budget_phys(imin_:imax_,jmin_:jmax_,0:nterms_Ep))
  E_budget_phys = 0.0_WP

  ! Momentum budgets - conditional
  allocate(U_budget_cond(imin_:imax_,jmid:jmax_,1:nbins_cond,1:2,0:nterms_Uc))
  U_budget_cond = 0.0_WP

  ! Scalar budgets - conditional
  ns_out_min = isc_OH
  ns_out_max = isc_H2O
  allocate(SC_budget_cond(imin_:imax_,jmid:jmax_,1:nbins_cond,ns_out_min:ns_out_max,1:nterms_Sc))
  SC_budget_cond = 0.0_WP

  ! TKE budgets - conditional
  allocate(E_budget_cond(imin_:imax_,jmid:jmax_,1:nbins_cond,0:nterms_Ec))
  E_budget_cond = 0.0_WP

  return
end subroutine volumeStats_budgets_init


subroutine volumeStats_budgets_compute
  use volumeStats_budgets
  use volumeStats_metric
  implicit none
  
  integer :: i, j, k, ibin, isc, ii, jj, st
  integer :: iplane, iunit
  character(len=str_short) :: of_name, indx_name
  real(WP), dimension(-st2:st1+1) :: tmpst
  real(WP), dimension(imin_:imax_,jmin_:jmax_) :: tmpxy
  real(WP), dimension(imin_:imax_,jmid:jmax_,1:nbins_cond) :: tmpxyc, tmpxyc_out


  ! Momentum budgets -- physical --------------------------------------------
  do ii=1,2
     ! Unsteady term
     U_budget_phys(:,:,ii,0) = U_budget_phys_cf(:,:,ii,0)

     ! Convective term
     do jj=1,3
        do j=jmin_,jmax_
           do i=imin_,imax_
              U_budget_phys(i,j,ii,1) = U_budget_phys(i,j,ii,1) - &
                   ( Umi(i,j,ii) * Umi(i,j,jj) * dRHOdxm(i,j,jj) &
                   + RHOm(i,j) * Umi(i,j,ii) * dUdxm(i,j,jj,jj) &
                   + RHOm(i,j) * Umi(i,j,jj) * dUdxm(i,j,ii,jj) )
           end do
        end do
     end do

     ! Turbulent stress
     do jj=1,3
        do j=jmin_,jmax_
           do i=imin_,imax_
              U_budget_phys(i,j,ii,2) = U_budget_phys(i,j,ii,2) - &
                   ( dRHOUiUjdxm(i,j,ii,jj) - &
                   ( Umi(i,j,ii)*Umi(i,j,jj)*dRHOdxm(i,j,jj) &
                   + RHOm(i,j) * Umi(i,j,ii) * dUdxm(i,j,jj,jj) &
                   + RHOm(i,j) * Umi(i,j,jj) * dUdxm(i,j,ii,jj) ) )
           end do
        end do
     end do

     ! Pressure gradient
     do j=jmin_,jmax_
        do i=imin_,imax_
           U_budget_phys(i,j,ii,3) = -dPdxm(i,j,ii)
        end do
     end do

     ! Viscous stress
     do jj=1,3
        do j=jmin_,jmax_
           do i=imin_,imax_
              U_budget_phys(i,j,ii,4) = U_budget_phys(i,j,ii,4) + dTAUdxm(i,j,ii,jj)
           end do
        end do
     end do
  end do


  ! Scalar budgets -- conditional ---------------------------------------

  if (combust) then
     do isc=ns_out_min, ns_out_max

        ! Unclosed sample space flux (1)
        !  -- tmp_SC_cz_* terms come from volumeStats_average_all
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 SC_budget_cond(i,j,ibin,isc,1) = 2.0_WP * tmp_SC_cz_1(i,j,ibin,isc) * pdf_Cz(i,j,ibin)
              end do
           end do
        end do
        if (irank.eq.iroot) print *, isc, ' conditional SC (1) done'

        ! Unclosed sample space flux (2-3)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 SC_budget_cond(i,j,ibin,isc,2) = - 0.5_WP * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * tmp_SC_cz_2(i,j,ibin-st2+1:ibin+st1+1,isc) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:) * tmp_SC_cz_2(i,j,ibin-st2  :ibin+st1  ,isc) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        ! Using regularized derivatives
        tmpxyc = tmp_SC_cz_2(imin_:imax_,jmid:jmax_,1:nbins_cond,isc) * pdf_Cz(imin_:imax_,jmid:jmax_,1:nbins_cond)
        if (irank.eq.iroot) print *, 'reg_deriv', isc
        call reg_deriv(tmpxyc, tmpxyc_out, 0)
        SC_budget_cond(imin_:imax_,jmid:jmax_,1:nbins_cond,isc,3) = - 0.5_WP * tmpxyc_out(imin_:imax_,jmid:jmax_,1:nbins_cond)
        if (irank.eq.iroot) print *, isc, ' conditional SC (2-3) done'

        ! Sample space flux with primary closure term 1 (4 and 5)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 SC_budget_cond(i,j,ibin,isc,4) = 0.5_WP * chi_cz(i,j,ibin) * pdf_Cz(i,j,ibin) * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * SCm_cz(i,j,ibin-st2+1:ibin+st1+1,isc) ) &
                      - sum(interp_c(ibin  ,:) * SCm_cz(i,j,ibin-st2  :ibin+st1  ,isc) ) )
              end do
           end do
        end do
        ! Using regularized derivatives
        tmpxyc = SCm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,isc)
        if (irank.eq.iroot) print *, 'reg_deriv', isc
        call reg_deriv(tmpxyc, tmpxyc_out, 0)
        SC_budget_cond(imin_:imax_,jmid:jmax_,1:nbins_cond,isc,5) = 0.5_WP &
             * chi_cz(imin_:imax_,jmid:jmax_,1:nbins_cond) &
             * pdf_Cz(imin_:imax_,jmid:jmax_,1:nbins_cond) &
             * tmpxyc_out(imin_:imax_,jmid:jmax_,1:nbins_cond)
        if (irank.eq.iroot) print *, isc, ' conditional SC (4-5) done'

        ! Sample space flux with primary closure term 2 (6 and 7)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 SC_budget_cond(i,j,ibin,isc,6) = - 0.5_WP * SCm_cz(i,j,ibin,isc) * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * chi_cz(i,j,ibin-st2+1:ibin+st1+1) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:) * chi_cz(i,j,ibin-st2  :ibin+st1  ) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        ! Using regularized derivatives
        tmpxyc = chi_cz(imin_:imax_,jmid:jmax_,1:nbins_cond) * pdf_Cz(imin_:imax_,jmid:jmax_,1:nbins_cond)
        if (irank.eq.iroot) print *, 'reg_deriv', isc
        call reg_deriv(tmpxyc, tmpxyc_out, 0)
        SC_budget_cond(imin_:imax_,jmid:jmax_,1:nbins_cond,isc,7) = -0.5_WP &
             * SCm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,isc) &
             * tmpxyc_out(imin_:imax_,jmid:jmax_,1:nbins_cond)
        if (irank.eq.iroot) print *, isc, ' conditional SC (6-7) done'

     end do
  end if


  ! Momentum budgets -- conditional ---------------------------------------
  if (combust) then
     do ii=1,2
        ! Unsteady term - U_cond_temp includes negative sign
        U_budget_cond(:,:,:,ii,0) = (U_cond_temp(:,:,:,ii) + Um_cz(:,:,:,ii)*RHO_cond_temp(:,:,:))/RHOm_cz(:,:,:)
        !U_budget_cond(:,:,:,ii,0) = dUdt_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii) !U_cond_temp(:,:,:,ii)

        ! Convective term (1)
        do jj=1,3
           do ibin=1,nbins_cond
              do j=jmid,jmax_
                 do i=imin_,imax_
                    U_budget_cond(i,j,ibin,ii,1) = U_budget_cond(i,j,ibin,ii,1) &
                         - Um_cz(i,j,ibin,jj) * dUdxm_cz(i,j,ibin,ii,jj)
                 end do
              end do
           end do
        end do
        do jj=1,2
           do ibin=1,nbins_cond
              do j=jmid,jmax_
                 do i=imin_,imax_
                    U_budget_cond(i,j,ibin,ii,10) = U_budget_cond(i,j,ibin,ii,10) &
                         - Um_cz(i,j,ibin,jj) * dUdxm_cz_reg(i,j,ibin,ii,jj)
                 end do
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (1) done'

        ! Turbulent stress (2)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 ! dk/dz
                 U_budget_cond(i,j,ibin,ii,2) = U_budget_cond(i,j,ibin,ii,2) - dxi(i) * &
                      ( sum(interp_sc_x(i+1,j,:)*RHOm_cz(i-st2+1:i+st1+1,j,ibin) * pdf_Cz(i-st2+1:i+st1+1,j,ibin) * &
                      ( UiUjm_cz(i-st2+1:i+st1+1,j,ibin,ii,1) - Um_cz(i-st2+1:i+st1+1,j,ibin,ii)*Um_cz(i-st2+1:i+st1+1,j,ibin,1) )) &
                      - sum(interp_sc_x(i  ,j,:)*RHOm_cz(i-st2  :i+st1  ,j,ibin) * pdf_Cz(i-st2  :i+st1  ,j,ibin) * &
                      ( UiUjm_cz(i-st2  :i+st1  ,j,ibin,ii,1) - Um_cz(i-st2  :i+st1  ,j,ibin,ii)*Um_cz(i-st2  :i+st1  ,j,ibin,1) )) )
                 ! dk/dy
                 U_budget_cond(i,j,ibin,ii,2) = U_budget_cond(i,j,ibin,ii,2) - dyi(j) * &
                      ( sum(interp_sc_y(i,j+1,:)*RHOm_cz(i,j-st2+1:j+st1+1,ibin) * pdf_Cz(i,j-st2+1:j+st1+1,ibin) * &
                      ( UiUjm_cz(i,j-st2+1:j+st1+1,ibin,ii,1) - Um_cz(i,j-st2+1:j+st1+1,ibin,ii)*Um_cz(i,j-st2+1:j+st1+1,ibin,1) )) &
                      - sum(interp_sc_y(i,j  ,:)*RHOm_cz(i,j-st2  :j+st1  ,ibin) * pdf_Cz(i,j-st2  :j+st1  ,ibin) * &
                      ( UiUjm_cz(i,j-st2  :j+st1  ,ibin,ii,1) - Um_cz(i,j-st2  :j+st1  ,ibin,ii)*Um_cz(i,j-st2  :j+st1  ,ibin,1) )) )
                 ! Complete the weighting
                 U_budget_cond(i,j,ibin,ii,2) = U_budget_cond(i,j,ibin,ii,2) / (RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin))
              end do
           end do
        end do
        ! Using regularized derivatives
        do jj=1,2
           tmpxyc = RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond) * pdf_Cz(imin_:imax_,jmid:jmax_,1:nbins_cond) * &
                ( UiUjm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,jj) &
                - Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii) * Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,jj) )
           if (irank.eq.iroot) print *, 'reg_deriv', ii, jj
           call reg_deriv(tmpxyc, tmpxyc_out, jj)
           U_budget_cond(:,:,:,ii,11) = U_budget_cond(:,:,:,ii,11) &
                - tmpxyc_out / (RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond)*pdf_Cz(imin_:imax_,jmid:jmax_,1:nbins_cond))
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (2) done'


        ! Diffusive transport in manifold space (3)
        !  -- tmp_U_cz_* terms come from volumeStats_average_all
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 ! d/dx
                 do st=-st2,st1+1
                    tmpst(st) = dci(ibin) * &  ! d/d lambda
                         ( sum(interp_c(ibin+1,:)*tmp_U_cz_1(i+st,j,ibin-st2+1:ibin+st1+1,ii,1) * pdf_Cz(i+st,j,ibin-st2+1:ibin+st1+1)) &
                         - sum(interp_c(ibin  ,:)*tmp_U_cz_1(i+st,j,ibin-st2  :ibin+st1  ,ii,1) * pdf_Cz(i+st,j,ibin-st2  :ibin+st1  )) )
                 end do
                 U_budget_cond(i,j,ibin,ii,3) = U_budget_cond(i,j,ibin,ii,3) - dxi(i) * &
                      ( sum(interp_sc_x(i+1,j,:)*tmpst(-st2+1:st1+1)) &
                      - sum(interp_sc_x(i  ,j,:)*tmpst(-st2  :st1  )) )
                 ! d/dy
                 do st=-st2,st1+1
                    tmpst(st) = dci(ibin) * &  ! d/d lambda
                         ( sum(interp_c(ibin+1,:)*tmp_U_cz_1(i,j+st,ibin-st2+1:ibin+st1+1,ii,2) * pdf_Cz(i,j+st,ibin-st2+1:ibin+st1+1)) &
                         - sum(interp_c(ibin  ,:)*tmp_U_cz_1(i,j+st,ibin-st2  :ibin+st1  ,ii,2) * pdf_Cz(i,j+st,ibin-st2  :ibin+st1  )) )
                 end do
                 U_budget_cond(i,j,ibin,ii,3) = U_budget_cond(i,j,ibin,ii,3) - dyi(j) * &
                      ( sum(interp_sc_y(i,j+1,:)*tmpst(-st2+1:st1+1)) &
                      - sum(interp_sc_y(i,j  ,:)*tmpst(-st2  :st1  )) )
                 ! Complete the weighting
                 U_budget_cond(i,j,ibin,ii,3) = U_budget_cond(i,j,ibin,ii,3) / (RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin))
              end do
           end do
        end do
        ! Using regularized derivatives
        do jj=1,2
           tmpxyc = 0.0_WP
           do ibin=1,nbins_cond
              do j=jmid,jmax_
                 do i=imin_,imax_
                    tmpxyc(i,j,ibin) = dci(ibin) * &  ! d/d lambda
                         ( sum(interp_c(ibin+1,:)*tmp_U_cz_1(i,j,ibin-st2+1:ibin+st1+1,ii,jj) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                         - sum(interp_c(ibin  ,:)*tmp_U_cz_1(i,j,ibin-st2  :ibin+st1  ,ii,jj) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
                 end do
              end do
           end do
           if (irank.eq.iroot) print *,'reg_deriv', ii, jj 
           call reg_deriv(tmpxyc, tmpxyc_out, jj)
           U_budget_cond(:,:,:,ii,12) = U_budget_cond(:,:,:,ii,12) &
                - tmpxyc_out / (RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond)*pdf_Cz(imin_:imax_,jmid:jmax_,1:nbins_cond))
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (3) done'

        ! Diffusion of progress variable in manifold space (4)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 ! d/dx
                 do st=-st2,st1+1
                    tmpst(st) = dci(ibin) * &  ! d/d lambda
                         ( sum(interp_c(ibin+1,:) * tmp_U_cz_2(i+st,j,ibin-st2+1:ibin+st1+1,1) * pdf_Cz(i+st,j,ibin-st2+1:ibin+st1+1)) &
                         - sum(interp_c(ibin  ,:) * tmp_U_cz_2(i+st,j,ibin-st2  :ibin+st1  ,1) * pdf_Cz(i+st,j,ibin-st2  :ibin+st1  )) )
                 end do
                 U_budget_cond(i,j,ibin,ii,4) = U_budget_cond(i,j,ibin,ii,4) + dxi(i) * &
                      ( sum(interp_sc_x(i+1,j,:)*tmpst(-st2+1:st1+1)) &
                      - sum(interp_sc_x(i  ,j,:)*tmpst(-st2  :st1  )) )
                 ! d/dy
                 do st=-st2,st1+1
                    tmpst(st) = dci(ibin) * &  ! d/d lambda
                         ( sum(interp_c(ibin+1,:) * tmp_U_cz_2(i,j+st,ibin-st2+1:ibin+st1+1,2) * pdf_Cz(i,j+st,ibin-st2+1:ibin+st1+1)) &
                         - sum(interp_c(ibin  ,:) * tmp_U_cz_2(i,j+st,ibin-st2  :ibin+st1  ,2) * pdf_Cz(i,j+st,ibin-st2  :ibin+st1  )) )
                 end do
                 U_budget_cond(i,j,ibin,ii,4) = U_budget_cond(i,j,ibin,ii,4) + dyi(j) * &
                      ( sum(interp_sc_y(i,j+1,:)*tmpst(-st2+1:st1+1)) &
                      - sum(interp_sc_y(i,j  ,:)*tmpst(-st2  :st1  )) )
                 ! Complete the weighting
                 U_budget_cond(i,j,ibin,ii,4) = U_budget_cond(i,j,ibin,ii,4) * &
                      Um_cz(i,j,ibin,ii) / (RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin))
              end do
           end do
        end do
        ! Using regularized derivatives
        do jj=1,2
           tmpxyc = 0.0_WP
           do ibin=1,nbins_cond
              do j=jmid,jmax_
                 do i=imin_,imax_
                    tmpxyc(i,j,ibin) = dci(ibin) * &  ! d/d lambda
                         ( sum(interp_c(ibin+1,:)*tmp_U_cz_2(i,j,ibin-st2+1:ibin+st1+1,jj) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                         - sum(interp_c(ibin  ,:)*tmp_U_cz_2(i,j,ibin-st2  :ibin+st1  ,jj) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
                 end do
              end do
           end do
           if (irank.eq.iroot) print *, 'reg_deriv', ii, jj
           call reg_deriv(tmpxyc, tmpxyc_out, jj)
           U_budget_cond(:,:,:,ii,13) = U_budget_cond(:,:,:,ii,13) + Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii) * tmpxyc_out / &
                (RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond)*pdf_Cz(imin_:imax_,jmid:jmax_,1:nbins_cond))
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (4) done'

        ! Pressure gradient (5)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 U_budget_cond(i,j,ibin,ii,5) = -dPdxm_cz(i,j,ibin,ii)/RHOm_cz(i,j,ibin)
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (5) done'

        ! Stress tensor (6)
        do jj=1,3
           do ibin=1,nbins_cond
              do j=jmid,jmax_
                 do i=imin_,imax_
                    U_budget_cond(i,j,ibin,ii,6) = dTAUdxm_cz(i,j,ibin,ii,jj)/RHOm_cz(i,j,ibin)
                 end do
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (6) done'

        ! Conditional gradient of velocity-prog. source term correlation (7)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 U_budget_cond(i,j,ibin,ii,7) = -1.0_WP/(RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin)) * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * tmp_U_cz_3(i,j,ibin-st2+1:ibin+st1+1,ii) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:) * tmp_U_cz_3(i,j,ibin-st2  :ibin+st1  ,ii) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (7) done'

        ! PDF equation prog. source term (8)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 U_budget_cond(i,j,ibin,ii,8) = Um_cz(i,j,ibin,ii) / (RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin)) * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * w_dot_cz(i,j,ibin-st2+1:ibin+st1+1) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:) * w_dot_cz(i,j,ibin-st2  :ibin+st1  ) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (8) done'

        ! Conditional gradient of prog. var. dissipation rate (9)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 U_budget_cond(i,j,ibin,ii,9) =  -0.5_WP/(RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin)) * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * chi_cz(i,j,ibin-st2+1:ibin+st1+1) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:) * chi_cz(i,j,ibin-st2  :ibin+st1  ) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) ) &
                      * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * Um_cz(i,j,ibin-st2+1:ibin+st1+1,ii) ) &
                      - sum(interp_c(ibin  ,:) * Um_cz(i,j,ibin-st2  :ibin+st1  ,ii) ) )
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (9) done'

        ! Unclosed sample space flux - term 1 (14)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 U_budget_cond(i,j,ibin,ii,14) = 2.0_WP * tmp_U_cz_4(i,j,ibin,ii) * pdf_Cz(i,j,ibin)
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (14) done'
        ! Take derivative of (14) wrt. lambda for term in budget equation (20)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 U_budget_cond(i,j,ibin,ii,20) = 2.0_WP / (RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin)) * dci(ibin)* &
                      ( sum(interp_c(ibin+1,:) * tmp_U_cz_4(i,j,ibin-st2+1:ibin+st1+1,ii) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:) * tmp_U_cz_4(i,j,ibin-st2  :ibin+st1  ,ii) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (20) done'
        

        ! Unclosed sample space flux - term 2 (15)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 tmpxyc(i,j,ibin) = - 0.5_WP * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * tmp_U_cz_5(i,j,ibin-st2+1:ibin+st1+1,ii) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:) * tmp_U_cz_5(i,j,ibin-st2  :ibin+st1  ,ii) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        U_budget_cond(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,15) = tmpxyc
        if (irank.eq.iroot) print *, ii, ' conditional U (15) done'
        ! Take derivative of (15) wrt. lambda for term in budget equation (21)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 U_budget_cond(i,j,ibin,ii,21) = dci(ibin) / (RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin)) * &
                      ( sum(interp_c(ibin+1,:) * tmpxyc(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:) * tmpxyc(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (21) done'

        ! Sample space flux with primary closure term 1 (16)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 U_budget_cond(i,j,ibin,ii,16) = 0.5_WP * chi_cz(i,j,ibin) * pdf_Cz(i,j,ibin) * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * Um_cz(i,j,ibin-st2+1:ibin+st1+1,ii) ) &
                      - sum(interp_c(ibin  ,:) * Um_cz(i,j,ibin-st2  :ibin+st1  ,ii) ) )
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (16) done'

        ! Sample space flux with primary closure term 2 (17)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 U_budget_cond(i,j,ibin,ii,17) = - 0.5_WP * Um_cz(i,j,ibin,ii) * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * chi_cz(i,j,ibin-st2+1:ibin+st1+1) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:) * chi_cz(i,j,ibin-st2  :ibin+st1  ) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (17) done'

        ! Conditional flux of velocity-prog. source term correlation (18)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 U_budget_cond(i,j,ibin,ii,18) = -1.0_WP/(RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin)) * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * (tmp_U_cz_3(i,j,ibin-st2+1:ibin+st1+1,ii) - w_dot_cz(i,j,ibin-st2+1:ibin+st1+1)*Um_cz(i,j,ibin-st2+1:ibin+st1+1,ii)) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:) * (tmp_U_cz_3(i,j,ibin-st2  :ibin+st1  ,ii) - w_dot_cz(i,j,ibin-st2  :ibin+st1  )*Um_cz(i,j,ibin-st2  :ibin+st1  ,ii)) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (18) done'

        ! Conditional flux of velocity-prog. source term correlation (19)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 U_budget_cond(i,j,ibin,ii,19) = -w_dot_cz(i,j,ibin)/RHOm_cz(i,j,ibin) * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * Um_cz(i,j,ibin-st2+1:ibin+st1+1,ii)) &
                      - sum(interp_c(ibin  ,:) * Um_cz(i,j,ibin-st2  :ibin+st1  ,ii)) )
              end do
           end do
        end do
        if (irank.eq.iroot) print *, ii, ' conditional U (19) done'

     end do
  end if

  ! TKE -------------------------------------------------------------------
  do jj=1,3
     do j=jmin_,jmax_
        do i=imin_,imax_
           TKE_phys(i,j) = TKE_phys(i,j) + R_stress(i,j,jj,jj)
        end do
     end do
  end do

  ! TKE budgets -- physical --------------------------------------------------
  ! Unsteady term
  do j=jmin_,jmax_
     do i=imin_,imax_
        E_budget_phys(i,j,0) = E_temp_m(i,j)
     end do
  end do
  do ii=1,2
     do j=jmin_,jmax_
        do i=imin_,imax_
           E_budget_phys(i,j,0) = E_budget_phys(i,j,0) - Umi(i,j,ii)*U_budget_phys(i,j,ii,0) - RHOm(i,j)*Umi(i,j,ii)*U_temp_m(i,j,ii)
           !E_budget_phys(i,j,0) = E_budget_phys(i,j,0) - 2.0_WP*Umi(i,j,ii)*U_budget_phys(i,j,ii,0)
        end do
     end do
  end do
  
  ! Convective term
  do ii=1,3
     do j=jmin_,jmax_
        do i=imin_,imax_
           E_budget_phys(i,j,1) = E_budget_phys(i,j,1) - ( RHOm(i,j)*(UiUjm(i,j,ii,ii) - Umi(i,j,ii)*Umi(i,j,ii)) )
        end do
     end do
  end do
  tmpxy = 0.0_WP
  do jj=1,3
     do j=jmin_,jmax_
        do i=imin_,imax_
           tmpxy(i,j) = tmpxy(i,j) + dUdxm(i,j,jj,jj)
        end do
     end do
  end do
  do j=jmin_,jmax_
     do i=imin_,imax_
        E_budget_phys(i,j,1) = E_budget_phys(i,j,1) * tmpxy(i,j)
     end do
  end do
  do jj=1,3
     tmpxy = 0.0_WP
     do ii=1,3
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmpxy(i,j) = tmpxy(i,j) + Umi(i,j,ii)*Umi(i,j,ii)*dRHOdxm(i,j,jj) + 2.0_WP*RHOm(i,j)*Umi(i,j,ii)*dUdxm(i,j,ii,jj)
           end do
        end do
     end do
     do j=jmin_,jmax_
        do i=imin_,imax_
           E_budget_phys(i,j,1) = E_budget_phys(i,j,1) - Umi(i,j,jj)*(dRHOUiUidxm(i,j,jj) - tmpxy(i,j))
        end do
     end do
  end do
     

!!$  do jj=1,3
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           E_budget_phys(i,j,1) = E_budget_phys(i,j,1) - &
!!$                ( dRHOdxm(i,j,jj)*Umi(i,j,jj) + RHOm(i,j)*dUdxm(i,j,jj,jj) )*TKE_phys(i,j) 
!!$        end do
!!$     end do
!!$  end do
!!$  do j=jmin_,jmax_
!!$     do i=imin_,imax_
!!$        ! dk/dx
!!$        E_budget_phys(i,j,1) = E_budget_phys(i,j,1) - RHOm(i,j)*Umi(i,j,1)*dxi(i) * &
!!$             ( sum(interp_sc_x(i+1,j,:)*TKE_phys(i-st2+1:i+st1+1,j)) &
!!$             - sum(interp_sc_x(i  ,j,:)*TKE_phys(i-st2  :i+st1  ,j)) )
!!$        ! dk/dy
!!$        E_budget_phys(i,j,1) = E_budget_phys(i,j,1) - RHOm(i,j)*Umi(i,j,2)*dyi(j) * &
!!$             ( sum(interp_sc_y(i,j+1,:)*TKE_phys(i,j-st2+1:j+st1+1)) &
!!$             - sum(interp_sc_y(i,j  ,:)*TKE_phys(i,j-st2  :j+st1  )) )
!!$     end do
!!$  end do

  ! Production term
  do jj=1,3
     do ii=1,3
        do j=jmin_,jmax_
           do i=imin_,imax_
              E_budget_phys(i,j,2) = E_budget_phys(i,j,2) + &
                   R_stress(i,j,ii,jj) * dUdxm(i,j,ii,jj)
           end do
        end do
     end do
  end do
  do j=jmin_,jmax_
     do i=imin_,imax_
        E_budget_phys(i,j,2) = -2.0_WP * RHOm(i,j) * E_budget_phys(i,j,2)
     end do
  end do

  ! Turbulent transport term
  tmpxy = 0.0_WP
  do ii=1,3
     do j=jmin_,jmax_
        do i=imin_,imax_
           tmpxy(i,j) = tmpxy(i,j) + UiUjm(i,j,ii,ii)
        end do
     end do
  end do
  do jj=1,3
     do j=jmin_,jmax_
        do i=imin_,imax_
           E_budget_phys(i,j,3) = E_budget_phys(i,j,3) - ( &
                + dRHOUiUiUjdxm(i,j,jj) - ( dRHOUiUidxm(i,j,jj)*Umi(i,j,jj) + RHOm(i,j)*tmpxy(i,j)*dUdxm(i,j,jj,jj) ) )
        end do
     end do
  end do
  do ii=1,3
     do jj=1,3
        do j=jmin_,jmax_
           do i=imin_,imax_
              E_budget_phys(i,j,3) = E_budget_phys(i,j,3) - ( &
                   - 2.0_WP*(dRHOUiUjdxm(i,j,ii,jj)*Umi(i,j,ii) + RHOm(i,j)*UiUjm(i,j,ii,jj)*dUdxm(i,j,ii,jj)) &
                   + 2.0_WP*( 2.0_WP*RHOm(i,j)*Umi(i,j,ii)*Umi(i,j,jj)*dUdxm(i,j,ii,jj) &
                   + RHOm(i,j)*Umi(i,j,ii)*Umi(i,j,ii)*dUdxm(i,j,jj,jj) &
                   + Umi(i,j,ii)*Umi(i,j,ii)*Umi(i,j,jj)*dRHOdxm(i,j,jj) ) )
           end do
        end do
     end do
  end do


!!$  do ii=1,3
!!$     do j=jmin_,jmax_
!!$        do i=imin_,imax_
!!$           ! d/dx
!!$           E_budget_phys(i,j,3) = E_budget_phys(i,j,3) - dxi(i) * &
!!$                ( sum(interp_sc_x(i+1,j,:)*RHOm(i-st2+1:i+st1+1,j)*vel_triple_corr(i-st2+1:i+st1+1,j,ii,1)) &
!!$                - sum(interp_sc_x(i  ,j,:)*RHOm(i-st2  :i+st1  ,j)*vel_triple_corr(i-st2  :i+st1  ,j,ii,1)) )
!!$           ! d/dy
!!$           E_budget_phys(i,j,3) = E_budget_phys(i,j,3) - dyi(j) * &
!!$                ( sum(interp_sc_y(i,j+1,:)*RHOm(i,j-st2+1:j+st1+1)*vel_triple_corr(i,j-st2+1:j+st1+1,ii,2)) &
!!$                - sum(interp_sc_y(i,  j,:)*RHOm(i,j-st2  :j+st1  )*vel_triple_corr(i,j-st2  :j+st1  ,ii,2)) )
!!$        end do
!!$     end do
!!$  end do

  ! Velocity-pressure gradient correlation
  do j=jmin_,jmax_
     do i=imin_,imax_
        E_budget_phys(i,j,4) = E_budget_phys(i,j,4) - 2.0_WP * &
             ( tmp_E_phys(i,j,1) - sum(Umi(i,j,:) * dPdxm(i,j,:)) )
     end do
  end do

  ! Viscous transport
  do j=jmin_,jmax_
     do i=imin_,imax_
        E_budget_phys(i,j,5) = 2.0_WP*tmp_E_phys(i,j,2)
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=jmin_,jmax_
           do i=imin_,imax_
              E_budget_phys(i,j,5) = E_budget_phys(i,j,5) - 2.0_WP * &
                   ( dUdxm(i,j,ii,jj)*TAUm(i,j,ii,jj) + Umi(i,j,ii)*dTAUdxm(i,j,ii,jj) )
           end do
        end do
     end do
  end do
  
  ! Dissipation
  do j=jmin_,jmax_
     do i=imin_,imax_
        E_budget_phys(i,j,6) = E_budget_phys(i,j,6) - 2.0_WP*tmp_E_phys(i,j,3)
     end do
  end do
  do jj=1,3
     do ii=1,3
        do j=jmin_,jmax_
           do i=imin_,imax_
              E_budget_phys(i,j,6) = E_budget_phys(i,j,6) + 2.0_WP * &
                   ( TAUm(i,j,ii,jj)*dUdxm(i,j,ii,jj) )
           end do
        end do
     end do
  end do


  ! TKE budgets -- conditional ----------------------------------------------
  if (combust) then
     ! Unsteady term (0) - E_cond_temp includes negative sign
     tmpxyc = 0.0_WP
     do ii=1,3
        tmpxyc = tmpxyc + UiUjm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,ii)
     end do
     E_budget_cond(:,:,:,0) = (E_cond_temp + tmpxyc*RHO_cond_temp)/RHOm_cz
     tmpxyc = 0.0_WP
     do ii=1,3
        tmpxyc = tmpxyc &
             + Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii) * &
             ( U_cond_temp(imin_:imax_,jmid:jmax_,1:nbins_cond,ii) &
             + Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii)*RHO_cond_temp(imin_:imax_,jmid:jmax_,1:nbins_cond) ) &
             / RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond)
!!$             * U_budget_cond(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,0)
     end do
     E_budget_cond(:,:,:,0) = E_budget_cond(:,:,:,0) - 2.0_WP*tmpxyc ! negative because U_budget_cond(0) has negative sign included

     ! Convection term (1)
     tmpxyc = 0.0_WP
     do ii=1,3
        tmpxyc = tmpxyc &
             + UiUjm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,ii) &
             - Um_cz   (imin_:imax_,jmid:jmax_,1:nbins_cond,ii) &
             * Um_cz   (imin_:imax_,jmid:jmax_,1:nbins_cond,ii)
     end do
     do jj=1,2
        call reg_deriv(tmpxyc, tmpxyc_out, jj)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 E_budget_cond(i,j,ibin,1) = E_budget_cond(i,j,ibin,1) &
                      - Um_cz(i,j,ibin,jj) * tmpxyc_out(i,j,ibin)
              end do
           end do
        end do
     end do
     if (irank.eq.iroot) print *, 'conditional TKE (1) done'

     ! Production term (2)
     do ii=1,3
        do jj=1,2
           do ibin=1,nbins_cond
              do j=jmid,jmax_
                 do i=imin_,imax_
                    E_budget_cond(i,j,ibin,2) = E_budget_cond(i,j,ibin,2) - 2.0_WP *  &
                         ( UiUjm_cz(i,j,ibin,ii,jj) - Um_cz(i,j,ibin,ii) * Um_cz(i,j,ibin,jj) ) &
                         * dUdxm_cz_reg(i,j,ibin,ii,jj)
                 end do
              end do
           end do
        end do
     end do
     if (irank.eq.iroot) print *, 'conditional TKE (2) done'

     ! Turbulent transport (3)
     do jj=1,2
        tmpxyc = 0.0_WP
        do ii=1,3
           tmpxyc = tmpxyc &
                + UiUiUjm_cz     (imin_:imax_,jmid:jmax_,1:nbins_cond,ii,jj) &
                - UiUjm_cz       (imin_:imax_,jmid:jmax_,1:nbins_cond,ii,ii)*Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,jj) &
                - 2.0_WP*UiUjm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,jj)*Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii) &
                + 2.0_WP*Um_cz   (imin_:imax_,jmid:jmax_,1:nbins_cond,ii)*Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii)*Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,jj)
        end do
        tmpxyc = tmpxyc &
             * RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond) &
             * pdf_Cz (imin_:imax_,jmid:jmax_,1:nbins_cond)
        if (irank.eq.iroot) print *, 'reg_deriv', jj
        call reg_deriv(tmpxyc, tmpxyc_out, jj)
        E_budget_cond(:,:,:,3) = E_budget_cond(:,:,:,3) &
             - tmpxyc_out / (RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond)*pdf_Cz(imin_:imax_,jmid:jmax_,1:nbins_cond))
     end do
     if (irank.eq.iroot) print *, 'conditional TKE (3) done'


     ! Diffusive transport in manifold space (4)
     ! Using regularized derivatives
     do jj=1,2
        tmpxyc = 0.0_WP
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 tmpxyc(i,j,ibin) = dci(ibin) * &  ! d/d lambda
                      ( sum(interp_c(ibin+1,:)*tmp_E_cz_1(i,j,ibin-st2+1:ibin+st1+1,jj) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:)*tmp_E_cz_1(i,j,ibin-st2  :ibin+st1  ,jj) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        if (irank.eq.iroot) print *,'reg_deriv', jj 
        call reg_deriv(tmpxyc, tmpxyc_out, jj)
        E_budget_cond(:,:,:,4) = E_budget_cond(:,:,:,4) &
             - tmpxyc_out / (RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond)*pdf_Cz(imin_:imax_,jmid:jmax_,1:nbins_cond))
     end do
     if (irank.eq.iroot) print *, 'conditional TKE (4) done'


     ! Diffusion of progress variable in manifold space (5)
     ! Using regularized derivatives
     do ii=1,3
        do jj=1,2
           tmpxyc = 0.0_WP
           do ibin=1,nbins_cond
              do j=jmid,jmax_
                 do i=imin_,imax_
                    tmpxyc(i,j,ibin) = dci(ibin) * &  ! d/d lambda
                         ( sum(interp_c(ibin+1,:)*tmp_U_cz_1(i,j,ibin-st2+1:ibin+st1+1,ii,jj) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                         - sum(interp_c(ibin  ,:)*tmp_U_cz_1(i,j,ibin-st2  :ibin+st1  ,ii,jj) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
                 end do
              end do
           end do
           if (irank.eq.iroot) print *, 'reg_deriv', ii, jj
           call reg_deriv(tmpxyc, tmpxyc_out, jj)
           E_budget_cond(:,:,:,5) = E_budget_cond(:,:,:,5) + 2.0_WP * Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii) * tmpxyc_out / &
                (RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond)*pdf_Cz(imin_:imax_,jmid:jmax_,1:nbins_cond))
        end do
     end do
     if (irank.eq.iroot) print *, 'conditional TKE (5) done'


     ! Diffusion of progress variable in manifold space (6)
     ! Need better names
     do jj=1,2
        tmpxyc = 0.0_WP
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 tmpxyc(i,j,ibin) = dci(ibin) * &  ! d/d lambda
                      ( sum(interp_c(ibin+1,:)*tmp_U_cz_2(i,j,ibin-st2+1:ibin+st1+1,jj) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:)*tmp_U_cz_2(i,j,ibin-st2  :ibin+st1  ,jj) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        if (irank.eq.iroot) print *, 'reg_deriv', jj
        call reg_deriv(tmpxyc, tmpxyc_out, jj)
        E_budget_cond(:,:,:,6) = E_budget_cond(:,:,:,6) + tmpxyc_out
     end do
     tmpxyc = 0.0_WP
     do ii=1,3
        tmpxyc = tmpxyc + &
             ( UiUjm_cz    (imin_:imax_,jmid:jmax_,1:nbins_cond,ii,ii) &
             - 2.0_WP*Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii)**2 )
     end do
     E_budget_cond(:,:,:,6) = E_budget_cond(:,:,:,6) * &
          tmpxyc / (RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond)*pdf_Cz(imin_:imax_,jmid:jmax_,1:nbins_cond))
     if (irank.eq.iroot) print *, 'conditional TKE (6) done'


     ! Conditional gradient of TKE-prog. source term correlation (7)
     do ibin=1,nbins_cond
        do j=jmid,jmax_
           do i=imin_,imax_
              E_budget_cond(i,j,ibin,7) = -1.0_WP/(RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin)) * dci(ibin) * &
                   ( sum(interp_c(ibin+1,:) * tmp_E_cz_2(i,j,ibin-st2+1:ibin+st1+1) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                   - sum(interp_c(ibin  ,:) * tmp_E_cz_2(i,j,ibin-st2  :ibin+st1  ) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
           end do
        end do
     end do
     if (irank.eq.iroot) print *, 'conditional TKE (7) done'


     ! Conditional gradient of vel-prog. source term correlation (8)
     do ii=1,3
        tmpxyc = 0.0_WP
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 tmpxyc(i,j,ibin) = 1.0_WP/(RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin)) * dci(ibin) * &
                      ( sum(interp_c(ibin+1,:) * tmp_U_cz_3(i,j,ibin-st2+1:ibin+st1+1,ii) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                      - sum(interp_c(ibin  ,:) * tmp_U_cz_3(i,j,ibin-st2  :ibin+st1  ,ii) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
              end do
           end do
        end do
        E_budget_cond(:,:,:,8) = E_budget_cond(:,:,:,8) &
             + 2.0_WP * Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii) * tmpxyc
     end do
     if (irank.eq.iroot) print *, 'conditional TKE (8) done'


     ! PDF equation prog. source term (9)
     do ibin=1,nbins_cond
        do j=jmid,jmax_
           do i=imin_,imax_
              E_budget_cond(i,j,ibin,9) = 1.0_WP/(RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin)) * dci(ibin) * &
                   ( sum(interp_c(ibin+1,:) * w_dot_cz(i,j,ibin-st2+1:ibin+st1+1) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                   - sum(interp_c(ibin  ,:) * w_dot_cz(i,j,ibin-st2  :ibin+st1  ) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) )
           end do
        end do
     end do
     tmpxyc = 0.0_WP
     do ii=1,3
        tmpxyc = tmpxyc + &
             ( UiUjm_cz    (imin_:imax_,jmid:jmax_,1:nbins_cond,ii,ii) &
             - 2.0_WP*Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii)**2 )
     end do
     E_budget_cond(:,:,:,9) = E_budget_cond(:,:,:,9) * tmpxyc
     if (irank.eq.iroot) print *, 'conditional TKE (9) done'


     ! Conditional gradient of prog. var. dissipation rate (10)
     tmpxyc = 0.0_WP
     do ii=1,3
        tmpxyc = tmpxyc + &
             ( UiUjm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,ii) &
             - Um_cz   (imin_:imax_,jmid:jmax_,1:nbins_cond,ii)**2 )
     end do
     do ibin=1,nbins_cond
        do j=jmid,jmax_
           do i=imin_,imax_
              E_budget_cond(i,j,ibin,10) =  -0.5_WP/(RHOm_cz(i,j,ibin)*pdf_Cz(i,j,ibin)) * dci(ibin) * &
                   ( sum(interp_c(ibin+1,:) * chi_cz(i,j,ibin-st2+1:ibin+st1+1) * pdf_Cz(i,j,ibin-st2+1:ibin+st1+1)) &
                   - sum(interp_c(ibin  ,:) * chi_cz(i,j,ibin-st2  :ibin+st1  ) * pdf_Cz(i,j,ibin-st2  :ibin+st1  )) ) * dci(ibin) * &
                   ( sum(interp_c(ibin+1,:) * tmpxyc(i,j,ibin-st2+1:ibin+st1+1) ) &
                   - sum(interp_c(ibin  ,:) * tmpxyc(i,j,ibin-st2  :ibin+st1  ) ) )
           end do
        end do
     end do
     if (irank.eq.iroot) print *, 'conditional TKE (10) done'


     ! Velocity-pressure gradient (11)
     tmpxyc = 0.0_WP
     do ii=1,3
        tmpxyc = tmpxyc + &
             ( Um_cz   (imin_:imax_,jmid:jmax_,1:nbins_cond,ii) &
             * dPdxm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii) )
     end do
     E_budget_cond(:,:,:,11) = &
          - 2.0_WP/RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond) * &
          ( tmp_E_cz_3 - tmpxyc )
     if (irank.eq.iroot) print *, 'conditional TKE (11) done'


     ! Viscous transport (12)


     ! Viscous dissipation (13)

     ! Combined viscous effects (12a)
     tmpxyc = 0.0_WP
     do ii=1,3
        do jj=1,3
           tmpxyc = tmpxyc + &
                ( Um_cz     (imin_:imax_,jmid:jmax_,1:nbins_cond,ii) &
                * dTAUdxm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,jj) )
        end do
     end do
     E_budget_cond(:,:,:,12) = &
          + 2.0_WP/RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond) * &
          ( tmp_E_cz_4 - tmpxyc )
     if (irank.eq.iroot) print *, 'conditional TKE (12a) done'


  end if


  ! Outputs -----------------------------------------------------------------
  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_budgets_U_cf_'//trim(of_name), action='write')
     i = pnindx(iplane)
     do j=jmin_,jmax_
        write (iunit,'(ES22.13)',advance='no'), y(j)
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*y_condm(i-1:i,j))
        do ii=1,3
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*Um(i-1:i,j,ii))
        end do
        do ii=1,2
           do jj=0,nterms_Up
              if (ii.eq.2.and.jj.eq.nterms_Up) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*U_budget_phys_cf(i-1:i,j,ii,jj))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*U_budget_phys_cf(i-1:i,j,ii,jj))
              end if
           end do
        end do
     end do
     close(iunit)
  end do

  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_budgets_U_'//trim(of_name), action='write')
     i = pnindx(iplane)
     do j=jmin_,jmax_
        write (iunit,'(ES22.13)',advance='no'), y(j)
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*y_condm(i-1:i,j))
        do ii=1,3
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*Umi(i-1:i,j,ii))
        end do
        do ii=1,2
           do jj=0,nterms_Up
              if (ii.eq.2.and.jj.eq.nterms_Up) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*U_budget_phys(i-1:i,j,ii,jj))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*U_budget_phys(i-1:i,j,ii,jj))
              end if
           end do
        end do
     end do
     close(iunit)
  end do

  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_budgets_TKE_'//trim(of_name), action='write')
     i = pnindx(iplane)
     do j=jmin_,jmax_
        write (iunit,'(ES22.13)',advance='no'), y(j)
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*y_condm(i-1:i,j))
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*TKE_phys(i-1:i,j))
        do ii=0,nterms_Ep
           if (ii.eq.nterms_Ep) then
              write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*E_budget_phys(i-1:i,j,ii))
           else
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*E_budget_phys(i-1:i,j,ii))
           end if
        end do
     end do
     close(iunit)
  end do

  if (combust) then
     do jj=1,nterms_Sc
        do isc=ns_out_min,ns_out_max
           write(indx_name,'(I2.2)'), jj
           do iplane=pnmin_,pnmax_
              call get_name(iplane,of_name)
              open(unit=iunit, file=trim(output_name)//'_budgets_SC_cond_'//trim(SC_name(isc))//'_'//trim(indx_name)//'_'//trim(of_name), action='write')
              i = pnindx(iplane)
              do ibin=1,nbins_cond
                 write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
                 do j=jmid,jmax_
                    if (j.eq.jmax_) then
                       write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*SC_budget_cond(i-1:i,j,ibin,isc,jj))
                    else
                       write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SC_budget_cond(i-1:i,j,ibin,isc,jj))
                    end if
                 end do
              end do
              close(iunit)
           end do
        end do
     end do

     do jj=0,nterms_Uc
        do ii=1,2
           write(indx_name,'(I1,I2.2)'), ii,jj
           do iplane=pnmin_,pnmax_
              call get_name(iplane,of_name)
              open(unit=iunit, file=trim(output_name)//'_budgets_U_cond_'//trim(indx_name)//'_'//trim(of_name), action='write')
              i = pnindx(iplane)
              do ibin=1,nbins_cond
                 write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
                 do j=jmid,jmax_
                    if (j.eq.jmax_) then
                       write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*U_budget_cond(i-1:i,j,ibin,ii,jj))
                    else
                       write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*U_budget_cond(i-1:i,j,ibin,ii,jj))
                    end if
                 end do
              end do
              close(iunit)
           end do
        end do
     end do

     do jj=0,nterms_Ec
        write(indx_name,'(I2.2)'), jj
        do iplane=pnmin_,pnmax_
           call get_name(iplane,of_name)
           open(unit=iunit, file=trim(output_name)//'_budgets_TKE_cond_'//trim(indx_name)//'_'//trim(of_name), action='write')
           i = pnindx(iplane)
           do ibin=1,nbins_cond
              write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
              do j=jmid,jmax_
                 if (j.eq.jmax_) then
                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*E_budget_cond(i-1:i,j,ibin,jj))
                 else
                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*E_budget_cond(i-1:i,j,ibin,jj))
                 end if
              end do
           end do
           close(iunit)
        end do
     end do

  end if

  print *, irank, 'budget output done'

  return
end subroutine volumeStats_budgets_compute
