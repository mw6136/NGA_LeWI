module monitor_hit
  use precision
  implicit none

  ! Homogeneous isotropic turbulence
  real(WP) :: TKE,eps,n_kolmogorov,keta,tau_eddy,L_int,Re_L,eps_favre
  real(WP), dimension(:), pointer :: variance,zmean,covariance,skewness
  real(WP), dimension(:), pointer :: scalar_stat
  real(WP) :: urms2,nu_favre
  real(WP) :: lambda,Re_lambda
  logical :: is_joint_def, is_skw_def
  integer :: njoint, nstat
  integer, dimension(:), pointer :: isc_joint1,isc_joint2

end module monitor_hit


! ================================= !
! Initialize the monitor/hit module !
! ================================= !
subroutine monitor_hit_init
  use monitor_hit
  use data
  use parser
  implicit none
  integer :: isc, ii
  character(len=str_medium), dimension(:), pointer :: vars_in_joint
  logical, dimension(:), pointer :: found1,found2

  ! Create a file to monitor at each timestep
  call monitor_create_file_step('hit',8)
  call monitor_set_header(1,'TKE','r')
  call monitor_set_header(2,'urms','r')
  call monitor_set_header(3,'epsilon','r')
  call monitor_set_header(4,'tau_eddy','r')
  call monitor_set_header(5,'eta','r')
  call monitor_set_header(6,'keta','r')
  call monitor_set_header(7,'Re_turb','r')
  call monitor_set_header(8,'Re_lambda','r')
  
  if (nscalar.ge.1) then
     ! Allocate
     allocate(zmean(nscalar))
     allocate(variance(nscalar))
     allocate(skewness(nscalar))
     call parser_is_defined('Joint scalar PDF variables',is_joint_def) ! Account for covariances
     if (is_joint_def) then     
        call parser_getsize('Joint scalar PDF variables',njoint)  
        if (mod(njoint,2) .ne. 0) call die('monitor_hit_init: must provide scalar pairs for joint pdfs')
        allocate(vars_in_joint(njoint))
        njoint = njoint/2
        allocate(covariance(njoint))
        allocate(isc_joint1(njoint))
        allocate(isc_joint2(njoint))
        allocate(found1(njoint))
        allocate(found2(njoint))
     else
        njoint=0   
     end if
     call parser_read('Save skewness data',is_skw_def,.false.)
     nstat = 2                        ! means and variances only (2 stats per var)
     if (is_skw_def) nstat = nstat+1  ! also save skewness (1 more stat per var)
     allocate(scalar_stat(nstat*nscalar+njoint))

     ! Create header
     call monitor_create_file_step('hit_scalar',nstat*nscalar+njoint)
     do isc=1,nscalar
        call monitor_set_header(nstat*(isc-1)+1,'avg_'//SC_name(isc),'r')
        call monitor_set_header(nstat*(isc-1)+2,'var_'//SC_name(isc),'r')
        if (is_skw_def) call monitor_set_header(nstat*(isc-1)+3,'skw_'//SC_name(isc),'r')
     end do

     ! Headers for covariances
     if (is_joint_def) then     
        found1 = .false.
        found2 = .false.
        call parser_read('Joint scalar PDF variables',vars_in_joint)
        do ii = 1,njoint ! find scalar indices
           do isc=1,nscalar
              if (trim(SC_name(isc)).eq.trim(vars_in_joint(ii*2-1))) then
                 isc_joint1(ii) = isc
                 found1(ii) = .true.
              end if
              if (trim(SC_name(isc)).eq.trim(vars_in_joint(ii*2))  ) then
                 isc_joint2(ii) = isc
                 found2(ii) = .true.
              end if
           end do
        end do
        do ii=1,njoint ! check that all were found
           if ((found1(ii) .eq. .false.) .or. (found2(ii) .eq. .false.)) then
              call die('monitor_hit_init: invalid scalar requested')
           end if
        end do
        deallocate(found1)
        deallocate(found2)
        do isc=1,njoint ! set header
           call monitor_set_header(nstat*nscalar+isc,'cov_'//trim(SC_name(isc_joint1(isc))(5:7)) //'_'//trim(SC_name(isc_joint2(isc))(5:7)),'r')
        end do

     end if
  end if
  
  return
end subroutine monitor_hit_init


! ========================================== !
! Compute the quantities relevant to monitor !
! ========================================== !
subroutine monitor_hit_compute
  use monitor_hit
  use data
  use metric_generic
  use memory
  use math
  implicit none
  
  real(WP) :: buf1,buf2
  real(WP) :: RHOmean,div1,div2,div3
  integer  :: isc,i,j,k
  
  ! Favre(urms^2) = <rhoU.U/3> / <RHO>
  buf1 = 0.0_WP
  buf2 = 0.0_WP
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           buf1 = buf1 + vol(i,j) * ( &
                sum(interp_u_xm(i,j,:)*U(i-st1:i+st2,j,k)*rhoU(i-st1:i+st2,j,k)) + &
                sum(interp_v_ym(i,j,:)*V(i,j-st1:j+st2,k)*rhoV(i,j-st1:j+st2,k)) + &
                sum(interp_w_zm(i,j,:)*W(i,j,k-st1:k+st2)*rhoW(i,j,k-st1:k+st2)) )
           buf2 = buf2 + RHO(i,j,k) * vol(i,j)
        end do
     end do
  end do
  call parallel_sum(buf1,urms2)
  urms2 = urms2/vol_total
  call parallel_sum(buf2,RHOmean)
  RHOmean = RHOmean/vol_total
  urms2 = urms2/RHOmean
  urms2 = urms2/3.0_WP
  
  ! TKE = 3/2 * urms^2
  TKE = 1.5_WP*urms2
  
  ! Calculate Mean. variance, and skewness of the scalars
  do isc=1,nscalar
     buf1 = 0.0_WP ! mean
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              buf1 = buf1 + RHO(i,j,k)*SC(i,j,k,isc) * vol(i,j)
           end do
        end do
     end do
     call parallel_sum(buf1,zmean(isc))
     zmean(isc) = zmean(isc)/vol_total
     zmean(isc) = zmean(isc)/RHOmean
     buf1 = 0.0_WP ! variance
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              buf1 = buf1 + RHO(i,j,k)*(SC(i,j,k,isc)-zmean(isc))**2 * vol(i,j)
           end do
        end do
     end do
     call parallel_sum(buf1,variance(isc))
     variance(isc) = variance(isc)/vol_total
     variance(isc) = variance(isc)/RHOmean
     scalar_stat(nstat*(isc-1)+1) = zmean(isc)
     scalar_stat(nstat*(isc-1)+2) = variance(isc)

     if (is_skw_def) then
        buf1 = 0.0_WP ! skewness
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 buf2 = buf2 + RHO(i,j,k)*(SC(i,j,k,isc)-zmean(isc))**3 * vol(i,j)
              end do
           end do
        end do
        call parallel_sum(buf1,skewness(isc))
        skewness(isc) = skewness(isc)/vol_total
        skewness(isc) = skewness(isc)/RHOmean
        skewness(isc) = skewness(isc)/variance(isc)**1.5
        scalar_stat(nstat*(isc-1)+3) = skewness(isc)
     end if
  end do

  ! Covariances
  do isc=1,njoint
     buf1 = 0.0_WP
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              buf1 = buf1 + RHO(i,j,k)*SC(i,j,k,isc_joint1(isc))*SC(i,j,k,isc_joint2(isc)) * vol(i,j)
           end do
        end do
     end do
     call parallel_sum(buf1,covariance(isc))
     covariance(isc) = covariance(isc)/vol_total
     covariance(isc) = covariance(isc)/RHOmean
     covariance(isc) = covariance(isc) - zmean(isc_joint1(isc))*zmean(isc_joint2(isc))
     scalar_stat(nstat*nscalar+isc) = covariance(isc)
  end do

  ! Favre average of nu = <mu>/<RHO>
  buf1 = 0.0_WP
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           buf1 = buf1 + VISC(i,j,k) * vol(i,j)
        end do
     end do
  end do
  call parallel_sum(buf1,nu_favre)
  nu_favre = nu_favre/vol_total
  nu_favre = nu_favre/RHOmean
  
  ! Interpolate velocities
  call interpolate_velocities
  
  ! Compute VISC.S
  !     ( 1 4 6 )
  ! S = ( 4 2 5 )
  !     ( 6 5 3 )
  call strainrate_compute(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
  tmp1 = 2.0_WP*VISC*tmp1
  tmp2 = 2.0_WP*VISC*tmp2
  tmp3 = 2.0_WP*VISC*tmp3
  tmp4 = 2.0_WP*VISC*tmp4
  tmp5 = 2.0_WP*VISC*tmp5
  tmp6 = 2.0_WP*VISC*tmp6
  
  ! Compute -U.(div(VISC.S))
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           ! div(VISC.S)
           div1=div_u(i,j,0)*(interp_sc_x(i,j,-1)*tmp1(i-1,j,k)+interp_sc_x(i,j,0)*tmp1(i,j,k)) + div_u(i,j,1)*(interp_sc_x(i+1,j,-1)*tmp1(i,j,k)+interp_sc_x(i+1,j,0)*tmp1(i+1,j,k)) + &
                div_v(i,j,0)*(interp_sc_y(i,j,-1)*tmp4(i,j-1,k)+interp_sc_y(i,j,0)*tmp4(i,j,k)) + div_v(i,j,1)*(interp_sc_y(i,j+1,-1)*tmp4(i,j,k)+interp_sc_y(i,j+1,0)*tmp4(i,j+1,k)) + &
                div_w(i,j,0)*(interp_sc_z(i,j,-1)*tmp6(i,j,k-1)+interp_sc_z(i,j,0)*tmp6(i,j,k)) + div_w(i,j,1)*(interp_sc_z(i,j,-1)  *tmp6(i,j,k)+interp_sc_z(i,j,0)  *tmp6(i,j,k+1))
           div2=div_u(i,j,0)*(interp_sc_x(i,j,-1)*tmp4(i-1,j,k)+interp_sc_x(i,j,0)*tmp4(i,j,k)) + div_u(i,j,1)*(interp_sc_x(i+1,j,-1)*tmp4(i,j,k)+interp_sc_x(i+1,j,0)*tmp4(i+1,j,k)) + &
                div_v(i,j,0)*(interp_sc_y(i,j,-1)*tmp2(i,j-1,k)+interp_sc_y(i,j,0)*tmp2(i,j,k)) + div_v(i,j,1)*(interp_sc_y(i,j+1,-1)*tmp2(i,j,k)+interp_sc_y(i,j+1,0)*tmp2(i,j+1,k)) + &
                div_w(i,j,0)*(interp_sc_z(i,j,-1)*tmp5(i,j,k-1)+interp_sc_z(i,j,0)*tmp5(i,j,k)) + div_w(i,j,1)*(interp_sc_z(i,j,-1)  *tmp5(i,j,k)+interp_sc_z(i,j,0)  *tmp5(i,j,k+1))
           div3=div_u(i,j,0)*(interp_sc_x(i,j,-1)*tmp6(i-1,j,k)+interp_sc_x(i,j,0)*tmp6(i,j,k)) + div_u(i,j,1)*(interp_sc_x(i+1,j,-1)*tmp6(i,j,k)+interp_sc_x(i+1,j,0)*tmp6(i+1,j,k)) + &
                div_v(i,j,0)*(interp_sc_y(i,j,-1)*tmp5(i,j-1,k)+interp_sc_y(i,j,0)*tmp5(i,j,k)) + div_v(i,j,1)*(interp_sc_y(i,j+1,-1)*tmp5(i,j,k)+interp_sc_y(i,j+1,0)*tmp5(i,j+1,k)) + &
                div_w(i,j,0)*(interp_sc_z(i,j,-1)*tmp3(i,j,k-1)+interp_sc_z(i,j,0)*tmp3(i,j,k)) + div_w(i,j,1)*(interp_sc_z(i,j,-1)  *tmp3(i,j,k)+interp_sc_z(i,j,0)  *tmp3(i,j,k+1))
           ! div(VISC.S)
           div1=+sum(div_u(i,j,:)*tmp1(i-st1:i+st2,j,k)) &
                +sum(div_v(i,j,:)*tmp4(i,j-st1:j+st2,k)) &
                +sum(div_w(i,j,:)*tmp6(i,j,k-st1:k+st2))
           div2=+sum(div_u(i,j,:)*tmp4(i-st1:i+st2,j,k)) &
                +sum(div_v(i,j,:)*tmp2(i,j-st1:j+st2,k)) &
                +sum(div_w(i,j,:)*tmp5(i,j,k-st1:k+st2))
           div3=+sum(div_u(i,j,:)*tmp6(i-st1:i+st2,j,k)) &
                +sum(div_v(i,j,:)*tmp5(i,j-st1:j+st2,k)) &
                +sum(div_w(i,j,:)*tmp3(i,j,k-st1:k+st2))
           ! -U.(div(VISC.S))
           tmp7(i,j,k) = &
                - sum(interp_u_xm(i,j,:)*U(i-st1:i+st2,j,k)) * div1 &
                - sum(interp_v_ym(i,j,:)*V(i,j-st1:j+st2,k)) * div2 &
                - sum(interp_w_zm(i,j,:)*W(i,j,k-st1:k+st2)) * div3
           
        end do
     end do
  end do
  
  ! Epsilon = -<U.(div(VISC.S))> 
  buf1 = 0.0_WP
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           buf1 = buf1 + tmp7(i,j,k) * vol(i,j)
        end do
     end do
  end do
  call parallel_sum(buf1,eps)
  eps = eps/vol_total
  eps_favre = eps/RHOmean
  
  if (eps_favre.eq.0.0_WP .or. nu_favre.eq.0.0_WP) then
     eps_favre = 1.0_WP
     nu_favre  = 1.0_WP
  end if
  
  ! n_kolmogorov
  n_kolmogorov = (nu_favre**3/eps_favre)**(0.25_WP)
  keta = Pi*n_kolmogorov*dzi
  
  ! L_int
  L_int = TKE**(1.5_WP)/eps_favre
  
  ! Re_L
  Re_L = TKE**2/(eps_favre*nu_favre)
  
  ! Lambda
  if (Re_L.le.0.0_WP) then
     lambda = 0.0_WP
  else
     lambda = L_int*sqrt(10.0_WP)*Re_L**(-0.5_WP)
  end if
  
  ! Re_lambda
  Re_lambda = urms2**(0.5_WP)*lambda/nu_favre
  
  ! tau_eddy
  tau_eddy = TKE/eps_favre
  
  ! Transfer values to monitor
  call monitor_select_file('hit')
  call monitor_set_single_value(1,TKE)
  call monitor_set_single_value(2,sqrt(urms2))
  call monitor_set_single_value(3,eps_favre)
  call monitor_set_single_value(4,tau_eddy)
  call monitor_set_single_value(5,n_kolmogorov)
  call monitor_set_single_value(6,keta)
  call monitor_set_single_value(7,Re_L)
  call monitor_set_single_value(8,Re_lambda)
  if (nscalar.ge.1) then
     call monitor_select_file('hit_scalar')
     call monitor_set_array_values(scalar_stat)
  end if
  
  return
end subroutine monitor_hit_compute

