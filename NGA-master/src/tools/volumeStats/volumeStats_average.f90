module volumeStats_average
  use volumeStats

  ! Mean conditional variable in physical space
  real(WP), dimension(:,:),     pointer :: y_condm

  ! Progress variable pdf
  real(WP), dimension(:,:),     pointer :: pdf_C
  real(WP), dimension(:,:,:),   pointer :: pdf_Cz

  ! Average over z-direction: *m
  ! Conditional average over y- and z-directions: *m_c
  ! Conditional average over z-direction: *m_cz
  real(WP), dimension(:,:),       pointer :: RHOm, RHOm_x, RHOm_y, RHOm_z
  real(WP), dimension(:,:),       pointer :: RHOm_c
  real(WP), dimension(:,:,:),     pointer :: RHOm_cz, RHOm_cz_old
  real(WP), dimension(:,:),       pointer :: NU
  real(WP), dimension(:,:),       pointer :: NU_c
  real(WP), dimension(:,:,:),     pointer :: NU_cz
  real(WP), dimension(:,:,:),     pointer :: Um, Umi
  real(WP), dimension(:,:,:),     pointer :: Um_c
  real(WP), dimension(:,:,:,:),   pointer :: Um_cz, Um_cz_old
  real(WP), dimension(:,:,:),     pointer :: SCm
  real(WP), dimension(:,:,:),     pointer :: SCm_c
  real(WP), dimension(:,:,:,:),   pointer :: SCm_cz
  real(WP), dimension(:,:,:),     pointer :: DIFFm
  real(WP), dimension(:,:,:),     pointer :: DIFFm_c
  real(WP), dimension(:,:,:,:),   pointer :: DIFFm_cz
  real(WP), dimension(:,:,:),     pointer :: src_SCm
  real(WP), dimension(:,:,:),     pointer :: src_SCm_c
  real(WP), dimension(:,:,:,:),   pointer :: src_SCm_cz
  real(WP), dimension(:,:,:,:),   pointer :: TAUm
  real(WP), dimension(:,:,:,:),   pointer :: TAUm_c
  real(WP), dimension(:,:,:,:,:), pointer :: TAUm_cz
  real(WP), dimension(:,:,:),     pointer :: w_dot_cz
  real(WP), dimension(:,:,:),     pointer :: chi_cz

  ! Mean gradients
  real(WP), dimension(:,:,:),     pointer :: dRHOdt_cz
  real(WP), dimension(:,:,:,:),   pointer :: dUdt_cz

  real(WP), dimension(:,:,:),     pointer :: dRHOdxm
  real(WP), dimension(:,:,:),     pointer :: dRHOdxm_c
  real(WP), dimension(:,:,:),     pointer :: dRHOdxm_cc
  real(WP), dimension(:,:,:),     pointer :: dPdxm
  real(WP), dimension(:,:,:),     pointer :: dPdxm_c
  real(WP), dimension(:,:,:,:),   pointer :: dPdxm_cz
  real(WP), dimension(:,:,:,:),   pointer :: dUdx_tmpv
  real(WP), dimension(:,:,:,:),   pointer :: dUdxm
  real(WP), dimension(:,:,:,:),   pointer :: dRHOUiUjdxm
  real(WP), dimension(:,:,:),     pointer :: dRHOUiUidxm
  real(WP), dimension(:,:,:),     pointer :: dRHOUiUiUjdxm
  real(WP), dimension(:,:,:,:),   pointer :: dUdxm_c
  real(WP), dimension(:,:,:,:,:), pointer :: dUdxm_cz
  real(WP), dimension(:,:,:,:,:), pointer :: dUdxm_cz_reg
  real(WP), dimension(:,:,:,:),   pointer :: dSCdxm
  real(WP), dimension(:,:,:,:),   pointer :: dSCdxm_c
  real(WP), dimension(:,:,:,:),   pointer :: dTAUdxm
  real(WP), dimension(:,:,:,:),   pointer :: dTAUdxm_2cd
  real(WP), dimension(:,:,:,:),   pointer :: dTAUdxm_c
  real(WP), dimension(:,:,:,:,:), pointer :: dTAUdxm_cz

  ! Double and triple correlations
  real(WP), dimension(:,:,:,:),   pointer :: UiUjm
  real(WP), dimension(:,:,:,:),   pointer :: UiUjm_c
  real(WP), dimension(:,:,:,:,:), pointer :: UiUjm_cz
  real(WP), dimension(:,:,:,:),   pointer :: UiSCm
  real(WP), dimension(:,:,:,:),   pointer :: UiSCm_c
  real(WP), dimension(:,:,:,:,:), pointer :: UiSCm_cz
  real(WP), dimension(:,:,:,:),   pointer :: UiUiUjm
  real(WP), dimension(:,:,:,:),   pointer :: UiUiUjm_c
  real(WP), dimension(:,:,:,:,:), pointer :: UiUiUjm_cz

  ! Velocity and scalar moments
  real(WP), dimension(:,:,:,:),   pointer :: R_stress
  real(WP), dimension(:,:,:,:),   pointer :: R_stress_c
  real(WP), dimension(:,:,:,:),   pointer :: SC_flux
  real(WP), dimension(:,:,:,:),   pointer :: SC_flux_c
  real(WP), dimension(:,:,:,:),   pointer :: vel_triple_corr
  real(WP), dimension(:,:,:,:),   pointer :: vel_triple_corr_c

  ! For budgets
  real(WP), dimension(:,:,:,:),   pointer :: U_save_temp
  real(WP), dimension(:,:,:,:),   pointer :: Uc_save_temp
  real(WP), dimension(:,:,:),     pointer :: E_save_temp
  real(WP), dimension(:,:,:),     pointer :: Ec_save_temp
  real(WP), dimension(:,:,:),     pointer :: RHO_save_temp
  real(WP), dimension(:,:,:,:),   pointer :: U_cond_temp
  real(WP), dimension(:,:,:),     pointer :: E_cond_temp
  real(WP), dimension(:,:,:),     pointer :: RHO_cond_temp
  real(WP), dimension(:,:,:),     pointer :: U_temp_m
  real(WP), dimension(:,:),       pointer :: E_temp_m
  real(WP), dimension(:,:,:),     pointer :: tmp_E_phys, tmp_E_cond
  real(WP), dimension(:,:,:,:),   pointer :: tmp_SC_cz_1
  real(WP), dimension(:,:,:,:),   pointer :: tmp_SC_cz_2
  real(WP), dimension(:,:,:,:,:), pointer :: tmp_U_cz_1
  real(WP), dimension(:,:,:,:),   pointer :: tmp_U_cz_2
  real(WP), dimension(:,:,:,:),   pointer :: tmp_U_cz_3
  real(WP), dimension(:,:,:,:),   pointer :: tmp_U_cz_4
  real(WP), dimension(:,:,:,:),   pointer :: tmp_U_cz_5
  real(WP), dimension(:,:,:,:),   pointer :: tmp_E_cz_1
  real(WP), dimension(:,:,:),     pointer :: tmp_E_cz_2
  real(WP), dimension(:,:,:),     pointer :: tmp_E_cz_3
  real(WP), dimension(:,:,:),     pointer :: tmp_E_cz_4

  ! Momentum budgets - physical space - at cell faces
  integer, parameter :: nterms_Up=4
  real(WP), dimension(:,:,:,:), pointer :: U_budget_phys_cf

contains
  
  ! ========================================================== !
  ! Compute the average in the z-direction                     !
  ! ========================================================== !
  subroutine phys_avg_z(sol,mean)
    implicit none
    real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_), intent(in) :: sol
    real(WP), dimension(imin_:imax_,jmin_:jmax_), intent(inout) :: mean
    real(WP) :: tmp_nzi
    integer :: i,j,k

    tmp_nzi = 1.0_WP/real(kmax_-stz-(kmin_+stz)+1,WP)

    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do j=jmin_,jmax_
       do i=imin_,imax_
          mean(i,j) = mean(i,j)*real(ntime_curr-1,WP)
       end do
    end do
    !$OMP END DO
    !$OMP DO PRIVATE(j,i) REDUCTION(+:mean)
    do k=kmin_+stz,kmax_-stz
       do j=jmin_,jmax_
          do i=imin_,imax_
             mean(i,j) = mean(i,j) + sol(i,j,k)*tmp_nzi
          end do
       end do
    end do
    !$OMP END DO
    !$OMP DO PRIVATE(i)
    do j=jmin_,jmax_
       do i=imin_,imax_
          mean(i,j) = mean(i,j)/real(ntime_curr,WP)
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    return
  end subroutine phys_avg_z
  

  ! ===================================================================== !
  ! Compute the conditional average at each xy-location                   !
  ! ===================================================================== !
  subroutine cond_avg_z(sol,mean,nsz_opt_in)
    implicit none
    integer, optional :: nsz_opt_in
    integer :: nsz_opt
    real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_), intent(in) :: sol
    real(WP), dimension(imin_:imax_,jmid:jmax_,1:nbins_cond), intent(inout) :: mean
    integer(8),  dimension(imin_:imax_,jmid:jmax_,1:nbins_cond) :: nsz
    integer :: i,j,k,ii

    if (.not.combust) return
    if (present(nsz_opt_in)) then
       nsz_opt = nsz_opt_in
    else
       nsz_opt = 1
    end if

    ! Initialize count
    if (nsz_opt.eq.1) then
       do j=jmid,jmax_
          do i=imin_,imax_
             nsz(i,j,:) = nsz_cond_total(i,j,:)
          end do
       end do
    else
       do j=jmid,jmax_
          do i=imin_,imax_
             nsz(i,j,:) = nsz_cond_temp_total(i,j,:)
          end do
       end do
    end if
       
    !$OMP PARALLEL
    !$OMP DO PRIVATE(j,i)
    do ii=1,nbins_cond
       do j=jmid,jmax_
          do i=imin_,imax_
             mean(i,j,ii) = mean(i,j,ii)*real(nsz(i,j,ii),WP)
          end do
       end do
    end do
    !$OMP END DO
!!$    !$OMP DO PRIVATE(j,i,ii) REDUCTION(+:nsz,mean)
!!$    ! bin_index needs be first set in cond_map
!!$    do k=kmin_,kmax_
!!$       do j=jmin_,jmid-1
!!$          do i=imin_,imax_
!!$             ii = bin_index(i,j,k)
!!$             nsz(i,jmax_-j-1,ii) = nsz(i,jmax_-j-1,ii) + 1
!!$             mean(i,jmax_-j-1,ii) = mean(i,jmax_-j-1,ii) + sol(i,j,k)
!!$          end do
!!$       end do
!!$    end do
    !$OMP DO PRIVATE(j,i,ii) REDUCTION(+:nsz,mean)
    do k=kmin_+stz,kmax_-stz
       do j=jmid,jmax_
          do i=imin_,imax_
             ii = bin_index(i,j,k)
             nsz(i,j,ii) = nsz(i,j,ii) + 1
             mean(i,j,ii) = mean(i,j,ii) + sol(i,j,k)
          end do
       end do
    end do
    !$OMP END DO
    !$OMP DO PRIVATE(i,j)
    do ii=1,nbins_cond
       do j=jmid,jmax_
          do i=imin_,imax_
             if (nsz(i,j,ii).gt.0) then
                mean(i,j,ii) = mean(i,j,ii)/real(nsz(i,j,ii),WP)
             else
                mean(i,j,ii) = 0.0_WP
             end if
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! Save for cumulative count
    if (nsz_opt.eq.1) then
       do j=jmid,jmax_
          do i=imin_,imax_
             nsz_cond(i,j,:) = nsz(i,j,:) - nsz_cond_total(i,j,:)
          end do
       end do
    else
       do j=jmid,jmax_
          do i=imin_,imax_
             nsz_cond_temp(i,j,:) = nsz(i,j,:) - nsz_cond_temp_total(i,j,:)
          end do
       end do
    end if

    return
  end subroutine cond_avg_z
  

  ! ==================================================================== !
  ! Compute the conditional average at each x-location                   !
  ! ==================================================================== !
  subroutine cond_avg_yz(sol,mean)
    implicit none
    real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_), intent(in) :: sol
    real(WP), dimension(imin_:imax_,1:nbins_cond), intent(inout) :: mean
    integer(8), dimension(imin_:imax_,1:nbins_cond) :: nsyz
    integer :: i,j,k,ii

    if (.not.combust) return

    ! Initialize count
    do i=imin_,imax_
       nsyz(i,:) = nsyz_cond_total(i,:)
    end do

    !$OMP PARALLEL
    !$OMP DO PRIVATE(i)
    do ii=1,nbins_cond
       do i=imin_,imax_
          mean(i,ii) = mean(i,ii)*real(nsyz(i,ii),WP)
       end do
    end do
    !$OMP END DO
    !$OMP DO PRIVATE(j,i,ii) REDUCTION(+:nsyz)
    ! bin_index needs be first set in cond_map
    do k=kmin_+stz,kmax_-stz
       do j=jmid,jmax_
          do i=imin_,imax_
             !ii = min(1,max(nbins_cond,bin_index(i,j,k)))
             ii = bin_index(i,j,k)
             nsyz(i,ii) = nsyz(i,ii) + 1
          end do
       end do
    end do
    !$OMP END DO
    !$OMP DO PRIVATE(j,i,ii) REDUCTION(+:mean)
    do k=kmin_+stz,kmax_-stz
       do j=jmid,jmax_
          do i=imin_,imax_
             !ii = min(1,max(nbins_cond,bin_index(i,j,k)))
             ii = bin_index(i,j,k)
             mean(i,ii) = mean(i,ii) + sol(i,j,k)
          end do
       end do
    end do
    !$OMP END DO
    !$OMP DO PRIVATE(i)
    do ii=1,nbins_cond
       do i=imin_,imax_
          if (nsyz(i,ii).gt.0) then
             mean(i,ii) = mean(i,ii)/real(nsyz(i,ii),WP)
          else
             mean(i,ii) = 0.0_WP
          end if
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! Save for cumulative count
    do i=imin_,imax_
       nsyz_cond(i,:) = nsyz(i,:) - nsyz_cond_total(i,:)
    end do

    return
  end subroutine cond_avg_yz
  

  ! ===================================================================== !
  ! Compute regularized derivatives                                       !
  ! ===================================================================== !
  subroutine reg_deriv(mean,grad,dir)
    use volumeStats_metric
    implicit none
    integer, intent(in) :: dir
    real(WP), dimension(imin_:imax_,jmid:jmax_,1:nbins_cond), intent(in)  :: mean
    real(WP), dimension(imin_:imax_,jmid:jmax_,1:nbins_cond), intent(out) :: grad
    real(WP), dimension(:,:), allocatable :: D_mat, A_mat, E_mat, L_mat, H_mat
    real(WP), dimension(:),   allocatable :: g_vec, s_vec, u_vec, f_vec
    integer,  dimension(:),   allocatable :: pivot
    integer  :: cmin_loc, cmax_loc, imin_loc, imax_loc, jmin_loc, jmax_loc
    real(WP) :: alpha, eps
    integer :: i,j,k,ip,it
    integer :: ii,jj,kk
    integer :: err

    ! Initialize output
    grad = 0.0_WP

    if (.not.combust) return

    ! Fixed parameters
    eps   = 1.0e-9_WP

    select case(dir)
    case(0)
       ! Set alpha higher for x
       alpha = 5.0e-8_WP

       ! Allocate arrays for x-derivative
       allocate(D_mat(1:nbins_cond,1:nbins_cond+1))
       allocate(A_mat(1:nbins_cond,1:nbins_cond+1))

       do ip=1,nbins_cond
          ! Differentiation operator
          D_mat(ip,ip  ) = -dci(ip)
          D_mat(ip,ip+1) = +dci(ip)
          ! Integration operator
          A_mat(ip:nbins_cond,ip  ) = 1.0_WP/dci(ip)
          A_mat(ip:nbins_cond,ip+1) = 1.0_WP/dci(ip) !added this late
       end do

       !! !$OMP PARALLEL COPYIN(E_mat,L_mat,H_mat,u_vec,f_vec,g_vec,s_vec,pivot) SHARED(mean,grad,D_mat,A_mat) PRIVATE(i,cmin_loc,cmax_loc,ip,it,ii,jj,kk,err,E_mat,L_mat,H_mat,u_vec,f_vec,g_vec,s_vec,pivot)
       allocate(E_mat(1:nbins_cond  ,1:nbins_cond  ))
       allocate(L_mat(1:nbins_cond+1,1:nbins_cond+1))
       allocate(H_mat(1:nbins_cond+1,1:nbins_cond+1))
       allocate(u_vec(1:nbins_cond+1))
       allocate(f_vec(1:nbins_cond+1))
       allocate(g_vec(1:nbins_cond+1))
       allocate(s_vec(1:nbins_cond+1))
       allocate(pivot(1:nbins_cond+1))

       !! !$OMP DO
       do j=jmid,jmax_
          loop0: do i=imin_,imax_
             ! Find the x-min and x-max points with non-NaN values of mean
             ! Use the conditional mean U-velocity
             ip = 1
             do while (((Um_cz(i,j,ip,1).ne.Um_cz(i,j,ip,1)).or. (abs(Um_cz(i,j,ip,1))>1.0e50_WP)).and.(ip.lt.nbins_cond))
                ip = ip+1
             end do
             cmin_loc = ip
             do while (((Um_cz(i,j,ip,1).eq.Um_cz(i,j,ip,1)).and.(abs(Um_cz(i,j,ip,1))<1.0e50_WP)).and.(ip.lt.nbins_cond))
                ip = ip+1
             end do
             cmax_loc = ip-1
             if (cmin_loc.ge.cmax_loc) then
                cycle loop0
             end if

             ! Initialize solution (f) and derivative (u) arrays
             f_vec(cmin_loc:cmax_loc) = mean(i,j,cmin_loc:cmax_loc) - mean(i,j,cmin_loc)
             u_vec(cmin_loc) = 0.0_WP
             do ip=cmin_loc,cmax_loc-1
                u_vec(ip+1) = (mean(i,j,ip+1)-mean(i,j,ip))*dci(ip)
             end do
             u_vec(cmax_loc+1) = 0.0_WP

             ! Gradient descent iteration to converge the Euler-Lagrange equation
             do it=1,50
                ! Diagonal matrix
                do ip=cmin_loc,cmax_loc
                   E_mat(ip,ip) = 1.0_WP/sqrt( (u_vec(ip+1)-u_vec(ip))**2 + eps )/dci(ip)
                end do
                ! L = D* E D
                ! C = A* B : c(i,j) = c(i,j) + a(k,i)*b(k,j)
                H_mat = 0.0_WP
                do kk=cmin_loc,cmax_loc
                   do jj=cmin_loc,cmax_loc
                      do ii=cmin_loc,cmax_loc+1
                         H_mat(ii,jj) = H_mat(ii,jj) + D_mat(kk,ii)*E_mat(kk,jj)
                      end do
                   end do
                end do
                ! C = A B  : c(i,j) = c(i,j) + a(i,k)*b(k,j)
                L_mat = 0.0_WP
                do kk=cmin_loc,cmax_loc
                   do jj=cmin_loc,cmax_loc+1
                      do ii=cmin_loc,cmax_loc+1
                         L_mat(ii,jj) = L_mat(ii,jj) + H_mat(ii,kk)*D_mat(kk,jj)
                      end do
                   end do
                end do
                ! H = A* A + alpha L
                H_mat = 0.0_WP
                do kk=cmin_loc,cmax_loc
                   do jj=cmin_loc,cmax_loc+1
                      do ii=cmin_loc,cmax_loc+1
                         H_mat(ii,jj) = H_mat(ii,jj) + A_mat(kk,ii)*A_mat(kk,jj)
                      end do
                   end do
                end do
                do jj=cmin_loc,cmax_loc+1
                   do ii=cmin_loc,cmax_loc+1
                      H_mat(ii,jj) = H_mat(ii,jj) + alpha*L_mat(ii,jj)
                   end do
                end do
                ! RHS: g = A* (A u - f) + alpha L u
                g_vec = 0.0_WP
                s_vec = 0.0_WP ! use as a temporary variable
                do jj=cmin_loc,cmax_loc+1
                   do ii=cmin_loc,cmax_loc
                      s_vec(ii) = s_vec(ii) + A_mat(ii,jj)*u_vec(jj) ! A u
                   end do
                end do
                do ii=cmin_loc,cmax_loc
                   s_vec(ii) = s_vec(ii) - f_vec(ii) ! A u - f
                end do
                do jj=cmin_loc,cmax_loc
                   do ii=cmin_loc,cmax_loc+1
                      g_vec(ii) = g_vec(ii) + A_mat(jj,ii)*s_vec(jj) ! A* (A u -f)
                   end do
                end do
                do jj=cmin_loc,cmax_loc+1
                   do ii=cmin_loc,cmax_loc+1
                      g_vec(ii) = g_vec(ii) + alpha*L_mat(ii,jj)*u_vec(jj) ! A* (A u -f) + alpha L u
                   end do
                end do

                ! Solve H s = -g using DGESV
                g_vec  = -g_vec
                ii = cmax_loc - cmin_loc + 2
                call DGESV(ii, 1, H_mat(cmin_loc:cmax_loc+1,cmin_loc:cmax_loc+1), &
                     ii, pivot(cmin_loc:cmax_loc+1), g_vec(cmin_loc:cmax_loc+1), ii, err)

                ! u_n+1 = s + u_n
                do ii=cmin_loc,cmax_loc+1
                   u_vec(ii) = g_vec(ii) + u_vec(ii)
                end do
             end do

             ! Save the derivative
             grad(i,j,cmin_loc:cmax_loc+1) = u_vec(cmin_loc:cmax_loc+1)
          end do loop0
       end do
       !! !$OMP END DO

       ! Clean up
       deallocate(E_mat, L_mat, H_mat)
       deallocate(u_vec, f_vec, g_vec, s_vec)
       deallocate(pivot)

       !! !$OMP END PARALLEL
       deallocate(D_mat, A_mat)


    case(1)
       ! Set alpha higher for x
       alpha = 5.0e-9_WP

       ! Allocate arrays for x-derivative
       allocate(D_mat(imin_:imax_  ,imin_:imax_+1))
       allocate(A_mat(imin_:imax_  ,imin_:imax_+1))

       do i=imin_,imax_
          ! Differentiation operator
          D_mat(i,i  ) = -dxmi(i)
          D_mat(i,i+1) = +dxmi(i)
          ! Integration operator
          A_mat(i:imax_,i  ) = 1.0_WP/dxmi(i)
          A_mat(i:imax_,i+1) = 1.0_WP/dxmi(i) !added this late
       end do

       !! !$OMP PARALLEL COPYIN(E_mat,L_mat,H_mat,u_vec,f_vec,g_vec,s_vec,pivot) SHARED(mean,grad,D_mat,A_mat) PRIVATE(i,imin_loc,imax_loc,ip,it,ii,jj,kk,err,E_mat,L_mat,H_mat,u_vec,f_vec,g_vec,s_vec,pivot)
       allocate(E_mat(imin_:imax_  ,imin_:imax_  ))
       allocate(L_mat(imin_:imax_+1,imin_:imax_+1))
       allocate(H_mat(imin_:imax_+1,imin_:imax_+1))
       allocate(u_vec(imin_:imax_+1))
       allocate(f_vec(imin_:imax_+1))
       allocate(g_vec(imin_:imax_+1))
       allocate(s_vec(imin_:imax_+1))
       allocate(pivot(imin_:imax_+1))

       !! !$OMP DO
       do j=jmid,jmax_
          loop1: do ip=1,nbins_cond
             ! Find the x-min and x-max points with non-NaN values of mean
             ! Use the conditional mean U-velocity
             i = imin_+4
!!$             loop11: do while (i.lt.imax_)
!!$                do while (Um_cz(i,j,ip,1).ne.Um_cz(i,j,ip,1))
!!$                   i = i+1
!!$                end do
!!$                imin_loc = i
!!$                do while (Um_cz(i,j,ip,1).eq.Um_cz(i,j,ip,1))
!!$                   i = i+1
!!$                end do
!!$                imax_loc = i-1
!!$                if ((imin_loc.ge.imax_loc).or.((imax_loc-imin_loc).lt.3)) then
!!$                   i = i + 1
!!$                   cycle loop11
!!$                else
!!$                   exit loop11
!!$                end if
!!$             end do loop11
!!$             if (i.eq.imax_.or.i.eq.(imax_-1)) then
!!$                print *, '(x) cycling', j, ip, i, imax_
!!$                cycle loop1
!!$             end if
             do while (((Um_cz(i,j,ip,1).ne.Um_cz(i,j,ip,1)).or. (abs(Um_cz(i,j,ip,1))>1.0e3_WP)).and.(i.lt.imax_))
                i = i+1
             end do
             imin_loc = i
             do while (((Um_cz(i,j,ip,1).eq.Um_cz(i,j,ip,1)).and.(abs(Um_cz(i,j,ip,1))<1.0e3_WP)).and.(i.lt.imax_))
                i = i+1
             end do
             imax_loc = i-1
             if (imin_loc.ge.imax_loc) then
                !print *, '(x) cycling', j, ip, imin_loc, imax_loc
                cycle loop1
             end if

             ! Initialize solution (f) and derivative (u) arrays
             f_vec(imin_loc:imax_loc) = mean(imin_loc:imax_loc,j,ip) - mean(imin_loc,j,ip)
             u_vec(imin_loc) = 0.0_WP
             do i=imin_loc,imax_loc-1
                u_vec(i+1) = (mean(i+1,j,ip)-mean(i,j,ip))/(xm(i+1)-xm(i))
             end do
             u_vec(imax_loc+1) = 0.0_WP

             ! Gradient descent iteration to converge the Euler-Lagrange equation
             do it=1,50
                ! Diagonal matrix
                do i=imin_loc,imax_loc
                   E_mat(i,i) = (xm(i+1)-xm(i))/sqrt( (u_vec(i+1)-u_vec(i))**2 + eps )
                end do
                ! L = D* E D
                ! C = A* B : c(i,j) = c(i,j) + a(k,i)*b(k,j)
                H_mat = 0.0_WP
                do kk=imin_loc,imax_loc
                   do jj=imin_loc,imax_loc
                      do ii=imin_loc,imax_loc+1
                         H_mat(ii,jj) = H_mat(ii,jj) + D_mat(kk,ii)*E_mat(kk,jj)
                      end do
                   end do
                end do
                ! C = A B  : c(i,j) = c(i,j) + a(i,k)*b(k,j)
                L_mat = 0.0_WP
                do kk=imin_loc,imax_loc
                   do jj=imin_loc,imax_loc+1
                      do ii=imin_loc,imax_loc+1
                         L_mat(ii,jj) = L_mat(ii,jj) + H_mat(ii,kk)*D_mat(kk,jj)
                      end do
                   end do
                end do
                ! H = A* A + alpha L
                H_mat = 0.0_WP
                do kk=imin_loc,imax_loc
                   do jj=imin_loc,imax_loc+1
                      do ii=imin_loc,imax_loc+1
                         H_mat(ii,jj) = H_mat(ii,jj) + A_mat(kk,ii)*A_mat(kk,jj)
                      end do
                   end do
                end do
                do jj=imin_loc,imax_loc+1
                   do ii=imin_loc,imax_loc+1
                      H_mat(ii,jj) = H_mat(ii,jj) + alpha*L_mat(ii,jj)
                   end do
                end do
                ! RHS: g = A* (A u - f) + alpha L u
                g_vec = 0.0_WP
                s_vec = 0.0_WP ! use as a temporary variable
                do jj=imin_loc,imax_loc+1
                   do ii=imin_loc,imax_loc
                      s_vec(ii) = s_vec(ii) + A_mat(ii,jj)*u_vec(jj) ! A u
                   end do
                end do
                do ii=imin_loc,imax_loc
                   s_vec(ii) = s_vec(ii) - f_vec(ii) ! A u - f
                end do
                do jj=imin_loc,imax_loc
                   do ii=imin_loc,imax_loc+1
                      g_vec(ii) = g_vec(ii) + A_mat(jj,ii)*s_vec(jj) ! A* (A u -f)
                   end do
                end do
                do jj=imin_loc,imax_loc+1
                   do ii=imin_loc,imax_loc+1
                      g_vec(ii) = g_vec(ii) + alpha*L_mat(ii,jj)*u_vec(jj) ! A* (A u -f) + alpha L u
                   end do
                end do

                ! Solve H s = -g using DGESV
                g_vec  = -g_vec
                ii = imax_loc - imin_loc + 2
                call DGESV(ii, 1, H_mat(imin_loc:imax_loc+1,imin_loc:imax_loc+1), &
                     ii, pivot(imin_loc:imax_loc+1), g_vec(imin_loc:imax_loc+1), ii, err)

                ! u_n+1 = s + u_n
                do ii=imin_loc,imax_loc+1
                   u_vec(ii) = g_vec(ii) + u_vec(ii)
                end do
             end do
!!$             if (abs(err).gt.0) then
!!$                print *, irank, ' Something went wrong with DGESV (x) at ', ip, j
!!$             end if

             ! Save the derivative
             grad(imin_loc:imax_loc+1,j,ip) = u_vec(imin_loc:imax_loc+1)
          end do loop1
       end do
       !! !$OMP END DO

       ! Clean up
       deallocate(E_mat, L_mat, H_mat)
       deallocate(u_vec, f_vec, g_vec, s_vec)
       deallocate(pivot)

       !! !$OMP END PARALLEL
       deallocate(D_mat, A_mat)

    case(2)
       ! Set alpha for y
       alpha = 1.0e-9_WP

       ! Allocate arrays for y-derivative
       allocate(D_mat(jmid:jmax_  ,jmid:jmax_+1))
       allocate(A_mat(jmid:jmax_  ,jmid:jmax_+1))

       do j=jmid,jmax_
          ! Differentiation operator
          D_mat(j,j  ) = -dymi(j)
          D_mat(j,j+1) = +dymi(j)
          ! Integration operator
          A_mat(j:jmax_,j) = 1.0_WP/dymi(j)
       end do

       !! !$OMP PARALLEL COPYIN(E_mat,L_mat,H_mat,u_vec,f_vec,g_vec,s_vec,pivot) SHARED(mean,grad,D_mat,A_mat) PRIVATE(j,jmin_loc,jmax_loc,ip,it,ii,jj,kk,err,E_mat,L_mat,H_mat,u_vec,f_vec,g_vec,s_vec,pivot)
       allocate(E_mat(jmid:jmax_  ,jmid:jmax_  ))
       allocate(L_mat(jmid:jmax_+1,jmid:jmax_+1))
       allocate(H_mat(jmid:jmax_+1,jmid:jmax_+1))
       allocate(u_vec(jmid:jmax_+1))
       allocate(f_vec(jmid:jmax_+1))
       allocate(g_vec(jmid:jmax_+1))
       allocate(s_vec(jmid:jmax_+1))
       allocate(pivot(jmid:jmax_+1))

       !! !$OMP DO
       do i=imin_,imax_
          loop2: do ip=1,nbins_cond
             ! Find the y-min and y-max points with non-NaN values of mean
             ! Use the conditional mean U-velocity
             j = jmid+1
             do while (((Um_cz(i,j,ip,1).ne.Um_cz(i,j,ip,1)).or. (abs(Um_cz(i,j,ip,1))>1.0e3_WP)).and.(j.lt.jmax_))
                j = j+1
             end do
             jmin_loc = j
             do while (((Um_cz(i,j,ip,1).eq.Um_cz(i,j,ip,1)).and.(abs(Um_cz(i,j,ip,1))<1.0e3_WP)).and.(j.lt.jmax_))
                j = j+1
             end do
             jmax_loc = j-1
             if (jmin_loc.ge.jmax_loc) then
                !print *, '(y) cycling', i, ip, jmin_loc, jmax_loc
                cycle loop2
             end if

             ! Initialize solution (f) and derivative (u) arrays
             f_vec(jmin_loc:jmax_loc) = mean(i,jmin_loc:jmax_loc,ip) - mean(i,jmin_loc,ip)
             u_vec(jmin_loc) = 0.0_WP
             do j=jmin_loc,jmax_loc-1
                u_vec(j+1) = (mean(i,j+1,ip)-mean(i,j,ip))*dymi(j) !/(ym(j+1)-ym(j))
             end do
             u_vec(jmax_loc+1) = 0.0_WP

             ! Gradient descent iteration to converge the Euler-Lagrange equation
             do it=1,50
                ! Diagonal matrix
                do j=jmin_loc,jmax_loc
                   !E_mat(j,j) = (ym(j+1)-ym(j))/sqrt( (u_vec(j+1)-u_vec(j))**2 + eps )
                   E_mat(j,j) = 1.0_WP/(dymi(j)*sqrt( (u_vec(j+1)-u_vec(j))**2 + eps ))
                end do
                ! L = D* E D
                ! C = A* B : c(i,j) = c(i,j) + a(k,i)*b(k,j)
                H_mat = 0.0_WP
                do kk=jmin_loc,jmax_loc
                   do jj=jmin_loc,jmax_loc
                      do ii=jmin_loc,jmax_loc+1
                         H_mat(ii,jj) = H_mat(ii,jj) + D_mat(kk,ii)*E_mat(kk,jj)
                      end do
                   end do
                end do
                ! C = A B  : c(i,j) = c(i,j) + a(i,k)*b(k,j)
                L_mat = 0.0_WP
                do kk=jmin_loc,jmax_loc
                   do jj=jmin_loc,jmax_loc+1
                      do ii=jmin_loc,jmax_loc+1
                         L_mat(ii,jj) = L_mat(ii,jj) + H_mat(ii,kk)*D_mat(kk,jj)
                      end do
                   end do
                end do
                ! H = A* A + alpha L
                H_mat = 0.0_WP
                do kk=jmin_loc,jmax_loc
                   do jj=jmin_loc,jmax_loc+1
                      do ii=jmin_loc,jmax_loc+1
                         H_mat(ii,jj) = H_mat(ii,jj) + A_mat(kk,ii)*A_mat(kk,jj)
                      end do
                   end do
                end do
                do jj=jmin_loc,jmax_loc+1
                   do ii=jmin_loc,jmax_loc+1
                      H_mat(ii,jj) = H_mat(ii,jj) + alpha*L_mat(ii,jj)
                   end do
                end do
                ! RHS: g = A* (A u - f) + alpha L u
                g_vec = 0.0_WP
                s_vec = 0.0_WP ! use as a temporary variable
                do jj=jmin_loc,jmax_loc+1
                   do ii=jmin_loc,jmax_loc
                      s_vec(ii) = s_vec(ii) + A_mat(ii,jj)*u_vec(jj) ! A u
                   end do
                end do
                do ii=jmin_loc,jmax_loc
                   s_vec(ii) = s_vec(ii) - f_vec(ii) ! A u - f
                end do
                do jj=jmin_loc,jmax_loc
                   do ii=jmin_loc,jmax_loc+1
                      g_vec(ii) = g_vec(ii) + A_mat(jj,ii)*s_vec(jj) ! A* (A u -f)
                   end do
                end do
                do jj=jmin_loc,jmax_loc+1
                   do ii=jmin_loc,jmax_loc+1
                      g_vec(ii) = g_vec(ii) + alpha*L_mat(ii,jj)*u_vec(jj) ! A* (A u -f) + alpha L u
                   end do
                end do

                ! Solve H s = -g using DGESV
                g_vec  = -g_vec
                ii = jmax_loc - jmin_loc + 2
                call DGESV(ii, 1, H_mat(jmin_loc:jmax_loc+1,jmin_loc:jmax_loc+1), &
                     ii, pivot(jmin_loc:jmax_loc+1), g_vec(jmin_loc:jmax_loc+1), ii, err)

                ! u_n+1 = s + u_n
                do ii=jmin_loc,jmax_loc+1
                   u_vec(ii) = g_vec(ii) + u_vec(ii)
                end do
             end do
!!$             if (abs(err).gt.0) then
!!$                print *, irank, ' Something went wrong with DGESV (y) at ', ip, i
!!$             end if

             ! Save the derivative
             grad(i,jmin_loc:jmax_loc+1,ip) = u_vec(jmin_loc:jmax_loc+1)
          end do loop2
       end do
       !! !$OMP END DO

       ! Clean up
       deallocate(E_mat, L_mat, H_mat)
       deallocate(u_vec, f_vec, g_vec, s_vec)
       deallocate(pivot)

       !! !$OMP END PARALLEL
       deallocate(D_mat, A_mat)

    end select

    return
  end subroutine reg_deriv

end module volumeStats_average

! ========================================================== !
! Initialize the module                                      !
! ========================================================== !
subroutine volumeStats_average_init
  use volumeStats_average
  implicit none

  ! For unsteady terms
  allocate(U_save_temp (imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3))
  allocate(Uc_save_temp(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:2))
  allocate(E_save_temp (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(Ec_save_temp(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(RHO_save_temp (imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(U_temp_m    (imin_:imax_,jmin_:jmax_,1:2))
  allocate(E_temp_m    (imin_:imax_,jmin_:jmax_))
  allocate(U_cond_temp (imin_:imax_,jmid:jmax_,1:nbins_cond,1:3))
  allocate(E_cond_temp (imin_:imax_,jmid:jmax_,1:nbins_cond))
  allocate(RHO_cond_temp (imin_:imax_,jmid:jmax_,1:nbins_cond))

  ! Momentum budgets - physical 
  allocate(U_budget_phys_cf(imin_:imax_,jmin_:jmax_,1:2,0:nterms_Up))
  U_budget_phys_cf = 0.0_WP

  ! Initialize velocity module
  call volumeStats_velocity_init

  ! Physical mean
  allocate(y_condm(imin_:imax_,jmin_:jmax_))
  allocate(RHOm   (imin_:imax_,jmin_:jmax_))
  allocate(RHOm_x (imin_:imax_,jmin_:jmax_))
  allocate(RHOm_y (imin_:imax_,jmin_:jmax_))
  allocate(RHOm_z (imin_:imax_,jmin_:jmax_))
  allocate(NU     (imin_:imax_,jmin_:jmax_))
  allocate(Um     (imin_:imax_,jmin_:jmax_,1:3)) ! At cell faces
  allocate(Umi    (imin_:imax_,jmin_:jmax_,1:3)) ! At cell centers
!!$  allocate(Um     (imino_:imaxo_,jmino_:jmaxo_,1:3)); Um = 0.0_WP
  allocate(SCm    (imin_:imax_,jmin_:jmax_,1:nscalar))
  allocate(DIFFm  (imin_:imax_,jmin_:jmax_,1:N_tot+1))
  allocate(src_SCm(imin_:imax_,jmin_:jmax_,1:N_tot+1))
  allocate(TAUm   (imin_:imax_,jmin_:jmax_,1:3,1:3)) ! i,j
  
  ! Physical mean gradients
  if (need_gradients) then
     allocate(dRHOdxm (imin_:imax_,jmin_:jmax_,1:3))
     allocate(dPdxm   (imin_:imax_,jmin_:jmax_,1:3))
     allocate(dUdx_tmpv(imin_:imax_,jmin_:jmax_,1:3,1:3)) ! i,j
     allocate(dUdxm   (imin_:imax_,jmin_:jmax_,1:3,1:3)) ! i,j
     allocate(dRHOUiUjdxm(imin_:imax_,jmin_:jmax_,1:3,1:3)) ! i,j
     allocate(dRHOUiUidxm(imin_:imax_,jmin_:jmax_,1:3)) ! i,j
     allocate(dRHOUiUiUjdxm(imin_:imax_,jmin_:jmax_,1:3)) ! i,j
     allocate(dTAUdxm (imin_:imax_,jmin_:jmax_,1:3,1:3)) ! i,j
     allocate(dTAUdxm_2cd(imin_:imax_,jmin_:jmax_,1:3,1:3)) ! i,j
     allocate(dSCdxm (imin_:imax_,jmin_:jmax_,1:3,1:nscalar)) ! j,k
  end if

  ! Physical space correlations
  allocate(UiUjm  (imin_:imax_,jmin_:jmax_,1:3,1:3))
  allocate(UiSCm  (imin_:imax_,jmin_:jmax_,1:3,1:nscalar))
  allocate(UiUiUjm(imin_:imax_,jmin_:jmax_,1:3,1:3))

  ! Physical space velocity and scalar moments
  allocate(R_stress(imin_:imax_,jmin_:jmax_,1:3,1:3))
  allocate(SC_flux (imin_:imax_,jmin_:jmax_,1:3,1:nscalar))
  allocate(vel_triple_corr(imin_:imax_,jmin_:jmax_,1:3,1:3))

  ! For physical space budgets
  allocate(tmp_E_phys(imin_:imax_,jmin_:jmax_,1:3))

  ! Conditional mean -- yz average
  allocate(RHOm_c   (imin_:imax_,1:nbins_cond))
  allocate(NU_c     (imin_:imax_,1:nbins_cond))
  allocate(Um_c     (imin_:imax_,1:nbins_cond,1:3))
  allocate(SCm_c    (imin_:imax_,1:nbins_cond,1:nscalar))
  allocate(DIFFm_c  (imin_:imax_,1:nbins_cond,1:N_tot+1))
  allocate(src_SCm_c(imin_:imax_,1:nbins_cond,1:N_tot+1))
  allocate(TAUm_c   (imin_:imax_,1:nbins_cond,1:3,1:3)) ! i,j
  
  ! Conditional mean gradients -- yz average
  if (need_gradients) then
     allocate(dRHOdxm_c (imin_:imax_,1:nbins_cond,1:3))
     allocate(dRHOdxm_cc(imin_:imax_,1:nbins_cond,1:3))
     allocate(dPdxm_c  (imin_:imax_,1:nbins_cond,1:3))
     allocate(dUdxm_c  (imin_:imax_,1:nbins_cond,1:3,1:3)) ! i,j
     allocate(dTAUdxm_c(imin_:imax_,1:nbins_cond,1:3,1:3)) ! i,j
     allocate(dSCdxm_c (imin_:imax_,1:nbins_cond,1:3,1:nscalar)) ! j,k
  end if

  ! Conditional correlations -- yz average
  allocate(UiUjm_c  (imin_:imax_,1:nbins_cond,1:3,1:3))
  allocate(UiSCm_c  (imin_:imax_,1:nbins_cond,1:3,1:nscalar))
  allocate(UiUiUjm_c(imin_:imax_,1:nbins_cond,1:3,1:3))

  ! Conditional velocity and scalar moments -- yz average
  allocate(R_stress_c(imin_:imax_,1:nbins_cond,1:3,1:3))
  allocate(SC_flux_c (imin_:imax_,1:nbins_cond,1:3,1:nscalar))
  allocate(vel_triple_corr_c(imin_:imax_,1:nbins_cond,1:3,1:3))

  ! For conditional budgets -- yz average
  allocate(tmp_E_cond(imin_:imax_,1:nbins_cond,1:3))

  ! Conditional mean -- z average
  !  -- Increase the stencil width to 0:nbins_cond+1 for conditional gradients
  allocate(RHOm_cz   (imin_:imax_,jmid:jmax_,0:nbins_cond+1))
  allocate(NU_cz     (imin_:imax_,jmid:jmax_,0:nbins_cond+1))
  allocate(Um_cz     (imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3))
!!$  allocate(Um_cz     (imino_:imaxo_,jmid-1:jmaxo_,0:nbins_cond+1,1:3)); Um_cz = 0.0_WP
  allocate(SCm_cz    (imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:nscalar))
  allocate(DIFFm_cz  (imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:N_tot+1))
  allocate(src_SCm_cz(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:N_tot+1))
  allocate(TAUm_cz   (imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3,1:3)) ! i,j
  allocate(w_dot_cz  (imin_:imax_,jmid:jmax_,0:nbins_cond+1))
  allocate(chi_cz    (imin_:imax_,jmid:jmax_,0:nbins_cond+1))

  ! For temporal derivatives
  allocate(RHOm_cz_old(imin_:imax_,jmid:jmax_,0:nbins_cond+1))
  allocate(Um_cz_old  (imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3))
  RHOm_cz_old = 0.0_WP
  Um_cz_old = 0.0_WP

  ! Conditional mean gradients -- z average
  if (need_gradients) then
     allocate(dRHOdt_cz (imin_:imax_,jmid:jmax_,0:nbins_cond+1))
     allocate(dUdt_cz   (imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3))

     allocate(dPdxm_cz  (imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3))
     allocate(dUdxm_cz  (imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3,1:3)) ! i,j
     allocate(dUdxm_cz_reg(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3,1:2)) ! i,j
     allocate(dTAUdxm_cz(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3,1:3)) ! i,j
  end if

  ! Conditional correlations -- z average
  allocate(UiUjm_cz  (imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3,1:3))
  allocate(UiSCm_cz  (imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3,1:nscalar))
  allocate(UiUiUjm_cz(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3,1:3))

  ! Conditional velocity and scalar moments -- z average
  !    !! not implemented

  ! For conditional budgets -- z average
  allocate(tmp_SC_cz_1(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:nscalar))
  allocate(tmp_SC_cz_2(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:nscalar))
  allocate(tmp_U_cz_1(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3,1:3))
  allocate(tmp_U_cz_2(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3))
  allocate(tmp_U_cz_3(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3))
  allocate(tmp_U_cz_4(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3))
  allocate(tmp_U_cz_5(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3))
  allocate(tmp_E_cz_1(imin_:imax_,jmid:jmax_,0:nbins_cond+1,1:3))
  allocate(tmp_E_cz_2(imin_:imax_,jmid:jmax_,0:nbins_cond+1))
  allocate(tmp_E_cz_3(imin_:imax_,jmid:jmax_,0:nbins_cond+1))
  allocate(tmp_E_cz_4(imin_:imax_,jmid:jmax_,0:nbins_cond+1))

  ! Progress variable pdf
  allocate(pdf_C(imin_:imax_,1:nbins_cond))
  allocate(pdf_Cz(imin_:imax_,jmid:jmax_,0:nbins_cond+1))
  pdf_C  = 0.0_WP
  pdf_Cz = 0.0_WP

  return
end subroutine volumeStats_average_init


! ========================================================== !
! Compute the averages of all variables in space (physical   !
!    and conditional) and time                               !
! ========================================================== !
subroutine volumeStats_average_all
  use volumeStats_average
  use volumeStats_metric
  implicit none

  integer :: i, j, k, ii, jj, isc, ibin
  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_) :: tmpxyz, tmpxyz2, tmp_tau
  real(WP) :: dti, ntime2, ntime1i

  dti = 1.0_WP/dt
  ntime2  = real(ntime_curr-2,WP)
  ntime1i = 1.0_WP/real(ntime_curr-1,WP)

  ! Compute the conditional mapping
  call cond_map(SC(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:nscalar),y_cond(imin_:imax_,jmin_:jmax_,kmin_:kmax_))

  ! Compute the conditional pdf
!!  !$OMP PARALLEL DO PRIVATE(j,i,ii) REDUCTION(+:pdf_C,pdf_Cz)
!!$  do k=kmin_,kmax_
!!$     do j=jmin_,jmid-1
!!$        do i=imin_,imax_
!!$           ii = bin_index(i,j,k)
!!$           pdf_C (i,ii)   = pdf_C (i,ii)   + 1.0_WP
!!$           pdf_Cz(i,jmax_-j-1,ii) = pdf_Cz(i,jmax_-j-1,ii) + 1.0_WP
!!$        end do
!!$     end do
!!$  end do
  ! Over all ny
  do k=kmin_+stz,kmax_-stz
     do j=jmin_,jmax_
        do i=imin_,imax_
           ii = bin_index(i,j,k)
           pdf_C (i,ii)   = pdf_C (i,ii) + 1.0_WP
        end do
     end do
  end do
  ! From jmid to ny -- density-weighted
  do k=kmin_+stz,kmax_-stz
     do j=jmid,jmax_
        do i=imin_,imax_
           ii = bin_index(i,j,k)
           pdf_Cz(i,j,ii) = pdf_Cz(i,j,ii) + RHO(i,j,k)
        end do
     end do
  end do
!!  !$OMP END PARALLEL DO

  ! Update terms in physical mean velocity budgets and unsteady terms
  call volumeStats_velocity_compute

  ! Mean progress variable
  call phys_avg_z(RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_)*y_cond, y_condm)

  ! Density
  call phys_avg_z (RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_), RHOm)
  call cond_avg_yz(RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_), RHOm_c)
  call cond_avg_z (RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_), RHOm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond))
  ! Unsteady density for conditional budgets
  if (ntime_curr.gt.1) then
     !$OMP PARALLEL DO PRIVATE(i,j)
     do ibin=1,nbins_cond
        do j=jmid,jmax_
           do i=imin_,imax_
              dRHOdt_cz(i,j,ibin) = (dRHOdt_cz(i,j,ibin)*ntime2 + (RHOm_cz(i,j,ibin) - RHOm_cz_old(i,j,ibin))*dti)*ntime1i
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  ! Save old conditional mean for next time step
  RHOm_cz_old = RHOm_cz

  ! Viscosity
  call phys_avg_z (VISC(imin_:imax_,jmin_:jmax_,kmin_:kmax_), NU)
  call cond_avg_yz(VISC(imin_:imax_,jmin_:jmax_,kmin_:kmax_), NU_c)
  call cond_avg_z (VISC(imin_:imax_,jmin_:jmax_,kmin_:kmax_), NU_cz(imin_:imax_,jmid:jmax_,1:nbins_cond))

  ! Mean density and velocity components at the cell faces
  tmpxyz = 0.0_WP
  !$OMP PARALLEL DO PRIVATE(j,i)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           tmpxyz(i,j,k) = sum(interp_sc_x(i,j,:)*RHO(i-st2:i+st1,j,k))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call phys_avg_z(tmpxyz, RHOm_x)
  tmpxyz = tmpxyz * U(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1)
  call phys_avg_z(tmpxyz, Um(:,:,1))

  !$OMP PARALLEL DO PRIVATE(j,i)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           tmpxyz(i,j,k) = sum(interp_sc_y(i,j,:)*RHO(i,j-st2:j+st1,k))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call phys_avg_z(tmpxyz, RHOm_y)
  tmpxyz = tmpxyz * U(imin_:imax_,jmin_:jmax_,kmin_:kmax_,2)
  call phys_avg_z(tmpxyz, Um(:,:,2))

  !$OMP PARALLEL DO PRIVATE(j,i)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           tmpxyz(i,j,k) = sum(interp_sc_z(i,j,:)*RHO(i,j,k-st2:k+st1))
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call phys_avg_z(tmpxyz, RHOm_z)
  tmpxyz = tmpxyz * U(imin_:imax_,jmin_:jmax_,kmin_:kmax_,3)
  call phys_avg_z(tmpxyz, Um(:,:,3))

  ! Mean velocity components at the cell centers
  do ii=1,3
     tmpxyz = RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_) * Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)
     call phys_avg_z (tmpxyz, Umi (:,:,ii))
     call cond_avg_yz(tmpxyz, Um_c(:,:,ii)) !! update these?
     call cond_avg_z (tmpxyz, Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii))
     ! Unsteady velocity for conditional budgets
     if (ntime_curr.gt.1) then
        !$OMP PARALLEL DO PRIVATE(i,j)
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 dUdt_cz(i,j,ibin,ii) = (dUdt_cz(i,j,ibin,ii)*ntime2 + (Um_cz(i,j,ibin,ii) - Um_cz_old(i,j,ibin,ii))*dti)*ntime1i
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end if
     ! Save old conditional mean for next time step
     Um_cz_old(:,:,:,ii) = Um_cz(:,:,:,ii)
  end do

  if (combust) then
     ! Scalars
     do isc=1,nscalar
        tmpxyz = RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_) * SC(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc)
        call phys_avg_z (tmpxyz, SCm  (:,:,isc))
        call cond_avg_yz(tmpxyz, SCm_c(:,:,isc))
        call cond_avg_z (tmpxyz, SCm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,isc))
     end do

     ! Diffusivity
     do isc=1,N_tot+1
        tmpxyz = DIFF(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc)
        call phys_avg_z (tmpxyz, DIFFm(:,:,isc))
        call cond_avg_yz(tmpxyz, DIFFm_c(:,:,isc))
        call cond_avg_z (tmpxyz, DIFFm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,isc))
     end do

     ! Chemical source terms
     do isc=1,N_tot+1
        tmpxyz = src_SC(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc)
        call phys_avg_z (tmpxyz, src_SCm  (:,:,isc))
        call cond_avg_yz(tmpxyz, src_SCm_c(:,:,isc))
        call cond_avg_z (tmpxyz, src_SCm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,isc))
     end do

     ! Progress variable volumetric source term
     call cond_avg_z(src_SC(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc_PROG)/(prog_lower-prog_upper), w_dot_cz(imin_:imax_,jmid:jmax_,1:nbins_cond))

     ! Velocity -- scalar moments
     do isc=1,nscalar
        do ii=1,3
           tmpxyz = RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_) * &
                Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii) * &
                SC(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc)
           call phys_avg_z (tmpxyz, UiSCm  (:,:,ii,isc))
           call cond_avg_yz(tmpxyz, UiSCm_c(:,:,ii,isc))
           call cond_avg_z (tmpxyz, UiSCm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,isc))
        end do
     end do
  end if

  ! Velocity -- velocity moments
  do jj=1,3
     do ii=1,3
        tmpxyz = RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_) * &
             Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii) * &
             Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj)
        call phys_avg_z (tmpxyz, UiUjm  (:,:,ii,jj))
        call cond_avg_yz(tmpxyz, UiUjm_c(:,:,ii,jj))
        call cond_avg_z (tmpxyz, UiUjm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,jj))
     end do
  end do

  do jj=1,3
     do ii=1,3
        tmpxyz = RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_) * &
             Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii) * &
             Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii) * &
             Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj)
        call phys_avg_z (tmpxyz, UiUiUjm  (:,:,ii,jj))
        call cond_avg_yz(tmpxyz, UiUiUjm_c(:,:,ii,jj))
        call cond_avg_z (tmpxyz, UiUiUjm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,jj))
     end do
  end do


  if (need_gradients) then
     ! Density gradient
     do ii=1,3
        tmpxyz = dRHOdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)
        call phys_avg_z (tmpxyz, dRHOdxm  (:,:,ii))
        call cond_avg_yz(tmpxyz, dRHOdxm_c(:,:,ii))
     end do

     ! For z-gradients of Favre-averaged velocity
     do jj=3,3 ! since we take the divergence manually for non-z components
        do ii=1,3
           tmpxyz = ( & 
                + RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_)*     dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,jj) &
                + Ui (imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)*dRHOdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj) )
           call phys_avg_z(tmpxyz, dUdx_tmpv(:,:,ii,jj))
        end do
     end do

     ! For turbulent stress term in momentum budgets
     ! dUdx are all at cell centers
     do jj=1,3
        do ii=1,3
           tmpxyz = RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_) * ( & 
                + Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)*dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj,jj) &
                + Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj)*dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,jj) ) &
                + Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)*Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj)*dRHOdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj)
           call phys_avg_z(tmpxyz, dRHOUiUjdxm(:,:,ii,jj))
        end do
     end do

     ! For convection term in TKE budgets
     tmpxyz2 = 0.0_WP
     do ii=1,3
        tmpxyz2 = tmpxyz2 + Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)*Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)
     end do
     ! dUdx are all at cell centers
     do jj=1,3
        tmpxyz = 0.0_WP
        do ii=1,3
           tmpxyz = tmpxyz + RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_) * ( & 
                + 2.0_WP*Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)*dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,jj) ) &
                + Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)**2 * dRHOdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj)
        end do
        call phys_avg_z(tmpxyz, dRHOUiUidxm(:,:,jj))

        ! For turbulent transport term
        tmpxyz = tmpxyz*Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj) &
             + RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_) * tmpxyz2 * dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj,jj)
        call phys_avg_z(tmpxyz, dRHOUiUiUjdxm(:,:,jj))
     end do

     if (combust) then
        do ii=1,3
           ! Account for definition of progress variable
           ! Finalize later
           call cond_avg_yz(RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_) * &
                dSCdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,isc_PROG)/(prog_lower-prog_upper), dRHOdxm_cc(:,:,ii))
        end do
     end if

     ! Stress tensor gradients
     do jj=1,3
        do ii=1,3
           tmpxyz = dTAUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,jj)
           call phys_avg_z (tmpxyz, dTAUdxm  (:,:,ii,jj))
           call cond_avg_yz(tmpxyz, dTAUdxm_c(:,:,ii,jj))
           call cond_avg_z (tmpxyz, dTAUdxm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,jj))
        end do
     end do

     ! Pressure gradient
     do jj=1,3
        tmpxyz = dPdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj)
        call phys_avg_z (tmpxyz, dPdxm  (:,:,jj))
        call cond_avg_yz(tmpxyz, dPdxm_c(:,:,jj))
        call cond_avg_z (tmpxyz, dPdxm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,jj))
     end do

     if (combust) then
        ! Scalar gradients

        ! Scalar dissipation rate
        tmpxyz = 0.0_WP
        do jj=1,3
           tmpxyz = tmpxyz + &
                ( dSCdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj,isc_PROG)/(prog_lower-prog_upper) )**2
        end do
        ! DIFF is RHO*D
        tmpxyz = tmpxyz * 2.0_WP*DIFF(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc_PROG)
        call cond_avg_z(tmpxyz, chi_cz(imin_:imax_,jmid:jmax_,1:nbins_cond))
        do isc=1,nscalar
           tmp_tau = tmpxyz * SC(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc)
           call cond_avg_z(tmp_tau, tmp_SC_cz_2(imin_:imax_,jmid:jmax_,1:nbins_cond,isc))
        end do
        do ii=1,2
           tmp_tau = tmpxyz * Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)
           call cond_avg_z(tmp_tau, tmp_U_cz_5(imin_:imax_,jmid:jmax_,1:nbins_cond,ii))
        end do
     end if

     ! For budgets
     tmpxyz = 0.0_WP
     do ii=1,3
        tmpxyz = tmpxyz + &
             Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii) * &
             dPdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)
     end do
     call phys_avg_z(tmpxyz, tmp_E_phys(:,:,1))
     call cond_avg_yz(tmpxyz, tmp_E_cond(:,:,1))

     tmpxyz = 0.0_WP
     do jj=1,3
        do ii=1,3
           ! Compute TAU and average
           if (ii.eq.jj) then
              tmp_tau = 2.0_WP*VISC(imin_:imax_,jmin_:jmax_,kmin_:kmax_) * &
                   ( dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,ii) &
                   - 1.0_WP/3.0_WP * &
                   ( dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1,1) &
                   + dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,2,2) &
                   + dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,3,3) ) )
           else
              tmp_tau = VISC(imin_:imax_,jmin_:jmax_,kmin_:kmax_) * &
                   ( dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,jj) &
                   + dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj,ii) )
           end if
           call phys_avg_z (tmp_tau, TAUm   (:,:,ii,jj))
           call cond_avg_yz(tmp_tau, TAUm_c (:,:,ii,jj))
           call cond_avg_z (tmp_tau, TAUm_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,jj))

           ! Component for budget
           tmpxyz = tmpxyz + dUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,jj) * tmp_tau
        end do
     end do
     call phys_avg_z(tmpxyz, tmp_E_phys(:,:,3))
     call cond_avg_yz(tmpxyz, tmp_E_cond(:,:,3))
     do jj=1,3
        do ii=1,3
           tmpxyz = tmpxyz &
                + Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii) * dTAUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,jj)
        end do
     end do
     call phys_avg_z(tmpxyz, tmp_E_phys(:,:,2))
     call cond_avg_yz(tmpxyz, tmp_E_cond(:,:,2))


     ! For conditional scalar budgets -- z average
     if (combust) then
        ! tmp1
        do isc=1,nscalar
           tmpxyz = 0.0_WP
           do jj=1,3
              tmpxyz = tmpxyz + &
                   DIFF (imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc_PROG) * &
                   dSCdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj,isc) * &
                   dSCdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj,isc_PROG)/(prog_lower-prog_upper)
           end do
           call cond_avg_z(tmpxyz, tmp_SC_cz_1(imin_:imax_,jmid:jmax_,1:nbins_cond,isc))
        end do

        ! For conditional velocity budgets -- z average
        ! tmp1
        tmpxyz = 0.0_WP
        do jj=1,3
           do ii=1,3
              tmpxyz = Ui    (imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii) * &
                   DIFF (imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc_PROG) * &
                   dSCdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj,isc_PROG)/(prog_lower-prog_upper)
              call cond_avg_z(tmpxyz, tmp_U_cz_1(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,jj))
           end do
        end do

        ! tmp2
        do jj=1,3
           tmpxyz = DIFF (imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc_PROG) * &
                dSCdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj,isc_PROG)/(prog_lower-prog_upper)
           call cond_avg_z(tmpxyz, tmp_U_cz_2(imin_:imax_,jmid:jmax_,1:nbins_cond,jj))
        end do

        ! tmp3 -- src_SC is mass source term; need volumetric here
        do jj=1,3
           tmpxyz = Ui     (imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj) * &
                src_SC(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc_PROG)/(prog_lower-prog_upper)
           call cond_avg_z(tmpxyz, tmp_U_cz_3(imin_:imax_,jmid:jmax_,1:nbins_cond,jj))
        end do

        ! tmp 4
        do ii=1,2
           tmpxyz = 0.0_WP
           do jj=1,3
              tmpxyz = tmpxyz + &
                   DIFF (imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc_PROG) * &
                   dUdx (imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,jj) * &
                   dSCdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj,isc_PROG)/(prog_lower-prog_upper)
           end do
           call cond_avg_z(tmpxyz, tmp_U_cz_4(imin_:imax_,jmid:jmax_,1:nbins_cond,ii))
        end do

        ! For conditional energy budgets -- z average
        ! tmp_E_cz_1
        ! tmpxyz is uiui
        tmpxyz = 0.0_WP
        do ii=1,3
           tmpxyz = tmpxyz + &
                Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii) * &
                Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)
        end do
        do jj=1,3
           tmp_tau = tmpxyz * &
                DIFF (imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc_PROG) * &
                dSCdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj,isc_PROG)/(prog_lower-prog_upper)
           call cond_avg_z(tmp_tau, tmp_E_cz_1(imin_:imax_,jmid:jmax_,1:nbins_cond,jj))
        end do

        ! tmp_E_cz_2
        tmp_tau = tmpxyz * &
             src_SC(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc_PROG)/(prog_lower-prog_upper)
        call cond_avg_z(tmp_tau, tmp_E_cz_2(imin_:imax_,jmid:jmax_,1:nbins_cond))

        ! tmp_E_cz_3
        tmp_tau = 0.0_WP
        do ii=1,3
           tmp_tau = tmp_tau + &
                Ui   (imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii) * &
                dPdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)
        end do
        call cond_avg_z(tmp_tau, tmp_E_cz_3(imin_:imax_,jmid:jmax_,1:nbins_cond))

        ! tmp_E_cz_4
        tmp_tau = 0.0_WP
        do ii=1,3
           do jj=1,3
              tmp_tau = tmp_tau + &
                   Ui    (imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii) * &
                   dTAUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,jj)
           end do
        end do
        call cond_avg_z(tmp_tau, tmp_E_cz_4(imin_:imax_,jmid:jmax_,1:nbins_cond))
        
     end if
  end if

  return
end subroutine volumeStats_average_all


subroutine volumeStats_average_finalize
  use volumeStats_average
  use volumeStats_metric
  implicit none

  integer :: i, j, k, ii, jj, isc, ibin
  integer :: iplane, iunit
  character(len=str_short) :: of_name, indx_name
  real(WP), dimension(imin_:imax_,1:nbins_cond) :: tmpxc, dPDFdc
  real(WP) :: tmp, w1, w2
  real(WP) :: r_half, U_cl

  ! Finalize the pdfs
  pdf_C  = pdf_C *real(nbins_cond,WP)/real(ny*(nz-2*stz)*ntime_curr,WP) ! be careful if we reduce the number of nz in average
  !$OMP PARALLEL DO
  do ibin=1,nbins_cond
     do j=jmid,jmax_
        do i=imin_,imax_
           pdf_Cz(i,j,ibin) = pdf_Cz(i,j,ibin)*real(nbins_cond,WP) &
                / (real((nz-2*stz)*ntime_curr,WP) * RHOm_cz(i,j,ibin))
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  ! Finalize the unsteady velocity terms
  ! Um_cz is still RHO*Ui
  do ii=1,3
     !$OMP PARALLEL DO
     do ibin=1,nbins_cond
        do j=jmid,jmax_
           do i=imin_,imax_
              dUdt_cz(i,j,ibin,ii) = dUdt_cz(i,j,ibin,ii)/RHOm_cz(i,j,ibin) &
                   - Um_cz(i,j,ibin,ii)/RHOm_cz(i,j,ibin)**2 * dRHOdt_cz(i,j,ibin)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end do

  ! Complete the density weighting
  !$OMP PARALLEL
  !$OMP DO
  do j=jmin_,jmax_
     do i=imin_,imax_
        y_condm(i,j) = y_condm(i,j)/RHOm(i,j)
     end do
  end do
  !$OMP END DO

  ! Viscosity in NGA is RHO*VISC
  !$OMP DO
  do j=jmin_,jmax_
     do i=imin_,imax_
        NU(i,j) = NU(i,j)/RHOm(i,j)
     end do
  end do
  !$OMP END DO
  !$OMP DO
  do j=1,nbins_cond
     do i=imin_,imax_
        NU_c(i,j) = NU_c(i,j)/RHOm_c(i,j)
     end do
  end do
  !$OMP END DO
  !$OMP DO
  do ibin=1,nbins_cond
     do j=jmid,jmax_
        do i=imin_,imax_
           NU_cz(i,j,ibin) = NU_cz(i,j,ibin)/RHOm_cz(i,j,ibin)
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Velocity components at cell faces
  do j=jmin_,jmax_
     do i=imin_,imax_
        Um(i,j,1) = Um(i,j,1)/RHOm_x(i,j)
        Um(i,j,2) = Um(i,j,2)/RHOm_y(i,j)
        Um(i,j,3) = Um(i,j,3)/RHOm_z(i,j)
     end do
  end do

  ! Velocity components at cell centers
  do jj=1,3
     !$OMP PARALLEL DO
     do j=jmin_,jmax_
        do i=imin_,imax_
           Umi(i,j,jj) = Umi(i,j,jj)/RHOm(i,j)
        end do
     end do
     !$OMP END PARALLEL DO
  end do
  do jj=1,3
     !$OMP PARALLEL DO
     do ibin=1,nbins_cond
        do i=imin_,imax_
           Um_c(i,ibin,jj) = Um_c(i,ibin,jj)/RHOm_c(i,ibin)
        end do
     end do
     !$OMP END PARALLEL DO
  end do
  do jj=1,3
     !$OMP PARALLEL DO
     do ibin=1,nbins_cond
        do j=jmid,jmax_
           do i=imin_,imax_
              Um_cz(i,j,ibin,jj) = Um_cz(i,j,ibin,jj)/RHOm_cz(i,j,ibin)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end do

  ! Scalars
  if (combust) then
     do isc=1,N_tot+1
        !$OMP PARALLEL DO
        do j=jmin_,jmax_
           do i=imin_,imax_
              SCm(i,j,isc) = SCm(i,j,isc)/RHOm(i,j)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
     do isc=1,N_tot+1
        !$OMP PARALLEL DO
        do ibin=1,nbins_cond
           do i=imin_,imax_
              SCm_c(i,ibin,isc) = SCm_c(i,ibin,isc)/RHOm_c(i,ibin)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
     do isc=1,N_tot+1
        !$OMP PARALLEL DO
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 SCm_cz(i,j,ibin,isc) = SCm_cz(i,j,ibin,isc)/RHOm_cz(i,j,ibin)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end if

  ! Velocity moments
  do jj=1,3
     do ii=1,3
        !$OMP PARALLEL DO
        do j=jmin_,jmax_
           do i=imin_,imax_
              UiUjm(i,j,ii,jj) = UiUjm(i,j,ii,jj)/RHOm(i,j)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end do
  do jj=1,3
     do ii=1,3
        !$OMP PARALLEL DO
        do ibin=1,nbins_cond
           do i=imin_,imax_
              UiUjm_c(i,ibin,ii,jj) = UiUjm_c(i,ibin,ii,jj)/RHOm_c(i,ibin)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end do
  do jj=1,3
     do ii=1,3
        !$OMP PARALLEL DO
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 UiUjm_cz(i,j,ibin,ii,jj) = UiUjm_cz(i,j,ibin,ii,jj)/RHOm_cz(i,j,ibin)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end do
  do jj=1,3
     do ii=1,3
        !$OMP PARALLEL DO
        do j=jmin_,jmax_
           do i=imin_,imax_
              UiUiUjm(i,j,ii,jj) = UiUiUjm(i,j,ii,jj)/RHOm(i,j)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end do
  do jj=1,3
     do ii=1,3
        !$OMP PARALLEL DO
        do ibin=1,nbins_cond
           do i=imin_,imax_
              UiUiUjm_c(i,ibin,ii,jj) = UiUiUjm_c(i,ibin,ii,jj)/RHOm_c(i,ibin)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end do
  do jj=1,3
     do ii=1,3
        !$OMP PARALLEL DO
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 UiUiUjm_cz(i,j,ibin,ii,jj) = UiUiUjm_cz(i,j,ibin,ii,jj)/RHOm_cz(i,j,ibin)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end do

  ! Scalar-velocity moments
  if (combust) then
     do isc=1,nscalar
        do ii=1,3
           !$OMP PARALLEL DO
           do j=jmin_,jmax_
              do i=imin_,imax_
                 UiSCm(i,j,ii,isc) = UiSCm(i,j,ii,isc)/RHOm(i,j)
              end do
           end do
           !$OMP END PARALLEL DO
        end do
     end do
     do isc=1,nscalar
        do ii=1,3
           !$OMP PARALLEL DO
           do ibin=1,nbins_cond
              do i=imin_,imax_
                 UiSCm_c(i,ibin,ii,isc) = UiSCm_c(i,ibin,ii,isc)/RHOm_c(i,ibin)
              end do
           end do
           !$OMP END PARALLEL DO
        end do
     end do
     do isc=1,nscalar
        do ii=1,3
           !$OMP PARALLEL DO
           do ibin=1,nbins_cond
              do j=jmid,jmax_
                 do i=imin_,imax_
                    UiSCm_cz(i,j,ibin,ii,isc) = UiSCm_cz(i,j,ibin,ii,isc)/RHOm_cz(i,j,ibin)
                 end do
              end do
           end do
           !$OMP END PARALLEL DO
        end do
     end do

     ! Scalar flux
     do isc=1,nscalar
        do ii=1,3
           !$OMP PARALLEL DO
           do j=jmin_,jmax_
              do i=imin_,imax_
                 SC_flux(i,j,ii,isc) = UiSCm(i,j,ii,isc) - Umi(i,j,ii)*SCm(i,j,isc)
              end do
           end do
           !$OMP END PARALLEL DO
        end do
     end do
     do isc=1,nscalar
        do ii=1,3
           !$OMP PARALLEL DO
           do ibin=1,nbins_cond
              do i=imin_,imax_
                 SC_flux_c(i,ibin,ii,isc) = UiSCm_c(i,ibin,ii,isc) - Um_c(i,ibin,ii)*SCm_c(i,ibin,isc)
              end do
           end do
           !$OMP END PARALLEL DO
        end do
     end do

     ! Diffusivity in NGA is RHO*DIFF
     do isc=1,N_tot+1
        !$OMP PARALLEL DO
        do j=jmin_,jmax_
           do i=imin_,imax_
              DIFFm(i,j,isc) = DIFFm(i,j,isc)/RHOm(i,j)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
     do isc=1,N_tot+1
        !$OMP PARALLEL DO
        do j=1,nbins_cond
           do i=imin_,imax_
              DIFFm_c(i,j,isc) = DIFFm_c(i,j,isc)/RHOm_c(i,j)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end if

  ! Reynolds stress
  do jj=1,3
     do ii=1,3
        !$OMP PARALLEL DO
        do j=jmin_,jmax_
           do i=imin_,imax_
              R_stress(i,j,ii,jj) = UiUjm(i,j,ii,jj) - Umi(i,j,ii)*Umi(i,j,jj)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end do
  do jj=1,3
     do ii=1,3
        !$OMP PARALLEL DO
        do j=1,nbins_cond
           do i=imin_,imax_
              R_stress_c(i,j,ii,jj) = UiUjm_c(i,j,ii,jj) - Um_c(i,j,ii)*Um_c(i,j,jj)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end do

  ! Velocity triple correlation
  do jj=1,3
     do ii=1,3
        !$OMP PARALLEL DO
        do j=jmin_,jmax_
           do i=imin_,imax_
              vel_triple_corr(i,j,ii,jj) = &
                   + UiUiUjm(i,j,ii,jj) &
                   - UiUjm  (i,j,ii,ii)*Umi(i,j,jj) &
                   - 2.0_WP*UiUjm(i,j,ii,jj)*Umi(i,j,ii) &
                   + 2.0_WP*Umi(i,j,ii)*Umi(i,j,ii)*Umi(i,j,jj)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end do
  do jj=1,3
     do ii=1,3
        !$OMP PARALLEL DO
        do j=1,nbins_cond
           do i=imin_,imax_
              vel_triple_corr_c(i,j,ii,jj) = &
                   + UiUiUjm_c(i,j,ii,jj) &
                   - UiUjm_c  (i,j,ii,ii)*Um_c(i,j,jj) &
                   - 2.0_WP*UiUjm_c(i,j,ii,jj)*Um_c(i,j,ii) &
                   + 2.0_WP*Um_c(i,j,ii)*Um_c(i,j,ii)*Um_c(i,j,jj)
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end do

  ! Mean velocity gradients
!!$  do j=jmin_,jmax_
!!$     do i=imin_,imax_
!!$        Um(i,j,1) = U(i,j,20,1)
!!$        Um(i,j,2) = U(i,j,20,2)
!!$        Um(i,j,3) = U(i,j,20,3)
!!$        dUdxm(i,j,1,1) = sum(div_u(i,j,:)*U(i-st1:i+st2,j,20,1)) !at cell centers
!!$        dUdxm(i,j,2,2) = sum(div_v(i,j,:)*U(i,j-st1:j+st2,20,2))
!!$        dUdxm(i,j,3,3) = sum(div_w(i,j,:)*U(i,j,20-st1:20+st2,3))
!!$     end do
!!$  end do

  ! Mean scalar gradients
  !  -- Compute this manually for the density weighting
  dSCdxm = 0.0_WP
  do isc=1,nscalar
     do j=jmin_+st2,jmax_-st2
        do i=imin_+st2,imax_-st2
           ! dSC/dx
           dSCdxm(i,j,1,isc) = sum(grad_x(i,j,:)*SCm(i-st2:i+st1,j,isc))
           ! dSC/dy
           dSCdxm(i,j,2,isc) = sum(grad_y(i,j,:)*SCm(i,j-st2:j+st1,isc))
        end do
     end do
  end do

  ! Mean velocity gradients
  !  -- Compute this manually for the density weighting
  dUdxm = 0.0_WP
  do j=jmin_+st2,jmax_-st2
     do i=imin_+st2,imax_-st2
        ! dU/dx
        dUdxm(i,j,1,1) = sum(div_u(i,j,:)*Um(i-st1:i+st2,j,1))
        ! dU/dy
        dUdxm(i,j,1,2) = dyi(j)*( +sum(interp_sc_y(i,j+1,:)*Umi(i,j-st2+1:j+st1+1,1)) &
                                  -sum(interp_sc_y(i,j,  :)*Umi(i,j-st2  :j+st1  ,1)) )
        ! dU/dz
        dUdxm(i,j,1,3) = dUdx_tmpv(i,j,1,3)/RHOm(i,j) - Umi(i,j,1)*dRHOdxm(i,j,3)/RHOm(i,j)
        ! dV/dx
        dUdxm(i,j,2,1) = dxi(i)*( +sum(interp_sc_x(i+1,j,:)*Umi(i-st2+1:i+st1+1,j,2)) &
                                  -sum(interp_sc_x(i  ,j,:)*Umi(i-st2  :i+st1  ,j,2)) )
        ! dV/dy
        dUdxm(i,j,2,2) = sum(div_v(i,j,:)*Um(i,j-st1:j+st2,2))
        ! dV/dz
        dUdxm(i,j,2,3) = dUdx_tmpv(i,j,2,3)/RHOm(i,j) - Umi(i,j,2)*dRHOdxm(i,j,3)/RHOm(i,j)
        ! dW/dx
        dUdxm(i,j,3,1) = dxi(i)*( +sum(interp_sc_x(i+1,j,:)*Umi(i-st2+1:i+st1+1,j,3)) &
                                  -sum(interp_sc_x(i  ,j,:)*Umi(i-st2  :i+st1  ,j,3)) )
        ! dW/dy
        dUdxm(i,j,3,2) = dyi(j)*( +sum(interp_sc_y(i,j+1,:)*Umi(i,j-st2+1:j+st1+1,3)) &
                                  -sum(interp_sc_y(i,j  ,:)*Umi(i,j-st2  :j+st1  ,3)) )
        ! dW/dz
        dUdxm(i,j,3,3) = dUdx_tmpv(i,j,3,3)/RHOm(i,j) - Umi(i,j,3)*dRHOdxm(i,j,3)/RHOm(i,j)
     end do
  end do

  do ii=1,3    !!!! THESE CAN BE BETTER
     do ibin=1,nbins_cond
        do j=jmid,jmax_
           do i=imin_,imax_
              ! dUmi/dx
              dUdxm_cz(i,j,ibin,ii,1) = dxi(i) * &
                   ( sum(interp_sc_x(i+1,j,:)*Um_cz(i-st2+1:i+st1+1,j,ibin,ii)) &
                   - sum(interp_sc_x(i  ,j,:)*Um_cz(i-st2  :i+st1  ,j,ibin,ii)) )
           end do
        end do
     end do
     if (ii.eq.2) then
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 ! dUmi/dy
                 dUdxm_cz(i,j,ibin,ii,2) = dyi(j) * &
                      ( sum(interp_vm_y(i,j+1,:)*Um_cz(i,j-st2+1:j+st1+1,ibin,ii)) &
                      - sum(interp_vm_y(i,j  ,:)*Um_cz(i,j-st2  :j+st1  ,ibin,ii)) )
              end do
           end do
        end do
     else
        do ibin=1,nbins_cond
           do j=jmid,jmax_
              do i=imin_,imax_
                 ! dUmi/dy
                 dUdxm_cz(i,j,ibin,ii,2) = dyi(j) * &
                      ( sum(interp_sc_y(i,j+1,:)*Um_cz(i,j-st2+1:j+st1+1,ibin,ii)) &
                      - sum(interp_sc_y(i,j  ,:)*Um_cz(i,j-st2  :j+st1  ,ibin,ii)) )
              end do
           end do
        end do
     end if
     ! dUmi/dz
     dUdxm_cz(:,:,ibin,ii,3) = 0.0_WP
  end do

  ! Conditional mean gradients from regularized derivatives
  do jj=1,2
     do ii=1,3
        if (irank.eq.iroot) print *, 'reg_deriv', ii,jj
        call reg_deriv(Um_cz(imin_:imax_,jmid:jmax_,1:nbins_cond,ii), &
             dUdxm_cz_reg(imin_:imax_,jmid:jmax_,1:nbins_cond,ii,jj), jj)
     end do
  end do
     

  ! Conditional mean gradients WITH pdf
  do ii=1,3
     tmpxc = dRHOdxm_cc(:,:,ii)
     do j=1,nbins_cond
        do i=imin_,imax_
           ! Gradient in sample space
           tmp = 0.5_WP*( tmpxc(i,j+1)*pdf_C(i,j+1) - tmpxc(i,j-1)*pdf_C(i,j-1) )/real(nbins_cond,WP)
           ! Neglects gradient of pdf
           dRHOdxm_cc(i,j,ii) = dRHOdxm_c(i,j,ii) - tmp/pdf_C(i,j)
        end do
     end do
  end do

  ! Gradient of pdf in sample space
  do j=1,nbins_cond
     do i=imin_,imax_
        dPDFdc(i,j) = 0.5_WP*( pdf_C(i,j+1) - pdf_C(i,j-1) )/real(nbins_cond,WP)
     end do
  end do

  ! Mean stress tensor gradients
  do ii=1,3
     do j=jmin_,jmax_
        do i=imin_,imax_
           ! dTAUm/dx
           dTAUdxm_2cd(i,j,ii,1) = dxi(i) * &
                ( sum(interp_sc_x(i+1,j,:)*TAUm(i-st2+1:i+st1+1,j,ii,1)) &
                - sum(interp_sc_x(i  ,j,:)*TAUm(i-st2  :i+st1  ,j,ii,1)) )
        end do
     end do
     if (ii.eq.2) then
        do j=jmin_,jmax_
           do i=imin_,imax_
              ! dTAUm/dy
              dTAUdxm_2cd(i,j,ii,2) = dyi(j) * &
                   ( sum(interp_vm_y(i,j+1,:)*TAUm(i,j-st2+1:j+st1+1,ii,2)) &
                   - sum(interp_vm_y(i,j  ,:)*TAUm(i,j-st2  :j+st1  ,ii,2)) )
           end do
        end do
     else
        do j=jmin_,jmax_
           do i=imin_,imax_
              ! dTAUm/dy
              dTAUdxm_2cd(i,j,ii,2) = dyi(j) * &
                   ( sum(interp_sc_y(i,j+1,:)*TAUm(i,j-st2+1:j+st1+1,ii,2)) &
                   - sum(interp_sc_y(i,j  ,:)*TAUm(i,j-st2  :j+st1  ,ii,2)) )
           end do
        end do
     end if
     ! dTAUm/dz
     dTAUdxm_2cd(i,j,ii,3) = 0.0_WP
  end do


  ! Compute centerline axial velocity, jet half width, and turbulence stats
  if (irank.eq.iroot) then
     print '(a12,ES15.5,a2)', 'Total time: ', time(ntime_curr)-time(1), ' s'
     print       '(3a15)' ,'x (m)', 'r1/2 (m)', 'U_cl (m/s)'!, 'u_rms (m/s)', 'TKE (m^2/s^2)', 'L_cl (mm)', 'tau_L (1/s)', 'Re_L', 'L times', 'eta_cl (um)', 'tau_eta (1/s)'
  end if

  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     i = pnindx(iplane)
     open(unit=iunit,  file=trim(output_name)//'_centerline_stats_'//trim(of_name), action='write')
     write(iunit,'(3a15)'),'x (m)', 'r1/2 (m)', 'U_cl (m/s)'!, 'u_rms (m/s)', 'TKE (m^2/s^2)', 'L_cl (mm)', 'tau_L (1/s)', 'Re_L', 'L times', 'eta_cl (um)', 'tau_eta (1/s)'
     
     r_half = -1.0_WP
     U_cl = sum(interp_pnxm(iplane,:)*Um(i-1:i,jmid,1))

     ! March upward from the centerline to find the half-width
     do j=jmid,ny-1
        if (sum(interp_pnxm(iplane,:)*Um(i-1:i,j,1)).gt.0.5_WP*U_cl .and. &
             sum(interp_pnxm(iplane,:)*Um(i-1:i,j+1,1)).le.0.5_WP*U_cl) then
           w1 = (0.5_WP*U_cl - sum(interp_pnxm(iplane,:)*Um(i-1:i,j+1,1))) &
                / (sum(interp_pnxm(iplane,:)*Um(i-1:i,j+1,1)) - sum(interp_pnxm(iplane,:)*Um(i-1:i,j,1)))
           w2 = 1.0_WP - w1
           r_half = w1*ym(j) + w2*ym(j+1)
        end if
     end do

     print       '(3ES15.6)' , sum(interp_pnxm(iplane,:)*xm(i-1:i)), r_half, U_cl
     write(iunit,'(3ES15.6)'), sum(interp_pnxm(iplane,:)*xm(i-1:i)), r_half, U_cl
     close(iunit)

     ! Compute velocity and turbulence stats at centerline
     
  end do

  ! Write some output
  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_vol_data_'//trim(of_name), action='write')
     i = pnindx(iplane)
     do j=jmin_,jmax_
        write (iunit,'(ES22.13)',advance='no'), y(j) ! y-coordinate
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*y_condm(i-1:i,j))
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*RHOm(i-1:i,j))
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*src_SCm(i-1:i,j,isc_H2O))
        do jj=1,3 ! 5
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*Um(i-1:i,j,jj))
        end do
        do ii=1,3 ! 8
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*UiUjm(i-1:i,j,ii,jj))
           end do
        end do
        do ii=1,3 ! 17
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*R_stress(i-1:i,j,ii,jj))
           end do
        end do
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*NU(i-1:i,j)) !26
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*DIFFm(i-1:i,j,isc_H2O)) !27
        do ii=1,3 !28
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*vel_triple_corr(i-1:i,j,ii,jj))
           end do
        end do
        do isc=isc_H2O,isc_H2O !37
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SCm(i-1:i,j,isc))
        end do
        do isc=isc_H2O,isc_H2O !38
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*UiSCm(i-1:i,j,jj,isc))
           end do
        end do
        do isc=isc_H2O,isc_H2O !41
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SC_flux(i-1:i,j,jj,isc))
           end do
        end do
        do ii=1,3 !44
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dTAUdxm(i-1:i,j,ii,jj))
           end do
        end do
        do ii=1,3 !53
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dTAUdxm_2cd(i-1:i,j,ii,jj))
           end do
        end do
        do ii=1,3 !62
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*TAUm(i-1:i,j,ii,jj))
           end do
        end do
        do jj=1,3 !71
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dPdxm(i-1:i,j,jj))
        end do
        do isc=isc_OH,isc_H2O !74
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SCm(i-1:i,j,isc))
        end do
        do isc=isc_OH,isc_H2O !77
           do jj=1,3
              if (isc.eq.isc_H2O.and.jj.eq.3) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*SC_flux(i-1:i,j,jj,isc))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SC_flux(i-1:i,j,jj,isc))
              end if
           end do
        end do
        
     end do
     close(iunit)
  end do


  ! Mean scalar gradients -- can extend this for d/dz gradients as for velocity
  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_dSCdxm_'//trim(of_name), action='write')
     i = pnindx(iplane)
     do j=jmin_,jmax_
        write (iunit,'(ES22.13)',advance='no'), y(j) ! y-coordinate
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*y_condm(i-1:i,j))
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*RHOm(i-1:i,j))
        do isc=1,nscalar
           do jj=1,2
              if (isc.eq.nscalar .and. jj.eq.2) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*dSCdxm(i-1:i,j,jj,isc))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dSCdxm(i-1:i,j,jj,isc))
              end if
           end do
        end do
     end do
     close(iunit)
  end do


  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_vol_data_cond_'//trim(of_name), action='write')
     i = pnindx(iplane)
     do j=1,nbins_cond
        write (iunit,'(ES22.13)',advance='no'), bins_cond(j) ! conditional variable
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*RHOm_c(i-1:i,j)) ! placeholder
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*RHOm_c(i-1:i,j))
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*src_SCm_c(i-1:i,j,isc_H2O))
        do jj=1,3
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*Um_c(i-1:i,j,jj))
        end do
        do ii=1,3
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*UiUjm_c(i-1:i,j,ii,jj))
           end do
        end do
        do ii=1,3
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*R_stress_c(i-1:i,j,ii,jj))
           end do
        end do
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*NU_c(i-1:i,j))
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*DIFFm_c(i-1:i,j,isc_H2O))
        do ii=1,3
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*vel_triple_corr_c(i-1:i,j,ii,jj))
           end do
        end do
        do isc=isc_H2O,isc_H2O
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SCm_c(i-1:i,j,isc))
        end do
        do isc=isc_H2O,isc_H2O
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*UiSCm_c(i-1:i,j,jj,isc))
           end do
        end do
        do isc=isc_H2O,isc_H2O
           do jj=1,3
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SC_flux_c(i-1:i,j,jj,isc))
           end do
        end do
        do ii=1,3
           do jj=1,3
              if (ii.eq.3.and.jj.eq.3) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*dTAUdxm_c(i-1:i,j,ii,jj))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dTAUdxm_c(i-1:i,j,ii,jj))
              end if
           end do
        end do
     end do
     close(iunit)
  end do

  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_pdf_C_'//trim(of_name), action='write')
     i = pnindx(iplane)
     do j=1,nbins_cond
        write (iunit,'(ES22.13)',advance='no'), bins_cond(j) ! conditional variable
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*pdf_C(i-1:i,j))
        write (iunit,'(ES22.13)',advance='no'), pdf_C(i,j)  
        write (iunit,'(ES22.13)'), dPDFdc(i,j)  
     end do
     close(iunit)
  end do

  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_grad_phys_'//trim(of_name), action='write')
     i = pnindx(iplane)
     do j=jmin_,jmax_
        write (iunit,'(ES22.13)',advance='no'), y(j) ! y-coordinate
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*y_condm(i-1:i,j))
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*RHOm(i-1:i,j))
        do jj=1,3
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dRHOdxm(i-1:i,j,jj))
        end do
        do ii=1,3
           do jj=1,3
              if (ii.eq.3.and.jj.eq.3) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*dUdxm(i-1:i,j,ii,jj))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dUdxm(i-1:i,j,ii,jj))
              end if
           end do
        end do
     end do
     close(iunit)
  end do

  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_grad_cond_'//trim(of_name), action='write')
     i = pnindx(iplane)
     do j=1,nbins_cond
        write (iunit,'(ES22.13)',advance='no'), bins_cond(j) ! conditional variable
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*RHOm_c(i-1:i,j)) ! placeholder
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*RHOm_c(i-1:i,j))
        do jj=1,3
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dRHOdxm_c(i-1:i,j,jj))
        end do
        do jj=1,3
           if (jj.eq.3) then
              write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*dRHOdxm_cc(i-1:i,j,jj))
           else
              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dRHOdxm_cc(i-1:i,j,jj))
           end if
        end do
     end do
     close(iunit)
  end do

  if (combust) then

     do iplane=pnmin_,pnmax_
        call get_name(iplane,of_name)
        open(unit=iunit, file=trim(output_name)//'_pdf_Cz_'//trim(of_name), action='write')
        i = pnindx(iplane)
        write (iunit,'(ES22.13)',advance='no'), 0.0_WP
        do j=jmid,jmax_
           if (j.eq.jmax_) then
              write (iunit,'(ES22.13)'), ym(j)
           else
              write (iunit,'(ES22.13)',advance='no'), ym(j)
           end if
        end do
        do ibin=1,nbins_cond
           write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
           do j=jmid,jmax_
              if (j.eq.jmax_) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*pdf_Cz(i-1:i,j,ibin))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*pdf_Cz(i-1:i,j,ibin))
              end if
           end do
        end do
        close(iunit)
     end do

     do iplane=pnmin_,pnmax_
        call get_name(iplane,of_name)
        open(unit=iunit, file=trim(output_name)//'_RHOm_cz_'//trim(of_name), action='write')
        i = pnindx(iplane)
        do ibin=1,nbins_cond
           write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
           do j=jmid,jmax_
              if (j.eq.jmax_) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*RHOm_cz(i-1:i,j,ibin))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*RHOm_cz(i-1:i,j,ibin))
              end if
           end do
        end do
        close(iunit)
     end do

     do jj=1,3
        write(indx_name,'(I1)'), jj
        do iplane=pnmin_,pnmax_
           call get_name(iplane,of_name)
           open(unit=iunit, file=trim(output_name)//'_U'//trim(indx_name)//'m_cz_'//trim(of_name), action='write')
           i = pnindx(iplane)
           do ibin=1,nbins_cond
              write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
              do j=jmid,jmax_
                 if (j.eq.jmax_) then
                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*Um_cz(i-1:i,j,ibin,jj))
                 else
                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*Um_cz(i-1:i,j,ibin,jj))
                 end if
              end do
           end do
           close(iunit)
        end do
     end do

     do isc=1,nscalar
        do iplane=pnmin_,pnmax_
           call get_name(iplane,of_name)
           open(unit=iunit, file=trim(output_name)//'_SCm_cz_'//trim(SC_name(isc))//'_'//trim(of_name), action='write')
           i = pnindx(iplane)
           do ibin=1,nbins_cond
              write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
              do j=jmid,jmax_
                 if (j.eq.jmax_) then
                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*SCm_cz(i-1:i,j,ibin,isc))
                 else
                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SCm_cz(i-1:i,j,ibin,isc))
                 end if
              end do
           end do
           close(iunit)
        end do
     end do

     do jj=1,3
        do ii=1,3
           write(indx_name,'(2I1)'), ii,jj
           do iplane=pnmin_,pnmax_
              call get_name(iplane,of_name)
              open(unit=iunit, file=trim(output_name)//'_dUdxm_cz_'//trim(indx_name)//'_'//trim(of_name), action='write')
              i = pnindx(iplane)
              do ibin=1,nbins_cond
                 write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
                 do j=jmid,jmax_
                    if (j.eq.jmax_) then
                       write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*dUdxm_cz(i-1:i,j,ibin,ii,jj))
                    else
                       write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dUdxm_cz(i-1:i,j,ibin,ii,jj))
                    end if
                 end do
              end do
              close(iunit)
           end do
        end do
     end do

     do jj=1,2
        do ii=1,2
           write(indx_name,'(2I1)'), ii,jj
           do iplane=pnmin_,pnmax_
              call get_name(iplane,of_name)
              open(unit=iunit, file=trim(output_name)//'_dUdxm_cz_reg_'//trim(indx_name)//'_'//trim(of_name), action='write')
              i = pnindx(iplane)
              do ibin=1,nbins_cond
                 write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
                 do j=jmid,jmax_
                    if (j.eq.jmax_) then
                       write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*dUdxm_cz_reg(i-1:i,j,ibin,ii,jj))
                    else
                       write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dUdxm_cz_reg(i-1:i,j,ibin,ii,jj))
                    end if
                 end do
              end do
              close(iunit)
           end do
        end do
     end do

     do jj=1,3
        do ii=1,3
           write(indx_name,'(2I1)'), ii,jj
           do iplane=pnmin_,pnmax_
              call get_name(iplane,of_name)
              open(unit=iunit, file=trim(output_name)//'_UiUj'//trim(indx_name)//'m_cz_'//trim(of_name), action='write')
              i = pnindx(iplane)
              do ibin=1,nbins_cond
                 write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
                 do j=jmid,jmax_
                    if (j.eq.jmax_) then
                       write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*UiUjm_cz(i-1:i,j,ibin,ii,jj))
                    else
                       write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*UiUjm_cz(i-1:i,j,ibin,ii,jj))
                    end if
                 end do
              end do
              close(iunit)
           end do
        end do
     end do

     do jj=1,3
        do ii=1,3
           write(indx_name,'(2I1)'), ii,jj
           do iplane=pnmin_,pnmax_
              call get_name(iplane,of_name)
              open(unit=iunit, file=trim(output_name)//'_Rstress_cz_'//trim(indx_name)//'_'//trim(of_name), action='write')
              i = pnindx(iplane)
              do ibin=1,nbins_cond
                 write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
                 do j=jmid,jmax_
                    if (j.eq.jmax_) then
                       write (iunit,'(ES22.13)'), &
                            sum(interp_pnxm(iplane,:)*(UiUjm_cz(i-1:i,j,ibin,ii,jj) - Um_cz(i-1:i,j,ibin,ii)*Um_cz(i-1:i,j,ibin,jj)))
                    else
                       write (iunit,'(ES22.13)',advance='no'), &
                            sum(interp_pnxm(iplane,:)*(UiUjm_cz(i-1:i,j,ibin,ii,jj) - Um_cz(i-1:i,j,ibin,ii)*Um_cz(i-1:i,j,ibin,jj)))
                    end if
                 end do
              end do
              close(iunit)
           end do
        end do
     end do

     do isc=1,nscalar
        do jj=1,3
           write(indx_name,'(I1)'), jj
           do iplane=pnmin_,pnmax_
              call get_name(iplane,of_name)
              open(unit=iunit, file=trim(output_name)//'_UiSCm_cz_'//trim(indx_name)//'_'//trim(SC_name(isc))//'_'//trim(of_name), action='write')
              i = pnindx(iplane)
              do ibin=1,nbins_cond
                 write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
                 do j=jmid,jmax_
                    if (j.eq.jmax_) then
                       write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*UiSCm_cz(i-1:i,j,ibin,jj,isc))
                    else
                       write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*UiSCm_cz(i-1:i,j,ibin,jj,isc))
                    end if
                 end do
              end do
              close(iunit)
           end do
        end do
     end do

     do isc=1,nscalar
        do jj=1,3
           write(indx_name,'(I1)'), jj
           do iplane=pnmin_,pnmax_
              call get_name(iplane,of_name)
              open(unit=iunit, file=trim(output_name)//'_SCflux_cz_'//trim(indx_name)//'_'//trim(SC_name(isc))//'_'//trim(of_name), action='write')
              i = pnindx(iplane)
              do ibin=1,nbins_cond
                 write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
                 do j=jmid,jmax_
                    if (j.eq.jmax_) then
                       write (iunit,'(ES22.13)'), &
                            sum(interp_pnxm(iplane,:)*(UiSCm_cz(i-1:i,j,ibin,jj,isc) - Um_cz(i-1:i,j,ibin,jj)*SCm_cz(i-1:i,j,ibin,isc)))
                    else
                       write (iunit,'(ES22.13)',advance='no'), &
                            sum(interp_pnxm(iplane,:)*(UiSCm_cz(i-1:i,j,ibin,jj,isc) - Um_cz(i-1:i,j,ibin,jj)*SCm_cz(i-1:i,j,ibin,isc)))
                    end if
                 end do
              end do
              close(iunit)
           end do
        end do
     end do



     ! Temporary arrays for conditional velocity budgets
     do jj=1,3
        do ii=1,3
           write(indx_name,'(2I1)'), ii,jj
           do iplane=pnmin_,pnmax_
              call get_name(iplane,of_name)
              open(unit=iunit, file=trim(output_name)//'_tmp_U_cz_1_'//trim(indx_name)//'_'//trim(of_name), action='write')
              i = pnindx(iplane)
              do ibin=1,nbins_cond
                 write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
                 do j=jmid,jmax_
                    if (j.eq.jmax_) then
                       write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*tmp_U_cz_1(i-1:i,j,ibin,ii,jj))
                    else
                       write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*tmp_U_cz_1(i-1:i,j,ibin,ii,jj))
                    end if
                 end do
              end do
              close(iunit)
           end do
        end do
     end do

     do jj=1,3
        write(indx_name,'(I1)'), jj
        do iplane=pnmin_,pnmax_
           call get_name(iplane,of_name)
           open(unit=iunit, file=trim(output_name)//'_tmp_U_cz_2_'//trim(indx_name)//'_'//trim(of_name), action='write')
           i = pnindx(iplane)
           do ibin=1,nbins_cond
              write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
              do j=jmid,jmax_
                 if (j.eq.jmax_) then
                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*tmp_U_cz_2(i-1:i,j,ibin,jj))
                 else
                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*tmp_U_cz_2(i-1:i,j,ibin,jj))
                 end if
              end do
           end do
           close(iunit)
        end do
     end do

     do jj=1,3
        write(indx_name,'(I1)'), jj
        do iplane=pnmin_,pnmax_
           call get_name(iplane,of_name)
           open(unit=iunit, file=trim(output_name)//'_tmp_U_cz_3_'//trim(indx_name)//'_'//trim(of_name), action='write')
           i = pnindx(iplane)
           do ibin=1,nbins_cond
              write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
              do j=jmid,jmax_
                 if (j.eq.jmax_) then
                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*tmp_U_cz_3(i-1:i,j,ibin,jj))
                 else
                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*tmp_U_cz_3(i-1:i,j,ibin,jj))
                 end if
              end do
           end do
           close(iunit)
        end do
     end do

     do jj=1,3
        write(indx_name,'(I1)'), jj
        do iplane=pnmin_,pnmax_
           call get_name(iplane,of_name)
           open(unit=iunit, file=trim(output_name)//'_dPdxm_cz_'//trim(indx_name)//'_'//trim(of_name), action='write')
           i = pnindx(iplane)
           do ibin=1,nbins_cond
              write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
              do j=jmid,jmax_
                 if (j.eq.jmax_) then
                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*dPdxm_cz(i-1:i,j,ibin,jj))
                 else
                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dPdxm_cz(i-1:i,j,ibin,jj))
                 end if
              end do
           end do
           close(iunit)
        end do
     end do

     do jj=1,3
        do ii=1,3
           write(indx_name,'(2I1)'), ii,jj
           do iplane=pnmin_,pnmax_
              call get_name(iplane,of_name)
              open(unit=iunit, file=trim(output_name)//'_dTAUdxm_cz_'//trim(indx_name)//'_'//trim(of_name), action='write')
              i = pnindx(iplane)
              do ibin=1,nbins_cond
                 write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
                 do j=jmid,jmax_
                    if (j.eq.jmax_) then
                       write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*dTAUdxm_cz(i-1:i,j,ibin,ii,jj))
                    else
                       write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dTAUdxm_cz(i-1:i,j,ibin,ii,jj))
                    end if
                 end do
              end do
              close(iunit)
           end do
        end do
     end do



     do iplane=pnmin_,pnmax_
        call get_name(iplane,of_name)
        open(unit=iunit, file=trim(output_name)//'_w_dot_cz_'//trim(of_name), action='write')
        i = pnindx(iplane)
        do ibin=1,nbins_cond
           write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
           do j=jmid,jmax_
              if (j.eq.jmax_) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*w_dot_cz(i-1:i,j,ibin))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*w_dot_cz(i-1:i,j,ibin))
              end if
           end do
        end do
        close(iunit)
     end do


     do iplane=pnmin_,pnmax_
        call get_name(iplane,of_name)
        open(unit=iunit, file=trim(output_name)//'_chi_cz_'//trim(of_name), action='write')
        i = pnindx(iplane)
        do ibin=1,nbins_cond
           write (iunit,'(ES22.13)',advance='no'), bins_cond(ibin) ! conditional variable
           do j=jmid,jmax_
              if (j.eq.jmax_) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*chi_cz(i-1:i,j,ibin))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*chi_cz(i-1:i,j,ibin))
              end if
           end do
        end do
        close(iunit)
     end do

  end if


!!$        if (iplane.eq.pnmax) then
!!$           write(iunit,'(3ES22.13)') &
!!$                sum(interp_pnxm(iplane,:)*VISC(i-1:i,j,nz/2)), &
!!$                sum(interp_pnxm(iplane,:)*DIFF(i-1:i,j,nz/2,isc_H2O)), &
!!$                sum(interp_pnxm(iplane,:)*src_SC(i-1:i,j,nz/2,isc_H2O))
!!$        else
!!$           write(iunit,'(3ES22.13)',advance='no') &
!!$                sum(interp_pnxm(iplane,:)*VISC(i-1:i,j,nz/2)), &
!!$                sum(interp_pnxm(iplane,:)*DIFF(i-1:i,j,nz/2,isc_H2O)), &
!!$                sum(interp_pnxm(iplane,:)*src_SC(i-1:i,j,nz/2,isc_H2O))
!!$        end if
!!$        write(iunit,'(ES22.13)',advance='no') sum(interp_pnxm(iplane,:)*RHO(i-1:i,j,nz/2))
!!$        do ii=1,3
!!$           write(iunit,'(ES22.13)',advance='no') sum(interp_pnxm(iplane,:)*U(i-1:i,j,nz/2,ii))
!!$        end do
!!$        do ii=1,3
!!$           do jj=1,3
!!$              write(iunit,'(ES22.13)',advance='no') sum(interp_pnxm(iplane,:)*dUdx(i-1:i,j,nz/2,ii,jj))
!!$           end do
!!$        end do
!!$        do jj=1,3
!!$           write(iunit,'(ES22.13)',advance='no') sum(interp_pnxm(iplane,:)*dPdx(i-1:i,j,nz/2,jj))
!!$        end do
!!$        do jj=1,3
!!$           write(iunit,'(ES22.13)',advance='no') sum(interp_pnxm(iplane,:)*dRHOdx(i-1:i,j,nz/2,jj))
!!$        end do
!!$        do isc=1,nscalar
!!$           do jj=1,3
!!$              if (iplane.eq.pnmax .and. isc.eq.nscalar .and. jj.eq.3) then
!!$                 write(iunit,'(ES22.13)') sum(interp_pnxm(iplane,:)*dSCdx(i-1:i,j,nz/2,jj,isc))
!!$              else
!!$                 write(iunit,'(ES22.13)',advance='no') sum(interp_pnxm(iplane,:)*dSCdx(i-1:i,j,nz/2,jj,isc))
!!$              end if
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$  close(iunit)
  print *, irank, 'avg output done'

  return
end subroutine volumeStats_average_finalize
