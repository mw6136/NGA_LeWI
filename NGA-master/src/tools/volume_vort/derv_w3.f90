module derv_w3
  use parallel
  use string
  use precision
  use geometry
  use vort_mod

  implicit none

contains

subroutine dx_w3(dir)
  implicit none
  integer :: m, l, dir
  integer :: i, j, k, ll, ii, maxx
  real(WP), dimension(5) :: xst, dst
  real(WP), dimension(2) :: beta, gamma, fp, fm, w, wt
  real(WP) :: eps, fluxp, fluxm
  real(WP), dimension(:), pointer :: xx

  if (dir == 1) then
     allocate(xx(nx1))
     xx(:) = x1(:)
     maxx = nx1
  end if
  if (dir == 2) then
     allocate(xx(ny1))
     xx(:) = y1(:)
     maxx = ny1
  end if
  if (dir == 3) then
     allocate(xx(nz1))
     xx(:) = z1(:)
     maxx = nz1
  end if

  ddata(:,:,:) = 0.0_WP
  eps = 1.0e-6_WP
  gamma(1) = 1.0_WP/3.0_WP; gamma(2) = 2.0_WP/3.0_WP
  do k = 1,nz1
     do j = 1,ny1
        do i =1,nx1
           if ((i==1 .and. dir==1) .or. (j==1 .and. dir==2) .or. (k==1 .and. dir==3)) then
              xp1 = xx(1) - dabs(xx(2)-xx(1))
              p1 = 0.0_WP
              do ll = 1,4
                 num = 1.0_WP; den = 1.0_WP
                 do ii = 1,4
                    if (ii == ll) cycle
                    num = num * (xp1-xx(ii))
                    den = den * (xx(ll) - xx(ii))
                 end do
              if (dir==1) p1 = p1 +(num/den)*data8(ll, j, k)
              if (dir==2) p1 = p1 +(num/den)*data8(i, ll, k)
              if (dir==3) p1 = p1 +(num/den)*data8(i, j, ll)
              end do
              xp2 = xx(1) - 2.0_WP*dabs(xx(2)-xx(1))
              xst(1) = xp1; xst(2:4) = xx(1:3)
              dst(1) = p1
              p2 = 0.0_WP
              if (dir==1) dst(2:4) = data8(1:3, j, k)
              if (dir==2) dst(2:4) = data8(i, 1:3, k)
              if (dir==3) dst(2:4) = data8(i, j, 1:3)
              do ll = 1,4
                 num = 1.0_WP; den = 1.0_WP
                 do ii = 1,4
                    if (ii == ll) cycle
                    num = num * (xp2-xst(ii))
                    den = den * (xst(ll) - xst(ii))
                 end do
              p2 = p2 +(num/den)*dst(ll)
              end do
           end if
           if ((i==2 .and. dir==1) .or. (j==2 .and. dir==2) .or. (k==2 .and. dir==3))  then
              xp1 = xx(1) - dabs(xx(2)-xx(1))
              p1 = 0.0_WP
              do ll = 1,4
                 num = 1.0_WP; den = 1.0_WP
                 do ii = 1,4
                    if (ii == ll) cycle
                    num = num * (xp1-xx(ii))
                    den = den * (xx(ll) - xx(ii))
                 end do
              if (dir==1) p1 = p1 +(num/den)*data8(ll, j, k)
              if (dir==2) p1 = p1 +(num/den)*data8(i, ll, k)
              if (dir==3) p1 = p1 +(num/den)*data8(i, j, ll)
              end do
           end if
           if ((i==maxx .and. dir==1) .or. (j==maxx .and. dir==2) .or. (k==maxx .and. dir==3)) then
              xp1 = xx(maxx) +(xx(maxx)-xx(maxx-1))
              p1 = 0.0_WP
              do ll = 1,4
                 num = 1.0_WP; den = 1.0_WP
                 do ii = 1,4
                    if (ii == ll) cycle
                    num = num * (xp1-xx(maxx-4+ii))
                    den = den * (xx(ll) - xx(maxx-4+ii))
                 end do
              if (dir==1) p1 = p1 +(num/den)*data8(maxx-4+ll, j, k)
              if (dir==2) p1 = p1 +(num/den)*data8(i, maxx-4+ll, k)
              if (dir==3) p1 = p1 +(num/den)*data8(i, j, maxx-4+ll)
              end do
              xp2 = xx(maxx) +2.0_WP*(xx(maxx)-xx(maxx-1))
              xst(1:3) = xx(maxx-2:maxx); xst(4) = xp1
              if (dir==1) dst(1:3) = data8(maxx-2:maxx, j, k)
              if (dir==2) dst(1:3) = data8(i, maxx-2:maxx, k)
              if (dir==3) dst(1:3) = data8(i, j, maxx-2:maxx)
              dst(4) = p1
              p2 = 0.0_WP
              do ll = 1,4
                 num = 1.0_WP; den = 1.0_WP
                 do ii = 1,4
                    if (ii == ll) cycle
                    num = num * (xp2-xst(ii))
                    den = den * (xst(ll) - xst(ii))
                 end do
              p2 = p2 +(num/den)*dst(ll)
              end do
           end if
           if ((i==maxx-1 .and. dir==1) .or. (j==maxx-1 .and. dir==2) .or. (k==maxx-1 .and. dir==3)) then
              xp1 = xx(maxx) +(xx(maxx)-xx(maxx-1))
              p1 = 0.0_WP
              do ll = 1,4
                 num = 1.0_WP; den = 1.0_WP
                 do ii = 1,4
                    if (ii == ll) cycle
                    num = num * (xp1-xx(maxx-4+ii))
                    den = den * (xx(ll) - xx(maxx-4+ii))
                 end do
              if (dir==1) p1 = (num/den)*data8(maxx-4+ll, j, k)
              if (dir==2) p1 = (num/den)*data8(i, maxx-4+ll, k)
              if (dir==3) p1 = (num/den)*data8(i, j, maxx-4+ll)
              end do
           end if

           if ((i==1 .and. dir==1) .or. (j==1 .and. dir==2) .or. (k==1.and. dir==3)) then
              xst(1) = xp2; xst(2) = xp1; xst(3:5) = xx(1:3)
              dst(1) = p2; dst(2) = p1
              if (dir==1) dst(3:5) = data8(1:3, j, k)
              if (dir==2) dst(3:5) = data8(i, 1:3, k)
              if (dir==3) dst(3:5) = data8(i, j, 1:3)
           else if ((i==2 .and. dir==1) .or. (j==2 .and. dir==2) .or. (k==2 .and. dir==3)) then
              xst(1) = xp1; xst(2:5) = xx(1:4)
              dst(1) = p1
              if (dir==1) dst(2:5) = data8(1:4, j, k)
              if (dir==2) dst(2:5) = data8(i, 1:4, k)
              if (dir==3) dst(2:5) = data8(i, j, 1:4) 
           else if ((i==maxx .and. dir==1) .or. (j==maxx .and. dir==2) .or. (k==maxx .and. dir==3)) then
              xst(1:3) = xx(maxx-2:maxx); xst(4) = xp1; xst(5) = xp2
              if (dir==1) dst(1:3) = data8(maxx-2:maxx, j, k)
              if (dir==2) dst(1:3) = data8(i, maxx-2:maxx, k)
              if (dir==3) dst(1:3) = data8(i, j, maxx-2:maxx)
              dst(4) = p1; dst(5) = p2
           else if ((i==maxx-1 .and. dir==1) .or. (j==maxx-1 .and. dir==2) .or. (k==maxx-1 .and. dir==3)) then
              xst(1:4) = xx(maxx-3:maxx); xst(5) = xp1
              if (dir==1) dst(1:4) = data8(maxx-3:maxx, j, k)
              if (dir==2) dst(1:4) = data8(i, maxx-3:maxx, k)
              if (dir==3) dst(1:4) = data8(i, j, maxx-3:maxx)
              dst(5) = p1
           else
              if (dir==1) then
                 xst(1:5) = xx(i-2:i+2)
                 dst(1:5) = data8(i-2:i+2, j, k)
              end if
              if (dir==2) then
                 xst(1:5) = xx(j-2:j+2)
                 dst(1:5) = data8(i, j-2:j+2, k)
              end if
              if (dir==3) then
                 xst(1:5) = xx(k-2:k+2)
                 dst(1:5) = data8(i, j, k-2:k+2)
              end if
           end if

           if (dst(3) .ge. 0.0_WP) then
              beta(1) = (dst(3)-dst(2))**2.0_WP
              beta(2) = (dst(4)-dst(3))**2.0_WP
              wt(1) = gamma(1)/((eps+beta(1))**2.0_WP)
              wt(2) = gamma(2)/((eps+beta(2))**2.0_WP)
              w(1) = wt(1)/(wt(1)+wt(2))
              w(2) = wt(2)/(wt(1)+wt(2))
              fp(1) = -0.5_WP*dst(2) +(3.0_WP/2.0_WP)*dst(3)
              fp(2) = 0.5_WP*dst(3) +0.5_WP*dst(4)
              fluxp = w(1)*fp(1) +w(2)*fp(2)
              beta(1) = (dst(2)-dst(1))**2.0_WP
              beta(2) = (dst(3)-dst(2))**2.0_WP
              wt(1) = gamma(1)/((eps+beta(1))**2.0_WP)
              wt(2) = gamma(2)/((eps+beta(2))**2.0_WP)
              w(1) = wt(1)/(wt(1)+wt(2))
              w(2) = wt(2)/(wt(1)+wt(2))
              fm(1) = -0.5_WP*dst(1) +(3.0_WP/2.0_WP)*dst(2)
              fm(2) = 0.5_WP*dst(2) +0.5_WP*dst(3)
              fluxm =  w(1)*fm(1) +w(2)*fm(2)
           else
              beta(1) = (dst(4)-dst(5))**2.0_WP
              beta(2) = (dst(3)-dst(4))**2.0_WP
              wt(1) = gamma(1)/((eps+beta(1))**2.0_WP)
              wt(2) = gamma(2)/((eps+beta(2))**2.0_WP)
              w(1) = wt(1)/(wt(1)+wt(2))
              w(2) = wt(2)/(wt(1)+wt(2))
              fp(1) = -0.5_WP*dst(5) +(3.0_WP/2.0_WP)*dst(4)
              fp(2) = 0.5_WP*dst(4) +0.5_WP*dst(3)
              fluxp = w(1)*fp(1) +w(2)*fp(2)
              beta(1) = (dst(3)-dst(4))**2.0_WP
              beta(2) = (dst(2)-dst(3))**2.0_WP
              wt(1) = gamma(1)/((eps+beta(1))**2.0_WP)
              wt(2) = gamma(2)/((eps+beta(2))**2.0_WP)
              w(1) = wt(1)/(wt(1)+wt(2))
              w(2) = wt(2)/(wt(1)+wt(2))
              fm(1) = -0.5_WP*dst(4) +(3.0_WP/2.0_WP)*dst(3)
              fm(2) = 0.5_WP*dst(3) +0.5_WP*dst(2)
              fluxm =  w(1)*fm(1) +w(2)*fm(2)
           endif

           ddata(i,j,k) = (fluxp-fluxm)/((xst(4)-xst(2))/2.0_WP)

        end do
     end do
  end do


end subroutine

end module derv_w3

