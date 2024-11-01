module volumeStats_metric
  use volumeStats
  implicit none

  ! Stencil lengths
  ! ---------------
  integer, parameter :: generic_order=2
  integer :: st1,st2
  integer :: stp
  integer :: stv1,stv2
  integer :: stvp
  integer :: stc1,stc2
  integer :: stcp
  
  ! Interpolation and derivation stencil
  ! ------------------------------------
  real(WP), dimension(:),   pointer :: coeff_deriv
  real(WP), dimension(:),   pointer :: coeff_interp  
  real(WP), dimension(:,:), pointer :: coeff_interp_2D
  
  ! Interpolation operators
  ! -----------------------
  ! Velocity interpolation
  real(WP), dimension(:,:,:), pointer :: interp_u_xm,interp_v_ym,interp_w_zm
  real(WP), dimension(:,:,:), pointer :: interp_uvw_x,interp_uvw_y,interp_uvw_z,interp_vm_y
  ! Scalar interpolation
  real(WP), dimension(:,:,:),   pointer :: interp_sc_x,interp_sc_y,interp_sc_z
  real(WP), dimension(:,:,:,:), pointer :: interp_sc_xy,interp_sc_yz,interp_sc_xz
  ! Velocity interpolators
  real(WP), dimension(:,:,:), pointer :: interp_Ju_y, interp_Ju_z
  real(WP), dimension(:,:,:), pointer :: interp_Jv_x, interp_Jv_z
  real(WP), dimension(:,:,:), pointer :: interp_Jw_x, interp_Jw_y
  
  ! Divergence operator
  ! -------------------
  real(WP), dimension(:,:,:), pointer :: div_u,div_v,div_w
  ! Divergence of the convective fluxes
  real(WP), dimension(:,:,:), pointer :: divc_xx,divc_xy,divc_xz
  real(WP), dimension(:,:,:), pointer :: divc_yx,divc_yy,divc_yz
  real(WP), dimension(:,:,:), pointer :: divc_zx,divc_zy,divc_zz
  ! Length of interpolation used in div
  integer, dimension(:,:,:), pointer :: interp_xx,interp_xy,interp_xz
  integer, dimension(:,:,:), pointer :: interp_yx,interp_yy,interp_yz
  integer, dimension(:,:,:), pointer :: interp_zx,interp_zy,interp_zz
  
  ! Gradient operators
  ! ------------------
  real(WP), dimension(:,:,:), pointer :: grad_x,grad_y,grad_z
  real(WP), dimension(:,:,:), pointer :: grad_xm,grad_ym,grad_zm
  
  ! Divergence operators
  ! --------------------
  ! Trace of velocity gradient tensor
  real(WP), dimension(:,:,:), pointer :: divv_u,divv_v,divv_w
  ! Divergence of the viscous fluxes
  real(WP), dimension(:,:,:), pointer :: divv_xx,divv_xy,divv_xz
  real(WP), dimension(:,:,:), pointer :: divv_yx,divv_yy,divv_yz
  real(WP), dimension(:,:,:), pointer :: divv_zx,divv_zy,divv_zz
  
  ! Gradient operators
  ! ------------------
  ! Velocity gradients
  real(WP), dimension(:,:,:), pointer :: grad_u_x,grad_u_y,grad_u_z
  real(WP), dimension(:,:,:), pointer :: grad_v_x,grad_v_y,grad_v_z
  real(WP), dimension(:,:,:), pointer :: grad_w_x,grad_w_y,grad_w_z

  ! Grid information
  real(WP), dimension(:), pointer :: dx, dxm
  real(WP), dimension(:), pointer :: dy, dym
  real(WP), dimension(:), pointer :: dxi, dyi
  real(WP), dimension(:), pointer :: dxmi, dymi
  real(WP) :: dz, dzi  

end module volumeStats_metric


! ========================================================== !
! Initialize the grid                                        !
! ========================================================== !
subroutine volumeStats_metric_grid
  use volumeStats_metric
  use parser
  use math
  implicit none

  integer :: i,j,k
  real(WP) :: dxm_loc, dym_loc

  ! Create the midpoints
  allocate(xm(imino:imaxo),ym(jmino:jmaxo),zm(kmino:kmaxo))
  do i=1,nx
     xm(i) = 0.5_WP*(x(i)+x(i+1))
  end do
!!$  dxm_loc = xm(2) - xm(1)
!!$  do i=1,imino
!!$     xm(i) = xm(2) - real(i,WP)*dxm_loc
!!$  end do
  xm(0) = 2.0_WP*xm(1)-xm(2)
  xm(nx+1) = 2.0_WP*xm(nx)-xm(nx-1)
  do j=1,ny
     ym(j) = 0.5_WP*(y(j)+y(j+1))
  end do
  ym(0) = 2.0_WP*ym(1)-ym(2)
  ym(ny+1) = 2.0_WP*ym(ny)-ym(ny-1)
  do k=1,nz
     zm(k) = 0.5_WP*(z(k)+z(k+1))
  end do
  zm(0) = 2.0_WP*zm(1)-zm(2)
  zm(nz+1) = 2.0_WP*zm(nz)-zm(nz-1)

  allocate(dx(imino:imaxo))
  allocate(dy(jmino:jmaxo))

  allocate(dxi(imino:imaxo))
  allocate(dyi(jmino:jmaxo))

  allocate(dxm(imino:imaxo-1))
  allocate(dym(jmino:jmaxo-1))

  allocate(dxmi(imino:imaxo-1))
  allocate(dymi(jmino:jmaxo-1))

  ! Compute short hand notations - x
  do i=imin-1,imax
     dx(i)   = x(i+1) - x(i)
     dxi(i)  = 1.0_WP/dx(i)
  end do
  do i=imin-1,imax
     dxm(i)  = xm(i+1) - xm(i)
     dxmi(i) = 1.0_WP/dxm(i)
  end do
  
  ! Compute short hand notations - y
  do j=jmin-1,jmax
     dy(j)   = y(j+1) - y(j)
     dyi(j)  = 1.0_WP/dy(j)
  end do
  do j=jmin-1,jmax
     dym(j)  = ym(j+1) - ym(j)
     dymi(j) = 1.0_WP/dym(j)
  end do

  ! Compute short hand notations - z
  dz = z(kmin+1) - z(kmin)
  dzi = 1.0_WP/dz

  ! Find the midpoint in the y-direction
  j = jmin
  do while (ym(j).lt.0.0_WP.and.j.lt.jmax)
     j = j+1
  end do
  jmid = j
  if (irank.eq.iroot) print *, 'jmid=',jmid, y(jmid), ym(jmid)

  return
end subroutine volumeStats_metric_grid


subroutine volumeStats_metric_grid_reset
  use volumeStats_metric
  implicit none

  deallocate(xm,ym,zm)
  deallocate(dx,dy)
  deallocate(dxi,dyi)
  deallocate(dxm,dym)
  deallocate(dxmi,dymi)

  call volumeStats_metric_grid

end subroutine volumeStats_metric_grid_reset

! ========================================================== !
! Initialize the module                                      !
! ========================================================== !
subroutine volumeStats_metric_init
  use volumeStats_metric
  use parser
  use math
  implicit none

  integer :: i,j,k,n,st,s1,s2
  real(WP), dimension(:,:), pointer :: A,B

  ! Set the stencil lengths
  st2 = generic_order/2
  st1 = st2 - 1
  stp = st1 + st2
  stv2 = generic_order/2
  stv1 = stv2 - 1
  stvp = stv1 + stv2
  stc2 = generic_order/2
  stc1 = stc2 - 1
  stcp = stc1 + stc2

  ! Compute the interpolation coefficients
  allocate(coeff_deriv(generic_order))
  allocate(coeff_interp(generic_order))
  allocate(coeff_interp_2D(generic_order,generic_order))
  coeff_interp = 0.5_WP
  coeff_interp_2D = 0.25_WP
  coeff_deriv(1) = -1.0_WP
  coeff_deriv(2) =  1.0_WP

  ! COMPUTE PRIMARY METRICS
  ! Allocate necessary arrays
  allocate(interp_u_xm(imin_-st2:imax_+st2,jmin_-st1-1:jmax_+st2,-st1:st2))
  allocate(interp_v_ym(imin_-st2:imax_+st2,jmin_-st1-1:jmax_+st2,-st1:st2))
  allocate(interp_w_zm(imin_-st2:imax_+st2,jmin_-st1-1:jmax_+st2,-st1:st2))

  ! Larger for strain rate calculation
  allocate(interp_uvw_x(imin_-st1:imax_+st2+1,jmin_-st2:jmax_+st2,-st2:st1))
  allocate(interp_uvw_y(imin_-st1:imax_+st2+1,jmin_-st2:jmax_+st2,-st2:st1))
  allocate(interp_uvw_z(imin_-st1:imax_+st2+1,jmin_-st2:jmax_+st2,-st2:st1))
  allocate(interp_vm_y (imin_-st1:imax_+st2+1,jmin_-st2:jmax_+st2,-st2:st1))
  
  allocate(interp_sc_x(imino_-st2:imaxo_+st2+1,jmino_-st2   :jmaxo_+st2   ,-st2:st1))
  allocate(interp_sc_y(imino_   :imaxo_   +1,jmin_-st2:jmax_+st2,-st2:st1))
  allocate(interp_sc_z(imino_   :imaxo_   +1,jmin_-st2:jmax_+st2,-st2:st1))

  ! Larger because of momentum fluxes
  allocate(interp_sc_xy(imino_+st2:imaxo_-st1,jmino_+st2:jmaxo_-st1,-st2:st1,-st2:st1))
  allocate(interp_sc_yz(imino_+st2:imaxo_-st1,jmino_+st2:jmaxo_-st1,-st2:st1,-st2:st1))
  allocate(interp_sc_xz(imino_+st2:imaxo_-st1,jmino_+st2:jmaxo_-st1,-st2:st1,-st2:st1))
  
  ! Larger because of momentum fluxes
  allocate(interp_Ju_y (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
  allocate(interp_Ju_z (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))

  allocate(interp_Jv_x (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
  allocate(interp_Jv_z (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))

  allocate(interp_Jw_x (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))
  allocate(interp_Jw_y (imin_-stc1:imax_+stc2,jmin_-stc1:jmax_+stc2,-stc2:stc1))

  ! Inside only
  allocate(div_u(imin_:imax_+1,jmin_-st1:jmax_+st2,-st1:st2))
  allocate(div_v(imin_:imax_+1,jmin_-st1:jmax_+st2,-st1:st2))
  allocate(div_w(imin_:imax_+1,jmin_-st1:jmax_+st2,-st1:st2))

  ! Only inside the domain
  allocate(divc_xx(imin_:imax_,jmin_:jmax_,-stc2:stc1))
  allocate(divc_xy(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(divc_xz(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(divc_yx(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(divc_yy(imin_:imax_,jmin_:jmax_,-stc2:stc1))
  allocate(divc_yz(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(divc_zx(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(divc_zy(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(divc_zz(imin_:imax_,jmin_:jmax_,-stc2:stc1))

  ! Only inside the domain
  allocate(interp_xx(imin_:imax_,jmin_:jmax_,-stc2:stc1))
  allocate(interp_xy(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(interp_xz(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(interp_yx(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(interp_yy(imin_:imax_,jmin_:jmax_,-stc2:stc1))
  allocate(interp_yz(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(interp_zx(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(interp_zy(imin_:imax_,jmin_:jmax_,-stc1:stc2))
  allocate(interp_zz(imin_:imax_,jmin_:jmax_,-stc2:stc1))

  ! Larger for diffusion term in scalar equation
  allocate(grad_x(imin_-st1:imax_+st2+1,jmin_-st1:jmax_+st2,-st2:st1))
  allocate(grad_y(imin_-st1:imax_+st2+1,jmin_-st1:jmax_+st2,-st2:st1))
  allocate(grad_z(imin_-st1:imax_+st2+1,jmin_-st1:jmax_+st2,-st2:st1))

  ! Only inside the domain
  allocate(grad_xm(imin_:imax_,jmin_:jmax_,-stp:stp))
  allocate(grad_ym(imin_:imax_,jmin_:jmax_,-stp:stp))
  allocate(grad_zm(imin_:imax_,jmin_:jmax_,-stp:stp))


  ! Interpolation of a scalar at the faces
  allocate(A(-st2:st1,generic_order))
  allocate(B(generic_order,-st2:st1))
  do i=imin_-st1,imax_+st2
     do st=-st2,st1
        do n=1,generic_order
           A(st,n) = (xm(i+st)-x(i))**(n-1)
        end do
     end do
     call inverse_matrix(A,B,generic_order)
     do j=jmino_,jmaxo_
        interp_sc_x(i,j,:) = B(1,:)
     end do
  end do
  do j=jmin_-st1,jmax_+st2
     do st=-st2,st1
        do n=1,generic_order
           A(st,n) = (ym(j+st)-y(j))**(n-1)
        end do
     end do
     call inverse_matrix(A,B,generic_order)
     do i=imino_,imaxo_
        interp_sc_y(i,j,:) = B(1,:)
     end do
  end do
  do j=jmin_-st1,jmax_+st2
     do i=imino_,imaxo_
        interp_sc_z(i,j,:) = coeff_interp
     end do
  end do

  ! Interpolation of the scalars at the corners
  do j=jmino_+st2,jmaxo_-st1
     do i=imino_+st2,imaxo_-st1
        interp_sc_xy(i,j,:,:) = coeff_interp_2D ! use interp_sc_x to make this better?
        interp_sc_yz(i,j,:,:) = coeff_interp_2D
        interp_sc_xz(i,j,:,:) = coeff_interp_2D
     end do
  end do

  ! Interpolation of the velocities at the center
  do j=jmin_-st1-1,jmax_+st2
     do i=imin_-st2,imax_+st2
        interp_u_xm(i,j,:) = coeff_interp
        interp_v_ym(i,j,:) = coeff_interp
        interp_w_zm(i,j,:) = coeff_interp
     end do
  end do

  ! Interpolation of the velocities at the corners
  do j=jmin_-stc1,jmax_+stc2
     do i=imin_-stc1,imax_+stc2
        interp_Ju_y(i,j,:) = coeff_interp * dy(j-stc2:j+stc1) * dymi(j-1)
        interp_Ju_z(i,j,:) = coeff_interp

        interp_Jv_x(i,j,:) = coeff_interp * dx(i-stc2:i+stc1) * dxmi(i-1)
        interp_Jv_z(i,j,:) = coeff_interp

        interp_Jw_x(i,j,:) = coeff_interp * dx(i-stc2:i+stc1) * dxmi(i-1)
        interp_Jw_y(i,j,:) = coeff_interp * dy(j-stc2:j+stc1) * dymi(j-1)
     end do
  end do

  ! Interpolate centered velocities at the faces
  ! THIS NEEDS TO BE BETTER
  do j=jmin_-st2,jmax_+st2
     do i=imin_-st1,imax_+st2+1
        interp_uvw_x(i,j,:)    = interp_sc_x(i,j,:)
        interp_uvw_y(i,j,:)    = interp_sc_y(i,j,:)
        interp_vm_y (i,j,:)    = interp_sc_y(i,j,:)
        interp_uvw_z(i,j,:)    = interp_sc_z(i,j,:)
     end do
  end do
  
  ! Divergence of a vector
  do j=jmin_-st1,jmax_+st2
     do i=imin_,imax_
        div_u(i,j,:) = coeff_deriv * dxi(i)
        div_v(i,j,:) = coeff_deriv * dyi(j)
        div_w(i,j,:) = coeff_deriv * dzi
     end do
  end do

  ! Divergence of a matrix
  do j=jmin_,jmax_
     do i=imin_,imax_
        divc_xx(i,j,:) = coeff_deriv * dxmi(i-1)
        divc_xy(i,j,:) = coeff_deriv * dyi(j)
        divc_xz(i,j,:) = coeff_deriv * dzi

        divc_yx(i,j,:) = coeff_deriv * dxi(i)
        divc_yy(i,j,:) = coeff_deriv * dymi(j-1)
        divc_yz(i,j,:) = coeff_deriv * dzi

        divc_zx(i,j,:) = coeff_deriv * dxi(i)
        divc_zy(i,j,:) = coeff_deriv * dyi(j)
        divc_zz(i,j,:) = coeff_deriv * dzi
     end do
  end do

  ! Length of interpolation for divergence of a matrix
  do j=jmin_,jmax_
     do i=imin_,imax_
        do st=-stc2,stc1
           interp_xx(i,j,st) = st
           interp_yy(i,j,st) = st
           interp_zz(i,j,st) = st
        end do
        do st=-stc1,stc2
           interp_xy(i,j,st) = st
           interp_xz(i,j,st) = st
           interp_yx(i,j,st) = st
           interp_yz(i,j,st) = st           
           interp_zx(i,j,st) = st
           interp_zy(i,j,st) = st
        end do
     end do
  end do

  ! Gradient of a scalar
  do j=jmin_-st1,jmax_+st2
     do i=imin_-st1,imax_+st2+1
        grad_x(i,j,:) = coeff_deriv * dxmi(i-1)
        grad_y(i,j,:) = coeff_deriv * dymi(j-1)
        grad_z(i,j,:) = coeff_deriv * dzi
     end do
  end do
  
  ! Centered Gradient operator
  do j=jmin_,jmax_
     do i=imin_,imax_
        ! Initialize
        grad_xm(i,j,:) = 0.0_WP
        grad_ym(i,j,:) = 0.0_WP
        grad_zm(i,j,:) = 0.0_WP
        ! Convolute
        do s1=-st1,st2
           do s2=-st2,st1
              grad_xm(i,j,s1+s2) = grad_xm(i,j,s1+s2) + &
                   interp_u_xm(i,j,s1)*grad_x(i+s1,j,s2)
              grad_ym(i,j,s1+s2) = grad_ym(i,j,s1+s2) + &
                   interp_v_ym(i,j,s1)*grad_y(i,j+s1,s2)
              grad_zm(i,j,s1+s2) = grad_zm(i,j,s1+s2) + &
                   interp_w_zm(i,j,s1)*grad_z(i,j,s2)
           end do
        end do
     end do
  end do

  ! Larger because of momentum fluxes
  allocate(divv_u(imin_-stv2:imax_+stv2,jmin_-stv2:jmax_+stv2,-stv1:stv2))
  allocate(divv_v(imin_-stv2:imax_+stv2,jmin_-stv2:jmax_+stv2,-stv1:stv2))
  allocate(divv_w(imin_-stv2:imax_+stv2,jmin_-stv2:jmax_+stv2,-stv1:stv2))

  ! Only inside the domain
  allocate(divv_xx(imin_-1:imax_+1,jmin_-1:jmax_+1,-stv2:stv1))
  allocate(divv_xy(imin_-1:imax_+1,jmin_-1:jmax_+1,-stv1:stv2))
  allocate(divv_xz(imin_-1:imax_+1,jmin_-1:jmax_+1,-stv1:stv2))
  allocate(divv_yx(imin_-1:imax_+1,jmin_-1:jmax_+1,-stv1:stv2))
  allocate(divv_yy(imin_-1:imax_+1,jmin_-1:jmax_+1,-stv2:stv1))
  allocate(divv_yz(imin_-1:imax_+1,jmin_-1:jmax_+1,-stv1:stv2))
  allocate(divv_zx(imin_-1:imax_+1,jmin_-1:jmax_+1,-stv1:stv2))
  allocate(divv_zy(imin_-1:imax_+1,jmin_-1:jmax_+1,-stv1:stv2))
  allocate(divv_zz(imin_-1:imax_+1,jmin_-1:jmax_+1,-stv2:stv1))

  ! Larger because of momentum fluxes
  allocate(grad_u_x(imin_-stv2:imax_+stv2  ,jmin_-stv2:jmax_+stv2,-stv1:stv2))
  allocate(grad_u_y(imin_-stv1:imax_+stv2+1,jmin_-stv1:jmax_+stv2,-stv2:stv1))
  allocate(grad_u_z(imin_-stv1:imax_+stv2+1,jmin_-stv1:jmax_+stv2,-stv2:stv1))

  ! Larger because of momentum fluxes
  allocate(grad_v_x(imin_-stv1:imax_+stv2+1,jmin_-stv1:jmax_+stv2,-stv2:stv1))
  allocate(grad_v_y(imin_-stv2:imax_+stv2  ,jmin_-stv2:jmax_+stv2,-stv1:stv2))
  allocate(grad_v_z(imin_-stv1:imax_+stv2+1,jmin_-stv1:jmax_+stv2,-stv2:stv1))

  ! Larger because of momentum fluxes
  allocate(grad_w_x(imin_-stv1:imax_+stv2+1,jmin_-stv1:jmax_+stv2,-stv2:stv1))
  allocate(grad_w_y(imin_-stv1:imax_+stv2+1,jmin_-stv1:jmax_+stv2,-stv2:stv1))
  allocate(grad_w_z(imin_-stv2:imax_+stv2  ,jmin_-stv2:jmax_+stv2,-stv1:stv2))

  ! Use first plane
  k = kmin_

  ! Divergence of a vector
  do j=jmin_-stv2,jmax_+stv2
     do i=imin_-stv2,imax_+stv2
        call hofdd(generic_order,x(i-stv1:i+stv2),xm(i),divv_u(i,j,:))
        call hofdd(generic_order,y(j-stv1:j+stv2),ym(j),divv_v(i,j,:))
        call hofdd(generic_order,z(k-stv1:k+stv2),zm(k),divv_w(i,j,:))
     end do
  end do

  ! Divergence of a matrix
  do j=jmin_-1,jmax_+1
     do i=imin_-1,imax_+1
        call hofdd(generic_order,xm(i-stv2:i+stv1),x (i),divv_xx(i,j,:))
        call hofdd(generic_order,y (j-stv1:j+stv2),ym(j),divv_xy(i,j,:))
        call hofdd(generic_order,z (k-stv1:k+stv2),zm(k),divv_xz(i,j,:))

        call hofdd(generic_order,x (i-stv1:i+stv2),xm(i),divv_yx(i,j,:))
        call hofdd(generic_order,ym(j-stv2:j+stv1),y (j),divv_yy(i,j,:))
        call hofdd(generic_order,z (k-stv1:k+stv2),zm(k),divv_yz(i,j,:))

        call hofdd(generic_order,x (i-stv1:i+stv2),xm(i),divv_zx(i,j,:))
        call hofdd(generic_order,y (j-stv1:j+stv2),ym(j),divv_zy(i,j,:))
        call hofdd(generic_order,zm(k-stv2:k+stv1),z (k),divv_zz(i,j,:))
     end do
  end do

  ! Gradient of a vector
  do j=jmin_-stv2,jmax_+stv2
     do i=imin_-stv2,imax_+stv2
        call hofdd(generic_order,x(i-stv1:i+stv2),xm(i),grad_u_x(i,j,:))
        call hofdd(generic_order,y(j-stv1:j+stv2),ym(j),grad_v_y(i,j,:))
        call hofdd(generic_order,z(k-stv1:k+stv2),zm(k),grad_w_z(i,j,:))
     end do
  end do
  do j=jmin_-stv1,jmax_+stv2
     do i=imin_-stv1,imax_+stv2+1
        call hofdd(generic_order,xm(i-stv2:i+stv1),x(i),grad_v_x(i,j,:))
        call hofdd(generic_order,xm(i-stv2:i+stv1),x(i),grad_w_x(i,j,:))
        call hofdd(generic_order,ym(j-stv2:j+stv1),y(j),grad_u_y(i,j,:))
        call hofdd(generic_order,ym(j-stv2:j+stv1),y(j),grad_w_y(i,j,:))
        call hofdd(generic_order,zm(k-stv2:k+stv1),z(k),grad_u_z(i,j,:))
        call hofdd(generic_order,zm(k-stv2:k+stv1),z(k),grad_v_z(i,j,:))
     end do
  end do

  return
end subroutine volumeStats_metric_init

! ========================================================== !
! Calculate the gradient of a cell-centered quantity         !
! ========================================================== !
subroutine volumeStats_metric_centered_gradient(dPfiltdx, Pfilt, skip)
  use volumeStats_metric
  implicit none

  integer :: i, j, k
  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3), intent(OUT) :: dPfiltdx
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(IN) :: Pfilt
  integer,  intent(IN) :: skip

  !$OMP PARALLEL
  !$OMP DO PRIVATE(j,i)
  DO k=kmin_,kmax_, skip
     DO j=jmin_,jmax_
        DO i=imin_,imax_
           dPfiltdx(i,j,k,1) = dxi(i)*( +sum(interp_sc_x(i+1,j,:)*Pfilt(i-st2+1:i+st1+1,j,k)) &
                                        -sum(interp_sc_x(i,j,:)  *Pfilt(i-st2  :i+st1  ,j,k)) )
           dPfiltdx(i,j,k,2) = dyi(j)*( +sum(interp_sc_y(i,j+1,:)*Pfilt(i,j-st2+1:j+st1+1,k)) &
                                        -sum(interp_sc_y(i,j,:)  *Pfilt(i,j-st2  :j+st1  ,k)) )
           dPfiltdx(i,j,k,3) = dzi   *( +sum(interp_sc_z(i,j,:)  *Pfilt(i,j,k-st2+1:k+st1+1)) &
                                        -sum(interp_sc_z(i,j,:)  *Pfilt(i,j,k-st2  :k+st1  )) )
        END DO
     END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  
end subroutine volumeStats_metric_centered_gradient


! ========================================================== !
! Initialize the module                                      !
! ========================================================== !
subroutine volumeStats_metric_gradients
  use volumeStats_metric
  implicit none

  integer :: i, j, k, isc, ii, jj, kk
  real(WP), dimension(-stv2:stv2):: FX, FY, FZ
  real(WP), dimension(1:2) :: tmpv

  ! Compute all gradients at cell centers
  ! Velocity
  !$OMP PARALLEL
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           ! dU/dx
           dUdx(i,j,k,1,1) = sum(div_u(i,j,:)*U(i-st1:i+st2,j,k,1))
           ! dU/dy
           dUdx(i,j,k,1,2) = dyi(j)*( +sum(interp_sc_y(i,j+1,:)*Ui(i,j-st2+1:j+st1+1,k,1)) &
                                      -sum(interp_sc_y(i,j,  :)*Ui(i,j-st2  :j+st1  ,k,1)) )
           ! dU/dz
           dUdx(i,j,k,1,3) = dzi   *( +sum(interp_sc_z(i,j,:)*Ui(i,j,k-st2+1:k+st1+1,1)) &
                                      -sum(interp_sc_z(i,j,:)*Ui(i,j,k-st2  :k+st1  ,1)) )
           ! dV/dx
           dUdx(i,j,k,2,1) = dxi(i)*( +sum(interp_sc_x(i+1,j,:)*Ui(i-st2+1:i+st1+1,j,k,2)) &
                                      -sum(interp_sc_x(i  ,j,:)*Ui(i-st2  :i+st1  ,j,k,2)) )
           ! dV/dy
           dUdx(i,j,k,2,2) = sum(div_v(i,j,:)*U(i,j-st1:j+st2,k,2))
           ! dV/dz
           dUdx(i,j,k,2,3) = dzi   *( +sum(interp_sc_z(i,j,:)*Ui(i,j,k-st2+1:k+st1+1,2)) &
                                      -sum(interp_sc_z(i,j,:)*Ui(i,j,k-st2  :k+st1  ,2)) )
           ! dW/dx
           dUdx(i,j,k,3,1) = dxi(i)*( +sum(interp_sc_x(i+1,j,:)*Ui(i-st2+1:i+st1+1,j,k,3)) &
                                      -sum(interp_sc_x(i  ,j,:)*Ui(i-st2  :i+st1  ,j,k,3)) )
           ! dW/dy
           dUdx(i,j,k,3,2) = dyi(j)*( +sum(interp_sc_y(i,j+1,:)*Ui(i,j-st2+1:j+st1+1,k,3)) &
                                      -sum(interp_sc_y(i,j  ,:)*Ui(i,j-st2  :j+st1  ,k,3)) )
           ! dW/dz
           dUdx(i,j,k,3,3) = sum(div_w(i,j,:)*U(i,j,k-st1:k+st2,3))
        end do
     end do
  end do
  !$OMP END DO
  ! Pressure
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           dPdx(i,j,k,1) = dxi(i)*( +sum(interp_sc_x(i+1,j,:)*P(i-st2+1:i+st1+1,j,k)) &
                                    -sum(interp_sc_x(i,j,:)  *P(i-st2:i+st1,j,k))     )
           dPdx(i,j,k,2) = dyi(j)*( +sum(interp_sc_y(i,j+1,:)*P(i,j-st2+1:j+st1+1,k)) &
                                    -sum(interp_sc_y(i,j,:)  *P(i,j-st2:j+st1,k))     )
           dPdx(i,j,k,3) = dzi   *( +sum(interp_sc_z(i,j,:)  *P(i,j,k-st2+1:k+st1+1)) &
                                    -sum(interp_sc_z(i,j,:)  *P(i,j,k-st2:k+st1))     )
        end do
     end do
  end do
  !$OMP END DO
  ! Density
  !$OMP DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           dRHOdx(i,j,k,1) = dxi(i)*( +sum(interp_sc_x(i+1,j,:)*RHO(i-st2+1:i+st1+1,j,k)) &
                                      -sum(interp_sc_x(i,j,:)  *RHO(i-st2:i+st1,j,k))     )
           dRHOdx(i,j,k,2) = dyi(j)*( +sum(interp_sc_y(i,j+1,:)*RHO(i,j-st2+1:j+st1+1,k)) &
                                      -sum(interp_sc_y(i,j,:)  *RHO(i,j-st2:j+st1,k))     )
           dRHOdx(i,j,k,3) = dzi   *( +sum(interp_sc_z(i,j,:)  *RHO(i,j,k-st2+1:k+st1+1)) &
                                      -sum(interp_sc_z(i,j,:)  *RHO(i,j,k-st2:k+st1))     )
        end do
     end do
  end do
  !$OMP END DO
  ! Scalars
  if (combust) then
     do isc=1,nscalar
        !$OMP DO
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 dSCdx(i,j,k,1,isc) = dxi(i)*( +sum(interp_sc_x(i+1,j,:)*SC(i-st2+1:i+st1+1,j,k,isc)) &
                      -sum(interp_sc_x(i  ,j,:)*SC(i-st2  :i+st1  ,j,k,isc)) )
                 dSCdx(i,j,k,2,isc) = dyi(j)*( +sum(interp_sc_y(i,j+1,:)*SC(i,j-st2+1:j+st1+1,k,isc)) &
                      -sum(interp_sc_y(i,j  ,:)*SC(i,j-st2  :j+st1  ,k,isc)) )
                 dSCdx(i,j,k,3,isc) = dzi   *( +sum(interp_sc_z(i,j,:)  *SC(i,j,k-st2+1:k+st1+1,isc)) &
                      -sum(interp_sc_z(i,j,:)  *SC(i,j,k-st2  :k+st1  ,isc)) )
              end do
           end do
        end do
        !$OMP END DO
     end do
  end if

  ! Stress tensor gradients at the cell midpoints
  !  -- Assumes velocity components U(:,:,ii) are at the cell faces
  !$OMP DO PRIVATE(j,i,ii,jj,kk,FX,FY,FZ,tmpv)
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           ! dTAUdx11
           FX = 0.0_WP
           do ii=i-stv2,i+stv2
              FX(ii-i) = &
                   + 2.0_WP*VISC(ii,j,k)*( &
                   + sum(grad_u_x(ii,j,:)*U(ii-stv1:ii+stv2,j,k,1)) &
                   - 1.0_WP/3.0_WP*( sum(divv_u(ii,j,:)*U(ii-stv1:ii+stv2,j,k,1))  &
                                   + sum(divv_v(ii,j,:)*U(ii, j-stv1:j+stv2,k,2))  &
                                   + sum(divv_w(ii,j,:)*U(ii, j,k-stv1:k+stv2,3))) )
           end do
           tmpv(1) = sum(divv_xx(i  ,j,:)*FX(-stv2:stv1))
           tmpv(2) = sum(divv_xx(i+1,j,:)*FX( stv1:stv2))
           dTAUdx(i,j,k,1,1) = sum(interp_u_xm(i,j,:)*tmpv)
           ! dTAUdx12
           FY = 0.0_WP
           do ii=i,i+1
              do jj=j-stv1,j+stv2
                 FY(jj-j) = &
                      + sum(interp_sc_xy(ii,jj,:,:)*VISC(ii-st2:ii+st1,jj-st2:jj+st1,k)) * &
                      ( sum(grad_u_y(ii,jj,:)*U(ii,jj-stv2:jj+stv1,k,1)) &
                      + sum(grad_v_x(ii,jj,:)*U(ii-stv2:ii+stv1,jj,k,2)) )
              end do
              tmpv(ii-i+1) = sum(divv_xy(ii,j,:)*FY(-stv1:stv2))
           end do
           dTAUdx(i,j,k,1,2) = sum(interp_u_xm(i,j,:)*tmpv)
           ! dTAUdx13
           FZ = 0.0_WP
           do ii=i,i+1
              do kk=k-stv1,k+stv2
                 FZ(kk-k) = &
                      + sum(interp_sc_xz(ii,j,:,:)*VISC(ii-st2:ii+st1,j,kk-st2:kk+st1)) * &
                      ( sum(grad_u_z(ii,j,:)*U(ii,j,kk-stv2:kk+stv1,1)) &
                      + sum(grad_w_x(ii,j,:)*U(ii-stv2:ii+stv1,j,kk,3)) )
              end do
              tmpv(ii-i+1) = sum(divv_xz(ii,j,:)*FZ(-stv1:stv2))
           end do
           dTAUdx(i,j,k,1,3) = sum(interp_u_xm(i,j,:)*tmpv)

           ! dTAUdx21
           FX = 0.0_WP
           do jj=j,j+1
              do ii=i-stv1,i+stv2
                 FX(ii-i) = &
                      + sum(interp_sc_xy(ii,jj,:,:)*VISC(ii-st2:ii+st1,jj-st2:jj+st1,k)) * &
                      ( sum(grad_u_y(ii,jj,:)*U(ii,jj-stv2:jj+stv1,k,1)) &
                      + sum(grad_v_x(ii,jj,:)*U(ii-stv2:ii+stv1,jj,k,2)) )
              end do
              tmpv(jj-j+1) = sum(divv_yx(i,jj,:)*FX(-stv1:stv2))
           end do
           dTAUdx(i,j,k,2,1) = sum(interp_v_ym(i,j,:)*tmpv)
           ! dTAUdx22
           FY = 0.0_WP
           do jj=j-stv2,j+stv2
              FY(jj-j) = &
                   + 2.0_WP*VISC(i,jj,k)*( &
                   + sum(grad_v_y(i,jj,:)*U(i,jj-stv1:jj+stv2,k,2)) &
                   - 1.0_WP/3.0_WP*( sum(divv_u(i,jj,:)*U(i-stv1:i+stv2,jj, k,1))  &
                                   + sum(divv_v(i,jj,:)*U(i,jj-stv1:jj+stv2,k,2))  &
                                   + sum(divv_w(i,jj,:)*U(i,jj, k-stv1:k+stv2,3))) )
           end do
           tmpv(1) = sum(divv_yy(i,j  ,:)*FY(-stv2:stv1))
           tmpv(2) = sum(divv_yy(i,j+1,:)*FY( stv1:stv2))
           dTAUdx(i,j,k,2,2) = sum(interp_v_ym(i,j,:)*tmpv)
           ! dTAUdx23
           FZ = 0.0_WP
           do jj=j,j+1
              do kk=k-stv1,k+stv2
                 FZ(kk-k) = &
                      + sum(interp_sc_yz(i,jj,:,:)*VISC(i,jj-st2:jj+st1,kk-st2:kk+st1)) * &
                      ( sum(grad_v_z(i,jj,:)*U(i,jj,kk-stv2:kk+stv1,2)) &
                      + sum(grad_w_y(i,jj,:)*U(i,jj-stv2:jj+stv1,kk,3)) )
              end do
              tmpv(jj-j+1) = sum(divv_yz(i,jj,:)*FZ(-stv1:stv2))
           end do
           dTAUdx(i,j,k,2,3) = sum(interp_v_ym(i,j,:)*tmpv)

           !dTAUdx31
           FX = 0.0_WP
           do kk=k,k+1
              do ii=i-stv1,i+stv2
                 FX(ii-i) = &
                      + sum(interp_sc_xz(ii,j,:,:)*VISC(ii-st2:ii+st1,j,kk-st2:kk+st1)) * &
                      ( sum(grad_u_z(ii,j,:)*U(ii,j,kk-stv2:kk+stv1,1)) &
                      + sum(grad_w_x(ii,j,:)*U(ii-stv2:ii+stv1,j,kk,3)) )
              end do
              tmpv(kk-k+1) = sum(divv_zx(i,j,:)*FX(-stv1:stv2))
           end do
           dTAUdx(i,j,k,3,1) = sum(interp_w_zm(i,j,:)*tmpv)
           ! dTAUdx32
           FY = 0.0_WP
           do kk=k,k+1
              do jj=j-stv1,j+stv2
                 FY(jj-j) = &
                      + sum(interp_sc_yz(i,jj,:,:)*VISC(i,jj-st2:jj+st1,kk-st2:kk+st1)) * &
                      ( sum(grad_v_z(i,jj,:)*U(i,jj,kk-stv2:kk+stv1,2)) &
                      + sum(grad_w_y(i,jj,:)*U(i,jj-stv2:jj+stv1,kk,3)) )
              end do
              tmpv(kk-k+1) = sum(divv_zy(i,j,:)*FY(-stv1:stv2))
           end do
           dTAUdx(i,j,k,3,2) = sum(interp_w_zm(i,j,:)*tmpv)
           ! dTAUdx33
           FZ = 0.0_WP
           do kk=k-stv2,k+stv2
              FZ(kk-k) = &
                   + 2.0_WP*VISC(i,j,kk)*( &
                   + sum(grad_w_z(i,j,:)*U(i,j,kk-stv1:kk+stv2,3)) &
                   - 1.0_WP/3.0_WP*( sum(divv_u(i,j,:)*U(i-stv1:i+stv2,j, kk,1))  &
                                   + sum(divv_v(i,j,:)*U(i,j-stv1:j+stv2, kk,2))  &
                                   + sum(divv_w(i,j,:)*U(i,j,kk-stv1:kk+stv2,3)) ))
           end do
           tmpv(1) = sum(divv_zz(i,j,:)*FZ(-stv2:stv1))
           tmpv(2) = sum(divv_zz(i,j,:)*FZ( stv1:stv2))
           dTAUdx(i,j,k,3,3) = sum(interp_w_zm(i,j,:)*tmpv)
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  return
end subroutine volumeStats_metric_gradients


! ========================================================== !
! Compute strain rate                                        !
! ========================================================== !
subroutine volumeStats_metric_strainrate(dUdx_in, Sij_out)
  use volumeStats_metric
  implicit none

  integer :: i, j, k, ii, jj
  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3,1:3), intent(OUT) :: Sij_out
  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3,1:3), intent(IN) :: dUdx_in

  !$OMP PARALLEL
  do jj=1,3
     do ii=jj,3
        !$OMP DO PRIVATE(j,i)
        do k=kmin_,kmax_
           do j=jmin_,jmax_
              do i=imin_,imax_
                 Sij_out(i,j,k,ii,jj) = 0.5_WP * ( dUdx_in(i,j,k,ii,jj) + dUdx_in(i,j,k,jj,ii) )
              end do
           end do
        end do
        !$OMP END DO
     end do
  end do
  !$OMP END PARALLEL

  Sij_out(:,:,:,1,2) = Sij_out(:,:,:,2,1)
  Sij_out(:,:,:,1,3) = Sij_out(:,:,:,3,1)
  Sij_out(:,:,:,2,3) = Sij_out(:,:,:,3,2)

end subroutine volumeStats_metric_strainrate
