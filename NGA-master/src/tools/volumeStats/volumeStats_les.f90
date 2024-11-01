! ============================================================ !
!                    volumeStats_les.f90                       !
!                                                              !
! Authors: Omkar B. Shende and Jonathan F. MacArt              !
! Date:   May -- July, 2018                                    !
! ============================================================ !
module volumeStats_les
  use volumeStats

  ! Filter values
  real(WP), dimension(:,:), allocatable :: xfilt, yfilt, zfilt  !NOT POINTERS
  integer,  dimension(:,:), allocatable :: xfbounds, yfbounds, zfbounds

  ! Filtered quantities
  real(WP), dimension(:,:,:,:),   pointer :: Uifilt
  real(WP), dimension(:,:,:,:,:), pointer :: UiUjfilt
  real(WP), dimension(:,:,:,:,:), pointer :: dUifiltdx

  real(WP), dimension(:,:,:),     pointer :: rhofilt
  real(WP), dimension(:,:,:),     pointer :: VISCfilt

  real(WP), dimension(:,:,:,:),   pointer :: SCfilt
  real(WP), dimension(:,:,:,:,:), pointer :: UiSCfilt
  real(WP), dimension(:,:,:,:,:), pointer :: dSCfiltdx
  !Progress Variable
  real(WP), dimension(:,:,:),     pointer :: PROGfilt
  real(WP), dimension(:,:,:,:),   pointer :: dPROGfiltdx

  real(WP), dimension(:,:,:),     pointer :: epsTA
  real(WP), dimension(:,:,:),     pointer :: eps
  real(WP), dimension(:,:,:),     pointer :: Dafilt


  !filtered density variable as well rhofilt

!!$  real(WP), dimension(:,:,:,:,:), pointer :: SFF        !Subfilter scalar flux
!!$  real(WP), dimension(:,:,:,:,:), pointer :: SFS        !Subfilter stress

  real(WP), dimension(:,:,:,:,:), pointer :: Sij
  real(WP), dimension(:,:,:,:,:), pointer :: Sijfilt
!!$  real(WP), dimension(:,:,:),     pointer :: ETAkol

!!$  real(WP), dimension(:,:,:,:,:), pointer :: SFFmodel

  real(WP), dimension(:,:), pointer :: SCval

  real(WP) :: Smag, Sc_t, sL, dF, ca, tau 
  integer  :: unburned, burned, k_skip, isc_CASE
       
  ! Data type for joint pdfs
  type PDF_VAR
     character(len=str_short) :: name
     integer*8, dimension(:,:,:), pointer :: pdf
!     real(WP),  dimension(:,:,:), pointer :: data
     real(WP) :: min_x, max_x
     real(WP) :: min_y, max_y
     logical :: init, j_flip
  end type PDF_VAR
  
  ! Joint pdf variables
  integer :: nvar_pdf
  type(PDF_VAR), dimension(:), pointer :: CPDF

contains
  ! ========================================================== !
  ! Filter kernel for LES                                      !
  ! ========================================================== !
  subroutine filter_kernel_offs(input,output,offset,d_weight)
    implicit none

    integer, intent(in) :: offset
    real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(in)    :: input
    real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(inout) :: output
    integer  :: a, b, c, i, j, k
    logical :: d_weight

    !$OMP PARALLEL
    do k = kmin_+offset, kmax_+offset, k_skip
       !$OMP DO PRIVATE(j, i, c, b, a)
       do j = jmin_, jmax_
          do i = imin_, imax_
             do c = zfbounds(k,1), zfbounds(k,2)
                do b = yfbounds(j,1), yfbounds(j,2)
                   do a = xfbounds(i,1), xfbounds(i,2)
                      output(i,j,k) = output(i,j,k) + xfilt(i,a)*yfilt(j,b)*zfilt(k,c)*input(a,b,c)
                   end do
                end do
             end do
          end do
       end do
       !$OMP END DO
    end do
    !$OMP END PARALLEL

    if (d_weight) then
       ! Complete the density weighting
       !$OMP PARALLEL
       do k = kmin_+offset, kmax_+offset, k_skip
          !$OMP DO PRIVATE(j, i)
          do j = jmin_, jmax_
             do i = imin_, imax_
                output(i,j,k) = output(i,j,k)/RHOfilt(i,j,k)
             end do
          end do
          !$OMP END DO
       end do
       !$OMP END PARALLEL
    end if
    
  end subroutine filter_kernel_offs

  ! ========================================================== !
  ! PDF kernel                                                 !
  ! ========================================================== !
  subroutine cond_pdf_kernel(PDF,data)
    implicit none

    real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_), intent(in) :: data
    type(PDF_VAR), intent(inout) :: PDF
    real(WP) :: tmp1
    integer :: i, j, k, ix, iy

!!$    ! Density weighting and y-flip if necessary
!!$    if (PDF%j_flip) then
!!$       do k=kmin_,kmax_,k_skip
!!$          do j=jmin_,jmid-1
!!$             do i=imin_,imax_
!!$                data(i,j,k) = -1.0_WP*data(i,j,k)
!!$             end do
!!$          end do
!!$       end do
!!$    end if

    ! Get the min/max if not initialized
    if (.not.PDF%init) then
       tmp1 = maxval(data)
       PDF%max_y = tmp1 + 0.05_WP*abs(tmp1)
       tmp1 = minval(data)
       PDF%min_y = tmp1 - 0.05_WP*abs(tmp1)
       PDF%init = .true.
    end if
    
    ! Update the PDF
    do k=kmin_+k_skip,kmax_-k_skip,k_skip
       do j=jmin_,jmax_
          do i=imin_,imax_
             ix   = bin_index(i,j,k) ! Set from PROGfilt
             tmp1 = (PDF%max_y - data(i,j,k))/(PDF%max_y - PDF%min_y)
             iy   = min(nbins_cond, max(1,floor(tmp1*nbins_cond)+1))
             PDF%pdf(i,ix,iy) = PDF%pdf(i,ix,iy) + 1
          end do
       end do
    end do
    
  end subroutine cond_pdf_kernel
  
end module volumeStats_les


! ========================================================== !
! Initialize the module                                      !
! ========================================================== !
subroutine volumeStats_les_init
  use volumeStats_les
  use parser
  implicit none

  integer :: ivar, isc, ii, jj
  real(WP), dimension(2) :: pdf_lim
  character(len=str_short) :: tmpstr, tmpstr2
  
  ! Filtered scalars
  allocate( SCfilt  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nscalar))
  allocate(dSCfiltdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3,1:nscalar))

  ! Filtered velocities
  allocate(Uifilt   (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:3))
  allocate(dUifiltdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3,1:3))
  allocate(UiUjfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_, 1:3, 1:3))
  allocate(UiSCfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_, 1:3, 1:nscalar))

  ! Filtered scalar quantities
  allocate(rhofilt(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(VISCfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(PROGfilt(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))

  ! Dissipation rates
  allocate(epsTA(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(eps(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(Dafilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_))

  ! Filtered gradients
  allocate(dPROGfiltdx(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:3))
  allocate(Sij(imin_:imax_,jmin_:jmax_,kmin_:kmax_, 1:3, 1:3))
  allocate(Sijfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_, 1:3, 1:3))

  
  ! Kolmogorov scale
!!$  allocate(ETAkol(imin_:imax_,jmin_:jmax_,kmin_:kmax_))

  ! Shouldn't need these parameters -- move to matlab
  
  Smag = 0.17*0.17                  !Smagorinsky constant
  Sc_t = 0.7                        !Turbulent Schmidt number
  sL = 1.195                        !Laminar flame speed
  dF = 0.000435                     !Flame thickness
  ca = 1/1.0                        !c_alpha from CTR proposal
  tflaminar = sL/dF                 !1/flame time scale
  tau = 6.7                         !Heat release factor, Kacr???
  burned = 1
  unburned = 2

  ! Skip in Z direction
  call parser_read('LES skip k',k_skip,1)
  isc_CASE = isc_H2O 

  ! Scalar burned/unburned values
  !   For 9-species hydrogen mechanism
  allocate(SCval(1:nscalar,1:2))
  SCval(isc_ZMIX,   burned) = 0.00000E+00
  SCval(isc_ZMIX, unburned) = 1.00000E+00
  SCval(isc_H,      burned) = 8.90741E-06
  SCval(isc_H,    unburned) = 0.00000E+00
  SCval(isc_O,      burned) = 5.89262E-05
  SCval(isc_O,    unburned) = 0.00000E+00
  SCval(isc_O2,     burned) = 4.11205E-03
  SCval(isc_O2,   unburned) = 1.69544E-01
  SCval(isc_OH,     burned) = 1.43446E-03
  SCval(isc_OH,   unburned) = 0.00000E+00
  SCval(isc_H2,     burned) = 2.66284E-04
  SCval(isc_H2,   unburned) = 2.13626E-02
  SCval(isc_HO2,    burned) = 3.16753E-07
  SCval(isc_HO2,  unburned) = 0.00000E+00
  SCval(isc_H2O,    burned) = 1.84539E-01
  SCval(isc_H2O,  unburned) = 0.00000E+00
  SCval(isc_H2O2,   burned) = 7.00111E-08
  SCval(isc_H2O2, unburned) = 0.00000E+00
  SCval(isc_T,      burned) = 2.04749E+03
  SCval(isc_T,    unburned) = 3.00000E+02


  ! Joint pdfs
  nvar_pdf = 1 + 3 + 8*nscalar + 3 + 2*6
  ! Allocate the pdf type
  allocate(CPDF(nvar_pdf))
  CPDF(:)%init = .false.
  ivar = 1

  ! Density
  CPDF(ivar)%name = 'RHO'
  call parser_read('PDF limits RHO',pdf_lim)
  CPDF(ivar)%min_y = pdf_lim(1); CPDF(ivar)%max_y = pdf_lim(2)
  CPDF(ivar)%init = .true.
  CPDF(ivar)%j_flip = .false.
  ivar = ivar + 1

  ! Velocity
  CPDF(ivar)%name = 'U'
  call parser_read('PDF limits U',pdf_lim)
  CPDF(ivar)%min_y = pdf_lim(1); CPDF(ivar)%max_y = pdf_lim(2)
  CPDF(ivar)%init = .true.
  CPDF(ivar)%j_flip = .false.
  ivar = ivar + 1
  CPDF(ivar)%name = 'V'
  call parser_read('PDF limits V',pdf_lim)
  CPDF(ivar)%min_y = pdf_lim(1); CPDF(ivar)%max_y = pdf_lim(2)
  CPDF(ivar)%init = .true.
  CPDF(ivar)%j_flip = .true.
  ivar = ivar + 1
  CPDF(ivar)%name = 'W'
  call parser_read('PDF limits W',pdf_lim)
  CPDF(ivar)%min_y = pdf_lim(1); CPDF(ivar)%max_y = pdf_lim(2)
  CPDF(ivar)%init = .true.
  CPDF(ivar)%j_flip = .false.
  ivar = ivar + 1

  ! Filtered strain-rate
  do ii=1,3
     do jj=ii,3
        write(tmpstr, '(I1)') ii
        write(tmpstr2,'(I1)') jj
        CPDF(ivar)%name = 'S'//trim(tmpstr)//trim(tmpstr2)
        CPDF(ivar)%j_flip = .false.     
        ivar = ivar + 1
     end do
  end do

  ! Subfilter stress
  do ii=1,3
     do jj=ii,3
        write(tmpstr, '(I1)') ii
        write(tmpstr2,'(I1)') jj
        CPDF(ivar)%name = 'SFS'//trim(tmpstr)//trim(tmpstr2)
        CPDF(ivar)%j_flip = .false.     
        ivar = ivar + 1
     end do
  end do
  
  ! Scalars
  do isc=isc_O2,isc_H2O
     CPDF(ivar)%name = trim(SC_name(isc))
     CPDF(ivar)%min_y = SCval(isc,2)
     CPDF(ivar)%max_y = SCval(isc,1)
     CPDF(ivar)%init = .true.
     CPDF(ivar)%j_flip = .false.     
     ivar = ivar + 1
  end do

  ! Scalar gradients
  do isc=isc_O2,isc_H2O
     do ii=1,3
        write(tmpstr,'(I1)') ii
        tmpstr = 'G'//trim(tmpstr)//'_'//trim(SC_name(isc))
        CPDF(ivar)%name = tmpstr
        CPDF(ivar)%j_flip = .false.     
        ivar = ivar + 1
     end do
  end do
  
  ! Subfilter scalar flux
  do isc=isc_O2,isc_H2O
     do ii=1,3
        write(tmpstr,'(I1)') ii
        tmpstr = 'F'//trim(tmpstr)//'_'//trim(SC_name(isc))
        CPDF(ivar)%name = tmpstr
!!$        call parser_read('PDF limits '//trim(tmpstr),pdf_lim)
!!$        CPDF(ivar)%min_y = pdf_lim(1)
!!$        CPDF(ivar)%max_y = pdf_lim(2)
        if (ii.eq.2) then
           CPDF(ivar)%j_flip = .true.
        else
           CPDF(ivar)%j_flip = .false.
        end if
        ivar = ivar + 1
     end do
  end do

  ! Scalar flux/gradient alignment
  do isc=isc_O2,isc_H2O
     CPDF(ivar)%name = 'ASF_'//trim(SC_name(isc))
     CPDF(ivar)%min_y = -1.0_WP
     CPDF(ivar)%max_y =  1.0_WP
     CPDF(ivar)%init = .true.
     CPDF(ivar)%j_flip = .false.     
     ivar = ivar + 1
  end do

  ! Anisotropic strain-rate/anisotropic subfilter stress alignment
  do ii=1,3
     write(tmpstr,'(I1)') ii
     CPDF(ivar)%name = 'Avec_'//tmpstr
     CPDF(ivar)%min_y = -1.0_WP
     CPDF(ivar)%max_y =  1.0_WP
     CPDF(ivar)%init = .true.
     CPDF(ivar)%j_flip = .false.     
     ivar = ivar + 1
  end do     
  
  ! Allocate the workspaces and indices
  do ivar=1,nvar_pdf
     allocate(CPDF(ivar)%pdf(imin_:imax_,1:nbins_cond,1:nbins_cond))
     CPDF(ivar)%pdf = 0
     !allocate(CPDF(ivar)%index_x(1:nbins_cond))
     !allocate(CPDF(ivar)%index_y(1:nbins_cond))
     CPDF(ivar)%min_x = 0.0_WP
     CPDF(ivar)%max_x = 1.0_WP
  end do
  
  return
end subroutine volumeStats_les_init


! ========================================================== !
! Initialize the filter -- constant width                    !
! ========================================================== !
subroutine volumeStats_les_filter_init
  use volumeStats_les
  use volumeStats_metric
  implicit none

  real(WP) :: tempsum, dist, dx_loc, dy_loc, dz_loc
  integer :: i,j,k, flag
  integer :: a,b,c

  allocate(xfilt(imin_:imax_, imino_:imaxo_))                     !Note xfilt is *NOT* just a pointer
  allocate(yfilt(jmin_:jmax_, jmino_:jmaxo_))
  allocate(zfilt(kmin_-1:kmax_+1, kmino_:kmaxo_))

  allocate(xfbounds(imin_:imax_, 1:2))
  allocate(yfbounds(jmin_:jmax_, 1:2))
  allocate(zfbounds(kmino_:kmaxo_, 1:2))
 
  ! X-filter
  do i=imin_,imax_
     tempsum = 0.0_WP
     flag = 0
     do a=imino_,imaxo_
        !Find |x(i,j,k)-x0|
        
        if (a .EQ. imino_) then
           dx_loc = ABS(xm(a+1)-xm(a))
        else if(a .EQ. imaxo_) then
           dx_loc = ABS(xm(a)-xm(a-1))
        else
           dx_loc = 0.5_WP*ABS(xm(a+1)-xm(a-1))
        endif
        
        select case(trim(filter_type))
        case('homogeneous')
           dist = ABS(xm(a)-xm(i))
           
        case ('inhomogeneous')
           dist = ABS(xm(a)-xm(i))*dxmi(i)
!           dx_loc = abs(dxm(a))

        case default
           stop 'filter type not implemented'
        end select
        
        if(dist .LE. filter_width) then
           ! Point is within the filter radius
           ! Compute the filter weight
           xfilt(i, a) = dx_loc*exp(-6.0_WP*(dist/filter_width)**2)
           ! Compute running total of filter value
           tempsum = tempsum + xfilt(i,a)
           
           if (flag .EQ. 0) then
              xfbounds(i, 1) = a
              flag = 1 
           end if
           xfbounds(i,2) = a
           
        else
           !If outside filter radius, no contribution
           xfilt(i, a) = 0.0_WP
        end if
        
     end do
     !Normalize filter weights*dx_loc
     xfilt(i,:) = xfilt(i,:)/tempsum
  end do
  print *, 'xfbounds = ', xfbounds(imin_,:)

!!$  print *, 'xfbounds'
!!$  print '(3a10)', 'ifmin', 'i', 'ifmax'
!!$  do i=imin_,imax_
!!$     print '(3i10)', xfbounds(i,1), i, xfbounds(i,2)
!!$  end do
  
  ! Y-filter
  do j = jmin_, jmax_
     tempsum = 0.0_WP
     flag = 0
     do b = jmin_,jmax_ !jmino_, jmaxo_

        if(b .EQ. jmino_) then
           dy_loc = ABS(ym(b+1)-ym(b))
        else if(b .EQ. jmaxo_) then
           dy_loc = ABS(ym(b)-ym(b-1))
        else
           dy_loc = 0.5_WP*ABS(ym(b+1)-ym(b-1))
        endif
        
        select case(trim(filter_type))
        case('homogeneous')
           dist = ABS(ym(b)-ym(j))
           
        case ('inhomogeneous')
           dist = ABS(ym(b)-ym(j))*dymi(j)
!           dy_loc = abs(dym(b))
        case default
           stop 'filter type not implemented'
        end select

        if(dist .LE. filter_width) then
           yfilt(j, b) = dy_loc*exp(-6.0_WP*(dist/filter_width)**2)
           
           tempsum = tempsum + yfilt(j, b)
           if (flag .EQ. 0) then
              yfbounds(j, 1) = b
              flag = 1 
           end if
           yfbounds(j,2) = b
        else
           yfilt(j, b) = 0.0_WP
        end if
     end do
     yfilt(j,:) = yfilt(j,:)/tempsum
  end do

!!$  print *, 'yfbounds'
!!$  print '(3a10)', 'jfmin', 'j', 'jfmax'
!!$  do j=jmin_,jmax_
!!$     print '(3i10)', yfbounds(j,1), j, yfbounds(j,2)
!!$  end do
  
  ! Z-filter
  do k = kmin_-1, kmax_+1
     tempsum = 0.0_WP
     flag = 0
     dz_loc = z(2) - z(1)
     
     do c = kmino_, kmaxo_
        select case(trim(filter_type))
        case('homogeneous')
           dist = abs(real(k-c,WP))*dz_loc
        case ('inhomogeneous')
           dist = abs(real(k-c,WP))
        case default
           stop 'filter type not implemented'
        end select
        
        if(dist .LE. filter_width) then
           
           zfilt(k, c) = dz_loc*exp(-6.0_WP*(dist/filter_width)**2)
           
           tempsum = tempsum + zfilt(k, c)
           if (flag .EQ. 0) then
              zfbounds(k, 1) = c
              flag = 1 
           end if
           zfbounds(k,2) = c
        else
           zfilt(k, c) = 0.0_WP
        end if
     end do
     zfilt(k,:) = zfilt(k,:)/tempsum
  end do

!!$  print *, 'zfbounds'
!!$  print '(3a10)', 'kfmin', 'k', 'kfmax'
!!$  do j=kmin_,kmax_
!!$     print '(3i10)', zfbounds(j,1), j, zfbounds(j,2)
!!$  end do
  
  return
end subroutine volumeStats_les_filter_init



! ========================================================== !
! Initialize the filter -- constant ratio                    !
! ========================================================== !
subroutine volumeStats_les_filter_init_inhomogeneous
  use volumeStats_les
  use volumeStats_metric
  implicit none

  real(WP) :: tempsum, dist
  integer :: i,j,k, flag
  integer :: a,b,c

  allocate(xfilt(imin_:imax_, imino_:imaxo_))                     !Note xfilt is *NOT* just a pointer
  allocate(yfilt(jmin_:jmax_, jmino_:jmaxo_))
  allocate(zfilt(kmin_-1:kmax_+1, kmino_:kmaxo_))

  allocate(xfbounds(imin_:imax_, 1:2))
  allocate(yfbounds(jmin_:jmax_, 1:2))
  allocate(zfbounds(kmino_:kmaxo_, 1:2))
 
  ! X-filter
  do i=imin_,imax_
     tempsum = 0.0_WP
     flag = 0
     do a=imino_,imaxo_
        !Find |x(i,j,k)-x0|
        dist = ABS(xm(a)-xm(i))*dxmi(i)
        
        if(dist .LE. filter_width) then
           ! Point is within the filter radius
           ! Compute the filter weight
           xfilt(i, a) = abs(dxm(a))*exp(-6.0_WP*(dist/filter_width)**2)
           ! Compute running total of filter value
           tempsum = tempsum + xfilt(i,a)
           
           if (flag .EQ. 0) then
              xfbounds(i, 1) = a
              flag = 1 
           end if
           xfbounds(i,2) = a
           
        else
           !If outside filter radius, no contribution
           xfilt(i, a) = 0.0_WP
        end if
        
     end do
     !Normalize filter weights*dx
     xfilt(i,:) = xfilt(i,:)/tempsum
  end do
  print *, 'xfbounds = ', xfbounds(imin_,:)

  ! Y-filter
  do j = jmin_, jmax_
     tempsum = 0.0_WP
     flag = 0
     do b = jmino_, jmaxo_
        dist = ABS(ym(b)-ym(j))*dymi(j)
        if(dist .LE. filter_width) then
           yfilt(j, b) = abs(dym(b))*exp(-6.0_WP*(dist/filter_width)**2)
           
           tempsum = tempsum + yfilt(j, b)
           if (flag .EQ. 0) then
              yfbounds(j, 1) = b
              flag = 1 
           end if
           yfbounds(j,2) = b
        else
           yfilt(j, b) = 0.0_WP
        end if
     end do
     yfilt(j,:) = yfilt(j,:)/tempsum
  end do
!  print *, yfilt(ny/2,yfbounds(ny/2,1):yfbounds(ny/2,2))
  
  ! Z-filter
  do k = kmin_-1, kmax_+1
     tempsum = 0.0_WP
     flag = 0
     
     do c = kmino_, kmaxo_
        dist = abs(real(k-c,WP))
        if(dist .LE. filter_width) then
           
           zfilt(k, c) = abs(dz)*exp(-6.0_WP*(dist/filter_width)**2)
           
           tempsum = tempsum + zfilt(k, c)
           if (flag .EQ. 0) then
              zfbounds(k, 1) = c
              flag = 1 
           end if
           zfbounds(k,2) = c
        else
           zfilt(k, c) = 0.0_WP
        end if
     end do
     zfilt(k,:) = zfilt(k,:)/tempsum
  end do
!  print *, 'zfbounds = ', zfbounds(nz/2,:)
  
  return
end subroutine volumeStats_les_filter_init_inhomogeneous



! ========================================================== !
! Calculate conditional filtered variable                    !
! ========================================================== !
subroutine volumeStats_les_cond
  use volumeStats_les
  implicit none
  
  integer :: i, j, k
  real(WP) :: tmp1

  if (trim(cond_var).eq.'PROG') then
     !$OMP PARALLEL DO PRIVATE(j,i)
     do k=kmin_,kmax_,k_skip
        do j=jmin_,jmax_
           do i=imin_,imax_
              ! Define the progress variable based on O2
              PROGfilt(i,j,k) = (prog_upper - SCfilt(i,j,k,isc_PROG))/(prog_upper-prog_lower)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  else
     stop ' volumeStats_les_cond: conditional variable not implemented'
  end if

end subroutine volumeStats_les_cond


! ========================================================== !
! Time averaging in 3D                                       !
! ========================================================== !
subroutine volumeStats_les_timeavg(input, output)
  use volumeStats_les
  implicit none

  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_), intent(in)    :: input
  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_), intent(inout) :: output
  real(WP) :: prev_ntime, inv_curr_ntime
  integer  :: i, j, k

  prev_ntime = real(ntime_curr-1, WP)
  inv_curr_ntime = 1.0_WP/real(ntime_curr,WP)

  !$OMP PARALLEL PRIVATE(j, i)
  !$OMP DO 
  do k=kmin_,kmax_,k_skip
    do j=jmin_,jmax_
      do i=imin_,imax_
        output(i,j,k) = output(i,j,k)*prev_ntime  !Undoes time averaging
      enddo
    enddo
  enddo
  !$OMP END DO

  !$OMP DO
  do k=kmin_,kmax_,k_skip
    do j=jmin_,jmax_
      do i=imin_,imax_
        output(i,j,k) = output(i,j,k) + input(i,j,k)
      end do
    end do 
  end do
  !$OMP END DO

  !$OMP DO 
  do k=kmin_,kmax_,k_skip
    do j=jmin_,jmax_
      do i=imin_,imax_
        output(i,j,k) = output(i,j,k)*inv_curr_ntime !Time averaging
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL  

  return
end subroutine volumeStats_les_timeavg


! ========================================================== !
! Apply the filter                                           !
! ========================================================== !
subroutine volumeStats_les_filter(input, output, d_weight)
  use volumeStats_les
  implicit none

  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(in)    :: input
  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_), intent(out) :: output
  logical :: d_weight
  real(WP) :: prev_ntime, inv_curr_ntime
  integer  :: a, b, c, i, j, k

  output = 0.0_WP

  ! Call the filter kernel
!  call filter_kernel(input, output, 0)
  
!!$  prev_ntime = real(ntime_curr-1, WP)
!!$  inv_curr_ntime = 1.0_WP/real(ntime_curr,WP)
  
  !$OMP PARALLEL 
!!$  do k = kmin_, kmax_, skip
!!$     !$OMP DO PRIVATE(j, i)
!!$     do j = jmin_, jmax_
!!$        do i = imin_, imax_
!!$           output(i,j,k) = output(i,j,k)*prev_ntime  !Undoes time averaging
!!$        end do
!!$     end do
!!$     !$OMP END DO
!!$  end do

  
  do k = kmin_, kmax_, k_skip
     !$OMP DO PRIVATE(j, i, c, b, a)
     do j = jmin_, jmax_
        do i = imin_, imax_
           do c = zfbounds(k,1), zfbounds(k,2)
              do b = yfbounds(j,1), yfbounds(j,2)
                 do a = xfbounds(i,1), xfbounds(i,2)
                    output(i,j,k) = output(i,j,k) + xfilt(i,a)*yfilt(j,b)*zfilt(k,c)*input(a,b,c)
                 end do
              end do
           end do
        end do
     end do
     !$OMP END DO
  end do
  
  !$OMP END PARALLEL
!!$  ! Call the filter kernel
!!$  call filter_kernel(input,output, 0)
!!$
!!$  ! Filter the field with necessary offsets for gradients
!!$  if (need_offs) then
!!$     call filter_kernel(input,output,-1)
!!$     call filter_kernel(input,output,+1)
!!$  end if
!!$  !$OMP PARALLEL

  if (d_weight) then
     !$OMP PARALLEL
     do k = kmin_, kmax_, k_skip
        !$OMP DO PRIVATE(j, i)
        do j = jmin_, jmax_
           do i = imin_, imax_
              ! Undo the density weighting
              output(i,j,k) = output(i,j,k)/RHOfilt(i,j,k)
!!$           output(i,j,k) = output(i,j,k)*inv_curr_ntime !Time averaging
           end do
        end do
        !$OMP END DO
     end do
     !$OMP END PARALLEL
  end if
  
  return
end subroutine volumeStats_les_filter


! ========================================================== !
! Apply the filter including offsets for gradients           !
! ========================================================== !
subroutine volumeStats_les_filter_offs(input, output, d_weight)
  use volumeStats_les
  implicit none

  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(in)  :: input
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(out) :: output
  logical :: d_weight
  real(WP) :: prev_ntime, inv_curr_ntime
  integer  :: a, b, c, i, j, k

  output = 0.0_WP

  ! Call the filter kernel
  call filter_kernel_offs(input, output, -1, d_weight)
  call filter_kernel_offs(input, output,  0, d_weight)
  call filter_kernel_offs(input, output, +1, d_weight)
  
  return
end subroutine volumeStats_les_filter_offs


! ========================================================== !
! Filter all quanties of interest                            !
! ========================================================== !
subroutine volumeStats_les_filter_all
  use volumeStats_les
  implicit none

  integer :: i, j, k, ii, jj, isc, ivar
  !real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: tmpxyzo
  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_) :: tmpxyz
  real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3) :: tmpxyz3
  real(WP), dimension(1:3,1:3) :: evec_SFS, evec_Sij
  real(WP), dimension(1:3) :: eval
  real(WP), dimension(1:ny) :: tmpv1
  real(WP) :: tmp, trace
  integer :: ierr

  ! Filtered density and viscosity
  call volumeStats_les_filter_offs( &
       RHO (imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_), &
       RHOfilt(imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_), .false. )

  call volumeStats_les_filter( &
       RHO (imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_)* &
       VISC(imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_), &
       VISCfilt(imin_:imax_, jmin_:jmax_, kmin_:kmax_), .true. )
  if (irank.eq.iroot) print *, 'Done RHO and VISC'

  ! Filtered velocity
  do ii=1,3
     call volumeStats_les_filter_offs( &
          RHO(imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_) * &
          Ui (imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_, ii), &
          Uifilt(imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_, ii), .true. )
     if (irank.eq.iroot) print *, 'Done velocity', ii
  end do

  ! Filtered velocity gradients
  do ii=1,3
     call volumeStats_mpi_update_volume_border(Uifilt(:,:,:,ii))
     call volumeStats_metric_centered_gradient(tmpxyz3, Uifilt(:,:,:,ii), k_skip)
     do jj=1,3
        do k=kmin_,kmax_,k_skip
           do j=jmin_,jmax_
              do i=imin_,imax_
                 dUifiltdx(i,j,k,ii,jj) = tmpxyz3(i,j,k,jj)
              end do
           end do
        end do
     end do
  end do
  
  ! Filtered Strain Rate
  call volumeStats_metric_strainrate(dUifiltdx, Sijfilt)
  
  ! Filtered velocity correlations
  do jj=1,3
     do ii=jj,3
        call volumeStats_les_filter( &
             Ui (imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_, ii) * &
             Ui (imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_, jj) * &
             RHO(imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_), &
             UiUjfilt(imin_:imax_, jmin_:jmax_, kmin_:kmax_, ii, jj), .true. )
        if (irank.eq.iroot) print *, 'Done velocity correlation', ii, jj
     end do
  end do

  ! Can be eliminated (symmetric)
  UiUjfilt(:,:,:,1,2) = UiUjfilt(:,:,:,2,1)
  UiUjfilt(:,:,:,1,3) = UiUjfilt(:,:,:,3,1)
  UiUjfilt(:,:,:,2,3) = UiUjfilt(:,:,:,3,2)
  
  ! Filtered scalars
  !   Enable flag for offset filtered quantities -- necessary for gradients
  do isc=isc_O2,isc_H2O
     call volumeStats_les_filter_offs( &
          SC (imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_, isc) * &
          RHO(imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_), &
          SCfilt(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,isc), .true. )
     if (irank.eq.iroot) print *, 'Done scalar ', trim(SC_name(isc))
  end do

  ! Filtered scalar gradients
  do isc=isc_O2,isc_H2O
     call volumeStats_mpi_update_volume_border(SCfilt(:,:,:,isc))
     call volumeStats_metric_centered_gradient(dSCfiltdx(:,:,:,:,isc), SCfilt(:,:,:,isc), k_skip)
  end do
  
  ! Filtered velocity-scalar correlation
  do isc=isc_O2,isc_H2O
     do ii=1,3
        call volumeStats_les_filter( &
             RHO(imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_) * &
             Ui (imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_, ii) * &
             SC (imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_, isc), &
             UiSCfilt(imin_:imax_, jmin_:jmax_, kmin_:kmax_,ii,isc), .true. )
     end do
  end do
  PRINT *, irank, 'Done filtered velocities and scalars'
  
  ! Unfiltered anisotropic strain rate
  call volumeStats_metric_strainrate(dUdx, Sij)
  ! Remove the trace
  do ii=1,3
     !$OMP PARALLEL DO PRIVATE(j,i)
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              trace = dUdx(i,j,k,1,1) + dUdx(i,j,k,2,2) + dUdx(i,j,k,3,3)
              Sij(i,j,k,ii,ii) = Sij(i,j,k,ii,ii) - trace
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end do

!!$  do jj=1,3
!!$     do ii=1,3 
!!$        call volumeStats_mpi_update_volume_border(Sij(:,:,:,ii,jj))
!!$     end do
!!$  end do
 
  ! Unfiltered dissipation rate -- verify this
  eps = 0.0_WP
  do jj=1,3
     do ii=1,3
        !$OMP PARALLEL DO PRIVATE(j,i)
        do k=kmin_,kmax_,k_skip
           do j=jmin_,jmax_
              do i=imin_,imax_
                 eps(i,j,k) = eps(i,j,k) &
                      + Sij(i,j,k,ii,jj)*Sij(i,j,k,ii,jj)
              end do
           end do
        end do
        !$OMP END PARALLEL DO
     end do
  end do
  !$OMP PARALLEL DO PRIVATE(j,i)
  do k=kmin_,kmax_,k_skip
     do j=jmin_,jmax_
        do i=imin_,imax_
           eps(i,j,k) = eps(i,j,k) * 2.0_WP*VISC(i,j,k)/RHO(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  call volumeStats_mpi_update_volume_border(eps)


  ! Filtered dissipation rate
  call volumeStats_les_timeavg( &
       eps  (imin_:imax_, jmin_:jmax_, kmin_:kmax_), &
       epsTA(imin_:imax_, jmin_:jmax_, kmin_:kmax_) ) 
  

!!$  do jj=1,3
!!$     do ii=jj,3 
!!$        call volumeStats_les_filter( &
!!$             RHO(imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_) * &
!!$             Sij(imino_:imaxo_, jmino_:jmaxo_, kmino_:kmaxo_, ii, jj), &
!!$             Sijfilt(imin_:imax_, jmin_:jmax_, kmin_:kmax_, ii, jj), .true., .false. ) 
!!$     end do
!!$  end do
!!$  
!!$  ! Can eliminate this
!!$  Sijfilt(:,:,:,1,2) = Sijfilt(:,:,:,2,1)
!!$  Sijfilt(:,:,:,1,3) = Sijfilt(:,:,:,3,1)
!!$  Sijfilt(:,:,:,2,3) = Sijfilt(:,:,:,3,2)

  
!!$  ! Filter Damkohler number
!!$  call volumeStats_les_timeavg( &
!!$       tflaminar * (filter_width*filter_width/eps(imin_:imax_,jmin_:jmax_,kmin_:kmax_))**(1.0_WP/3.0_WP), &
!!$       Dafilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_) )
  
  PRINT *, irank, 'Filtering step done'

  ! Update the filtered progress variable and bin_index
  call cond_map(SCfilt  (imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:nscalar), &
                PROGfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  
  ! Update the density-weighted pdfs
  ivar = 1
  call cond_pdf_kernel(CPDF(ivar), RHOfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  ivar = ivar + 1
  
  ! Velocity
  do ii=1,3
     tmpxyz = Uifilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)
     call cond_pdf_kernel(CPDF(ivar), tmpxyz)
     ivar = ivar + 1
  end do

  ! Filtered strain-rate
  do ii=1,3
     do jj=ii,3
        tmpxyz = Sijfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,jj)
        call cond_pdf_kernel(CPDF(ivar), tmpxyz)
        ivar = ivar + 1
     end do
  end do
  
  ! Subfilter stress
  do ii=1,3
     do jj=ii,3
        tmpxyz = UiUjfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,jj) &
             - Uifilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)*Uifilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,jj)
        call cond_pdf_kernel(CPDF(ivar), tmpxyz)
        ivar = ivar + 1
     end do
  end do
  
  ! Scalars
  do isc=isc_O2,isc_H2O
     call cond_pdf_kernel(CPDF(ivar), SCfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc))
     ivar = ivar + 1
  end do

  ! Scalar gradients
  do isc=isc_O2,isc_H2O
     do ii=1,3
        call cond_pdf_kernel(CPDF(ivar), dSCfiltdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,isc))
        ivar = ivar + 1
     end do
  end do

  ! Subfilter scalar flux
  do isc=isc_O2,isc_H2O
     do ii=1,3
        tmpxyz = UiSCfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,isc) - &
             Uifilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)*SCfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc)
        call cond_pdf_kernel(CPDF(ivar), tmpxyz)
        ivar = ivar + 1
     end do
  end do

  ! Scalar flux/gradient alignment
  do isc=isc_O2,isc_H2O
     do ii=1,3
        tmpxyz3(:,:,:,ii) = UiSCfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii,isc) - &
             Uifilt  (imin_:imax_,jmin_:jmax_,kmin_:kmax_,ii)*SCfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc)
     end do
     tmpxyz = &
          sum( tmpxyz3 * dSCfiltdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3,isc), dim=4 ) &
             / ( sqrt(sum( tmpxyz3**2, dim=4)) &
               * sqrt(sum(dSCfiltdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3,isc)**2, dim=4)) )
     call cond_pdf_kernel(CPDF(ivar), tmpxyz)
     ivar = ivar + 1
  end do
  
  if (irank.eq.iroot) print *, '     -- Done velocity and scalar pdfs '


  ! Subfilter stress and strain-rate eigenvalues
  tmpxyz3 = 0.0_WP
  
  do k=kmin_,kmax_,k_skip
     do j=jmin_,jmax_
        do i=imin_,imax_
           
           ! Eigenvectors of anisotropic subfilter stress
           do ii=1,3
              do jj=ii,3
                 evec_SFS(ii,jj) = UiUjfilt(i,j,k,ii,jj) - Uifilt(i,j,k,ii)*Uifilt(i,j,k,jj)
              end do
           end do
           trace = evec_SFS(1,1) + evec_SFS(2,2) + evec_SFS(3,3)
           do ii=1,3
              evec_SFS(ii,ii) = evec_SFS(ii,ii) - trace
           end do
           call DSYEV('V', 'U', 3, evec_SFS, 3, eval, tmpv1, ny, ierr)
           if (ierr.ne.0) print *, 'Something went wrong with DSYEV -- SFF'
           
           ! Eigenvectors of anisotropic filtered strain-rate
           do ii=1,3
              do jj=ii,3
                 evec_Sij(ii,jj) = Sijfilt(i,j,k,ii,jj)
              end do
           end do
           trace = evec_Sij(1,1) + evec_Sij(2,2) + evec_Sij(3,3)
           do ii=1,3
              evec_Sij(ii,ii) = evec_Sij(ii,ii) - trace
           end do
           call DSYEV('V', 'U', 3, evec_Sij, 3, eval, tmpv1, ny, ierr)
           if (ierr.ne.0) print *, 'Something went wrong with DSYEV -- Sij'

           ! Compute alignment
           do ii=1,3
              tmp = 0.0_WP
              do jj=1,3
                 tmp = tmp + evec_SFS(ii,jj) * evec_Sij(4-ii,jj)
              end do
              tmpxyz3(i,j,k,ii) = tmp
           end do

        end do
     end do
  end do

  ! Compute stress/strain-rate alignment pdfs
  do ii=1,3
     call cond_pdf_kernel(CPDF(ivar), tmpxyz3(:,:,:,ii))
     ivar = ivar + 1
  end do

  if (irank.eq.iroot) print *, '     -- Done eigenvector alignment pdfs '
  
  return
end subroutine volumeStats_les_filter_all


! ========================================================== !
! Calculate final quantities, output them                    !
! ========================================================== !
subroutine volumeStats_les_finalize
  use volumeStats_les
  use parser
  implicit none

  integer  :: i, j, k, isc, ii, jj, iunit, iplane, ivar
  real(WP) :: nvec, Sijmag, const
  character(len=str_short) :: of_name
 

!!$  !Calculate Modeled SFF in y (2) direction
!!$  !$OMP PARALLEL PRIVATE(j, i, nvec, isc, const, Sijmag, ii, jj)
!!$  !$OMP DO 
!!$  do k = kmin_, kmax_,k_skip
!!$     do j = jmin_, jmax_
!!$        do i = imin_, imax_
!!$           !nvec = dSCfiltdx(i,j,k,isc_T,1)/SQRT(dSCfiltdx(i,j,k,isc_T,1)*dSCfiltdx(i,j,k,isc_T,1) &
!!$           !                                   + dSCfiltdx(i,j,k,isc_T,2)*dSCfiltdx(i,j,k,isc_T,2) &
!!$           !                                   + dSCfiltdx(i,j,k,isc_T,3)*dSCfiltdx(i,j,k,isc_T,3) )
!!$           !const = SQRT(VISCfilt(i,j,k)*Dafilt(i,j,k)*tflaminar/RHOfilt(i,j,k))
!!$
!!$           Sijmag = 0.0_WP
!!$           do jj = 1,3
!!$              do ii=1,3
!!$                 Sijmag = Sijmag + 2.0_WP*Sijfilt(i,j,k,ii,jj)*Sijfilt(i,j,k,ii,jj)
!!$              end do
!!$           end do
!!$           Sijmag = SQRT(Sijmag)
!!$
!!$           do isc =1,nscalar
!!$              SFFmodel(i,j,k,isc,1) = (SCfilt(i,j,k,isc)-SCval(isc,unburned))*(SCval(isc,burned)-SCfilt(i,j,k,isc)) &
!!$                   /(SCval(isc,burned)-SCval(isc,unburned))
!!$              SFFmodel(i,j,k,isc,2) = Sijmag
!!$           end do
!!$
!!$        end do
!!$     end do
!!$  end do
!!$  !$OMP END DO
!!$  !$OMP END PARALLEL

  !Filter Damkohler number
!!$OMP PARALLEL PRIVATE(j, i)
!!$OMP DO 
  !do k = kmin_, kmax_,k_skip
  !  do j = jmin_, jmax_
  !    do i = imin_, imax_
  !      Dafilt(i,j,k) = tflaminar * (filter_width*filter_width/epsfilt(i,j,k))**(1.0_WP/3.0_WP)
  !    enddo
  !  enddo
  !end do
!!$OMP END DO
!!$OMP END PARALLEL

  !Calculate Kolmogorov Length Scale
!!$OMP PARALLEL PRIVATE(j, i)
!!$OMP DO 
  !  do k = kmin_, kmax_,k_skip
  !    do j = jmin_, jmax_
  !      do i = imin_, imax_
  !        ETAkol(i,j,k) = SQRT(SQRT(VISC(i,j,k)*VISC(i,j,k)*VISC(i,j,k)/(eps(i,j,k)*RHO(i,j,k)*RHO(i,j,k)*RHO(i,j,k))))
  !      enddo
  !    enddo
  !  end do
!!$OMP END DO
!!$OMP END PARALLEL


!!$  !Compute UiYj_(bar) -  Ui_(bar)*Yj_(bar)
!!$  do isc = 1, nscalar
!!$     do ii = 1, 3
!!$        !$OMP PARALLEL PRIVATE(j, i)
!!$        !$OMP DO 
!!$        do k = kmin_, kmax_,k_skip
!!$           do j = jmin_, jmax_
!!$              do i = imin_, imax_
!!$                 SFF(i,j,k,ii,isc) = UiSCfilt(i,j,k,ii,isc) - Uifilt(i,j,k,ii)*SCfilt(i,j,k,isc) 
!!$              end do
!!$           end do
!!$        end do
!!$        !$OMP END DO
!!$        !$OMP END PARALLEL
!!$     enddo
!!$  enddo

  ! do isc = 1, nscalar
  !     !$OMP PARALLEL PRIVATE(j, i, nvec)
  !     !$OMP DO 
  !     do k = kmin_, kmax_,k_skip
  !       do j = jmin_, jmax_
  !         do i = imin_, imax_
  !           nvec = 1/SQRT(dSCfiltdx(i,j,k,isc_T,1)*dSCfiltdx(i,j,k,isc_T,1) &
  !                       + dSCfiltdx(i,j,k,isc_T,2)*dSCfiltdx(i,j,k,isc_T,2) &
  !                       + dSCfiltdx(i,j,k,isc_T,3)*dSCfiltdx(i,j,k,isc_T,3) )
  !           SFFnorm(i,j,k,isc) =( SFF(i,j,k,1,isc)*dSCfiltdx(i,j,k,isc_T,1) &
  !                               + SFF(i,j,k,2,isc)*dSCfiltdx(i,j,k,isc_T,2) &
  !                               + SFF(i,j,k,3,isc)*dSCfiltdx(i,j,k,isc_T,3) )*nvec
  !         end do
  !       end do
  !     end do
  !     !$OMP END DO
  !     !$OMP END PARALLEL
  ! enddo

!!$  !Compute UiUj_(bar) -  Ui_(bar)*Uj_(bar)
!!$  do jj = 1, 3
!!$     do ii = 1, 3
!!$        !$OMP PARALLEL PRIVATE(j, i)
!!$        !$OMP DO 
!!$        do k = kmin_, kmax_,k_skip
!!$           do j = jmin_, jmax_
!!$              do i = imin_, imax_
!!$                 SFS(i,j,k,ii,jj) = UiUjfilt(i,j,k,ii,jj) - Uifilt(i,j,k,ii)*Uifilt(i,j,k,jj)
!!$              end do
!!$           end do
!!$        end do
!!$        !$OMP END DO
!!$        !$OMP END PARALLEL
!!$     enddo
!!$  enddo

  iunit = 10
  !Dump data
  do iplane = pnmin_, pnmax_

     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_les_basic_data_'//trim(of_name), action='write')
     i = pnindx(iplane)

     do j = jmin_, jmax_
        do k = kmin_+k_skip, kmax_-k_skip,k_skip
           write (iunit,'(ES22.13)',advance='no'), ym(j) !y-coordinate
           write (iunit,'(ES22.13)',advance='no'), zm(k) ! z-coordinate
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*PROGfilt(i-1:i,j,k))
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*epsTA(i-1:i,j,k))
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*VISCfilt(i-1:i,j,k))
           write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*rhofilt(i-1:i,j,k))

           do ii=1,3
              if(ii .eq. 3) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*Uifilt(i-1:i,j,k,ii))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*Uifilt(i-1:i,j,k,ii))
              endif
           end do
        end do
     end do
     close(iunit)
  end do



  ! Filtered scalars
  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_les_SC_data_'//trim(of_name), action='write')
     i = pnindx(iplane)

     do j = jmin_, jmax_
        do k = kmin_+k_skip, kmax_-k_skip,k_skip
           write (iunit,'(ES22.13)',advance='no'), ym(j) !y-coordinate
           write (iunit,'(ES22.13)',advance='no'), zm(k) ! z-coordinate
           do isc=1,nscalar
              if(isc .eq. nscalar) then
                 write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*SCfilt(i-1:i,j,k,isc))
              else
                 write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SCfilt(i-1:i,j,k,isc))
              endif
           end do
        end do
     end do
     close(iunit)
  end do

  
  ! Filtered scalar gradients
  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_les_dSCdx_data_'//trim(of_name), action='write')
     i = pnindx(iplane)

     do j = jmin_, jmax_
        do k = kmin_+k_skip, kmax_-k_skip,k_skip
           write (iunit,'(ES22.13)',advance='no'), ym(j) !y-coordinate
           write (iunit,'(ES22.13)',advance='no'), zm(k) ! z-coordinate
           do isc=1,nscalar
              do ii=1,3
                 if (isc.eq.nscalar .and. ii.eq.3) then
                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*dSCfiltdx(i-1:i,j,k,ii,isc))
                 else
                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dSCfiltdx(i-1:i,j,k,ii,isc))
                 end if
              end do
           end do
        end do
     end do
     close(iunit)
  end do
  
  ! Filtered velocity-scalar correlation
  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_les_UiSC_data_'//trim(of_name), action='write')
     i = pnindx(iplane)

     do j = jmin_, jmax_
        do k = kmin_+k_skip, kmax_-k_skip, k_skip
           write (iunit,'(ES22.13)',advance='no'), ym(j) !y-coordinate
           write (iunit,'(ES22.13)',advance='no'), zm(k) ! z-coordinate
           do isc=1,nscalar
              do ii=1,3
                 if((ii .eq. 3) .and. (isc.eq.nscalar)) then
                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*UiSCfilt(i-1:i,j,k,ii,isc))
                 else
                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*UiSCfilt(i-1:i,j,k,ii,isc))
                 endif
              end do
           end do
        end do
     end do
     close(iunit)
  end do
  
  ! Filtered velocity gradient
  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_les_dUdx_data_'//trim(of_name), action='write')
     i = pnindx(iplane)

     do j = jmin_, jmax_
        do k = kmin_+k_skip, kmax_-k_skip,k_skip
           write (iunit,'(ES22.13)',advance='no'), ym(j) !y-coordinate
           write (iunit,'(ES22.13)',advance='no'), zm(k) ! z-coordinate
           do jj=1,3
              do ii=1,3
                 if((ii .eq. 3) .and. (jj .eq. 3)) then
                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*dUifiltdx(i-1:i,j,k,ii,jj))
                 else
                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*dUifiltdx(i-1:i,j,k,ii,jj))
                 endif
              end do
           end do
        end do
     end do
     close(iunit)
  end do
  
  ! Filtered strain rate
  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_les_Sij_data_'//trim(of_name), action='write')
     i = pnindx(iplane)

     do j = jmin_, jmax_
        do k = kmin_+k_skip, kmax_-k_skip,k_skip
           write (iunit,'(ES22.13)',advance='no'), ym(j) !y-coordinate
           write (iunit,'(ES22.13)',advance='no'), zm(k) ! z-coordinate
           do jj=1,3
              do ii=1,3
                 if((ii .eq. 3) .and. (jj .eq. 3)) then
                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*Sijfilt(i-1:i,j,k,ii,jj))
                 else
                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*Sijfilt(i-1:i,j,k,ii,jj))
                 endif
              end do
           end do
        end do
     end do
     close(iunit)
  end do
  
  ! Filtered velocity double correlation
  do iplane=pnmin_,pnmax_
     call get_name(iplane,of_name)
     open(unit=iunit, file=trim(output_name)//'_les_UiUj_data_'//trim(of_name), action='write')
     i = pnindx(iplane)

     do j = jmin_, jmax_
        do k = kmin_+k_skip, kmax_-k_skip, k_skip
           write (iunit,'(ES22.13)',advance='no'), ym(j) !y-coordinate
           write (iunit,'(ES22.13)',advance='no'), zm(k) ! z-coordinate
           do jj=1,3
              do ii=1,3
                 if((ii .eq. 3) .and. (jj .eq. 3)) then
                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*UiUjfilt(i-1:i,j,k,ii, jj))
                 else
                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*UiUjfilt(i-1:i,j,k,ii, jj))
                 endif
              end do
           end do
        end do
     end do
     close(iunit)
  end do


  ! Write conditional pdfs
  do ivar=1,nvar_pdf
     
     do iplane=pnmin_,pnmax_
        call get_name(iplane,of_name)
        open(unit=iunit, file=trim(output_name)//'_les_pdf_'//trim(CPDF(ivar)%name)//'_'//trim(of_name), &
             action='write')
        i = pnindx(iplane)

        ! Write header
        write (iunit,'(ES22.13)',advance='no'), 0.0_WP
        do jj=1,nbins_cond
           if (jj.eq.nbins_cond) then
              write (iunit,'(ES22.13)'), &
                   real(jj-1,WP)*(CPDF(ivar)%max_y - CPDF(ivar)%min_y)/real(nbins_cond,WP) + CPDF(ivar)%min_y
           else
              write (iunit,'(ES22.13)',advance='no'), &
                   real(jj-1,WP)*(CPDF(ivar)%max_y - CPDF(ivar)%min_y)/real(nbins_cond,WP) + CPDF(ivar)%min_y
           end if
        end do

        ! Write data
        do ii=1,nbins_cond
           write (iunit,'(ES22.13)',advance='no'), bins_cond(ii) ! conditional variable
           do jj=1,nbins_cond
              if (jj.eq.nbins_cond) then
                 write (iunit,'(ES22.13)'), real(CPDF(ivar)%pdf(i,ii,jj),WP)
              else
                 write (iunit,'(ES22.13)',advance='no'), real(CPDF(ivar)%pdf(i,ii,jj),WP)
              end if
           end do
        end do
        
        close(iunit)
     end do
     
  end do
  
  
!!$  do iplane=pnmin_,pnmax_
!!$
!!$     call get_name(iplane,of_name)
!!$     open(unit=iunit, file=trim(output_name)//'_les_SFS_data_'//trim(of_name), action='write')
!!$     i = pnindx(iplane)
!!$
!!$     do j = jmin_, jmax_
!!$        do k = kmin_+k_skip, kmax_-k_skip,k_skip
!!$           write (iunit,'(ES22.13)',advance='no'), ym(j) !y-coordinate
!!$           write (iunit,'(ES22.13)',advance='no'), zm(k) ! z-coordinate
!!$           do jj=1,3
!!$              do ii=1,3
!!$                 if((ii .eq. 3) .and. (jj .eq. 3)) then
!!$                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*SFS(i-1:i,j,k,ii,jj))
!!$                 else
!!$                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SFS(i-1:i,j,k,ii,jj))
!!$                 endif
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$     close(iunit)
!!$  end do


!!$  do iplane=pnmin_,pnmax_
!!$
!!$     call get_name(iplane,of_name)
!!$     open(unit=iunit, file=trim(output_name)//'_les_SFF_data_'//trim(of_name), action='write')
!!$     i = pnindx(iplane)
!!$
!!$     do j = jmin_, jmax_
!!$        do k = kmin_+k_skip, kmax_-k_skip,k_skip
!!$           write (iunit,'(ES22.13)',advance='no'), ym(j) !y-coordinate
!!$           write (iunit,'(ES22.13)',advance='no'), zm(k) ! z-coordinate
!!$           do isc=isc_CASE,isc_CASE
!!$              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SFFmodel(i-1:i,j,k,isc,1))
!!$              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SFFmodel(i-1:i,j,k,isc,2))
!!$              write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*Dafilt(i-1:i,j,k))
!!$              do ii=1,3
!!$                 if((isc .eq. isc_CASE) .and. (ii .eq. 3)) then
!!$                    write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*SFF(i-1:i,j,k,ii,isc))
!!$                 else
!!$                    write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*SFF(i-1:i,j,k,ii,isc))
!!$                 endif
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$     close(iunit)
!!$  end do

  print *, irank, 'Done output'
  
end subroutine volumeStats_les_finalize


! ========================================================== !
! Output parallel data decomp for Ensight visualiziation     !
! ========================================================== !

subroutine volumeStats_les_ensight
  use volumeStats_les
  use fileio
  use parser
  implicit none

  integer :: i, j, k, isc
  integer :: iunit, ierr
  character(len=str_medium) :: fname, folder, filename
  integer :: nx1, ny1, nz1, nvarsT, ntime1

  !Ensight variables
  real(SP), dimension(:), pointer :: xm1, ym1, zm1
  integer, dimension(:,:,:), pointer :: iblank
  real(SP), dimension(:,:,:), pointer :: var4
  real(WP) :: ddum
  integer, dimension(:), pointer :: var_list
  integer :: idum, ipart
  character(len=80) :: cdum, buffer
  logical :: file_is_there

  !MPI variables
  integer, dimension(3) :: gsizes,lsizes,start
  integer :: fileview,datasize,gdatasize
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: disp

  gsizes(1) = nx
  gsizes(2) = ny
  gsizes(3) = nz
  lsizes(1) = nx_
  lsizes(2) = ny_
  lsizes(3) = nz_
  start(1) = imin_-imin
  start(2) = jmin_-jmin
  start(3) = kmin_-kmin
  datasize = lsizes(1)*lsizes(2)*lsizes(3)
  gdatasize= gsizes(1)*gsizes(2)*gsizes(3)

  call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,start,MPI_ORDER_FORTRAN,MPI_REAL_SP,fileview,ierr)
  call MPI_TYPE_COMMIT(fileview,ierr)

  nx1 = nx_
  ny1 = ny_
  nz1 = nz_

  if (irank .eq. iroot) print *, 'Starting ensight output'
  ! JFM: might have to change MPI_COMM_WORLD to comm
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call parser_read('Ensight output', folder)
  folder = trim(folder) // '/'
  call system("mkdir -p " // trim(folder))
  allocate(var4(nx1,ny1,nz1))

  write(fname, '(''vol_data.1_'', I8.8)') 
  filename = "volume_data/"//trim(fname)
  print*, irank, trim(filename)

  !++++++++++++++++++++++++++++++++++++ CASE FILE++++++++++++++++++++++++++++++++++++++!
  if (irank .eq. iroot) then
     ! ** Write the case file **
     iunit = 2
     fname = 'arts.case'
     filename = trim(folder)//trim(fname)        
     open(iunit,file=trim(filename),form="formatted",iostat=ierr,status="REPLACE")
     !print*, 'File opened:', filename
     ! Write the case
     cdum = 'FORMAT'
     write(iunit,'(a80)') cdum
     cdum = 'type: ensight gold'
     write(iunit,'(a80)') cdum
     cdum = 'GEOMETRY'
     write(iunit,'(a80)') cdum
     cdum = 'model: 1 1 geometry'
     write(iunit,'(a80)') cdum
     cdum = 'VARIABLE'
     write(iunit,'(a80)') cdum 

     cdum = 'vector per node: V_filt V_filt'
     write(iunit,'(a80)') cdum

     do i = 4, nvars
        select case(trim(names(i)))
        case ('P')
           buffer = trim(adjustl('C')) // '_filt' !Conditional variable
        case default
           buffer = trim(adjustl(names(i))) // '_filt'
        end select
        cdum = 'scalar per node: '//trim(buffer)//' '//trim(buffer)
        write(iunit,'(a80)') cdum
     end do

     cdum = 'TIME'
     write(iunit,'(a80)') cdum 
     cdum = 'time set: 1'
     write(iunit,'(a80)') cdum 
     cdum = 'number of steps: 1'
     write(iunit,'(a80)') cdum         
     cdum = 'filename start number: 1'
     write(iunit,'(a80)') cdum        
     cdum = 'filename increment: 1'
     write(iunit,'(a80)') cdum
     cdum = 'time values:'
     ddum = time(nfiles)
     write(iunit, '(a12,ES12.5)') cdum,ddum

     ! Close the file
     close(iclose(iunit))
     !print *, irank, 'arts.case completed' 

     !+++++++++++++++++++++++++++++++++++++++++GEOMETRY+++++++++++++++++++++++++++++++++++++++++++!

     fname = 'geometry'
     filename = trim(folder)//trim(fname)
     iunit = 10
     ! ** Open the grid file to write **
     call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
     !print*, 'File opened:', filename
     ! Write the geometry
     buffer = 'C Binary'
     call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
     buffer = 'Ensight Gold Geometry File'
     call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
     buffer = 'Structured Geometry'
     call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
     buffer = 'node id off'
     call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
     buffer = 'element id off'
     call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
     ! Cell centers
     buffer = 'part'
     call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
     ipart = 1
     call BINARY_FILE_WRITE(iunit,ipart,1,kind(ipart),ierr)
     buffer = 'Complete geometry'
     call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
     buffer = 'block rectilinear iblanked'
     call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
     !print*, 'Header completed'

     call BINARY_FILE_WRITE(iunit,nx,1,kind(nx1),ierr)
     call BINARY_FILE_WRITE(iunit,ny,1,kind(ny1),ierr)
     call BINARY_FILE_WRITE(iunit,nz,1,kind(nz1),ierr)
     !print*, 'N points completed'

     allocate(xm1(1:nx),ym1(1:ny),zm1(1:nz))
     allocate(iblank(nx,ny,nz))
     iblank = 1
     xm1(1:nx) = x(imin:imax)
     ym1(1:ny) = y(jmin:jmax)
     zm1(1:nz) = z(kmin:kmax)

     call BINARY_FILE_WRITE(iunit,xm1,nx,kind(xm1),ierr)
     call BINARY_FILE_WRITE(iunit,ym1,ny,kind(ym1),ierr)
     call BINARY_FILE_WRITE(iunit,zm1,nz,kind(zm1),ierr)
     call BINARY_FILE_WRITE(iunit,iblank,nx*ny*nz,kind(iblank),ierr)

     call BINARY_FILE_CLOSE(iunit,ierr)
     !print *, irank, 'geometry completed'
     deallocate(xm1, ym1, zm1, iblank)

  end if

  !++++++++++++++++++++++++++++++++++++ DATA FILES+++++++++++++++++++++++++++++++++++++!

  do i = 1, 3
     select case(trim(names(i)))
     case ('U')
        fname = trim(folder)//'V_filt' !vector of velocity

        inquire(file=trim(fname),exist=file_is_there)

        if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(trim(fname),mpi_info,ierr)
        call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_WRONLY+MPI_MODE_CREATE,mpi_info,iunit,ierr)

        if(irank .eq. iroot) then
           !call BINARY_FILE_OPEN(iunit,trim(fname),"w",ierr)
           !print *, irank, iroot, '573'
           cdum = 'velocity'
           !call BINARY_FILE_WRITE(iunit,cdum,80,kind(cdum),ierr)
           call MPI_FILE_WRITE(iunit,cdum,80,MPI_CHARACTER,status,ierr)
           cdum = 'part'
           !call BINARY_FILE_WRITE(iunit,cdum,80,kind(cdum),ierr)
           call MPI_FILE_WRITE(iunit,cdum,80,MPI_CHARACTER,status,ierr)
           idum = 1
           !call BINARY_FILE_WRITE(iunit,idum,1,kind(idum),ierr)
           call MPI_FILE_WRITE(iunit,idum,1,MPI_INTEGER,status,ierr)
           cdum = 'block'
           !call BINARY_FILE_WRITE(iunit,cdum,80,kind(cdum),ierr)
           call MPI_FILE_WRITE(iunit,cdum,80,MPI_CHARACTER,status,ierr)

        end if

     case default
     end select

     select case(trim(names(i)))
     case ('U')
        var4(1:nx1,1:ny1,1:nz1) = Uifilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1)
     case ('V')
        var4(1:nx1,1:ny1,1:nz1) = Uifilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,2)
     case ('W')
        var4(1:nx1,1:ny1,1:nz1) = Uifilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,3)
     case default
     end select

     disp = 3*80+(i-1)*gdatasize*4+4
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview,"native",mpi_info,ierr)
     call MPI_FILE_WRITE_ALL(iunit,var4,datasize,MPI_REAL_SP,status,ierr)

     !call BINARY_FILE_WRITE(iunit,var4,nx1*ny1*nz1,kind(var4),ierr)
  end do

  call MPI_FILE_CLOSE(iunit,ierr)
  !call BINARY_FILE_CLOSE(iunit,ierr)  

  isc = 1
  do i = 4,nvars
     select case(trim(names(i)))
     case ('P')
        fname = trim(folder)//trim(adjustl('C')) // '_filt' !Conditional variable
     case default
        fname = trim(folder)//trim(adjustl(names(i))) // '_filt'
     end select

     inquire(file=trim(fname),exist=file_is_there)
     if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(trim(fname),mpi_info,ierr)

     call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_WRONLY+MPI_MODE_CREATE,mpi_info,iunit,ierr)

     if(irank .eq. iroot) then 
        cdum = trim(adjustl(names(i)))
        !call BINARY_FILE_WRITE(iunit,cdum,80,kind(cdum),ierr)
        call MPI_FILE_WRITE(iunit,cdum,80,MPI_CHARACTER,status,ierr)
        cdum = 'part'
        !call BINARY_FILE_WRITE(iunit,cdum,80,kind(cdum),ierr)
        call MPI_FILE_WRITE(iunit,cdum,80,MPI_CHARACTER,status,ierr)
        idum = 1
        !call BINARY_FILE_WRITE(iunit,idum,1,kind(idum),ierr)
        call MPI_FILE_WRITE(iunit,idum,1,MPI_INTEGER,status,ierr)
        cdum = 'block'
        !call BINARY_FILE_WRITE(iunit,cdum,80,kind(cdum),ierr)
        call MPI_FILE_WRITE(iunit,cdum,80,MPI_CHARACTER,status,ierr)
     end if

     select case(trim(names(i)))
     case ('P')
        var4(1:nx1,1:ny1,1:nz1) = PROGfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
     case ('RHO')
        var4(1:nx1,1:ny1,1:nz1) = rhofilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
     case('HR')
        ! Do nothing
     case default
        ! Read scalars here
        var4(1:nx1,1:ny1,1:nz1) = SCfilt(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc)
        isc = isc + 1
     end select

     disp = 3*80+4
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview,"native",mpi_info,ierr)
     call MPI_FILE_WRITE_ALL(iunit,var4,datasize,MPI_REAL_SP,status,ierr)
     call MPI_FILE_CLOSE(iunit,ierr)

     !call BINARY_FILE_WRITE(iunit,var4,nx1*ny1*nz1,kind(var4),ierr)
     !call BINARY_FILE_CLOSE(iunit,ierr)   

  end do

  print *, irank, 'Completed ensight dump'

end subroutine volumeStats_les_ensight
