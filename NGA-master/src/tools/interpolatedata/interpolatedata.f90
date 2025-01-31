module interpolatedata
  use precision
  use string
  use fileio
  implicit none
  
  ! Sizes
  integer :: nx1,ny1,nz1
  integer :: nx2,ny2,nz2
  integer :: xper1,yper1,zper1
  integer :: xper2,yper2,zper2
  integer :: icyl1
  integer :: icyl2
  integer :: nvar
  real(WP) :: dt,time
  character(len=str_short), dimension(:), pointer :: names
  character(len=str_medium) :: config1
  character(len=str_medium) :: config2
  integer, dimension(:,:), pointer :: mask1
  integer, dimension(:,:), pointer :: mask2
  
  ! Data
  real(WP), dimension(:,:,:), pointer :: data1
  real(WP), dimension(:,:,:), pointer :: data2
  
  ! Mesh
  real(WP), dimension(:), pointer :: x1,xm1
  real(WP), dimension(:), pointer :: y1,ym1
  real(WP), dimension(:), pointer :: z1,zm1
  real(WP), dimension(:), pointer :: x2,xm2
  real(WP), dimension(:), pointer :: y2,ym2
  real(WP), dimension(:), pointer :: z2,zm2
  
end module interpolatedata


program interpolator
  use interpolatedata
  implicit none
  
  character(len=str_medium) :: fconfig1
  character(len=str_medium) :: fconfig2
  character(len=str_medium) :: fdata1
  character(len=str_medium) :: fdata2
  integer :: iunit,iunit1,iunit2,var
  integer :: i,j,k,ierr
  
  ! Get file names from user
  print*,'================================'
  print*,'| ARTS - Data file interpolator |'
  print*,'================================'
  print*
  print "(a28,$)", " Enter source config file : "
  read "(a)", fconfig1
  print "(a33,$)", " Enter destination config file : "
  read "(a)", fconfig2
  print "(a26,$)", " Enter source data file : "
  read "(a)", fdata1
  print "(a31,$)", " Enter destination data file : "
  read "(a)", fdata2
  
  ! Read the source config file
  call BINARY_FILE_OPEN(iunit,trim(fconfig1),"r",ierr)
  call BINARY_FILE_READ(iunit,config1,str_medium,kind(config1),ierr)
  call BINARY_FILE_READ(iunit,icyl1,1,kind(icyl1),ierr)
  call BINARY_FILE_READ(iunit,xper1,1,kind(xper1),ierr)
  call BINARY_FILE_READ(iunit,yper1,1,kind(yper1),ierr)
  call BINARY_FILE_READ(iunit,zper1,1,kind(zper1),ierr)
  call BINARY_FILE_READ(iunit,nx1,1,kind(nx1),ierr)
  call BINARY_FILE_READ(iunit,ny1,1,kind(ny1),ierr)
  call BINARY_FILE_READ(iunit,nz1,1,kind(nz1),ierr)
  print*,'Source file'
  print*,'Config : ',config1
  print*,'Grid :',nx1,'x',ny1,'x',nz1
  allocate(x1(nx1+1),y1(ny1+1),z1(nz1+1))
  allocate(mask1(nx1,ny1))
  call BINARY_FILE_READ(iunit,x1,nx1+1,kind(x1),ierr)
  call BINARY_FILE_READ(iunit,y1,ny1+1,kind(y1),ierr)
  call BINARY_FILE_READ(iunit,z1,nz1+1,kind(z1),ierr)
  call BINARY_FILE_READ(iunit,mask1,nx1*ny1,kind(mask1),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! Read the destination config file
  call BINARY_FILE_OPEN(iunit,trim(fconfig2),"r",ierr)
  call BINARY_FILE_READ(iunit,config2,str_medium,kind(config2),ierr)
  call BINARY_FILE_READ(iunit,icyl2,1,kind(icyl2),ierr)
  call BINARY_FILE_READ(iunit,xper2,1,kind(xper2),ierr)
  call BINARY_FILE_READ(iunit,yper2,1,kind(yper2),ierr)
  call BINARY_FILE_READ(iunit,zper2,1,kind(zper2),ierr)
  call BINARY_FILE_READ(iunit,nx2,1,kind(nx2),ierr)
  call BINARY_FILE_READ(iunit,ny2,1,kind(ny2),ierr)
  call BINARY_FILE_READ(iunit,nz2,1,kind(nz2),ierr)
  print*,'Destination file'
  print*,'Config : ',config2
  print*,'Grid :',nx2,'x',ny2,'x',nz2
  allocate(x2(nx2+1),y2(ny2+1),z2(nz2+1))
  allocate(mask2(nx2,ny2))
  call BINARY_FILE_READ(iunit,x2,nx2+1,kind(x2),ierr)
  call BINARY_FILE_READ(iunit,y2,ny2+1,kind(y2),ierr)
  call BINARY_FILE_READ(iunit,z2,nz2+1,kind(z2),ierr)
  call BINARY_FILE_READ(iunit,mask2,nx2*ny2,kind(mask2),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! Finish to create the meshes
  allocate(xm1(nx1),ym1(ny1),zm1(nz1))
  do i=1,nx1
     xm1(i) = 0.5_WP*(x1(i)+x1(i+1))
  end do
  do j=1,ny1
     ym1(j) = 0.5_WP*(y1(j)+y1(j+1))
  end do
  do k=1,nz1
     zm1(k) = 0.5_WP*(z1(k)+z1(k+1))
  end do
  allocate(xm2(nx2),ym2(ny2),zm2(nz2))
  do i=1,nx2
     xm2(i) = 0.5_WP*(x2(i)+x2(i+1))
  end do
  do j=1,ny2
     ym2(j) = 0.5_WP*(y2(j)+y2(j+1))
  end do
  do k=1,nz2
     zm2(k) = 0.5_WP*(z2(k)+z2(k+1))
  end do
  
  ! Read the source data file
  call BINARY_FILE_OPEN(iunit1,trim(fdata1),"r",ierr)
  call BINARY_FILE_READ(iunit1,nx1, 1,kind(nx1), ierr)
  call BINARY_FILE_READ(iunit1,ny1, 1,kind(ny1), ierr)
  call BINARY_FILE_READ(iunit1,nz1, 1,kind(nz1), ierr)
  call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
  call BINARY_FILE_READ(iunit1,dt,  1,kind(dt),  ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables in source data file: ',names
  
  ! Dump the new data file
  call BINARY_FILE_OPEN(iunit2,trim(fdata2),"w",ierr)
  call BINARY_FILE_WRITE(iunit2,nx2, 1,kind(nx2), ierr)
  call BINARY_FILE_WRITE(iunit2,ny2, 1,kind(ny2), ierr)
  call BINARY_FILE_WRITE(iunit2,nz2, 1,kind(nz2), ierr)
  call BINARY_FILE_WRITE(iunit2,nvar,1,kind(nvar),ierr)
  call BINARY_FILE_WRITE(iunit2,dt,  1,kind(dt),  ierr)
  call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
  do var=1,nvar
     call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
  end do
  
  ! Prepare for interpolation
  allocate(data1(nx1,ny1,nz1))
  allocate(data2(nx2,ny2,nz2))
  
  ! Interpolate
  do var=1,nvar
     ! Read variable
     call BINARY_FILE_READ(iunit1,data1,nx1*ny1*nz1,kind(data1),ierr)
     ! Perform interpolation
     select case (trim(adjustl(names(var))))
     case ('U')
        call interp_data(data1,data2,'U')
     case ('rhoU')
        call interp_data(data1,data2,'U')
     case ('V')
        call interp_data(data1,data2,'V')
     case ('rhoV')
        call interp_data(data1,data2,'V')
     case ('W')
        call interp_data(data1,data2,'W')
     case ('rhoW')
        call interp_data(data1,data2,'W')
     case default
        call interp_data(data1,data2,'SC')
     end select
     ! Write variable
     call BINARY_FILE_WRITE(iunit2,data2,nx2*ny2*nz2,kind(data2),ierr)
     ! Dump status
     print*,trim(adjustl(names(var))),' done'
  end do
  
  ! Close files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  
end program interpolator


! --------------------------------------- !
! Interpolate the data between two arrays !
! For a given staggering 'dir'            !
! --------------------------------------- !
subroutine interp_data(A1,A2,dir)
  use interpolatedata
  implicit none
  
  real(WP), dimension(nx1,ny1,nz1), intent(in)  :: A1
  real(WP), dimension(nx2,ny2,nz2), intent(out) :: A2
  character(len=*) :: dir

  real(WP), dimension(:), pointer :: x,y,z
  integer  :: nx,ny,nz
  integer  :: i2,j2,k2
  integer  :: i1,j1,k1
  integer  :: ni,nj,nk
  real(WP) :: xd,yd,zd
  
  real(WP), dimension(0:1,0:1,0:1) :: ww
  real(WP) :: wx1,wx2,wy1,wy2,wz1,wz2

  ! dir = U  => x  - ym - zm
  ! dir = V  => xm - y  - zm
  ! dir = W  => xm - ym - z
  ! dir = SC => xm - ym - zm
  
  !$OMP PARALLEL DO PRIVATE(k2,j2,i2,i1,j1,k1,ni,nj,nk,xd,yd,zd,x,y,z,nx,ny,nz,wx1,wy1,wz1,wx2,wy2,wz2,ww)
  do k2=1,nz2
     do j2=1,ny2
        do i2=1,nx2
         ! Find the position for interpolation
           select case (trim(adjustl(dir)))
           case ('U')
              xd = x2 (i2)
              yd = ym2(j2)
              zd = zm2(k2)
              x => x1;  nx = nx1+1
              y => ym1; ny = ny1
              z => zm1; nz = nz1
           case ('V')
              xd = xm2(i2)
              yd = y2 (j2)
              zd = zm2(k2)
              x => xm1; nx = nx1
              y => y1;  ny = ny1+1
              z => zm1; nz = nz1
           case ('W')
              xd = xm2(i2)
              yd = ym2(j2)
              zd = z2 (k2)
              x => xm1; nx = nx1
              y => ym1; ny = ny1
              z => z1;  nz = nz1+1
           case ('SC')
              xd = xm2(i2)
              yd = ym2(j2)
              zd = zm2(k2)
              x => xm1; nx = nx1
              y => ym1; ny = ny1
              z => zm1; nz = nz1
           end select
           
           ! Find the nearest points in mesh1
           call bisection(xd,i1,x,nx)
           call bisection(yd,j1,y,ny)
           call bisection(zd,k1,z,nz)
           i1 = max(1,min(i1,nx1-1))
           j1 = max(1,min(j1,ny1-1))
           k1 = max(1,min(k1,nz1-1))
           
           ! Interpolate the point
           wx1 = 1.0_WP
           wy1 = 1.0_WP
           wz1 = 1.0_WP
           if (nx.ne.1) wx1 = (x(i1+1)-xd   )/(x(i1+1)-x(i1))
           if (ny.ne.1) wy1 = (y(j1+1)-yd   )/(y(j1+1)-y(j1))
           if (nz.ne.1) wz1 = (z(k1+1)-zd   )/(z(k1+1)-z(k1))
           
           ! Points outside the domain
           if (xd.lt.x(1)) then
              i1  = 1
              wx1 = 1.0_WP
           end if
           if (xd.gt.x(nx1)) then
              i1  = nx1-1
              wx1 = 0.0_WP
           end if
           if (yd.lt.y(1)) then
              j1  = 1
              wy1 = 1.0_WP
           end if
           if (yd.gt.y(ny1)) then
              j1  = ny1-1
              wy1 = 0.0_WP
           end if
           if (zd.lt.z(1)) then
              k1  = 1
              wz1 = 1.0_WP
           end if
           if (zd.gt.z(nz1)) then
              k1  = nz1-1
              wz1 = 0.0_WP
           end if

           ! Masks
           ! JFM 8/18/17 - This can be better for scalars, but fixes mass conservation issues
!!$           if (mask1(i1,j1)  .gt.0) then
!!$              wx1 = 0.0_WP
!!$              wy1 = 0.0_WP
!!$              wz1 = 0.0_WP
!!$           end if
           if (trim(adjustl(dir)).eq.'SC') then
              if (mask1(i1,j1).gt.0.and.mask1(i1,j1-1).gt.0) wy1 = 0.0_WP
              if (mask1(i1,j1).gt.0.and.mask1(i1-1,j1).gt.0) wx1 = 0.0_WP
              if (mask1(i1+1,j1).gt.0) wx1 = 1.0_WP
              if (mask1(i1,j1+1).gt.0) wy1 = 1.0_WP
           end if
           
           ! Compute the other coefficients
           wx2 = 1.0_WP-wx1
           wy2 = 1.0_WP-wy1
           wz2 = 1.0_WP-wz1
           
           ! Combine the interpolation coefficients to form a tri-linear interpolation
           ww(0,0,0)=wx1*wy1*wz1
           ww(1,0,0)=wx2*wy1*wz1
           ww(0,1,0)=wx1*wy2*wz1
           ww(1,1,0)=wx2*wy2*wz1
           ww(0,0,1)=wx1*wy1*wz2
           ww(1,0,1)=wx2*wy1*wz2
           ww(0,1,1)=wx1*wy2*wz2
           ww(1,1,1)=wx2*wy2*wz2
           
           ! Perform the actual interpolation on A
           A2(i2,j2,k2) = 0.0_WP
           do nk=0,min(k1+1,nz1)-k1
              do nj=0,min(j1+1,ny1)-j1
                 do ni=0,min(i1+1,nx1)-i1
                    A2(i2,j2,k2) = A2(i2,j2,k2) + ww(ni,nj,nk)*A1(i1+ni,j1+nj,k1+nk)
                 end do
              end do
           end do

        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  return
end subroutine interp_data



! -------------------------------------------------- !
! Bisection routine                                  !
! Gets an array and its size as well as a position x !
! Returns the index between 1 and nx                 !
! x(iloc)<=xloc<x(iloc+1)                            !
! Assuming xarray is monotonically increasing        !
! -------------------------------------------------- !
subroutine bisection(xloc,iloc,x,nx)
  use precision
  implicit none
  
  real(WP), intent(in)  :: xloc
  integer,  intent(out) :: iloc
  integer,  intent(in)  :: nx
  real(WP), dimension(1:nx), intent(in) :: x
  
  integer :: il,im,iu
  
  ! Take care of outside points
  if (xloc.lt.x(1)) then
     iloc = 1
  else if (xloc.ge.x(nx)) then
     iloc = nx-1
  else
     
     ! Initialize lower and upper limits
     il=1
     iu=nx
     
     ! While not done
     do while (iu-il.gt.1)
        ! Compute a mid-point
        im=(iu+il)/2
        ! Replace lower of upper limit as appropriate
        if (xloc.ge.x(im)) then
           il=im
        else
           iu=im
        end if
     end do
     
     ! Return
     iloc = il
  end if
  
  return
end subroutine bisection
