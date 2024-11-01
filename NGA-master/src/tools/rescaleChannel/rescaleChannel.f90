! -------------------------------------------------------------------------- !
!                          rescaleChannel.F90                                !
!     Purpose     : Rescales a channel data file, keeping z-direction        !
!     Date        : October 25, 2016                                         !
!     Written by  : Jonathan F. MacArt                                       !
! -------------------------------------------------------------------------- !


program rescaleChannel
  use precision
  use string
  use fileio
  use parser
  use cli_reader
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  character(len=str_medium) :: input_name
  character(len=str_long) :: fconfig1,fconfig2,fdata1,fdata2
  character(len=str_short), dimension(:), pointer :: names
  character(len=str_medium) :: config
  integer, dimension(:,:), pointer :: mask

  ! Sizes
  integer :: nx,ny,nz1,nz2
  integer :: icyl
  integer :: xper,yper,zper

  ! Mesh
  real(WP), dimension(:), pointer :: x1
  real(WP), dimension(:), pointer :: y1
  real(WP), dimension(:), pointer :: z1
  real(WP), dimension(:), pointer :: x2
  real(WP), dimension(:), pointer :: y2
  real(WP), dimension(:), pointer :: z2
  integer :: nvar,k,kk
  integer :: iunit,iunit1,iunit2,ierr,var
  real(WP) :: dt,time
  real(WP) :: h_scale,v_scale,vert_offs
  real(WP) :: Lx1,Lx2,Ly1,Ly2,Lz,Lz1
  real(WP) :: dz1,dz2,dz_tmp,zloc1,zloc2,wz1,wz2
  integer :: nrep

  ! Data
  real(WP), dimension(:,:,:), pointer :: data

  ! Parse the input file
  call get_command_argument(1,input_name)
  call parser_init
  call parser_parsefile(input_name)

  ! Get the file names
  call parser_read('Source config file',     fconfig1)
  call parser_read('Destination config file',fconfig2)
  call parser_read('Source data file',     fdata1)
  call parser_read('Destination data file',fdata2)

  ! Read the source config file
  call BINARY_FILE_OPEN(iunit,trim(fconfig1),"r",ierr)
  call BINARY_FILE_READ(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_READ(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_READ(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_READ(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz1,1,kind(nz1),ierr)
  print*,'Source file'
  print*,'Config : ',config
  allocate(x1(nx+1),y1(ny+1),z1(nz1+1))
  allocate(mask(nx,ny))
  call BINARY_FILE_READ(iunit,x1,nx+1,kind(x1),ierr)
  call BINARY_FILE_READ(iunit,y1,ny+1,kind(y1),ierr)
  call BINARY_FILE_READ(iunit,z1,nz1+1,kind(z1),ierr)
  call BINARY_FILE_READ(iunit,mask,nx*ny,kind(mask),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Compute the scaling
  call parser_read('Scaling distance',h_scale)
  call parser_read('Scaling velocity',v_scale)
  call parser_read('Destination nz',nz2,nz1)
  call parser_read('Vertical offset',vert_offs,0.0_WP)
  allocate(x2(nx+1),y2(ny+1),z2(nz2+1))
  x2 = x1/h_scale
  y2 = y1-vert_offs
  y2 = y2/h_scale+vert_offs
  Lx1 = x1(nx+1)-x1(1)
  Lx2 = x2(nx+1)-x2(1)
  Ly1 = y1(ny+1)-y1(1)
  Ly2 = y2(ny+1)-y2(1)  
  Lz  = z1(nz1+1)-z1(1)
  dz1 = Lz/real(nz1,WP)/h_scale
  dz2 = Lz/real(nz2,WP)
  do k=1,nz2+1
     z2(k) = real(k-1,WP)/real(nz2)*Lz - 0.5_WP*Lz
  end do
  print*,'Source grid :',nx,'x',ny,'x',nz1
  print*,'Dest.  grid :',nx,'x',ny,'x',nz2
  print '(a,3ES13.5)', 'Source dims :', Lx1, Ly1, Lz
  print '(a,3ES13.5)', 'Dest.  dims :', Lx2, Ly2, Lz
  
  ! Write the destination config file
  call BINARY_FILE_OPEN(iunit,trim(fconfig2),"w",ierr)
  call BINARY_FILE_WRITE(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_WRITE(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_WRITE(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_WRITE(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_WRITE(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit,nz2,1,kind(nz2),ierr)
  call BINARY_FILE_WRITE(iunit,x2,nx+1,kind(x2),ierr)
  call BINARY_FILE_WRITE(iunit,y2,ny+1,kind(y2),ierr)
  call BINARY_FILE_WRITE(iunit,z2,nz2+1,kind(z2),ierr)
  call BINARY_FILE_WRITE(iunit,mask,nx*ny,kind(mask),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! Read the source data file
  call BINARY_FILE_OPEN(iunit1,trim(fdata1),"r",ierr)
  call BINARY_FILE_READ(iunit1,nx, 1,kind(nx), ierr)
  call BINARY_FILE_READ(iunit1,ny, 1,kind(ny), ierr)
  call BINARY_FILE_READ(iunit1,nz1, 1,kind(nz1), ierr)
  call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
  call BINARY_FILE_READ(iunit1,dt,  1,kind(dt),  ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables in source data file: ',names
  
  ! Start the new data file
  call BINARY_FILE_OPEN(iunit2,trim(fdata2),"w",ierr)
  call BINARY_FILE_WRITE(iunit2,nx, 1,kind(nx), ierr)
  call BINARY_FILE_WRITE(iunit2,ny, 1,kind(ny), ierr)
  call BINARY_FILE_WRITE(iunit2,nz2, 1,kind(nz2), ierr)
  call BINARY_FILE_WRITE(iunit2,nvar,1,kind(nvar),ierr)
  call BINARY_FILE_WRITE(iunit2,dt,  1,kind(dt),  ierr)
  call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
  do var=1,nvar
     call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
  end do

  ! Allocate arrays
  allocate(data(nx,ny,nz1))

  ! Shift z to make things easier
  Lz1 = Lz/h_scale
  z1 = z1 + 0.5_WP*Lz
  z1 = z1/h_scale
  z2 = z2 + 0.5_WP*Lz

  ! Scale the channel's length, height, and velocity, matching the width
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,data,nx*ny*nz1,kind(data),ierr)
     ! Scale velocity
     if (trim(names(var)).eq.'U'.or.trim(names(var)).eq.'V'.or.trim(names(var)).eq.'W') then
        data = data/v_scale
        print *, 'Scaled ', trim(names(var))
     end if
     ! If h_scale > 1 : H decreases; need more data in z -- append
     ! If h_scale < 1 : H increases; need less data in z -- interpolate
     print '(2a5,6a13)', 'k','kk','z1','z1+dz1','corr z2','z2','wz1','wz2'
     k  = 1
     nrep = 0
     do kk=1,nz2
        ! If necessary, increment the source file
        do while (z2(kk)-real(nrep,WP)*(Lz1-dz1).ge.z1(k+1)) 
           k = k + 1
        end do
        ! If necessary, start from the beginning of the source file
        if (k+1.gt.nz1) then
           k = 1
           nrep = nrep + 1
        end if

        ! Interpolate in z
        wz2 = (z2(kk)-real(nrep,WP)*(Lz1-dz1)-z1(k))/dz1
        wz1 = 1.0_WP-wz2
        print '(2i5,6ES13.5)', k,kk,z1(k),z1(k+1),z2(kk)-real(nrep,WP)*(Lz1-dz1),z2(kk),wz1,wz2
        call BINARY_FILE_WRITE(iunit2,(wz1*data(:,:,k)+wz2*data(:,:,k+1)),nx*ny,kind(data),ierr)
     end do
  end do
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  
end program rescaleChannel
