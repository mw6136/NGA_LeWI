program flreg
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz
  integer :: xper,yper,zper
  integer :: icyl
  real(WP), dimension(:), pointer :: x,y,z
  integer, dimension(:,:), pointer :: mask
  !integer, dimension(:,:,:), pointer :: iblank
  integer :: iunit,ierr
  character(len=str_medium) :: filename1,filename2,filename3
  character(len=str_medium) :: config
  integer :: i
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:,:), pointer :: data
  real(WP), dimension(:,:,:), pointer :: C
  integer :: var,nvar
  real(WP) :: L
  real(WP) :: dt,time
  real(WP) :: epsilon

  ! Read file name from standard input
  print*,'=================================='
  print*,'| ARTS - Flame Regimes evaluator |'
  print*,'=================================='
  print*
  print "(a,$)", " config file : "
  read "(a)", filename1
  print "(a,$)", " data file : "
  read "(a)", filename2
  print "(a,$)", " spectrum file : "
  read "(a)", filename3

  ! ** Open the config file to read **
  call BINARY_FILE_OPEN(iunit,trim(filename1),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_READ(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_READ(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_READ(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)
  print*,'Config : ',config
  print*,'Grid :',nx,'x',ny,'x',nz

  ! Read grid field
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(mask(nx,ny))
  call BINARY_FILE_READ(iunit,x,nx+1,kind(x),ierr)
  call BINARY_FILE_READ(iunit,y,ny+1,kind(y),ierr)
  call BINARY_FILE_READ(iunit,z,nz+1,kind(z),ierr)
  call BINARY_FILE_READ(iunit,mask,nx*ny,kind(mask),ierr)

  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Get domain size
  L = x(nx+1) - x(1)
  print*,'Domain size:',L

  ! ** Open the data file to read **
  call BINARY_FILE_OPEN(iunit,trim(filename2),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit,nvar,1,kind(nvar),ierr)

  ! Read additional stuff
  call BINARY_FILE_READ(iunit,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit,time,1,kind(time),ierr)
  print*,'Data file at time :',time

  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit,names(var),str_short,kind(names),ierr)
  end do
  print*,'Data name :'
  do var=1,nvar
     print*, var, ' - ',names(var)
  end do

  ! Read data field
  allocate(data(nx,ny,nz,nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit,data(:,:,:,var),nx*ny*nz,kind(data),ierr)
  end do
  call BINARY_FILE_CLOSE(iunit,ierr)

  !Progress variable C=Y_H2O+Y_H2+Y_CO+YCO2 
  allocate(data(nx,ny,nz))
  C(:,:,:) = data(:,:,:,12)+data(:,:,:,13)+data(:,:,:,16)+data(:,:,:,17)



  ! Output it
  iunit=iopen()
  open(iunit,file=filename3,form='formatted')
  !write(11,*) 
  close(iclose(iunit))

  ! Deallocation
  deallocate(x,y,z)
  deallocate(mask)
  deallocate(data)

end program flreg
