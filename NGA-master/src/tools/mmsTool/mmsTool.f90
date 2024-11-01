program mmsTool
  use precision
  use string
  use fileio
  use cli_reader
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz,nvar,i,j,k,ii
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:), pointer :: x,y,z
  real(WP), dimension(:,:,:), pointer :: data
  integer :: iunit1,iunit2,ierr,var,icyl,xper,yper,zper
  character(len=str_medium) :: topfolder,filename,simu_name
  character(len=str_short) :: varname,varname2
  character(len=str_short), dimension(7) :: gridfolders
  real(WP) :: dt,time,value,xm
  real(WP) :: Lx, tau, DIFF, RHO, U, Z_ex, err_L2
  real(WP), parameter :: pi = 4.0_WP*atan(1.0_WP)

  ! Read file name from standard input
  call get_command_argument(1,topfolder)
  !print*,'=============================='
  !print*,'| ARTS - mms comparison tool |'
  !print*,'=============================='
  !print*
  !print "(a15,$)", " folder name : "
  !read "(a)", topfolder

  ! Simulation parameters
  Lx   = 0.025_WP
  tau  = 4.0e-5_WP
  gridfolders(1) = '64'
  gridfolders(2) = '128'
  gridfolders(3) = '256'
  gridfolders(4) = '512'
  gridfolders(5) = '1024'
  gridfolders(6) = '2048'
  gridfolders(7) = '4096'

  do ii=1,7
     ! ** Open the config file **
     filename = trim(topfolder)//'/'//trim(gridfolders(ii))//'/'//'config-n'//trim(gridfolders(ii))
     call BINARY_FILE_OPEN(iunit1,trim(filename),"r",ierr)
     call BINARY_FILE_READ(iunit1,simu_name,str_medium,kind(simu_name),ierr)
     call BINARY_FILE_READ(iunit1,icyl,1,kind(icyl),ierr)
     call BINARY_FILE_READ(iunit1,xper,1,kind(xper),ierr)
     call BINARY_FILE_READ(iunit1,yper,1,kind(yper),ierr)
     call BINARY_FILE_READ(iunit1,zper,1,kind(zper),ierr)
     call BINARY_FILE_READ(iunit1,nx,1,kind(nx),ierr)
     call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
     call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
     allocate(x(1:nx+1))
     allocate(y(1:ny+1))
     allocate(z(1:nz+1))
     call BINARY_FILE_READ(iunit1,x,nx+1,kind(x),ierr)
     call BINARY_FILE_READ(iunit1,y,ny+1,kind(y),ierr)
     call BINARY_FILE_READ(iunit1,z,nz+1,kind(z),ierr)
     call BINARY_FILE_CLOSE(iunit1)

     ! ** Open the data file to read **
     filename = trim(topfolder)//'/'//trim(gridfolders(ii))//'/'//'data.1'
     !print *, filename
     call BINARY_FILE_OPEN(iunit1,trim(filename),"r",ierr)
     
     ! Read sizes
     call BINARY_FILE_READ(iunit1,nx,1,kind(nx),ierr)
     call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
     call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
     call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
     !print*,'Grid :',nx,'x',ny,'x',nz
     
     ! Read additional stuff
     call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
     call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
     !print*,'Timestep size     :',dt
     !print*,'Data file at time :',time
     
     ! Read variable names
     allocate(names(nvar))
     do var=1,nvar
        call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
     end do
     !print*,'Variables : ',names
     !print*,'There are ',nvar,' variables.'

     ! Allocate arrays
     allocate(data(nx,ny,nz))
     
     do var=1,nvar
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        !print*,"min: ",minval(data)," - max: ",maxval(data)
        !print*,maxloc(data)
     end do
     ! We now have the mixture fraction field (var=nvar)

     ! Compute L2 norm vs exact solution
     err_L2 = 0.0_WP
     do k=1,nz
        do j=1,ny
           do i=1,nx
              xm     = 0.5_WP*(x(i)+x(i+1))
              Z_ex   = 0.5_WP*sin(2.0_WP*pi*xm/Lx)*sin(2.0_WP*pi*time/tau) + 0.5_WP
              err_L2 = err_L2 + (data(i,j,k)-Z_ex)**2
              !print '(ES22.13,ES22.13,ES22.13)', xm, Z_ex, data(i,j,k)
           end do
        end do
     end do
     err_L2 = sqrt(err_L2/nx)
     print '(A6,ES22.13)', trim(gridfolders(ii)), err_L2
     
     ! Close the file
     call BINARY_FILE_CLOSE(iunit1,ierr)
     deallocate(x)
     deallocate(y)
     deallocate(z)
     deallocate(data)
  end do
  
end program mmsTool
