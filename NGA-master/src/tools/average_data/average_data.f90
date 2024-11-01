! ============================================================ !
!                      average_data.f90                        !
! Program to perform an average across periodic directions     !
!    and time. Requires NGA config file and data files.        !
!                                                              !
! Author: Jonathan F. MacArt                                   !
! Date:   October 10, 2017                                     !
! ============================================================ !
module average_data
  use parallel
  use precision
  use string
  implicit none
  
  ! Sizes
  integer :: nx,ny,nz
  integer :: nx2,ny2,nz2
  integer :: xper,yper,zper
  integer :: icyl
  integer :: nvar
  real(WP) :: dt,dt_file,time,time_prev
  character(len=str_short), dimension(:), pointer :: names
  character(len=str_medium) :: config
  integer, dimension(:,:), pointer :: mask
  
  ! Data
  real(WP), dimension(:,:,:), pointer :: data
  real(WP), dimension(:,:,:), pointer :: data_avg
  
  ! Mesh
  real(WP), dimension(:), pointer :: x,x2
  real(WP), dimension(:), pointer :: y,y2
  real(WP), dimension(:), pointer :: z,z2

end module average_data


! ========================================================== !
! MAIN ROUTINE                                               !
! ========================================================== !
program average_data_main
  use average_data
  use parser
  use fileio
  use cli_reader
  implicit none
  
  character(len=str_medium) :: input_name
  character(len=str_medium) :: fconfig, fconfig2
  character(len=str_medium) :: data_dir
  character(len=str_medium) :: fdata, output_file
  character(len=str_medium), dimension(:), pointer :: data_files
  integer :: Nfiles,Nfiles_max,ifile
  integer :: iunit,iunit2,var
  integer :: i,j,k,ierr,ii,jj,i2
  real(WP) :: tmp
  real(WP), dimension(:,:), pointer :: tmpxy

!!$  ! Initialize the parallel module
!!$  call parallel_init
!!$  call parallel_init_filesystem

  ! Parse the input file
!!$  if (irank.eq.iroot) then
     call get_command_argument(1,input_name)
!!$  end if
!!$  call MPI_BCAST(input_name, str_medium, MPI_CHARACTER, iroot-1, MPI_COMM_WORLD, ierr)
  call parser_init
  call parser_parsefile(input_name)

  ! Read the source config file
  call parser_read('Data root dir',data_dir,'.')
  call parser_read('Configuration file',fconfig)
  fconfig2 = trim(data_dir)//'/'//trim(fconfig)
  call BINARY_FILE_OPEN(iunit,trim(fconfig2),"r",ierr)
  call BINARY_FILE_READ(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_READ(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_READ(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_READ(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(mask(nx,ny))
  call BINARY_FILE_READ(iunit,x,nx+1,kind(x),ierr)
  call BINARY_FILE_READ(iunit,y,ny+1,kind(y),ierr)
  call BINARY_FILE_READ(iunit,z,nz+1,kind(z),ierr)
  call BINARY_FILE_READ(iunit,mask,nx*ny,kind(mask),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Account for periodicity
  nx2 = nx
  ny2 = ny
  nz2 = nz
  if (xper.eq.1) nx2 = 1
  if (yper.eq.1) ny2 = 1
  if (zper.eq.1) nz2 = 1
  allocate(x2(nx2+1))
  allocate(y2(ny2+1))
  allocate(z2(nz2+1))
  x2 = x
  y2 = y
  z2 = z
  if (xper.eq.1) x2 = 0.0_WP
  if (yper.eq.1) y2 = 0.0_WP
  if (zper.eq.1) then
     do k=1,nz2+1
        z2(k) = real(k-1,WP)*(z(nz+1)-z(1))/real(nz2,WP)
     end do
  end if

  ! Write the config file for the average data
  fconfig2 = trim(fconfig)//'_avg'
  call BINARY_FILE_OPEN(iunit,trim(fconfig2),"w",ierr)
  call BINARY_FILE_WRITE(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_WRITE(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_WRITE(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_WRITE(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_WRITE(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_WRITE(iunit,nx2,1,kind(nx2),ierr)
  call BINARY_FILE_WRITE(iunit,ny2,1,kind(ny2),ierr)
  call BINARY_FILE_WRITE(iunit,nz2,1,kind(nz2),ierr)
  call BINARY_FILE_WRITE(iunit,x2,nx2+1,kind(x2),ierr)
  call BINARY_FILE_WRITE(iunit,y2,ny2+1,kind(y2),ierr)
  call BINARY_FILE_WRITE(iunit,z2,nz2+1,kind(z2),ierr)
  call BINARY_FILE_WRITE(iunit,mask,nx*ny,kind(mask),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Get the list of data files to read
  call parser_getsize('Data files',Nfiles)
  allocate(data_files(Nfiles))
  call parser_read('Data files',data_files)
  call parser_read('Max Nfiles',Nfiles_max,Nfiles)
  call parser_read('Output file',output_file)
  if (Nfiles_max.gt.Nfiles) Nfiles_max = Nfiles

  ! Allocate data file to read and space for averaged data
  allocate(data(nx,ny,nz))
  if (zper.eq.1) then
     allocate(tmpxy(nx,ny))
  else
     print *, 'not implemented'
  end if
     
  data_avg = 0.0_WP

  ! Read the files and compute statistics
  time_prev = 0.0_WP
  do ifile=1,Nfiles_max
     fdata = trim(data_dir)//'/'//trim(data_files(ifile))
     call BINARY_FILE_OPEN(iunit,trim(fdata),"r",ierr)
     call BINARY_FILE_READ(iunit,nx,  1,kind(nx),  ierr)
     call BINARY_FILE_READ(iunit,ny,  1,kind(ny),  ierr)
     call BINARY_FILE_READ(iunit,nz,  1,kind(nz),  ierr)
     call BINARY_FILE_READ(iunit,nvar,1,kind(nvar),ierr)
     call BINARY_FILE_READ(iunit,dt,  1,kind(dt),  ierr)
     call BINARY_FILE_READ(iunit,time,1,kind(time),ierr)
     if (ifile.eq.1) then
        allocate(names(nvar))
     end if
     if (ifile.eq.1.and.zper.eq.1) then
        allocate(data_avg(nx2,ny2,nvar))
     end if
     do var=1,nvar
        call BINARY_FILE_READ(iunit,names(var),str_short,kind(names),ierr)
     end do

     dt_file = time-time_prev
     print '(a8,i3,a8,ES12.5,a8,ES12.5)','  ifile=',ifile,'   time=',time,'     dt=',dt_file

     do var=1,nvar
        ! Read the data
        call BINARY_FILE_READ(iunit,data,nx*ny*nz,kind(data),ierr)
        
        ! Update the temporal average
        if (zper.eq.1) then
           tmpxy = 0.0_WP
           do k=1,nz
              do j=1,ny
                 do i=1,nx
                    tmpxy(i,j) = tmpxy(i,j) + data(i,j,k)
                 end do
              end do
           end do
           tmpxy = tmpxy/real(nz,WP)

           do j=1,ny
              do i=1,nx
                 data_avg(i,j,var) = (real(ifile-1,WP)*data_avg(i,j,var) + tmpxy(i,j))/real(ifile,WP)
              end do
           end do
        else
           print *, 'not implemented'
        end if
   
        print *, trim(names(var)),' done'
     end do
     call BINARY_FILE_CLOSE(iunit,ierr)

     ! Save data for next time step
     time_prev = time
  end do

  ! Write the averaged data file
  call BINARY_FILE_OPEN(iunit,trim(output_file),"w",ierr)
  call BINARY_FILE_WRITE(iunit,nx2,  1,kind(nx2),  ierr)
  call BINARY_FILE_WRITE(iunit,ny2,  1,kind(ny2),  ierr)
  call BINARY_FILE_WRITE(iunit,nz2,  1,kind(nz2),  ierr)
  call BINARY_FILE_WRITE(iunit,nvar,1,kind(nvar),ierr)
  call BINARY_FILE_WRITE(iunit,dt,  1,kind(dt),  ierr)
  call BINARY_FILE_WRITE(iunit,time,1,kind(time),ierr)
  do var=1,nvar
     call BINARY_FILE_WRITE(iunit,names(var),str_short,kind(names),ierr)
  end do
  do var=1,nvar
     call BINARY_FILE_WRITE(iunit,data_avg(:,:,var),nx2*ny2,kind(data_avg),ierr)
  end do


!!$  ! Initialize the MPI module and read data in parallel
!!$  call average_data_MPI_init
  ! Use parallel_sum
  ! Could use parallel_gather_dir to communicate along lines  

end program average_data_main
