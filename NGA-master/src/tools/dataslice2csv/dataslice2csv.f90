program dataslice2csv
  use precision
  use string
  use fileio
  use cli_reader
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz,nvar, i, j, kk
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:), pointer :: data
  integer :: iunit1,iunit2,ierr,var
  character(len=str_medium) :: filename1,filename2
  character(len=str_short) :: varname
  real(WP) :: dt,time,value

  ! Read file name from the command line
  call get_command_argument(1,filename1)

  ! ** Open the data file to read **
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit1,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
  print*,'Grid :',nx,'x',ny,'x',nz
    
  ! Read additional stuff
  call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  print*,'Timestep size     :',dt
  print*,'Data file at time :',time
  
  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables : ',names
  print*,'There are ',nvar,' variables.'

  ! Open the output CSV file and print the header
  ! Order: nx, ny
  open(unit=iunit2, file='zslice_'//trim(filename1)//'.csv', action='write')
  write (iunit2,'(I10)') nx
  write (iunit2,'(I10)') ny


  ! Allocate arrays
  allocate(data(nx,ny,nz))
  kk = nz/2
  do var=1,nvar
     ! Read the data
     call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)

     ! Write the slice at z=0 to the CSV file
     do j=1,ny-1
        do i=1,nx
           write (iunit2,'(ES22.13)',advance='no'), data(i,j,kk)
        end do
     end do
     do i=1,nx-1
        write (iunit2,'(ES22.13)',advance='no'), data(i,ny,kk)
     end do
     write (iunit2,'(ES22.13)'), data(nx,ny,kk)
  end do
     
  ! Close the files
  close(iunit2)
  call BINARY_FILE_CLOSE(iunit1,ierr)
  
end program dataslice2csv
