program editInflow
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: ntime,ny,nz,nvar,newvar,n_del,icyl
  character(len=str_short), dimension(:), pointer :: names,names2
  character(len=str_short), dimension(10) :: names_del
  real(WP), dimension(:,:), pointer :: data,avg_data,line_data
  real(WP), dimension(:), pointer :: y,z
  integer :: iunit1,iunit2,ierr,var,iunit3,var2,itime,i,j,k,choice
  character(len=str_medium) :: filename1,filename2
  character(len=str_short) :: varname,varname2
  real(WP) :: dt,time,value

  ! Read file name from standard input
  print*,'========================'
  print*,'| ARTS - inflow Editor |'
  print*,'========================'
  print*
  print "(a20,$)", " inflow file : "
  read "(a)", filename1
  print "(a17,$)", " file to print: "
  read "(a)", filename2

  ! ** Open the data file to read **
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit1,ntime,1,kind(ntime),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
  print*,'Grid :',ntime,'ntime',ny,'ny',nz,'nz'
 
  ! Read additional stuff
  call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  print*,'Inflow file at time :',time
  
  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables : ',names
  print*,'There are ',nvar,' variables.'

  ! Allocate arrays, read spatial coordinates
  allocate(data(ny,nz))
  allocate(y(ny+1))
  allocate(z(nz+1))
  call BINARY_FILE_READ(iunit1,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit1,y,ny+1,kind(y),ierr)
  call BINARY_FILE_READ(iunit1,z,nz+1,kind(z),ierr)

 ! ** Ask what to do **
  print*
  print*, "1. Take average"
  print*, "2. Show a point through time"
  print "(a9,$)", "Choice : "
  read "(i1)", choice
  
  ! Case dependent operation
  select case(choice)
  case(1) ! Average inflow quantities for entire file
     ! Give info
     print*
     print*, "Averaging inflow quantities..."

     ! Averaging operations
     ! Allocate arrays for averaged quantities
     newvar=3*nvar
     allocate(avg_data(ny,newvar))
     avg_data = 0

     ! Gather and sum the data
     do itime=1,ntime
        do var=1,nvar
           call BINARY_FILE_READ(iunit1,data(:,:),ny*nz,kind(data),ierr)
           do k=1,nz
              do j=1,ny
                 avg_data(j,var) = avg_data(j,var) + data(j,k)
              end do
           end do
           ! Squares
           do k=1,nz
              do j=1,ny
                 avg_data(j,var+nvar) = avg_data(j,var+nvar) + data(j,k)**2
              end do
           end do
           ! Cubes
           do k=1,nz
              do j=1,ny
                 avg_data(j,var+nvar+nvar) = avg_data(j,var+nvar+nvar) + data(j,k)**3
              end do
           end do
        end do
     end do

     ! Close the file
     call BINARY_FILE_CLOSE(iunit1,ierr)

     ! Actually take the averages
     do var=1,newvar
        do j=1,ny
           avg_data(j,var) = avg_data(j,var)/(nz*ntime)
        end do
     end do

     ! Print to file
     open (iunit2, file=filename2)
     write(iunit2,'(10000a20)') 'y', (names(var), var=1,nvar), (names(var), var=1,nvar),"^2", (names(var), var=1,nvar),"^3"
     do j=1,ny
        write(iunit2,'(10000ES20.12)') y(j), (avg_data(j,var) , var=1,newvar)
     end do
     close(iclose(iunit2))

  case (2) ! For a single point, track data in time
     allocate(line_data(ntime,nvar))
     line_data = 0

     ! Gather the data
     do itime=1,ntime
        do var=1,nvar
           call BINARY_FILE_READ(iunit1,data(:,:),ny*nz,kind(data),ierr)
           line_data(itime,var) = data(ny,10)
        end do
     end do

     ! Close the file
     call BINARY_FILE_CLOSE(iunit1,ierr)

     ! Print to file
     open (iunit2, file=filename2)
     write(iunit2,'(10000a20)') 't/T', (names(var), var=1,nvar)
     do itime=1,ntime
        write(iunit2,'(10000ES20.12)') itime, (line_data(itime,var) , var=1,nvar)
     end do
     close(iclose(iunit2))
  end select

end program editInflow
