program editVolumeData
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: ntime,nx,ny,nz,nvar,n_del
  real(WP), dimension(:), pointer :: x,y,z
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:), pointer :: data
  integer :: i,j,k
  integer :: iunit1,iunit2,ierr,var,choice
  character(len=str_medium_2) :: filename1,filename2
  character(len=str_short) :: varname,varname2
  character(len=str_short), dimension(10) :: names_del
  real(WP) :: dt,time,value,data_value

  ! Read file name from standard input
  print*,'============================='
  print*,'| ARTS - volume data Editor |'
  print*,'============================='
  print*
  print "(a28,$)", " data file before editing : "
  read "(a)", filename1
  print "(a27,$)", " data file after editing : "
  read "(a)", filename2

  ! ** Open the data file to read **
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit1,ntime,1,kind(ntime),ierr)
  call BINARY_FILE_READ(iunit1,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
  print*,'Grid :',nx,'x',ny,'x',nz

  ! Read the grid
  allocate(x(nx))
  allocate(y(ny))
  allocate(z(nz))
  call BINARY_FILE_READ(iunit1,x,nx,kind(x),ierr)
  call BINARY_FILE_READ(iunit1,y,ny,kind(y),ierr)
  call BINARY_FILE_READ(iunit1,z,nz,kind(z),ierr)
  
  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables : ',names
  print*,'There are ',nvar,' variables.'
    
  ! Read additional stuff
  call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  print*,'Timestep size     :',dt
  print*,'Data file at time :',time

  ! Allocate arrays
  allocate(data(nx,ny,nz))

  ! ** Ask what to do **
  print*
  print*, "1. Print Min/Max of variable"
  print*, "2. Add variable -- NOT READY"
  print*, "3. Delete variable -- NOT READY"
  print*, "4. Reset variable"
  print "(a9,$)", "Choice : "
  read "(i1)", choice
  
  ! Case dependent operation
  select case(choice)
     
  case(1) ! Print min/max of all variables
     do var=1,nvar
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        print*,trim(names(var))
        print*,"min: ",minval(data)," - max: ",maxval(data)
        print*,maxloc(data)
     end do
     
  case (2) ! Add variable
     !  NEEDS TO BE UPDATED FOR VOLUME FILE FORMAT
     print "(a16,$)", "Variable name : "
     read "(a)", varname
     print "(a16,$)", "Default value : "
     read(*,*) value
     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar+1,1,kind(nvar),ierr)
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
     do var=1,nvar
        call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
     end do
     call BINARY_FILE_WRITE(iunit2,varname,str_short,kind(varname),ierr)
     do var=1,nvar
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do
     data = value
     call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     call BINARY_FILE_CLOSE(iunit2,ierr)
     
  case (3) ! Delete variables 
     !  NEEDS TO BE UPDATED FOR VOLUME FILE FORMAT
     print*,'You can delete at most 10 variables, press q to exit'
     n_del = 0
     do while (n_del .le. 10)
        print "(a16,$)", "Variable name : "
        read "(a)", varname
        if (varname .ne. 'q') then
           n_del = n_del + 1
           names_del(n_del) = varname
        else
           exit   
        end if
     end do

     print*,'how many are deleted? ',n_del
     print*,'They are :',names_del(1:n_del)

     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar-n_del,1,kind(nvar),ierr)
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)

     do var=1,nvar
        if (.not. any(names_del.eq.(trim(adjustl(names(var)))))) &
             call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
     end do
     do var=1,nvar
        call BINARY_FILE_READ(iunit1,data,nx*ny*nz,kind(data),ierr)
        if (.not. any(names_del.eq.(trim(adjustl(names(var)))))) &
             call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
     end do
     call BINARY_FILE_CLOSE(iunit2,ierr)
          
  case (4) ! Reset variable
     print "(a26,$)", " variable name to reset : "
     read "(a)", varname2
     print "(a13,$)", " new value : "
     read(*,*), data_value
     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,ntime,1,kind(ntime),ierr)
     call BINARY_FILE_WRITE(iunit2,nx,1,kind(nx),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar,1,kind(nvar),ierr)

     call BINARY_FILE_WRITE(iunit2,x,nx,kind(x),ierr)
     call BINARY_FILE_WRITE(iunit2,y,ny,kind(y),ierr)
     call BINARY_FILE_WRITE(iunit2,z,nz,kind(z),ierr)
     
     do var=1,nvar
        call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
     end do
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
     
     do var=1,nvar
        if (trim(names(var)).eq.trim(varname2)) then
           call BINARY_FILE_READ (iunit1,data,nx*ny*nz,kind(data),ierr)
           data = data_value
           call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
        else
           call BINARY_FILE_READ (iunit1,data,nx*ny*nz,kind(data),ierr)
           call BINARY_FILE_WRITE(iunit2,data,nx*ny*nz,kind(data),ierr)
        end if
     end do
     call BINARY_FILE_CLOSE(iunit2,ierr)
      
 case default
     stop "Unknown choice"
  end select
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  
end program editVolumeData
