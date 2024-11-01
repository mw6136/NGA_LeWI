program data2ensight
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz,nvar
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:), pointer :: data8
  real(SP), dimension(:,:,:), pointer :: data4
  integer :: iunit,iunit2,ierr
  integer :: var
  character(len=str_medium) :: filename1,filename2, directory
  real(WP) :: dt,time

  character(len=str_short), dimension(:), pointer :: new_names
  character(len=80) :: buffer
  integer :: var1,index1,ibuffer
  logical :: vel_found

  ! Read file name from standard input
  print*,'========================================='
  print*,'| ARTS - data to ENSIGHT GOLD converter |'
  print*,'========================================='
  print*
  print "(a13,$)", " data file : "
  read "(a)", filename1
  print "(a21,$)", " ensight directory : "
  read "(a)", directory

  !call CREATE_FOLDER(trim(directory))
  call system("mkdir -p "//trim(directory))

  ! ** Open the data file to read **
  call BINARY_FILE_OPEN(iunit,trim(filename1),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit,nvar,1,kind(nvar),ierr)
  print*,'Grid :',nx,'x',ny,'x',nz
  
  ! Read additional stuff
  call BINARY_FILE_READ(iunit,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit,time,1,kind(time),ierr)
  print*,'Data file at time :',time
  
  ! Read variable names
  allocate(names(nvar))
  vel_found = .false.
  do var=1,nvar
     call BINARY_FILE_READ(iunit,names(var),str_short,kind(names),ierr)
     if (trim(names(var)).eq.'U') vel_found = .true.
  end do
  print*,'Variables : ',names
  if (vel_found) then ! Data file
     index1 = 4
     var1 = 4
  else ! Optdata file
     index1 = 1
     var1 = 1
  end if

  ! Read data field
  allocate(data8(nx,ny,nz))
  allocate(data4(nx,ny,nz))
  allocate(new_names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit,data8,nx*ny*nz,kind(data8),ierr)
     print*,maxval(abs(data8)), ' at ', maxloc(abs(data8))

     ! ** Convert the data **
     data4(:,:,:) = data8(:,:,:)
     select case(trim(names(var)))
     case ('U')
        new_names(1) = 'U'
     case ('V')
        new_names(2) = 'V'
     case ('W')
        new_names(3) = 'W'
     case default
        new_names(index1) = names(var)
        index1 = index1+1
     end select

     ! ** Open the data file to write **
     ! Vector: velocity
     if (vel_found.and.var.lt.var1) then
        if (var.eq.1) then
           filename2 = trim(directory) // '/V'
           call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
           buffer = 'velocity'
           call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
           buffer = 'part'
           call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
           ibuffer = 1
           call BINARY_FILE_WRITE(iunit2,ibuffer,1,kind(ibuffer),ierr)
           buffer = 'block'
           call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        end if
        call BINARY_FILE_WRITE(iunit2,data4(:,:,:),nx*ny*nz,kind(data4),ierr)
        if (var.eq.3) then
           call BINARY_FILE_CLOSE(iunit2,ierr)
        end if
     end if

     ! Scalar: pressure, ... and all variables from Opdata
     if (var.ge.var1) then
        filename2 = trim(directory) // '/' // trim(new_names(var))
        call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
        buffer = new_names(var)
        call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        buffer = 'part'
        call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        ibuffer = 1
        call BINARY_FILE_WRITE(iunit2,ibuffer,1,kind(ibuffer),ierr)
        buffer = 'block'
        call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        call BINARY_FILE_WRITE(iunit2,data4(:,:,:),nx*ny*nz,kind(data4),ierr)
        call BINARY_FILE_CLOSE(iunit2,ierr)
     end if

  end do
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! ** Write the case file **
  iunit = iopen()
  filename2 = trim(directory) // '/arts.case'
  open(iunit,file=trim(filename2),form="formatted",iostat=ierr,status="REPLACE")
  ! Write the case
  write(iunit,'(a)') 'FORMAT'
  write(iunit,'(a)') 'type: ensight gold'
  write(iunit,'(a)') 'GEOMETRY'
  write(iunit,'(a)') 'model: 1 1 geometry'

  write(iunit,'(a)') 'VARIABLE'
  if (vel_found) write(iunit,'(a)') 'vector per node: Velocity V'
  do var=var1,nvar
     write(iunit,'(4a)') 'scalar per node: ', trim(new_names(var)),' ',trim(new_names(var))
  end do

  write(iunit,'(a)') 'TIME'
  write(iunit,'(a)') 'time set: 1'
  write(iunit,'(a)') 'number of steps: 1'
  write(iunit,'(a)') 'filename start number: 1'
  write(iunit,'(a)') 'filename increment: 1'
  write(iunit,'(a,ES12.5)') 'time values:',time
  ! Close the file
  close(iclose(iunit))
  
end program data2ensight
