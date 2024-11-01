program mergeInflow
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: ntime1,ny,nz,ntime2,ny2,nz2,nvar1,nvar2,icyl,icyl2,ntime3,n,nvar3
  character(len=str_short), dimension(:), pointer :: names,names2
  real(WP), dimension(:,:), pointer :: data
  real(WP), dimension(:), pointer :: y,z,y2,z2
  integer :: iunit1,iunit2,iunit3,ierr,var,var2,itime
  character(len=str_medium) :: filename1,filename2,filename3
  character(len=str_short) :: varname,varname2
  real(WP) :: dt,time,dt2,time2,value

  ! Read file name from standard input
  print*,'======================'
  print*,'| ARTS - inflow Merger |'
  print*,'======================'
  print*
  print "(a17,$)", " inflow file 1 : "
  read "(a)", filename1
  print "(a17,$)", " inflow file 2 : "
  read "(a)", filename2
  print "(a27,$)", " inflow file when merged : "
  read "(a)", filename3

  ! ** Open the data file 1 to read **
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit1,ntime1,1,kind(ntime1),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar1,1,kind(nvar1),ierr)
  print*,'Grid :',ntime1,'time',ny,'x',nz
  
  ! Read additional stuff
  call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  print*,'Inflow file at time :',time
  print*,'Timestep size :',dt

  ! Read variables
  allocate(names(nvar1))
  do var=1,nvar1
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables : ',names
  
  ! ** Open the data file 2 to read **
  call BINARY_FILE_OPEN(iunit2,trim(filename2),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit2,ntime2,1,kind(ntime2),ierr)
  call BINARY_FILE_READ(iunit2,ny2,1,kind(ny2),ierr)
  call BINARY_FILE_READ(iunit2,nz2,1,kind(nz2),ierr)
  call BINARY_FILE_READ(iunit2,nvar2,1,kind(nvar2),ierr)
  if (ny.ne.ny2 .or. nz.ne.nz2) then
     print *, "WARNING: sizes not the same between files, errors will exist"
  end if
  print*,'Grid :',ntime2,'time',ny2,'x',nz2
  
  ! Read additional stuff
  call BINARY_FILE_READ(iunit2,dt2,1,kind(dt2),ierr)
  call BINARY_FILE_READ(iunit2,time2,1,kind(time2),ierr)
  if (dt.ne.dt2) then
     print *, "WARNING: timesteps not the same between files, errors will exist"
  end if
  print*,'Inflow file at time :',time2
  print*,'Timestep size :',dt2

  ! Read variables
  allocate(names2(nvar2))
  if (nvar1.ne.nvar2) then
     print *, "WARNING: nvar not the same between files, errors may exist"
  end if
  do var=1,nvar2
     call BINARY_FILE_READ(iunit2,names2(var),str_short,kind(names2),ierr)
  end do
  
  ! Add times in preparation for merging
  ntime3 = ntime1+ntime2
  
  ! Pull grid information
  allocate(y(ny+1))
  allocate(z(nz+1))
  call BINARY_FILE_READ(iunit1,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit1,y,ny+1,kind(y),ierr)
  call BINARY_FILE_READ(iunit1,z,nz+1,kind(z),ierr)
  allocate(y2(ny+1))
  allocate(z2(nz+1))
  call BINARY_FILE_READ(iunit2,icyl2,1,kind(icyl2),ierr)
  call BINARY_FILE_READ(iunit2,y2,ny+1,kind(y2),ierr)
  call BINARY_FILE_READ(iunit2,z2,nz+1,kind(z2),ierr)
  
  ! Write new header
  time = 0.0_WP
  call BINARY_FILE_OPEN(iunit3,trim(filename3),"w",ierr)
  call BINARY_FILE_WRITE(iunit3,ntime3,1,kind(ntime3),ierr)
  call BINARY_FILE_WRITE(iunit3,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit3,nz,1,kind(nz),ierr)
  call BINARY_FILE_WRITE(iunit3,nvar1,1,kind(nvar1),ierr)
  call BINARY_FILE_WRITE(iunit3,dt,1,kind(dt),ierr)
  call BINARY_FILE_WRITE(iunit3,time,1,kind(time),ierr)
  do var=1,nvar1
     call BINARY_FILE_WRITE(iunit3,names(var),str_short,kind(names),ierr)
  end do
  call BINARY_FILE_WRITE(iunit3,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_WRITE(iunit3,y,ny+1,kind(y),ierr)
  call BINARY_FILE_WRITE(iunit3,z,nz+1,kind(z),ierr)

  ! Read and write the data field
  allocate(data(ny,nz))
  do n=1,ntime1
     do var=1,nvar1
        data = 0.0_WP
        call BINARY_FILE_READ (iunit1,data(:,:),ny*nz,kind(data),ierr)
        call BINARY_FILE_WRITE(iunit3,data(:,:),ny*nz,kind(data),ierr)
     end do
  end do
  do n=ntime1+1,ntime3
     do var=1,nvar1
        data = 0.0_WP
        call BINARY_FILE_READ (iunit2,data(:,:),ny*nz,kind(data),ierr)
        call BINARY_FILE_WRITE(iunit3,data(:,:),ny*nz,kind(data),ierr)
     end do
  end do
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  call BINARY_FILE_CLOSE(iunit3,ierr)
  
end program mergeInflow
