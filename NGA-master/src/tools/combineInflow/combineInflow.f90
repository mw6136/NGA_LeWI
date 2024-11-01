! -------------------------------------------------------------------------- !
!                           COMBINEINFLOW.F90                                !
!     Purpose     : Combines two inflow files into one                       !
!     Options     : Regular jet (two inflow channels)                        !
!                   Coflowing slot jet (jetcart2 - three inflow channels)    !
!     Date        : February 29, 2016                                        !
!     Modified by : Jonathan F. MacArt                                       !
! -------------------------------------------------------------------------- !

program combineInflow
  use precision
  use string
  use fileio
  use parser
  use cli_reader
  implicit none

  character(len=str_medium) :: input_name
  character(len=str_long) :: filename,filename_pipe,filename_coflow
  integer  :: ierr,var,n,nn,j,nscalar
  integer  :: iunit,iunit_pipe,iunit_coflow
  integer  :: ntime,ntime_pipe,ntime_coflow
  integer  :: ny,ny_pipe,ny_coflow
  integer  :: nz,nz_pipe,nz_coflow
  integer  :: nvar,nvar_pipe,nvar_coflow
  integer  :: icyl
  real(WP) :: dt,dt_pipe,dt_coflow,dt_ratio
  real(WP) :: time,time_pipe,time_coflow,wtime1,wtime2
  character(len=str_short), dimension(:), pointer :: names,names_pipe,names_coflow
  real(WP), dimension(:,:), pointer :: inflow,inflow_pipe,inflow_coflow
  real(WP), dimension(:,:,:), pointer :: inflow_prev, inflow_next
  real(WP), dimension(:), pointer :: y,y_pipe,y_coflow
  real(WP), dimension(:), pointer :: z,z_pipe,z_coflow

  logical :: jetcart2

  ! Parse the input file
  call get_command_argument(1,input_name)
  call parser_init
  call parser_parsefile(input_name)

  ! Pipe or concentric cartesian jet?
  call parser_read('Concentric Cartesian jet',jetcart2,.false.)

  ! Open the pipe inflow
  call parser_read('Pipe inflow file',filename_pipe)
  call BINARY_FILE_OPEN(iunit_pipe,trim(filename_pipe),"r",ierr)
  if (ierr.ne.0) stop "combineInflow: cannot open pipe file"
  call BINARY_FILE_READ(iunit_pipe,ntime_pipe,1,kind(ntime_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,ny_pipe,1,kind(ny_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,nz_pipe,1,kind(nz_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,nvar_pipe,1,kind(nvar_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,dt_pipe,1,kind(dt_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,time_pipe,1,kind(time_pipe),ierr)
  allocate(names_pipe(nvar_pipe))
  do var=1,nvar_pipe
     call BINARY_FILE_READ(iunit_pipe,names_pipe(var),str_short,kind(names_pipe),ierr)
  end do

  ! Open the coflow inflow
  call parser_read('Coflow inflow file',filename_coflow)
  call BINARY_FILE_OPEN(iunit_coflow,trim(filename_coflow),"r",ierr)
  if (ierr.ne.0) stop "combineInflow: cannot open coflow file"
  call BINARY_FILE_READ(iunit_coflow,ntime_coflow,1,kind(ntime_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,ny_coflow,1,kind(ny_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,nz_coflow,1,kind(nz_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,nvar_coflow,1,kind(nvar_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,dt_coflow,1,kind(dt_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,time_coflow,1,kind(time_coflow),ierr)
  allocate(names_coflow(nvar_coflow))
  do var=1,nvar_coflow
     call BINARY_FILE_READ(iunit_coflow,names_coflow(var),str_short,kind(names_coflow),ierr)
  end do

  ! Check the consistency of the two inflows
  if (nz_pipe.ne.nz_coflow) stop "combineInflow: two inflows should have same nz"
  if (nvar_pipe.ne.nvar_coflow) stop "combineInflow: two inflows should have same variables"
!!$  if (dt_pipe.ne.dt_coflow) stop "combineInflow: two inflows should have same dt"
  do var=1,nvar_pipe
     if (trim(names_pipe(var)).ne.trim(names_coflow(var))) &
          stop "combineInflow: two inflows should have same variables"
  end do

  ! Grid parameters
  if (jetcart2) then
     ny = ny_pipe + (ny_coflow - 1)*2
  else
     ny = ny_pipe + ny_coflow + 1
  end if
  nz = nz_pipe
  nvar = nvar_pipe

  ! Time step matching
  dt_ratio = dt_pipe/dt_coflow
  print *, 'Timestep ratio',dt_ratio
  if (dt_ratio.gt.1.0_WP) then
     ntime_pipe = floor(real(ntime_pipe,WP)*dt_ratio)
  else
     ntime_coflow = floor(real(ntime_coflow,WP)/dt_ratio)
  end if
  ntime = min(ntime_pipe,ntime_coflow)
  dt = min(dt_pipe,dt_coflow)
  time = 0.0_WP
  allocate(names(nvar))
  names = names_pipe

  ! Print some stats
  print*
  print*,'Common parameters:'
  print*,'-> nz',nz
  print*
  print*,'Updated parameters:'
  write(*,'(a9,3a13)') ' ','Pipe','Coflow','Combined'
  write(*,'(a9,3ES13.5)') ' -> dt   ',dt_pipe,dt_coflow,dt
  write(*,'(a9,3i13)')    ' -> ntime',ntime_pipe,ntime_coflow,ntime
  write(*,'(a9,3i13)')    ' -> ny   ',ny_pipe,ny_coflow,ny
  print*
  print*,'Variables '
  print*,'-> ',names
  print*

  ! Read the grids
  allocate(y_pipe(ny_pipe+1))
  allocate(z_pipe(nz_pipe+1))
  call BINARY_FILE_READ(iunit_pipe,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit_pipe,y_pipe,ny_pipe+1,kind(y_pipe),ierr)
  call BINARY_FILE_READ(iunit_pipe,z_pipe,nz_pipe+1,kind(z_pipe),ierr)
  allocate(y_coflow(ny_coflow+1))
  allocate(z_coflow(nz_coflow+1))
  call BINARY_FILE_READ(iunit_coflow,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit_coflow,y_coflow,ny_coflow+1,kind(y_coflow),ierr)
  call BINARY_FILE_READ(iunit_coflow,z_coflow,nz_coflow+1,kind(z_coflow),ierr)

  ! Create the new grid
  allocate(y(ny+1))
  if (jetcart2) then
     y(1:ny_coflow) = y_coflow(1:ny_coflow)
     y(ny_coflow+1:ny_coflow+ny_pipe-1) = y_pipe(2:ny_pipe-1)
     do j=1,ny_coflow
        y(ny_coflow+ny_pipe+j-1) = -y_coflow(ny_coflow+1-j)
     end do
  else
     y(1:ny_pipe+1) = y_pipe
     y(ny_pipe+2:ny_pipe+ny_coflow+2) = y_coflow
  end if
!!$  do j=1,ny+1
!!$     print *, j, y(j)
!!$  end do
  allocate(z(nz+1))
  z = z_pipe

  ! Write the new header
  call parser_read('Combined inflow file',filename)
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
  call BINARY_FILE_WRITE(iunit,ntime,1,kind(ntime),ierr)
  call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
  call BINARY_FILE_WRITE(iunit,nvar,1,kind(nvar),ierr)
  call BINARY_FILE_WRITE(iunit,dt,1,kind(dt),ierr)
  call BINARY_FILE_WRITE(iunit,time,1,kind(time),ierr)
  do var=1,nvar
     call BINARY_FILE_WRITE(iunit,names(var),str_short,kind(names),ierr)
  end do
  call BINARY_FILE_WRITE(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_WRITE(iunit,y,ny+1,kind(y),ierr)
  call BINARY_FILE_WRITE(iunit,z,nz+1,kind(z),ierr)

  ! Read and write data field
  allocate(inflow(ny,nz))
  allocate(inflow_pipe(ny_pipe,nz_pipe))
  allocate(inflow_coflow(ny_coflow,nz_coflow))
  if (dt_ratio.gt.1.0_WP) then
     allocate(inflow_prev(ny_pipe,nz_pipe,nvar))
     allocate(inflow_next(ny_pipe,nz_pipe,nvar))
  else if (dt_ratio.lt.1.0_WP) then
     allocate(inflow_prev(ny_coflow,nz_coflow,nvar))
     allocate(inflow_next(ny_coflow,nz_coflow,nvar))
  end if
  time_pipe   = 0.0_WP
  time_coflow = 0.0_WP
  wtime1 = 1.0_WP
  wtime2 = 0.0_WP
  ! If interpolating, prep the interpolation steps
  if (dt_ratio.gt.1.0_WP) then
     ! Pipe has larger dt
     call BINARY_FILE_READ(iunit_pipe,inflow_prev(:,:,:),ny_pipe*nz_pipe*nvar,kind(inflow_prev),ierr)
     call BINARY_FILE_READ(iunit_pipe,inflow_next(:,:,:),ny_pipe*nz_pipe*nvar,kind(inflow_next),ierr)
     time_pipe = time_pipe + dt_pipe
  else if (dt_ratio.lt.1.0_WP) then
     ! Coflow has larger dt
     call BINARY_FILE_READ(iunit_coflow,inflow_prev(:,:,:),ny_coflow*nz_coflow*nvar,kind(inflow_prev),ierr)
     call BINARY_FILE_READ(iunit_coflow,inflow_next(:,:,:),ny_coflow*nz_coflow*nvar,kind(inflow_next),ierr)
     time_coflow = time_coflow + dt_coflow
  end if
  ! Write the header
  if (dt_ratio.ne.1.0_WP) then
     write(*,'(a6,5a13)') 'n','t_pipe','t_coflow','t_interp','wtime1','wtime2'
  else
     write(*,'(a6,2a13)') 'n','t_pipe','t_coflow'
  end if

  ! Loop over time
  do n=1,ntime
     if (dt_ratio.gt.1.0_WP) then
        ! Pipe has larger dt
        write(*,'(i6,5ES13.5)') n, time_pipe, time_coflow, wtime1*(time_pipe-dt_pipe)+wtime2*time_pipe, wtime1, wtime2
     else if (dt_ratio.lt.1.0_WP) then
        ! Coflow has larger dt
        write(*,'(i6,5ES13.5)') n, time_pipe, time_coflow, wtime1*(time_coflow-dt_coflow)+wtime2*time_coflow, wtime1, wtime2
     else
        write(*,'(i6,5ES13.5)') n, time_pipe, time_coflow
     end if
     do var=1,nvar
        inflow = 0.0_WP
        ! Read, interpolating in time if necessary
        if (dt_ratio.gt.1.0_WP) then
           ! Pipe has larger dt, interpolate it and read the coflow
           inflow_pipe(:,:) = wtime1*inflow_prev(:,:,var) + wtime2*inflow_next(:,:,var)
           call BINARY_FILE_READ(iunit_coflow,inflow_coflow(:,:),ny_coflow*nz_coflow,kind(inflow_coflow),ierr)
        else if (dt_ratio.lt.1.0_WP) then
           ! Coflow has larger dt, interpolate it and read the pipe
           call BINARY_FILE_READ(iunit_pipe  ,inflow_pipe  (:,:),ny_pipe  *nz_pipe  ,kind(inflow_pipe  ),ierr)
           inflow_coflow(:,:) = wtime1*inflow_prev(:,:,var) + wtime2*inflow_next(:,:,var)
        else
           call BINARY_FILE_READ(iunit_pipe  ,inflow_pipe  (:,:),ny_pipe  *nz_pipe  ,kind(inflow_pipe  ),ierr)
           call BINARY_FILE_READ(iunit_coflow,inflow_coflow(:,:),ny_coflow*nz_coflow,kind(inflow_coflow),ierr)
        end if
        ! Convert
        if (jetcart2) then
           inflow(1:ny_coflow,:) = inflow_coflow
           inflow(ny_coflow+1:ny_coflow+ny_pipe-1,:) = inflow_pipe(2:ny_pipe,:)
           inflow(ny_coflow+ny_pipe:ny_coflow+ny_pipe+ny_coflow-2,:) = inflow_coflow(2:ny_coflow,:)
        else
           inflow(1:ny_pipe,:) = inflow_pipe
           inflow(ny_pipe+2:ny_pipe+ny_coflow+1,:) = inflow_coflow
        end if
        ! Write
        call BINARY_FILE_WRITE(iunit,inflow(:,:),ny*nz,kind(inflow),ierr)
     end do

     ! Advance the time pointer
     if (dt_ratio.gt.1.0_WP) then
        time_coflow = time_coflow + dt_coflow
     else if (dt_ratio.lt.1.0_WP) then
        time_pipe = time_pipe + dt_pipe
     else
        time_pipe   = time_pipe   + dt_pipe
        time_coflow = time_coflow + dt_coflow
     end if

     ! Interpolate in time if necessary
     if (dt_ratio.gt.1.0_WP) then
        ! Pipe has larger dt; interpolate it to the coflow time
        if (time_coflow.gt.time_pipe) then
           inflow_prev = inflow_next
           call BINARY_FILE_READ(iunit_pipe,inflow_next(:,:,:),ny_pipe*nz_pipe*nvar,kind(inflow_next),ierr)
           time_pipe = time_pipe + dt_pipe
        end if
        wtime1 = (time_pipe-time_coflow)/dt_pipe
        wtime2 = 1.0_WP - wtime1
     else if (dt_ratio.lt.1.0_WP) then
        ! Coflow has larger dt, interpolate it to the pipe time
        if (time_pipe.gt.time_coflow) then
           inflow_prev = inflow_next
           call BINARY_FILE_READ(iunit_coflow,inflow_next(:,:,:),ny_coflow*nz_coflow*nvar,kind(inflow_next),ierr)
           time_coflow = time_coflow + dt_coflow
        end if
        wtime1 = (time_coflow-time_pipe)/dt_coflow
        wtime2 = 1.0_WP-wtime1
     end if
  end do

  ! Close the files
  call BINARY_FILE_CLOSE(iunit,ierr)
  call BINARY_FILE_CLOSE(iunit_pipe,ierr)
  call BINARY_FILE_CLOSE(iunit_coflow,ierr)

end program combineInflow
