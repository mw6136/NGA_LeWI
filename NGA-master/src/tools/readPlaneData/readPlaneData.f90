! -------------------------------------------------------------------------- !
!                           READPLANEDATA.F90                                !
!     Purpose     : Reads header information from binary plane data files    !
!     Date        : June 14, 2016                                            !
!     Modified by : Jonathan F. MacArt                                       !
! -------------------------------------------------------------------------- !
program readPlaneData
  use precision
  use string
  use fileio
  use cli_reader
  implicit none

  character(len=str_medium) :: filename
  integer :: iunit, ierr
  integer :: ny, nz, nplanes, nvars, ntime
  integer :: it,n,i,j,k
  real(WP) :: dt, time, Lz
  real(WP), dimension(:),       pointer :: x
  real(WP), dimension(:),       pointer :: y
  real(WP), dimension(:,:,:,:), pointer :: data
  character(len=str_short), dimension(:), pointer :: names

  call get_command_argument(1,filename)

  print *, "  plane file : ", trim(filename)

  ! Open the binary file
  call BINARY_FILE_OPEN(iunit, trim(filename), "r", ierr)

  ! Read data sizes
  call BINARY_FILE_READ(iunit, ntime, 1, kind(ntime), ierr)
  call BINARY_FILE_READ(iunit, ny, 1, kind(ny), ierr)
  call BINARY_FILE_READ(iunit, nz, 1, kind(nz), ierr)
  call BINARY_FILE_READ(iunit, nplanes, 1, kind(nplanes), ierr)
  call BINARY_FILE_READ(iunit, nvars, 1, kind(nvars), ierr)

  print *, 'ntime :   ', ntime
  print *, 'ny :      ', ny
  print *, 'nz :      ', nz
  print *, 'nplanes : ', nplanes
  print *, 'nvars :   ', nvars

  ! Read y-locations
  allocate(y(ny))
  do j=1,ny
     call BINARY_FILE_READ(iunit, y(j), 1, kind(y), ierr)
  end do
  print *, 'y = ', y

  ! Read Lz
  call BINARY_FILE_READ(iunit, Lz, 1, kind(Lz), ierr)  
  print *, 'Lz = ', Lz

  ! Read x-locations of planes
  allocate(x(nplanes))
  do n=1,nplanes
     call BINARY_FILE_READ(iunit, x(n), 1, kind(x), ierr)
  end do
  print *, 'planes at x = ', x

  ! Read variable names
  allocate(names(nvars))
  do n=1,nvars
     call BINARY_FILE_READ(iunit, names(n), str_short, kind(names), ierr)
  end do
  print *, 'Variables : ', names

  ! Read the data
  allocate(data(nvars,nplanes,ny,nz))
  do it=1,ntime
     ! Read time info
     call BINARY_FILE_READ(iunit, dt, 1, kind(dt), ierr)
     call BINARY_FILE_READ(iunit, time, 1, kind(time), ierr)
     print *, '****************************'
     print *, 'ntime : ', it
     print *, 'dt    : ', dt
     print *, 'time  : ', time
     do n=1,nvars
        print *, 'var=',n,names(n)
        do k=1,nz
           do j=1,ny
              do i=1,nplanes
                 call BINARY_FILE_READ(iunit, data(n,i,j,k), 1, kind(data), ierr)
              end do
           end do
        end do
!!$        do i=1,nplanes
!!$           print *, 'plane=',i,x(i)
!!$           do j=1,ny
!!$              print *, 'j=',j
!!$              do k=1,nz
!!$                 print *, data(n,i,j,k)
!!$              end do
!!$           end do
!!$        end do
     end do
  end do

  ! Close the file
  call BINARY_FILE_CLOSE(iunit, ierr)

end program readPlaneData
