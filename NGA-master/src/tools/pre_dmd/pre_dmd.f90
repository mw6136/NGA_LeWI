! ============================================================ !
!                         dmd.f90                              !
!                                                              !
! Author: Temistocle Grenga                                    !
! Date:   March 21, 2017                                       !
! ============================================================ !

program pre_dmd
  use parallel
  use string
  use precision
  use geometry
  use partition
  use masks
  use fileio
  use data
  use masks

  implicit none
  character(len=str_medium) :: folder, filename1, filename_base, fname
  integer :: i, j, k, l, m, l_var
  integer :: iunit1, ierr, iunit2
  integer :: n_snap, start_file, step_file, nvars, nvars1, ntime1, nx1, ny1, nz1
  integer, dimension(:), pointer :: var_list
  real(WP), dimension(:), pointer :: datarms, datamean, datamagn
  real(WP), dimension(:,:,:), pointer :: data8
  real(WP), dimension(:), pointer :: x1, y1, z1, x2, y2, z2
  character(len=str_short) :: cdum
  integer :: idum
  real(WP) :: ddum
  character(len=str_short), dimension(:), pointer :: names
!  integer :: minnx, maxnx, minny, maxny, minnz, maxnz
!  real(WP) :: min_x, min_y, max_y, min_z, max_z

  !Initializa MPI enironment
  call parallel_init()

  start_file = 6465
  step_file = 1
  n_snap = 400
  nvars = 16

!  min_x = 0.00432_WP*5.0_WP
!  min_y = -0.00432_WP
!  max_y = 0.00432_WP
!  min_z = -0.00432_WP
!  max_z = 0.00432_WP

  allocate(var_list(nvars))
  var_list = (/ 1:16 /)
  !var_list = (/ 1, 2, 3, 4, 10, 12 /)

  write(folder, '(''Ured2_'',I4.4,''_'', I2.2,''/'')') n_snap, step_file

  filename_base = 'vol_data.1'
  l = start_file +step_file*n_snap
  write(fname, '(''_'', I8.8)') l
  fname = trim(filename_base)//trim(fname)
  filename1 = "volume_data/"//trim(fname)
  print*, filename1

  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)
  call BINARY_FILE_READ(iunit1, ntime1, 1, kind(ntime1), ierr)
  call BINARY_FILE_READ(iunit1, nx1, 1,    kind(nx1), ierr)
  call BINARY_FILE_READ(iunit1, ny1, 1,    kind(ny1), ierr)
  call BINARY_FILE_READ(iunit1, nz1, 1,    kind(nz1), ierr)
  call BINARY_FILE_READ(iunit1, nvars1, 1, kind(nvars1), ierr)
  allocate(x1(0:nx1+1), y1(0:ny1+1), z1(0:nz1+1))
  allocate(x2(nx1), y2(ny1), z2(nz1))
  call BINARY_FILE_READ(iunit1, x2, nx1, kind(x1), ierr)
  call BINARY_FILE_READ(iunit1, y2, ny1, kind(y1), ierr)
  call BINARY_FILE_READ(iunit1, z2, nz1, kind(z1), ierr)
  call BINARY_FILE_CLOSE(iunit1,ierr)

  print*, minval(x2), maxval(x2)
  print*, minval(y2), maxval(y2)
  print*, minval(z2), maxval(z2)

  x1(1:nx1) = x2(1:nx1)
  y1(1:ny1) = y2(1:ny1)
  z1(1:ny1) = z2(1:ny1)
  x1(0) = x1(1) - (x1(2) - x1(1))
  y1(0) = y1(1) - (y1(2) - y1(1))
  z1(0) = z1(1) - (z1(2) - z1(1))
  x1(nx1+1) = x1(nx1) + (x1(nx1) - x1(nx1-1))
  y1(ny1+1) = y1(ny1) + (y1(ny1) - y1(ny1-1))
  z1(nz1+1) = z1(nz1) + (z1(nz1) - z1(nz1-1))

  print*, x1(0), x1(1), x1(2)
  print*, x1(nx1-1), x1(nx1), x1(nx1+1)
  print*, y1(0), y1(1), y1(2)
  print*, y1(ny1-1), y1(ny1), y1(ny1+1)
  print*, z1(0), z1(1), z1(2)
  print*, z1(nz1-1), z1(nz1), z1(nz1+1)

!  do i= 1,nx1
!     minnx = i
!     if (x1(i) .ge. min_x) EXIT
!  end do
!  maxnx = nx1-1
!  do i= 1,ny1
!     minny = i
!     if (y1(i) .ge. min_y) EXIT
!  end do
!  do i= ny1,1,-1
!     maxny = i
!     if (y1(i) .le. max_y) EXIT
!  end do
!  do i= 1,nz1
!     minnz = i
!     if (z1(i) .ge. min_z) EXIT
!  end do
!  minnz = max(minnz,2)
!  max_z = min(-z1(minnz),max_z)
!  do i= nz1,1,-1
!     maxnz = i
!     if (z1(i) .le. max_z) EXIT
!  end do
!  maxnz = min(maxnz,nz1-1)

  allocate(datamean(nvars), datamagn(nvars), data8(nx1, ny1, nz1))
  allocate(datarms(nvars), names(nvars1))
  datarms(:) = 0.0_WP
  l_var = nx1*ny1*nz1

  print*, nx1, ny1, nz1, l_var

  datamean(:) = 0.0_WP; datamagn = 0.0_WP
  do i = 1,n_snap
     l = start_file +step_file*(i-1)
     write(fname, '(''_'', I8.8)') l
     fname = trim(filename_base)//trim(fname)
     filename1 = "volume_data/"//trim(fname)
     print*, filename1

     call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

     do j = 1, 5
        call BINARY_FILE_READ(iunit1, idum, 1, kind(idum), ierr)
     end do
     call BINARY_FILE_READ(iunit1, x2, nx1, kind(x1), ierr)
     call BINARY_FILE_READ(iunit1, y2, ny1, kind(y1), ierr)
     call BINARY_FILE_READ(iunit1, z2, nz1, kind(z1), ierr)

     do j = 1,nvars1
        call BINARY_FILE_READ(iunit1, names(j), str_short, kind(names), ierr)
     end do
     !print*, names

     do j = 1, 2
        call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
        !print*, ddum
     end do

     do j = 1, var_list(1)-1
        call BINARY_FILE_READ(iunit1,data8,l_var,kind(data8),ierr)
     end do

     do l = 1,nvars
        data8(:,:,:) = 0.0_WP
        call BINARY_FILE_READ(iunit1,data8,l_var,kind(data8),ierr)
        !print*, var_list(l), minval(data8(:,:,:)), maxval(data8(:,:,:))
        do k = 1,nz1
           do j = 1,ny1
              do m = 1,nx1
                 datamean(l) =  datamean(l) +data8(m,j,k)*(0.5_WP*DABS(x1(m+1)-x1(m-1))*0.5_WP*DABS(y1(j+1)-y1(j-1))*0.5_WP*DABS(z1(k+1)-z1(k-1)))
              end do
           end do
        end do

        do k = 1,nz1
           do j = 1,ny1
              do m = 1,nx1
                 datamagn(l) =  datamagn(l) +data8(m,j,k)*data8(m,j,k)*(0.5_WP*DABS(x1(m+1)-x1(m-1))*0.5_WP*DABS(y1(j+1)-y1(j-1))*0.5_WP*DABS(z1(k+1)-z1(k-1)))
              end do
           end do
        end do
   
        if ( l == nvars)  GO TO 211
        m = var_list(l+1) - var_list(l) -1
        do j = 1,m
           call BINARY_FILE_READ(iunit1,data8,l_var,kind(data8),ierr)
        end do
     end do
   211 CONTINUE

     call BINARY_FILE_CLOSE(iunit1, ierr)

  end do

  do l = 1,nvars
     datamean(l) = datamean(l)/(REAL(n_snap,WP)* &
                DABS(0.5_WP*(x1(nx1+1)+x1(nx1))-0.5_WP*(x1(0)+x1(1)))* &
                DABS(0.5_WP*(y1(ny1+1)+y1(ny1))-0.5_WP*(y1(0)+y1(1)))* &
                DABS(0.5_WP*(z1(nz1+1)+z1(nz1))-0.5_WP*(z1(0)+z1(1))))
     datamagn(l) = datamagn(l)/(REAL(n_snap,WP)* &
                DABS(0.5_WP*(x1(nx1+1)+x1(nx1))-0.5_WP*(x1(0)+x1(1)))* &
                DABS(0.5_WP*(y1(ny1+1)+y1(ny1))-0.5_WP*(y1(0)+y1(1)))* &
                DABS(0.5_WP*(z1(nz1+1)+z1(nz1))-0.5_WP*(z1(0)+z1(1))))
     datamagn(l) = DSQRT(datamagn(l))
  end do

  print*, 'mean'
  print*, datamean
  print*, 'magn'
  print*, datamagn


!  do i = 1,n_snap
!     l = start_file +step_file*(i-1)
!     write(fname, '(''_'', I8.8)') l
!     fname = trim(filename_base)//trim(fname)
!     filename1 = "volume_data/"//trim(fname)
!     print*, filename1
!
!     call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)
!
!     do j = 1, 5
!        call BINARY_FILE_READ(iunit1, idum, 1, kind(idum), ierr)
!     end do
!     do j = 1, nx1+ny1+nz1
!        call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
!     end do
!
!     do j = 1,nvars1
!        call BINARY_FILE_READ(iunit1, names(j), str_short, kind(cdum), ierr)
!     end do
!
!     do j = 1, 2
!        call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
!     end do
!
!     do j = 1, var_list(1)-1
!        call BINARY_FILE_READ(iunit1,data8,l_var,kind(data8),ierr)
!     end do
!     
!     do l = 1,nvars
!        data8(:,:,:) = 0.0_WP
!        call BINARY_FILE_READ(iunit1,data8,l_var,kind(data8),ierr)
!        do k = 1,nz1
!           do j = 1,ny1
!              do m = 1,nx1
!                 datarms(l) =  datarms(l)+(data8(m,j,k)-datamean(l))**2.0_WP
!              end do
!           end do
!        end do
!        if ( l == nvars)  GO TO 212
!        m = var_list(l+1) - var_list(l) -1
!        do j = 1,m
!           call BINARY_FILE_READ(iunit1,data8,l_var,kind(data8),ierr)
!        end do
!     end do
!   212 CONTINUE
!     call BINARY_FILE_CLOSE(iunit1, ierr)
!  end do

!  do l = 1,nvars
!     datarms(l) = datarms(l)/REAL(l_var*n_snap,WP)
!     datarms(l) = DSQRT(datarms(l))
!  end do

!  do l = 1,nvars
!     fname = trim(adjustl(names(var_list(l)))) // "_mean"
!     filename1 = trim(folder)//trim(fname)
!     call BINARY_FILE_OPEN(iunit1,trim(filename1),"w",ierr)
!     call BINARY_FILE_WRITE(iunit1,datamean(l),nx1*ny1,kind(datamean),ierr)
!     call BINARY_FILE_CLOSE(iunit1,ierr)
!  end do
 
  fname = 'var_mean'
  filename1 = trim(folder)//trim(fname)
  open(iunit1,file=trim(filename1),form="formatted",iostat=ierr,status="REPLACE")
  do i = 1,nvars
     write(iunit1,'(E24.12E4)') datamean(i)
  end do
  close(iclose(iunit1))

  fname = 'var_magn'
  filename1 = trim(folder)//trim(fname)
  open(iunit2,file=trim(filename1),form="formatted",iostat=ierr,status="REPLACE")
  do i = 1,nvars
     write(iunit2,'(E24.12E4)') datamagn(i)
  end do
  close(iclose(iunit2))
   
!  fname = 'var_rms'
!  filename1 = trim(folder)//trim(fname)
!  open(iunit2,file=trim(filename1),form="formatted",iostat=ierr,status="REPLACE")
!  do i = 1,nvars
!     write(iunit2,'(E24.12E4)') datarms(i)
!  end do
!  close(iclose(iunit2))  

  ! Finalize the parallel environment
  call parallel_final()
  print *, irank, 'Completed post_dmd'

end program pre_dmd
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
