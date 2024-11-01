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
  integer :: iunit1, ierr
  integer :: n_snap, start_file, step_file, nvars, nvars1, ntime1, nx1, ny1, nz1
  integer, dimension(:), pointer :: var_list
  real(WP), dimension(:), pointer :: datarms
  real(WP), dimension(:,:,:), pointer :: datamean, data8
  character(len=str_short) :: cdum
  integer :: idum
  real(WP) :: ddum
  character(len=str_short), dimension(:), pointer :: names
 
  !Initializa MPI enironment
  call parallel_init()

  start_file = 10640
  step_file = 1
  n_snap = 100
  nvars = 6

  allocate(var_list(nvars))
  var_list = (/ 1, 2, 3, 4, 10, 12 /)
  !var_list = (/ 1, 2, 3, 4, 10, 12 /)

  write(folder, '(''mean_wS_'',I4.4,''_'', I3.3,''/'')') n_snap, step_file

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
  call BINARY_FILE_CLOSE(iunit1,ierr)

  allocate(datamean(nx1, ny1, nvars), data8(nx1, ny1, nz1))
  l_var = nx1*ny1*nz1

  datamean(:,:,:) = 0.0_WP
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
     do j = 1, nx1+ny1+nz1
        call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
     end do

     do j = 1,nvars1
        call BINARY_FILE_READ(iunit1, cdum, str_short, kind(cdum), ierr)
     end do

     do j = 1, 2
        call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
     end do

     do j = 1, var_list(1)-1
        call BINARY_FILE_READ(iunit1,data8,l_var,kind(data8),ierr)
     end do

     do l = 1,nvars
        data8(:,:,:) = 0.0_WP
        call BINARY_FILE_READ(iunit1,data8,l_var,kind(data8),ierr)
        do j = 1,nz1
           datamean(:,:,l) =  datamean(:,:,l) +data8(:,:,j)/REAL(nz1*n_snap,WP)
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

  allocate(datarms(nvars), names(nvars1))
  datarms(:) = 0.0_WP

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
     do j = 1, nx1+ny1+nz1
        call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
     end do

     do j = 1,nvars1
        call BINARY_FILE_READ(iunit1, names(j), str_short, kind(cdum), ierr)
     end do

     do j = 1, 2
        call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
     end do

     do j = 1, var_list(1)-1
        call BINARY_FILE_READ(iunit1,data8,l_var,kind(data8),ierr)
     end do
     
     do l = 1,nvars
        data8(:,:,:) = 0.0_WP
        call BINARY_FILE_READ(iunit1,data8,l_var,kind(data8),ierr)
        do k = 1,nz1
           do j = 1,ny1
              do m = 1,nx1
                 datarms(l) =  datarms(l)+(data8(m,j,k)-datamean(m,j,l))**2.0_WP
              end do
           end do
        end do
        if ( l == nvars)  GO TO 212
        m = var_list(l+1) - var_list(l) -1
        do j = 1,m
           call BINARY_FILE_READ(iunit1,data8,l_var,kind(data8),ierr)
        end do
     end do
   212 CONTINUE
     call BINARY_FILE_CLOSE(iunit1, ierr)
  end do

  do l = 1,nvars
     datarms(l) = datarms(l)/REAL(l_var,WP)
     datarms(l) = DSQRT(datarms(l))
  end do

  do l = 1,nvars
     fname = trim(adjustl(names(var_list(l)))) // "_mean"
     filename1 = trim(folder)//trim(fname)
     call BINARY_FILE_OPEN(iunit1,trim(filename1),"w",ierr)
     call BINARY_FILE_WRITE(iunit1,datamean(:,:,l),nx1*ny1,kind(datamean),ierr)
     call BINARY_FILE_CLOSE(iunit1,ierr)
  end do
     
  fname = 'var_rms'
  filename1 = trim(folder)//trim(fname)
  open(iunit1,file=trim(filename1),form="formatted",iostat=ierr,status="REPLACE")
  do i = 1,nvars
     write(iunit1,'(E24.12E4)') datarms(i)
  end do
  close(iclose(iunit1))  

  ! Finalize the parallel environment
  call parallel_final()
  print *, irank, 'Completed post_dmd'

end program pre_dmd
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
