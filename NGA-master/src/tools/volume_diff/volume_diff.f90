! ============================================================ !
!                   volume_diff.f90                            !
!                                                              !
! Author: Temistocle Grenga                                    !
! Date:   June, 2017                                           !
! ============================================================ !

program volume_diff
  use parallel
  use string
  use precision
  use geometry
  use partition
  use masks
  use fileio
  use data
  use masks
  use diff_mod
  use derv_w3_diff

  implicit none
  character(len=str_medium) :: filename1, folder, filename_base, fname
  integer :: iunit1, ierr, iunit2
  integer :: i, j, k, l, m, ll
  integer :: n_snap, start_file, step_file, nvars, w_snap
  integer :: var_list, diff_list
  integer :: idum, nvars1
  real(WP) :: ddum
  character(len=str_short) :: cdum
  integer, dimension(2) :: derv, dir
  character(len=80) :: buffer
  real(SP), dimension(:,:,:), pointer :: data4
  real(WP), dimension(:,:,:), pointer :: data8b
  character(len=str_short), dimension(:), pointer :: names

  !Initializa MPI enironment
  call parallel_init()

  start_file = 6465 !8830 !6465 !15687 !21525
  step_file = 1
  n_snap = 400
  nvars = 1

  var_list = 16
  diff_list = 28

  filename_base = 'vol_data.1'

  if (nproc /= n_snap+1) then
     if (irank == 1) print*, 'nproc .neq. nsnap+1'
        call parallel_final()
        GO TO 110
  end if

  w_snap = MOD(irank,n_snap+1)

  l = start_file +step_file*w_snap
  write(fname, '(''_'', I8.8)') l
  fname = trim(filename_base)//trim(fname)
  filename1 = "volume_data/"//trim(fname)
  print*, filename1
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

  call BINARY_FILE_READ(iunit1, idum, 1, kind(idum), ierr)
  call BINARY_FILE_READ(iunit1, nx1, 1, kind(idum), ierr)
  call BINARY_FILE_READ(iunit1, ny1, 1, kind(idum), ierr)
  call BINARY_FILE_READ(iunit1, nz1, 1, kind(idum), ierr)
  call BINARY_FILE_READ(iunit1, nvars1, 1, kind(nvars1), ierr)
  allocate (x1(nx1), y1(ny1), z1(nz1))
  if (irank == 1) print*, nx1, ny1, nz1
  call BINARY_FILE_READ(iunit1, x1, nx1, kind(x1), ierr)
  call BINARY_FILE_READ(iunit1, y1, ny1, kind(y1), ierr)
  call BINARY_FILE_READ(iunit1, z1, nz1, kind(z1), ierr)
  allocate(names(nvars1))
  do i = 1,nvars1
     call BINARY_FILE_READ(iunit1, names(i), str_short, kind(cdum), ierr)
  end do
  do i = 1, 2
     call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
     print*, ddum
  end do

  allocate(data8(nx1,ny1,nz1), ddata(nx1,ny1,nz1))
  allocate(datadx(nx1,ny1,nz1), datady(nx1,ny1,nz1), datadz(nx1,ny1,nz1))
  data8(:,:,:) = 0.0_WP
  do ll = 1, var_list-1
     call BINARY_FILE_READ(iunit1, data8, nx1*ny1*nz1,kind(data8),ierr)
  end do
  data8(:,:,:) = 0.0_WP
  call BINARY_FILE_READ(iunit1, data8, nx1*ny1*nz1, kind(data8), ierr)

  ddata(:,:,:) = 0.0_WP
  call dx_w3(1)
  datadx(1:nx1,1:ny1,1:nz1) = ddata(1:nx1,1:ny1,1:nz1)

  ddata(:,:,:) = 0.0_WP
  call dx_w3(2)
  datady(1:nx1,1:ny1,1:nz1) = ddata(1:nx1,1:ny1,1:nz1)

  ddata(:,:,:) = 0.0_WP
  call dx_w3(3)
  datadz(1:nx1,1:ny1,1:nz1) = ddata(1:nx1,1:ny1,1:nz1)

  allocate(data8b(nx1,ny1,nz1))
  data8b(:,:,:) = 0.0_WP

  allocate(data4(nx1,ny1,nz1))

  if (names(var_list) == 'ZMIX') then
     print*, 'Into ZMIX'

     do k = 1,nz1
        do j = 1,ny1
           do i = 1,nx1
              data8b(i,j,k) = (datadx(i,j,k)**2.0_WP) +(datady(i,j,k)**2.0_WP) +(datadz(i,j,k)**2.0_WP)
           end do
        end do
     end do
     do ll = 1, (diff_list -var_list -1)
        call BINARY_FILE_READ(iunit1, data8, nx1*ny1*nz1,kind(data8),ierr)
     end do
     data8(:,:,:) = 0.0_WP
     call BINARY_FILE_READ(iunit1, data8, nx1*ny1*nz1, kind(data8), ierr)
     data8b(1:nx1,1:ny1,1:nz1) = 2.0_WP*data8(1:nx1,1:ny1,1:nz1)*data8b(1:nx1,1:ny1,1:nz1)
     data4(1:nx1,1:ny1,1:nz1) = REAL(data8b(1:nx1,1:ny1,1:nz1),SP)

  else

     do ll = 1, (diff_list -var_list -1)
        call BINARY_FILE_READ(iunit1, data8, nx1*ny1*nz1,kind(data8),ierr)
     end do
     data8b(:,:,:) = 0.0_WP
     call BINARY_FILE_READ(iunit1, data8b, nx1*ny1*nz1, kind(data8), ierr)

     data8(1:nx1,1:ny1,1:nz1) = data8b(1:nx1,1:ny1,1:nz1)*datadx(1:nx1,1:ny1,1:nz1)
     call dx_w3(1)
     datadx(:,:,:) = 0.0_WP
     datadx(1:nx1,1:ny1,1:nz1) = ddata(1:nx1,1:ny1,1:nz1)

     data8(:,:,:) = 0.0_WP
     data8(1:nx1,1:ny1,1:nz1) = data8b(1:nx1,1:ny1,1:nz1)*datady(1:nx1,1:ny1,1:nz1)
     call dx_w3(2)
     datady(:,:,:) = 0.0_WP
     datady(1:nx1,1:ny1,1:nz1) = ddata(1:nx1,1:ny1,1:nz1)

     data8(:,:,:) = 0.0_WP
     data8(1:nx1,1:ny1,1:nz1) = data8b(1:nx1,1:ny1,1:nz1)*datadz(1:nx1,1:ny1,1:nz1)
     call dx_w3(3)
     datadz(:,:,:) = 0.0_WP
     datadz(1:nx1,1:ny1,1:nz1) = ddata(1:nx1,1:ny1,1:nz1)

     data4(1:nx1,1:ny1,1:nz1) = REAL(datadx(1:nx1,1:ny1,1:nz1),SP) +REAL(datady(1:nx1,1:ny1,1:nz1),SP) +REAL(datadz(1:nx1,1:ny1,1:nz1),SP)

  end if

  call BINARY_FILE_CLOSE(iunit1,ierr)
  deallocate(data8, data8b, ddata)

  folder = 'diff_snap/'
  filename1 = trim(adjustl(names(var_list))) // "_"
  write(filename1(len_trim(filename1)+1:len_trim(filename1)+6),'(i6.6)') start_file+w_snap
  filename1 = trim(folder)//trim(filename1)
  call BINARY_FILE_OPEN(iunit2,trim(filename1),"w",ierr)
  call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
  buffer = 'part'
  call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
  idum = 1
  call BINARY_FILE_WRITE(iunit2,idum,1,kind(idum),ierr)
  buffer = 'block'
  call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
  call BINARY_FILE_WRITE(iunit2,data4,nx1*ny1*nz1,kind(data4),ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  print*, irank, filename1, 'written'

  deallocate(datadx, datady, datadz)
  deallocate(data4)

  call parallel_final()

 110 Continue

  print*, irank, 'Complete'

end program volume_diff
