! ============================================================ !
!                   volume_vort.f90                            !
!                                                              !
! Author: Temistocle Grenga                                    !
! Date:   April 27, 2017                                       !
! ============================================================ !

program volume_vort
  use parallel
  use string
  use precision
  use geometry
  use partition
  use masks
  use fileio
  use data
  use masks
  use vort_mod
  use derv_w3

  implicit none
  character(len=str_medium) :: filename1, folder, filename_base, fname
  integer :: iunit1, ierr, iunit2
  integer :: i, j, k, l, m
  integer :: n_snap, start_file, step_file, nvars, w_snap
  integer, dimension(:), pointer :: var_list
  integer :: idum, nvars1
  real(WP) :: ddum
  character(len=str_short) :: cdum
  integer, dimension(2) :: derv, dir
  character(len=80) :: buffer
  real(SP), dimension(:,:,:), pointer :: data4

  !Initializa MPI enironment
  call parallel_init()

  start_file = 6465 !8830 !15687 !21525
  step_file = 1
  n_snap = 400
  nvars = 3

  allocate(var_list(nvars))
  var_list = (/ 1, 2, 3 /)

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
  do i = 1,nvars1
     call BINARY_FILE_READ(iunit1, cdum, str_short, kind(cdum), ierr)
  end do
  do i = 1, 2
     call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
     print*, ddum
  end do

  allocate(data8(nx1,ny1,nz1), ddata(nx1,ny1,nz1), vort_x(nx1,ny1,nz1), vort_y(nx1,ny1,nz1), vort_z(nx1,ny1,nz1))
  data8(:,:,:) = 0.0_WP; vort_x(:,:,:) = 0.0_WP; vort_y(:,:,:) = 0.0_WP; vort_z(:,:,:) = 0.0_WP
  do m = 1, 3
     call BINARY_FILE_READ(iunit1, data8, nx1*ny1*nz1, kind(data8), ierr)
     do l = 1,2
        if (m == 1) then
           if (l == 1) then
              call dx_w3(2)
              vort_z(:,:,:) = vort_z(:,:,:) +ddata(:,:,:)
           else
              call dx_w3(3)
              vort_y(:,:,:) = vort_y(:,:,:) -ddata(:,:,:)
           end if
        end if
        if (m == 2) then
           if (l == 1) then
              call dx_w3(3)
              vort_x(:,:,:) = vort_x(:,:,:) +ddata(:,:,:)
           else
              call dx_w3(1)
              vort_z(:,:,:) = vort_z(:,:,:) -ddata(:,:,:)
           end if
        end if
        if (m == 3) then
           if (l == 1) then
              call dx_w3(1)
              vort_y(:,:,:) = vort_y(:,:,:) +ddata(:,:,:)
           else
              call dx_w3(2)
              vort_x(:,:,:) = vort_x(:,:,:) -ddata(:,:,:)
           end if
        end if     
     end do
     data8(:,:,:) = 0.0_WP
  end do
  print*, 'Vorticity evaluation completed'

  call BINARY_FILE_CLOSE(iunit1,ierr)
  deallocate(data8)
  allocate(data4(nx1,ny1,nz1))

  folder = 'vort_snap/'
  filename1 = 'vort_x_'
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
  data4(1:nx1,1:ny1,1:nz1) =  REAL(vort_x(1:nx1,1:ny1,1:nz1),SP)
  call BINARY_FILE_WRITE(iunit2,data4,nx1*ny1*nz1,kind(data4),ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  print*, irank, filename1, 'written'

  filename1 = 'vort_y_'
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
  data4(1:nx1,1:ny1,1:nz1) =  REAL(vort_y(1:nx1,1:ny1,1:nz1),SP)
  call BINARY_FILE_WRITE(iunit2,data4,nx1*ny1*nz1,kind(data4),ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  print*, irank, filename1, 'written'

  filename1 = 'vort_z_'
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
  data4(1:nx1,1:ny1,1:nz1) =  REAL(vort_z(1:nx1,1:ny1,1:nz1),SP)
  call BINARY_FILE_WRITE(iunit2,data4,nx1*ny1*nz1,kind(data4),ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  print*, irank, filename1, 'written'

  allocate(vort(nx1,ny1,nz1))
  do k = 1,nz1
     do j = 1,ny1
        do i =1,nx1
           vort(i,j,k) = DSQRT(vort_x(i,j,k)**2.0_WP+vort_y(i,j,k)**2.0_WP+vort_z(i,j,k)**2.0_WP)  
        end do
     end do
  end do
  filename1 = 'vort_'
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
  data4(1:nx1,1:ny1,1:nz1) =  REAL(vort(1:nx1,1:ny1,1:nz1),SP)
  call BINARY_FILE_WRITE(iunit2,data4,nx1*ny1*nz1,kind(data4),ierr)
  call BINARY_FILE_CLOSE(iunit2,ierr)
  print*, irank, filename1, 'written'
  deallocate(vort)

  deallocate(vort_x, vort_y, vort_z)

  call parallel_final()

 110 Continue

  print*, irank, 'Complete'

end program volume_vort
