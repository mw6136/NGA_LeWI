! ============================================================ !
!                         dmd.f90                              !
!                                                              !
! Author: Temistocle Grenga                                    !
! Date:   October 24, 2016                                     !
! ============================================================ !

program dmd_3D_p3
  use parallel
  use string
  use precision
  use geometry
  use partition
  use masks
  use fileio
  use data
  use masks
  !use dump_volume
  implicit none
  integer :: i, j, k, l, ii 
  integer :: iunit1, iunit2, ierr
  integer :: n_snap, start_file, step_file
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: nx1, ny1, nz1, nvars, nvars1, ntime1, np
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nd1_MOK, nd2_MOK, nd3_MOK, ii_MOK, ll_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK, str_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK, NTIME_MOK, idx_MOK
  real(WP) :: xx(3), varv, ddum
  integer :: size, ibuffer, ipart, ll
  character(len=str_short), dimension(:), pointer :: names
  character(len=str_medium) :: fname, folder, filename1
  integer, dimension(:), pointer :: ifile
  logical :: file_is_there
  real(WP), dimension(:), pointer :: x1, y1, z1
  real(SP), dimension(:), pointer :: xm1, ym1, zm1
  real(SP) :: phi4
  character(len=80) :: buffer
  integer, dimension(:,:,:), pointer :: iblank
  integer, dimension(:), pointer :: var_list
  real(WP), dimension(:,:,:), pointer :: data8
  real(SP), dimension(:,:,:), pointer :: data4
  real(WP) :: min_x, max_x, min_y, max_y, min_z, max_z
  integer :: minnx, maxnx, minny, maxny, minnz, maxnz
  logical :: subdomain

  if (irank == 1) print *, 'Starting dmd_3D Part 3'

  !Initializa MPI enironment
  call parallel_init()

  start_file = 8830 !6465 !8830 !15286 !21500
  step_file = 2
  n_snap = 200
  nvars = 1

  subdomain = .false.
  min_x = 0.00432_WP*7.0_WP
  max_x = 0.00432_WP*8.0_WP
  min_y = 0.00432_WP*-2.00_WP
  max_y = 0.00432_WP*2.00_WP
  min_z = 0.00432_WP*-0.5_WP
  max_z = 0.00432_WP*0.5_WP

  allocate(var_list(nvars))
  var_list = (/ 1/)
  !var_list = (/ 1, 2, 3, 4, 10, 12, 13, 14, 15 /)
  !var_list = (/ 1, 2, 3, 4, 10, 12 /)

  print*, 'Read Variables'
  write(folder, '(''U_'',I4.4,''_'', I2.2,''_tail/'')') n_snap, step_file
  l = start_file +step_file*n_snap
  write(fname, '(''vol_data.1_'', I8.8)') l
  filename1 = "volume_data/"//trim(fname)
  print*, trim(filename1)
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)
  call BINARY_FILE_READ(iunit1, ntime1, 1, kind(ntime1), ierr)
  call BINARY_FILE_READ(iunit1, nx1, 1,    kind(nx1), ierr)
  call BINARY_FILE_READ(iunit1, ny1, 1,    kind(ny1), ierr)
  call BINARY_FILE_READ(iunit1, nz1, 1,    kind(nz1), ierr)
  call BINARY_FILE_READ(iunit1, nvars1, 1, kind(nvars1), ierr)
  np = nx1*ny1*nz1
  if (irank == 1) print*, 'n var', nvars1
  if (irank == 1) print*, 'n x', nx1
  if (irank == 1) print*, 'n y', ny1
  if (irank == 1) print*, 'n z', nz1

  allocate(x1(nx1), y1(ny1), z1(nz1))
  do i=1,nx1
     call BINARY_FILE_READ(iunit1, x1(i), 1, kind(x1), ierr)
  end do
  do i=1,ny1
     call BINARY_FILE_READ(iunit1, y1(i), 1, kind(y1), ierr)
  end do
  do i=1,nz1
     call BINARY_FILE_READ(iunit1, z1(i), 1, kind(z1), ierr)
  end do

  if(.not.subdomain) then
     min_x = minval(x1); minnx = 1
     max_x = maxval(x1); maxnx = nx1
     min_y = minval(y1); minny = 1
     max_y = maxval(y1); maxny = ny1
     min_z = minval(z1); minnz = 1
     max_z = maxval(z1); maxnz = nz1
  else
     do i = 1,nx1
        minnx = i
        if (x1(i) .ge. min_x) EXIT
     end do
     do i = nx1,1,-1
        maxnx = i
        if (x1(i) .le. max_x) EXIT
     end do
     do i= 1,ny1
        minny = i
        if (y1(i) .ge. min_y) EXIT
     end do
     do i= ny1,1,-1
        maxny = i
        if (y1(i) .le. max_y) EXIT
     end do
     do i= 1,nz1
        minnz = i
        if (z1(i) .ge. min_z) EXIT
     end do
     do i= nz1,1,-1
        maxnz = i
        if (z1(i) .le. max_z) EXIT
     end do
     if (irank == 1 ) then
        print *, 'Subdomain'
        print *, 'x_min = ', x1(minnx); print *, 'x_max = ', x1(maxnx)
        print *, 'y_min = ', y1(minny); print *, 'y_max = ', y1(maxny)
        print *, 'z_min = ', z1(minnz); print *, 'z_max = ', z1(maxnz)
        print *, 'x_min = ', minnx
        print *, 'x_max = ', maxnx, nx1
        print *, 'y_min = ', minny
        print *, 'y_max = ', maxny, ny1
        print *, 'z_min = ', minnz
        print *, 'z_max = ', maxnz, nz1
     end if
  end if

  allocate(names(nvars1),ifile(nvars1))
  do i = 1,nvars1
     call BINARY_FILE_READ(iunit1, names(i), str_short, kind(names), ierr)
  end do
  do i = 1, 2
     call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
  end do
  if (irank == 1) print*, 'names', names

        ! ** Write the case file **
        fname = 'arts.case'
        filename1 = trim(folder)//trim(fname)
        open(iunit2,file=trim(filename1),form="formatted",iostat=ierr) !,status="REPLACE")
        if (irank == 1) print*, 'File opened ', filename1
        ! Write the case
        write(iunit2,'(a)') 'FORMAT'
        write(iunit2,'(a)') 'type: ensight gold'
        write(iunit2,'(a)') 'GEOMETRY'
        write(iunit2,'(a)') 'model: 1 1 geometry'
        write(iunit2,'(a)') 'VARIABLE'
        do i=1,nvars
           do j = 1,n_snap
              buffer = trim(adjustl(names(var_list(i)))) // "_"
              write(buffer(len_trim(buffer)+1:len_trim(buffer)+3),'(i3.3)') j
              write(iunit2,'(6a)') 'scalar per node: ', trim(buffer),' ',trim(buffer)
           end do
           write(iunit2,'(6a)') 'scalar per node: ', trim(names(var_list(i))),' ',trim(names(var_list(i)))
        end do
        write(iunit2,'(a)') 'TIME'
        write(iunit2,'(a)') 'time set: 1'
        write(iunit2,'(a)') 'number of steps: 1'
        write(iunit2,'(a)') 'filename start number: 1'
        write(iunit2,'(a)') 'filename increment: 1'
        print*, ddum
        write(iunit2,'(a,ES12.5)') 'time values:',ddum
        ! Close the file
        close(iclose(iunit2))
        print *, irank, 'arts.case completed'


        fname = 'geometry'
        filename1 = trim(folder)//trim(fname)
        ! ** Open the grid file to write **
        call BINARY_FILE_OPEN(iunit2,trim(filename1),"w",ierr)
        if (irank == 1) print*, 'File opened', filename1
        ! Write the geometry
        buffer = 'C Binary'
        call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        buffer = 'Ensight Gold Geometry File'
        call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        buffer = 'Structured Geometry'
        call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        buffer = 'node id off'
        call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        buffer = 'element id off'
        call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        ! Cell centers
        buffer = 'part'
        call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        ipart = 1
        call BINARY_FILE_WRITE(iunit2,ipart,1,kind(ipart),ierr)
        buffer = 'Complete geometry'
        call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        buffer = 'block rectilinear iblanked'
        call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
        print*, 'Header completed'

        call BINARY_FILE_WRITE(iunit2,maxnx-minnx+1,1,kind(nx),ierr)
        call BINARY_FILE_WRITE(iunit2,maxny-minny+1,1,kind(ny),ierr)
        call BINARY_FILE_WRITE(iunit2,maxnz-minnz+1,1,kind(nz),ierr)
        print*, 'N points completed'
  
        allocate(iblank(maxnx-minnx+1,maxny-minny+1,maxnz-minnz+1), xm1(maxnx-minnx+1), ym1(maxny-minny+1), zm1(maxny-minny+1))
        iblank = 1
        xm1(1:maxnx-minnx+1) = REAL(x1(minnx:maxnx),SP)
        ym1(1:maxny-minny+1) = REAL(y1(minny:maxny),SP)
        zm1(1:maxnz-minnz+1) = REAL(z1(minnz:maxnz),SP)

        call BINARY_FILE_WRITE(iunit2,xm1,maxnx-minnx+1,kind(xm1),ierr)
        call BINARY_FILE_WRITE(iunit2,ym1,maxny-minny+1,kind(ym1),ierr)
        call BINARY_FILE_WRITE(iunit2,zm1,maxnz-minnz+1,kind(zm1),ierr)
        call BINARY_FILE_WRITE(iunit2,iblank,(maxnx-minnx+1)*(maxny-minny+1)*(maxnz-minnz+1),kind(iblank),ierr)

        call BINARY_FILE_CLOSE(iunit2,ierr)
        print *, irank, 'geometry completed'
        deallocate(x1, y1, z1, xm1, ym1, zm1, iblank)

  allocate(data8(nx1,ny1,nz1), data4(maxnx-minnx+1,maxny-minny+1,maxnz-minnz+1))
  do ll = 1, var_list(1)-1
     call BINARY_FILE_READ(iunit1,data8,nx1*ny1*nz1,kind(data8),ierr)
  end do
  do i = 1,nvars
     fname = trim(adjustl(names(var_list(i))))
     filename1 = trim(folder)//trim(fname)
     inquire(file=filename1,exist=file_is_there)
     if (.not.(file_is_there)) &
        call BINARY_FILE_OPEN(ifile(i),trim(filename1),"w",ierr)
     if (irank == 1) print*, 'File opened', filename1

     buffer = trim(adjustl(names(var_list(i))))
     print*, buffer
     call BINARY_FILE_WRITE(ifile(i),buffer,80,kind(buffer),ierr)
     buffer = 'part'
     call BINARY_FILE_WRITE(ifile(i),buffer,80,kind(buffer),ierr)
     ibuffer = 1
     call BINARY_FILE_WRITE(ifile(i),ibuffer,1,kind(ibuffer),ierr)
     buffer = 'block'
     call BINARY_FILE_WRITE(ifile(i),buffer,80,kind(buffer),ierr)
     if (irank == 1) print*, 'Header completed'

     data8(:,:,:) = 0.0_WP
     call BINARY_FILE_READ(iunit1,data8,nx1*ny1*nz1,kind(data8),ierr)

     data4(1:maxnx-minnx+1, 1:maxny-minny+1, 1:maxnz-minnz+1) = REAL(data8(minnx:maxnx, minny:maxny, minnz:maxnz),SP)

     call BINARY_FILE_WRITE(ifile(i),data4,(maxnx-minnx+1)*(maxny-minny+1)*(maxnz-minnz+1),kind(data4),ierr)

     call BINARY_FILE_CLOSE(ifile(i),ierr)

     if ( i == nvars)  GO TO 213
     ll = var_list(i+1) - var_list(i) -1
     do j = 1,ll
        call BINARY_FILE_READ(iunit1,data8,nx1*ny1*nz1,kind(data8),ierr)
     end do
 213 CONTINUE
     print *, irank, names(var_list(i)), 'Completed'
  end do

  call MPI_FILE_CLOSE(iunit1, ierr)

  ! Finalize the parallel environment
  call parallel_final()
  print *, irank, 'Completed dmd_3D'

end program dmd_3D_p3
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
