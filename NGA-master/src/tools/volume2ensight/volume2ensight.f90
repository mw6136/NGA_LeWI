! ============================================================ !
!                   volume2ensight.f90                         !
!                                                              !
! Author: Temistocle Grenga                                    !
! Date:   March 10, 2017                                       !
! ============================================================ !

program volume2ensight
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
  integer :: i, j, k, l, m
  integer :: iunit1, iunit2, ierr
  integer :: n_snap, start_file, step_file, nvars
  character(len=str_medium) :: fname, folder, filename1
  integer :: nx1, ny1, nz1, nvarsT, ntime1
  real(WP), dimension(:), pointer :: x1, y1, z1
  real(SP), dimension(:), pointer :: xm1, ym1, zm1
  integer, dimension(:,:,:), pointer :: iblank
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:), pointer :: var8
  real(SP), dimension(:,:,:), pointer :: var4
  real(WP) :: ddum
  integer, dimension(:), pointer :: var_list
  integer :: idum, ipart, avg, ll
  character(len=80) :: cdum, buffer
  logical :: file_is_there

  if (irank == 1) print *, 'converting dump volume to ensight'

  !Initializa MPI enironment
  call parallel_init()

  start_file = 8805
  step_file = 1
  n_snap = 430
  nvars = 10

  avg = FLOOR(REAL(n_snap,WP)/REAL(nproc,WP))

  if (avg*nproc .ne. n_snap) then
     if (irank == 1) print *, 'ERROR'
     if (irank == 1) print *, 'avg*nproc /= n_snap'
     if (irank == 1) print *, 'EXIT'
     call parallel_final()
     STOP
  end if
  if (irank == 1) print *, 'avg = ',avg

  allocate(var_list(nvars))
  var_list = (/ 1, 5, 7, 8, 10, 11, 12, 13, 14, 15 /)

  print*, 'Read Variables'
  write(folder, '(''ensight/'')')
 
  do m = 1, avg
     l = start_file +step_file*avg*(irank-1)+step_file*(m-1)
     write(fname, '(''vol_data.1_'', I8.8)') l
     filename1 = "volume_data/"//trim(fname)
     print*, irank, trim(filename1)

     ! Open the binary files
     call BINARY_FILE_OPEN(iunit1, trim(filename1), "r", ierr)

     ! Read data sizes
     call BINARY_FILE_READ(iunit1, ntime1, 1, kind(ntime1), ierr)
     call BINARY_FILE_READ(iunit1, nx1, 1,    kind(nx1), ierr)
     call BINARY_FILE_READ(iunit1, ny1, 1,    kind(ny1), ierr)
     call BINARY_FILE_READ(iunit1, nz1, 1,    kind(nz1), ierr)
     call BINARY_FILE_READ(iunit1, nvarsT, 1, kind(nvarsT), ierr)

     if (m == 1) allocate(x1(nx1), y1(ny1), z1(nz1), names(nvarsT))
     do i=1,nx1
        call BINARY_FILE_READ(iunit1, x1(i), 1, kind(x1), ierr)
     end do
     do i=1,ny1
        call BINARY_FILE_READ(iunit1, y1(i), 1, kind(y1), ierr)
     end do
     do i=1,nz1
        call BINARY_FILE_READ(iunit1, z1(i), 1, kind(y1), ierr)
     end do

     ! Read variable names
     do i = 1,nvarsT
        call BINARY_FILE_READ(iunit1, names(i), str_short, kind(names), ierr)
     end do
     if(irank == 1 ) print *, names

     do i = 1, 2
        call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
     end do

     allocate(var8(nx1,ny1,nz1), var4(nx1,ny1,nz1))

     do i = 1, var_list(1)-1
        call BINARY_FILE_READ(iunit1,var8,nx1*ny1*nz1,kind(var8),ierr)
     end do

     do i = 1,nvars
        if(irank == 1 ) print *, names(var_list(i))
        call BINARY_FILE_READ(iunit1,var8,nx1*ny1*nz1,kind(var8),ierr)
        if (irank == 1 ) then
           print *, 'minimum ', minval(var8)
           print *, 'maximum ', maxval(var8)
        end if
        fname = trim(folder)//trim(adjustl(names(var_list(i)))) // "_"
        write(fname(len_trim(fname)+1:len_trim(fname)+6),'(i6.6)') l
        print*,fname
        call BINARY_FILE_OPEN(iunit2,trim(fname),"w",ierr)
        cdum = trim(adjustl(names(i))) // "_1_"
        write(cdum(len_trim(cdum)+1:len_trim(cdum)+6),'(i6.6)') l
        call BINARY_FILE_WRITE(iunit2,cdum,80,kind(cdum),ierr)
        cdum = 'part'
        call BINARY_FILE_WRITE(iunit2,cdum,80,kind(cdum),ierr)
        idum = 1
        call BINARY_FILE_WRITE(iunit2,idum,1,kind(idum),ierr)
        cdum = 'block'
        call BINARY_FILE_WRITE(iunit2,cdum,80,kind(cdum),ierr)
        var4(1:nx1,1:ny1,1:nz1) = var8(1:nx1,1:ny1,1:nz1)
        call BINARY_FILE_WRITE(iunit2,var4,nx1*ny1*nz1,kind(var4),ierr)
        call BINARY_FILE_CLOSE(iunit2,ierr)   
        if ( l == nvars)  GO TO 213
        ll = var_list(i+1) - var_list(i) -1
        do j = 1,ll
           call BINARY_FILE_READ(iunit1,var8,nx1*ny1*nz1,kind(var8),ierr)
        end do
     end do

 213 CONTINUE
     ! Close the file  
     call BINARY_FILE_CLOSE(iunit1, ierr)

  end do

     if (irank == 1) then
        ! ** Write the case file **
        fname = 'arts.case'
        filename1 = trim(folder)//trim(fname)
        inquire(file=filename1,exist=file_is_there)
        if (file_is_there) go to 111
        open(iunit2,file=trim(filename1),form="formatted",iostat=ierr,status="REPLACE")
        if (irank == 1) print*, 'File opened', filename1
        ! Write the case
        write(iunit2,'(a)') 'FORMAT'
        write(iunit2,'(a)') 'type: ensight gold'
        write(iunit2,'(a)') 'GEOMETRY'
        write(iunit2,'(a)') 'model: 1 1 geometry'
        write(iunit2,'(a)') 'VARIABLE'
        do i=1,nvars
           do j = 1,n_snap
              l = start_file +step_file*(j-1)
              buffer = trim(adjustl(names(var_list(i)))) // "_"
              write(buffer(len_trim(buffer)+1:len_trim(buffer)+6),'(i6.6)') l
              write(iunit2,'(6a)') 'scalar per node: ', trim(buffer),' ',trim(buffer)
           end do
        end do
        write(iunit2,'(a)') 'TIME'
        write(iunit2,'(a)') 'time set: 1'
        write(iunit2,'(a)') 'number of steps: 1'
        write(iunit2,'(a)') 'filename start number: 1'
        write(iunit2,'(a)') 'filename increment: 1'
        ddum = 1.94615E-02
        write(iunit2,'(a,ES12.5)') 'time values:',ddum
        ! Close the file
        close(iclose(iunit2))
        print *, irank, 'arts.case completed'
        go to 115

 111 continue

        fname = 'add_to_arts.case'
        filename1 = trim(folder)//trim(fname)
        open(iunit2,file=trim(filename1),form="formatted",iostat=ierr,status="REPLACE")
        if (irank == 1) print*, 'File opened', filename1
        do i=1,nvars
           do j = 1,n_snap
              l = start_file +step_file*(j-1)
              buffer = trim(adjustl(names(var_list(i)))) // "_"
              write(buffer(len_trim(buffer)+1:len_trim(buffer)+6),'(i6.6)') l
              write(iunit2,'(6a)') 'scalar per node: ', trim(buffer),' ',trim(buffer)
           end do
        end do
        close(iclose(iunit2))
        print *, irank, 'arts.case completed'
     end if

 115 continue

     if (irank == nproc) then
        fname = 'geometry'
        filename1 = trim(folder)//trim(fname)
        inquire(file=filename1,exist=file_is_there)
        if (file_is_there) go to 112
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

        call BINARY_FILE_WRITE(iunit2,nx1,1,kind(nx),ierr)
        call BINARY_FILE_WRITE(iunit2,ny1,1,kind(ny),ierr)
        call BINARY_FILE_WRITE(iunit2,nz1,1,kind(nz),ierr)
        print*, 'N points completed'

        allocate(xm1(1:nx1),ym1(1:ny1),zm1(1:nz1))
        allocate(iblank(nx1,ny1,nz1))
        iblank = 1
        xm1(1:nx1) = x1(1:nx1)
        ym1(1:ny1) = y1(1:ny1)
        zm1(1:nz1) = z1(1:nz1)

        call BINARY_FILE_WRITE(iunit2,xm1,nx1,kind(xm1),ierr)
        call BINARY_FILE_WRITE(iunit2,ym1,ny1,kind(ym1),ierr)
        call BINARY_FILE_WRITE(iunit2,zm1,nz1,kind(zm1),ierr)
        call BINARY_FILE_WRITE(iunit2,iblank,nx1*ny1*nz1,kind(iblank),ierr)

        call BINARY_FILE_CLOSE(iunit2,ierr)
        print *, irank, 'geometry completed'
        deallocate(xm1, ym1, zm1, iblank)
     end if

 112 continue

  deallocate(x1, y1, z1)
  deallocate(names)

  ! Finalize the parallel environment
  call parallel_final()
  print *, irank, 'Completed volume2ensight'


end program volume2ensight
