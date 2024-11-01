! ============================================================ !
!                         dmd.f90                              !
!                                                              !
! Author: Temistocle Grenga                                    !
! Date:   September 7, 2016                                    !
! ============================================================ !

program dmd_3D_p1
  use parallel
  use string
  use precision
  use geometry
  use partition
  use masks
  use fileio
  use data
  !use dump_volume
  implicit none
  character(len=str_long) :: strname
  integer :: i, j, k, l, m, ll ,pr
  integer :: iunit1, iunit2, ierr, iunit3
  integer :: nx1, ny1, nz1, nvarsT, nvars1, ntime1
  integer :: n_snap, l_vec, l_var
  integer, dimension(:), pointer :: var_list, nelc
  integer :: H_nel, lvec
  integer, dimension(:), pointer :: H_lim, i_idx, j_idx
  integer, dimension(:,:,:), pointer :: H_idx
  real(WP) :: avg, U_bulk(3), T_bulk, P_bulk, ddum
  real(WP) :: min_x, max_x, min_y, max_y, min_z, max_z
  integer :: minnx, maxnx, minny, maxny, minnz, maxnz
  real(WP), dimension(:), pointer :: x1, y1, z1, U0
  real(WP), dimension(:,:), pointer :: H, H1, H_recv
  real(WP), dimension(:,:,:), pointer :: xi, xj
  character(len=str_short), dimension(:), pointer :: names
  character(len=str_short) :: cdum
  integer, dimension(MPI_STATUS_SIZE) :: status
  real(WP) :: curr_dt, curr_time
  logical :: dimensional, subdomain
  character(len=str_medium) :: filename1, fname, folder, arg
  integer :: start_file, step_file, idum

  !if(irank == 1 ) print *, 'Starting dmd_3D'

  !Initialize MPI environment
  call parallel_init()

  dimensional = .false.
  U_bulk(1:3) = 23.35_WP !34.86_WP NP !23.35_WP Low Ka !93.44_WP High Ka
  !U_bulk(1) = 21.1357928603924_WP
  !U_bulk(2) = 5.15000756525904_WP
  !U_bulk(3) = 2.65841626242270_WP
  T_bulk = 300.0_WP
  P_bulk = 101325_WP

  subdomain = .true.
  min_x = 0.00432_WP*1.0_WP
  max_x = 0.00432_WP*2.0_WP
  min_y = 0.00432_WP*-0.75_WP
  max_y = 0.00432_WP*0.75_WP
  min_z = 0.00432_WP*-0.5_WP
  max_z = 0.00432_WP*0.5_WP

  start_file = 8830 !6465 !8830 !NP NON-REACT 15686  REACTIVE 21500
  step_file = 2
  n_snap = 200
  nvars1 = 3
  allocate(var_list(nvars1))
  var_list = (/ 1, 2, 3 /) 

  write(folder, '(''U_'',I4.4,''_'', I2.2,''_head/'')') n_snap, step_file

  write(fname, '(''vol_data.1_'', I8.8)') start_file
  filename1 = "volume_data/"//trim(fname)
  if(irank == 1 ) print *, filename1

  ! Open the binary files
  call BINARY_FILE_OPEN(iunit1, trim(filename1), "r", ierr)

  ! Read data sizes
  call BINARY_FILE_READ(iunit1, ntime1, 1, kind(ntime1), ierr)
  call BINARY_FILE_READ(iunit1, nx1, 1,    kind(nx1), ierr)
  call BINARY_FILE_READ(iunit1, ny1, 1,    kind(ny1), ierr)
  call BINARY_FILE_READ(iunit1, nz1, 1,    kind(nz1), ierr)
  call BINARY_FILE_READ(iunit1, nvarsT, 1, kind(nvars1), ierr)
  if(irank == 1 ) then
     print *, 'ntime :   ', ntime1
     print *, 'nx :      ', nx1
     print *, 'ny :      ', ny1
     print *, 'nz :      ', nz1
     print *, 'nvars :   ', nvarsT
  end if
  if(irank == 1 ) print *, 'nvars :   ', nvars1, nvarsT

  ! Read locations in all directions
  allocate(x1(nx1), y1(ny1), z1(nz1))
  do i=1,nx1
     call BINARY_FILE_READ(iunit1, x1(i), 1, kind(x1), ierr)
  end do
  do i=1,ny1
     call BINARY_FILE_READ(iunit1, y1(i), 1, kind(y1), ierr)
  end do
  do i=1,nz1
     call BINARY_FILE_READ(iunit1, z1(i), 1, kind(y1), ierr)
  end do
  if(irank == 1 ) then
     print *, 'DATA'
     print *, 'x_min = ', minval(x1); print *, 'x_max = ', maxval(x1)
     print *, 'y_min = ', minval(y1); print *, 'y_max = ', maxval(y1)
     print *, 'z_min = ', minval(z1); print *, 'z_max = ', maxval(z1)
  end if

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

  deallocate(x1, y1, z1)

  ! Read variable names
  allocate(names(nvarsT))
  do i = 1,nvarsT
        call BINARY_FILE_READ(iunit1, names(i), str_short, kind(names), ierr)
  end do
  if (irank == 1 ) then
     print *, names
     do i = 1,nvars1
        print *, names(var_list(i))
     end do
  end if

!  ! Close the file  
  call BINARY_FILE_CLOSE(iunit1, ierr)

  allocate(U0(nvars1))
  if (dimensional) then
     U0 = 1.0D0
  else
     do i = 1,nvars1
        if (names(var_list(i)) == 'U') then
           U0(i) = U_bulk(1)
        else if (names(var_list(i)) == 'V') then
           U0(i) = U_bulk(2)
        else if (names(var_list(i)) == 'W') then
           U0(i) = U_bulk(3)
        else if (names(var_list(i)) == 'T') then
           U0(i) = T_bulk
        else if (names(var_list(i)) == 'P') then
           U0(i) = P_bulk
        else
           U0(i) = 1.0D0
        end if
     end do
  end if
  if(irank == 1 ) print *, 'U0 = ', U0

  if(irank == 1 ) print *, 'n_snap = ', n_snap

  l_vec = nx1*ny1*nz1*nvars1
  l_var = nx1*ny1*nz1
  !if(irank == 1 ) print *, 'vector length = ', l_vec 

  !if(irank == 1 ) print *, 'number of cores = ', nproc

  H_nel = 0
  do i=n_snap,1,-1
     H_nel = H_nel+i
  end do
  !if(irank == 1 ) print *, 'number of elements in H (symm) matrix = ', H_nel
  H_nel = H_nel +n_snap
  !if(irank == 1 ) print *, 'Total number of elements = ', H_nel

  avg = REAL(H_nel,WP)/REAL(nproc,WP)
  if (avg.lt.1.0D0) then
     print*,'WARNING: nproc > H_nel'
     go to 110
  end if
  allocate(H_lim(nproc+1))
  H_lim(1) = 1
  H_lim(nproc+1) = H_nel+1
  do i=2,nproc
     H_lim(i) = 1 + FLOOR(avg*(i-1))
  end do
  allocate(nelc(nproc))
  do i=1,nproc
     nelc(i) = H_lim(i+1) - H_lim(i)
  end do
  !print*, irank, 'number of elements = ', nelc(irank)
  !if (nelc(irank) == 0) GO TO 111

  allocate(i_idx(nelc(irank)), j_idx(nelc(irank)))
  k = 0; l = 0
  do i = 1,n_snap
     do j = i,n_snap+1
        k = k +1
        if (k == H_lim(irank+1)) go to 112
        if (k .lt. H_lim(irank)) CYCLE
        l = l+1
        i_idx(l) = i
        j_idx(l) = j
     end do
  end do
 112 CONTINUE
  !do i = 1,nelc(irank)
  !   print*, irank, i_idx(i), j_idx(i)
  !end do

  allocate(H(n_snap,n_snap+1))
  H = 0.0D0

  allocate(xi(nx1,ny1,nz1), xj(nx1,ny1,nz1))

  do i = 1,nelc(irank)
     print*, irank, i_idx(i), j_idx(i)

     if (i_idx(i) == j_idx(i)) then
        ! Open the binary files
        l = start_file +step_file*(i_idx(i)-1)
        write(fname, '(''vol_data.1_'', I8.8)') l
        filename1 = "volume_data/"//trim(fname)
        print*, irank, i_idx(i)-1, filename1
        call BINARY_FILE_OPEN(iunit1, trim(filename1), "r", ierr)

        do j = 1, 5
           call BINARY_FILE_READ(iunit1, idum, 1, kind(idum), ierr)
        end do
        do j = 1, nx1+ny1+nz1
           call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
        end do
        do j = 1,nvarsT
           call BINARY_FILE_READ(iunit1, cdum, str_short, kind(cdum), ierr)
        end do
        do j = 1, 2
           call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
        end do
        do ll = 1, var_list(1)-1
           call BINARY_FILE_READ(iunit1,xi,l_var,kind(xi),ierr)
        end do

        do l = 1,nvars1
           xi(:,:,:) = 0.0_WP
           call BINARY_FILE_READ(iunit1,xi,l_var,kind(xi),ierr)
           if (names(var_list(l)) == 'P') xi(:,:,:) = xi(:,:,:) +P_bulk
           !print*, xi(1,1,1), xi(nx1,ny1,nz1)
           do k = minnz,maxnz
              do j = minny,maxny
                 do m = minnx,maxnx
                    H(i_idx(i),j_idx(i)) = H(i_idx(i),j_idx(i))+(xi(m,j,k)/U0(l))**2.0_WP
                 end do
              end do
           end do
           if ( l == nvars1)  GO TO 212
           ll = var_list(l+1) - var_list(l) -1
           do j = 1,ll
              call BINARY_FILE_READ(iunit1,xi,l_var,kind(xi),ierr)
           end do
        end do
 212 CONTINUE

        call BINARY_FILE_CLOSE(iunit1, ierr)

!     else if (i_idx(i) == i_idx(i-1)) then
!        ! Open the binary files
!        l = start_file +step_file*(j_idx(i)-1)
!
!        write(fname, '(''vol_data.1_'', I8.8)') l
!        filename1 = "volume_data/"//trim(fname)
!        print*, irank, j_idx(i)-1, filename1
!        call BINARY_FILE_OPEN(iunit2, trim(filename1), "r", ierr)
!
!        do j = 1, 5
!           call BINARY_FILE_READ(iunit2, idum, 1, kind(idum), ierr)
!        end do
!        do j = 1, nx1+ny1+nz1
!           call BINARY_FILE_READ(iunit2, ddum, 1, kind(ddum), ierr)
!        end do
!        do j = 1,nvarsT
!           call BINARY_FILE_READ(iunit2, cdum, str_short, kind(cdum), ierr)
!        end do
!        do j = 1, 2
!           call BINARY_FILE_READ(iunit2, ddum, 1, kind(ddum), ierr)
!        end do
!        do ll = 1, var_list(1)-1
!           call BINARY_FILE_READ(iunit2,xj,l_var,kind(xj),ierr)
!        end do
!
!        do l = 1,nvars1
!           xj(:,:,:) = 0.0_WP
!           call BINARY_FILE_READ(iunit2,xj,l_var,kind(xj),ierr)
!           do k = minnz,maxnz
!              do j = minny,maxny
!                 do m = minnx,maxnx
!                    H(i_idx(i),j_idx(i)) = H(i_idx(i),j_idx(i))+(xi(m,j,k)/U0(l))*(xj(m,j,k)/U0(l))
!                 end do
!              end do
!           end do
!           if ( l == nvars1)  GO TO 211
!           ll = var_list(l+1) - var_list(l) -1
!           do j = 1,ll
!              call BINARY_FILE_READ(iunit2,xj,l_var,kind(xj),ierr)
!           end do
!        end do
! 211 CONTINUE
!
!        call BINARY_FILE_CLOSE(iunit2, ierr)
     else
        ! Open the binary files
        l = start_file +step_file*(i_idx(i)-1)
        write(fname, '(''vol_data.1_'', I8.8)') l
        filename1 = "volume_data/"//trim(fname)
        print*, irank, i_idx(i)-1, filename1
        call BINARY_FILE_OPEN(iunit1, trim(filename1), "r", ierr)

        l = start_file +step_file*(j_idx(i)-1)
        write(fname, '(''vol_data.1_'', I8.8)') l
        filename1 = "volume_data/"//trim(fname)
        !print*, irank, j_idx(i)-1, filename1
        call BINARY_FILE_OPEN(iunit2, trim(filename1), "r", ierr)

        do j = 1, 5
           call BINARY_FILE_READ(iunit1, idum, 1, kind(idum), ierr)
           call BINARY_FILE_READ(iunit2, idum, 1, kind(idum), ierr)
        end do
        do j = 1, nx1+ny1+nz1
           call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
           call BINARY_FILE_READ(iunit2, ddum, 1, kind(ddum), ierr)
        end do
        do j = 1,nvarsT
           call BINARY_FILE_READ(iunit1, cdum, str_short, kind(cdum), ierr)
           call BINARY_FILE_READ(iunit2, cdum, str_short, kind(cdum), ierr)
        end do
        do j = 1, 2
           call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
           call BINARY_FILE_READ(iunit2, ddum, 1, kind(ddum), ierr)
        end do
        do ll = 1, var_list(1)-1
           call BINARY_FILE_READ(iunit1,xi,l_var,kind(xi),ierr)
           call BINARY_FILE_READ(iunit2,xj,l_var,kind(xj),ierr)
        end do

        do l = 1,nvars1
           xi(:,:,:) = 0.0_WP; xj(:,:,:) = 0.0_WP
           call BINARY_FILE_READ(iunit1,xi,l_var,kind(xi),ierr)
           call BINARY_FILE_READ(iunit2,xj,l_var,kind(xj),ierr)
           if (names(var_list(l)) == 'P') then
              xi(:,:,:) = xi(:,:,:) +P_bulk
              xj(:,:,:) = xj(:,:,:) +P_bulk
           end if
           !print*, xi(1,1,1), xi(nx1,ny1,nz1)
           do k = minnz,maxnz
              do j = minny,maxny
                 do m = minnx,maxnx
                    H(i_idx(i),j_idx(i)) = H(i_idx(i),j_idx(i)) +(xi(m,j,k)/U0(l))*(xj(m,j,k)/U0(l))
                 end do
              end do
           end do
           if ( l == nvars1)  GO TO 213
           ll = var_list(l+1) - var_list(l) -1
           do j = 1,ll
              call BINARY_FILE_READ(iunit1,xi,l_var,kind(xi),ierr)
              call BINARY_FILE_READ(iunit2,xj,l_var,kind(xj),ierr)
           end do
        end do
 213 CONTINUE

        call BINARY_FILE_CLOSE(iunit1, ierr)
        call BINARY_FILE_CLOSE(iunit2, ierr)
     end if

  end do
  deallocate(xi, xj)

  allocate(H_recv(n_snap,n_snap+1))  
  call MPI_ALLREDUCE(H, H_recv, n_snap*(n_snap+1), MPI_REAL_WP, MPI_SUM, MPI_COMM_WORLD, ierr)

  write( filename1, '(''/H_data_'',I4.4,''_'', I3.3, ''_head.dat'')') n_snap, step_file
  filename1 = trim(folder)//trim(filename1)
  if(irank == 1 ) print*, irank, filename1
  call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filename1), MPI_MODE_WRONLY+MPI_MODE_CREATE, mpi_info, iunit2, ierr)

  if(irank == 1 ) then
     call MPI_FILE_WRITE(iunit2, n_snap, 1, MPI_INTEGER, status, ierr)
     !if(irank == 1 ) print*, irank, 'writing'
     do i = 1,n_snap
        do j = i,n_snap+1
           !print *, i, j, H_recv(i,j)
           !write(678,*) i, j, H_recv(i,j)
           call MPI_FILE_WRITE(iunit2, H_recv(i,j), 1, MPI_REAL_WP, status, ierr)
        end do
     end do
  end if

  call MPI_FILE_CLOSE(iunit2, ierr)

  deallocate(i_idx, j_idx)
  deallocate(H)
  deallocate(H_recv)
 111 CONTINUE
  deallocate(H_lim)

 110 CONTINUE
  print *, irank, 'Completed dmd_3D'

  ! Finalize the parallel environment
  call parallel_final()

end program dmd_3D_p1
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
