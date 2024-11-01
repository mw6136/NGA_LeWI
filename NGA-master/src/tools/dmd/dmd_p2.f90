! ============================================================ !
!                         dmd.f90                              !
!                                                              !
! Author: Temistocle Grenga                                    !
! Date:   September 27, 2016                                   !
! ============================================================ !

program dmd_3D_p2
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
  character(len=str_medium) :: filename1
  integer :: i, j, k, l, ll, o, pr
  integer :: iunit1, iunit2, ierr, iunit5, iunit6
  real(WP), dimension(:,:), pointer :: H, U1, H1, M, Lambda, VL, VR, VI
  real(WP), dimension(:,:), pointer :: dummyarr, dummyarr2
  real(WP), dimension(:), pointer :: H2, S
  integer :: n_snap, start_file, step_file, name_snap
  integer, dimension(MPI_STATUS_SIZE) :: status
  CHARACTER*1 :: JOBZ, UPLO, JOBVL, JOBVR
  INTEGER :: LWORK, LIWORK, INFO, nb
  real(WP), dimension(:), pointer :: WORK, dummy, dummy2
  integer, dimension(:), pointer :: IWORK, IPIV1
  integer, dimension(:,:), pointer :: IWKN, IPIV2
  integer :: nx1, ny1, nz1, nvars, nvars1, ntime1, np
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nd1_MOK, nd2_MOK, nd3_MOK, ii_MOK, ll_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK, str_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK, NTIME_MOK, idx_MOK, snap_MOK
  CHARACTER*100:: finame
  LOGICAL :: cc, lprint
  LOGICAL, dimension(:), pointer :: log_compl
  COMPLEX(WP), dimension(:,:), pointer :: VC, VCstar, dummyarrC, dummyarrC2, T, test
  COMPLEX(WP), dimension(:), pointer :: LambdaC, dummy3, AP, WORKc1, D, SS
  integer :: w_snap, w_var, size, ibuffer, ipart, idum
  character(len=str_short), dimension(:), pointer :: names
  character(len=str_medium) :: file, folder, filename_base, fname
  integer :: ifile
  logical :: file_is_there
  real(WP), dimension(:), pointer :: x1, y1, z1
  real(SP), dimension(:), pointer :: xm1, ym1, zm1
  real(SP) :: phi4
  character(len=80) :: buffer
  integer, dimension(:,:,:), pointer :: iblank
  real(WP), dimension(:,:,:), pointer :: data8
  real(SP), dimension(:,:,:), pointer :: data4
  real(WP) :: dt1, time1, ddum, rsnap
  real(WP) :: min_x, max_x, min_y, max_y, min_z, max_z
  integer :: minnx, maxnx, minny, maxny, minnz, maxnz
  integer, dimension(:), pointer :: var_list
  logical :: subdomain, allmode
  integer :: nmode
  if (irank == 1) print *, 'Starting dmd_3D Part 2'

  !Initialize MPI enironment
  call parallel_init()

  start_file = 8830 !6465 !8830 !15286 !21500
  step_file = 2
  name_snap = 200
  nvars = 1
  allmode = .false.
  nmode = 10

  subdomain = .false.
  min_x = 0.00432_WP*7.0_WP
  max_x = 0.00432_WP*8.0_WP
  min_y = 0.00432_WP*-2.00_WP
  max_y = 0.00432_WP*2.00_WP
  min_z = 0.00432_WP*-0.5_WP
  max_z = 0.00432_WP*0.5_WP

  allocate(var_list(nvars))
  var_list = (/ 1 /)
  !var_list = (/ 1, 2, 3, 4, 10, 12, 13, 14, 15 /)
  rsnap = REAL(name_snap,WP)

  write(folder, '(''U_'',I4.4,''_'',I2.2,''_tail/'')') name_snap, step_file
  write( filename1, '(''H_data_'',I4.4,''_'', I3.3, ''_tail.dat'')') name_snap, step_file
  filename1 = trim(folder)//trim(filename1)
  if (irank == 1) print*, irank, filename1
  call MPI_FILE_OPEN(MPI_COMM_WORLD, filename1, MPI_MODE_RDONLY, mpi_info, iunit2, ierr)
  call MPI_FILE_READ(iunit2, n_snap, 1, MPI_INTEGER, status, ierr)
  if (irank == 1) print*, 'snapshots', n_snap
  if (name_snap /= n_snap) then
     if (irank == 1) print*, 'WARNING: name_snap /= n_snap'
     if (irank == 1) print*, 'name snap ', name_snap
     if (irank == 1) print*, 'n_snap ', n_snap
     !if (irank == 1) print*, 'EXIT'
     !go to 110
     n_snap = name_snap
  end if

  allocate(H(n_snap,n_snap), H2(n_snap))
  H = 0.0D0; H2 = 0.0D0
  do i = 1,n_snap
     do j = i,n_snap
        call MPI_FILE_READ(iunit2, H(i,j), 1, MPI_REAL_WP, status, ierr)
     end do
     call MPI_FILE_READ(iunit2, H2(i), 1, MPI_REAL_WP, status, ierr)
  end do
  call MPI_FILE_CLOSE(iunit2, ierr)

  !if (irank == 1) then
  !   do i = 1,n_snap
  !      do j = i,n_snap
  !         write(456,'(2E24.12E4)') H(i,j)
  !      end do
  !   end do
  !end if

  allocate(U1(n_snap,n_snap), dummy(n_snap), dummyarr(n_snap,n_snap), dummyarr2(n_snap,n_snap))
  k = 0; l = 0
  do i = 1,n_snap
     do j = i,n_snap
        dummyarr(i, j) = H(i, j)
     end do
     dummy(i) = REAL(i,WP)
  end do

  JOBZ = 'V'
  UPLO = 'U'
  LWORK = (n_snap+2)*n_snap !for DSYEV
  allocate (WORK(LWORK), S(n_snap))
  call DSYEV( JOBZ, UPLO, n_snap, dummyarr, n_snap, S, WORK, LWORK, INFO )
  deallocate(WORK)
  if (irank == 1) then
     print*, 'H eigenvalues'
     do i = 1,n_snap
        print*,i, S(i)
     end do
  endif

  call DSORT (S, dummy, n_snap, -2)
  do i = 1,n_snap
     U1(:,i) = dummyarr(:,INT(dummy(i)))
  end do
  deallocate(dummy)

  !allocate(SS(n_snap))
  do i = 1,n_snap
     do j = i+1,n_snap
        H(j, i) = H(i, j)  !Fill symmetric matrix
     end do
     !SS(i) = DSQRT(S(i))
     !S(i) = DSQRT(DABS(S(i)))    !sqrt of eigenvalues
     S(i) = DSQRT(S(i))
  end do
  if (irank == 1) then
     !print*, 'sqrt complex H eigenvalues'
     !do i = 1,n_snap
     !   print*,i, SS(i)
     !end do
     !print*, 'H eigenvector'
     !do i = 1,n_snap
     !   print*,i, U1(i,:)
     !end do
  endif

  !call parallel_final()
  !print *, irank, 'Stop'
  !STOP


  allocate(H1(n_snap,n_snap))
  do i = 1,n_snap
     do j = 1,n_snap-1
        H1(i,j) = H(i,j+1)
     end do
     H1(i,n_snap) = H2(i)
     S(i) = 1.0D0/S(i)
  end do
  deallocate(H2)

  dummyarr = 0.0D0
  do i = 1,n_snap
     dummyarr(:,i) = U1(:,i)*S(i)
  end do
  !if (irank == 1) then
  !   print*, 'US'
  !   do i = 1,n_snap
  !      print*, dummyarr(i,:)
  !   end do
  !endif

  do i = 1,n_snap
     do j = 1,n_snap
        dummyarr2(i,j) = 0.0D0
        do k = 1,n_snap
           dummyarr2(i,j) = dummyarr2(i,j)+ H1(i,k)*dummyarr(k,j)
        end do
     end do
  end do
  !if (irank == 1) then
  !   print*, 'HUS'
  !   do i = 1,n_snap
  !      print*, dummyarr2(i,:)
  !   end do
  !endif

  do i = 1,n_snap
     do j = 1,n_snap
        dummyarr(i,j) = 0.0D0
        do k = 1,n_snap 
           dummyarr(i,j) = dummyarr(i,j)+ U1(k,i)*dummyarr2(k,j)
        end do
     end do
  end do
  !if (irank == 1) then
  !   print*, 'UHUS'
  !   do i = 1,n_snap
  !      print*, dummyarr(i,:)
  !   end do
  !endif
  allocate(M(n_snap,n_snap))
  M = 0.0D0
  do i = 1,n_snap
     M(i,:) = S(i)*dummyarr(i,:)
  end do
  !if (irank == 1) then
  !   print*, 'M'
  !   do i = 1,n_snap
  !      print*, M(i,:)
  !   end do
  !end if

  JOBVL = 'N'
  JOBVR = 'V'
  LWORK = 4*n_snap
  allocate(Lambda(n_snap,2), VL(n_snap,n_snap), VR(n_snap,n_snap), WORK(LWORK) )
  call DGEEV( JOBVL, JOBVR, n_snap, M, n_snap, Lambda(:,1),  Lambda(:,2), VL, n_snap, VR, n_snap, WORK, LWORK, INFO )
  deallocate(WORK)

  allocate(log_compl(n_snap), LambdaC(n_snap), VC(n_snap,n_snap), VCstar(n_snap,n_snap))
  log_compl = .false.
  do i = 1,n_snap
     if (DABS(Lambda(i,2)) .gt. 1.0D-13) then
        if (irank == 1) print*,i, 'Complex eigenvalue'
        log_compl(i) = .true.
     end if
  end do

  cc = .false.
  do i = 1,n_snap
     if (log_compl(i)) then
        if (cc) then
           cc = .false.
        else
           cc = .true.
           LambdaC(i) = cmplx(Lambda(i,1), Lambda(i,2))
           LambdaC(i+1) = cmplx(Lambda(i+1,1), Lambda(i+1,2))
           do j = 1,n_snap
              VC(j,i) = cmplx(VR(j,i), VR(j,i+1))
              VC(j,i+1) = cmplx(VR(j,i), -VR(j,i+1))
              VCstar(i,j) = cmplx(VR(j,i), -VR(j,i+1))
              VCstar(i+1,j) = cmplx(VR(j,i), VR(j,i+1))
           end do
        end if
     else
        LambdaC(i) = cmplx(Lambda(i,1), 0.0D0)
        do j = 1,n_snap
           VC(j,i) = cmplx(VR(j,i), 0.0D0)
           VCstar(i,j) = cmplx(VR(j,i), 0.0D0)
        end do
     end if
  end do

  allocate(dummy(n_snap), dummy2(n_snap))
  do i = 1,n_snap
     dummy(i) = 0.0D0
     do k = 1,n_snap
        dummy(i) = dummy(i)+ U1(k,i)*H(k,1)
     end do
  end do
  do i = 1,n_snap
     dummy(i) = dummy(i)*S(i)
  end do
  allocate(dummy3(n_snap))
  do i = 1,n_snap
     dummy3(i) = 0.0D0
     do k = 1,n_snap
        dummy3(i) = dummy3(i)+ VCstar(i,k)*dummy(k)
     end do
  end do

  allocate(dummyarrC(n_snap,n_snap))
  do i = 1,n_snap
     do j = 1,n_snap
        dummyarrC(i,j) = 0.0D0
        do k = 1,n_snap
           dummyarrC(i,j) = dummyarrC(i,j)+ VCstar(i,k)*VC(k,j)
        end do
     end do
  end do

  ! Complex general matrix
  LWORK = n_snap*n_snap
  allocate(IPIV1(n_snap), WORKc1(LWORK))
  call ZGETRF(n_snap, n_snap, dummyarrC, n_snap, IPIV1, INFO)
  call ZGETRI(n_snap, dummyarrC, n_snap, IPIV1, WORKc1, LWORK, INFO)
  deallocate(IPIV1, WORKc1)

  allocate(D(n_snap))
  do i = 1,n_snap
     D(i) = 0.0D0
     do k = 1,n_snap
        D(i) = D(i)+ dummyarrC(i,k)*dummy3(k)
     end do
  end do

  dummyarrC = (0.0D0, 0.0D0)
  do i = 1,n_snap
     dummyarrC(:,i) = VC(:,i)*D(i)
  end do

  allocate(dummyarrC2(n_snap, n_snap))
  dummyarrC2 = (0.0D0, 0.0D0)
  do i = 1,n_snap
     dummyarrC2(i,:) = dummyarrC(i,:)*S(i)
  end do

  allocate(T(n_snap, n_snap))
  do i = 1,n_snap
     do j = 1,n_snap
        T(i,j) = 0.0D0
        do k = 1,n_snap
           T(i,j) = T(i,j)+ U1(i,k)*dummyarrC2(k,j)
        end do
     end do
  end do

  lprint = .true.
  if (lprint) then
  if (irank == 1) then
     print*,'T'
     do i = 1,n_snap
        print*,T(i,:)
     end do
     print*,'Lambda'
     do i = 1,n_snap
        print*,Lambda(i,1), Lambda(i,2)
        write(777,'(2E24.12E4)') Lambda(i,1), Lambda(i,2)
     end do

     fname = 'Tarray'
     filename1 = trim(folder)//trim(fname)
     open(iunit6,file=trim(filename1),form="formatted",iostat=ierr)
     do i = 1,n_snap
        write(iunit6,'(500E24.12E4)') T(i,1:n_snap)
     end do
     close(iclose(iunit6))

     fname = 'Lambda'
     filename1 = trim(folder)//trim(fname)
     open(iunit5,file=trim(filename1),form="formatted",iostat=ierr) 
     do i = 1,n_snap
        write(iunit5,'(2E24.12E4)') Lambda(i,1), Lambda(i,2)
     end do
     close(iclose(iunit5))
  end if
  end if

  deallocate(dummy, dummy2, dummyarr, dummyarr2)
  deallocate(H, H1, U1, S, D)
  deallocate(VC, VCstar, dummyarrC, dummyarrC2, dummy3)

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

  !np = nx1*ny1*nz1
  !if (irank == 1) print*, 'n var', nvars

  ! Resize some integers so MPI can write even the biggest files
  WP_MOK    = int(WP,        MPI_Offset_kind)
  str_MOK   = int(str_short, MPI_Offset_kind)
  NVARS_MOK = int(nvars1,    MPI_Offset_kind)
  NTIME_MOK = int(ntime1,    MPI_Offset_kind)
  nd1_MOK   = int(nx1,       MPI_Offset_kind)
  nd2_MOK   = int(ny1,       MPI_Offset_kind)
  nd3_MOK   = int(nz1,       MPI_Offset_kind)

  if (allmode) then 
     if (nproc /= n_snap*nvars) then
        if (irank == 1) print*, 'nproc .neq. n_snap*nvars1'
           call parallel_final()
           GO TO 110
     end if
     w_snap = MOD(irank,n_snap)
     w_var = CEILING(REAL(irank,WP)/REAL(n_snap,WP))
     if (w_snap == 0) w_snap = n_snap
  else
     w_snap = MOD(irank,nmode)
     if (w_snap == 0) w_snap = nmode
     w_snap = w_snap !+nmode
     w_var = 1
     !w_var = CEILING(REAL(irank,WP)/REAL(100,WP))
  end if

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

  deallocate(x1, y1, z1)

  allocate(names(nvars1))
  do i = 1,nvars1
     call BINARY_FILE_READ(iunit1, names(i), str_short, kind(names), ierr)
  end do
  if (irank == 1) print*, 'vars', names
  call BINARY_FILE_CLOSE(iunit1,ierr)

  fname = trim(folder)//trim(adjustl(names(var_list(w_var)))) // "_"
  write(fname(len_trim(fname)+1:len_trim(fname)+3),'(i3.3)') w_snap
  print*, irank, fname
  call BINARY_FILE_OPEN(ifile,trim(fname),"w",ierr)
!  print*, irank, fname


  buffer = trim(adjustl(names(var_list(w_var)))) // "_"
  write(buffer(len_trim(buffer)+1:len_trim(buffer)+3),'(i3.3)') w_snap
  call BINARY_FILE_WRITE(ifile,buffer,80,kind(buffer),ierr)
  buffer = 'part'
  call BINARY_FILE_WRITE(ifile,buffer,80,kind(buffer),ierr)
  ibuffer = 1
  call BINARY_FILE_WRITE(ifile,ibuffer,1,kind(ibuffer),ierr)
  buffer = 'block'
  call BINARY_FILE_WRITE(ifile,buffer,80,kind(buffer),ierr)

  allocate(data8(nx1,ny1,nz1))
  allocate(data4(nx1,ny1,nz1))
  data4(:,:,:) = 0.0_SP

  do l = 1,n_snap

     i = start_file +(l-1)
     write(fname, '(''_'', I8.8)') i
     fname = trim(filename_base)//trim(fname)
     filename1 = "volume_data/"//trim(fname)

     call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(filename1), MPI_MODE_RDONLY, mpi_info, iunit1, ierr) 
     !call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

     !do i = 1, 5
     !   call BINARY_FILE_READ(iunit1, idum, 1, kind(idum), ierr)
     !end do
     !do i = 1, nx1+ny1+nz1
     !   call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
     !end do

     !do i = 1,nvars1
     !   call BINARY_FILE_READ(iunit1, names(i), str_short, kind(names), ierr)
     !end do

     !do i = 1, 2
     !   call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
     !end do

     !do ll = 1, var_list(w_var)-1
     !   call BINARY_FILE_READ(iunit1,data8,nx1*ny1*nz1,kind(data8),ierr)
     !end do
     data8(:,:,:) = 0.0_WP

     ll_MOK = int(var_list(w_var)-1,MPI_Offset_kind)
     disp = 5*4 + (nd1_MOK+nd2_MOK+nd3_MOK)*WP_MOK &
            + str_MOK*NVARS_MOK &
            + 2*WP_MOK &
            + ll_MOK*nd1_MOK*nd2_MOK*nd3_MOK*WP_MOK
     call MPI_FILE_SET_VIEW(iunit1, disp, MPI_INTEGER, MPI_INTEGER, "native",  mpi_info, ierr)
     call MPI_FILE_READ(iunit1, data8, nx1*ny1*nz1, MPI_REAL_WP, status, ierr)
     !call BINARY_FILE_READ(iunit1,data8,nx1*ny1*nz1,kind(data8),ierr)
     
     data4(1:maxnx-minnx+1,1:maxny-minny+1,1:maxnz-minnz+1) = REAL(data8(minnx:maxnx,minny:maxny,minnz:maxnz),SP)
     !data4(1:maxnx-minnx+1,1:maxny-minny+1,1:maxnz-minnz+1) = data4(1:maxnx-minnx+1,1:maxny-minny+1,1:maxnz-minnz+1)+ REAL(data8(minnx:maxnx,minny:maxny,minnz:maxnz)*T(l,w_snap)*(LambdaC(w_snap)**rsnap),SP)

     call MPI_FILE_CLOSE(iunit1, ierr)
     !call BINARY_FILE_CLOSE(iunit1,ierr)
  end do
  call BINARY_FILE_WRITE(ifile,data4,(maxnx-minnx+1)*(maxny-minny+1)*(maxnz-minnz+1),kind(data4),ierr)

  call BINARY_FILE_CLOSE(ifile,ierr)

!  if (irank == 1) then
!     fname = 'Lambda'
!     filename1 = trim(folder)//trim(fname)
!     open(iunit5,file=trim(filename1),form="formatted",iostat=ierr,status="new")
!     do i = 1,n_snap
!        write(iunit5,'(2E24.12E4)') Lambda(i,1), Lambda(i,2)
!     end do
!     close(iclose(iunit5))
!  end if


 110 continue

  ! Finalize the parallel environment
  call parallel_final()
  print *, irank, 'Completed dmd_3D'

end program dmd_3D_p2
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

