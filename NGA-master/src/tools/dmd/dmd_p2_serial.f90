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
  !use dump_volume
  implicit none
  character(len=str_long) :: filename1
  integer :: i, j, k, l, ll, o, pr
  integer :: iunit1, iunit2, ierr
  real(WP), dimension(:,:), pointer :: H, U1, H1, M, Lambda, VL, VR, VI
  real(WP), dimension(:,:), pointer :: dummyarr, dummyarr2
  real(WP), dimension(:), pointer :: H2, S
  integer :: n_snap
  integer, dimension(MPI_STATUS_SIZE) :: status
  CHARACTER*1 :: JOBZ, UPLO, JOBVL, JOBVR
  INTEGER :: LWORK, LIWORK, INFO, nb
  real(WP), dimension(:), pointer :: WORK, dummy, dummy2
  integer, dimension(:), pointer :: IWORK, IPIV1
  integer, dimension(:,:), pointer :: IWKN, IPIV2
  integer :: nx1, ny1, nz1, nvars1, ntime1, np, jump_var
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nd1_MOK, nd2_MOK, nd3_MOK, ii_MOK, jj_MOK, kk_MOK, ll_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK, str_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK, NTIME_MOK, idx_MOK, snap_MOK
  real(WP) :: xx(3)
  real(WP), dimension(:), pointer :: vel
  CHARACTER*100:: finame
  complex(WP), dimension(:,:), pointer :: phi
  LOGICAL :: cc, debug1, debug2, debug3
  LOGICAL, dimension(:), pointer :: log_compl
  COMPLEX(WP), dimension(:,:), pointer :: VC, VCstar, dummyarrC, dummyarrC2, T, test
  COMPLEX(WP), dimension(:), pointer :: LambdaC, dummy3, AP, WORKc1, D

  print *, 'Starting dmd_3D Part 2'

  !Initializa MPI enironment
  call parallel_init()

        debug1 = .false.
        debug2 = .false.
        debug3 = .false.

        if (debug3) GO TO 202
        if (debug2) GO TO 201

        if (debug1) then
                print*, 'Running TOY-Problem'
                n_snap = 3
                allocate(H(n_snap,n_snap))
                H(1,:) = (/1, 4, 3/)
                H(2,:) = (/4, 1, 0/)
                H(3,:) = (/3, 0, 1/)

                !print*, H(:,1)
                !print*, H(:,2)
                !print*, H(:,3)
        else

  filename1 = "./H_data.42.40"
  call MPI_FILE_OPEN(MPI_COMM_WORLD, filename1, MPI_MODE_RDONLY, mpi_info, iunit2, ierr)
  !call BINARY_FILE_OPEN(iunit2, trim(filename1), "r", ierr)
  call MPI_FILE_READ(iunit2, n_snap, 1, MPI_INTEGER, status, ierr)
  !call BINARY_FILE_READ(iunit2, n_snap, 1, kind(n_snap), ierr)
  print*, n_snap
  allocate(H(n_snap,n_snap), H2(n_snap))
  H = 0.0D0; H2 = 0.0D0
  do i = 1,n_snap
     do j = i,n_snap
        call MPI_FILE_READ(iunit2, H(i,j), 1, MPI_REAL_WP, status, ierr)
        !call BINARY_FILE_READ(iunit2, H(i,j), 1, kind(H), ierr)
     end do
     call MPI_FILE_READ(iunit2, H2(i), 1, MPI_REAL_WP, status, ierr)
     !call BINARY_FILE_READ(iunit2, H2(i), 1, kind(H2), ierr)
  end do
  
  !do i = 1,n_snap
  !   do j = i,n_snap
  !      print*, irank, 'H(',i,',',j,') = ',H(i,j)
  !   end do
  !end do
  !do i = 1,n_snap
  !   print*, irank, 'H2(',i,') = ',H2(i)
  !end do

  call MPI_FILE_CLOSE(iunit2, ierr)
  !call BINARY_FILE_CLOSE(iunit2, ierr)

        end if !debug1

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
  !print*,'DSYEV'
  call DSYEV( JOBZ, UPLO, n_snap, dummyarr, n_snap, S, WORK, LWORK, INFO )
  !print*,'Eigenvectors'
  !do i = 1,n_snap
  !   print*, i,'--',S(i)
  !   print*, dummyarr(:,i)
  !end do
  deallocate(WORK)

  call DSORT (S, dummy, n_snap, -2)
  do i = 1,n_snap
     U1(:,i) = dummyarr(:,INT(dummy(i)))
  end do
  deallocate(dummy)

  !print*,'Ordered Eigenvalues and Eigenvectors'
  !do i = 1,n_snap
  !   print*, i,'--',S(i)
  !   print*, U1(:,i)
  !end do

  do i = 1,n_snap
     do j = i+1,n_snap
        H(j, i) = H(i, j)  !Fill symmetric matrix
     end do
     S(i) = DSQRT(S(i))    !sqrt of eigenvalues
  end do

  !*** Verify that HU=US ***
  !dummyarr = 0.0D0
  !dummyarr2 = 0.0D0
  !do i = 1,n_snap
  !   do j = 1,n_snap
  !      do k = 1,n_snap
  !         dummyarr(i,j) = dummyarr(i,j)+ H(i,k)*U1(k,j)
  !      end do
  !   end do
  !end do
  !do i = 1,n_snap
  !   dummyarr2(:,i) = U1(:,i)*S(i)  !S has to be the eigenvalues, not their sqrt
  !end do
  !print*,'HU'
  !do i = 1,n_snap
  !   print*, dummyarr(i,:)
  !end do
  !print*,'US'
  !do i = 1,n_snap
  !   print*, dummyarr2(i,:)
  !end do

  allocate(H1(n_snap,n_snap))
  do i = 1,n_snap
     do j = 1,n_snap-1
        H1(i,j) = H(i,j+1)
     end do
     H1(i,n_snap) = H2(i)
     S(i) = 1.0D0/S(i)
  end do
  deallocate(H2)

 201 CONTINUE !debug2

        if(debug2) then
                print*, 'Running TOY-Problem - 2'
                n_snap = 3
                allocate(dummyarr(n_snap,n_snap), dummyarr2(n_snap,n_snap))
                allocate(H1(n_snap,n_snap), U1(n_snap,n_snap), S(n_snap))
                U1(1,:) = (/1, 4, 3/)
                U1(2,:) = (/4, 1, 0/)
                U1(3,:) = (/3, 0, 1/)
                H1(1,:) = (/1, 2, 3/)
                H1(2,:) = (/1, 1, 1/)
                H1(3,:) = (/3, 4, 1/)
                S(:) = (/2, 2, 2/)
        end if

  dummyarr = 0.0D0
  do i = 1,n_snap
     dummyarr(:,i) = U1(:,i)*S(i)
  end do
  do i = 1,n_snap
     do j = 1,n_snap
        dummyarr2(i,j) = 0.0D0
        do k = 1,n_snap
           dummyarr2(i,j) = dummyarr2(i,j)+ H1(i,k)*dummyarr(k,j)
        end do
     end do
  end do
  do i = 1,n_snap
     do j = 1,n_snap
        dummyarr(i,j) = 0.0D0
        do k = 1,n_snap 
           dummyarr(i,j) = dummyarr(i,j)+ U1(k,i)*dummyarr2(k,j)
        end do
     end do
  end do
  allocate(M(n_snap,n_snap))
  M = 0.0D0
  do i = 1,n_snap
     M(i,:) = S(i)*dummyarr(i,:)
  end do
  !print*,'M'
  !do i = 1,n_snap
  !   print*, M(i,:)
  !end do

        if(debug2) then
                print*, 'Running TOY-Problem - 2'
                M(1,:) = (/1, 4, 3/)
                M(2,:) = (/4, 1, 0/)
                M(3,:) = (/3, 0, 1/)
                print*,'M'
                do i = 1,n_snap
                        print*, M(i,:)
                end do
        end if

 202 CONTINUE
        if(debug3) then
                print*, 'Running TOY-Problem - 3'
                n_snap = 3
                allocate(dummyarr(n_snap,n_snap), dummyarr2(n_snap,n_snap))
                allocate(H(n_snap,n_snap), U1(n_snap,n_snap), S(n_snap), M(n_snap,n_snap))
                U1(1,:) = (/1, 4, 3/)
                U1(2,:) = (/4, 1, 0/)
                U1(3,:) = (/3, 0, 1/)
                H(1,:) = (/1, 2, 3/)
                H(2,:) = (/1, 1, 1/)
                H(3,:) = (/3, 4, 1/)
                S(:) = (/2, 2, 2/)
                M(1,:) = (/-2, 2, -1/)
                M(2,:) = (/-3, 1, 5/)
                M(3,:) = (/-1, 2, 0/)
        end if

  JOBVL = 'N'
  JOBVR = 'V'
  LWORK = 4*n_snap
  allocate(Lambda(n_snap,2), VL(n_snap,n_snap), VR(n_snap,n_snap), WORK(LWORK) )
  call DGEEV( JOBVL, JOBVR, n_snap, M, n_snap, Lambda(:,1),  Lambda(:,2), VL, n_snap, VR, n_snap, WORK, LWORK, INFO )
  deallocate(WORK)

  !print*,'Lambda'
  !do i = 1,n_snap
  !   print*,Lambda(i,1),  Lambda(i,2)
  !end do
  !print*,'Eigenvectors'
  !do i = 1,n_snap
  !   print*, i,'--',Lambda(i,1),  Lambda(i,2)
  !   print*, VR(:,i)
  !end do
  !*** Verify that MV=VL ***
  !Need to redefine M
  !dummyarr = 0.0D0
  !dummyarr2 = 0.0D0
  !do i = 1,n_snap
  !   do j = 1,n_snap
  !      do k = 1,n_snap
  !         dummyarr(i,j) = dummyarr(i,j)+ M(i,k)*VR(k,j)
  !      end do
  !   end do
  !end do
  !do i = 1,n_snap
  !   do j = 1,n_snap
  !      dummyarr2(j,i) = VR(j,i)*Lambda(i,1)
  !   end do
  !end do
  !print*,'MV'
  !do i = 1,n_snap
  !   print*, dummyarr(i,:)
  !end do
  !print*,'VL'
  !do i = 1,n_snap
  !   print*, dummyarr2(i,:)
  !end do

  allocate(log_compl(n_snap), LambdaC(n_snap), VC(n_snap,n_snap), VCstar(n_snap,n_snap))
  log_compl = .false.
  do i = 1,n_snap
     if (DABS(Lambda(i,2)) .gt. 1.0D-13) then
        print*,i, 'Complex eigenvalue'
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

  !print*, 'VC'
  !do i = 1,n_snap
  !   print*, i, VC(:,i)
  !end do
  !print*, 'VCstar'
  !do i = 1,n_snap
  !   print*, i, VCstar(i,:)
  !end do

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

  !print*, 'VCstarSUH'
  !print*, dummy3(:)

  allocate(dummyarrC(n_snap,n_snap))
  do i = 1,n_snap
     do j = 1,n_snap
        dummyarrC(i,j) = 0.0D0
        do k = 1,n_snap
           dummyarrC(i,j) = dummyarrC(i,j)+ VCstar(i,k)*VC(k,j)
        end do
     end do
  end do

  !allocate(test(n_snap, n_snap))
  !test = dummyarrC
  !print*, 'VV'
  !do i = 1,n_snap
  !   print*, dummyarrC(i,:)
  !end do

  ! Complex Hermitian Matrix - Bunch-Kaufman
  !UPLO = 'U'
  !nb = ILAENV( 1, 'CHETRF', UPLO, n_snap, n_snap, n_snap, n_snap )
  !print*, 'nb', nb
  !LWORK = n_snap*n_snap
  !allocate(IPIV1(n_snap), WORKc1(LWORK))
  !call CHETRF( UPLO, n_snap, dummyarrC, n_snap, IPIV1, WORKc1, LWORK, INFO )
  !deallocate(WORKc1)
  !allocate(WORKc1(n_snap))
  !call CHETRI( UPLO, n_snap, dummyarrC, n_snap, IPIV1, WORKc1, INFO )
  !deallocate(IPIV1, WORKc1)

  ! Complex Hermitian Matrix - Bunch-Kaufman
  !UPLO = 'U'
  !allocate(IPIV1(n_snap), WORKc1(n_snap), AP(n_snap*(n_snap+1)/2))
  !do i = 1,n_snap
  !   do j = i,n_snap !Upper triangle
  !      AP(i + (j-1)*j/2) = dummyarrC(i,j)
  !   end do
  !end do
  !call CHPTRF( UPLO, n_snap, AP, IPIV1, INFO )
  !call CHPTRI( UPLO, n_snap, AP, IPIV1, WORKc1, INFO )
  !do i = 1,n_snap
  !   do j = i,n_snap
  !      dummyarrC(i,j) = AP(i + (j-1)*j/2)
  !   end do
  !end do
  !deallocate(IPIV1, WORKc1, AP)

  ! Complex general matrix
  LWORK = n_snap*n_snap
  allocate(IPIV1(n_snap), WORKc1(LWORK))
  call ZGETRF(n_snap, n_snap, dummyarrC, n_snap, IPIV1, INFO)
  call ZGETRI(n_snap, dummyarrC, n_snap, IPIV1, WORKc1, LWORK, INFO)
  deallocate(IPIV1, WORKc1)

  ! Complex general matrix
  !LWORK = n_snap*n_snap
  !allocate(IPIV2(n_snap, n_snap), WORKc1(LWORK))
  !call CGETRF(n_snap, n_snap, dummyarrC, n_snap, IPIV2, INFO)
  !call CGETRI(n_snap, dummyarrC, n_snap, IPIV2, WORKc1, LWORK, INFO)
  !deallocate(IPIV2, WORKc1)

  !Real general matrix
  !LWORK = n_snap*n_snap
  !allocate(IWKN(n_snap, n_snap), WORK(LWORK))
  !call DGETRF(n_snap, n_snap, dummyarr, n_snap, IWKN, INFO)
  !call DGETRI(n_snap, dummyarr, n_snap, IWKN, WORK, LWORK, INFO)
  !deallocate(IWKN, WORK)

  !print*, 'Info', INFO
  !print*, 'Inverse VV'
  !do i = 1,n_snap
  !   print*, dummyarrC(i,:)
  !end do

  !check Inverse
  !allocate(dummyarrC2(n_snap, n_snap))
  !do i = 1,n_snap
  !   do j = 1,n_snap
  !      dummyarrC2(i,j) = 0.0D0
  !      do k = 1,n_snap
  !         dummyarrC2(i,j) = dummyarrC2(i,j) + dummyarrC(i,k)*test(k,i)
  !      end do
  !   end do
  !end do
  !print*, 'Test Inverse VV'
  !do i = 1,n_snap
  !   print*, dummyarrC2(i,:)
  !end do

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

  print*,'T'
  do i = 1,n_snap
     print*,T(i,:)
  end do
  print*,'Lambda'
  do i = 1,n_snap
       print*,Lambda(i,1), Lambda(i,2)
  end do

  open(9000,file='Lambda.dat',status="unknown")
  do i = 1,n_snap
     write(9000,'(2E24.12E4)') Lambda(i,1), Lambda(i,2)
  end do
  
  deallocate(dummy, dummy2, dummyarr, dummyarr2)
  deallocate(H, H1, U1, S, D)
  deallocate(VC, VCstar, dummyarrC, dummyarrC2, dummy3)

  filename1 = "./vol_data.42"
  call MPI_FILE_OPEN(MPI_COMM_WORLD, filename1, MPI_MODE_RDONLY, mpi_info, iunit1, ierr)

  call MPI_FILE_READ(iunit1, ntime1, 1, MPI_INTEGER, status, ierr)
  call MPI_FILE_READ(iunit1, nx1, 1, MPI_INTEGER, status, ierr)
  call MPI_FILE_READ(iunit1, ny1, 1, MPI_INTEGER, status, ierr)
  call MPI_FILE_READ(iunit1, nz1, 1, MPI_INTEGER, status, ierr)
  call MPI_FILE_READ(iunit1, nvars1, 1, MPI_INTEGER, status, ierr)
  np = nx1*ny1*nz1
  !nvars1 = 1
  jump_var = 0
  print*, nvars1

  ! Resize some integers so MPI can write even the biggest files
  WP_MOK    = int(WP,        MPI_Offset_kind)
  str_MOK   = int(str_short, MPI_Offset_kind)
  NVARS_MOK = int(nvars1,    MPI_Offset_kind)
  NTIME_MOK = int(ntime1,    MPI_Offset_kind)
  nd1_MOK   = int(nx1,       MPI_Offset_kind)
  nd2_MOK   = int(ny1,       MPI_Offset_kind)
  nd3_MOK   = int(nz1,       MPI_Offset_kind)

!  open(9001,file='Umode.dat',status="unknown")
!  open(9002,file='Vmode.dat',status="unknown")
!  open(9003,file='Wmode.dat',status="unknown")

  do i = 1,nvars1
     write( finame, '(''U'', I2.2,''mode.dat'')') i+jump_var
     open(9000+i, file=finame, status="unknown")
  end do
  open(8001,file='U.dat',status="unknown")

  allocate(phi(3,n_snap), vel(nvars1))
  do i = 1,nx1
  ii_MOK = int(i,MPI_Offset_kind)
  do j = 1,ny1
  jj_MOK = int(j,MPI_Offset_kind)
  do k = 1,nz1
     kk_MOK = int(k,MPI_Offset_kind)
     !idx_MOK = int((i-1)*ny1*nz1+(j-1)*nz1+k,MPI_Offset_kind)
     idx_MOK = int((k-1)*nx1*ny1+(j-1)*nx1+i,MPI_Offset_kind)
     disp = 5*4 + (ii_MOK-1)*WP_MOK
     call MPI_FILE_SET_VIEW(iunit1, disp, MPI_INTEGER, MPI_INTEGER, "native", mpi_info, ierr)
     call MPI_FILE_READ(iunit1, xx(1), 1, MPI_REAL_WP, status, ierr)
     disp = 5*4 + (nd1_MOK+jj_MOK-1)*WP_MOK
     call MPI_FILE_SET_VIEW(iunit1, disp, MPI_INTEGER, MPI_INTEGER, "native", mpi_info, ierr)
     call MPI_FILE_READ(iunit1, xx(2), 1, MPI_REAL_WP, status, ierr)
     disp = 5*4 + (nd1_MOK+nd2_MOK+kk_MOK-1)*WP_MOK
     call MPI_FILE_SET_VIEW(iunit1, disp, MPI_INTEGER, MPI_INTEGER, "native", mpi_info, ierr)
     call MPI_FILE_READ(iunit1, xx(3), 1, MPI_REAL_WP, status, ierr)
     phi = 0.0D0
     do l = 1,n_snap
        snap_MOK = int(l,MPI_Offset_kind)
        do ll = 1,nvars1
           ll_MOK = int(ll+jump_var,MPI_Offset_kind)
           disp = 5*4 + (nd1_MOK+nd2_MOK+nd3_MOK)*WP_MOK + str_MOK*NVARS_MOK &
                  + 2*(snap_MOK)*WP_MOK &
                  + nd1_MOK*nd2_MOK*nd3_MOK*WP_MOK*NVARS_MOK*(snap_MOK-1) &
                  + ((ll_MOK-1)*(nd1_MOK*nd2_MOK*nd3_MOK)+idx_MOK-1)*WP_MOK
           call MPI_FILE_SET_VIEW(iunit1, disp, MPI_INTEGER, MPI_INTEGER, "native", mpi_info, ierr)
           call MPI_FILE_READ(iunit1, vel(ll), 1, MPI_REAL_WP, status, ierr)
           do o = 1,n_snap
              phi(ll,o) = phi(ll,o) + vel(ll)*T(l,o)*LambdaC(o)
           end do
        end do
     end do
     do ll = 1,nvars1
        write(9000+ll,'(1024E24.12E4)') xx(1), xx(2), xx(3),(REAL(phi(ll,l),WP),l=1,n_snap)
     end do
     write(8001,'(1024E24.12E4)') xx(1), xx(2), xx(3), (vel(l),l=1,nvars1)

!        disp = 5*4 + (nd1_MOK+nd2_MOK+nd3_MOK)*WP_MOK + str_MOK*NVARS_MOK &
!               + 2*(snap_MOK)*WP_MOK &
!               + nd1_MOK*nd2_MOK*nd3_MOK*WP_MOK*NVARS_MOK*(snap_MOK-1) &
!               + (idx_MOK-1)*WP_MOK
!        call MPI_FILE_SET_VIEW(iunit1, disp, MPI_INTEGER, MPI_INTEGER, "native", mpi_info, ierr)
!        call MPI_FILE_READ(iunit1, vel(1), 1, MPI_REAL_WP, status, ierr)
!        disp = 5*4 + (nd1_MOK+nd2_MOK+nd3_MOK)*WP_MOK + str_MOK*NVARS_MOK &
!               + 2*(snap_MOK)*WP_MOK &
!               + nd1_MOK*nd2_MOK*nd3_MOK*WP_MOK*NVARS_MOK*(snap_MOK-1) &
!               + ((nd1_MOK*nd2_MOK*nd3_MOK)+idx_MOK-1)*WP_MOK
!        call MPI_FILE_SET_VIEW(iunit1, disp, MPI_INTEGER, MPI_INTEGER, "native", mpi_info, ierr)
!        call MPI_FILE_READ(iunit1, vel(2), 1, MPI_REAL_WP, status, ierr)
!        disp = 5*4 + (nd1_MOK+nd2_MOK+nd3_MOK)*WP_MOK + str_MOK*NVARS_MOK &
!               + 2*(snap_MOK)*WP_MOK &
!               + nd1_MOK*nd2_MOK*nd3_MOK*WP_MOK*NVARS_MOK*(snap_MOK-1) &
!               + (2*(nd1_MOK*nd2_MOK*nd3_MOK)+idx_MOK-1)*WP_MOK
!        call MPI_FILE_SET_VIEW(iunit1, disp, MPI_INTEGER, MPI_INTEGER, "native", mpi_info, ierr)
!        call MPI_FILE_READ(iunit1, vel(3), 1, MPI_REAL_WP, status, ierr)
!        do o = 1,n_snap
!           phi(1,o) = phi(1,o) + vel(1)*T(l,o)*LambdaC(o)
!           phi(2,o) = phi(2,o) + vel(2)*T(l,o)*LambdaC(o)
!           phi(3,o) = phi(3,o) + vel(3)*T(l,o)*LambdaC(o)
!        end do
!     end do
!     write(9001,'(1024E24.12E4)') xx(1), xx(2), xx(3), (phi(1,l),l=1,n_snap)
!     write(9002,'(1024E24.12E4)') xx(1), xx(2), xx(3), (phi(2,l),l=1,n_snap)
!     write(9003,'(1024E24.12E4)') xx(1), xx(2), xx(3), (phi(3,l),l=1,n_snap)
!     write(8001,'(1024E24.12E4)') xx(1), xx(2), xx(3), (vel(l),l=1,3) 
  end do  
  end do
  end do

  call MPI_FILE_CLOSE(iunit1, ierr)

  print *, 'Completed dmd_3D'

end program dmd_3D_p2
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

