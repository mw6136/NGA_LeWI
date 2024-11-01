! ========================================================== !
! Read number of records from binary files using mpi         !
! ========================================================== !
subroutine volumeStats_mpi_read_nrec
  use volumeStats
  implicit none

  integer :: iunit, ierr, i, j, n
  integer :: ifile
  integer :: nvars_
  real(WP), dimension(:), pointer :: x_,y_,z_
  character(len=str_long) :: filename

  ! Print the output header
  if (irank.eq.iroot) print '(2a10,a15)', 'ifile','nseq','time'
  
  ! Loop over the input files
  do ifile=1,nfiles
     ! Open the binary file
     call BINARY_FILE_OPEN(iunit, trim(data_files(ifile)), "r", ierr)

     ! Read data sizes
     call BINARY_FILE_READ(iunit, ntime(ifile), 1, kind(ntime), ierr)
     call BINARY_FILE_READ(iunit, nx_,1, kind(nx_), ierr)
     call BINARY_FILE_READ(iunit, ny_,1, kind(ny_), ierr)
     call BINARY_FILE_READ(iunit, nz_,1, kind(nz_), ierr)
     call BINARY_FILE_READ(iunit, nvars_,1, kind(nvars_), ierr)

     allocate(x_(nx_))
     allocate(y_(ny_))
     allocate(z_(nz_))
     do j=1,nx_
        call BINARY_FILE_READ(iunit, x_(j), 1, kind(x_), ierr)
     end do
     do j=1,ny_
        call BINARY_FILE_READ(iunit, y_(j), 1, kind(y_), ierr)
     end do
     do n=1,nz_
        call BINARY_FILE_READ(iunit, z_(n), 1, kind(z_), ierr)
     end do

     ! Account for the new data and check consistency
     if (ifile.eq.1) then
        nx = nx_
        ny = ny_
        nz = nz_
        nzi = 1.0_WP/real(nz,WP)
        ! Grid in the domain
        allocate(x(1-noverx:nx+noverx))
        allocate(y(1-novery:ny+novery))
        allocate(z(1-nover :nz+nover ))
        x(1:nx) = x_
        y(1:ny) = y_
        z(1:nz) = z_
        ! Grid in the ghost cells
        do i=1-noverx,0
           x(i) = 0.0_WP + real(i-1,WP)*(x(2)-x(1))
        end do
        do i=nx+1,nx+noverx
           x(i) = x(nx) + real(i-nx,WP)*(x(nx)-x(nx-1))
        end do
        do i=1-novery,0
           y(i) = 0.0_WP + real(i-1,WP)*(y(2)-y(1))
        end do
        do i=ny+1,ny+novery
           y(i) = y(ny) + real(i-ny,WP)*(y(ny)-y(ny-1))
        end do
        do i=1-nover,0
           z(i) = 0.0_WP + real(i-1,WP)*(z(2)-z(1))
        end do
        do i=nz+1,nz+nover
           z(i) = z(nz) + real(i-nz,WP)*(z(nz)-z(nz-1))
        end do
        nvars = nvars_
     else
        if (nx.ne.nx_) then
           print *, 'Inconsistent nx in file ', trim(data_files(ifile))
           print *, nx, nx_
           return
        end if
        if (ny.ne.ny_) then
           print *, 'Inconsistent ny in file ', trim(data_files(ifile))
           print *, ny, ny_
           return
        end if
        if (nz.ne.nz_) then
           print *, 'Inconsistent nz in file ', trim(data_files(ifile))
           print *, nz, nz_
           return
        end if
        do i=1,nx
           if (x_(i).ne.x(i)) then
              print *, 'Inconsistency in x in file ', trim(data_files(ifile))
              return
           end if
        end do
        do i=1,ny
           if (y_(i).ne.y(i)) then
              print *, 'Inconsistency in y in file ', trim(data_files(ifile))
              return
           end if
        end do
        do i=1,nz
           if (z_(i).ne.z(i)) then
              print *, 'Inconsistency in z in file ', trim(data_files(ifile))
              return
           end if
        end do
     end if

     ! Read variable names
     if (ifile.eq.1) allocate(names(nvars))
     do n=1,nvars
        call BINARY_FILE_READ(iunit, names(n), str_short, kind(names), ierr)
     end do

     ! Read the timestamp
     call BINARY_FILE_READ(iunit, dt, 1, kind(dt), ierr)
     call BINARY_FILE_READ(iunit, time(ifile), 1, kind(time), ierr)

     ! Print some info
     if (irank.eq.iroot) print '(2i10,ES15.5)', ifile,ntime(ifile),time(ifile)
     !if (irank.eq.iroot) print *, names

     ! Close the file
     call BINARY_FILE_CLOSE(iunit, ierr)
     deallocate(x_)
     deallocate(y_)
     deallocate(z_)
  end do

  return
end subroutine volumeStats_mpi_read_nrec


! ========================================================== !
! Initialize the MPI decomposition                           !
! ========================================================== !
subroutine volumeStats_mpi_init
  use volumeStats
  implicit none

  integer :: var, ierr, isc
  integer, dimension(3) :: gsizes, lsizes, start
  integer :: imins_, imaxs_, kmins_, kmaxs_
  
  ! Allocate the MPI type
  STAT_DATA_NVARS = nvars
  allocate(STAT_DATA(STAT_DATA_NVARS))
  
  ! Set the variable names from the data file
  do var=1,STAT_DATA_NVARS
     STAT_DATA(var)%name = names(var)
  end do

  ! Allocate all arrays to read
  allocate(U  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:3))
  allocate(Ui (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:3))
  allocate(RHO(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(P  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  U   = 0.0_WP
  Ui  = 0.0_WP
  RHO = 0.0_WP
  P   = 0.0_WP
  if (combust) then
     allocate(SC (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nscalar))
     allocate(SC_name(1:nscalar))
     SC  = 0.0_WP
  end if

  ! Link the pointers
  if (use_planes) then
     imins_ = imino_
     imaxs_ = imaxo_
     kmins_ = kmin_
     kmaxs_ = kmax_
  else
     imins_ = imin_
     imaxs_ = imax_
     kmins_ = kmin_
     kmaxs_ = kmax_
  end if
     
  isc = 1
  do var=1,STAT_DATA_NVARS
     select case(trim(STAT_DATA(var)%name))
     case ('U')
        STAT_DATA(var)%var => U(imins_:imaxs_,jmin_:jmax_,kmins_:kmaxs_,1)
     case ('V')
        STAT_DATA(var)%var => U(imins_:imaxs_,jmin_:jmax_,kmins_:kmaxs_,2)
     case ('W')
        STAT_DATA(var)%var => U(imins_:imaxs_,jmin_:jmax_,kmins_:kmaxs_,3)
     case ('P')
        STAT_DATA(var)%var => P(imins_:imaxs_,jmin_:jmax_,kmins_:kmaxs_)
     case ('RHO')
        STAT_DATA(var)%var => RHO(imins_:imaxs_,jmin_:jmax_,kmins_:kmaxs_)
     case('HR')
        ! Don't read HR
        nvars = nvars-1
     case default
        ! Read scalars here
        STAT_DATA(var)%var => SC(imins_:imaxs_,jmin_:jmax_,kmins_:kmaxs_,isc)
        SC_name(isc) = trim(STAT_DATA(var)%name)
        if (trim(STAT_DATA(var)%name).eq.'ZMIX') isc_ZMIX = isc
        if (trim(STAT_DATA(var)%name).eq.'H')    isc_H    = isc
        if (trim(STAT_DATA(var)%name).eq.'O2') then
           isc_PROG = isc
           isc_O2   = isc
        end if
        if (trim(STAT_DATA(var)%name).eq.'O')    isc_O    = isc
        if (trim(STAT_DATA(var)%name).eq.'OH')   isc_OH   = isc
        if (trim(STAT_DATA(var)%name).eq.'H2')   isc_H2   = isc
        if (trim(STAT_DATA(var)%name).eq.'H2O')  then
           !isc_PROG = isc
           isc_H2O  = isc
        end if
        if (trim(STAT_DATA(var)%name).eq.'HO2')  isc_HO2  = isc
        if (trim(STAT_DATA(var)%name).eq.'H2O2') isc_H2O2 = isc
        if (trim(STAT_DATA(var)%name).eq.'T')    isc_T    = isc
        isc = isc + 1
     end select
  end do

  ! Global sizes
  gsizes(1) = nx
  gsizes(2) = ny
  gsizes(3) = nz

  if (use_planes) then
     ! Plane decomposition: read overlap cells
     ! Local sizes
     lsizes(1) = nxo_
     lsizes(2) = ny_
     lsizes(3) = nz_

     ! Starting points
     start(1) = imino_-imin
     start(2) = jmin_-jmin
     start(3) = kmin_-kmin
  else
     ! Volume decomposition: ignore overlap cells
     ! Local sizes
     lsizes(1) = nx_
     lsizes(2) = ny_
     lsizes(3) = nz_

     ! Starting points
     start(1) = imin_-imin
     start(2) = jmin_-jmin
     start(3) = kmin_-kmin
  end if
  
  ! Define the view for each variable
  do var=1,STAT_DATA_NVARS
     call MPI_TYPE_CREATE_SUBARRAY(3, gsizes, lsizes, start, &
          MPI_ORDER_FORTRAN, MPI_REAL_WP, STAT_DATA(var)%view, ierr)
     call MPI_TYPE_COMMIT(STAT_DATA(var)%view, ierr)
  end do

  ! Allocate computed variables
  allocate(VISC  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(DIFF  (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:N_tot+1))
  allocate(src_SC(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:N_tot+1))
  VISC = 0.0_WP
  DIFF = 0.0_WP
  src_SC = 0.0_WP

  if (.not.combust) then
     RHO  = rho_in
     VISC = visc_in
     DIFF = diff_in
     VISC = visc_in
     nscalar = 1
  end if

  ! Allocate gradients
  allocate(dRHOdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3))
  allocate(dPdx  (imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3))
  allocate(dUdx  (imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3,1:3))
  allocate(dSCdx (imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3,1:nscalar))
  allocate(dTAUdx(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3,1:3))

  return
end subroutine volumeStats_mpi_init


! ========================================================== !
! Read input data from binary files using mpi                !
! ========================================================== !
subroutine volumeStats_mpi_read(option)
  use volumeStats
  use volumeStats_metric
  implicit none

  integer, intent(in) :: option
  integer :: iunit, ierr, var, data_size
  integer :: ifile, ifile_prev, irec, j, k, i, isc, ii, jj
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK, str_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK
  character(len=str_long) :: filename
  integer :: isgood, isgood_r

  ! Resize some integers so MPI can read even the biggest files
  nx_MOK  = int(nx,       MPI_Offset_kind)
  ny_MOK  = int(ny,       MPI_Offset_kind)
  nz_MOK  = int(nz,       MPI_Offset_kind)
  WP_MOK  = int(WP,       MPI_Offset_kind)
  str_MOK = int(str_short,MPI_Offset_kind)
  NVARS_MOK = int(STAT_DATA_NVARS,MPI_Offset_kind)

  if (irank.eq.iroot) print '(3a10,2a15)', 'step', 'irec', 'ifile', 'time', 'dt'
  ntime_curr = 0
  ifile_prev = 1

  ! Number of I/O devices to be used for file striping
  if (use_lfs_stripes) then
     call MPI_Info_set(mpi_info, "striping_factor", "16", ierr )
     call MPI_Info_set(mpi_info, "cb_nodes", "16", ierr )
     call MPI_INFO_SET(mpi_info,"romio_cb_read","enable",ierr)
     call MPI_INFO_SET(mpi_info,"romio_ds_read","automatic",ierr)
  end if

  ! Loop over the files
  do irec=1,nrec
     ifile = irec*nfile_step
     ! Get the file name to read in
     filename = trim(data_files(ifile))
     if (use_mpiiofs) filename = trim(mpiiofs) // ":" // trim(filename)
     isgood = 1

     ! Open the file
     call MPI_FILE_OPEN(comm, filename, MPI_MODE_RDONLY, mpi_info, iunit, ierr)

     ! Size of local arrays
     if (use_planes) then
        data_size = nxo_*ny_*nz_
     else
        data_size = nx_*ny_*nz_
     end if
     ntime_curr = ntime_curr + 1
     dt = time(ifile) - time(ifile_prev)

     ! Print time info
     if (irank.eq.iroot) print '(3i10,2ES15.5)', &
          option, irec, ifile, time(ifile), dt

     ! Read data in
     do var=1,nvars
        var_MOK = int(var, MPI_Offset_kind)
        disp = 5*4 + (nx_MOK+ny_MOK+nz_MOK)*WP_MOK + NVARS_MOK*str_MOK + 2*WP_MOK &
             + nx_MOK*ny_MOK*nz_MOK*(var_MOK-1)*WP_MOK
        call MPI_FILE_SET_VIEW(iunit, disp, MPI_REAL_WP, STAT_DATA(var)%view, &
             "native", mpi_info, ierr)
        call MPI_FILE_READ_ALL(iunit, STAT_DATA(var)%var, data_size, &
             MPI_REAL_WP, status, ierr)
     end do

     ! Check data
     do var=1,nvars
        if (trim(STAT_DATA(var)%name).eq.'RHO') then
           if (minval(STAT_DATA(var)%var).lt.1e-10) then
              print *, irank, 'min RHO=', minval(STAT_DATA(var)%var)
              ! It's not good
              isgood = 0
           end if
        end if
     end do

     ! See if there's a failure
     call MPI_REDUCE(isgood, isgood_r, 1, MPI_INTEGER, MPI_MIN, iroot-1, comm, ierr)
     ! Tell our friends if the data is good or not
     if (irank.eq.iroot) isgood = isgood_r
     call MPI_BCAST(isgood, 1, MPI_INTEGER, iroot-1, comm, ierr)

     if (isgood.eq.1) then
        ! Zero the conditional counters for this timestep
        nsyz_cond = 0
        nsz_cond  = 0
        nsz_cond_temp = 0
        
        ! Update borders
        select case(option)
        case(1)
           ! RANS-type budgets
           do ii=1,3
              call volumeStats_mpi_update_border(U(:,:,:,ii))
           end do
           call volumeStats_mpi_update_border(P)
           call volumeStats_mpi_update_border(RHO)
           if (combust) then
              do isc=1,nscalar
                 call volumeStats_mpi_update_border(SC(:,:,:,isc))
              end do
           end if
           
        case(2)
           ! LES analysis
           do ii=1,3
              call volumeStats_mpi_update_volume_border(U(:,:,:,ii))
           end do
           call volumeStats_mpi_update_volume_border(P)
           call volumeStats_mpi_update_volume_border(RHO)
           if (combust) then
              do isc=1,nscalar
                 call volumeStats_mpi_update_volume_border(SC(:,:,:,isc))
              end do
           end if

        case default
           stop 'volumeStats_mpi_read: option not implemented'
        end select

        ! Compute interpolated velocities at cell centers
        if (use_planes) then
           !$OMP PARALLEL DO
           do k=kmino_,kmaxo_
              do j=jmino_,jmaxo_
                 do i=imino_,imaxo_
                    Ui(i,j,k,1) = 0.5_WP*sum(U(i-st1:i+st2,j,k,1))
                    Ui(i,j,k,2) = 0.5_WP*sum(U(i,j-st1:j+st2,k,2))
                    Ui(i,j,k,3) = 0.5_WP*sum(U(i,j,k-st1:k+st2,3))
                 end do
              end do
           end do
           !$OMP END PARALLEL DO
        else
           !$OMP PARALLEL DO
           do k=kmin_-st2,kmax_+st2
              do j=jmin_-st2,jmax_+st2
                 do i=imin_-st2,imax_+st2
                    Ui(i,j,k,1) = sum(interp_u_xm(i,j,:) * U(i-st1:i+st2,j,k,1))
                    Ui(i,j,k,2) = sum(interp_v_ym(i,j,:) * U(i,j-st1:j+st2,k,2))
                    Ui(i,j,k,3) = sum(interp_w_zm(i,j,:) * U(i,j,k-st1:k+st2,3))
                 end do
              end do
           end do
           !$OMP END PARALLEL DO
        end if
        
        ! Compute derived quantities - VISC, DIFF, src_SC, diff_flux
        if (combust) then
           call volumeStats_finitechem_viscosity
           call volumeStats_finitechem_diffusivity
           ! call volumeStats_finitechem_diffusion !! NEED TO IMPLEMENT
           call volumeStats_finitechem_source
        end if

        ! Compute all gradients
        !  -- Must be called after interpolating velocities to Ui, Vi, Wi
        !  -- Must be called after computing viscosity
        !  -- Must be called before saving interpolated velocities to U
        if (need_gradients) call volumeStats_metric_gradients
        
        ! Determine what to do with the data
        select case(option)
        case(1)
           ! Update RANS-averaged quantities
           call volumeStats_average_all
        case(2)
           ! Update LES-filtered quantities
           call volumeStats_les_filter_all
        end select

        ! Update the conditional averaging counters
        nsyz_cond_total = nsyz_cond_total + nsyz_cond
        nsz_cond_total  = nsz_cond_total  + nsz_cond
        nsz_cond_temp_total  = nsz_cond_temp_total  + nsz_cond_temp
     end if

     ! Close the file
     call MPI_FILE_CLOSE(iunit, ierr)
     ifile_prev = ifile
  end do

  return
end subroutine volumeStats_mpi_read



! ========================================================== !
! Update borders for the data fields                         !
! ========================================================== !
subroutine volumeStats_mpi_update_border(A)
  use volumeStats
  implicit none
  integer :: i,j,k
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(inout) :: A

  if (nover.eq.0) return

  ! Update -z overlap cells
  do k=kmino_,kmin_-1
     !$OMP PARALLEL DO
     do j=jmin_,jmax_
        do i=imin_,imax_
           A(i,j,k) = A(i,j,k+kmax_)
        end do
     end do
     !$OMP END PARALLEL DO
  end do
  ! Update +z overlap cells
  do k=kmax_+1,kmaxo_
     !$OMP PARALLEL DO
     do j=jmin_,jmax_
        do i=imin_,imax_
           A(i,j,k) = A(i,j,k-kmax_)
        end do
     end do
     !$OMP END PARALLEL DO
  end do

  ! Update -y overlap cells
  ! NEED TO COMMUNICATE FOR VOLUME AVERAGING/FILTERING
  do k=kmin_,kmax_
     do j=jmino_,jmin_-1
        do i=imin_,imax_
           A(i,j,k) = A(i,jmin,k)
        end do
     end do
  end do
  ! Update +y overlap cells
  do k=kmin_,kmax_
     do j=jmax_+1,jmaxo_
        do i=imin_,imax_
           A(i,j,k) = A(i,jmax,k)
        end do
     end do
  end do

  ! Update -x overlap cells
  ! NEED TO COMMUNICATE FOR VOLUME AVERAGING/FILTERING
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imino_,imin_-1
           A(i,j,k) = A(imin_,j,k)
        end do
     end do
  end do
  ! Update +x overlap cells
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imax_+1,imaxo_
           A(i,j,k) = A(imax_,j,k)
        end do
     end do
  end do

  return
end subroutine volumeStats_mpi_update_border


! ========================================================== !
! Update volume data borders for the data fields             !
! ========================================================== !
subroutine volumeStats_mpi_update_volume_border(A)
  use volumeStats
  implicit none

  integer :: i,j,k
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(inout) :: A
  real(WP), dimension(:,:,:), pointer :: buf1
  real(WP), dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr

  if (nover.eq.0) return
  if (noverx.eq.0) return

  ! Update -z overlap cells
  do k=kmino_,kmin_-1
     !$OMP PARALLEL DO
     do j=jmin_,jmax_
        do i=imin_,imax_
           A(i,j,k) = A(i,j,k+kmax_)
        end do
     end do
     !$OMP END PARALLEL DO
  end do
  ! Update +z overlap cells
  do k=kmax_+1,kmaxo_
     !$OMP PARALLEL DO
     do j=jmin_,jmax_
        do i=imin_,imax_
           A(i,j,k) = A(i,j,k-kmax_)
        end do
     end do
     !$OMP END PARALLEL DO
  end do

  if (.not.use_planes) then
     call volumeStats_mpi_comm_border_x(A,nxo_,nyo_,nzo_,noverx)
     call volumeStats_mpi_comm_border_y(A,nxo_,nyo_,nzo_,novery)
  end if

  return
end subroutine volumeStats_mpi_update_volume_border



!=============================================================================!

subroutine volumeStats_mpi_comm_border_x(A,n1,n2,n3,no)
  use volumeStats
  implicit none
  ! Input parameters
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  real(WP), dimension(n1,n2,n3), intent(inout) :: A
  ! Local variables
  real(WP), dimension(:,:,:), pointer :: buf1
  real(WP), dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr
  
  ! Initialize buffer
  allocate(buf1(no,n2,n3))
  allocate(buf2(no,n2,n3))
  icount = no*n2*n3
  
  ! Copy left buffer
  buf1(:,:,:) = A(1+no:1+no+no-1,:,:)
  
  ! Send left buffer to left neighbour
  call MPI_CART_SHIFT(comm,0,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste left buffer to right
  if (isource.NE.MPI_PROC_NULL) A(n1-no+1:n1,:,:) = buf2(:,:,:)
  
  ! Copy right buffer
  buf1(:,:,:) = A(n1-no-no+1:n1-no,:,:)
  
  ! Send right buffer to right neighbour
  call MPI_CART_SHIFT(comm,0,1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste right buffer to left
  if (isource.NE.MPI_PROC_NULL) A(1:no,:,:) = buf2(:,:,:)
  
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine volumeStats_mpi_comm_border_x

!=============================================================================!

subroutine volumeStats_mpi_comm_border_y(A,n1,n2,n3,no)
  use volumeStats
  implicit none
  ! Input parameters
  integer, intent(in) :: n1
  integer, intent(in) :: n2
  integer, intent(in) :: n3
  integer, intent(in) :: no
  real(WP), dimension(n1,n2,n3), intent(inout) :: A
  ! Local variables
  real(WP), dimension(:,:,:), pointer :: buf1
  real(WP), dimension(:,:,:), pointer :: buf2
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: icount
  integer :: isource
  integer :: idest
  integer :: ierr

  ! Initialize buffer
  allocate(buf1(n1,no,n3))
  allocate(buf2(n1,no,n3))
  icount = n1*no*n3
  
  ! Copy left buffer
  buf1(:,:,:) = A(:,1+no:1+no+no-1,:)
  
  ! Send lower buffer to lower neighbour
  call MPI_CART_SHIFT(comm,1,-1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste lower buffer to top
  if (isource.NE.MPI_PROC_NULL) A(:,n2-no+1:n2,:) = buf2(:,:,:)
  
  ! Copy right buffer
  buf1(:,:,:) = A(:,n2-no-no+1:n2-no,:)
  
  ! Send right buffer to upper neighbour
  call MPI_CART_SHIFT(comm,1,1,isource,idest,ierr)
  call MPI_SENDRECV(buf1,icount,MPI_REAL_WP,idest,0,buf2,icount,MPI_REAL_WP,isource,0,comm,istatus,ierr)
  
  ! Paste right buffer to left
  if (isource.NE.MPI_PROC_NULL) A(:,1:no,:) = buf2(:,:,:)
  
  ! Deallocate
  deallocate(buf1)
  deallocate(buf2)
  nullify(buf1)
  nullify(buf2)
  
  return
end subroutine volumeStats_mpi_comm_border_y

!=============================================================================!



subroutine volumeStats_mpi_topology_init
  use volumeStats
  use parser
  implicit none
  integer :: ierr
  integer, dimension(3) :: dims
  logical, dimension(3) :: isper
  logical, dimension(3) :: dir
  logical :: reorder
  integer :: ndims
  integer, dimension(3) :: coords
  
  ! Read topology from input file
  !  Default is to decompose in x only
  !  Don't decompose in z
  call parser_read('Processors along X',npx,nproc)
  call parser_read('Processors along Y',npy,1)
  !call parser_read('Processors along Z',npz,1)
  
  ! Set MPI topology
  ndims = 3
  dims(1) = npx
  dims(2) = npy
  dims(3) = 1

  isper(1) = .false.
  isper(2) = .false.
  isper(3) = .false.
  reorder = .true.

  call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,isper,reorder,comm,ierr)
  call MPI_COMM_RANK(comm,irank,ierr)
  call MPI_CART_COORDS(comm,irank,ndims,coords,ierr)
  irank = irank + 1
  iproc = coords(1) + 1
  jproc = coords(2) + 1
  kproc = coords(3) + 1
  
  ! Define a root processor at coordinates (1,1,1)
  dims = 0
  call MPI_CART_RANK(comm,dims,iroot,ierr)
  iroot = iroot + 1

  return
end subroutine volumeStats_mpi_topology_init
