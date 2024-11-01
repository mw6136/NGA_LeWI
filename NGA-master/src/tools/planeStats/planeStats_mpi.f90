! ========================================================== !
! Read number of records from binary files using mpi         !
! ========================================================== !
subroutine planeStats_mpi_read_nrec
  use planeStats
  implicit none

  integer :: iunit, ierr, i, j, n, it
  integer :: ifile
  integer :: nvars_
  real(WP), dimension(:), pointer :: d1_, d2_, dp_
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: ndp_MOK, nd1_MOK, nd2_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK, str_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK, NTIME_MOK
  character(len=str_long) :: filename
  real(WP) :: curr_time, stop_time_prev

  stop_time_prev = 0.0_WP
  allocate(nvars(nfiles))

  ! Loop over the input files
  do ifile=1,nfiles
     if (irank.eq.iroot) print *, trim(data_files(ifile))

     ! Open the binary file
     call BINARY_FILE_OPEN(iunit, trim(data_files(ifile)), "r", ierr)

     ! Read data sizes
     call BINARY_FILE_READ(iunit, ntime(ifile), 1, kind(ntime), ierr)
     call BINARY_FILE_READ(iunit, nd1_,    1, kind(nd1_), ierr)
     call BINARY_FILE_READ(iunit, nd2_,    1, kind(nd2_), ierr)
     call BINARY_FILE_READ(iunit, nplanes_,1, kind(nplanes_), ierr)
     call BINARY_FILE_READ(iunit, dplanes, 1, kind(dplanes), ierr)
     call BINARY_FILE_READ(iunit, nvars_,  1, kind(nvars_), ierr)

     allocate(d1_(nd1_))
     allocate(d2_(nd2_))
     allocate(dp_(nplanes_))
     do j=1,nd1_
        call BINARY_FILE_READ(iunit, d1_(j), 1, kind(d1_), ierr)
     end do
     do j=1,nd2_
        call BINARY_FILE_READ(iunit, d2_(j), 1, kind(d2_), ierr)
     end do
     do n=1,nplanes_
        call BINARY_FILE_READ(iunit, dp_(n), 1, kind(dp_), ierr)
     end do

     if (dplanes.eq.1) then
        ! xy planes
        nx = nd1_
        ny = nd2_
        if (irank.eq.iroot) print *, 'planes at z=', dp_
     elseif (dplanes.eq.2) then
        ! xz planes
        nx = nd1_
        nz = nd2_
        if (irank.eq.iroot) print *, 'planes at y=', dp_
     elseif (dplanes.eq.3) then
        ! yz planes
        ny = nd1_
        nz = nd2_
        if (irank.eq.iroot) print *, 'planes at x=', dp_
     end if

     allocate(names(nvars_))
     nvars(ifile) = nvars_
     if (irank.eq.iroot) print *, 'ifile=',ifile,', nvars=',nvars_

     ! Account for the new data and check consistency
     if (ifile.eq.1) then
        allocate(d1(nd1_))
        allocate(d2(nd2_))
        allocate(dp(nplanes_))
        nd1      = nd1_
        nd2      = nd2_
        nplanes  = nplanes_
        d1       = d1_
        d2       = d2_
        dp       = dp_
        if (dplanes.eq.1) then
           ! For AGN
           allocate(x(nx))
           allocate(y(ny))
           x = d1
           y = d2
        elseif (dplanes.eq.3) then
           ! For JFM budgets
           allocate(x(nplanes))
           allocate(y(ny))
           x = dp
           y = d1
        end if
     else
        if (nd1.ne.nd1_) then
           print *, 'Inconsistent nd1 in file ', trim(data_files(ifile))
           print *, nd1, nd1_
           return
        end if
        if (nd2.ne.nd2_) then
           print *, 'Inconsistent nd2 in file ', trim(data_files(ifile))
           print *, nd2, nd2_
           return
        end if
        if (nplanes.ne.nplanes_) then
           print *, 'Inconsistent nplanes in file ', trim(data_files(ifile))
           print *, nplanes, nplanes_
           return
        end if
        do i=1,nd1
           if (d1_(i).ne.d1(i)) then
              print *, 'Inconsistency in d1 in file ', trim(data_files(ifile))
              return
           end if
        end do
        do i=1,nd2
           if (d2_(i).ne.d2(i)) then
              print *, 'Inconsistency in d2 in file ', trim(data_files(ifile))
              return
           end if
        end do
        do i=1,nplanes
           if (dp_(i).ne.dp(i)) then
              print *, 'Inconsistency in dp in file ', trim(data_files(ifile))
              return
           end if
        end do
     end if

     ! Read variable names
     do n=1,nvars_
        call BINARY_FILE_READ(iunit, names(n), str_short, kind(names), ierr)
     end do
     if (irank.eq.iroot) print *, 'ifile=',ifile, names

     ! Read the first timestamp
     call BINARY_FILE_READ(iunit, dt, 1, kind(dt), ierr)
     call BINARY_FILE_READ(iunit, curr_time, 1, kind(curr_time), ierr)

     ! Close the file
     call BINARY_FILE_CLOSE(iunit, ierr)
     deallocate(d1_)
     deallocate(d2_)
     deallocate(dp_)
     if (ifile.lt.nfiles) deallocate(names)

     ! Look for the first unique timestamp in order to avoid duplicate records
     filename = trim(data_files(ifile))
     if (use_mpiiofs) filename = trim(mpiiofs) // ":" // trim(filename)
     call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, mpi_info, iunit, ierr)
     ! Resize some integers so MPI can read even the biggest files
     ndp_MOK   = int(nplanes,       MPI_Offset_kind)
     nd1_MOK   = int(nd1,           MPI_Offset_kind)
     nd2_MOK   = int(nd2,           MPI_Offset_kind)
     WP_MOK    = int(WP,            MPI_Offset_kind)
     str_MOK   = int(str_short,     MPI_Offset_kind)
     NVARS_MOK = int(nvars(ifile),  MPI_Offset_kind)

     it = 1
     do while (it.lt.ntime(ifile) .and. curr_time.le.stop_time_prev)
        NTIME_MOK = int(it, MPI_Offset_kind)
        ! Read the timestamp
        disp = 6*4 + (ndp_MOK+nd1_MOK+nd2_MOK)*WP_MOK + NVARS_MOK*str_MOK + 2*(NTIME_MOK-1)*WP_MOK &
             + ndp_MOK*nd1_MOK*nd2_MOK*NVARS_MOK*(NTIME_MOK-1)*WP_MOK + 1*WP_MOK
        call MPI_FILE_SET_VIEW(iunit, disp, MPI_INTEGER, MPI_INTEGER, &
             "native", mpi_info, ierr)
        call MPI_FILE_READ(iunit, curr_time, 1, MPI_REAL_WP, status, ierr)
        it = it + 1
     end do
     
     ! Set the file's start time
     ntime_start(ifile) = it

     ! Save the file's ending time
     NTIME_MOK = int(ntime(ifile), MPI_Offset_kind)
     disp = 6*4 + (ndp_MOK+nd1_MOK+nd2_MOK)*WP_MOK + NVARS_MOK*str_MOK + 2*(NTIME_MOK-1)*WP_MOK &
          + ndp_MOK*nd1_MOK*nd2_MOK*NVARS_MOK*(NTIME_MOK-1)*WP_MOK + 1*WP_MOK
     call MPI_FILE_SET_VIEW(iunit, disp, MPI_INTEGER, MPI_INTEGER, &
          "native", mpi_info, ierr)
     call MPI_FILE_READ(iunit, stop_time_prev, 1, MPI_REAL_WP, status, ierr)

     ! Close the file
     call MPI_FILE_CLOSE(iunit, ierr)
     if (irank.eq.iroot) print *, ' '
  end do

  return
end subroutine planeStats_mpi_read_nrec


! ========================================================== !
! Initialize the MPI decomposition                           !
! ========================================================== !
subroutine planeStats_mpi_init
  use planeStats
  implicit none

  integer :: var, ierr, isc, isc_dx1, isc_dx2, isc_dx3, isc_src
  integer :: isc_diff, isc_q1, isc_q2, isc_q3, isc_dq
  integer, dimension(3) :: gsizes, lsizes, start
  character(len=str_short) :: sel
  
  ! Allocate the MPI type
  STAT_DATA_NVARS = maxval(nvars)
  allocate(STAT_DATA(STAT_DATA_NVARS))
  
  ! Set the variable names from the data file
  do var=1,STAT_DATA_NVARS
     STAT_DATA(var)%name = names(var)
  end do
  !if (irank.eq.iroot) print *, 'Variables to read: ', STAT_DATA(:)%name

  ! UPDATE THESE COMMENTS
  ! Workspace data structure
  !   (1) Budgets in progress variable space:
  !         pnmin_:pnmax_ = plane indices (major axis)
  !         iprogmin:iprogmax = progress variable bins
  !         (third index unused)
  !   (2) Budgets in mixture fraction space:
  !         pnmin_:pnmax_ = plane indices (major axis)
  !         izmixmin:izmixmax = mixture fraction bins
  !         (third index unused)
  !   (3) Budgets in spectral space:
  !         ipmin_:ipmax_ = probe indices (major axis, y-location)
  !         1:nrec = series of temporal records
  !         1:nz = points in the periodic direction

  ! Allocate all arrays to read
  allocate(SC      (1:nd2,1:nd1,pnmin_:pnmax_,1:nscalar))
  allocate(RHO     (1:nd2,1:nd1,pnmin_:pnmax_))
  allocate(dRHOdx  (1:nd2,1:nd1,pnmin_:pnmax_,1:3))
  allocate(VISC    (1:nd2,1:nd1,pnmin_:pnmax_))
  allocate(P       (1:nd2,1:nd1,pnmin_:pnmax_))
  allocate(U       (1:nd2,1:nd1,pnmin_:pnmax_,1:3))
  allocate(dUdx    (1:nd2,1:nd1,pnmin_:pnmax_,1:3,1:3))
  allocate(dPdx    (1:nd2,1:nd1,pnmin_:pnmax_,1:3))
  allocate(dTAUdx  (1:nd2,1:nd1,pnmin_:pnmax_,1:3,1:3))
  allocate(dSCdx   (1:nd2,1:nd1,pnmin_:pnmax_,1:nscalar,1:3))
  allocate(src_SC  (1:nd2,1:nd1,pnmin_:pnmax_,1:nscalar-1))
  allocate(DIFF    (1:nd2,1:nd1,pnmin_:pnmax_,1:nscalar))
  allocate(diff_flux(1:nd2,1:nd1,pnmin_:pnmax_,1:nscalar,1:3))
  allocate(diff_flux_div(1:nd2,1:nd1,pnmin_:pnmax_,1:nscalar))

  allocate(SC_name(1:nscalar))

  ! Link the pointers
  ! Copy the data while reading in MPI (faster processing),
  !   but keep select case below to count isc
  isc = 1
  isc_dx1 = 1
  isc_dx2 = 1
  isc_dx3 = 1
  isc_src = 1
  isc_diff = 1
  isc_q1 = 1
  isc_q2 = 1
  isc_q3 = 1
  isc_dq = 1
  do var=1,STAT_DATA_NVARS
     sel = STAT_DATA(var)%name(1:3)
     if (trim(sel).eq.'') then
        if (irank.eq.iroot) print *, 'Found a blank variable; skipping'
        exit
     end if
     if (sel.ne.'dx1'.and.sel.ne.'dx2'.and.sel.ne.'dx3'.and.sel.ne.'src'.and.sel.ne.'DIF'.and.sel.ne.'q1_'.and.sel.ne.'q2_'.and.sel.ne.'q3_'.and.sel.ne.'dq_') then
        sel = STAT_DATA(var)%name
     end if
     select case(trim(sel))
     case ('U')
        !STAT_DATA(var)%var => U(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('V')
        !STAT_DATA(var)%var => V(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('W')
        !STAT_DATA(var)%var => W(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('P')
        !STAT_DATA(var)%var => P(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('RHO')
        !STAT_DATA(var)%var => RHO(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('VISC')
        !STAT_DATA(var)%var => VISC(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dUdx11')
        !STAT_DATA(var)%var => dUdx11(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dUdx12')
        !STAT_DATA(var)%var => dUdx12(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dUdx13')
        !STAT_DATA(var)%var => dUdx13(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dUdx21')
        !STAT_DATA(var)%var => dUdx21(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dUdx22')
        !STAT_DATA(var)%var => dUdx22(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dUdx23')
        !STAT_DATA(var)%var => dUdx23(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dUdx31')
        !STAT_DATA(var)%var => dUdx31(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dUdx32')
        !STAT_DATA(var)%var => dUdx32(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dUdx33')
        !STAT_DATA(var)%var => dUdx33(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dRHOdx1')
        !STAT_DATA(var)%var => dRHOdx1(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dRHOdx2')
        !STAT_DATA(var)%var => dRHOdx2(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dRHOdx3')
        !STAT_DATA(var)%var => dRHOdx3(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dPdx1')
        !STAT_DATA(var)%var => dPdx1(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dPdx2')
        !STAT_DATA(var)%var => dPdx2(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dPdx3')
        !STAT_DATA(var)%var => dPdx3(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dTAUdx11')
        !STAT_DATA(var)%var => dTAUdx11(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dTAUdx12')
        !STAT_DATA(var)%var => dTAUdx12(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dTAUdx13')
        !STAT_DATA(var)%var => dTAUdx13(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dTAUdx21')
        !STAT_DATA(var)%var => dTAUdx21(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dTAUdx22')
        !STAT_DATA(var)%var => dTAUdx22(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dTAUdx23')
        !STAT_DATA(var)%var => dTAUdx23(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dTAUdx31')
        !STAT_DATA(var)%var => dTAUdx31(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dTAUdx32')
        !STAT_DATA(var)%var => dTAUdx32(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dTAUdx33')
        !STAT_DATA(var)%var => dTAUdx33(i1min:i1max,i2min:i2max,i3min:i3max)
     case ('dx1')
        isc_dx1 = isc_dx1 + 1
     case ('dx2')
        isc_dx2 = isc_dx2 + 1
     case ('dx3')
        isc_dx3 = isc_dx3 + 1
     case ('src')
        isc_src = isc_src + 1
     case ('DIF')
        isc_DIFF = isc_DIFF + 1
     case ('q1_')
        isc_q1 = isc_q1 + 1
     case ('q2_')
        isc_q2 = isc_q2 + 1
     case ('q3_')
        isc_q3 = isc_q3 + 1
     case ('dq_')
        isc_dq = isc_dq + 1
     case default
        ! Read scalars here
        SC_name(isc) = trim(STAT_DATA(var)%name)
        !STAT_DATA(var)%var => SC(i1min:i1max,i2min:i2max,i3min:i3max,isc)
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
  if (irank.eq.iroot) print *, 'isc=', isc-1, isc_dx1-1, isc_dx2-1, isc_dx3-1, isc_src-1
  if (irank.eq.iroot) print *, '    ', isc_DIFF-1, isc_q1-1, isc_q2-1, isc_q3-1, isc_dq-1

  ! Set all initial values to zero
  do var=1,STAT_DATA_NVARS
     allocate(STAT_DATA(var)%var(pnmin_:pnmax_,1:nd1,1:nd2))
     STAT_DATA(var)%var = 0.0_WP
  end do

  ! Global sizes
  gsizes(1) = nplanes
  gsizes(2) = nd1
  gsizes(3) = nd2

  ! Local sizes
  lsizes(1) = nplanes_
  lsizes(2) = nd1
  lsizes(3) = nd2

  ! Starting points
  start(1) = pnmin_-pnmin
  start(2) = 0
  start(3) = 0
  
  ! Define the view for each variable
  do var=1,STAT_DATA_NVARS
     call MPI_TYPE_CREATE_SUBARRAY(3, gsizes, lsizes, start, &
          MPI_ORDER_FORTRAN, MPI_REAL_WP, STAT_DATA(var)%view, ierr)
     call MPI_TYPE_COMMIT(STAT_DATA(var)%view, ierr)
  end do

  return
end subroutine planeStats_mpi_init


! ========================================================== !
! Read input data from binary files using mpi                !
! ========================================================== !
subroutine planeStats_mpi_read(option)
  use planeStats
  implicit none

  integer, intent(in) :: option
  integer :: vmin,vmax
  integer :: iunit, ierr, var, data_size
  integer :: ifile, it, j, k, i, ii, jj, isc, isc_dx1, isc_dx2, isc_dx3, isc_src
  integer :: isc_diff, isc_q1, isc_q2, isc_q3, isc_dq
  integer :: ntime_prev
  integer :: iplane
  real(WP) :: tmp1, tmp2, time_prev
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nd2_MOK, ndp_MOK, nd1_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK, str_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK, NTIME_MOK
  integer(kind=MPI_Offset_kind) :: PROBE_I_MOK, PROBE_J_MOK
  character(len=str_long) :: filename
  character(len=str_short) :: sel
  integer :: isgood, isgood_r

  ! Default not to use dSC -- update after reading data
  use_dSC = .false.
  use_qSC = .false.

  ! Resize some integers so MPI can read even the biggest files
  nd1_MOK   = int(nd1,           MPI_Offset_kind)
  nd2_MOK   = int(nd2,           MPI_Offset_kind)
  ndp_MOK   = int(nplanes,       MPI_Offset_kind)
  WP_MOK    = int(WP,            MPI_Offset_kind)
  str_MOK   = int(str_short,     MPI_Offset_kind)

  ! Set the initial time index
  ntime_curr = 0
  ntime_dSC  = 0
  ntime_qSC  = 0
  dt = 0.0_WP

  if (irank.eq.iroot) print '(a10,4a7,2a14)', 'step', 'rec', 'nrec', 'f_rec', 'f_nrec', 'time', 'dt'

  ! Loop over the files
  do ifile=1,nfiles

     ! Select the variables to read
     ! In all cases, need scalars to compute the conditional mapping
     select case(option)
     case (1)
        ! This is kind of cheating, but should work
        ! U, V, W, P, N2, ...
        vmin = 5
        vmax = 5+nscalar+13
     case (2)
        ! Conventional budgets - averages
        !! MIGHT NEED TO UPDATE THIS!!
        vmin = 1
        vmax = nvars(ifile) !30+4*nscalar
     case (3)
        ! Conventional budgets
        vmin = 1
        vmax = nvars(ifile)
     case (4)
        ! Scatter data
        vmin = 1
        vmax = 6 + nscalar !nvars(ifile)
     case (5)
        ! Conditional budgets
        vmin = 1
        vmax = nvars(ifile) !6+2*9+2*3+4*nscalar !nvars(ifile)
     end select

     !if (irank.eq.iroot) print *, STAT_DATA(vmin:vmax)%name
     
     ! Account for more variables in subsequent data files
     NVARS_MOK = int(nvars(ifile),MPI_Offset_kind)
     ! Get the file name to read in
     filename = trim(data_files(ifile))
     if (use_mpiiofs) filename = trim(mpiiofs) // ":" // trim(filename)
     
     ! Open the file
     call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, mpi_info, iunit, ierr)

     ! Set the start and end points
     ntime_prev = ntime_curr

     ! Size of local arrays
     data_size = nplanes_*nd1_*nd2_
        
     ! Set flags if scalar gradients and species diffusive fluxes are available in this file
     use_dSC = .false.
     use_qSC = .false.

     ! Loop over time
     it     = ntime_start(ifile)
     isgood = 1
     dt     = 0.0_WP
     do while (it.le.ntime(ifile).and.isgood.eq.1)
        ntime_curr = ntime_prev + 1 + (it-ntime_start(ifile))/stepsize
        NTIME_MOK = int(it, MPI_Offset_kind)
       
        ! Read the time
        disp = 6*4 + (ndp_MOK+nd1_MOK+nd2_MOK)*WP_MOK + NVARS_MOK*str_MOK + (2*(NTIME_MOK-1)+1)*WP_MOK &
             + ndp_MOK*nd1_MOK*nd2_MOK*NVARS_MOK*(NTIME_MOK-1)*WP_MOK
        call MPI_FILE_SET_VIEW(iunit, disp, MPI_INTEGER, MPI_INTEGER, &
             "native", mpi_info, ierr)
        call MPI_FILE_READ(iunit, time(ntime_curr), 1, MPI_REAL_WP, status, ierr)
        if (ntime_curr.gt.1) dt = time(ntime_curr) - time_prev

        ! Print time info
        if (irank.eq.iroot) print '(i10,4i7,2ES14.5)', &
             option, ntime_curr, nrec, it, ntime(ifile), time(ntime_curr), dt

        ! Read data in
        do var=vmin,vmax
           var_MOK = int(var, MPI_Offset_kind)
           disp = 6*4 + (ndp_MOK+nd1_MOK+nd2_MOK)*WP_MOK + NVARS_MOK*str_MOK + 2*NTIME_MOK*WP_MOK &
                + ndp_MOK*nd1_MOK*nd2_MOK*NVARS_MOK*(NTIME_MOK-1)*WP_MOK &
                + ndp_MOK*nd1_MOK*nd2_MOK*(var_MOK-1)*WP_MOK
           call MPI_FILE_SET_VIEW(iunit, disp, MPI_REAL_WP, STAT_DATA(var)%view, &
                "native", mpi_info, ierr)
           call MPI_FILE_READ(iunit, STAT_DATA(var)%var, data_size, &
                MPI_REAL_WP, status, ierr)
        end do

        ! Check data
        do var=vmin,vmax
!!$           do k=1,nd2
!!$              do j=1,nd1
!!$                 do i=pnmin_,pnmax_
!!$                    if (STAT_DATA(var)%var(i,j,k).ne.STAT_DATA(var)%var(i,j,k)) then
!!$                       print *, irank, 'NaN in ', STAT_DATA(var)%name
!!$                    end if
!!$                 end do
!!$              end do
!!$           end do
!!$           if (maxval(abs(STAT_DATA(var)%var)).gt.1e30) then
!!$              print *, irank, 'Inf in ', STAT_DATA(var)%name
!!$              isgood = 0
!!$           end if
           if (trim(STAT_DATA(var)%name).eq.'RHO') then
              if (minval(STAT_DATA(var)%var).lt.1e-10) then
                 print *, irank, 'min RHO=', minval(STAT_DATA(var)%var)
                 ! It's not good
                 isgood = 0
              end if
           end if
        end do

        ! See if there's a failure
        call MPI_REDUCE(isgood, isgood_r, 1, MPI_INTEGER, MPI_MIN, iroot-1, MPI_COMM_WORLD, ierr)
        ! Tell our friends if the data is good or not
        if (irank.eq.iroot) isgood = isgood_r
        call MPI_BCAST(isgood, 1, MPI_INTEGER, iroot-1, MPI_COMM_WORLD, ierr)

        if (isgood.eq.1) then
           ! Zero the conditional counters for this timestep
           ns_cond = 0

           ! Re-order the arrays
           isc = 1
           isc_dx1 = 1
           isc_dx2 = 1
           isc_dx3 = 1
           isc_src = 1
           isc_diff = 1
           isc_q1 = 1
           isc_q2 = 1
           isc_q3 = 1
           isc_dq = 1
           do var=vmin,vmax
              sel = STAT_DATA(var)%name(1:3)
              if (trim(sel).eq.'') then
                 !if (irank.eq.iroot) print *, 'Found a blank variable; skipping'
                 exit
              end if
              if (sel.eq.'dx1'.or.sel.eq.'dx2'.or.sel.eq.'dx3'.or.sel.eq.'src') then
                 ! Compute scalar gradients
                 use_dSC = .true.
              else if (sel.eq.'DIF'.or.sel.eq.'q1_'.or.sel.eq.'q2_'.or.sel.eq.'q3_'.or.sel.eq.'dq_') then
                 ! Compute diffusive fluxes
                 use_qSC = .true.
              else
                 sel = STAT_DATA(var)%name
              end if
              select case(trim(sel))
              case ('U')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          U(k,j,i,1) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('V')
                 if (dplanes.eq.1) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          ! Flip negative side of jet to positive side
                          do k=1,nd2/2-1
                             U(k,j,i,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=1,nd1
                          do k=nd2/2,nd2
                             U(k,j,i,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 elseif (dplanes.eq.3) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1/2-1
                          ! Flip negative side of jet to positive side
                          do k=1,nd2
                             U(k,j,i,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=nd1/2,nd1
                          do k=1,nd2
                             U(k,j,i,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          do k=1,nd2
                             U(k,j,i,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 end if
              case ('W')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          U(k,j,i,3) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('P')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          P(k,j,i) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('RHO')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          RHO(k,j,i) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('VISC')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          VISC(k,j,i) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dUdx11')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dUdx(k,j,i,1,1) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dUdx12')
                 if (dplanes.eq.1) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          ! Flip negative side of jet to positive side
                          do k=1,nd2/2-1
                             dUdx(k,j,i,1,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=1,nd1
                          do k=nd2/2,nd2
                             dUdx(k,j,i,1,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 elseif (dplanes.eq.3) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1/2-1
                          do k=1,nd2
                             dUdx(k,j,i,1,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=nd1/2,nd1
                          do k=1,nd2
                             dUdx(k,j,i,1,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          do k=1,nd2
                             dUdx(k,j,i,1,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 end if
              case ('dUdx13')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dUdx(k,j,i,1,3) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dUdx21')
                 if (dplanes.eq.1) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          ! Flip negative side of jet to positive side
                          do k=1,nd2/2-1
                             dUdx(k,j,i,2,1) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=1,nd1
                          do k=nd2/2,nd2
                             dUdx(k,j,i,2,1) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 elseif (dplanes.eq.3) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1/2-1
                          do k=1,nd2
                             dUdx(k,j,i,2,1) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=nd1/2,nd1
                          do k=1,nd2
                             dUdx(k,j,i,2,1) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          do k=1,nd2
                             dUdx(k,j,i,2,1) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 end if
              case ('dUdx22')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dUdx(k,j,i,2,2) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dUdx23')
                 if (dplanes.eq.1) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          ! Flip negative side of jet to positive side
                          do k=1,nd2/2-1
                             dUdx(k,j,i,2,3) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=1,nd1
                          do k=nd2/2,nd2
                             dUdx(k,j,i,2,3) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 elseif (dplanes.eq.3) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1/2-1
                          do k=1,nd2
                             dUdx(k,j,i,2,3) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=nd1/2,nd1
                          do k=1,nd2
                             dUdx(k,j,i,2,3) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          do k=1,nd2
                             dUdx(k,j,i,2,3) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 end if
              case ('dUdx31')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dUdx(k,j,i,3,1) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dUdx32')
                 if (dplanes.eq.1) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          ! Flip negative side of jet to positive side
                          do k=1,nd2/2-1
                             dUdx(k,j,i,3,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=1,nd1
                          do k=nd2/2,nd2
                             dUdx(k,j,i,3,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 elseif (dplanes.eq.3) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1/2-1
                          do k=1,nd2
                             dUdx(k,j,i,3,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=nd1/2,nd1
                          do k=1,nd2
                             dUdx(k,j,i,3,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          do k=1,nd2
                             dUdx(k,j,i,3,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 end if
              case ('dUdx33')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dUdx(k,j,i,3,3) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dRHOdx1')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dRHOdx(k,j,i,1) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dRHOdx2')
                 if (dplanes.eq.1) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          ! Flip negative side of jet to positive side
                          do k=1,nd2/2-1
                             dRHOdx(k,j,i,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=1,nd1
                          do k=nd2/2,nd2
                             dRHOdx(k,j,i,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 elseif (dplanes.eq.3) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1/2-1
                          do k=1,nd2
                             dRHOdx(k,j,i,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=nd1/2,nd1
                          do k=1,nd2
                             dRHOdx(k,j,i,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          do k=1,nd2
                             dRHOdx(k,j,i,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 end if
              case ('dRHOdx3')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dRHOdx(k,j,i,3) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dPdx1')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dPdx(k,j,i,1) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dPdx2')
                 if (dplanes.eq.1) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          ! Flip negative side of jet to positive side
                          do k=1,nd2/2-1
                             dPdx(k,j,i,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=1,nd1
                          do k=nd2/2,nd2
                             dPdx(k,j,i,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 elseif (dplanes.eq.3) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1/2-1
                          do k=1,nd2
                             dPdx(k,j,i,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=nd1/2,nd1
                          do k=1,nd2
                             dPdx(k,j,i,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          do k=1,nd2
                             dPdx(k,j,i,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 end if
              case ('dPdx3')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dPdx(k,j,i,3) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dTAUdx11')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dTAUdx(k,j,i,1,1) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dTAUdx12')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dTAUdx(k,j,i,1,2) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dTAUdx13')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dTAUdx(k,j,i,1,3) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dTAUdx21')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dTAUdx(k,j,i,2,1) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dTAUdx22')
                 if (dplanes.eq.1) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          ! Flip negative side of jet to positive side
                          do k=1,nd2/2-1
                             dTAUdx(k,j,i,2,2) = -2.0_WP*STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=1,nd1
                          do k=nd2/2,nd2
                             dTAUdx(k,j,i,2,2) = 2.0_WP*STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 elseif (dplanes.eq.3) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1/2-1
                          do k=1,nd2
                             dTAUdx(k,j,i,2,2) = -2.0_WP*STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=nd1/2,nd1
                          do k=1,nd2
                             dTAUdx(k,j,i,2,2) = 2.0_WP*STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          do k=1,nd2
                             dTAUdx(k,j,i,2,2) = 2.0_WP*STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 end if
              case ('dTAUdx23')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dTAUdx(k,j,i,2,3) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dTAUdx31')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dTAUdx(k,j,i,3,1) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dTAUdx32')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dTAUdx(k,j,i,3,2) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dTAUdx33')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dTAUdx(k,j,i,3,3) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
              case ('dx1')
                 ! Scalar gradients x
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dSCdx(k,j,i,isc_dx1,1) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
                 isc_dx1 = isc_dx1 + 1
              case ('dx2')
                 ! Scalar gradients y
                 if (dplanes.eq.1) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          ! Flip negative side of jet to positive side
                          do k=1,nd2/2-1
                             dSCdx(k,j,i,isc_dx2,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=1,nd1
                          do k=nd2/2,nd2
                             dSCdx(k,j,i,isc_dx2,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 elseif (dplanes.eq.3) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1/2-1
                          do k=1,nd2
                             dSCdx(k,j,i,isc_dx2,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=nd1/2,nd1
                          do k=1,nd2
                             dSCdx(k,j,i,isc_dx2,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          do k=1,nd2
                             dSCdx(k,j,i,isc_dx2,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 end if
                 isc_dx2 = isc_dx2 + 1
              case ('dx3')
                 ! Scalar gradients z
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          dSCdx(k,j,i,isc_dx3,3) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
                 isc_dx3 = isc_dx3 + 1
              case ('src')
                 ! Scalar source terms
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          src_SC(k,j,i,isc_src) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
                 isc_src = isc_src + 1
              case ('DIF')
                 ! Species diffusivities
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          DIFF(k,j,i,isc_DIFF) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
                 isc_DIFF = isc_DIFF + 1
              case ('q1_')
                 ! Species diffusive fluxes (1,2,3)
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          diff_flux(k,j,i,isc_q1,1) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
                 isc_q1 = isc_q1 + 1
              case ('q2_')
                 if (dplanes.eq.1) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          ! Flip negative side of jet to positive side
                          do k=1,nd2/2-1
                             diff_flux(k,j,i,isc_q2,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=1,nd1
                          do k=nd2/2,nd2
                             diff_flux(k,j,i,isc_q2,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 elseif (dplanes.eq.3) then
                    do i=pnmin_,pnmax_
                       do j=1,nd1/2-1
                          do k=1,nd2
                             diff_flux(k,j,i,isc_q2,2) = -STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                       do j=nd1/2,nd1
                          do k=1,nd2
                             diff_flux(k,j,i,isc_q2,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 else
                    do i=pnmin_,pnmax_
                       do j=1,nd1
                          do k=1,nd2
                             diff_flux(k,j,i,isc_q2,2) = STAT_DATA(var)%var(i,j,k)
                          end do
                       end do
                    end do
                 end if
                 isc_q2 = isc_q2 + 1
              case ('q3_')
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          diff_flux(k,j,i,isc_q3,3) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
                 isc_q3 = isc_q3 + 1
              case ('dq_')
                 ! Species diffusive flux divergence
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          diff_flux_div(k,j,i,isc_dq) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
                 isc_dq = isc_dq + 1
              case default
                 ! Read scalars here
                 do i=pnmin_,pnmax_
                    do j=1,nd1
                       do k=1,nd2
                          SC(k,j,i,isc) = STAT_DATA(var)%var(i,j,k)
                       end do
                    end do
                 end do
                 isc = isc + 1
              end select
           end do

           if (use_dSC) ntime_dSC = ntime_dSC + 1
           if (use_qSC) ntime_qSC = ntime_qSC + 1

           ! Determine what to do with the data
           select case(option)
           case(2)
              ! Update averages
              call planeStats_average_all
           case(3)
              ! Compute the budgets
              call planeStats_budgets_compute
           case(4)
              ! Write scatterplot data
              call planeStats_scatter_write
           case(5)
              ! Conditional budgets
              call planeStats_cond_compute
           end select
           

           ! Update the conditional averaging counters
           ns_cond_total = ns_cond_total + ns_cond
        end if

        ! Increment the step
        it = it + stepsize
        time_prev = time(ntime_curr)
     end do

     ! Close the file
     call MPI_FILE_CLOSE(iunit, ierr)
  end do

  ! Write some output
  open(unit=iunit, file='plane_data', action='write')
  do j=1,ny
     do iplane=pnmin_,pnmax_
        if (iplane.eq.pnmin_) then
           write(iunit,'(ES22.13)',advance='no') y(j) ! y-coordinate
        end if
        do jj=1,3
           write(iunit,'(ES22.13)',advance='no') Um(j,iplane,jj)
        end do
        do ii=1,3
           do jj=1,3
              if (iplane.eq.pnmax_ .and. ii.eq.3.and.jj.eq.3) then
                 write(iunit,'(ES22.13)') dTAUdxm(j,iplane,ii,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') dTAUdxm(j,iplane,ii,jj)
              end if
           end do
        end do
!!$        if (iplane.eq.pnmax_) then
!!$           write(iunit,'(3ES22.13)') VISC(nz/2,j,iplane), DIFF(nz/2,j,iplane,isc_H2O), src_SC(nz/2,j,iplane,isc_H2O)
!!$        else
!!$           write(iunit,'(3ES22.13)',advance='no') VISC(nz/2,j,iplane), DIFF(nz/2,j,iplane,isc_H2O), src_SC(nz/2,j,iplane,isc_H2O)
!!$        end if
!!$        write(iunit,'(ES22.13)',advance='no') RHO(nz/2,j,iplane)
!!$        do ii=1,3
!!$           write(iunit,'(ES22.13)',advance='no') U(nz/2,j,iplane,ii)
!!$        end do
!!$        do ii=1,3
!!$           do jj=1,3
!!$              write(iunit,'(ES22.13)',advance='no') dUdx(nz/2,j,iplane,ii,jj)
!!$           end do
!!$        end do
!!$        do jj=1,3
!!$           write(iunit,'(ES22.13)',advance='no') dPdx(nz/2,j,iplane,jj)
!!$        end do
!!$        do jj=1,3
!!$           write(iunit,'(ES22.13)',advance='no') dRHOdx(nz/2,j,iplane,jj)
!!$        end do
!!$        do isc=1,nscalar
!!$           do jj=1,3
!!$              if (iplane.eq.pnmax_ .and. isc.eq.nscalar .and. jj.eq.3) then
!!$                 write(iunit,'(ES22.13)') dSCdx(nz/2,j,iplane,isc,jj)
!!$              else
!!$                 write(iunit,'(ES22.13)',advance='no') dSCdx(nz/2,j,iplane,isc,jj)
!!$              end if
!!$           end do
!!$        end do
     end do
  end do
  close(iunit)
  print *, 'output done'

  return
end subroutine planeStats_mpi_read
