! ========================================================== !
!                       dump_volume.f90                      !
! Dumps some variables in 3D subdomain of specified size     !
!                                                            !
! Author: T. Grenga                                          !
! Date:   August 17, 2016                                    !
! ========================================================== !
module dump_volume
  use parallel
  use precision
  use geometry
  use partition
  use masks
  use fileio
  use data
  use dump_data
  implicit none

  logical :: in_domain

  ! Indices of edges locations
  integer :: vxmin, vxmax, vxmin_, vxmax_
  integer :: vymin, vymax, vymin_, vymax_
  integer :: vzmin, vzmax, vzmin_, vzmax_

  ! MPI IO variables
  integer :: VOLUME_IO_COMM, VOLUME_NON_COMM
  integer :: irank_p, nsteps_out
  character(len=str_medium_2) :: filename_base, folder_base
  integer :: VOLUME_IO_NVARS
  type(MPI_IO_VAR), dimension(:), pointer :: VOLUME_IO_DATA

end module dump_volume


! ========================================================== !
! Initialization                                             !
! ========================================================== !
subroutine dump_volume_init
  use dump_volume
  use parser
  use time_info
  implicit none
  integer :: i,j,k,n,var,ierr,isc
  integer, dimension(3) :: gsizes, lsizes, start
  real(WP) :: vol_xmin, vol_xmax, vol_ymin, vol_ymax, vol_zmin, vol_zmax, tmp
  character(len=str_medium) :: cmd
  real(WP), dimension(:,:,:), pointer :: global_arr
  character(len=str_short) :: buf
  integer :: lfs_stripe

  ! Create and start the timer
  call timing_create('dump_volI')
  call timing_create('dump_volC')
  call timing_start ('dump_volI')

  call parser_read('Volume data folder',folder_base,'volume_data')
  if (irank.eq.iroot) call system("mkdir -p "//trim(folder_base))

  ! Increase the LUSTRE stripe size to speed file writes
  if (trim(mpiiofs).eq.'lustre'.and.irank.eq.iroot) then
     call parser_read('Volume folder lfs stripes',lfs_stripe,16)
     write(buf, '(I3)') lfs_stripe
     print *, "lfs setstripe -c "//buf//" $PWD/"//trim(folder_base)
     call system("lfs setstripe -c "//buf//" $PWD/"//trim(folder_base))
  end if

  ! Read locations of volume's edges
  ! Determine if belong to the current process
  call parser_read('Volume xmin',vol_xmin)
  call parser_read('Volume xmax',vol_xmax)
  call parser_read('Volume ymin',vol_ymin)
  call parser_read('Volume ymax',vol_ymax)
  call parser_read('Volume zmin',vol_zmin)
  call parser_read('Volume zmax',vol_zmax)

!!$  print *, vol_xmin, x(imin)
!!$  print *, vol_xmax, x(imax+1)
!!$  print *, vol_ymin, y(jmin)
!!$  print *, vol_ymax, y(jmax+1)
!!$  print *, vol_zmin, z(kmin)
!!$  print *, vol_zmax, z(kmax+1)

  if (vol_xmin.gt.vol_xmax) &
     call die('dump_volume_init: xmin must be less than xmax')
  if (vol_ymin.gt.vol_ymax) &
     call die('dump_volume_init: ymin must be less than ymax')
  if (vol_zmin.gt.vol_zmax) &
     call die('dump_volume_init: zmin must be less than zmax')
  if (vol_xmin.lt.x(imin) .or. vol_ymin.lt.y(jmin) .or. vol_zmin.lt.z(kmin) .or. &
      vol_xmax.gt.x(imax+1) .or. vol_ymax.gt.y(jmax+1) .or. vol_zmax.gt.z(kmax+1)) &
     call die('dump_volume_init: volume not in domain')

  if ((vol_xmin.le.x(imax_) .and. vol_xmax.gt.x(imin_)) .and. &
      (vol_ymin.le.y(jmax_) .and. vol_ymax.gt.y(jmin_)) .and. &
      (vol_zmin.le.z(kmax_) .and. vol_zmax.gt.z(kmin_))) &
     in_domain = .true.

  if (in_domain) then
     ! Indices in x
     if (vol_xmin.ge.x(imin_)) then 
        j = imin_
        do while (x(j).lt.vol_xmin)
           j = j+1
        end do
        vxmin_ = j
     else 
        vxmin_  = imin_
     end if
     if (vol_xmax.le.x(imax_)) then 
        j = imin_
        do while (x(j).lt.vol_xmax)
           j = j+1
        end do
        vxmax_ = j
     else 
        vxmax_  = imax_
     end if

     ! Indices in y
     if (vol_ymin.ge.y(jmin_)) then 
        j = jmin_
        do while (y(j).lt.vol_ymin)
           j = j+1
        end do
        vymin_ = j
     else 
        vymin_  = jmin_
     end if
     if (vol_ymax.le.y(jmax_)) then
        j = jmin_
        do while (y(j).lt.vol_ymax)
           j = j+1
        end do
        vymax_ = j
     else 
        vymax_  = jmax_
     end if

     ! Indices in z
     if (vol_zmin.ge.z(kmin_)) then
        j = kmin_
        do while (z(j).lt.vol_zmin)
           j = j+1
        end do
        vzmin_ = j
     else
        vzmin_  = kmin_
     end if
     if (vol_zmax.le.z(kmax_)) then
        j = kmin_
        do while (z(j).lt.vol_zmax)
           j = j+1
        end do
        vzmax_ = j
     else
        vzmax_  = kmax_
     end if
  end if

  ! Create a communicator for processes with planes to dump
  irank_p = -1
  if (in_domain) then
     call MPI_COMM_SPLIT(MPI_COMM_WORLD, 1, irank, VOLUME_IO_COMM, ierr)
     call MPI_COMM_RANK(VOLUME_IO_COMM, irank_p, ierr)
  else
     call MPI_COMM_SPLIT(MPI_COMM_WORLD, 0, irank, VOLUME_NON_COMM, ierr)
  end if

  ! Get information on all the planes from the communicator
  if (in_domain) then
     ! Get the total nx, ny, nz
     call MPI_ALLREDUCE(vxmin_, vxmin, 1, MPI_INTEGER, MPI_MIN, VOLUME_IO_COMM, ierr)
     call MPI_ALLREDUCE(vxmax_, vxmax, 1, MPI_INTEGER, MPI_MAX, VOLUME_IO_COMM, ierr)
     call MPI_ALLREDUCE(vymin_, vymin, 1, MPI_INTEGER, MPI_MIN, VOLUME_IO_COMM, ierr)
     call MPI_ALLREDUCE(vymax_, vymax, 1, MPI_INTEGER, MPI_MAX, VOLUME_IO_COMM, ierr)
     call MPI_ALLREDUCE(vzmin_, vzmin, 1, MPI_INTEGER, MPI_MIN, VOLUME_IO_COMM, ierr)
     call MPI_ALLREDUCE(vzmax_, vzmax, 1, MPI_INTEGER, MPI_MAX, VOLUME_IO_COMM, ierr)

     ! Get the name of the output file
     call parser_read('Volume data file to write',filename_base)
     nsteps_out = 0

     ! Set the number of output variables
     if (combust) then
        VOLUME_IO_NVARS = 6 + nscalar
     else
        VOLUME_IO_NVARS = 4
     end if
     allocate(VOLUME_IO_DATA(VOLUME_IO_NVARS))

     ! Name the output variables
     VOLUME_IO_DATA(1)%name = 'U'
     VOLUME_IO_DATA(2)%name = 'V'
     VOLUME_IO_DATA(3)%name = 'W'
     VOLUME_IO_DATA(4)%name = 'P'
     if (combust) then
        VOLUME_IO_DATA(5)%name = 'RHO'
        do isc=1,nscalar
           VOLUME_IO_DATA(5+isc)%name = trim(SC_name(isc))
        end do
        VOLUME_IO_DATA(6+nscalar)%name = 'HR'
     end if

     ! Link the variables
     do var=1,VOLUME_IO_NVARS
        select case(trim(VOLUME_IO_DATA(var)%name))
        case ('U')
           VOLUME_IO_DATA(var)%var => U(vxmin_:vxmax_,vymin_:vymax_,vzmin_:vzmax_)
        case ('V')
           VOLUME_IO_DATA(var)%var => V(vxmin_:vxmax_,vymin_:vymax_,vzmin_:vzmax_)
        case ('W')
           VOLUME_IO_DATA(var)%var => W(vxmin_:vxmax_,vymin_:vymax_,vzmin_:vzmax_)
        case ('P')
           VOLUME_IO_DATA(var)%var => P(vxmin_:vxmax_,vymin_:vymax_,vzmin_:vzmax_)
        case default
           ! Need to interpolate density and scalars in time
           allocate(VOLUME_IO_DATA(var)%var(vxmin_:vxmax_,vymin_:vymax_,vzmin_:vzmax_))
        end select
     end do

     ! Global sizes
     gsizes(1) = vxmax - vxmin + 1
     gsizes(2) = vymax - vymin + 1
     gsizes(3) = vzmax - vzmin + 1

     ! Local sizes
     lsizes(1) = vxmax_ - vxmin_ + 1
     lsizes(2) = vymax_ - vymin_ + 1
     lsizes(3) = vzmax_ - vzmin_ + 1

     ! Starting points
     start(1) = vxmin_ - vxmin
     start(2) = vymin_ - vymin
     start(3) = vzmin_ - vzmin

     ! Define the view for each variable
     do var=1,VOLUME_IO_NVARS
        call MPI_TYPE_CREATE_SUBARRAY(3, gsizes, lsizes, start, &
             MPI_ORDER_FORTRAN, MPI_REAL_WP, VOLUME_IO_DATA(var)%view, ierr)
        call MPI_TYPE_COMMIT(VOLUME_IO_DATA(var)%view, ierr)
     end do

  end if

  ! Stop the timer
  call timing_stop('dump_volI')

  call monitor_log("dump_volume_init COMPLETED")

  return
end subroutine dump_volume_init


! ========================================================== !
! Dump all data on the volume                                !
! ========================================================== !
subroutine dump_volume_3D
  use dump_volume
  implicit none

  integer :: n, j
  integer :: ifile, ierr, var, data_size
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(5) :: dims
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nd1_MOK, nd2_MOK, nd3_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK, str_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK, NTIME_MOK
  integer :: vd1, vd2, vd3, vd1_, vd2_, vd3_
  character(len=str_medium) :: buffer, cmd
  character(len=str_medium_2) :: filename
  logical :: file_is_there
  real(WP) :: tmp

  ! Start the timer
  call timing_start('dump_volC')

  if (in_domain) then
     write(filename, '(''_'', I8.8)') output_save(isave_vol)
     filename = trim(folder_base)//'/'//trim(filename_base)//trim(filename)
     if (use_mpiiofs) filename = trim(mpiiofs)//":" // trim(filename)

     ! Increment the number of output steps
     nsteps_out = nsteps_out + 1

     ! Resize some integers so MPI can write even the biggest files
     WP_MOK    = int(WP,            MPI_Offset_kind)
     str_MOK   = int(str_short,     MPI_Offset_kind)
     NVARS_MOK = int(VOLUME_IO_NVARS,MPI_Offset_kind)
     NTIME_MOK = int(nsteps_out,    MPI_Offset_kind)

     vd1 = vxmax - vxmin + 1
     vd2 = vymax - vymin + 1
     vd3 = vzmax - vzmin +1
     vd1_ = vxmax_ - vxmin_ + 1
     vd2_ = vymax_ - vymin_ + 1
     vd3_ = vzmax_ - vzmin_ + 1
     nd1_MOK = int(vd1, MPI_Offset_kind)
     nd2_MOK = int(vd2, MPI_Offset_kind)
     nd3_MOK = int(vd3, MPI_Offset_kind)

     ! Interpolate data in x and compute derivatives
     if (combust) call dump_volume_compute 

     ! Open the file to write
     call MPI_FILE_OPEN(VOLUME_IO_COMM, trim(filename), MPI_MODE_WRONLY+MPI_MODE_CREATE, mpi_info, ifile, ierr)

     ! Write the header
     if (irank_p.eq.0) then
        ! Write dimensions
        dims(1) = nsteps_out
        dims(2) = vd1
        dims(3) = vd2
        dims(4) = vd3
        dims(5) = VOLUME_IO_NVARS
        call MPI_FILE_WRITE(ifile, dims, 5, MPI_INTEGER, status, ierr)

        ! Write locations directions
        do j=vxmin,vxmax
           call MPI_FILE_WRITE(ifile, x(j), 1, MPI_REAL_WP, status, ierr)
        end do
        do j=vymin,vymax
           call MPI_FILE_WRITE(ifile, y(j), 1, MPI_REAL_WP, status, ierr)
        end do
        do j=vzmin,vzmax
           call MPI_FILE_WRITE(ifile, z(j), 1, MPI_REAL_WP, status, ierr)
        end do
        ! Write variable names
        do var=1,VOLUME_IO_NVARS
           call MPI_FILE_WRITE(ifile,VOLUME_IO_DATA(var)%name,str_short,MPI_CHARACTER,status,ierr)
        end do
        ! Write time info
        call MPI_FILE_WRITE(ifile,dt,  1,MPI_REAL_WP,status,ierr)
        call MPI_FILE_WRITE(ifile,time,1,MPI_REAL_WP,status,ierr)
     end if

     ! Size of local arrays
     data_size = vd1_*vd2_*vd3_

     ! Write the data for each variable
     do var=1,VOLUME_IO_NVARS
        var_MOK = int(var,MPI_Offset_kind)
        disp = 5*4 + (nd1_MOK+nd2_MOK+nd3_MOK)*WP_MOK + str_MOK*NVARS_MOK + 2*WP_MOK &
               + nd1_MOK*nd2_MOK*nd3_MOK*WP_MOK*(var_MOK-1)
        call MPI_FILE_SET_VIEW(ifile, disp, MPI_REAL_WP, VOLUME_IO_DATA(var)%view, &
             "native", mpi_info, ierr)
        call MPI_FILE_WRITE_ALL(ifile, VOLUME_IO_DATA(var)%var, data_size, &
             MPI_REAL_WP, status, ierr)
     end do

     ! Close the file
     call MPI_FILE_CLOSE(ifile, ierr)

  end if

  ! Stop the timer
  call timing_stop('dump_volC')

  return
end subroutine dump_volume_3D


! ========================================================== !
! Compute all output data on the planes                      !
! ========================================================== !
subroutine dump_volume_compute
  use dump_volume
  use data
  use scalar
  implicit none
  integer :: i, j, k, isc

  !$OMP PARALLEL
  ! Interpolate RHO and SC from n+3/2 to n+1
  !$OMP DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           RHO_s(i,j,k) = RHO(i,j,k)
        end do
     end do
  end do
  !$OMP END DO
  !$OMP DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           RHO(i,j,k) = 0.5_WP*(RHOold(i,j,k) + RHO(i,j,k))
        end do
     end do
  end do
  !$OMP END DO
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SC_s(i,j,k,isc) = SC(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SC(i,j,k,isc) = 0.5_WP*(SCold(i,j,k,isc) + SC(i,j,k,isc))
           end do
        end do
     end do
     !$OMP END DO
  end do
  !$OMP END PARALLEL

  do k=vzmin_,vzmax_
     do j=vymin_,vymax_
        do i=vxmin_,vxmax_
           VOLUME_IO_DATA(5)%var(i,j,k) = RHO(i,j,k) ! Density
           do isc=1,nscalar
              VOLUME_IO_DATA(5+isc)%var(i,j,k) = SC(i,j,k,isc) ! Scalars
           end do
           VOLUME_IO_DATA(6+nscalar)%var(i,j,k) = HR_gp(i,j,k) ! Heat release rate
        end do
     end do
  end do

  !$OMP PARALLEL
  ! Return the saved values of RHO and SC at n+3/2
  !$OMP DO
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           RHO(i,j,k) = RHO_s(i,j,k)
        end do
     end do
  end do
  !$OMP END DO
  do isc=1,nscalar
     !$OMP DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              SC(i,j,k,isc) = SC_s(i,j,k,isc)
           end do
        end do
     end do
     !$OMP END DO
  end do
  !$OMP END PARALLEL

  return
end subroutine dump_volume_compute


! ========================================================== !
! Finalize the communicators                                 !
! ========================================================== !
subroutine dump_volume_finalize
  use dump_volume
  implicit none
  integer :: ierr

  if (in_domain) then
     call MPI_COMM_FREE(VOLUME_IO_COMM, ierr)
  else
     call MPI_COMM_FREE(VOLUME_NON_COMM, ierr)
  end if

  return
end subroutine dump_volume_finalize
                                   

