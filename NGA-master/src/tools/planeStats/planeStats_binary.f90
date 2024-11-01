! ========================================================== !
! Read number of records from binary files                   !
! ========================================================== !
subroutine spectrum_span_binary_read_nrec
  use spectrum_span
  implicit none

  integer :: iunit, ierr, i, j, n
  integer :: ifile, ntime, ny, nvars
  integer :: ny_, nz_, nplanes, nvars_
  real(WP) :: Lz_
  real(WP), dimension(:), pointer :: x, y, x_, y_

  n_records = 0

  ! Loop over the input files
  do ifile=1,n_files
     ! Open the binary file
     call BINARY_FILE_OPEN(iunit, trim(data_files(ifile)), "r", ierr)

     ! Read data sizes
     call BINARY_FILE_READ(iunit, ntime, 1, kind(ntime), ierr)
     call BINARY_FILE_READ(iunit, ny_,   1, kind(ny), ierr)
     call BINARY_FILE_READ(iunit, nz_,   1, kind(nz), ierr)
     call BINARY_FILE_READ(iunit, nplanes, 1, kind(nplanes), ierr)
     call BINARY_FILE_READ(iunit, nvars_,  1, kind(nvars), ierr)
     allocate(y_(ny_))
     allocate(x_(nplanes))
     ! Read y-locations
     do j=1,ny_
        call BINARY_FILE_READ(iunit, y_(j), 1, kind(y_), ierr)
     end do
     ! Read Lz
     call BINARY_FILE_READ(iunit, Lz_, 1, kind(Lz_), ierr)
     ! Read x-locations of planes
     do n=1,nplanes
        call BINARY_FILE_READ(iunit, x_(n), 1, kind(x_), ierr)
     end do

     ! Account for the new data and check consistency
     if (ifile.eq.1) then
        allocate(y(ny_))
        allocate(x(nplanes))
        ny       = ny_
        nz       = nz_
        n_probes = nplanes
        nvars    = nvars_
        Lz       = Lz_
        y        = y_
        x        = x_
     else
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
        if (n_probes.ne.nplanes) then ! This one isn't necessary -- need to indicate which we care about
           print *, 'Inconsistent nplanes in file ', trim(data_files(ifile))
           print *, n_probes, nplanes
           return
        end if
        if (nvars.ne.nvars_) then
           print *, 'Inconsistent nvars in file ', trim(data_files(ifile))
           print *, nvars, nvars_
           return
        end if
        if (Lz.ne.Lz_) then
           print *, 'Inconsistent Lz in file ', trim(data_files(ifile))
           print *, Lz, Lz_
           return
        end if
        do i=1,ny
           if (y_(i).ne.y(i)) then
              print *, 'Inconsistency in y in file ', trim(data_files(ifile))
              return
           end if
        end do
        do i=1,nplanes
           if (x_(i).ne.x(i)) then
              print *, 'Inconsistency in x in file ', trim(data_files(ifile))
              return
           end if
        end do
     end if

     ! Update the number of records
     n_records = n_records + ntime

     ! Close the file
     call BINARY_FILE_CLOSE(iunit, ierr)
     deallocate(x_)
     deallocate(y_)
  end do

  deallocate(x)
  deallocate(y)

  return
end subroutine spectrum_span_binary_read_nrec


! ========================================================== !
! Read input data from binary files                          !
! ========================================================== !
subroutine spectrum_span_binary_read
  use spectrum_span
  implicit none

  integer :: iunit, ierr
  integer :: ifile, it, n, i, j, k, istart, ii, jj
  integer :: nvars, nplanes, ntime, ny
  integer :: ntime_prev, ntime_curr
  real(WP) :: dt, time
  real(WP), dimension(:), pointer :: x, y
  character(len=str_short), dimension(:), pointer :: names
  integer :: iprobe, iplane
  real(WP) :: tmp
  real(WP), dimension(:,:,:), pointer :: PROG, tmpv
  integer, dimension(:,:,:), pointer :: plane_jndx_

  ! First figure out where the progress variable is
  allocate(plane_jndx_(n_records,n_plane_xloc,n_PROG))
  plane_jndx_ = 0

  ntime_curr = 0
  do ifile=1,n_files
     ! Open the file and read its header
     call BINARY_FILE_OPEN(iunit, trim(data_files(ifile)), "r", ierr)

     ! Read data sizes
     call BINARY_FILE_READ(iunit, ntime,   1, kind(ntime), ierr)
     call BINARY_FILE_READ(iunit, ny,      1, kind(ny), ierr)
     call BINARY_FILE_READ(iunit, nz,      1, kind(nz), ierr)
     call BINARY_FILE_READ(iunit, nplanes, 1, kind(nplanes), ierr)
     call BINARY_FILE_READ(iunit, nvars,   1, kind(nvars), ierr)
     allocate(y(ny))
     allocate(x(nplanes))
     allocate(names(nvars))
     ! Read y-locations
     do j=1,ny
        call BINARY_FILE_READ(iunit, y(j), 1, kind(y), ierr)
     end do
     ! Read Lz
     call BINARY_FILE_READ(iunit, Lz, 1, kind(Lz), ierr)
     ! Read x-locations of planes
     do i=1,nplanes
        call BINARY_FILE_READ(iunit, x(i), 1, kind(x), ierr)
     end do
     ! Read variable names
     do n=1,nvars
        call BINARY_FILE_READ(iunit, names(n), str_short, kind(names), ierr)
     end do
     print *, 'file ', ifile, ' read the header, computing PROG'

     ! Find the plane indices to match plane_xloc
     if (ifile.eq.1) then
        do i=1,n_plane_xloc
           do ii=1,nplanes
              if (plane_xloc(i).eq.x(ii)) then
                 plane_indx(i) = ii
              end if
           end do
        end do
     end if

     ! Set the start and end points
     ntime_prev = ntime_curr
     ntime_curr = ntime_curr + ntime

     ! Allocate the arrayw
     allocate(PROG(nplanes,ny,nz))
     allocate(tmpv(nplanes,ny,nz))

     ! Read the data
     do it=ntime_prev+1,ntime_curr
        ! Read time info
        call BINARY_FILE_READ(iunit, dt,   1, kind(dt), ierr)
        call BINARY_FILE_READ(iunit, time, 1, kind(time), ierr)
        print *, it, time
        ! Select the variables
        do n=1,nvars
           call BINARY_FILE_READ(iunit, tmpv, nplanes*ny*nz, kind(tmpv), ierr)
           if (trim(names(n)).eq.'O2') then
              PROG = (1.695440e-01-tmpv)/(1.695440e-01-4.112050e-03)
           end if
!           do k=1,nz
!              do j=1,ny
!                 do i=1,nplanes
!                    call BINARY_FILE_READ(iunit, tmp, 1, kind(tmp), ierr)
!                    if (trim(names(n)).eq.'O2') then ! implement arbitrary definition of PROG
!                       PROG(i,j,k) = (1.695440e-01-tmp)/(1.695440e-01-4.112050e-03)
!                    end if
!                 end do
!              end do
!           end do
        end do
        ! Find the prog
        do i=1,n_plane_xloc
           ii = plane_indx(i)
           do j=ny/2,ny ! needs to be updated!!
              tmp = sum(PROG(ii,j,:))/nz
              do jj=1,n_PROG
                 if (tmp.ge.plane_PROG(jj).and.plane_jndx_(it,i,jj).eq.0) then
                    plane_jndx_(it,i,jj) = j
!                    print *, 'found prog=', tmp, j
                 end if
              end do
           end do
        end do
     end do
    
     ! Close the file
     call BINARY_FILE_CLOSE(iunit, ierr)

     deallocate(x)
     deallocate(y)
     deallocate(names)
     deallocate(PROG)
     deallocate(tmpv)
  end do

  ! Average the progress variable index
  jj = 1
  do i=1,n_plane_xloc
     do j=1,n_PROG
        plane_jndx(jj) = sum(plane_jndx_(:,i,j))/n_records
        jj = jj + 1
     end do
  end do


  ! Might want to just read the data into memory...... but might not be possible for arbitrary file size
  ! Now store the data
  ntime_curr = 0
  do ifile=1,n_files
     ! Open the file and read its header
     call BINARY_FILE_OPEN(iunit, trim(data_files(ifile)), "r", ierr)

     ! Read data sizes
     call BINARY_FILE_READ(iunit, ntime,   1, kind(ntime), ierr)
     call BINARY_FILE_READ(iunit, ny,      1, kind(ny), ierr)
     call BINARY_FILE_READ(iunit, nz,      1, kind(nz), ierr)
     call BINARY_FILE_READ(iunit, nplanes, 1, kind(nplanes), ierr)
     call BINARY_FILE_READ(iunit, nvars,   1, kind(nvars), ierr)
     allocate(y(ny))
     allocate(x(nplanes))
     allocate(names(nvars))
     ! Read y-locations
     do j=1,ny
        call BINARY_FILE_READ(iunit, y(j), 1, kind(y), ierr)
     end do
     ! Read Lz
     call BINARY_FILE_READ(iunit, Lz, 1, kind(Lz), ierr)
     ! Read x-locations of planes
     do i=1,nplanes
        call BINARY_FILE_READ(iunit, x(i), 1, kind(x), ierr)
     end do
     ! Read variable names
     do n=1,nvars
        call BINARY_FILE_READ(iunit, names(n), str_short, kind(names), ierr)
     end do

     ! Set the start and end points
     ntime_prev = ntime_curr
     ntime_curr = ntime_curr + ntime

     ! Read the data
     allocate(tmpv(nplanes,ny,nz))
     ! This is super slow. Either decrease data output frequency or go to MPI.
     do it=ntime_prev+1,ntime_curr
        ! Read time info
        call BINARY_FILE_READ(iunit, dt,   1, kind(dt), ierr)
        call BINARY_FILE_READ(iunit, time, 1, kind(time), ierr)
        print *, 'reading data at itime=', it
        ! Select the variables
        do n=1,nvars
           call BINARY_FILE_READ(iunit, tmpv, nplanes*ny*nz, kind(tmpv), ierr)
           do i=1,n_plane_xloc
              do jj=1+(i-1)*n_PROG,i*n_PROG
                 if (n.eq. 1        ) U(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq. 2        ) V(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq. 3        ) W(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq. 4        ) P(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq. 5+nscalar) RHO(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq. 6+nscalar) VISC(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq. 7+nscalar) dUdx11(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq. 8+nscalar) dUdx12(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq. 9+nscalar) dUdx13(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.10+nscalar) dUdx21(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.11+nscalar) dUdx22(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.12+nscalar) dUdx23(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.13+nscalar) dUdx31(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.14+nscalar) dUdx32(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.15+nscalar) dUdx33(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.16+nscalar) dRHOdx1(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.17+nscalar) dRHOdx2(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.18+nscalar) dRHOdx3(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.19+nscalar) dPdx1(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.20+nscalar) dPdx2(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.21+nscalar) dPdx3(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.22+nscalar) dTAUdx11(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.23+nscalar) dTAUdx12(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.24+nscalar) dTAUdx13(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.25+nscalar) dTAUdx21(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.26+nscalar) dTAUdx22(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.27+nscalar) dTAUdx23(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.28+nscalar) dTAUdx31(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.29+nscalar) dTAUdx32(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.30+nscalar) dTAUdx33(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
                 if (n.eq.31+nscalar) HR(jj,it,:) = tmpv(plane_indx(i),plane_jndx(jj),:)
              end do
           end do
        end do
!!$           do k=1,nz
!!$              jj = 1
!!$              do j=1,ny
!!$                 ii = 1
!!$                 do i=1,nplanes
!!$                    !call BINARY_FILE_READ(iunit, tmp, 1, kind(tmp), ierr)
!!$                    if (j.eq.plane_jndx(jj) .and. i.eq.plane_indx(ii)) then
!!$                       if (n.eq. 1        ) U(jj,it,k) = tmp
!!$                       if (n.eq. 2        ) V(jj,it,k) = tmp
!!$                       if (n.eq. 3        ) W(jj,it,k) = tmp
!!$                       if (n.eq. 4        ) P(jj,it,k) = tmp
!!$                       if (n.eq. 5+nscalar) RHO(jj,it,k) = tmp
!!$                       if (n.eq. 6+nscalar) VISC(jj,it,k) = tmp
!!$                       if (n.eq. 7+nscalar) dUdx11(jj,it,k) = tmp
!!$                       if (n.eq. 8+nscalar) dUdx12(jj,it,k) = tmp
!!$                       if (n.eq. 9+nscalar) dUdx13(jj,it,k) = tmp
!!$                       if (n.eq.10+nscalar) dUdx21(jj,it,k) = tmp
!!$                       if (n.eq.11+nscalar) dUdx22(jj,it,k) = tmp
!!$                       if (n.eq.12+nscalar) dUdx23(jj,it,k) = tmp
!!$                       if (n.eq.13+nscalar) dUdx31(jj,it,k) = tmp
!!$                       if (n.eq.14+nscalar) dUdx32(jj,it,k) = tmp
!!$                       if (n.eq.15+nscalar) dUdx33(jj,it,k) = tmp
!!$                       if (n.eq.16+nscalar) dRHOdx1(jj,it,k) = tmp
!!$                       if (n.eq.17+nscalar) dRHOdx2(jj,it,k) = tmp
!!$                       if (n.eq.18+nscalar) dRHOdx3(jj,it,k) = tmp
!!$                       if (n.eq.19+nscalar) dPdx1(jj,it,k) = tmp
!!$                       if (n.eq.20+nscalar) dPdx2(jj,it,k) = tmp
!!$                       if (n.eq.21+nscalar) dPdx3(jj,it,k) = tmp
!!$                       if (n.eq.22+nscalar) dTAUdx11(jj,it,k) = tmp
!!$                       if (n.eq.23+nscalar) dTAUdx12(jj,it,k) = tmp
!!$                       if (n.eq.24+nscalar) dTAUdx13(jj,it,k) = tmp
!!$                       if (n.eq.25+nscalar) dTAUdx21(jj,it,k) = tmp
!!$                       if (n.eq.26+nscalar) dTAUdx22(jj,it,k) = tmp
!!$                       if (n.eq.27+nscalar) dTAUdx23(jj,it,k) = tmp
!!$                       if (n.eq.28+nscalar) dTAUdx31(jj,it,k) = tmp
!!$                       if (n.eq.29+nscalar) dTAUdx32(jj,it,k) = tmp
!!$                       if (n.eq.30+nscalar) dTAUdx33(jj,it,k) = tmp
!!$                       if (n.eq.31+nscalar) HR(jj,it,k) = tmp
!!$                       if (ii.le.n_plane_xloc) ii = ii + 1
!!$                       if (jj.le.n_plane_xloc*n_PROG) jj = jj + 1
!!$                    end if
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
     end do

     ! Close the file
     call BINARY_FILE_CLOSE(iunit, ierr)

     deallocate(x)
     deallocate(y)
     deallocate(names)
     deallocate(tmpv)
  end do

  return
end subroutine spectrum_span_binary_read
