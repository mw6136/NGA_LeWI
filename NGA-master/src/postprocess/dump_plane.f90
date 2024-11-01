! ========================================================== !
!                       dump_plane.f90                       !
! Dumps all variables (and some gradients) on planes in the  !
!     x-y, x-z and y-z plane at specified locations.         !
!                                                            !
! Author: Jonathan F. MacArt, T. Grenga                      !
! Date:   August 26, 2016                                    !
! ========================================================== !
module dump_plane
  use parallel
  use precision
  use geometry
  use partition
  use masks
  use fileio
  use data
  use scalar
  implicit none

  ! Number of planes (global, local)
  integer :: nplanes, nplanes_
  integer :: pnmin, pnmax, pnmin_, pnmax_
  character(len=2) :: plane_l

  ! Indices of plane locations
  integer :: pxmin, pxmax, pxmin_, pxmax_
  integer :: pymin, pymax, pymin_, pymax_
  integer :: pzmin, pzmax, pzmin_, pzmax_
  integer :: pc1min, pc1max, pc1min_, pc1max_
  integer :: pc2min, pc2max, pc2min_, pc2max_
  integer,  dimension(:), pointer :: index_x,index_xm
  real(WP), dimension(:), pointer :: plane_x
  
  ! Interpolation in x
  real(WP), dimension(:), pointer :: interp_x, interp_xm

  ! MPI IO variables
  integer :: PLANE_IO_COMM, PLANE_NON_COMM
  integer :: irank_p, nsteps_out
  character(len=str_medium) :: filename  
  integer :: PLANE_IO_NVARS
  type(MPI_IO_VAR), dimension(:), pointer :: PLANE_IO_DATA

  real(WP), dimension(:,:,:,:), pointer :: rhs_fc_out

end module dump_plane


! ========================================================== !
! Initialization                                             !
! ========================================================== !
subroutine dump_plane_init
  use dump_plane
  use parser
  use time_info
  implicit none

  integer :: i,j,k,n,var,ierr,isc
  integer, dimension(3) :: gsizes, lsizes, start
  logical :: isdefx, isdefl
  real(WP) :: plane_xmin, plane_xmax, plane_ymin, plane_ymax, plane_zmin, plane_zmax, tmp

  ! Create and start the timer
  call timing_create('dump_plane')
!  call timing_create('dump_planeC')
  call timing_start ('dump_plane')

  ! Read the input file
  call parser_is_defined('Plane orientation',isdefl)
  call parser_is_defined('Plane locations',isdefx)
  if (isdefl .and. isdefx) then
     call parser_getsize('Plane locations',nplanes)
  else
     call die('dump_plane_init: must specify *Plane orientation * and *Plane locations * in the input file')
  end if
  pnmin = 1
  pnmax = nplanes

  ! Allocate arrays
  allocate(plane_x(nplanes))

  ! Read orientation, locations and dimensions
  ! Determine which planes, if any, belong to the current process
  call parser_read('Plane orientation',plane_l)
  call parser_read('Plane locations',plane_x)

  select case(plane_l)
     case ('xy')
       call parser_read('Plane xmin',plane_xmin)
       call parser_read('Plane xmax',plane_xmax)
       call parser_read('Plane ymin',plane_ymin)
       call parser_read('Plane ymax',plane_ymax)
         if (plane_xmin.gt.plane_xmax) &
            call die('dump_plane_init: xmin must be less than or equal to xmax')
         if (plane_ymin.gt.plane_ymax) &
            call die('dump_plane_init: ymin must be less than or equal to ymax')
         if (plane_xmin.lt.x(imin) .or. plane_ymin.lt.y(jmin) .or. &
             plane_xmax.ge.x(imax+1) .or. plane_ymax.ge.y(jmax+1)) &
               call die('dump_plane_init: plane not in domain')

       nplanes_ = 0
       do n=1,nplanes
          if (plane_x(n).lt.z(kmin) .or. plane_x(n).ge.z(kmax+1)) &
             call die('dump_plane_init: plane z-location not in domain')
          if (plane_x(n).ge.z(kmin_) .and. plane_x(n).lt.z(kmax_+1) .and. &
              plane_xmax.ge.x(imin_) .and. plane_xmin.lt.x(imax_+1) .and. &
              plane_ymax.ge.y(jmin_) .and. plane_ymin.lt.y(jmax_+1)) then
             ! Set the local start and end indices
             if (nplanes_.eq.0) pnmin_ = n
             pnmax_   = n
             nplanes_ = nplanes_ + 1
          end if
       end do

     case ('xz')
       call parser_read('Plane xmin',plane_xmin)
       call parser_read('Plane xmax',plane_xmax)
       call parser_read('Plane zmin',plane_zmin)
       call parser_read('Plane zmax',plane_zmax)
         if (plane_xmin.gt.plane_xmax) &
            call die('dump_plane_init: xmin must be less than or equal to xmax')
         if (plane_zmin.gt.plane_zmax) &
            call die('dump_plane_init: zmin must be less than or equal to zmax')
         if (plane_xmin.lt.x(jmin) .or. plane_zmin.lt.z(kmin) .or. &
             plane_xmax.gt.x(jmax+1) .or. plane_zmax.gt.z(kmax+1)) &
               call die('dump_plane_init: plane not in domain')

         if (plane_zmax.eq.z(kmax+1)) plane_zmax = z(kmax)

       nplanes_ = 0
       do n=1,nplanes
          if (plane_x(n).lt.y(jmin) .or. plane_x(n).ge.y(jmax+1)) &
             call die('dump_plane_init: plane y-location not in domain')
          if (plane_x(n).ge.y(jmin_) .and. plane_x(n).lt.y(jmax_+1) .and. &
              plane_xmax.ge.x(imin_) .and. plane_xmin.lt.x(imax_+1) .and. &
              plane_zmax.ge.z(kmin_) .and. plane_zmin.lt.z(kmax_+1)) then
             ! Set the local start and end indices
             if (nplanes_.eq.0) pnmin_ = n
             pnmax_   = n
             nplanes_ = nplanes_ + 1
          end if
       end do

     case ('yz')
       call parser_read('Plane ymin',plane_ymin)
       call parser_read('Plane ymax',plane_ymax)
       call parser_read('Plane zmin',plane_zmin)
       call parser_read('Plane zmax',plane_zmax)
         if (plane_ymin.gt.plane_ymax) &
            call die('dump_plane_init: ymin must be less than or equal to ymax')
         if (plane_zmin.gt.plane_zmax) &
            call die('dump_plane_init: zmin must be less than or equal to zmax')
         if (plane_ymin.lt.y(jmin) .or. plane_zmin.lt.z(kmin) .or. &
             plane_ymax.gt.y(jmax+1) .or. plane_zmax.gt.z(kmax+1)) &
               call die('dump_plane_init: plane not in domain')

       nplanes_ = 0
       do n=1,nplanes
          if (plane_x(n).lt.x(imin) .or. plane_x(n).ge.x(imax+1)) &
             call die('dump_plane_init: plane x-location not in domain')
          if (plane_x(n).ge.x(imin_) .and. plane_x(n).lt.x(imax_+1) .and. &
              plane_ymax.ge.y(jmin_) .and. plane_ymin.lt.y(jmax_+1) .and. &
              plane_zmax.ge.z(kmin_) .and. plane_zmin.lt.z(kmax_+1)) then
             ! Set the local start and end indices
             if (nplanes_.eq.0) pnmin_ = n
             pnmax_   = n
             nplanes_ = nplanes_ + 1
          end if
       end do

  end select

  ! Allocate local arrays
  allocate(interp_x (nplanes)); interp_x  = 0.0_WP
  allocate(interp_xm(nplanes)); interp_xm = 0.0_WP
  allocate(index_x  (nplanes)); index_x   = 0
  allocate(index_xm (nplanes)); index_xm  = 0

  ! Get indices and interpolation
  if (nplanes_.gt.0) then
     ! Indices and interpolation
     index_x   = 0
     index_xm  = 0
     interp_x  = 0.0_WP
     interp_xm = 0.0_WP
     select case(plane_l)
        case ('xy')
           do n=pnmin_,pnmax_
              i = kmin_
              do while (z(i+1).le.plane_x(n))
                 i = i+1
              end do
              index_x(n)  = i
              interp_x(n) = (plane_x(n)-z(i))/(z(i+1)-z(i))
              i = kmino_
              do while (zm(i+1).le.plane_x(n))
                 i = i+1
              end do
              index_xm(n)  = i
              interp_xm(n) = (plane_x(n)-zm(i))/(zm(i+1)-zm(i))
              i = index_x(n)
           end do

           ! Indices in x
           if (plane_xmin.ge.x(imin_)) then ! plane_xmin is on this CPU
              j = imin_
              do while (x(j).lt.plane_xmin)
                 j = j+1
              end do
              pxmin_ = j
           else ! Plane starts below this CPU
              pxmin_  = imin_
           end if
           if (plane_xmax.le.x(imax_)) then ! plane_xmax is on this CPU
              j = imin_
              do while (x(j).lt.plane_xmax)
                 j = j+1
              end do
              pxmax_ = j
           else ! Plane ends above this CPU
              pxmax_  = imax_
           end if
           pc1min_ = pxmin_
           pc1max_ = pxmax_

           ! Indices in y
           if (plane_ymin.ge.y(jmin_)) then ! plane_ymin is on this CPU
              j = jmin_
              do while (y(j).lt.plane_ymin)
                 j = j+1
              end do
              pymin_ = j
           else ! Plane starts below this CPU
              pymin_  = jmin_
           end if
           if (plane_ymax.le.y(jmax_)) then ! plane_ymax is on this CPU
              j = jmin_
              do while (y(j).lt.plane_ymax)
                 j = j+1
              end do
              pymax_ = j
           else ! Plane ends above this CPU
              pymax_  = jmax_
           end if
           pc2min_ = pymin_
           pc2max_ = pymax_

        case ('xz')
           do n=pnmin_,pnmax_
              i = jmin_
              do while (y(i+1).le.plane_x(n))
                 i = i+1
              end do
              index_x(n)  = i
              interp_x(n) = (plane_x(n)-y(i))/(y(i+1)-y(i))
              i = jmino_
              do while (ym(i+1).le.plane_x(n))
                 i = i+1
              end do
              index_xm(n)  = i
              interp_xm(n) = (plane_x(n)-ym(i))/(ym(i+1)-ym(i))
              i = index_x(n)
           end do

           ! Indices in x
           if (plane_xmin.ge.x(imin_)) then ! plane_xmin is on this CPU
              j = imin_
              do while (x(j).lt.plane_xmin)
                 j = j+1
              end do
              pxmin_ = j
           else ! Plane starts below this CPU
              pxmin_  = imin_
           end if
           if (plane_xmax.lt.x(imax_)) then ! plane_xmax is on this CPU
              j = imin_
              do while (x(j).le.plane_xmax)
                 j = j+1
              end do
              pxmax_ = j
           else ! Plane ends above this CPU
              pxmax_  = imax_
           end if
           pc1min_ = pxmin_
           pc1max_ = pxmax_

           ! Indices in z
           if (plane_zmin.ge.z(kmin_)) then ! plane_zmin is on this CPU
              j = kmin_
              do while (z(j).lt.plane_zmin)
                 j = j+1
              end do
              pzmin_ = j
           else ! Plane starts below this CPU
              pzmin_  = kmin_
           end if
           if (plane_zmax.le.z(kmax_)) then ! plane_zmax is on this CPU
              j = kmin_
              do while (z(j).le.plane_zmax)
                 j = j+1
              end do
              pzmax_ = j
           else ! Plane ends above this CPU
              pzmax_  = kmax_
           end if
           pc2min_ = pzmin_
           pc2max_ = pzmax_

        case ('yz')
           do n=pnmin_,pnmax_
              i = imin_
              do while (x(i+1).le.plane_x(n))
                 i = i+1
              end do
              index_x(n)  = i
              interp_x(n) = (plane_x(n)-x(i))/(x(i+1)-x(i))
              i = imino_
              do while (xm(i+1).le.plane_x(n))
                 i = i+1
              end do
              index_xm(n)  = i
              interp_xm(n) = (plane_x(n)-xm(i))/(xm(i+1)-xm(i))
           end do
     
           ! Indices in y
           if (plane_ymin.ge.y(jmin_)) then ! plane_ymin is on this CPU
              j = jmin_
              do while (y(j+1).le.plane_ymin)
                 j = j+1
              end do
              pymin_ = j
           else ! Plane starts below this CPU
              pymin_  = jmin_
           end if
           if (plane_ymax.le.y(jmax_)) then ! plane_ymax is on this CPU
              j = jmin_
              do while (y(j+1).le.plane_ymax)
                 j = j+1
              end do
              pymax_ = j
           else ! Plane ends above this CPU
              pymax_  = jmax_
           end if
           pc1min_ = pymin_
           pc1max_ = pymax_
           !print '(i6,i6,i6,ES12.5,ES12.5)', irank, pymin_, pymax_, y(pymin_), y(pymax_)

           ! Indices in z
           if (plane_zmin.ge.z(kmin_)) then ! plane_zmin is on this CPU
              j = kmin_
              do while (z(j+1).le.plane_zmin)
                 j = j+1
              end do
              pzmin_ = j
           else ! Plane starts below this CPU
              pzmin_  = kmin_
           end if
           if (plane_zmax.le.z(kmax_)) then ! plane_zmax is on this CPU
              j = kmin_
              do while (z(j+1).le.plane_zmax)
                 j = j+1
              end do
              pzmax_ = j
           else ! Plane ends above this CPU
              pzmax_  = kmax_
           end if
           pc2min_ = pzmin_
           pc2max_ = pzmax_

     end select
  end if

  ! Create a communicator for processes with planes to dump
  irank_p = -1
  if (nplanes_.gt.0) then
     call MPI_COMM_SPLIT(MPI_COMM_WORLD, 1, irank, PLANE_IO_COMM, ierr)
     call MPI_COMM_RANK(PLANE_IO_COMM, irank_p, ierr)
  else
     call MPI_COMM_SPLIT(MPI_COMM_WORLD, 0, irank, PLANE_NON_COMM, ierr)
  end if

  ! Get information on all the planes from the communicator
  if (nplanes_.gt.0) then
     ! Get x-index and interpolation of each plane for output
     do n=1,nplanes
        call MPI_ALLREDUCE(MPI_IN_PLACE, index_x(n),  1, MPI_INTEGER, MPI_MAX, PLANE_IO_COMM, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, interp_x(n), 1, MPI_REAL_WP, MPI_MAX, PLANE_IO_COMM, ierr)
     end do     
     ! Get the total nx, ny, nz
     call MPI_ALLREDUCE(pc1min_, pc1min, 1, MPI_INTEGER, MPI_MIN, PLANE_IO_COMM, ierr)
     call MPI_ALLREDUCE(pc1max_, pc1max, 1, MPI_INTEGER, MPI_MAX, PLANE_IO_COMM, ierr)
     call MPI_ALLREDUCE(pc2min_, pc2min, 1, MPI_INTEGER, MPI_MIN, PLANE_IO_COMM, ierr)
     call MPI_ALLREDUCE(pc2max_, pc2max, 1, MPI_INTEGER, MPI_MAX, PLANE_IO_COMM, ierr)
  end if
  
  ! Get the name of the output file
  call parser_read('Plane data file to write',filename)
  if (use_mpiiofs) filename = trim(mpiiofs)//":"//trim(filename)
  
  nsteps_out = 0
  
  ! Set the number of output variables in 3D (2D)
  ! 3x (2x) Velocity components
  ! 1x      Pressure
  ! 1x      Density
  ! 1x      Viscosity (mu)
  ! 9x (4x) Velocity gradients
  ! 3x (2x) Density gradient
  ! 3x (2x) Pressure gradient
  ! 9x (4x) Stress tensor gradients
  ! 4x (3x) nscalar:
  !           1x      SC
  !           3x (2x) SC gradient
!!$  if (nz.eq.1) then
!!$     PLANE_IO_NVARS = 17 + nscalar*3
!!$  else
     PLANE_IO_NVARS = 30 + nscalar*4
!!$  end if

  ! If using finitechem, add further output  variables
  if (use_HR) then
     ! Add to number of output variables:
     ! N_nons+1 : Source terms
     ! nscalar  : Diffusivities
     ! 4*nscalar: 3x (2x) diffusive fluxes, 1x diffusive flux divergence
!!$     if (nz.eq.1) then
!!$        PLANE_IO_NVARS = PLANE_IO_NVARS + N_nons+1 + nscalar + 3*nscalar
!!$     else
        PLANE_IO_NVARS = PLANE_IO_NVARS + N_nons+1 + nscalar + 4*nscalar
!!$     end if

     if (nplanes_.gt.0) then
        ! Allocate arrays in core/data.f90
        allocate(diff_flux(1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nscalar))
        allocate(diff_flux_div(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:nscalar))
     end if
  end if

  if (nplanes_.gt.0) then
     ! Allocate array for output data
     allocate(PLANE_IO_DATA(PLANE_IO_NVARS))
     
     ! Name the output variables
     PLANE_IO_DATA(1)%name = 'U'
     PLANE_IO_DATA(2)%name = 'V'
     PLANE_IO_DATA(3)%name = 'W'
     PLANE_IO_DATA(4)%name = 'P'
     do isc=1,nscalar
        PLANE_IO_DATA(4+isc)%name = trim(SC_name(isc))
     end do
     PLANE_IO_DATA( 5+nscalar)%name = 'RHO'
     PLANE_IO_DATA( 6+nscalar)%name = 'VISC'
     ! Velocity gradients
     PLANE_IO_DATA( 7+nscalar)%name = 'dUdx11'
     PLANE_IO_DATA( 8+nscalar)%name = 'dUdx12'
     PLANE_IO_DATA( 9+nscalar)%name = 'dUdx13'
     PLANE_IO_DATA(10+nscalar)%name = 'dUdx21'
     PLANE_IO_DATA(11+nscalar)%name = 'dUdx22'
     PLANE_IO_DATA(12+nscalar)%name = 'dUdx23'
     PLANE_IO_DATA(13+nscalar)%name = 'dUdx31'
     PLANE_IO_DATA(14+nscalar)%name = 'dUdx32'
     PLANE_IO_DATA(15+nscalar)%name = 'dUdx33'
     ! Density gradient
     PLANE_IO_DATA(16+nscalar)%name = 'dRHOdx1'
     PLANE_IO_DATA(17+nscalar)%name = 'dRHOdx2'
     PLANE_IO_DATA(18+nscalar)%name = 'dRHOdx3'
     ! Pressure gradient
     PLANE_IO_DATA(19+nscalar)%name = 'dPdx1'
     PLANE_IO_DATA(20+nscalar)%name = 'dPdx2'
     PLANE_IO_DATA(21+nscalar)%name = 'dPdx3'
     ! Stress tensor gradients
     PLANE_IO_DATA(22+nscalar)%name = 'dTAUdx11'
     PLANE_IO_DATA(23+nscalar)%name = 'dTAUdx12'
     PLANE_IO_DATA(24+nscalar)%name = 'dTAUdx13'
     PLANE_IO_DATA(25+nscalar)%name = 'dTAUdx21'
     PLANE_IO_DATA(26+nscalar)%name = 'dTAUdx22'
     PLANE_IO_DATA(27+nscalar)%name = 'dTAUdx23'
     PLANE_IO_DATA(28+nscalar)%name = 'dTAUdx31'
     PLANE_IO_DATA(29+nscalar)%name = 'dTAUdx32'
     PLANE_IO_DATA(30+nscalar)%name = 'dTAUdx33'
     i = 30+nscalar
     ! Scalar gradients
     do isc=1,nscalar
        PLANE_IO_DATA(i+1+3*(isc-1))%name = 'dx1_'//trim(SC_name(isc))
        PLANE_IO_DATA(i+2+3*(isc-1))%name = 'dx2_'//trim(SC_name(isc))
        PLANE_IO_DATA(i+3+3*(isc-1))%name = 'dx3_'//trim(SC_name(isc))
     end do
     i = i+3*nscalar

     ! Further variables added when using finitechem
     if (use_HR) then
        do isc=1,N_nons+1
           ! Chemical source terms
           ! Assumes ordering of species [Y_k, T]
           PLANE_IO_DATA(i+isc)%name = 'src_'//trim(SC_name(isc))
        end do
        i = i+N_nons+1

        do isc=1,nscalar
           ! Species diffusivities
           PLANE_IO_DATA(i+isc)%name = 'DIFF_'//trim(SC_name(isc))
        end do
        i = i+nscalar

        do isc=1,nscalar
           ! Species diffusive fluxes
           PLANE_IO_DATA(i+1+3*(isc-1))%name = 'q1_'//trim(SC_name(isc))
           PLANE_IO_DATA(i+2+3*(isc-1))%name = 'q2_'//trim(SC_name(isc))
           PLANE_IO_DATA(i+3+3*(isc-1))%name = 'q3_'//trim(SC_name(isc))
        end do
        i = i+3*nscalar

        do isc=1,nscalar
           ! Divergence of species diffusive flux
           PLANE_IO_DATA(i+isc)%name = 'dq_'//trim(SC_name(isc))
        end do

        ! Allocate finitechem arrays
        allocate(rhs_fc_out(pnmin_:pnmax_,pc1min_:pc1max_,pc2min_:pc2max_,1:N_nons+1))
     end if

     ! Allocate the variables
     do var=1,PLANE_IO_NVARS
        allocate(PLANE_IO_DATA(var)%var(pnmin_:pnmax_,pc1min_:pc1max_,pc2min_:pc2max_))
        PLANE_IO_DATA(var)%var = 0.0_WP
     end do
     
     ! Global sizes
     gsizes(1) = nplanes
     gsizes(2) = pc1max  - pc1min  + 1
     gsizes(3) = pc2max  - pc2min  + 1

     ! Local sizes
     lsizes(1) = nplanes_
     lsizes(2) = pc1max_ - pc1min_ + 1
     lsizes(3) = pc2max_ - pc2min_ + 1

     ! Starting points
     start(1) = pnmin_-pnmin
     start(2) = pc1min_-pc1min
     start(3) = pc2min_-pc2min

     ! Define the view for each variable
     do var=1,PLANE_IO_NVARS
        call MPI_TYPE_CREATE_SUBARRAY(3, gsizes, lsizes, start, &
             MPI_ORDER_FORTRAN, MPI_REAL_WP, PLANE_IO_DATA(var)%view, ierr)
        call MPI_TYPE_COMMIT(PLANE_IO_DATA(var)%view, ierr)
     end do

  end if

  ! Stop the timer
  call timing_stop('dump_plane')

  call monitor_log("dump_plane_init COMPLETED")

  return
end subroutine dump_plane_init


! ========================================================== !
! Dump all data on the planes                                !
! ========================================================== !
subroutine dump_plane_3D
  use dump_plane
  implicit none

  integer :: n, j
  integer :: ifile, ierr, var, data_size
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer, dimension(6) :: dims
  integer(kind=MPI_Offset_kind) :: disp
  integer(kind=MPI_Offset_kind) :: nd1_MOK, nd2_MOK, np_MOK
  integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK, str_MOK
  integer(kind=MPI_Offset_kind) :: NVARS_MOK, NTIME_MOK
  integer :: pd1, pd2, pd1_, pd2_
  character(len=str_medium) :: buffer
  logical :: file_is_there
  real(WP) :: tmp

  ! Start the timer
  call timing_start('dump_plane')

  ! Increment the number of output steps
  nsteps_out = nsteps_out + 1

  if (nplanes_.gt.0) then
     ! Resize some integers so MPI can write even the biggest files
     np_MOK    = int(nplanes,       MPI_Offset_kind)
     WP_MOK    = int(WP,            MPI_Offset_kind)
     str_MOK   = int(str_short,     MPI_Offset_kind)
     NVARS_MOK = int(PLANE_IO_NVARS,MPI_Offset_kind)
     NTIME_MOK = int(nsteps_out,    MPI_Offset_kind)
     
     pd1     = pc1max  - pc1min  + 1
     pd2     = pc2max  - pc2min  + 1
     pd1_    = pc1max_ - pc1min_ + 1
     pd2_    = pc2max_ - pc2min_ + 1
     nd1_MOK = int(pd1, MPI_Offset_kind)
     nd2_MOK = int(pd2, MPI_Offset_kind)
     
     select case(plane_l)
     case ('xy')
        dims(5) = 1
     case ('xz')
        dims(5) = 2
     case ('yz')
        dims(5) = 3
     end select

     ! Interpolate data in x and compute derivatives
     call dump_plane_compute

     ! Open the file to write
     call MPI_FILE_OPEN(PLANE_IO_COMM, trim(filename), MPI_MODE_WRONLY+MPI_MODE_CREATE, mpi_info, ifile, ierr)

     ! Write the header
     if (nsteps_out.eq.1) then
        ! First step -- write the full header
        if (irank_p.eq.0) then
           ! Write dimensions
           dims(1) = nsteps_out
           dims(2) = pd1
           dims(3) = pd2
           dims(4) = nplanes
           dims(6) = PLANE_IO_NVARS
           call MPI_FILE_WRITE(ifile, dims, 6, MPI_INTEGER, status, ierr)
           ! Write locations directions
           select case(plane_l)
              case ('xy')
                 do j=pc1min,pc1max
                    call MPI_FILE_WRITE(ifile, x(j), 1, MPI_REAL_WP, status, ierr)
                 end do
                 do j=pc2min,pc2max
                    call MPI_FILE_WRITE(ifile, y(j), 1, MPI_REAL_WP, status, ierr)
                 end do
                 do n=1,nplanes
                    tmp = interp_x(n)*(z(index_x(n)+1)-z(index_x(n))) + z(index_x(n))
                    call MPI_FILE_WRITE(ifile, tmp, 1, MPI_REAL_WP, status, ierr)
                 end do
              case ('xz')
                 do j=pc1min,pc1max
                    call MPI_FILE_WRITE(ifile, x(j), 1, MPI_REAL_WP, status, ierr)
                 end do
                 do j=pc2min,pc2max
                    call MPI_FILE_WRITE(ifile, z(j), 1, MPI_REAL_WP, status, ierr)
                 end do
                 do n=1,nplanes
                    tmp = interp_x(n)*(y(index_x(n)+1)-y(index_x(n))) + y(index_x(n))
                    call MPI_FILE_WRITE(ifile, tmp, 1, MPI_REAL_WP, status, ierr)
                 end do
              case ('yz')
                 do j=pc1min,pc1max
                    call MPI_FILE_WRITE(ifile, y(j), 1, MPI_REAL_WP, status, ierr)
                 end do
                 do j=pc2min,pc2max
                    call MPI_FILE_WRITE(ifile, z(j), 1, MPI_REAL_WP, status, ierr)
                 end do
                 do n=1,nplanes
                    tmp = interp_x(n)*(x(index_x(n)+1)-x(index_x(n))) + x(index_x(n))
                    call MPI_FILE_WRITE(ifile, tmp, 1, MPI_REAL_WP, status, ierr)
                 end do
           end select
           ! Write variable names
           do var=1,PLANE_IO_NVARS
              call MPI_FILE_WRITE(ifile,PLANE_IO_DATA(var)%name,str_short,MPI_CHARACTER,status,ierr)
           end do
           ! Write time info
           call MPI_FILE_WRITE(ifile,dt,  1,MPI_REAL_WP,status,ierr)
           call MPI_FILE_WRITE(ifile,time,1,MPI_REAL_WP,status,ierr)
        end if
     else
        ! Subsequent steps -- just update the step count and append dt and time
        disp = 0
        call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, MPI_INTEGER, &
             "native", mpi_info, ierr)
        if (irank_p.eq.0) then
           ! Update the step count
           call MPI_FILE_WRITE(ifile, nsteps_out, 1, MPI_INTEGER, status, ierr)
        end if
        disp = 6*4 + (nd1_MOK+nd2_MOK+np_MOK)*WP_MOK + str_MOK*NVARS_MOK + 2*(NTIME_MOK-1)*WP_MOK &
             + np_MOK*nd1_MOK*nd2_MOK*WP_MOK*NVARS_MOK*(NTIME_MOK-1)
        call MPI_FILE_SET_VIEW(ifile, disp, MPI_INTEGER, MPI_INTEGER, &
             "native", mpi_info, ierr)
        if (irank_p.eq.0) then
           ! Write time info
           call MPI_FILE_WRITE(ifile,dt,  1,MPI_REAL_WP,status,ierr)
           call MPI_FILE_WRITE(ifile,time,1,MPI_REAL_WP,status,ierr)
        end if
     end if

     ! Size of local arrays
     data_size = nplanes_*pd1_*pd2_

     ! Write the data for each variable
     do var=1,PLANE_IO_NVARS
        var_MOK = int(var,MPI_Offset_kind)
        disp = 6*4 + (nd1_MOK+nd2_MOK+np_MOK)*WP_MOK + str_MOK*NVARS_MOK + 2*NTIME_MOK*WP_MOK &
             + np_MOK*nd1_MOK*nd2_MOK*WP_MOK*NVARS_MOK*(NTIME_MOK-1) + np_MOK*nd1_MOK*nd2_MOK*WP_MOK*(var_MOK-1)
        call MPI_FILE_SET_VIEW(ifile, disp, MPI_REAL_WP, PLANE_IO_DATA(var)%view, &
             "native", mpi_info, ierr)
        call MPI_FILE_WRITE_ALL(ifile, PLANE_IO_DATA(var)%var, data_size, &
             MPI_REAL_WP, status, ierr)
     end do

     ! Close the file
     call MPI_FILE_CLOSE(ifile, ierr)

  end if

  ! Stop the timer
  call timing_stop('dump_plane')

  return
end subroutine dump_plane_3D


! ========================================================== !
! Compute all output data on the planes                      !
! ========================================================== !
subroutine dump_plane_compute
  use dump_plane
  use interpolate
  use strainrate
  use metric_velocity_visc
  implicit none

  integer :: m, n, isc, i, im, j, k, ivar
  integer :: ii, jj, kk
  real(WP) :: wx1,  wx2
  real(WP) :: wx1m, wx2m
  real(WP), dimension(3,3,2) :: dUdx
  real(WP), dimension(2) :: dSCdx
  !real(WP), dimension(3,2) :: dSCdxs
  !real(WP), dimension(3,3,0:1,-stv2:stv2) :: dUdxs
  !real(WP), dimension(2,-stv1:stv2) :: VISCi
  real(WP), dimension(-st1:st2) :: tmpv
  real(WP), dimension(-stv2:stv2) :: FX, FY, FZ
  real(WP), dimension(N_nons+1) :: sol

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

  if (use_HR) then
     diff_output = .true.
     call finitechem_diffusivity
     call finitechem_viscosity
     call finitechem_diffusion
     diff_output = .false.

!!$     rhs_fc = 0.0_WP
     rhs_fc_out = 0.0_WP
     ! Get rhs from finitechem
     do k=pc2min_,pc2max_!kmino_,kmaxo_
        do j=pc1min_,pc1max_!jmino_,jmaxo_
           do n=pnmin_,pnmax_
              im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
              ! This is correct for yz-planes. Need to update!
              do isc=1,N_nons
                 sol(isc) = wx1m*SC(im,j,k,isc_sc-1+isc) + wx2m*SC(im+1,j,k,isc_sc-1+isc)
              end do
              sol(N_nons+1) = wx1m*SC(im,j,k,isc_T) + wx2m*SC(im+1,j,k,isc_T)
              call finitechem_mono_rhs(sol,rhs_fc_out(n,j,k,1:N_nons+1))
!!$           do i=imino_,imaxo_
!!$              if (SC(i,j,k,isc_T).gt.0.0_WP) then
!!$                 if (.not.compressible .and. xper.ne.1) then
!!$                    sol(1:N_nons) = SC(i,j,k,isc_sc:isc_sc-1+N_nons)
!!$                    sol(N_nons+1) = SC(i,j,k,isc_T)
!!$                 else
!!$                    sol(1:N_nons) = 0.5_WP*(RHO(i,j,k)+RHOold(i,j,k))*SC(i,j,k,isc_sc:isc_sc-1+N_nons)
!!$                    sol(N_nons+1) = 0.5_WP*(RHO(i,j,k)+RHOold(i,j,k))*SC(i,j,k,isc_T)
!!$                 end if
!!$                 call finitechem_mono_rhs(sol,rhs_fc(i,j,k,1:N_nons+1))
!!$              end if
           end do
        end do
     end do
  end if

  ! Compute local values on the planes
  select case(plane_l)
  case ('xy')
     ! x->j, y->k, z->i
     do n=pnmin_,pnmax_
        i  = index_x(n);  wx2  = interp_x(n);  wx1  = 1.0_WP-wx2
        im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
        do j=pc1min_,pc1max_
           do k=pc2min_,pc2max_
              ! All values below are at cell centers
              PLANE_IO_DATA(1)%var(n,j,k) = wx1m*Ui(j,k,im) + wx2m*Ui(j,k,im+1) ! Velocity - U
              PLANE_IO_DATA(2)%var(n,j,k) = wx1m*Vi(j,k,im) + wx2m*Vi(j,k,im+1) ! Velocity - V
              PLANE_IO_DATA(3)%var(n,j,k) = wx1 *W (j,k,i ) + wx2 *W (j,k,i +1) ! Velocity - W
              PLANE_IO_DATA(4)%var(n,j,k) = wx1m*P (j,k,im) + wx2m*P (j,k,im+1) ! Pressure
              do isc=1,nscalar
                 PLANE_IO_DATA(4+isc)%var(n,j,k) = wx1m*SC(j,k,im,isc) + wx2m*SC(j,k,im+1,isc) ! Scalars
              end do
              PLANE_IO_DATA(5+nscalar)%var(n,j,k) = wx1m*RHO (j,k,im) + wx2m*RHO (j,k,im+1) ! Density
              PLANE_IO_DATA(6+nscalar)%var(n,j,k) = wx1m*VISC(j,k,im) + wx2m*VISC(j,k,im+1) ! Viscosity
              call vel_grad_local(j,k,i  , dUdx(:,:,1))
              call vel_grad_local(j,k,i+1, dUdx(:,:,2))
              m = 7
              do ii=1,3
                 do jj=1,3
                    PLANE_IO_DATA(m+nscalar)%var(n,j,k) = wx1m*dUdx(ii,jj,1) + wx2m*dUdx(ii,jj,2) ! Velocity gradients
                    m = m + 1
                 end do
              end do
              ! Density gradient x
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_x(j+jj,k,:)*RHO(j+jj-st2:j+jj+st1,k,im  ))
              end do
              dSCdx(1) = sum(interp_u_xm(j,k,:)*tmpv)
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_x(j+jj,k,:)*RHO(j+jj-st2:j+jj+st1,k,im+1))
              end do
              dSCdx(2) = sum(interp_u_xm(j,k,:)*tmpv)
              PLANE_IO_DATA(16+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              ! Density gradient y
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_y(j,k+kk,:)*RHO(j,k+kk-st2:k+kk+st1,im  ))
              end do
              dSCdx(1) = sum(interp_v_ym(j,k,:)*tmpv)
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_y(j,k+kk,:)*RHO(j,k+kk-st2:k+kk+st1,im+1))
              end do
              dSCdx(2) = sum(interp_v_ym(j,k,:)*tmpv)
              PLANE_IO_DATA(17+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              ! Density gradient z
              dSCdx(1) = sum(grad_z(j,k,:)*RHO(j,k,i  -st2:i  +st1))
              dSCdx(2) = sum(grad_z(j,k,:)*RHO(j,k,i+1-st2:i+1+st1))
              PLANE_IO_DATA(18+nscalar)%var(n,j,k) = wx1 *dSCdx(1) + wx2 *dSCdx(2)
              ! Pressure gradient x
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_x(j+jj,k,:)*P(j+jj-st2:j+jj+st1,k,im  ))
              end do
              dSCdx(1) = sum(interp_u_xm(j,k,:)*tmpv)
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_x(j+jj,k,:)*P(j+jj-st2:j+jj+st1,k,im+1))
              end do
              dSCdx(2) = sum(interp_u_xm(j,k,:)*tmpv)
              PLANE_IO_DATA(19+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              ! Pressure gradient y
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_y(j,k+kk,:)*P(j,k+kk-st2:k+kk+st1,im  ))
              end do
              dSCdx(1) = sum(interp_v_ym(j,k,:)*tmpv)
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_y(j,k+kk,:)*P(j,k+kk-st2:k+kk+st1,im+1))
              end do
              dSCdx(2) = sum(interp_v_ym(j,k,:)*tmpv)
              PLANE_IO_DATA(20+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              ! Pressure gradient z
              dSCdx(1) = sum(grad_z(j,k,:)*P(j,k,i  -st2:i  +st1))
              dSCdx(2) = sum(grad_z(j,k,:)*P(j,k,i+1-st2:i+1+st1))
              PLANE_IO_DATA(21+nscalar)%var(n,j,k) = wx1 *dSCdx(1) + wx2 *dSCdx(2)
           end do
        end do
     end do

     ivar = 30+nscalar
     do isc=1,nscalar
        do j=pc1min_,pc1max_
           do k=pc2min_,pc2max_
              do n=pnmin_,pnmax_
                 i  = index_x(n);  wx2  = interp_x(n);  wx1  = 1.0_WP-wx2
                 im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                 ! Scalar gradients
                 ! dx
                 do jj=-st1,st2
                    tmpv(jj) = sum(grad_x(j+jj,k,:)*SC(j+jj-st2:j+jj+st1,k,im  ,isc))
                 end do
                 dSCdx(1) = sum(interp_u_xm(j,k,:)*tmpv)
                 do jj=-st1,st2
                    tmpv(jj) = sum(grad_x(j+jj,k,:)*SC(j+jj-st2:j+jj+st1,k,im+1,isc))
                 end do
                 dSCdx(2) = sum(interp_u_xm(j,k,:)*tmpv)
                 PLANE_IO_DATA(ivar+1+3*(isc-1))%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
                 ! dy
                 do kk=-st1,st2
                    tmpv(kk) = sum(grad_y(j,k+kk,:)*SC(j,k+kk-st2:k+kk+st1,im  ,isc))
                 end do
                 dSCdx(1) = sum(interp_v_ym(j,k,:)*tmpv)
                 do kk=-st1,st2
                    tmpv(kk) = sum(grad_y(j,k+kk,:)*SC(j,k+kk-st2:k+kk+st1,im+1,isc))
                 end do
                 dSCdx(2) = sum(interp_v_ym(j,k,:)*tmpv)
                 PLANE_IO_DATA(ivar+2+3*(isc-1))%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
                 ! dz
                 dSCdx(1) = sum(grad_z(j,k,:)*SC(j,k,i  -st2:i  +st1,isc))
                 dSCdx(2) = sum(grad_z(j,k,:)*SC(j,k,i+1-st2:i+1+st1,isc))
                 PLANE_IO_DATA(ivar+3+3*(isc-1))%var(n,j,k) = wx1 *dSCdx(1) + wx2 *dSCdx(2)
              end do
           end do
        end do
     end do
     ivar = ivar+3*nscalar

     ! Finitechem outputs
     if (use_HR) then
        do isc=1,N_nons+1
           do k=pc2min_,pc2max_
              do j=pc1min_,pc1max_
                 do n=pnmin_,pnmax_
                    im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                    ! Scalar source terms
                    PLANE_IO_DATA(ivar+isc)%var(n,j,k) = wx1m*rhs_fc(j,k,im,isc) + wx2m*rhs_fc(j,k,im+1,isc)
                 end do
              end do
           end do
        end do
        ivar = ivar+N_nons+1

        do isc=1,nscalar
           do k=pc2min_,pc2max_
              do j=pc1min_,pc1max_
                 do n=pnmin_,pnmax_
                    im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                    ! Species diffusivities
                    PLANE_IO_DATA(ivar+isc)%var(n,j,k) = wx1m*DIFF(j,k,im,isc) + wx2m*DIFF(j,k,im+1,isc)
                 end do
              end do
           end do
        end do
        ivar = ivar+nscalar

        do isc=1,nscalar
           do k=pc2min_,pc2max_
              do j=pc1min_,pc1max_
                 do n=pnmin_,pnmax_
                    i  = index_x(n);  wx2  = interp_x(n);  wx1  = 1.0_WP-wx2
                    im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                    ! Species diffusive fluxes
                    dSCdx(1) = sum(interp_u_xm(j,k,:)*diff_flux(1,j-st1:j+st2,k,im  ,isc))
                    dSCdx(2) = sum(interp_u_xm(j,k,:)*diff_flux(1,j-st1:j+st2,k,im+1,isc))
                    PLANE_IO_DATA(ivar+1+3*(isc-1))%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
                    dSCdx(1) = sum(interp_v_ym(j,k,:)*diff_flux(2,j,k-st1:k+st2,im  ,isc))
                    dSCdx(2) = sum(interp_v_ym(j,k,:)*diff_flux(2,j,k-st1:k+st2,im+1,isc))
                    PLANE_IO_DATA(ivar+2+3*(isc-1))%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
                    PLANE_IO_DATA(ivar+3+3*(isc-1))%var(n,j,k) = wx1*diff_flux(3,j,k,i,isc) + wx2*diff_flux(3,j,k,i+1,isc)
                 end do
              end do
           end do
        end do
        ivar = ivar+3*nscalar

        do isc=1,nscalar
           do k=pc2min_,pc2max_
              do j=pc1min_,pc1max_
                 do n=pnmin_,pnmax_
                    im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                    ! Species diffusive flux divergence
                    PLANE_IO_DATA(ivar+isc)%var(n,j,k) = wx1m*diff_flux_div(j,k,im,isc) + wx2m*diff_flux_div(j,k,im+1,isc)
                 end do
              end do
           end do
        end do
        ivar = ivar+nscalar
     end if

     !! NEED STRESS TENSOR GRADIENTS

  case ('xz')
     ! x->j, y->, z->k
     do n=pnmin_,pnmax_
        i  = index_x(n);  wx2  = interp_x(n);  wx1  = 1.0_WP-wx2
        im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
        do j=pc1min_,pc1max_
           do k=pc2min_,pc2max_
              ! All values below are at cell centers
              PLANE_IO_DATA(1)%var(n,j,k) = wx1m*Ui(j,im,k) + wx2m*Ui(j,im+1,k) ! Velocity - U
              PLANE_IO_DATA(2)%var(n,j,k) = wx1 *V (j,i ,k) + wx2 *V (j,i +1,k) ! Velocity - V
              PLANE_IO_DATA(3)%var(n,j,k) = wx1m*Wi(j,im,k) + wx2m*Wi(j,im+1,k) ! Velocity - W
              PLANE_IO_DATA(4)%var(n,j,k) = wx1m*P (j,im,k) + wx2m*P (j,im+1,k) ! Pressure
              do isc=1,nscalar
                 PLANE_IO_DATA(4+isc)%var(n,j,k) = wx1m*SC(j,im,k,isc) + wx2m*SC(j,im+1,k,isc) ! Scalars
              end do
              PLANE_IO_DATA(5+nscalar)%var(n,j,k) = wx1m*RHO (j,im,k) + wx2m*RHO (j,im+1,k) ! Density
              PLANE_IO_DATA(6+nscalar)%var(n,j,k) = wx1m*VISC(j,im,k) + wx2m*VISC(j,im+1,k) ! Viscosity
              call vel_grad_local(j,i  ,k, dUdx(:,:,1))
              call vel_grad_local(j,i+1,k, dUdx(:,:,2))
              m = 7
              do ii=1,3
                 do jj=1,3
                    PLANE_IO_DATA(m+nscalar)%var(n,j,k) = wx1m*dUdx(ii,jj,1) + wx2m*dUdx(ii,jj,2) ! Velocity gradients
                    m = m + 1
                 end do
              end do
              ! Density gradient x
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_x(j+jj,im  ,:)*RHO(j+jj-st2:j+jj+st1,im  ,k))
              end do
              dSCdx(1) = sum(interp_u_xm(j,im  ,:)*tmpv)
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_x(j+jj,im+1,:)*RHO(j+jj-st2:j+jj+st1,im+1,k))
              end do
              dSCdx(2) = sum(interp_u_xm(j,im+1,:)*tmpv)
              PLANE_IO_DATA(16+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              ! Density gradient y
              dSCdx(1) = sum(grad_y(j,i  ,:)*RHO(j,i  -st2:i  +st1,k))
              dSCdx(2) = sum(grad_y(j,i+1,:)*RHO(j,i+1-st2:i+1+st1,k))
              PLANE_IO_DATA(17+nscalar)%var(n,j,k) = wx1 *dSCdx(1) + wx2 *dSCdx(2)
              ! Density gradient z
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_z(j,im  ,:)*RHO(j,im  ,k+kk-st2:k+kk+st1))
              end do
              dSCdx(1) = sum(interp_w_zm(j,im  ,:)*tmpv)
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_z(j,im+1,:)*RHO(j,im+1,k+kk-st2:k+kk+st1))
              end do
              dSCdx(2) = sum(interp_w_zm(j,im+1,:)*tmpv)
              PLANE_IO_DATA(18+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              ! Pressure gradient x
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_x(j+jj,im  ,:)*P(j+jj-st2:j+jj+st1,im  ,k))
              end do
              dSCdx(1) = sum(interp_u_xm(j,im  ,:)*tmpv)
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_x(j+jj,im+1,:)*P(j+jj-st2:j+jj+st1,im+1,k))
              end do
              dSCdx(2) = sum(interp_u_xm(j,im+1,:)*tmpv)
              PLANE_IO_DATA(19+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              ! Pressure gradient y
              dSCdx(1) = sum(grad_y(j,i  ,:)*P(j,i  -st2:i  +st1,k))
              dSCdx(2) = sum(grad_y(j,i+1,:)*P(j,i+1-st2:i+1+st1,k))
              PLANE_IO_DATA(20+nscalar)%var(n,j,k) = wx1 *dSCdx(1) + wx2 *dSCdx(2)
              ! Pressure gradient z
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_z(j,im  ,:)*P(j,im  ,k+kk-st2:k+kk+st1))
              end do
              dSCdx(1) = sum(interp_w_zm(j,im  ,:)*tmpv)
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_z(j,im+1,:)*P(j,im+1,k+kk-st2:k+kk+st1))
              end do
              dSCdx(2) = sum(interp_w_zm(j,im+1,:)*tmpv)
              PLANE_IO_DATA(21+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              if (use_HR) then
                 PLANE_IO_DATA(31+nscalar)%var(n,j,k) = wx1m*HR_gp(j,im,k) + wx2m*HR_gp(j,im+1,k) ! Heat release rate
              end if
           end do
        end do
     end do
     
  case ('yz')
     ! x->i, y->j, z->k
     ! All values below are at cell centers
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              i  = index_x(n);  wx2  = interp_x(n);  wx1  = 1.0_WP-wx2
              PLANE_IO_DATA(1)%var(n,j,k) = wx1 *U (i ,j,k) + wx2 *U (i +1,j,k) ! Velocity - U
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
              PLANE_IO_DATA(2)%var(n,j,k) = wx1m*Vi(im,j,k) + wx2m*Vi(im+1,j,k) ! Velocity - V
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
              PLANE_IO_DATA(3)%var(n,j,k) = wx1m*Wi(im,j,k) + wx2m*Wi(im+1,j,k) ! Velocity - W
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
              PLANE_IO_DATA(4)%var(n,j,k) = wx1m*P (im,j,k) + wx2m*P (im+1,j,k) ! Pressure
           end do
        end do
     end do
     do isc=1,nscalar
        do k=pc2min_,pc2max_
           do j=pc1min_,pc1max_
              do n=pnmin_,pnmax_
                 im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                 PLANE_IO_DATA(4+isc)%var(n,j,k) = wx1m*SC(im,j,k,isc) + wx2m*SC(im+1,j,k,isc) ! Scalars
              end do
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
              PLANE_IO_DATA(5+nscalar)%var(n,j,k) = wx1m*RHO (im,j,k) + wx2m*RHO (im+1,j,k) ! Density
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
              PLANE_IO_DATA(6+nscalar)%var(n,j,k) = wx1m*VISC(im,j,k) + wx2m*VISC(im+1,j,k) ! Viscosity
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
              call vel_grad_local(im  ,j,k, dUdx(:,:,1))
              call vel_grad_local(im+1,j,k, dUdx(:,:,2))
              m = 7
              do ii=1,3
                 do jj=1,3
                    PLANE_IO_DATA(m+nscalar)%var(n,j,k) = wx1m*dUdx(ii,jj,1) + wx2m*dUdx(ii,jj,2) ! Velocity gradients
                    m = m + 1
                 end do
              end do
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              i  = index_x(n);  wx2  = interp_x(n);  wx1  = 1.0_WP-wx2
              ! Density gradient x
              dSCdx(1) = sum(grad_x(i  ,j,:)*RHO(i  -st2:i  +st1,j,k))
              dSCdx(2) = sum(grad_x(i+1,j,:)*RHO(i+1-st2:i+1+st1,j,k))
              PLANE_IO_DATA(16+nscalar)%var(n,j,k) = wx1 *dSCdx(1) + wx2 *dSCdx(2)
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
              ! Density gradient y
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_y(im  ,j+jj,:)*RHO(im  ,j+jj-st2:j+jj+st1,k))
              end do
              dSCdx(1) = sum(interp_v_ym(im  ,j,:)*tmpv)
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_y(im+1,j+jj,:)*RHO(im+1,j+jj-st2:j+jj+st1,k))
              end do
              dSCdx(2) = sum(interp_v_ym(im+1,j,:)*tmpv)
              PLANE_IO_DATA(17+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
              ! Density gradient z
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_z(im  ,j,:)*RHO(im  ,j,k+kk-st2:k+kk+st1))
              end do
              dSCdx(1) = sum(interp_w_zm(im  ,j,:)*tmpv)
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_z(im+1,j,:)*RHO(im+1,j,k+kk-st2:k+kk+st1))
              end do
              dSCdx(2) = sum(interp_w_zm(im+1,j,:)*tmpv)
              PLANE_IO_DATA(18+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              i  = index_x(n);  wx2  = interp_x(n);  wx1  = 1.0_WP-wx2
              ! Pressure gradient x
              dSCdx(1) = sum(grad_x(i  ,j,:)*P(i  -st2:i  +st1,j,k))
              dSCdx(2) = sum(grad_x(i+1,j,:)*P(i+1-st2:i+1+st1,j,k))
              PLANE_IO_DATA(19+nscalar)%var(n,j,k) = wx1 *dSCdx(1) + wx2 *dSCdx(2)
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
              ! Pressure gradient y
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_y(im  ,j+jj,:)*P(im  ,j+jj-st2:j+jj+st1,k))
              end do
              dSCdx(1) = sum(interp_v_ym(im  ,j,:)*tmpv)
              do jj=-st1,st2
                 tmpv(jj) = sum(grad_y(im+1,j+jj,:)*P(im+1,j+jj-st2:j+jj+st1,k))
              end do
              dSCdx(2) = sum(interp_v_ym(im+1,j,:)*tmpv)
              PLANE_IO_DATA(20+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
           end do
        end do
     end do
     do k=pc2min_,pc2max_
        do j=pc1min_,pc1max_
           do n=pnmin_,pnmax_
              im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
              ! Pressure gradient z
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_z(im  ,j,:)*P(im  ,j,k+kk-st2:k+kk+st1))
              end do
              dSCdx(1) = sum(interp_w_zm(im  ,j,:)*tmpv)
              do kk=-st1,st2
                 tmpv(kk) = sum(grad_z(im+1,j,:)*P(im+1,j,k+kk-st2:k+kk+st1))
              end do
              dSCdx(2) = sum(interp_w_zm(im+1,j,:)*tmpv)
              PLANE_IO_DATA(21+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
           end do
        end do
     end do
     ivar = 30+nscalar
     do isc=1,nscalar
        do k=pc2min_,pc2max_
           do j=pc1min_,pc1max_
              do n=pnmin_,pnmax_
                 i  = index_x(n);  wx2  = interp_x(n);  wx1  = 1.0_WP-wx2
                 ! Scalar gradients
                 ! dx
                 dSCdx(1) = sum(grad_x(i  ,j,:)*SC(i  -st2:i  +st1,j,k,isc))
                 dSCdx(2) = sum(grad_x(i+1,j,:)*SC(i+1-st2:i+1+st1,j,k,isc))
                 PLANE_IO_DATA(ivar+1+3*(isc-1))%var(n,j,k) = wx1 *dSCdx(1) + wx2 *dSCdx(2)
              end do
           end do
        end do
        do k=pc2min_,pc2max_
           do j=pc1min_,pc1max_
              do n=pnmin_,pnmax_
                 im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                 ! dy
                 do jj=-st1,st2
                    tmpv(jj) = sum(grad_y(im  ,j+jj,:)*SC(im  ,j+jj-st2:j+jj+st1,k,isc))
                 end do
                 dSCdx(1) = sum(interp_v_ym(im  ,j,:)*tmpv)
                 do jj=-st1,st2
                    tmpv(jj) = sum(grad_y(im+1,j+jj,:)*SC(im+1,j+jj-st2:j+jj+st1,k,isc))
                 end do
                 dSCdx(2) = sum(interp_v_ym(im+1,j,:)*tmpv)
                 PLANE_IO_DATA(ivar+2+3*(isc-1))%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              end do
           end do
        end do
        do k=pc2min_,pc2max_
           do j=pc1min_,pc1max_
              do n=pnmin_,pnmax_
                 im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                 ! dz
                 do kk=-st1,st2
                    tmpv(kk) = sum(grad_z(im  ,j,:)*SC(im  ,j,k+kk-st2:k+kk+st1,isc))
                 end do
                 dSCdx(1) = sum(interp_w_zm(im  ,j,:)*tmpv)
                 do kk=-st1,st2
                    tmpv(kk) = sum(grad_z(im+1,j,:)*SC(im+1,j,k+kk-st2:k+kk+st1,isc))
                 end do
                 dSCdx(2) = sum(interp_w_zm(im+1,j,:)*tmpv)
                 PLANE_IO_DATA(ivar+3+3*(isc-1))%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              end do
           end do
        end do
     end do
     ivar = ivar+3*nscalar

     ! Finitechem outputs
     if (use_HR) then
        do isc=1,N_nons+1
           do k=pc2min_,pc2max_
              do j=pc1min_,pc1max_
                 do n=pnmin_,pnmax_
                    im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                    ! Scalar source terms
                    PLANE_IO_DATA(ivar+isc)%var(n,j,k) = rhs_fc_out(n,j,k,isc) !wx1m*rhs_fc(im,j,k,isc) + wx2m*rhs_fc(im+1,j,k,isc)
                 end do
              end do
           end do
        end do
        ivar = ivar+N_nons+1

        do isc=1,nscalar
           do k=pc2min_,pc2max_
              do j=pc1min_,pc1max_
                 do n=pnmin_,pnmax_
                    im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                    ! Species diffusivities
                    PLANE_IO_DATA(ivar+isc)%var(n,j,k) = wx1m*DIFF(im,j,k,isc) + wx2m*DIFF(im+1,j,k,isc)
                 end do
              end do
           end do
        end do
        ivar = ivar+nscalar

        do isc=1,nscalar
           do k=pc2min_,pc2max_
              do j=pc1min_,pc1max_
                 do n=pnmin_,pnmax_
                    i  = index_x(n);  wx2  = interp_x(n);  wx1  = 1.0_WP-wx2
                    im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                    ! Species diffusive fluxes
                    PLANE_IO_DATA(ivar+1+3*(isc-1))%var(n,j,k) = wx1*diff_flux(1,i,j,k,isc) + wx2*diff_flux(1,i+1,j,k,isc)
                    dSCdx(1) = sum(interp_v_ym(im  ,j,:)*diff_flux(2,im  ,j-st1:j+st2,k,isc))
                    dSCdx(2) = sum(interp_v_ym(im+1,j,:)*diff_flux(2,im+1,j-st1:j+st2,k,isc))
                    PLANE_IO_DATA(ivar+2+3*(isc-1))%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
                    dSCdx(1) = sum(interp_w_zm(im  ,j,:)*diff_flux(3,im  ,j,k-st1:k+st2,isc))
                    dSCdx(2) = sum(interp_w_zm(im+1,j,:)*diff_flux(3,im+1,j,k-st1:k+st2,isc))
                    PLANE_IO_DATA(ivar+3+3*(isc-1))%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
                 end do
              end do
           end do
        end do
        ivar = ivar+3*nscalar

        do isc=1,nscalar
           do k=pc2min_,pc2max_
              do j=pc1min_,pc1max_
                 do n=pnmin_,pnmax_
                    im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
                    ! Species diffusive flux divergence
                    PLANE_IO_DATA(ivar+isc)%var(n,j,k) = wx1m*diff_flux_div(im,j,k,isc) + wx2m*diff_flux_div(im+1,j,k,isc)
                 end do
              end do
           end do
        end do
        ivar = ivar+nscalar
     end if


     ! Stress tensor gradients
     do n=pnmin_,pnmax_
        i  = index_x(n);  wx2  = interp_x(n);  wx1  = 1.0_WP-wx2
        im = index_xm(n); wx2m = interp_xm(n); wx1m = 1.0_WP-wx2m
        do j=pc1min_,pc1max_
           do k=pc2min_,pc2max_
              ! dTAUdx11
              FX = 0.0_WP
              do ii=i-1,i+1
                 FX(ii-i) = &
                      + 2.0_WP*VISC(ii,j,k)*( &
                      + sum(grad_u_x(ii,j,:)*U(ii-stv1:ii+stv2,j,k)) &
                      - 1.0_WP/3.0_WP*( sum(divv_u(ii,j,:)*U(ii-stv1:ii+stv2,j,k)) &
                                      + sum(divv_v(ii,j,:)*V(ii,j-stv1:j+stv2,k))  &
                                      + sum(divv_w(ii,j,:)*W(ii,j,k-stv1:k+stv2))) )
              end do
              dSCdx(1) = sum(divv_xx(i  ,j,:)*FX(-stv2:stv1))
              dSCdx(2) = sum(divv_xx(i+1,j,:)*FX( stv1:stv2))
              PLANE_IO_DATA(22+nscalar)%var(n,j,k) = wx1*dSCdx(1) + wx2*dSCdx(2)
              ! dTAUdx12
              FY = 0.0_WP
              do ii=i,i+1
                 do jj=j-stv1,j+stv2
                    FY(jj-j) = &
                         + sum(interp_sc_xy(ii,jj,:,:)*VISC(ii-st2:ii+st1,jj-st2:jj+st1,k)) * &
                         ( sum(grad_u_y(ii,jj,:)*U(ii,jj-stv2:jj+stv1,k)) &
                         + sum(grad_v_x(ii,jj,:)*V(ii-stv2:ii+stv1,jj,k)) )
                 end do
                 dSCdx(ii-i+1) = sum(divv_xy(ii,j,:)*FY(-stv1:stv2))
              end do
              PLANE_IO_DATA(23+nscalar)%var(n,j,k) = wx1*dSCdx(1) + wx2*dSCdx(2)
              ! dTAUdx13
              FZ = 0.0_WP
              do ii=i,i+1
                 do kk=k-stv1,k+stv2
                    FZ(kk-k) = &
                         + sum(interp_sc_xz(ii,j,:,:)*VISC(ii-st2:ii+st1,j,kk-st2:kk+st1)) * &
                         ( sum(grad_u_z(ii,j,:)*U(ii,j,kk-stv2:kk+stv1)) &
                         + sum(grad_w_x(ii,j,:)*W(ii-stv2:ii+stv1,j,kk)) )
                 end do
                 dSCdx(ii-i+1) = sum(divv_xz(ii,j,:)*FZ(-stv1:stv2))
              end do
              PLANE_IO_DATA(24+nscalar)%var(n,j,k) = wx1*dSCdx(1) + wx2*dSCdx(2)

              ! dTAUdx21
              FX = 0.0_WP
              do jj=j,j+1
                 do ii=im-stv1,im+stv2+1
                    FX(ii-im-1) = &
                         + sum(interp_sc_xy(ii,jj,:,:)*VISC(ii-st2:ii+st1,jj-st2:jj+st1,k)) * &
                         ( sum(grad_u_y(ii,jj,:)*U(ii,jj-stv2:jj+stv1,k)) &
                         + sum(grad_v_x(ii,jj,:)*V(ii-stv2:ii+stv1,jj,k)) )
                 end do
                 dSCdx(jj-j+1) = wx1m*sum(divv_yx(im,jj,:)*FX(-stv2:stv1)) + wx2m*sum(divv_yx(im+1,jj,:)*FX(stv1:stv2))
              end do
              PLANE_IO_DATA(25+nscalar)%var(n,j,k) = sum(interp_v_ym(im,j,:)*dSCdx)
              ! dTAUdx22
              FY = 0.0_WP
              do ii=im,im+1
                 do jj=j-stv2,j+stv2
                    FY(jj-j) = &
                         + 2.0_WP*VISC(ii,jj,k)*( &
                         + sum(grad_v_y(ii,jj,:)*V(ii,jj-stv1:jj+stv2,k)) &
                         - 1.0_WP/3.0_WP*( sum(divv_u(ii,jj,:)*U(ii-stv1:ii+stv2,jj,k)) &
                                         + sum(divv_v(ii,jj,:)*V(ii,jj-stv1:jj+stv2,k)) &
                                         + sum(divv_w(ii,jj,:)*W(ii,jj,k-stv1:k+stv2))) )
                 end do
                 tmpv(1) = sum(divv_yy(ii,j  ,:)*FY(-stv2:stv1))
                 tmpv(2) = sum(divv_yy(ii,j+1,:)*FY( stv1:stv2))
                 dSCdx(ii-im+1) = sum(interp_v_ym(ii,j,:)*tmpv)
                 !dSCdx(ii-i+1) = sum(interp_v_ym(ii,j,:)*tmpv)
              end do
              PLANE_IO_DATA(26+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              !PLANE_IO_DATA(26+nscalar)%var(n,j,k) = wx1*dSCdx(1) + wx2*dSCdx(2)
              ! dTAUdx23
              FZ = 0.0_WP
              do ii=im,im+1
                 do jj=j,j+1
                    do kk=k-stv1,k+stv2
                       FZ(kk-k) = &
                            + sum(interp_sc_yz(ii,jj,:,:)*VISC(ii,jj-st2:jj+st1,kk-st2:kk+st1)) * &
                            ( sum(grad_v_z(ii,jj,:)*V(ii,jj,kk-stv2:kk+stv1)) &
                            + sum(grad_w_y(ii,jj,:)*W(ii,jj-stv2:jj+stv1,kk)) )
                    end do
                    tmpv(jj-j) = sum(divv_yz(ii,jj,:)*FZ(-stv1:stv2))
                 end do
                 dSCdx(ii-im+1) = sum(interp_v_ym(ii,j,:)*tmpv)
              end do
              PLANE_IO_DATA(27+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)

              !dTAUdx31
              FX = 0.0_WP
              do kk=k,k+1
                 do ii=im-stv1,im+stv2+1
                    FX(ii-im-1) = &
                         + sum(interp_sc_xz(ii,j,:,:)*VISC(ii-st2:ii+st1,j,kk-st2:kk+st1)) * &
                         ( sum(grad_u_z(ii,j,:)*U(ii,j,kk-stv2:kk+stv1)) &
                         + sum(grad_w_x(ii,j,:)*W(ii-stv2:ii+stv1,j,kk)) )
                 end do
                 dSCdx(kk-k+1) = wx1m*sum(divv_zx(im,j,:)*FX(-stv2:stv1)) + wx2m*sum(divv_zx(im+1,j,:)*FX(stv1:stv2))
              end do
              PLANE_IO_DATA(28+nscalar)%var(n,j,k) = sum(interp_w_zm(im,j,:)*dSCdx)
              ! dTAUdx32
              FY = 0.0_WP
              do ii=im,im+1
                 do kk=k,k+1
                    do jj=j-stv1,j+stv2
                       FY(jj-j) = &
                            + sum(interp_sc_yz(ii,jj,:,:)*VISC(ii,jj-st2:jj+st1,kk-st2:kk+st1)) * &
                            ( sum(grad_v_z(ii,jj,:)*V(ii,jj,kk-stv2:kk+stv1)) &
                            + sum(grad_w_y(ii,jj,:)*W(ii,jj-stv2:jj+stv1,kk)) )
                    end do
                    tmpv(kk-k) = sum(divv_zy(ii,j,:)*FY(-stv1:stv2))
                 end do
                 dSCdx(ii-im+1) = sum(interp_w_zm(ii,j,:)*tmpv)
              end do
              PLANE_IO_DATA(29+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
              ! dTAUdx33
              FZ = 0.0_WP
              do ii=im,im+1
                 do kk=k-stv2,k+stv2
                    FZ(kk-k) = &
                         + 2.0_WP*VISC(ii,j,kk)*( &
                         + sum(grad_w_z(ii,j,:)*W(ii,j,kk-stv1:kk+stv2)) &
                         - 1.0_WP/3.0_WP*( sum(divv_u(ii,j,:)*U(ii-stv1:ii+stv2,j,kk))  &
                                         + sum(divv_v(ii,j,:)*V(ii,j-stv1:j+stv2, kk))  &
                                         + sum(divv_w(ii,j,:)*W(ii,j,kk-stv1:kk+stv2)) ))
                 end do
                 tmpv(1) = sum(divv_zz(ii,j,:)*FZ(-stv2:stv1))
                 tmpv(2) = sum(divv_zz(ii,j,:)*FZ( stv1:stv2))
                 dSCdx(ii-im+1) = sum(interp_w_zm(ii,j,:)*tmpv)
              end do
              PLANE_IO_DATA(30+nscalar)%var(n,j,k) = wx1m*dSCdx(1) + wx2m*dSCdx(2)
           end do
        end do
     end do
  end select

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
end subroutine dump_plane_compute


! ========================================================== !
! Finalize the communicators                                 !
! ========================================================== !
subroutine dump_plane_finalize
  use dump_plane
  implicit none
  integer :: ierr

  if (nplanes_.gt.0) then
     call MPI_COMM_FREE(PLANE_IO_COMM,  ierr)
  else
     call MPI_COMM_FREE(PLANE_NON_COMM, ierr)
  end if

  return
end subroutine dump_plane_finalize
