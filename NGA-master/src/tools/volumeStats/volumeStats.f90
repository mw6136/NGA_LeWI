! ============================================================ !
!                       volumeStats.f90                        !
! Program to compute statistics in physical and conditional    !
!    space using volume data from an inhomogeneous             ! 
!    configuration. Requires input data from the routine in    !
!    postprocess/dump_volume.f90.                              !
!                                                              !
! Author: Jonathan F. MacArt                                   !
! Date:   October 16, 2017                                     !
! ============================================================ !
module volumeStats
  use parallel
  use precision
  use string
  implicit none

  ! Data file names
  character(len=str_long) :: data_dir, data_base
  character(len=str_long), dimension(:), pointer :: data_files
  integer :: nfiles, nfile_step
  integer, dimension(:), pointer :: ifile_start, ifile_end
  character(len=str_medium) :: output_name, scaling_type, cond_var

  ! Indices for planar decomposition
  logical :: use_planes
  logical :: need_gradients
  real(WP) :: H_cj
  integer :: pnmin, pnmin_, pnmax, pnmax_, nplanes_
  real(WP), dimension(:), pointer :: pnloc
  integer,  dimension(:), pointer :: pnindx
  real(WP), dimension(:,:), pointer :: interp_pnxm
  
  ! Parallel decomposition
  integer :: imin, imax, imino, imaxo
  integer :: jmin, jmax, jmino, jmaxo
  integer :: kmin, kmax, kmino, kmaxo
  integer :: imin_, imax_, imino_, imaxo_
  integer :: jmin_, jmax_, jmino_, jmaxo_
  integer :: kmin_, kmax_, kmino_, kmaxo_
  integer :: jmid
  integer :: nx, ny, nz, nvars
  integer :: nx_, ny_, nz_
  integer :: nxo_, nyo_, nzo_
  integer :: nover, noverx, novery
  real(WP) :: nzi
  real(WP), dimension(:), pointer :: x, y, z
  real(WP), dimension(:), pointer :: xm, ym, zm

  ! Conditional statistics
  integer, parameter :: stz = 4
  integer :: nbins_cond
  integer :: icondmin, icondmax
  real(WP) :: condmin, condmax, prog_upper, prog_lower
  real(WP), dimension(:),     pointer :: dci
  real(WP), dimension(:,:),   pointer :: interp_c
  real(WP), dimension(:),     pointer :: bins_cond
  integer,  dimension(:,:),   pointer :: nsyz_cond, nsyz_cond_total
  integer,  dimension(:,:,:), pointer :: nsz_cond, nsz_cond_total
  integer,  dimension(:,:,:), pointer :: nsz_cond_temp, nsz_cond_temp_total
  integer,  dimension(:,:,:), pointer :: bin_index
  real(WP), dimension(:,:,:), pointer :: y_cond

  ! Output options
  logical :: output_scatter, output_pdfs

  ! LES filter
  character(len=str_medium) :: filter_type
  real(WP) :: filter_width, tflaminar

  ! ================================== !
  ! Type definitions for MPI IO        !
  ! for a given variable :             !
  !   - name : name of the variable    !
  !   - var  : pointer to data array   !
  !   - view : MPI IO view on the file !
  ! ================================== !
  type MPI_IO_VAR
     character(len=str_short) :: name
     real(WP), dimension(:,:,:), pointer :: var
     integer :: view
  end type MPI_IO_VAR
  
  ! Arrays with all variables
  integer :: STAT_DATA_NVARS
  type(MPI_IO_VAR), dimension(:), pointer :: STAT_DATA
  
  ! Scalar information
  integer :: nscalar, N_tot
  real(WP), dimension(:), pointer :: W_sp
  character(len=str_short), dimension(:), pointer :: SC_name
  integer :: isc_ZMIX, isc_PROG, isc_H, isc_O2, isc_O, isc_OH, isc_H2, isc_H2O, isc_HO2, isc_H2O2, isc_T

  ! Time information
  real(WP) :: dt
  integer :: ntime_curr
  integer :: nrec, nrec_
  integer, dimension(:), pointer :: ntime
  real(WP), dimension(:), pointer :: time

  ! Data fields to read
  logical :: combust, use_lfs_stripes
  real(WP) :: rho_in, visc_in, temp_in, diff_in
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:),   pointer :: RHO
  real(WP), dimension(:,:,:),   pointer :: P
  real(WP), dimension(:,:,:,:), pointer :: U, Ui
  real(WP), dimension(:,:,:,:), pointer :: SC

  ! Transport properties
  real(WP), dimension(:,:,:),   pointer :: VISC
  real(WP), dimension(:,:,:,:), pointer :: DIFF

  ! Species source terms
  real(WP), dimension(:,:,:,:), pointer :: src_SC

  ! Species diffusive fluxes and divergence
  real(WP), dimension(:,:,:,:,:), pointer :: diff_flux
  real(WP), dimension(:,:,:,:),   pointer :: diff_flux_div

  ! Physical space gradients
  real(WP), dimension(:,:,:,:),   pointer :: dRHOdx
  real(WP), dimension(:,:,:,:),   pointer :: dPdx
  real(WP), dimension(:,:,:,:,:), pointer :: dUdx
  real(WP), dimension(:,:,:,:,:), pointer :: dSCdx
  real(WP), dimension(:,:,:,:,:), pointer :: dTAUdx

  ! Conditional gradients
  real(WP), dimension(:,:,:,:),   pointer :: dRHOdx_c
  real(WP), dimension(:,:,:,:),   pointer :: dPdx_c
  real(WP), dimension(:,:,:,:,:), pointer :: dUdx_c
  real(WP), dimension(:,:,:,:,:), pointer :: dSCdx_c
  real(WP), dimension(:,:,:,:,:), pointer :: dTAUdx_c

contains

  ! ========================================================== !
  ! Compute the conditional mapping                            !
  ! ========================================================== !
  subroutine cond_map(SC_in,cond_out)
    implicit none

    integer :: i, j, k
    real(WP) :: tmp1, tmp2, ZH, ZO
    real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:nscalar), intent(in) :: SC_in
    real(WP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_), intent(out) :: cond_out

    if (.not.combust) return
    
    if (trim(cond_var).eq.'ZMIX') then
       ! Based on mixing between N2-diluted center H2 jet and air
       tmp2 = 0.077942_WP/(2.0_WP*W_sp(isc_H)) + 2.0_WP*0.232_WP/(0.5_WP*W_sp(isc_O2))
       !$OMP PARALLEL DO PRIVATE(j,i,ZH,ZO,tmp1)
       do k=kmin_,kmax_
          do j=jmin_,jmax_
             do i=imin_,imax_
                ! Bilger (1988) mixture fraction (Peters, Turbulent Combustion, p. 175)
                ZH = + 1.0_WP*W_sp(isc_H)/W_sp(isc_H   )*SC_in(i,j,k,isc_H   ) &
                     + 2.0_WP*W_sp(isc_H)/W_sp(isc_H2  )*SC_in(i,j,k,isc_H2  ) &
                     + 1.0_WP*W_sp(isc_H)/W_sp(isc_OH  )*SC_in(i,j,k,isc_OH  ) &
                     + 2.0_WP*W_sp(isc_H)/W_sp(isc_H2O )*SC_in(i,j,k,isc_H2O ) &
                     + 1.0_WP*W_sp(isc_H)/W_sp(isc_HO2 )*SC_in(i,j,k,isc_HO2 ) &
                     + 2.0_WP*W_sp(isc_H)/W_sp(isc_H2O2)*SC_in(i,j,k,isc_H2O2)
                ZO = + 1.0_WP*W_sp(isc_O)/W_sp(isc_O   )*SC_in(i,j,k,isc_O   ) &
                     + 2.0_WP*W_sp(isc_O)/W_sp(isc_O2  )*SC_in(i,j,k,isc_O2  ) &
                     + 1.0_WP*W_sp(isc_O)/W_sp(isc_OH  )*SC_in(i,j,k,isc_OH  ) &
                     + 1.0_WP*W_sp(isc_O)/W_sp(isc_H2O )*SC_in(i,j,k,isc_H2O ) &
                     + 2.0_WP*W_sp(isc_O)/W_sp(isc_HO2 )*SC_in(i,j,k,isc_HO2 ) &
                     + 2.0_WP*W_sp(isc_O)/W_sp(isc_H2O2)*SC_in(i,j,k,isc_H2O2)
                tmp1 = ( ZH/(2.0_WP*W_sp(isc_H)) - 2.0_WP*(ZO-0.232_WP)/(0.5_WP*W_sp(isc_O2)) )/tmp2
                bin_index(i,j,k) = min(nbins_cond, max(1,floor(tmp1*nbins_cond)+1))
                cond_out(i,j,k) = tmp1
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    elseif (trim(cond_var).eq.'PROG') then
       !$OMP PARALLEL DO PRIVATE(j,i,tmp1)
       do k=kmin_,kmax_
          do j=jmin_,jmax_
             do i=imin_,imax_
                ! Define the progress variable based on O2
                tmp1 = (prog_upper - SC_in(i,j,k,isc_PROG))/(prog_upper-prog_lower)
                bin_index(i,j,k) = min(nbins_cond, max(1,floor(tmp1*nbins_cond)+1))
                cond_out(i,j,k) = tmp1
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine cond_map

  ! ========================================================== !
  ! Get the filename header                                    !
  ! ========================================================== !
  subroutine get_name(iplane,of_name)
    implicit none
    integer, intent(in) :: iplane
    character(len=str_short) :: of_name

    if (iplane<10) then
       write(of_name,'(A3,I1)') 'xH0', nint(pnloc(iplane)/H_cj)
    else
       write(of_name,'(A2,I2)') 'xH', nint(pnloc(iplane)/H_cj)
    end if
    
!!$    if (iplane<10) then
!!$       write(of_name,'(A1,I1)') '0', iplane
!!$    else
!!$       write(of_name,'(I2)') iplane
!!$    end if

    return
  end subroutine get_name

end module volumeStats


! ========================================================== !
! MAIN ROUTINE                                               !
! ========================================================== !
program volumeStats_main
  use volumeStats
  use parser
  use fileio
  use cli_reader
  implicit none

  character(len=str_medium) :: input_name, output_type
  integer :: ierr, iunit, ifile, irec, iext
  integer :: i, j, k, n, iplane, nseq, ifile_offs
  integer, dimension(:), pointer :: list
  logical :: isdef
  integer :: q, r
  integer :: tmpi1, tmpi2
  real(WP) :: tmp1, tmp2, tmp3

  ! Initialize the parallel environment
  call parallel_init

  ! Parse the input file
  if (irank.eq.iroot) then
     call get_command_argument(1,input_name)
  end if
  call MPI_BCAST(input_name, str_medium, MPI_CHARACTER, iroot-1, MPI_COMM_WORLD, ierr)
  call parser_init
  call parser_parsefile(input_name)

  ! Initialize filesystem
  call parallel_init_filesystem
  call volumeStats_mpi_topology_init

  ! Output type
  call parser_read('Output type', output_type)

  ! User provides a header for the output files
  !   N.B. -- Generates many files, so put these in a folder
  call parser_read('Output file', output_name)

  ! Option to set LFS stripe size in volumeStats_mpi_read
  call parser_read('Use LFS stripes',use_lfs_stripes,.false.)

  ! Binary data files
  ! User must provide the following in the input file:
  !    Data root dir :    Relative or absolute path to directory containing the volume data files
  !    Data file base :   Filename before the numeric sequence
  !    Data file range :  The start and end values of the data file range
  !                           Can be specified as a single range, e.g.,
  !                              1101,1105
  !                           or multiple ranges, e.g.,
  !                              1101,1105
  !                              2501,2505
  call parser_read('Data root dir', data_dir)
  call parser_read('Data file base',data_base)
  call parser_getsize('Data file range',nseq)
  if (mod(nseq,2).ne.0) stop 'Must specify start and end for each data file sequence'
  allocate(list(nseq))
  call parser_read('Data file range',list)
  allocate(ifile_start(nseq/2))
  allocate(ifile_end(nseq/2))
  nfiles = 0
  do n=1,nseq/2
     ifile_start(n) = list(1+(n-1)*2)
     ifile_end(n)   = list(2+(n-1)*2)
     nfiles = nfiles + ifile_end(n) - ifile_start(n) + 1
  end do
  if (irank.eq.iroot) print *, ifile_start
  if (irank.eq.iroot) print *, ifile_end
  if (irank.eq.iroot) print *, nfiles

  ! Create the data file names
  allocate(data_files(nfiles))
  ifile_offs = 0
  do n=1,nseq/2
     do ifile=1,ifile_end(n)-ifile_start(n)+1
        write(data_files(ifile+ifile_offs),'(I8.8)') ifile_start(n)+ifile-1
        data_files(ifile+ifile_offs) = trim(data_dir)//'/'//trim(data_base)//'_'//trim(data_files(ifile+ifile_offs))
        if (irank.eq.iroot) print *, trim(data_files(ifile+ifile_offs))
     end do
     ifile_offs = ifile_offs + ifile_end(n)-ifile_start(n)+1
  end do

  ! Allocate time indices
  allocate(ntime(nfiles))
  allocate(time(nfiles))

  ! Number of overlap cells - need for mpi_read_nrec
  nover = 2
  noverx = nover
  novery = nover
  call volumeStats_mpi_read_nrec

  ! For debugging 
  !   -- Use this to limit the number of records to read before generating output files
  call parser_is_defined('Max records',isdef)
  if (isdef) then
     call parser_read('Max records',nrec)
     ntime = nrec
     nfiles = nrec
  else
     nrec = nfiles
  end if

  ! Enable to skip records
  call parser_read('Stepsize',nfile_step,1)
  nrec = nrec/nfile_step
  if (irank.eq.iroot) print *, 'Number of records :', nrec

  ! Global cells include the entire domain
  imin = 1;  jmin = 1;  kmin = 1;
  imax = nx; jmax = ny; kmax = nz;
  ! Overlap cells with room for gradients
  imino = imin-noverx; imaxo = imax+noverx
  jmino = jmin-novery; jmaxo = jmax+novery
  kmino = kmin-nover;  kmaxo = kmax+nover


  ! Initialize the grid
  call volumeStats_metric_grid


  ! Set the local indices
  ! Two options:
  !    1. If planes are defined, restrict the grid to just the necessary neighboring cells
  !    2. If planes are not defined, partition the grid evenly among the processors
  call parser_is_defined('Plane x locations',isdef)
  if (isdef) then
     use_planes = .true.
     ! Get central jet height
     call parser_read('Central jet height',H_cj)
     pnmin = 1
     ! Get plane locations
     call parser_getsize('Plane x locations',pnmax)
     allocate(pnloc(pnmax))
     call parser_read('Plane x locations',pnloc)
     allocate(pnindx(pnmax))
     allocate(interp_pnxm(pnmax,2))
     pnindx = 0
     interp_pnxm = 0.0_WP

!!$     select case(trim(output_type))
!!$     case('les')
!!$        ! Find the planes on each processor
!!$        pnmax_ = 0
!!$        nplanes_ = 0
!!$        do iplane=pnmin,pnmax
!!$           if (pnloc(iplane).gt.xm(nx-1)) stop 'plane x-index not in domain'
!!$           if (pnloc(iplane).ge.x(imin_) .and. pnloc(iplane).lt.x(imax_)) then
!!$              if (nplanes_.eq.0) pnmin_ = iplane
!!$              pnmax_ = iplane
!!$              nplanes_ = nplanes_ + 1
!!$              
!!$              ! Compute the interpolation to the planes
!!$              i = imin_
!!$              do while (xm(i).lt.pnloc(iplane))
!!$                 i = i+1
!!$              end do
!!$              interp_pnxm(iplane,1) = (xm(i)-pnloc(iplane))/(xm(i)-xm(i-1))
!!$              interp_pnxm(iplane,2) = 1.0_WP-interp_pnxm(iplane,1)
!!$              pnindx(iplane) = i
!!$              pnloc (iplane) = interp_pnxm(iplane,1)*xm(i-1) + interp_pnxm(iplane,2)*xm(i)
!!$           end if
!!$        end do

!!$     case('budgets')
     ! RANS-type budgets
     ! Adjust the decomposition based on the amount of data we need to read
     !  -- Matching the number of processors to the number of planes
     !        will minimize the amount of data to read in each step
     ! Need to divide planes evenly on the processes
     if (pnmax.ne.nproc) stop 'Must decompose planes evenly on processors'

     ! Just one plane per processor
     pnmin_ = irank
     pnmax_ = irank
     nplanes_ = 1
     iplane = irank
     if (pnloc(iplane).gt.xm(nx-1)) stop 'plane x-index not in domain'

     ! Compute the interpolation to the planes
     i = imin
     do while (xm(i).lt.pnloc(iplane))
        i = i+1
     end do
     interp_pnxm(iplane,1) = (xm(i)-pnloc(iplane))/(xm(i)-xm(i-1))
     interp_pnxm(iplane,2) = 1.0_WP-interp_pnxm(iplane,1)
     pnindx(iplane) = i
     pnloc (iplane) = interp_pnxm(iplane,1)*xm(i-1) + interp_pnxm(iplane,2)*xm(i)

     ! Set the local indices from the plane location
     iext   = 4
     imin_  = pnindx(pnmin_)-iext
     imax_  = pnindx(pnmax_)+iext
     jmin_  = jmin
     jmax_  = ny

     print '(i5,a10,10i15)',    irank, ' pnindx= ', pnindx(pnmin_:pnmax_)
     print '(i5,a10,10ES15.5)', irank, ' pnloc=  ', pnloc(pnmin_:pnmax_)

  else
     ! Planes not defined; decompose the grid evenly
     use_planes = .false.

     ! Local indices - x
     q = nx/npx
     r = mod(nx,npx)
     if (iproc<=r) then
        nx_   = q+1
        imin_ = imin + (iproc-1)*(q+1)
     else
        nx_   = q
        imin_ = imin + r*(q+1) + (iproc-r-1)*q
     end if
     imax_  = imin_ + nx_ - 1

     ! Local indices - y
     !   Only decompose in y for LES
     q = ny/npy
     r = mod(ny,npy)
     if (jproc<=r) then
        ny_   = q+1
        jmin_ = jmin + (jproc-1)*(q+1)
     else
        ny_   = q
        jmin_ = jmin + r*(q+1) + (jproc-r-1)*q
     end if
     jmax_  = jmin_ + ny_ - 1

  end if

  ! Local indices - k
  kmin_ = kmin
  kmax_ = kmax


  ! For LES filtering: adjust the number of x- and y-overlap cells to accomodate
  !    the number of cells required for the filter width

  ! NEED TO COMMUNICATE IF CELLS OVERLAP -- GENERALIZE
  ! e.g., planes need only some cells from neighboring planes. give start/end of range of comm overlap cells.
  ! update mpi_volume_update_border

  select case(trim(output_type))
  case('les')
     ! Get the constant filter width
     call parser_read('Filter type', filter_type)
     call parser_read('Filter width', filter_width)
     print*, 'Filter width =',filter_width

     select case(trim(filter_type))
     case('homogeneous')
        ! Minimum x-grid spacing at the overlap
        if (iproc.eq.1) then
           tmp1 = x(imax_)-x(imax_-1)
        else if (iproc.eq.npx) then
           tmp1 = x(imin_+1)-x(imin_)
        else
           tmp1 = min(x(imin_+1)-x(imin_), x(imax_)-x(imax_-1))
        end if

        ! Minimum y-grid spacing at the overlap
        if (jproc.eq.1) then
           tmp2 = y(jmax_)-y(jmax_-1)
        else if (jproc.eq.npy) then
           tmp2 = y(jmin_+1)-y(jmin_)
        else
           tmp2 = min(y(jmin_+1)-y(jmin_), y(jmax_)-y(jmax_-1))
        end if

        ! Adjust the z-grid spacing and the limits for MPI file read
        tmp3 = z(kmin_+1) - z(kmin_)

        ! Compute the required number of overlap cells
        noverx = ceiling(filter_width/tmp1)
        novery = ceiling(filter_width/tmp2)
        
        if (noverx<nover) noverx = nover
        if (novery<nover) novery = nover

        nover  = ceiling(filter_width/tmp3)

     case('inhomogeneous')        
        ! changed for inhomogeneous filter
        noverx = ceiling(filter_width)+2
        novery = ceiling(filter_width)+2
        nover  = ceiling(filter_width)+2 
        
     case default
        stop 'filter type not implemented'
     end select
     
     ! Print from each process
     if (irank.eq.iroot) print '(a10,3a12,3a10)','irank','min_dxo','min_dyo','min_dz','noverx','novery','nover'
     print '(I10,3ES12.3,3I10)', irank, tmp1, tmp2, tmp3, noverx, novery, nover

     if (.not.use_planes) then
        ! Get the max noverx and novery from all processes
        tmpi1 = noverx
        tmpi2 = novery
        call MPI_ALLREDUCE(tmpi1, noverx, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
        call MPI_ALLREDUCE(tmpi2, novery, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
     end if

     ! Check if overlap cells are going to collide between processes
     !   -- A problem with planes, since MPI will read from imino_ to imaxo_
!!$     if (use_planes .and. noverx.gt.nx_) stop 'volumeStats: noverx is too large (reduce filter width)'

  case('budgets')
     if (npy.gt.1) stop 'volumeStats: cannot decompose domain in y for RANS budgets'

  case default
     ! Not LES filtering -- do not change noverx and novery
  end select


  ! Local data sizes
  nx_ = imax_-imin_+1
  ny_ = jmax_-jmin_+1
  nz_ = kmax_-kmin_+1
  nrec_    = nrec

  if (irank.eq.iroot) print '(10a10)', 'irank','imin_','imax_','jmin_','jmax_','kmin_','kmax_','nx_','ny_','nz_'
  print '(10i10)', irank,imin_,imax_,jmin_,jmax_,kmin_,kmax_,nx_,ny_,nz_
  if (irank.eq.iroot) print *, 'nx=',nx_,'ny=',ny_,'nz=', nz_
  call MPI_BARRIER(comm, ierr)

  
  ! Finalize the overlap cells
  imino_ = imin_-noverx
  imaxo_ = imax_+noverx
  jmino_ = jmin_-novery
  jmaxo_ = jmax_+novery
  kmino_ = kmin_-nover
  kmaxo_ = kmax_+nover
  if(imino_ .LT. 0)    imino_ = imin
  if(imaxo_ .GT. nx+1) imaxo_ = nx+1
  if(jmino_ .LT. 0)    jmino_ = jmin
  if(jmaxo_ .GT. ny+1) jmaxo_ = ny+1
  select case(trim(output_type))
  case('les')
     ! Don't clip kmino_ and kmaxo_
  case default
     if(kmino_ .LT. 0)    kmino_ = kmin
     if(kmaxo_ .GT. nz+1) kmaxo_ = nz+1
  end select
  nxo_ = imaxo_-imino_+1
  nyo_ = jmaxo_-jmino_+1
  nzo_ = kmaxo_-kmino_+1
  imino = imino_
  imaxo = imaxo_
  jmino = jmino_
  jmaxo = jmaxo_
  kmino = kmino_
  kmaxo = kmaxo_


!!$  select case(trim(output_type))
!!$  case('les')
!!$     call volumeStats_metric_grid_reset
!!$  end select

  ! Initialize the metrics -- must be done after decomposition
  call volumeStats_metric_init

  ! Read the chemistry type
  call parser_read('Combustion',combust,.true.)
  if (combust) then
     ! Get the number of scalars
     call GETNSPECIES(N_tot)
     allocate(W_sp (1:N_tot))
     W_sp = 0.0_WP
     call GETMOLARMASS(W_sp)
     nscalar = N_tot + 2
     call volumeStats_finitechem_init
  else
     N_tot = 0; nscalar = 0;
     W_sp = 0.0_WP
     call parser_read('Density',rho_in)
     call parser_read('Viscosity',visc_in)
     call parser_read('Temperature',temp_in)
     call parser_read('Diffusivity',diff_in)
  end if

  ! Allocate and decompose the arrays to read
  call volumeStats_mpi_init

  ! Set up the conditional averaging
  call parser_read('Cond var', cond_var,'none')
  call parser_read('Min cond var', condmin,0.0_WP)
  call parser_read('Max cond var', condmax,1.0_WP)

  ! Always read prog_upper and lower, since we're defining this for conditional averages
  call parser_read('Prog upper', prog_upper)
  call parser_read('Prog lower', prog_lower)
  call parser_read('Bin count', nbins_cond,100)

  ! Allocate the conditional metrics
  allocate(dci(1:nbins_cond+1))
  allocate(interp_c(1:nbins_cond+1,-1:0))

  ! All conditional bins
  ! dci is inverse of conditional bin spacing!
  dci      = real(nbins_cond,WP)/(condmax-condmin)
  interp_c = 0.5_WP
  ! First bin
  interp_c(1,-1) = 0.0_WP
  interp_c(1, 0) = 1.0_WP
  dci(1) = 2.0_WP*dci(1)
  ! Last bin
  interp_c(nbins_cond+1,-1) = 1.0_WP
  interp_c(nbins_cond+1, 0) = 0.0_WP
  dci(nbins_cond+1) = 2.0_WP*dci(nbins_cond+1)

  ! Set up the bins
  allocate(bins_cond(1:nbins_cond+1))
  do i=1,nbins_cond+1
     bins_cond(i) = (real(i-1,WP))*(condmax-condmin)/real(nbins_cond,WP) + condmin
  end do
  do i=1,nbins_cond
     bins_cond(i) = 0.5_WP*(bins_cond(i)+bins_cond(i+1))
  end do
  allocate(bin_index(imin_:imax_,jmin_:jmax_,kmin_:kmax_))

  allocate(nsyz_cond      (imin_:imax_,1:nbins_cond))
  allocate(nsyz_cond_total(imin_:imax_,1:nbins_cond))
  ! These go from jmid to jmax_ -- can average the symmetric -y and +y regions
  allocate(nsz_cond       (imin_:imax_,jmid:jmax_,1:nbins_cond))
  allocate(nsz_cond_total (imin_:imax_,jmid:jmax_,1:nbins_cond))
  allocate(nsz_cond_temp       (imin_:imax_,jmid:jmax_,1:nbins_cond))
  allocate(nsz_cond_temp_total (imin_:imax_,jmid:jmax_,1:nbins_cond))
  ! Conditional variable
  allocate(y_cond(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  y_cond = 0.0_WP

  ! Only compute gradients if necessary
  need_gradients = .false.

  ! Decide what to do
  select case(trim(output_type))
  case('budgets')
     nsyz_cond_total = 0
     nsz_cond_total  = 0
     nsz_cond_temp_total  = 0
     need_gradients  = .true.
     call volumeStats_average_init
     call volumeStats_mpi_read(1)
     call volumeStats_average_finalize
     call volumeStats_budgets_init
     call volumeStats_budgets_compute

  case('les')
     need_gradients  = .true.
     call volumeStats_les_init
     call volumeStats_les_filter_init
     call volumeStats_mpi_read(2)
     if (use_planes) then
        call volumeStats_les_finalize
     else
        call volumeStats_les_ensight
     end if

  case('verification')
     nsyz_cond_total = 0
     nsz_cond_total  = 0
     nsz_cond_temp_total  = 0
     need_gradients  = .true.
     call volumeStats_average_init
     call volumeStats_mpi_read(1)
     call volumeStats_average_finalize
     call volumeStats_budgets_init
     call volumeStats_budgets_compute

  case default
     print *, 'output type not implemented'
  end select


  ! Finalize MPI
  call MPI_FINALIZE(ierr)

end program volumeStats_main

