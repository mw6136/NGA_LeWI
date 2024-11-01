! ============================================================ !
!                        planeStats.f90                        !
! Program to compute conditional statistics in physical and    !
!    spectral space using planar data from an inhomogeneous    ! 
!    configuration. Requires input data from the routine in    !
!    postprocess/dump_plane.f90.                               !
!                                                              !
! Author: Jonathan F. MacArt                                   !
! Date:   August 29, 2016                                      !
! ============================================================ !
module planeStats
  use parallel
  use precision
  use string
  implicit none

  ! Data file names
  character(len=str_long) :: data_dir
  character(len=str_long), dimension(:), pointer :: data_files
  integer :: nfiles
  character(len=str_medium) :: output_name, scaling_type, cond_var
  integer, dimension(:), pointer :: ntime_start

  ! Indices for the probes
  integer :: nPROG
  real(WP), dimension(:),   pointer :: plane_PROG
  integer,  dimension(:,:), pointer :: plane_jndx
  
  ! Parallel decomposition
  integer :: ipmin, ipmin_, ipmax, ipmax_ ! Probes
  integer :: pnmin, pnmin_, pnmax, pnmax_ ! Planes
  integer :: nprobes, nprobes_
  integer :: nplanes, nplanes_
  integer :: dplanes !planes orientation: 1=xy, 2=xz, 3=yz
  integer :: nd1, nd1_, nd2, nd2_, nx, ny, nz

  ! Conditional statistics
  integer :: nbins_cond
  integer :: icondmin, icondmax
  real(WP) :: condmin, condmax, flameloc_cond, prog_upper, prog_lower
  real(WP), dimension(:), pointer :: bins_cond
  integer,  dimension(:,:,:), pointer :: bin_index
  integer,  dimension(:,:,:), pointer :: ns_cond, ns_cond_total
  real(WP), dimension(:,:,:), pointer :: y_cond

  ! Output options
  logical :: output_scatter, output_pdfs

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
  integer, dimension(:), pointer :: nvars
  logical :: use_dSC, use_qSC
  integer :: icount_all, icount_dSC, icount_qSC
  
  ! Generalized indices for STAT_DATA decomposition
!!$  integer :: i1min, i1max, i2min, i2max, i3min, i3max
  integer :: isc_ZMIX, isc_PROG, isc_H, isc_O2, isc_O, isc_OH, isc_H2, isc_H2O, isc_HO2, isc_H2O2, isc_T

  ! Data fields to read
  real(WP) :: dt
  integer :: ntime_curr, ntime_dSC, ntime_qSC
  integer :: nrec, nscalar, stepsize
  integer :: nrec_
  integer :: N_tot
  real(WP), dimension(:), pointer :: W_sp
  real(WP), dimension(:), pointer :: d1, d2, dp, x, y
  character(len=str_short), dimension(:), pointer :: names, SC_name
  integer, dimension(:), pointer :: ntime
  real(WP), dimension(:), pointer :: time
  ! Mean in physical space: *m
  ! Mean in conditional space: *m_c
  real(WP), dimension(:,:,:),     pointer :: RHO
  real(WP), dimension(:,:,:,:),   pointer :: dRHOdx
  real(WP), dimension(:,:,:),     pointer :: VISC, VISC_fc
  real(WP), dimension(:,:,:),     pointer :: P, Pf
  ! Physical space
  real(WP), dimension(:,:),       pointer :: RHOm, Pm
  real(WP), dimension(:,:,:),     pointer :: dRHOdxm
  real(WP), dimension(:,:,:,:),   pointer :: U, Ui
  real(WP), dimension(:,:,:,:,:), pointer :: dUdx
  ! Conditional space
  real(WP), dimension(:,:),       pointer :: RHOm_c, Pm_c
  real(WP), dimension(:,:,:),     pointer :: dRHOdxm_c
  real(WP), dimension(:,:,:,:),   pointer :: U_c
  real(WP), dimension(:,:,:,:,:), pointer :: dUdx_c
  ! Stay instantaneous
  real(WP), dimension(:,:,:,:),   pointer :: dPdx
  real(WP), dimension(:,:,:,:,:), pointer :: dTAUdx
  real(WP), dimension(:,:,:,:,:), pointer :: dUdxi
  real(WP) :: U_bulk, U_coflow, H_jet, P_mean

  !    Um : U_tilde (Favre averaged)
  !    _c : in conditional space
  !    Instantaneous quantities become fluctuating later on
  real(WP), dimension(:,:),     pointer :: NU, NU_c
  real(WP), dimension(:,:,:),   pointer :: Um
  real(WP), dimension(:,:,:),   pointer :: Um_c
  real(WP), dimension(:,:,:,:), pointer :: dUdxm
  real(WP), dimension(:,:,:,:), pointer :: dUdxm_c

  real(WP), dimension(:,:,:,:),  pointer :: dTAUdxm, dTAUdxm_c
  real(WP), dimension(:,:,:,:,:),pointer :: dTAUdxp, dTAUdxp_c
  real(WP), dimension(:,:,:),    pointer :: dPdxm, dPdxm_c
  real(WP), dimension(:,:,:,:),  pointer :: dPdxp, dPdxp_c

  ! Scalars
  real(WP), dimension(:,:,:,:), pointer :: SC, SC_c, SCi
  real(WP), dimension(:,:,:),   pointer :: SCm, SCm_c
  ! Scalar gradients
  real(WP), dimension(:,:,:,:,:), pointer :: dSCdx, dSCdxi
  real(WP), dimension(:,:,:,:,:), pointer :: dSCdx_c
  real(WP), dimension(:,:,:,:),   pointer :: dSCdxm
  real(WP), dimension(:,:,:,:),   pointer :: dSCdxm_c
!!$  real(WP), dimension(:,:,:,:,:), pointer :: dSCdxp
!!$  real(WP), dimension(:,:,:,:,:), pointer :: dSCdxp_c
  ! Source terms
  real(WP), dimension(:,:,:,:), pointer :: src_SC, src_SC_fc
  ! Species diffusivities
  real(WP), dimension(:,:,:,:),   pointer :: DIFF, DIFF_fc
  ! Species diffusive fluxes and divergence
  real(WP), dimension(:,:,:,:,:), pointer :: diff_flux
  real(WP), dimension(:,:,:,:),   pointer :: diff_flux_div

contains
  ! ========================================================== !
  ! Compute the conditional mapping                            !
  ! ========================================================== !
  subroutine cond_map
    implicit none
    integer :: iplane, j, k
    real(WP) :: tmp1, tmp2, ZH, ZO
    
    if (trim(cond_var).eq.'ZMIX') then
       ! Based on mixing between N2-diluted center H2 jet and air
       !tmp2 = 0.077942_WP/1.008_WP + 2.0_WP*0.232_WP/16.0_WP
       tmp2 = 0.077942_WP/(2.0_WP*W_sp(isc_H)) + 2.0_WP*0.232_WP/(0.5_WP*W_sp(isc_O2))
       do iplane=pnmin_,pnmax_
          do j=1,ny
             do k=1,nz
!!$              ! Conserved scalar mixture fraction
!!$              bin_index(k,j,iplane) = min(nbins_cond, max(1,floor(SC(k,j,iplane,isc_ZMIX)*nbins_cond)+1))

                ! Bilger (1988) mixture fraction (Peters, Turbulent Combustion, p. 175)
                !tmp1 = (SC(k,j,iplane,isc_H2)/1.008_WP + 2.0_WP*(0.232_WP-SC(k,j,iplane,isc_O2))/16.0_WP)/tmp2
                ZH = + 1.0_WP*W_sp(isc_H)/W_sp(isc_H   )*SC(k,j,iplane,isc_H   ) &
                     + 2.0_WP*W_sp(isc_H)/W_sp(isc_H2  )*SC(k,j,iplane,isc_H2  ) &
                     + 1.0_WP*W_sp(isc_H)/W_sp(isc_OH  )*SC(k,j,iplane,isc_OH  ) &
                     + 2.0_WP*W_sp(isc_H)/W_sp(isc_H2O )*SC(k,j,iplane,isc_H2O ) &
                     + 1.0_WP*W_sp(isc_H)/W_sp(isc_HO2 )*SC(k,j,iplane,isc_HO2 ) &
                     + 2.0_WP*W_sp(isc_H)/W_sp(isc_H2O2)*SC(k,j,iplane,isc_H2O2)
                ZO = + 1.0_WP*W_sp(isc_O)/W_sp(isc_O   )*SC(k,j,iplane,isc_O   ) &
                     + 2.0_WP*W_sp(isc_O)/W_sp(isc_O2  )*SC(k,j,iplane,isc_O2  ) &
                     + 1.0_WP*W_sp(isc_O)/W_sp(isc_OH  )*SC(k,j,iplane,isc_OH  ) &
                     + 1.0_WP*W_sp(isc_O)/W_sp(isc_H2O )*SC(k,j,iplane,isc_H2O ) &
                     + 2.0_WP*W_sp(isc_O)/W_sp(isc_HO2 )*SC(k,j,iplane,isc_HO2 ) &
                     + 2.0_WP*W_sp(isc_O)/W_sp(isc_H2O2)*SC(k,j,iplane,isc_H2O2)
                tmp1 = ( ZH/(2.0_WP*W_sp(isc_H)) - 2.0_WP*(ZO-0.232_WP)/(0.5_WP*W_sp(isc_O2)) )/tmp2
                bin_index(k,j,iplane) = min(nbins_cond, max(1,floor(tmp1*nbins_cond)+1))
                y_cond(k,j,iplane) = tmp1
             end do
          end do
       end do
    elseif (trim(cond_var).eq.'PROG') then
       do iplane=pnmin_,pnmax_
          do j=1,ny
             do k=1,nz
                ! Define the progress variable based on O2
                tmp1 = (prog_upper - SC(k,j,iplane,isc_PROG))/(prog_upper-prog_lower)
                bin_index(k,j,iplane) = min(nbins_cond, max(1,floor(tmp1*nbins_cond)+1))
                y_cond(k,j,iplane) = tmp1
             end do
          end do
       end do
    end if
  end subroutine cond_map
  
  ! ========================================================== !
  ! Compute the mean of a vector in the z-direction            !
  ! ========================================================== !
  subroutine zmean(sol,mean)
    implicit none
    real(WP), dimension(nz,ny), intent(in) :: sol
    real(WP), dimension(ny), intent(out) :: mean
    integer :: j

    do j=1,ny
       mean(j) = sum(sol(:,j))/real(nz,WP)
    end do
  end subroutine zmean

  ! ========================================================== !
  ! Compute the mean of a vector over the conditional space    !
  ! ========================================================== !
  subroutine cmean(iplane,sol,mean,icount_cond)
    implicit none
    integer, intent(in) :: iplane, icount_cond
    real(WP), dimension(1:nz,1:ny),    intent(in)  :: sol
    real(WP), dimension(1:nbins_cond), intent(inout) :: mean
    integer,  dimension(1:nbins_cond) :: ns
    integer :: ii, j, k

    ! Initialize count
    ns = ns_cond_total(:,iplane,icount_cond)

    ! bin_index needs to be set first in planeStats_budgets_compute
    do j=1,ny
       do k=1,nz
          ii = bin_index(k,j,iplane)
          ns(ii) = ns(ii) + 1
          mean(ii) = ( real(ns(ii)-1,WP)*mean(ii) + sol(k,j) ) / real(ns(ii),WP)
       end do
    end do

    ! Save for cumulative count
    ns_cond(:,iplane,icount_cond) = ns - ns_cond_total(:,iplane,icount_cond)

  end subroutine cmean

  ! ========================================================== !
  ! Compute the mean of a vector in both spaces                !
  ! ========================================================== !
  subroutine condmean(iplane,sol,mean_z,mean_c,icount_cond)
    implicit none
    integer, intent(in) :: iplane, icount_cond
    real(WP), dimension(1:nz,1:ny),    intent(in)    :: sol
    real(WP), dimension(1:ny),         intent(inout) :: mean_z
    real(WP), dimension(1:nbins_cond), intent(inout) :: mean_c
    integer,  dimension(1:nbins_cond) :: ns
    integer :: ii, j, k, ntime_avg

    select case(icount_cond)
    case(1)
       ntime_avg = ntime_curr
    case(2)
       ntime_avg = ntime_dSC
    case(3)
       ntime_avg = ntime_qSC
    end select

    ! Mean in z
    do j=1,ny
       mean_z(j) = ( real(ntime_avg-1,WP)*mean_z(j) + sum(sol(:,j))/real(nz,WP) ) / real(ntime_avg,WP)
    end do

    ! Initialize count
    ns = ns_cond_total(:,iplane,icount_cond)

    ! bin_index needs to be set first in planeStats_budgets_compute
    do j=1,ny
       do k=1,nz
          ii = bin_index(k,j,iplane)
          ns(ii) = ns(ii) + 1
          mean_c(ii) = ( real(ns(ii)-1,WP)*mean_c(ii) + sol(k,j) ) / real(ns(ii),WP)
       end do
    end do

    ! Save for cumulative count
    ns_cond(:,iplane,icount_cond) = ns - ns_cond_total(:,iplane,icount_cond)

  end subroutine condmean


!!$  ! ========================================================== !
!!$  ! Update a histogram of sol with bounds sol_min,sol_max      !
!!$  ! ========================================================== !
!!$  subroutine hist(sol,prob,sol_min,sol_max)
!!$    implicit none
!!$    real(WP), dimension(1:nz,1:ny), intent(in) :: sol
!!$    real(WP), dimension(1:nbins_cond), intent(out) :: prob
!!$    real(WP), intent(in) :: sol_min,sol_max
!!$    real(WP) :: tmp1
!!$    integer  :: j,k,bin
!!$
!!$    ! Zero the probability
!!$    prob = 0.0_WP
!!$
!!$    ! Count the distribution of sol into prob
!!$    do j=1,ny
!!$       do k=1,nz
!!$          tmp1 = (sol(k,j)-sol_min)/(sol_max-sol_min)
!!$          bin = min(nbins_cond, max(1,floor(tmp1*nbins_cond)+1))
!!$          prob(bin) = prob(bin) + 1.0_WP
!!$       end do
!!$    end do    
!!$  end subroutine hist

end module planeStats


! ========================================================== !
! MAIN ROUTINE                                               !
! ========================================================== !
program planeStats_main
  use planeStats
  use parser
  use fileio
  use cli_reader
  implicit none

  character(len=str_medium) :: input_name, output_type
  integer :: ierr, iunit, ifile, iprobe, irec
  integer :: i,j, k, n
  logical :: isdef
  real(WP), dimension(:),   pointer :: plane_PROG1

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

  ! Output type
  call parser_read('Output type', output_type)

  ! Size of each step for analysis
  call parser_read('Stepsize', stepsize, 1)

  ! Get the number of scalars as well as output parameters
  call parser_read('Number of scalars', nscalar)
  call parser_read('Output file', output_name)
  call GETNSPECIES(N_tot)
  allocate(W_sp (1:N_tot))
  W_sp = 0.0_WP
  call GETMOLARMASS(W_sp)
  
  ! Binary data files
  call parser_read('Data root dir', data_dir)
  call parser_getsize('Data files', nfiles)
  allocate(data_files(nfiles))
  call parser_read('Data files', data_files)
  do ifile=1,nfiles
     data_files(ifile) = trim(data_dir)//'/'//trim(data_files(ifile))
  end do

  ! Get scaling info
  call parser_read('Centerline velocity', U_bulk)
!!$  call parser_read('Coflow velocity', U_coflow)
  call parser_read('Jet height', H_jet)
  call parser_read('Cond flame loc', flameloc_cond)

  ! Get the number of records and basic parameters
  allocate(ntime(nfiles))
  allocate(ntime_start(nfiles))
  call planeStats_mpi_read_nrec
  
  ! Check that nplanes matches nproc - this isn't actually doing anything
  if (nproc.gt.nplanes) then
     print *, 'Error: nproc must not be greater than nplanes=', nplanes
  end if

  ! Synchronize some parameters between the processes
  nrec = 0
  do i=1,nfiles
     nrec = nrec + (ntime(i) - ntime_start(i) + 1)/stepsize + 1
  end do

  ! For debugging
  call parser_is_defined('Max records',isdef)
  if (isdef) then
     call parser_read('Max records',nrec)
     ntime = nrec*stepsize
     nfiles = 1
  end if
  if (irank.eq.iroot) print *, 'Number of records :', nrec
  allocate(time(1:nrec))

  ! Decompose the planes evenly among the processors
  ! Global indices
  pnmin = 1
  pnmax = nplanes
  ! Local indices -- decompose planes only
  pnmin_ = (irank-1)*pnmax/nproc + 1
  pnmax_ = irank*pnmax/nproc
  ! Local data sizes
  nplanes_ = pnmax_-pnmin_+1
  nrec_    = nrec
  print '(a6,i3,a8,i3,a8,i3,a10,i3)', 'irank=', irank, ' pnmin_=', pnmin_, ' pnmax_=', pnmax_, ' nplanes_=', nplanes_

  ! Initialize the metrics
!  call planeStats_metric_init

  ! Initialize the MPI arrays
  call planeStats_mpi_init

  ! Set up the conditional averaging
  call parser_read('Cond var', cond_var,'none')
  call parser_read('Min cond var', condmin,0.0_WP)
  call parser_read('Max cond var', condmax,1.0_WP)
  if (trim(cond_var).eq.'PROG') then
     call parser_read('Prog upper', prog_upper)
     call parser_read('Prog lower', prog_lower)
  else
     ! For mixture fraction
     prog_upper = 1.0_WP
     prog_lower = 0.0_WP
  end if
  call parser_read('Bin count', nbins_cond,100)
  ! Set up the bins
  allocate(bins_cond(1:nbins_cond+1))
  do i=1,nbins_cond+1
     bins_cond(i) = (real(i-1,WP))*(condmax-condmin)/real(nbins_cond,WP) + condmin
  end do
  allocate(bin_index(1:nz,1:ny,pnmin_:pnmax_))

  allocate(ns_cond      (1:nbins_cond,pnmin_:pnmax_,1:3))
  allocate(ns_cond_total(1:nbins_cond,pnmin_:pnmax_,1:3))
  allocate(y_cond(1:nd2,1:nd1,pnmin_:pnmax_))
  y_cond = 0.0_WP

  ! Indices to the conditional averaging counters
  icount_all = 1
  icount_dSC = 2
  icount_qSC = 3

  ! Decide what to do
  select case(trim(output_type))
  case('budgets')
     ! Initialize the averaging and budget modules
     call planeStats_average_init
     call planeStats_budgets_init

     ! Read and compute average quantities for budgets
     ns_cond_total = 0
     call planeStats_mpi_read(2)
     call planeStats_average_finalize
     ! Output alignment statistics
     call planeStats_average_output

     ! Read and compute budgets
     ns_cond_total = 0
     call planeStats_mpi_read(3)
     call planeStats_budgets_output

  case('scatter')
     ! Output data for scatter plots
     if (irank.eq.iroot) print *, 'ny,nz=', ny, nz
     call parser_read('Write scatter',output_scatter,.true.)
     call parser_read('Write pdfs',output_pdfs,.false.)
     call planeStats_scatter_init
     call planeStats_mpi_read(4)

  case('cond budgets')
     ! Compute conditional (CMC) budgets
     ns_cond_total = 0
     call planeStats_cond_init
     call planeStats_mpi_read(5)
     call planeStats_cond_finalize
     call planeStats_cond_output

  case('verification')
     ns_cond_total = 0
     call planeStats_average_init
     call planeStats_mpi_read(2)
     call planeStats_average_finalize
     call planeStats_average_output

  case default
     print *, 'output type not implemented'
  end select


  ! Finalize MPI
  call MPI_FINALIZE(ierr)

end program planeStats_main

