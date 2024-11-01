! ============================================================ !
!                       dnsbox_corr.f90                        !
! Program to compute two-point, two-time correlations from DNS !
!    of isotropic turbulence. Requires NGA data files in the   !
!    prescribed format.                                        !
!                                                              !
! Author: Jonathan F. MacArt                                   !
! Date:   September 30, 2017                                   !
! ============================================================ !
module dnsbox_corr
  use precision
  use string
  use, intrinsic :: ISO_C_binding
  implicit none

  ! OpenMP
  integer :: nthreads, void
  integer, external :: OMP_GET_MAX_THREADS
  
  ! Sizes
  integer :: nx,ny,nz
  integer*8 :: N_total
  integer :: xper,yper,zper
  integer :: icyl
  integer :: nvar
  real(WP) :: dt,dt_file,time,time_prev
  real(WP) :: dt1, dt2, time1, time2
  real(WP) :: L
  character(len=str_short), dimension(:), pointer :: names
  character(len=str_medium) :: config
  integer, dimension(:,:), pointer :: mask
  
  ! Data
  real(WP), dimension(:,:,:),   pointer :: data
  real(WP), dimension(:,:,:,:), pointer :: Um
  real(WP), dimension(:), pointer :: mom_U1
  real(WP), dimension(:), pointer :: mom_U2
  real(WP), dimension(:), pointer :: mom_U2c
  
  ! Mesh
  real(WP), dimension(:), pointer :: x,xm
  real(WP), dimension(:), pointer :: y,ym
  real(WP), dimension(:), pointer :: z,zm

  ! Switches for different computation methods
  integer :: temp_corr_step
  character(len=str_short) :: method
  logical :: use_FFT, use_BF

  ! Fourier transform of velocity
  integer, parameter :: kind_cpx=16
  complex*16, dimension(:,:,:),   pointer :: in, out
  complex*16, dimension(:,:,:),   pointer :: u_hat, u_hat2
  complex*16, dimension(:,:,:,:), pointer :: u_hat1, UU_hat1
  complex*16, dimension(:,:,:),   pointer :: UU_hat

  ! Two-point longitudinal velocity correlation
  real(WP), dimension(:), pointer :: r, temp_corr_time
  real(WP), dimension(:,:,:), pointer :: Bij, Fii, Kii
  real(WP), dimension(:,:), pointer :: Biij

  ! Energy spectrum
  real(WP) :: u_rms, epsilon, visc, RHO, TKE, lambda
  real(WP), dimension(:,:), pointer :: spectrum, spectrum2

end module dnsbox_corr


! ========================================================== !
! MAIN ROUTINE                                               !
! ========================================================== !
program dnsbox_corr_main
  use dnsbox_corr
  use parser
  use fileio
  use cli_reader
  implicit none
  include "fftw3.f03"
  
  character(len=str_medium) :: input_name
  character(len=str_medium) :: data_dir, out_dir
  character(len=str_long) :: fconfig
  character(len=str_long) :: fdata, fdata_p
  character(len=str_medium), dimension(:), pointer :: ftime
  integer :: Nfiles,Nfiles_max,ifile
  integer :: ifile1, ifile2, Nfiles_max1, Nfiles_max2, Nfiles_max_temp, start_file1
  integer :: iunit,iunit1,iunit2,var
  integer :: i,j,k,ierr,ii,jj,i2,iimax
  real(WP) :: tmp, tmp1, tmp2, dx, nxi
  real(WP) :: L_11, Re_L11, Re_lambda
  logical :: save_FFT

  ! Initialize FFTW OpenMP threads
  void = fftw_init_threads()
  !$OMP PARALLEL
  nthreads = OMP_GET_MAX_THREADS()
  !$OMP END PARALLEL
  call fftw_plan_with_nthreads(int(nthreads,c_int))

!!$  ! Initialize the parallel module
!!$  call parallel_init
!!$  call parallel_init_filesystem

  ! Parse the input file
!!$  if (irank.eq.iroot) then
     call get_command_argument(1,input_name)
!!$  end if
!!$  call MPI_BCAST(input_name, str_medium, MPI_CHARACTER, iroot-1, MPI_COMM_WORLD, ierr)
  call parser_init
  call parser_parsefile(input_name)

  ! Read the source config file
  call parser_read('Root dir',data_dir,'.')
  call parser_read('Configuration file',fconfig)
  fconfig = trim(data_dir)//'/'//trim(fconfig)
  call BINARY_FILE_OPEN(iunit,trim(fconfig),"r",ierr)
  call BINARY_FILE_READ(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_READ(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_READ(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_READ(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(mask(nx,ny))
  call BINARY_FILE_READ(iunit,x,nx+1,kind(x),ierr)
  call BINARY_FILE_READ(iunit,y,ny+1,kind(y),ierr)
  call BINARY_FILE_READ(iunit,z,nz+1,kind(z),ierr)
  call BINARY_FILE_READ(iunit,mask,nx*ny,kind(mask),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! Finish creating the meshes
  allocate(xm(nx),ym(ny),zm(nz))
  do i=1,nx
     xm(i) = 0.5_WP*(x(i)+x(i+1))
  end do
  do j=1,ny
     ym(j) = 0.5_WP*(y(j)+y(j+1))
  end do
  do k=1,nz
     zm(k) = 0.5_WP*(z(k)+z(k+1))
  end do

  ! Get domain size
  L = x(nx+1) - x(1)
  dx = L/real(nx,WP)
  print*,'Domain size:',L
  
  ! Save total number of grid points in a 64-bit integer
  N_total = nx*ny*nz

  ! Get the list of data files to read
  call parser_read('Data file prefix',fdata_p)
  call parser_getsize('Data files at time',Nfiles)
  print *, 'Nfiles=',Nfiles
  allocate(ftime(Nfiles))
  call parser_read('Data files at time',ftime)
  call parser_read('Max Nfiles',Nfiles_max,Nfiles)
  if (Nfiles_max.gt.Nfiles) Nfiles_max = Nfiles
  call parser_read('Start file index',start_file1,1)

  ! Properties
  call parser_read('Viscosity',visc)
  call parser_read('Density',RHO)
  visc = visc/RHO
  
  ! Choose computation method
  call parser_read('Method',method)
  if ((trim(method).eq.'FFT').or.(trim(method).eq.'both')) then
     use_FFT = .true.
     allocate(in (nx,nx,nx))
     allocate(out(nx,nx,nx))
     call parser_read('Save Fourier coefficients',save_FFT,.true.)
     print *, 'Save Fourier coefficients? ', save_FFT
  end if
  if ((trim(method).eq.'BF').or.(trim(method).eq.'both')) then
     use_BF = .true.
  end if

  ! Choose the output type
  call parser_read('Temporal correlation step',temp_corr_step)

  select case (temp_corr_step)
  case (1)
     ! FIRST STEP IN TEMPORAL CORRELATIONS ############################################################
     !   -- Compute FFT of velocity components
     !   -- Save single-snapshot Fourier coefficients for selected snapshots
     !   -- Output single-time, two- and three-point correlations

     ! Just look at three components of velocity
     allocate(data(nx,ny,nz))
     if (use_BF) allocate(Um(nx,ny,nz,3))

     ! Allocate moments
     allocate(mom_U1 (3))
     allocate(mom_U2 (3)) !6))
     allocate(mom_U2c(3)) !6))

     ! Allocate spectra
     if (use_FFT) then
        ! u_hat is complex (large memory footprint!)
        allocate(u_hat (nx,nx,nx))
        allocate(UU_hat(nx,nx,nx))
        allocate(spectrum(2,nx+1))
        allocate(spectrum2(2,nx+1))
     end if

     ! Allocate two-point correlations
     allocate(r(nx))
     allocate(Bij(nx,3,3))
     allocate(Biij(nx,3))
     do i=1,nx
        r(i) = dx*(i-1)
     end do

     ! Write headers for output
     print '(a7,12a15)', 'step', 'time', 'TKE (spect)', 'TKE (corr)', 'eps (int)', 'eps (avg)', 'u_rms', 'L_t', 'Re_t', 'L_11', 'Re_L11', 'lambda', 'Re_lambda'

!!$     open(unit=iunit, file='velocity_'//trim(fdata_p), form='formatted', action='write')
!!$     write (iunit, '(a7,4a15)'), 'step', 'time', '<u1u1>',  '<u2u2>', '<u3u3>'
!!$     !write (iunit, '(a7,7a15)'), 'step', 'time', '<u1u1>', '<u1u2>', '<u1u3>', '<u2u2>', '<u2u3>', '<u3u3>'
!!$     close(iunit)

     ! Read the files and compute statistics
     time_prev = 0.0_WP
     Nfiles_max1 = min(start_file1+Nfiles_max-1, Nfiles)
     do ifile=start_file1, Nfiles_max1
        fdata = trim(data_dir)//'/'//trim(fdata_p)//'_'//trim(ftime(ifile))
        !print *, trim(ftime(ifile))
        !print *, trim(fdata)
        call BINARY_FILE_OPEN(iunit,trim(fdata),"r",ierr)
        call BINARY_FILE_READ(iunit,nx,  1,kind(nx),  ierr)
        call BINARY_FILE_READ(iunit,ny,  1,kind(ny),  ierr)
        call BINARY_FILE_READ(iunit,nz,  1,kind(nz),  ierr)
        call BINARY_FILE_READ(iunit,nvar,1,kind(nvar),ierr)
        call BINARY_FILE_READ(iunit,dt,  1,kind(dt),  ierr)
        call BINARY_FILE_READ(iunit,time,1,kind(time),ierr)
        allocate(names(nvar))
        do var=1,nvar
           call BINARY_FILE_READ(iunit,names(var),str_short,kind(names),ierr)
        end do

        dt_file = time-time_prev

        ! Zero the spectrum for this timestep
        spectrum = 0.0_WP
        spectrum2 = 0.0_WP
        epsilon  = 0.0_WP

        ! Start the spectrum files to save
        if (save_FFT) then
           ! u_hat
           fdata = 'u_hat_'//trim(fdata_p)//'_'//trim(ftime(ifile))
           call BINARY_FILE_OPEN (iunit1,trim(fdata),"w",ierr)
           call BINARY_FILE_WRITE(iunit1,nx,  1,kind(nx),  ierr)
           call BINARY_FILE_WRITE(iunit1,3,   1,kind(nvar),ierr)
           call BINARY_FILE_WRITE(iunit1,dt,  1,kind(dt),  ierr)
           call BINARY_FILE_WRITE(iunit1,time,1,kind(time),ierr)

           ! uu_hat
           fdata = 'uu_hat_'//trim(fdata_p)//'_'//trim(ftime(ifile))
           !print *, fdata
           call BINARY_FILE_OPEN (iunit2,trim(fdata),"w",ierr)
           call BINARY_FILE_WRITE(iunit2,nx,  1,kind(nx),  ierr)
           call BINARY_FILE_WRITE(iunit2,3,   1,kind(nvar),ierr)
           call BINARY_FILE_WRITE(iunit2,dt,  1,kind(dt),  ierr)
           call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
        end if

        ! Read in the velocity components -------------------------------------------------------------
        do ii=1,3
           call BINARY_FILE_READ_LONG(iunit,data,N_total,kind(data),ierr)

           ! Get the spectral velocity
           if (use_FFT) then
              call get_spectrum(1,ii,data,u_hat,spectrum)

              ! Get the three-point, one-time longitudinal velocity correlations (Biii) using FFT
              !$OMP PARALLEL DO PRIVATE(i,j)
              do k=1,nz
                 do j=1,ny
                    do i=1,nx
                       data(i,j,k) = data(i,j,k) * data(i,j,k)
                    end do
                 end do
              end do
              !$OMP END PARALLEL DO
              call get_spectrum(2,ii,data,UU_hat,spectrum2)

              ! Get the two-point, one-time longitudinal velocity correlations (Bii) using FFT
              call get_long_correlation(ii,u_hat,u_hat,Bij(:,ii,ii))
              !print *, 'Done correlation Bii', ii

              call get_long_correlation(ii,UU_hat,u_hat,Biij(:,ii))
              !print *, 'Done correlation Biii', ii

              ! Save the Fourier components u_hat and uu_hat
              if (save_FFT) then
                 call BINARY_FILE_WRITE(iunit1,u_hat ,N_total,kind_cpx,ierr)
                 call BINARY_FILE_WRITE(iunit2,UU_hat,N_total,kind_cpx,ierr)
              end if
           end if

           !print *, trim(names(ii)),' done'
        end do
        call BINARY_FILE_CLOSE(iunit,ierr)


        ! FFT outputs ---------------------------------------------------------------------
        if (use_FFT) then
           if (save_FFT) then
              ! Close the binary files
              call BINARY_FILE_CLOSE(iunit1,ierr)
              call BINARY_FILE_CLOSE(iunit2,ierr)
           end if

           ! Output the spectrum
           open(unit=iunit, file='spectrum_'//trim(fdata_p)//'_'//trim(ftime(ifile)), action='write')
           do i=1,nx+1
              if ((spectrum(1,i).ne.0.0_WP).and.(spectrum(2,i).ne.0.0_WP)) then
                 write(iunit,'(2ES22.13)') spectrum(1,i), spectrum(2,i)
              end if
           end do
           close(iunit)

           ! Output the two-point correlations
           open(unit=iunit, file='Bij_FFT_'//trim(fdata_p)//'_'//trim(ftime(ifile)), action='write')
           do i=1,nx
              write(iunit,'(4ES22.13)') r(i), Bij(i,1,1), Bij(i,2,2), Bij(i,3,3)
           end do
           close(iunit)

           ! Output the three-point correlations
           open(unit=iunit, file='Biij_FFT_'//trim(fdata_p)//'_'//trim(ftime(ifile)), action='write')
           do i=1,nx
              write(iunit,'(4ES22.13)') r(i), Biij(i,1), Biij(i,2), Biij(i,3)
           end do
           close(iunit)
        end if


        ! Brute-force outputs -----------------------------------------------------------------
        if (use_BF) then
           ! Save the velocity field for brute-force outputs
           !$OMP PARALLEL DO PRIVATE(j,i)
           do k=1,nz
              do j=1,ny
                 do i=1,nx
                    Um(i,j,k,ii) = data(i,j,k)
                 end do
              end do
           end do
           !$OMP END PARALLEL DO

!!$        ! Interpolate velocities to cell centers -- only needed if multiplying different velocity components
!!$        ! Assumes uniform grid spacing and periodicity in x,y,z
!!$        select case(ii)
!!$        case(1)
!!$           ! Interpolate U
!!$           !$OMP PARALLEL DO
!!$           do k=1,nz
!!$              do j=1,ny
!!$                 do i=1,nx-1
!!$                    Um(i,j,k,1) = 0.5_WP*(data(i,j,k) + data(i+1,j,k))
!!$                 end do
!!$                 Um(nx,j,k,1) = 0.5_WP*(data(nx,j,k) + data(1,j,k))
!!$              end do
!!$           end do
!!$           !$OMP END PARALLEL DO
!!$        case(2)
!!$           !$OMP PARALLEL DO
!!$           do k=1,nz
!!$              do j=1,ny-1
!!$                 do i=1,nx
!!$                    Um(i,j,k,2) = 0.5_WP*(data(i,j,k) + data(i,j+1,k))
!!$                 end do
!!$              end do
!!$              do i=1,nx
!!$                 Um(i,ny,k,2) = 0.5_WP*(data(i,ny,k) + data(i,1,k))
!!$              end do
!!$           end do
!!$           !$OMP END PARALLEL DO
!!$        case(3)
!!$           !$OMP PARALLEL DO
!!$           do k=1,nz-1
!!$              do j=1,ny
!!$                 do i=1,nx
!!$                    Um(i,j,k,3) = 0.5_WP*(data(i,j,k) + data(i,j,k+1))
!!$                 end do
!!$              end do
!!$           end do
!!$           !$OMP END PARALLEL DO
!!$           do j=1,ny
!!$              do i=1,nx
!!$                 Um(i,j,nz,3) = 0.5_WP*(data(i,j,nz) + data(i,j,1))
!!$              end do
!!$           end do
!!$        end select

           ! Get the two-point, one-time longitudinal velocity correlations (Bii) using brute force
           ! B11
           nxi = 1.0_WP/real(N_total,WP)
           Bij = 0.0_WP
           do ii=1,3
              do i2=1,nx
                 tmp1 = 0.0_WP
                 tmp2 = 0.0_WP
                 data = cshift(Um(:,:,:,ii), shift=(i2-1), dim=ii) !! would probably be much more efficient just to operate on small regions (data parallelism)
                 !$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(+:tmp1,tmp2)
                 do k=1,nz
                    do j=1,ny
                       do i=1,nx
                          ! Compute the correlation
                          tmp1 = tmp1 + Um(i,j,k,ii)*data(i,j,k)
                          tmp2 = tmp2 + Um(i,j,k,ii)*Um(i,j,k,ii)*data(i,j,k)
                       end do
                    end do
                 end do
                 !$OMP END PARALLEL DO
                 Bij(i2,ii,ii) = tmp1*nxi
                 Biij(i2,ii) = tmp2*nxi
              end do
              print *, 'Done BF', ii
           end do

           ! Output the two-point correlations
           open(unit=iunit, file='Bij_BF_'//trim(fdata_p)//'_'//trim(ftime(ifile)), action='write')
           do i=1,nx
              write(iunit,'(4ES22.13)') r(i), Bij(i,1,1), Bij(i,2,2), Bij(i,3,3)
           end do
           close(iunit)

           ! Output the three-point correlations
           open(unit=iunit, file='Biij_BF_'//trim(fdata_p)//'_'//trim(ftime(ifile)), action='write')
           do i=1,nx
              write(iunit,'(4ES22.13)') r(i), Biij(i,1), Biij(i,2), Biij(i,3)
           end do
           close(iunit)

        end if

!!$        ! Compute the velocity moments
!!$        ! First moment (mean)
!!$        do jj=1,3
!!$           tmp = 0.0_WP
!!$           !$OMP PARALLEL DO REDUCTION(+:tmp)
!!$           do k=1,nz
!!$              do j=1,ny
!!$                 do i=1,nx
!!$                    tmp = tmp + Um(i,j,k,jj)
!!$                 end do
!!$              end do
!!$           end do
!!$           !$OMP END PARALLEL DO
!!$           mom_U1(jj) = tmp
!!$        end do
!!$        mom_U1 = mom_U1/real(N_total,WP)
!!$
!!$        ! Second moment (variance)
!!$        do i2=1,3
!!$           tmp = 0.0_WP
!!$           !$OMP PARALLEL DO REDUCTION(+:tmp)
!!$           do k=1,nz
!!$              do j=1,ny
!!$                 do i=1,nx
!!$                    tmp = tmp + Um(i,j,k,i2)*Um(i,j,k,i2)
!!$                 end do
!!$              end do
!!$           end do
!!$           !$OMP END PARALLEL DO
!!$           mom_U2(i2) = tmp
!!$        end do
!!$        mom_U2 = mom_U2/real(N_total,WP)
!!$
!!$        ! Second centralized moment
!!$        do i2=1,3
!!$           mom_U2c(i2) = mom_U2(i2) - mom_U1(i2)*mom_U1(i2)
!!$        end do

        ! Get and print some mean quantities
!!$     print '(4a20)', ' ', '1', '2', '3'
!!$     print '(a20,3ES20.7)',  '1st moment:         ', mom_U1
!!$     print '(7a20)', ' ', '1,1', '1,2', '1,3', '2,2', '2,3', '3,3'
!!$     print '(a20,6ES20.7)', '2nd moment:         ', mom_U2
!!$     print '(a20,6ES20.7)', '2nd central moment: ', mom_U2c
!!$     print '(a20,ES20.7)' , 'TKE (moments):      ' ,0.5_WP*(mom_U2c(1)+mom_U2c(4)+mom_U2c(6))

        TKE     = sum(spectrum(2,:))
        epsilon = visc*epsilon/real(N_total,WP)
        u_rms   = sqrt(2.0_WP*TKE/3.0_WP)
!!$        L_11    = sum(Bij(1:nx/2,1,1))*dx/mom_U2c(1)
        L_11    = sum(Bij(1:nx/2,1,1))*dx/Bij(1,1,1)
        Re_L11  = u_rms*L_11/visc
!!$        lambda  = sqrt(2.0_WP*mom_U2c(1)/(lambda/real(nx**3,WP)))
        lambda  = sqrt(2.0_WP*Bij(1,1,1)/(lambda/real(N_total,WP)))
        Re_lambda = u_rms*lambda/visc

!!$        print '(i7,12ES15.5)', ifile, time, TKE, 0.5_WP*(mom_U2c(1)+mom_U2c(4)+mom_U2c(6)), visc*2.0_WP*sum(spectrum(1,:)**2*spectrum(2,:)), &
        print '(i7,12ES15.5)', ifile, time, TKE, 0.5_WP*(Bij(1,1,1)+Bij(1,2,2)+Bij(1,3,3)), visc*2.0_WP*sum(spectrum(1,:)**2*spectrum(2,:)), &
             epsilon, u_rms, u_rms**3/epsilon, u_rms**4/visc/epsilon, L_11, Re_L11, lambda, Re_lambda

!!$        open(unit=iunit, file='velocity_'//trim(fdata_p), form='formatted', status='old', position='append', action='write')
!!$        write(iunit, '(i7,4ES15.5)') ifile, time, mom_U2c
!!$        close(iunit)

        ! Save data for next time step
        time_prev = time
     end do

  case (2)
     ! SECOND STEP IN TEMPORAL CORRELATIONS ############################################################
     !   -- Read in single-snapshot Fourier coefficients
     !   -- Compute two-time correlations using IFFT
     !   -- Output two- and three-point, two-time correlations

     ! Allocate data arrays
     ! Fourier coefficients for first and second snapshots
     iimax = 1  ! JUST LOOK AT 1-COMPONENT to save memory
     allocate(u_hat1(nx,nx,nx,iimax))
     allocate(u_hat2(nx,nx,nx))

     ! Workspace array for three-point correlations
     allocate(UU_hat1(nx,nx,nx,iimax))

     call parser_read('Temporal corr max Nfiles2',Nfiles_max_temp,Nfiles)

     ! Folder for output files
     call parser_read('Output dir',out_dir)

     ! Allocate two-point correlations
     allocate(temp_corr_time(Nfiles_max_temp))
     allocate(Fii(Nfiles_max_temp,nx,iimax)); Fii = 0.0_WP
     allocate(Kii(Nfiles_max_temp,nx,iimax)); Kii = 0.0_WP
     allocate(r(nx))
     do i=1,nx
        r(i) = dx*(i-1)
     end do

     ! Print header for terminal output
     print '(a10,a15,a10,a15)', 'ifile1', 'time1', 'ifile2', 'time2'

     ! Loop over temporal combinations
     time_prev = 0.0_WP
     Nfiles_max1 = min(start_file1+Nfiles_max-1, Nfiles)
     do ifile1=start_file1, Nfiles_max1
        ! Read the base data files
        ! u_hat
        fdata = 'u_hat_'//trim(fdata_p)//'_'//trim(ftime(ifile1))
        call BINARY_FILE_OPEN(iunit,trim(fdata),"r",ierr)
        call BINARY_FILE_READ(iunit,nx,   1,kind(nx),  ierr)
        call BINARY_FILE_READ(iunit,nvar, 1,kind(nvar),ierr)
        call BINARY_FILE_READ(iunit,dt1,  1,kind(dt1),  ierr)
        call BINARY_FILE_READ(iunit,time1,1,kind(time1),ierr)
        do var=1,iimax
           call BINARY_FILE_READ_LONG(iunit,u_hat1(:,:,:,var),N_total,kind_cpx,ierr)
        end do
        call BINARY_FILE_CLOSE(iunit,ierr)

        ! uu_hat
        fdata = 'uu_hat_'//trim(fdata_p)//'_'//trim(ftime(ifile1))
        call BINARY_FILE_OPEN(iunit,trim(fdata),"r",ierr)
        call BINARY_FILE_READ(iunit,nx,   1,kind(nx),  ierr)
        call BINARY_FILE_READ(iunit,nvar, 1,kind(nvar),ierr)
        call BINARY_FILE_READ(iunit,dt1,  1,kind(dt1),  ierr)
        call BINARY_FILE_READ(iunit,time1,1,kind(time1),ierr)
        do var=1,iimax
           call BINARY_FILE_READ_LONG(iunit,uu_hat1(:,:,:,var),N_total,kind_cpx,ierr) ! does this memory usage need to be reduced?
        end do
        call BINARY_FILE_CLOSE(iunit,ierr)


        ! Write headers for the output files
        open(unit=iunit, file=trim(out_dir)//'/Fii_'//trim(fdata_p)//'_'//trim(ftime(ifile1)), form='formatted', action='write')
        write(iunit, '(a22)', advance='no'), 'time'
        do ii=1,iimax
           do i=1,nx
              if (i.eq.nx .and. ii.eq.iimax) then
                 write(iunit, '(ES22.13)') r(i)
              else
                 write(iunit, '(ES22.13)', advance='no') r(i)
              end if
           end do
        end do
        close(iunit)

        open(unit=iunit, file=trim(out_dir)//'/Kii_'//trim(fdata_p)//'_'//trim(ftime(ifile1)), form='formatted', action='write')
        write(iunit, '(a22)', advance='no'), 'time'
        do ii=1,iimax
           do i=1,nx
              if (i.eq.nx .and. ii.eq.iimax) then
                 write(iunit, '(ES22.13)') r(i)
              else
                 write(iunit, '(ES22.13)', advance='no') r(i)
              end if
           end do
        end do
        close(iunit)

        ! Loop over second data files
        Nfiles_max2 = min(ifile1+Nfiles_max_temp-1, Nfiles)
        do ifile2=ifile1, Nfiles_max2
           ! Read the second data file if not equal to the first
           if (ifile2.gt.ifile1) then
              fdata = 'u_hat_'//trim(fdata_p)//'_'//trim(ftime(ifile2))
              call BINARY_FILE_OPEN(iunit,trim(fdata),"r",ierr)
              call BINARY_FILE_READ(iunit,nx,   1,kind(nx),  ierr)
              call BINARY_FILE_READ(iunit,nvar, 1,kind(nvar),ierr)
              call BINARY_FILE_READ(iunit,dt2,  1,kind(dt2),  ierr)
              call BINARY_FILE_READ(iunit,time2,1,kind(time2),ierr)
           else
              time2 = time1
           end if

           ! Save time information
           temp_corr_time(ifile2-ifile1+1) = time2 - time1
           
           do ii=1,iimax
              if (ifile2.gt.ifile1) then
                 call BINARY_FILE_READ_LONG(iunit,u_hat2,N_total,kind_cpx,ierr)
              else
                 u_hat2 = u_hat1(:,:,:,ii)
              end if

              ! Get the two-point, two-time longitudinal velocity correlations (Fii) using FFT -- IS THIS THE CORRECT DEF OF FII?
              call get_long_correlation(ii,u_hat1(:,:,:,ii),u_hat2,Fii(ifile2-ifile1+1,:,ii))
              !print *, 'Done correlation Bii ', ii, ifile1, ifile2
              
              ! Get the three-point, two-time longitudinal velocity correlations (Kii) using FFT
              call get_long_correlation(ii,uu_hat1(:,:,:,ii),u_hat2,Kii(ifile2-ifile1+1,:,ii))
              !print *, 'Done correlation Biii', ii, ifile1, ifile2
           end do

           if (ifile2.gt.ifile1) call BINARY_FILE_CLOSE(iunit,ierr) 

           ! Write to terminal
           print '(i10,ES15.5,i10,ES15.5)', ifile1, time1, ifile2, time2

           ! Output two-time correlations
           open(unit=iunit, file=trim(out_dir)//'/Fii_'//trim(fdata_p)//'_'//trim(ftime(ifile1)), form='formatted', status='old', position='append', action='write')
           write(iunit, '(ES22.13)', advance='no') temp_corr_time(ifile2-ifile1+1)
           do ii=1,iimax
              do i=1,nx
                 if (i.eq.nx .and. ii.eq.iimax) then
                    write(iunit, '(ES22.13)') Fii(ifile2-ifile1+1,i,ii)
                 else
                    write(iunit, '(ES22.13)', advance='no') Fii(ifile2-ifile1+1,i,ii)
                 end if
              end do
           end do
           close(iunit)

           open(unit=iunit, file=trim(out_dir)//'/Kii_'//trim(fdata_p)//'_'//trim(ftime(ifile1)), form='formatted', status='old', position='append', action='write')
           write(iunit, '(ES22.13)', advance='no') temp_corr_time(ifile2-ifile1+1)
           do ii=1,iimax
              do i=1,nx
                 if (i.eq.nx .and. ii.eq.iimax) then
                    write(iunit, '(ES22.13)') Kii(ifile2-ifile1+1,i,ii)
                 else
                    write(iunit, '(ES22.13)', advance='no') Kii(ifile2-ifile1+1,i,ii)
                 end if
              end do
           end do
           close(iunit)
        end do ! ifile2
     end do ! ifile1

  end select



!!$  ! Initialize the MPI module and read data in parallel
!!$  call dnsbox_corr_MPI_init
  ! Use parallel_sum
  ! Could use parallel_gather_dir to communicate along lines  

end program dnsbox_corr_main


! ============================================================ !
! Compute the velocity spectra                                 !
! ============================================================ !
subroutine get_spectrum(opt,ivar,A,vel,S)
  use dnsbox_corr
  use precision
  implicit none
  include 'fftw3.f03'  
  integer, intent(in) :: opt, ivar
  real(WP), dimension(nx,nx,nx), intent(in) :: A
  complex*16, dimension(nx,nx,nx), intent(out) :: vel
  real(WP), dimension(2,nx+1), intent(out) :: S
  real(WP), dimension(:,:,:), pointer :: B
  real(WP) :: pi,dk,kc,eps,kx,ky,kz,kk
  integer(KIND=8) :: plan
  integer :: i,j,k,ik
  complex*16, parameter :: ii =(0.0_WP,1.0_WP)
  
  ! Create the plan
  call dfftw_plan_dft_3d(plan,nx,nx,nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
  
  ! FFT
  in = A
  call dfftw_execute(plan)
  out = out/real(N_total,WP)
!  print *, ivar, tmp_int, minval(real(in)), maxval(real(in)), minval(real(out)), maxval(real(out))
  
  ! Save velocity
  !$OMP PARALLEL DO PRIVATE(j,i)
  do k=1,nx
     do j=1,nx
        do i=1,nx
           vel(i,j,k) = out(i,j,k)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
  if (opt.eq.1) then
     ! Compute the spectrum
     allocate(B(nx,nx,nx))
     !$OMP PARALLEL DO PRIVATE(j,i)
     do k=1,nx
        do j=1,nx
           do i=1,nx
              B(i,j,k) = sqrt(real(out(i,j,k)*conjg(out(i,j,k))))
           end do
        end do
     end do
     !$OMP END PARALLEL DO

     pi=acos(-1.0_WP)
     dk=2.0_WP*pi/L
     kc=pi*nx/L
     eps=kc/1000000.0_WP
     !$OMP PARALLEL DO PRIVATE(j,i,kx,ky,kz,kk,ik) REDUCTION(+:S)
     do k=1,nx
        do j=1,nx
           do i=1,nx
              ! Wavenumbers
              kx=real(i-1,WP)*dk
              if (i.gt.(nx/2+1)) kx=-real(nx-i+1,WP)*dk
              ky=real(j-1,WP)*dk
              if (j.gt.(nx/2+1)) ky=-real(nx-j+1,WP)*dk
              kz=real(k-1,WP)*dk
              if (k.gt.(nx/2+1)) kz=-real(nx-k+1,WP)*dk
              kk=sqrt(kx**2+ky**2+kz**2)
              ! Spectrum
              ik=1+idint(kk/dk+0.5_WP)
              if ((kk.gt.eps).and.(kk.le.kc)) then
                 S(2,ik)=S(2,ik)+0.5_WP*(B(i,j,k)**2)
              end if
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     do ik=1,nx+1
        S(1,ik)=dk*(ik-1)
     end do

     deallocate(B)

     ! Create the plan
     call dfftw_plan_dft_3d(plan,nx,nx,nx,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)

     ! Compute the derivatives - X
     !$OMP PARALLEL DO PRIVATE(j,i,kx)
     do k=1,nx
        do j=1,nx
           do i=1,nx
              kx=real(i-1,WP)*dk
              if (i.gt.(nx/2+1)) kx=-real(nx-i+1,WP)*dk
              in(i,j,k) = ii*kx*vel(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call dfftw_execute(plan)
     !$OMP PARALLEL DO PRIVATE(j,i) REDUCTION(+:epsilon)
     do k=1,nx
        do j=1,nx
           do i=1,nx
              epsilon = epsilon + real(out(i,j,k),WP)**2
           end do
        end do
     end do
     !$OMP END PARALLEL DO

     ! Save (du1/dx1)^2*nx^3 for Taylor scale
     if (ivar.eq.1) lambda = epsilon

     ! Compute the derivatives - Y
     !$OMP PARALLEL DO PRIVATE(j,i,ky)
     do k=1,nx
        do j=1,nx
           do i=1,nx
              ky=real(j-1,WP)*dk
              if (j.gt.(nx/2+1)) ky=-real(nx-j+1,WP)*dk
              in(i,j,k) = ii*ky*vel(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call dfftw_execute(plan)
     !$OMP PARALLEL DO PRIVATE(j,i) REDUCTION(+:epsilon)
     do k=1,nx
        do j=1,nx
           do i=1,nx
              epsilon = epsilon + real(out(i,j,k),WP)**2
           end do
        end do
     end do
     !$OMP END PARALLEL DO

     ! Compute the derivatives - Z
     !$OMP PARALLEL DO PRIVATE(j,i,kz)
     do k=1,nx
        do j=1,nx
           do i=1,nx
              kz=real(k-1,WP)*dk
              if (k.gt.(nx/2+1)) kz=-real(nx-k+1,WP)*dk
              in(i,j,k) = ii*kz*vel(i,j,k)
           end do
        end do
     end do
     !$OMP END PARALLEL DO
     call dfftw_execute(plan)
     !$OMP PARALLEL DO PRIVATE(j,i) REDUCTION(+:epsilon)
     do k=1,nx
        do j=1,nx
           do i=1,nx
              epsilon = epsilon + real(out(i,j,k),WP)**2
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if
  
  ! Clean up
  call dfftw_destroy_plan(plan)
  
  return
end subroutine get_spectrum


! ============================================================ !
! Compute the two-point longitudinal velocity correlation      !
!   along coordinate direction ii                              !
! ============================================================ !
subroutine get_long_correlation(ii,A1,A2,S)
  use dnsbox_corr
  use precision
  implicit none
  include 'fftw3.f03'  
  integer, intent(in) :: ii
  complex*16, dimension(nx,nx,nx), intent(in) :: A1, A2
  real(WP), dimension(nx), intent(out) :: S
  integer(KIND=8) :: plan
  integer :: i,j,k
  !complex*16, parameter :: ii =(0.0_WP,1.0_WP)
  
  ! Create the plan
  call dfftw_plan_dft_3d(plan,nx,nx,nx,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
  
  ! Compute the inverse transform
  !$OMP PARALLEL DO PRIVATE(j,i)
  do k=1,nx
     do j=1,nx
        do i=1,nx
           ! Velocity cross-spectrum tensor
           in(i,j,k) = A1(i,j,k)*conjg(A2(i,j,k))
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  call dfftw_execute(plan)

  ! Two-point velocity correlation
  select case(ii)
  case(1)
     S = real(out(:,1,1),WP)
  case(2)
     S = real(out(1,:,1),WP)
  case(3)
     S = real(out(1,1,:),WP)
  end select
  
  ! Clean up
  call dfftw_destroy_plan(plan)

  return
end subroutine get_long_correlation
