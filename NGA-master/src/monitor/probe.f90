module probe
  use precision
  use geometry
  use partition
  implicit none
  
  ! Number of probes and variables to probe
  integer :: nprobefiles
  integer :: nprobes
  integer :: nvars
  integer :: ncols
  character(len=str_medium), dimension(:), allocatable :: varnames
  character(len=str_medium), dimension(:), pointer :: probe_filename
  integer, dimension(:), pointer :: lengths ! number of probes in each file
  real(WP), dimension(:,:,:,:), allocatable :: varsdata 

  ! Probe locations
  real(WP), dimension(:), pointer :: probe_x,probe_y,probe_z
  logical,  dimension(:), pointer :: probe_local
  
  ! Interpolation
  real(WP), dimension(:), pointer :: interp_x,interp_xm
  real(WP), dimension(:), pointer :: interp_y,interp_ym
  real(WP), dimension(:), pointer :: interp_z,interp_zm
  integer,  dimension(:), pointer :: index_x,index_xm
  integer,  dimension(:), pointer :: index_y,index_ym
  integer,  dimension(:), pointer :: index_z,index_zm
  
  ! Values to monitor
  real(WP), dimension(:), pointer :: mval
  logical :: probe_spectra

  ! Information for z-line probe output
  logical :: probe_line
  ! DOESN'T WORK YET FOR DECOMPOSITION IN Z!!
  character(len=str_medium), dimension(:), pointer :: probe_line_foldernames
  character(len=str_medium), dimension(:,:), pointer :: probe_line_filenames
  real(WP), dimension(:,:), pointer :: mval_line
  integer, dimension(:), pointer :: indx_pfolder
  integer, dimension(:), pointer :: indx_pfile

end module probe


! ===================== !
! Initialize the probes !
! ===================== !
subroutine probe_init
  use probe
  use parser
  use data
  implicit none
  integer :: nargs,n,i,j,k,isc,m,mm,ipfile,ipfolder,n_indx
  logical :: isdef
  logical :: istrue = .false.
  real(WP), dimension(:,:), pointer :: list
  integer ::  max_probe_files = 64
  real(WP), dimension(:),   pointer :: xvals 
  real(WP), dimension(:,:), pointer :: yvals 
  character(len=str_medium) :: cmd, yloc

  ! Any probes ?
  call parser_is_defined('Probe locations',isdef)
  if (isdef) then
     call parser_getsize('Probe locations',nargs)
     if (mod(nargs,3).ne.0) then
        call die('probe_init: incorrect number of coordinates for probes locations')
     else
        nprobes = nargs/3
     end if
  else
     nprobes = 0
     return
  end if

  ! Additional variable to probe
  call parser_is_defined('Other probe variables',isdef)
  if (isdef) then
     call parser_getsize('Other probe variables',nvars)
     allocate(varnames(nvars))
     call parser_read('Other probe variables',varnames)
     allocate(varsdata(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nvars))
  else
     nvars = 0
  end if

  ! Allocate arrays
  allocate(list(3,nprobes))
  allocate(probe_x(nprobes))
  allocate(probe_y(nprobes))
  allocate(probe_z(nprobes))
  allocate(probe_local(nprobes))
  allocate(interp_x (nprobes))
  allocate(interp_xm(nprobes))
  allocate(interp_y (nprobes))
  allocate(interp_ym(nprobes))
  allocate(interp_z (nprobes))
  allocate(interp_zm(nprobes))
  allocate(index_x (nprobes))
  allocate(index_xm(nprobes))
  allocate(index_y (nprobes))
  allocate(index_ym(nprobes))
  allocate(index_z (nprobes))
  allocate(index_zm(nprobes))
  allocate(lengths(max_probe_files))
  allocate(xvals(max_probe_files))

  ! Read the locations
  call parser_read('Probe locations',list)
  probe_x = list(1,:)
  probe_y = list(2,:)
  probe_z = list(3,:)

  ! Output strain rate?
  call parser_read('Probe spectra',probe_spectra,.false.)

  ! Output z-line probes?
  call parser_read('Probe z-lines',probe_line,.false.)
   
  ! Separate probe files 
  ! NOTE: In order to combine probes on x-location in files, need to list
  !       probes with the same x-location sequentially in the input file
  ! (This could be better -- add search function to allow arbitrary ordering)
  nprobefiles = 1
  lengths     = 0
  lengths(1)  = nprobes
  xvals(1)    = 0.0_WP
  call parser_is_defined('Separate probe files',isdef)
  if (isdef) then
     call parser_read('Separate probe files',istrue)
     if (istrue) then
        lengths(1) = 0  ! Number of probes saved in this file
        xvals(1) = list(1,1)
        do n=1,nprobes
           if (list(1,n).ne.xvals(nprobefiles)) then
              nprobefiles = nprobefiles + 1
              if (nprobefiles.gt.max_probe_files) then
                 call die('probe_init: too many probe files (x-locations)')
              end if
              xvals(nprobefiles) = list(1,n)
              lengths(nprobefiles) = 1
           else
              lengths(nprobefiles) = lengths(nprobefiles) + 1
           end if 
        end do
     end if
  end if

  ! Gather y-locations
  if (probe_line) then
     if (npz.gt.1) then
        call die('probe_init: probe_line not yet implemented for decomposition in z!')
     end if
     allocate(yvals(max_probe_files,max_probe_files))
     allocate(indx_pfolder(nprobes))
     allocate(indx_pfile  (nprobes))
     n_indx = 0
     do i=1,nprobefiles
        ! Loop over the number of probe files (x-locations)
        do j=1,lengths(i)
           ! Loop over the y-locations for this x-location
           !    Assumes properly formatted list (see above!)
           n_indx = n_indx + 1
           yvals(i,j) = list(2,n_indx)

           ! Save the folder/file indexing
           indx_pfolder(n_indx) = i
           indx_pfile  (n_indx) = j
        end do
     end do
  end if

  ! Get interpolation
  do n=1,nprobes
     if (probe_x(n).lt.x(imin) .or. probe_x(n).ge.x(imax+1)) &
          call die('probe_init: probe x location not in domain')
     if (probe_y(n).lt.y(jmin) .or. probe_y(n).ge.y(jmax+1)) &
          call die('probe_init: probe y location not in domain')
     if (probe_z(n).lt.z(kmin) .or. probe_z(n).ge.z(kmax+1)) &
          call die('probe_init: probe z location not in domain')
     
     ! Determine if a probe is on the current CPU
     if ( probe_x(n).lt.x(imin_) .or. probe_x(n).ge.x(imax_+1) .or. &
          probe_y(n).lt.y(jmin_) .or. probe_y(n).ge.y(jmax_+1) .or. &
          probe_z(n).lt.z(kmin_) .or. probe_z(n).ge.z(kmax_+1)) then
        probe_local(n) = .false.
     else
        probe_local(n) = .true.

        ! Interpolation in x
        i = imin_
        do while(x(i+1).le.probe_x(n)) 
           i = i+1
        end do
        index_x(n)  = i
        interp_x(n) = (probe_x(n)-x(i))/(x(i+1)-x(i))
        i = imino_
        do while(xm(i+1).le.probe_x(n)) 
           i = i+1
        end do
        index_xm(n)  = i
        interp_xm(n) = (probe_x(n)-xm(i))/(xm(i+1)-xm(i))
        
        ! Interpolation in y
        j = jmin_
        do while(y(j+1).le.probe_y(n)) 
           j = j+1
        end do
        index_y(n)  = j
        interp_y(n) = (probe_y(n)-y(j))/(y(j+1)-y(j))
        j = jmino_
        do while(ym(j+1).le.probe_y(n)) 
           j = j+1
        end do
        index_ym(n)  = j
        interp_ym(n) = (probe_y(n)-ym(j))/(ym(j+1)-ym(j))
        
        ! Interpolation in z
        k = kmin_
        do while(z(k+1).le.probe_z(n)) 
           k = k+1
        end do
        index_z(n)  = k
        interp_z(n) = (probe_z(n)-z(k))/(z(k+1)-z(k))
        k = kmino_
        do while(zm(k+1).le.probe_z(n)) 
           k = k+1
        end do
        index_zm(n)  = k
        interp_zm(n) = (probe_z(n)-zm(k))/(zm(k+1)-zm(k))
     end if
  end do
  
  ! Create files to monitor - single-point probes
  allocate(probe_filename(nprobefiles))
  do ipfile=1,nprobefiles
     write(probe_filename(ipfile),'(a9,es12.6e2)') 'probes-x-',xvals(ipfile)
     if ( (ipfile.eq.1) .and. (.not.(istrue.and.isdef)) ) probe_filename(1)='probes'
     ! Count the number of columns per probe
     if (probe_spectra) then
        ! 3x velocity components
        ! Pressure
        ! Scalars
        ! Additional vars
        ! Density
        ! Viscosity
        ! 9x velocity field gradients
        ! 3x density gradient
        ! 3x hydrodynamic pressure gradient
        ! 9x stress tensor gradients
        ! Heat release rate
        ncols = 4 + nscalar + nvars + 2 + 9 + 3 + 3 + 9 + 1
     else
        ncols = 4 + nscalar + nvars
     end if
     call monitor_create_file_step(probe_filename(ipfile),ncols*lengths(ipfile))
     do n=1,lengths(ipfile)
        call monitor_set_header(1+(n-1)*ncols,'U','d')
        call monitor_set_header(2+(n-1)*ncols,'V','d')
        call monitor_set_header(3+(n-1)*ncols,'W','d')
        call monitor_set_header(4+(n-1)*ncols,'P','d')
        do isc=1,nscalar
           call monitor_set_header(4+isc+(n-1)*ncols,trim(SC_name(isc)),'d')
        end do
        do m=1,nvars
           call monitor_set_header(4+nscalar+m+(n-1)*ncols,trim(varnames(m)),'d')
        end do
        ! Additional information required to compute spectra
        if (probe_spectra) then
           call monitor_set_header(5 +nvars+nscalar+(n-1)*ncols,'RHO', 'd')
           call monitor_set_header(6 +nvars+nscalar+(n-1)*ncols,'VISC','d')
           call monitor_set_header( 7+nvars+nscalar+(n-1)*ncols,'dUdx11','d')
           call monitor_set_header( 8+nvars+nscalar+(n-1)*ncols,'dUdx12','d')
           call monitor_set_header( 9+nvars+nscalar+(n-1)*ncols,'dUdx13','d')
           call monitor_set_header(10+nvars+nscalar+(n-1)*ncols,'dUdx21','d')
           call monitor_set_header(11+nvars+nscalar+(n-1)*ncols,'dUdx22','d')
           call monitor_set_header(12+nvars+nscalar+(n-1)*ncols,'dUdx23','d')
           call monitor_set_header(13+nvars+nscalar+(n-1)*ncols,'dUdx31','d')
           call monitor_set_header(14+nvars+nscalar+(n-1)*ncols,'dUdx32','d')
           call monitor_set_header(15+nvars+nscalar+(n-1)*ncols,'dUdx33','d')
           call monitor_set_header(16+nvars+nscalar+(n-1)*ncols,'dRHOdx1','d')
           call monitor_set_header(17+nvars+nscalar+(n-1)*ncols,'dRHOdx2','d')
           call monitor_set_header(18+nvars+nscalar+(n-1)*ncols,'dRHOdx3','d')
           call monitor_set_header(19+nvars+nscalar+(n-1)*ncols,'dPdx1','d')
           call monitor_set_header(20+nvars+nscalar+(n-1)*ncols,'dPdx2','d')
           call monitor_set_header(21+nvars+nscalar+(n-1)*ncols,'dPdx3','d')
           call monitor_set_header(22+nvars+nscalar+(n-1)*ncols,'dTAUdx11','d')
           call monitor_set_header(23+nvars+nscalar+(n-1)*ncols,'dTAUdx12','d')
           call monitor_set_header(24+nvars+nscalar+(n-1)*ncols,'dTAUdx13','d')
           call monitor_set_header(25+nvars+nscalar+(n-1)*ncols,'dTAUdx21','d')
           call monitor_set_header(26+nvars+nscalar+(n-1)*ncols,'dTAUdx22','d')
           call monitor_set_header(27+nvars+nscalar+(n-1)*ncols,'dTAUdx23','d')
           call monitor_set_header(28+nvars+nscalar+(n-1)*ncols,'dTAUdx31','d')
           call monitor_set_header(29+nvars+nscalar+(n-1)*ncols,'dTAUdx32','d')
           call monitor_set_header(30+nvars+nscalar+(n-1)*ncols,'dTAUdx33','d')
           call monitor_set_header(31+nvars+nscalar+(n-1)*ncols,'HR','d')
        end if
     end do
  end do
  allocate(mval(ncols*nprobes))


  ! Create files to monitor - z-line probes
  if (probe_line) then
     allocate(probe_line_foldernames(nprobefiles))
     allocate(probe_line_filenames(nprobefiles,ncols*maxval(lengths(:))))
     
     do ipfolder=1,nprobefiles
        ! One folder for each x-location
        ! Generate the folder name
        write(probe_line_foldernames(ipfolder), '(a10,es12.6e2)') 'zprobes-x-',xvals(ipfolder)
        ! Make the folder
        if (irank.eq.iroot) then
           write(cmd, '(a27,es12.6e2)') 'mkdir -p monitor/zprobes-x-',xvals(ipfolder)
           call system( cmd )
        end if

        do ipfile=1,lengths(ipfolder)
           ! One file for each variable at each y-location
           write(yloc, '(es12.6e2)') yvals(ipfolder,ipfile)

           ! Initialize the files
           ! U
           probe_line_filenames(ipfolder, 1+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-U'
           call monitor_create_file_step(probe_line_filenames(ipfolder, 1 +(ipfile-1)*ncols), nz)
           do k=1,nz
              write(cmd, '(a3,i8)') 'nz-',k
              call monitor_set_header(k, cmd, 'd')
           end do
           ! V
           probe_line_filenames(ipfolder, 2+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-V'
           call monitor_create_file_step(probe_line_filenames(ipfolder, 2 +(ipfile-1)*ncols), nz)
           do k=1,nz
              write(cmd, '(a3,i8)') 'nz-',k
              call monitor_set_header(k, cmd, 'd')
           end do
           ! W
           probe_line_filenames(ipfolder, 3+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-W'
           call monitor_create_file_step(probe_line_filenames(ipfolder, 3 +(ipfile-1)*ncols), nz)
           do k=1,nz
              write(cmd, '(a3,i8)') 'nz-',k
              call monitor_set_header(k, cmd, 'd')
           end do
           ! P
           probe_line_filenames(ipfolder, 4+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-P'
           call monitor_create_file_step(probe_line_filenames(ipfolder, 4 +(ipfile-1)*ncols), nz)
           do k=1,nz
              write(cmd, '(a3,i8)') 'nz-',k
              call monitor_set_header(k, cmd, 'd')
           end do
           ! SC
           do isc=1,nscalar
              probe_line_filenames(ipfolder, 4+isc+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-' // trim(SC_name(isc))
              call monitor_create_file_step(probe_line_filenames(ipfolder, 4+isc+(ipfile-1)*ncols), nz)
              do k=1,nz
                 write(cmd, '(a3,i8)') 'nz-',k
                 call monitor_set_header(k, cmd, 'd')
              end do
           end do
           ! NVARS
           do m=1,nvars
              probe_line_filenames(ipfolder, 4+nscalar+m+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-' // trim(varnames(m))
              call monitor_create_file_step(probe_line_filenames(ipfolder, 4+nscalar+m+(ipfile-1)*ncols), nz)
              do k=1,nz
                 write(cmd, '(a3,i8)') 'nz-',k
                 call monitor_set_header(k, cmd, 'd')
              end do
           end do
           
           ! Data for spectra
           if (probe_spectra) then
              ! RHO
              probe_line_filenames(ipfolder, 5+nvars+nscalar+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-RHO'
              call monitor_create_file_step(probe_line_filenames(ipfolder, 5+nvars+nscalar+(ipfile-1)*ncols), nz)
              do k=1,nz
                 write(cmd, '(a3,i8)') 'nz-',k
                 call monitor_set_header(k, cmd, 'd')
              end do
              ! VISC
              probe_line_filenames(ipfolder, 6+nvars+nscalar+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-VISC'
              call monitor_create_file_step(probe_line_filenames(ipfolder, 6+nvars+nscalar+(ipfile-1)*ncols), nz)
              do k=1,nz
                 write(cmd, '(a3,i8)') 'nz-',k
                 call monitor_set_header(k, cmd, 'd')
              end do
              ! Velocity gradients
              mm = 1
              do m=1,3
                 do n=1,3
                    write(cmd, '(i1,i1)') m, n
                    probe_line_filenames(ipfolder, 6+mm+nvars+nscalar+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-dUdx' // trim(cmd)
                    call monitor_create_file_step(probe_line_filenames(ipfolder, 6+mm+nvars+nscalar+(ipfile-1)*ncols), nz)
                    do k=1,nz
                       write(cmd, '(a3,i8)') 'nz-',k
                       call monitor_set_header(k, cmd, 'd')
                    end do
                    mm = mm + 1
                 end do
              end do
              ! Density gradient
              do m=1,3
                 write(cmd, '(i1)') m
                 probe_line_filenames(ipfolder, 15+m+nvars+nscalar+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-dRHOdx' // trim(cmd)
                 call monitor_create_file_step(probe_line_filenames(ipfolder, 15+m+nvars+nscalar+(ipfile-1)*ncols), nz)
                 do k=1,nz
                    write(cmd, '(a3,i8)') 'nz-',k
                    call monitor_set_header(k, cmd, 'd')
                 end do
              end do
              ! Pressure gradient
              do m=1,3
                 write(cmd, '(i1)') m
                 probe_line_filenames(ipfolder, 18+m+nvars+nscalar+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-dPdx' // trim(cmd)
                 call monitor_create_file_step(probe_line_filenames(ipfolder, 18+m+nvars+nscalar+(ipfile-1)*ncols), nz)
                 do k=1,nz
                    write(cmd, '(a3,i8)') 'nz-',k
                    call monitor_set_header(k, cmd, 'd')
                 end do
              end do
              ! Stress tensor gradients
              mm = 1
              do m=1,3
                 do n=1,3
                    write(cmd, '(i1,i1)') m, n
                    probe_line_filenames(ipfolder, 21+mm+nvars+nscalar+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-dTAUdx' // trim(cmd)
                    call monitor_create_file_step(probe_line_filenames(ipfolder, 21+mm+nvars+nscalar+(ipfile-1)*ncols), nz)
                    do k=1,nz
                       write(cmd, '(a3,i8)') 'nz-',k
                       call monitor_set_header(k, cmd, 'd')
                    end do
                    mm = mm + 1
                 end do
              end do
              ! RHO
              probe_line_filenames(ipfolder, 31+nvars+nscalar+(ipfile-1)*ncols) = trim(probe_line_foldernames(ipfolder)) // '/y-' // trim(yloc) // '-HR'
              call monitor_create_file_step(probe_line_filenames(ipfolder, 31+nvars+nscalar+(ipfile-1)*ncols), nz)
              do k=1,nz
                 write(cmd, '(a3,i8)') 'nz-',k
                 call monitor_set_header(k, cmd, 'd')
              end do

           end if
           ! Done creating files
        end do
     end do

     ! mval_line: 1st index = variables for each probe
     !            2nd index = z-location
     allocate(mval_line(ncols,kmin:kmax))
  end if

  return
end subroutine probe_init


! ================== !
! Monitor the probes !
! ================== !
subroutine probe_monitor
  use probe
  use data
  use strainrate
  implicit none
  integer  :: n,isc,i,j,k,m,ipfile,nprev,nk,ii,jj,kk
  real(WP) :: wx1,wy1,wz1,wx2,wy2,wz2,tmp
  real(WP), dimension(3,3,8) :: dUdx
  real(WP), dimension(3,8)   :: dSCdx
  real(WP), dimension(3,3,-st2:1+st1,-st2:1+st1,-st2:1+st1) :: dUdxs
  
  ! Nothing to do if no probes
  if (nprobes.eq.0) return

  ! Get chemtable data for relevant vairables 
  if (nvars.ne.0) then
     do n=1,nvars
        call chemtable_lookup(trim(varnames(n)), varsdata(:,:,:,n))
     end do
  end if

  ! Compute local values at the probes
  mval = 0.0_WP
  do n=1,nprobes
     if (probe_local(n)) then
        ! Velocity - U
        i = index_x (n); wx2 = interp_x (n); wx1 = 1.0_WP-wx2
        j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
        k = index_zm(n); wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
        mval(1+(n-1)*(ncols)) = &
             + wx1 * ( wy1 * (wz1*U(i  ,j  ,k)+wz2*U(i  ,j  ,k+1)) &
                     + wy2 * (wz1*U(i  ,j+1,k)+wz2*U(i  ,j+1,k+1)) )&
             + wx2 * ( wy1 * (wz1*U(i+1,j  ,k)+wz2*U(i+1,j  ,k+1)) &
                     + wy2 * (wz1*U(i+1,j+1,k)+wz2*U(i+1,j+1,k+1)) )
        ! Velocity - V
        i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
        j = index_y (n); wy2 = interp_y (n); wy1 = 1.0_WP-wy2
        k = index_zm(n); wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
        mval(2+(n-1)*(ncols)) = &
             + wx1 * ( wy1 * (wz1*V(i  ,j  ,k)+wz2*V(i  ,j  ,k+1)) &
                     + wy2 * (wz1*V(i  ,j+1,k)+wz2*V(i  ,j+1,k+1)) )&
             + wx2 * ( wy1 * (wz1*V(i+1,j  ,k)+wz2*V(i+1,j  ,k+1)) &
                     + wy2 * (wz1*V(i+1,j+1,k)+wz2*V(i+1,j+1,k+1)) )
        ! Velocity - W
        i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
        j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
        k = index_z (n); wz2 = interp_z (n); wz1 = 1.0_WP-wz2
        mval(3+(n-1)*(ncols)) = &
             + wx1 * ( wy1 * (wz1*W(i  ,j  ,k)+wz2*W(i  ,j  ,k+1)) &
                     + wy2 * (wz1*W(i  ,j+1,k)+wz2*W(i  ,j+1,k+1)) )&
             + wx2 * ( wy1 * (wz1*W(i+1,j  ,k)+wz2*W(i+1,j  ,k+1)) &
                     + wy2 * (wz1*W(i+1,j+1,k)+wz2*W(i+1,j+1,k+1)) )
        ! Pressure
        i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
        j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
        k = index_zm(n); wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
        mval(4+(n-1)*(ncols)) = &
             + wx1 * ( wy1 * (wz1*P(i  ,j  ,k)+wz2*P(i  ,j  ,k+1)) &
                     + wy2 * (wz1*P(i  ,j+1,k)+wz2*P(i  ,j+1,k+1)) )&
             + wx2 * ( wy1 * (wz1*P(i+1,j  ,k)+wz2*P(i+1,j  ,k+1)) &
                     + wy2 * (wz1*P(i+1,j+1,k)+wz2*P(i+1,j+1,k+1)) )
        ! Scalars
        i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
        j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
        k = index_zm(n); wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
        do isc=1,nscalar
           mval(4+isc+(n-1)*(ncols)) = &
                + wx1 * ( wy1 * (wz1*SC(i  ,j  ,k,isc)+wz2*SC(i  ,j  ,k+1,isc)) &
                        + wy2 * (wz1*SC(i  ,j+1,k,isc)+wz2*SC(i  ,j+1,k+1,isc)) )&
                + wx2 * ( wy1 * (wz1*SC(i+1,j  ,k,isc)+wz2*SC(i+1,j  ,k+1,isc)) &
                        + wy2 * (wz1*SC(i+1,j+1,k,isc)+wz2*SC(i+1,j+1,k+1,isc)) )
        end do
        ! Other Variables
        do m=1,nvars
           mval(4+nscalar+m+(n-1)*(ncols)) = &
                + wx1 * ( wy1 * (wz1*varsdata(i  ,j  ,k,m)+wz2*varsdata(i  ,j  ,k+1,m)) &
                        + wy2 * (wz1*varsdata(i  ,j+1,k,m)+wz2*varsdata(i  ,j+1,k+1,m)) )&
                + wx2 * ( wy1 * (wz1*varsdata(i+1,j  ,k,m)+wz2*varsdata(i+1,j  ,k+1,m)) &
                        + wy2 * (wz1*varsdata(i+1,j+1,k,m)+wz2*varsdata(i+1,j+1,k+1,m)) )
        end do

        ! Information for spectra
        if (probe_spectra) then
           ! RHO(i,j,k)-> xm(i),ym(j),zm(k) at t=n+1/2
           i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
           j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
           k = index_zm(n); wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
           mval(5 +nvars+nscalar+(n-1)*(ncols)) = &
                + wx1 * ( wy1 * (wz1*RHO(i  ,j  ,k)+wz2*RHO(i  ,j  ,k+1)) &
                        + wy2 * (wz1*RHO(i  ,j+1,k)+wz2*RHO(i  ,j+1,k+1)) )&
                + wx2 * ( wy1 * (wz1*RHO(i+1,j  ,k)+wz2*RHO(i+1,j  ,k+1)) &
                        + wy2 * (wz1*RHO(i+1,j+1,k)+wz2*RHO(i+1,j+1,k+1)) )

           ! Viscosity
           mval(6 +nvars+nscalar+(n-1)*ncols) = &
                + wx1 * ( wy1 * (wz1*VISC(i  ,j  ,k)+wz2*VISC(i  ,j  ,k+1)) &
                        + wy2 * (wz1*VISC(i  ,j+1,k)+wz2*VISC(i  ,j+1,k+1)) )&
                + wx2 * ( wy1 * (wz1*VISC(i+1,j  ,k)+wz2*VISC(i+1,j  ,k+1)) &
                        + wy2 * (wz1*VISC(i+1,j+1,k)+wz2*VISC(i+1,j+1,k+1)) )

           ! Velocity gradients
           ! dUdx computed at cell centers using i,j,k indices at cell faces
           i = index_x(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
           j = index_y(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
           k = index_z(n); wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
           call vel_grad_local(i  ,j  ,k  , dUdx(:,:,1))
           call vel_grad_local(i  ,j  ,k+1, dUdx(:,:,2))
           call vel_grad_local(i  ,j+1,k  , dUdx(:,:,3))
           call vel_grad_local(i  ,j+1,k+1, dUdx(:,:,4))
           call vel_grad_local(i+1,j  ,k  , dUdx(:,:,5))
           call vel_grad_local(i+1,j  ,k+1, dUdx(:,:,6))
           call vel_grad_local(i+1,j+1,k  , dUdx(:,:,7))
           call vel_grad_local(i+1,j+1,k+1, dUdx(:,:,8))
           m = 7
           do ii=1,3
              do jj=1,3
                 mval(m +nvars+nscalar+(n-1)*ncols) = &
                      + wx1 * ( wy1 * (wz1*dUdx(ii,jj,1)+wz2*dUdx(ii,jj,2)) &
                              + wy2 * (wz1*dUdx(ii,jj,3)+wz2*dUdx(ii,jj,4)) )&
                      + wx2 * ( wy1 * (wz1*dUdx(ii,jj,5)+wz2*dUdx(ii,jj,6)) &
                              + wy2 * (wz1*dUdx(ii,jj,7)+wz2*dUdx(ii,jj,8)) )
                 m = m + 1
              end do
           end do

           ! Density gradient
           ! dRHOdx computed at cell faces using i,j,k indices at cell centers
           i = index_xm(n)
           j = index_ym(n)
           k = index_zm(n)
           dSCdx(1,1) = sum(grad_x(i  ,j  ,:)*RHO(i  -st2:i  +st1,j  ,k  ))
           dSCdx(1,2) = sum(grad_x(i  ,j  ,:)*RHO(i  -st2:i  +st1,j  ,k+1))
           dSCdx(1,3) = sum(grad_x(i  ,j+1,:)*RHO(i  -st2:i  +st1,j+1,k  ))
           dSCdx(1,4) = sum(grad_x(i  ,j+1,:)*RHO(i  -st2:i  +st1,j+1,k+1))
           dSCdx(1,5) = sum(grad_x(i+1,j  ,:)*RHO(i+1-st2:i+1+st1,j  ,k  ))
           dSCdx(1,6) = sum(grad_x(i+1,j  ,:)*RHO(i+1-st2:i+1+st1,j  ,k+1))
           dSCdx(1,7) = sum(grad_x(i+1,j+1,:)*RHO(i+1-st2:i+1+st1,j+1,k  ))
           dSCdx(1,8) = sum(grad_x(i+1,j+1,:)*RHO(i+1-st2:i+1+st1,j+1,k+1))
           dSCdx(2,1) = sum(grad_y(i  ,j  ,:)*RHO(i  ,j  -st2:j  +st1,k  ))
           dSCdx(2,2) = sum(grad_y(i  ,j  ,:)*RHO(i  ,j  -st2:j  +st1,k+1))
           dSCdx(2,3) = sum(grad_y(i  ,j+1,:)*RHO(i  ,j+1-st2:j+1+st1,k  ))
           dSCdx(2,4) = sum(grad_y(i  ,j+1,:)*RHO(i  ,j+1-st2:j+1+st1,k+1))
           dSCdx(2,5) = sum(grad_y(i+1,j  ,:)*RHO(i+1,j  -st2:j  +st1,k  ))
           dSCdx(2,6) = sum(grad_y(i+1,j  ,:)*RHO(i+1,j  -st2:j  +st1,k+1))
           dSCdx(2,7) = sum(grad_y(i+1,j+1,:)*RHO(i+1,j+1-st2:j+1+st1,k  ))
           dSCdx(2,8) = sum(grad_y(i+1,j+1,:)*RHO(i+1,j+1-st2:j+1+st1,k+1))
           dSCdx(3,1) = sum(grad_z(i  ,j  ,:)*RHO(i  ,j  ,k  -st2:k  +st1))
           dSCdx(3,2) = sum(grad_z(i  ,j  ,:)*RHO(i  ,j  ,k+1-st2:k+1+st1))
           dSCdx(3,3) = sum(grad_z(i  ,j+1,:)*RHO(i  ,j+1,k  -st2:k  +st1))
           dSCdx(3,4) = sum(grad_z(i  ,j+1,:)*RHO(i  ,j+1,k+1-st2:k+1+st1))
           dSCdx(3,5) = sum(grad_z(i+1,j  ,:)*RHO(i+1,j  ,k  -st2:k  +st1))
           dSCdx(3,6) = sum(grad_z(i+1,j  ,:)*RHO(i+1,j  ,k+1-st2:k+1+st1))
           dSCdx(3,7) = sum(grad_z(i+1,j+1,:)*RHO(i+1,j+1,k  -st2:k  +st1))
           dSCdx(3,8) = sum(grad_z(i+1,j+1,:)*RHO(i+1,j+1,k+1-st2:k+1+st1))
           m = 16
           do ii=1,3
              if (ii.eq.1) then
                 wx2 = interp_x(n);  wx1 = 1.0_WP-wx2
                 wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                 wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
              elseif (ii.eq.2) then
                 wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                 wy2 = interp_y(n);  wy1 = 1.0_WP-wy2
                 wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
              elseif (ii.eq.3) then
                 wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                 wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                 wz2 = interp_z(n);  wz1 = 1.0_WP-wz2
              end if
              mval(m +nvars+nscalar+(n-1)*(ncols)) = &
                   + wx1 * ( wy1 * (wz1*dSCdx(ii,1)+wz2*dSCdx(ii,2)) &
                           + wy2 * (wz1*dSCdx(ii,3)+wz2*dSCdx(ii,4)) )&
                   + wx2 * ( wy1 * (wz1*dSCdx(ii,5)+wz2*dSCdx(ii,6)) &
                           + wy2 * (wz1*dSCdx(ii,7)+wz2*dSCdx(ii,8)) )
              m = m + 1
           end do

           ! Pressure gradient
           ! dPdx computed at cell faces using i,j,k indices at cell centers
           dSCdx(1,1) = sum(grad_x(i  ,j  ,:)*P(i  -st2:i  +st1,j  ,k  ))
           dSCdx(1,2) = sum(grad_x(i  ,j  ,:)*P(i  -st2:i  +st1,j  ,k+1))
           dSCdx(1,3) = sum(grad_x(i  ,j+1,:)*P(i  -st2:i  +st1,j+1,k  ))
           dSCdx(1,4) = sum(grad_x(i  ,j+1,:)*P(i  -st2:i  +st1,j+1,k+1))
           dSCdx(1,5) = sum(grad_x(i+1,j  ,:)*P(i+1-st2:i+1+st1,j  ,k  ))
           dSCdx(1,6) = sum(grad_x(i+1,j  ,:)*P(i+1-st2:i+1+st1,j  ,k+1))
           dSCdx(1,7) = sum(grad_x(i+1,j+1,:)*P(i+1-st2:i+1+st1,j+1,k  ))
           dSCdx(1,8) = sum(grad_x(i+1,j+1,:)*P(i+1-st2:i+1+st1,j+1,k+1))
           dSCdx(2,1) = sum(grad_y(i  ,j  ,:)*P(i  ,j  -st2:j  +st1,k  ))
           dSCdx(2,2) = sum(grad_y(i  ,j  ,:)*P(i  ,j  -st2:j  +st1,k+1))
           dSCdx(2,3) = sum(grad_y(i  ,j+1,:)*P(i  ,j+1-st2:j+1+st1,k  ))
           dSCdx(2,4) = sum(grad_y(i  ,j+1,:)*P(i  ,j+1-st2:j+1+st1,k+1))
           dSCdx(2,5) = sum(grad_y(i+1,j  ,:)*P(i+1,j  -st2:j  +st1,k  ))
           dSCdx(2,6) = sum(grad_y(i+1,j  ,:)*P(i+1,j  -st2:j  +st1,k+1))
           dSCdx(2,7) = sum(grad_y(i+1,j+1,:)*P(i+1,j+1-st2:j+1+st1,k  ))
           dSCdx(2,8) = sum(grad_y(i+1,j+1,:)*P(i+1,j+1-st2:j+1+st1,k+1))
           dSCdx(3,1) = sum(grad_z(i  ,j  ,:)*P(i  ,j  ,k  -st2:k  +st1))
           dSCdx(3,2) = sum(grad_z(i  ,j  ,:)*P(i  ,j  ,k+1-st2:k+1+st1))
           dSCdx(3,3) = sum(grad_z(i  ,j+1,:)*P(i  ,j+1,k  -st2:k  +st1))
           dSCdx(3,4) = sum(grad_z(i  ,j+1,:)*P(i  ,j+1,k+1-st2:k+1+st1))
           dSCdx(3,5) = sum(grad_z(i+1,j  ,:)*P(i+1,j  ,k  -st2:k  +st1))
           dSCdx(3,6) = sum(grad_z(i+1,j  ,:)*P(i+1,j  ,k+1-st2:k+1+st1))
           dSCdx(3,7) = sum(grad_z(i+1,j+1,:)*P(i+1,j+1,k  -st2:k  +st1))
           dSCdx(3,8) = sum(grad_z(i+1,j+1,:)*P(i+1,j+1,k+1-st2:k+1+st1))
           m = 19
           do ii=1,3
              if (ii.eq.1) then
                 wx2 = interp_x(n);  wx1 = 1.0_WP-wx2
                 wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                 wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
              elseif (ii.eq.2) then
                 wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                 wy2 = interp_y(n);  wy1 = 1.0_WP-wy2
                 wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
              elseif (ii.eq.3) then
                 wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                 wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                 wz2 = interp_z(n);  wz1 = 1.0_WP-wz2
              end if
              mval(m +nvars+nscalar+(n-1)*(ncols)) = &
                   + wx1 * ( wy1 * (wz1*dSCdx(ii,1)+wz2*dSCdx(ii,2)) &
                           + wy2 * (wz1*dSCdx(ii,3)+wz2*dSCdx(ii,4)) )&
                   + wx2 * ( wy1 * (wz1*dSCdx(ii,5)+wz2*dSCdx(ii,6)) &
                           + wy2 * (wz1*dSCdx(ii,7)+wz2*dSCdx(ii,8)) )
              m = m + 1
           end do

           ! Stress tensor gradients
           ! dUdx computed at cell centers using i,j,k indices at cell faces
           m = 22
           i = index_x(n)
           j = index_y(n)
           k = index_z(n)
           do ii=-st2,1+st1
              do jj=-st2,1+st1
                 do kk=-st2,1+st1
                    call vel_grad_local(i+ii,j+jj,k+kk, dUdxs(:,:,ii,jj,kk))
                 end do
              end do
           end do
           i = index_xm(n)
           j = index_ym(n)
           k = index_zm(n)
           do ii=1,3
              if (ii.eq.1) then
                 dSCdx(1,1) = sum(grad_x(i  ,j  ,:)*VISC(i  -st2:i  +st1,j  ,k  )*(dUdxs(ii,1, -st2: +st1,0,0) + dUdxs(1,ii, -st2: +st1,0,0) - 2.0_WP/3.0_WP*dUdxs(1,1, -st2: +st1,0,0)))
                 dSCdx(1,2) = sum(grad_x(i  ,j  ,:)*VISC(i  -st2:i  +st1,j  ,k+1)*(dUdxs(ii,1, -st2: +st1,0,1) + dUdxs(1,ii, -st2: +st1,0,1) - 2.0_WP/3.0_WP*dUdxs(1,1, -st2: +st1,0,1)))
                 dSCdx(1,3) = sum(grad_x(i  ,j+1,:)*VISC(i  -st2:i  +st1,j+1,k  )*(dUdxs(ii,1, -st2: +st1,1,0) + dUdxs(1,ii, -st2: +st1,1,0) - 2.0_WP/3.0_WP*dUdxs(1,1, -st2: +st1,1,0)))
                 dSCdx(1,4) = sum(grad_x(i  ,j+1,:)*VISC(i  -st2:i  +st1,j+1,k+1)*(dUdxs(ii,1, -st2: +st1,1,1) + dUdxs(1,ii, -st2: +st1,1,1) - 2.0_WP/3.0_WP*dUdxs(1,1, -st2: +st1,1,1)))
                 dSCdx(1,5) = sum(grad_x(i+1,j  ,:)*VISC(i+1-st2:i+1+st1,j  ,k  )*(dUdxs(ii,1,1-st2:1+st1,0,0) + dUdxs(1,ii,1-st2:1+st1,0,0) - 2.0_WP/3.0_WP*dUdxs(1,1,1-st2:1+st1,0,0)))
                 dSCdx(1,6) = sum(grad_x(i+1,j  ,:)*VISC(i+1-st2:i+1+st1,j  ,k+1)*(dUdxs(ii,1,1-st2:1+st1,0,1) + dUdxs(1,ii,1-st2:1+st1,0,1) - 2.0_WP/3.0_WP*dUdxs(1,1,1-st2:1+st1,0,1)))
                 dSCdx(1,7) = sum(grad_x(i+1,j+1,:)*VISC(i+1-st2:i+1+st1,j+1,k  )*(dUdxs(ii,1,1-st2:1+st1,1,0) + dUdxs(1,ii,1-st2:1+st1,1,0) - 2.0_WP/3.0_WP*dUdxs(1,1,1-st2:1+st1,1,0)))
                 dSCdx(1,8) = sum(grad_x(i+1,j+1,:)*VISC(i+1-st2:i+1+st1,j+1,k+1)*(dUdxs(ii,1,1-st2:1+st1,1,1) + dUdxs(1,ii,1-st2:1+st1,1,1) - 2.0_WP/3.0_WP*dUdxs(1,1,1-st2:1+st1,1,1)))
              else
                 dSCdx(1,1) = sum(grad_x(i  ,j  ,:)*VISC(i  -st2:i  +st1,j  ,k  )*(dUdxs(ii,1, -st2: +st1,0,0) + dUdxs(1,ii, -st2: +st1,0,0)))
                 dSCdx(1,2) = sum(grad_x(i  ,j  ,:)*VISC(i  -st2:i  +st1,j  ,k+1)*(dUdxs(ii,1, -st2: +st1,0,1) + dUdxs(1,ii, -st2: +st1,0,1)))
                 dSCdx(1,3) = sum(grad_x(i  ,j+1,:)*VISC(i  -st2:i  +st1,j+1,k  )*(dUdxs(ii,1, -st2: +st1,1,0) + dUdxs(1,ii, -st2: +st1,1,0)))
                 dSCdx(1,4) = sum(grad_x(i  ,j+1,:)*VISC(i  -st2:i  +st1,j+1,k+1)*(dUdxs(ii,1, -st2: +st1,1,1) + dUdxs(1,ii, -st2: +st1,1,1)))
                 dSCdx(1,5) = sum(grad_x(i+1,j  ,:)*VISC(i+1-st2:i+1+st1,j  ,k  )*(dUdxs(ii,1,1-st2:1+st1,0,0) + dUdxs(1,ii,1-st2:1+st1,0,0)))
                 dSCdx(1,6) = sum(grad_x(i+1,j  ,:)*VISC(i+1-st2:i+1+st1,j  ,k+1)*(dUdxs(ii,1,1-st2:1+st1,0,1) + dUdxs(1,ii,1-st2:1+st1,0,1)))
                 dSCdx(1,7) = sum(grad_x(i+1,j+1,:)*VISC(i+1-st2:i+1+st1,j+1,k  )*(dUdxs(ii,1,1-st2:1+st1,1,0) + dUdxs(1,ii,1-st2:1+st1,1,0)))
                 dSCdx(1,8) = sum(grad_x(i+1,j+1,:)*VISC(i+1-st2:i+1+st1,j+1,k+1)*(dUdxs(ii,1,1-st2:1+st1,1,1) + dUdxs(1,ii,1-st2:1+st1,1,1)))
              end if
              if (ii.eq.2) then
                 dSCdx(2,1) = sum(grad_y(i  ,j  ,:)*VISC(i  ,j  -st2:j  +st1,k  )*(dUdxs(ii,2,0, -st2: +st1,0) + dUdxs(2,ii,0, -st2: +st1,0) - 2.0_WP/3.0_WP*dUdxs(2,2,0, -st2: +st1,0)))
                 dSCdx(2,2) = sum(grad_y(i  ,j  ,:)*VISC(i  ,j  -st2:j  +st1,k+1)*(dUdxs(ii,2,0, -st2: +st1,1) + dUdxs(2,ii,0, -st2: +st1,1) - 2.0_WP/3.0_WP*dUdxs(2,2,0, -st2: +st1,1)))
                 dSCdx(2,3) = sum(grad_y(i  ,j+1,:)*VISC(i  ,j+1-st2:j+1+st1,k  )*(dUdxs(ii,2,0,1-st2:1+st1,0) + dUdxs(2,ii,0,1-st2:1+st1,0) - 2.0_WP/3.0_WP*dUdxs(2,2,0,1-st2:1+st1,0)))
                 dSCdx(2,4) = sum(grad_y(i  ,j+1,:)*VISC(i  ,j+1-st2:j+1+st1,k+1)*(dUdxs(ii,2,0,1-st2:1+st1,1) + dUdxs(2,ii,0,1-st2:1+st1,1) - 2.0_WP/3.0_WP*dUdxs(2,2,0,1-st2:1+st1,1)))
                 dSCdx(2,5) = sum(grad_y(i+1,j  ,:)*VISC(i+1,j  -st2:j  +st1,k  )*(dUdxs(ii,2,1, -st2: +st1,0) + dUdxs(2,ii,1, -st2: +st1,0) - 2.0_WP/3.0_WP*dUdxs(2,2,1, -st2: +st1,0)))
                 dSCdx(2,6) = sum(grad_y(i+1,j  ,:)*VISC(i+1,j  -st2:j  +st1,k+1)*(dUdxs(ii,2,1, -st2: +st1,1) + dUdxs(2,ii,1, -st2: +st1,1) - 2.0_WP/3.0_WP*dUdxs(2,2,1, -st2: +st1,1)))
                 dSCdx(2,7) = sum(grad_y(i+1,j+1,:)*VISC(i+1,j+1-st2:j+1+st1,k  )*(dUdxs(ii,2,1,1-st2:1+st1,0) + dUdxs(2,ii,1,1-st2:1+st1,0) - 2.0_WP/3.0_WP*dUdxs(2,2,1,1-st2:1+st1,0)))
                 dSCdx(2,8) = sum(grad_y(i+1,j+1,:)*VISC(i+1,j+1-st2:j+1+st1,k+1)*(dUdxs(ii,2,1,1-st2:1+st1,1) + dUdxs(2,ii,1,1-st2:1+st1,1) - 2.0_WP/3.0_WP*dUdxs(2,2,1,1-st2:1+st1,1)))
              else
                 dSCdx(2,1) = sum(grad_y(i  ,j  ,:)*VISC(i  ,j  -st2:j  +st1,k  )*(dUdxs(ii,2,0, -st2: +st1,0) + dUdxs(2,ii,0, -st2: +st1,0)))
                 dSCdx(2,2) = sum(grad_y(i  ,j  ,:)*VISC(i  ,j  -st2:j  +st1,k+1)*(dUdxs(ii,2,0, -st2: +st1,1) + dUdxs(2,ii,0, -st2: +st1,1)))
                 dSCdx(2,3) = sum(grad_y(i  ,j+1,:)*VISC(i  ,j+1-st2:j+1+st1,k  )*(dUdxs(ii,2,0,1-st2:1+st1,0) + dUdxs(2,ii,0,1-st2:1+st1,0)))
                 dSCdx(2,4) = sum(grad_y(i  ,j+1,:)*VISC(i  ,j+1-st2:j+1+st1,k+1)*(dUdxs(ii,2,0,1-st2:1+st1,1) + dUdxs(2,ii,0,1-st2:1+st1,1)))
                 dSCdx(2,5) = sum(grad_y(i+1,j  ,:)*VISC(i+1,j  -st2:j  +st1,k  )*(dUdxs(ii,2,1, -st2: +st1,0) + dUdxs(2,ii,1, -st2: +st1,0)))
                 dSCdx(2,6) = sum(grad_y(i+1,j  ,:)*VISC(i+1,j  -st2:j  +st1,k+1)*(dUdxs(ii,2,1, -st2: +st1,1) + dUdxs(2,ii,1, -st2: +st1,1)))
                 dSCdx(2,7) = sum(grad_y(i+1,j+1,:)*VISC(i+1,j+1-st2:j+1+st1,k  )*(dUdxs(ii,2,1,1-st2:1+st1,0) + dUdxs(2,ii,1,1-st2:1+st1,0)))
                 dSCdx(2,8) = sum(grad_y(i+1,j+1,:)*VISC(i+1,j+1-st2:j+1+st1,k+1)*(dUdxs(ii,2,1,1-st2:1+st1,1) + dUdxs(2,ii,1,1-st2:1+st1,1)))
              end if
              if (ii.eq.3) then
                 dSCdx(3,1) = sum(grad_z(i  ,j  ,:)*VISC(i  ,j  ,k  -st2:k  +st1)*(dUdxs(ii,3,0,0, -st2: +st1) + dUdxs(3,ii,0,0, -st2: +st1) - 2.0_WP/3.0_WP*dUdxs(3,3,0,0, -st2: +st1)))
                 dSCdx(3,2) = sum(grad_z(i  ,j  ,:)*VISC(i  ,j  ,k+1-st2:k+1+st1)*(dUdxs(ii,3,0,0,1-st2:1+st1) + dUdxs(3,ii,0,0,1-st2:1+st1) - 2.0_WP/3.0_WP*dUdxs(3,3,0,0,1-st2:1+st1)))
                 dSCdx(3,3) = sum(grad_z(i  ,j+1,:)*VISC(i  ,j+1,k  -st2:k  +st1)*(dUdxs(ii,3,0,1, -st2: +st1) + dUdxs(3,ii,0,1, -st2: +st1) - 2.0_WP/3.0_WP*dUdxs(3,3,0,1, -st2: +st1)))
                 dSCdx(3,4) = sum(grad_z(i  ,j+1,:)*VISC(i  ,j+1,k+1-st2:k+1+st1)*(dUdxs(ii,3,0,1,1-st2:1+st1) + dUdxs(3,ii,0,1,1-st2:1+st1) - 2.0_WP/3.0_WP*dUdxs(3,3,0,1,1-st2:1+st1)))
                 dSCdx(3,5) = sum(grad_z(i+1,j  ,:)*VISC(i+1,j  ,k  -st2:k  +st1)*(dUdxs(ii,3,1,0, -st2: +st1) + dUdxs(3,ii,1,0, -st2: +st1) - 2.0_WP/3.0_WP*dUdxs(3,3,1,0, -st2: +st1)))
                 dSCdx(3,6) = sum(grad_z(i+1,j  ,:)*VISC(i+1,j  ,k+1-st2:k+1+st1)*(dUdxs(ii,3,1,0,1-st2:1+st1) + dUdxs(3,ii,1,0,1-st2:1+st1) - 2.0_WP/3.0_WP*dUdxs(3,3,1,0,1-st2:1+st1)))
                 dSCdx(3,7) = sum(grad_z(i+1,j+1,:)*VISC(i+1,j+1,k  -st2:k  +st1)*(dUdxs(ii,3,1,1, -st2: +st1) + dUdxs(3,ii,1,1, -st2: +st1) - 2.0_WP/3.0_WP*dUdxs(3,3,1,1, -st2: +st1)))
                 dSCdx(3,8) = sum(grad_z(i+1,j+1,:)*VISC(i+1,j+1,k+1-st2:k+1+st1)*(dUdxs(ii,3,1,1,1-st2:1+st1) + dUdxs(3,ii,1,1,1-st2:1+st1) - 2.0_WP/3.0_WP*dUdxs(3,3,1,1,1-st2:1+st1)))
              else
                 dSCdx(3,1) = sum(grad_z(i  ,j  ,:)*VISC(i  ,j  ,k  -st2:k  +st1)*(dUdxs(ii,3,0,0, -st2: +st1) + dUdxs(3,ii,0,0, -st2: +st1)))
                 dSCdx(3,2) = sum(grad_z(i  ,j  ,:)*VISC(i  ,j  ,k+1-st2:k+1+st1)*(dUdxs(ii,3,0,0,1-st2:1+st1) + dUdxs(3,ii,0,0,1-st2:1+st1)))
                 dSCdx(3,3) = sum(grad_z(i  ,j+1,:)*VISC(i  ,j+1,k  -st2:k  +st1)*(dUdxs(ii,3,0,1, -st2: +st1) + dUdxs(3,ii,0,1, -st2: +st1)))
                 dSCdx(3,4) = sum(grad_z(i  ,j+1,:)*VISC(i  ,j+1,k+1-st2:k+1+st1)*(dUdxs(ii,3,0,1,1-st2:1+st1) + dUdxs(3,ii,0,1,1-st2:1+st1)))
                 dSCdx(3,5) = sum(grad_z(i+1,j  ,:)*VISC(i+1,j  ,k  -st2:k  +st1)*(dUdxs(ii,3,1,0, -st2: +st1) + dUdxs(3,ii,1,0, -st2: +st1)))
                 dSCdx(3,6) = sum(grad_z(i+1,j  ,:)*VISC(i+1,j  ,k+1-st2:k+1+st1)*(dUdxs(ii,3,1,0,1-st2:1+st1) + dUdxs(3,ii,1,0,1-st2:1+st1)))
                 dSCdx(3,7) = sum(grad_z(i+1,j+1,:)*VISC(i+1,j+1,k  -st2:k  +st1)*(dUdxs(ii,3,1,1, -st2: +st1) + dUdxs(3,ii,1,1, -st2: +st1)))
                 dSCdx(3,8) = sum(grad_z(i+1,j+1,:)*VISC(i+1,j+1,k+1-st2:k+1+st1)*(dUdxs(ii,3,1,1,1-st2:1+st1) + dUdxs(3,ii,1,1,1-st2:1+st1)))
              end if
              do jj=1,3
                 if (jj.eq.1) then
                    wx2 = interp_x(n);  wx1 = 1.0_WP-wx2
                    wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                    wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
                 elseif (jj.eq.2) then
                    wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                    wy2 = interp_y(n);  wy1 = 1.0_WP-wy2
                    wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
                 elseif (jj.eq.3) then
                    wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                    wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                    wz2 = interp_z(n);  wz1 = 1.0_WP-wz2
                 end if
                 mval(m +nvars+nscalar+(n-1)*ncols) = &
                      + wx1 * ( wy1 * (wz1*dSCdx(jj,1)+wz2*dSCdx(jj,2)) &
                              + wy2 * (wz1*dSCdx(jj,3)+wz2*dSCdx(jj,4)) )&
                      + wx2 * ( wy1 * (wz1*dSCdx(jj,5)+wz2*dSCdx(jj,6)) &
                              + wy2 * (wz1*dSCdx(jj,7)+wz2*dSCdx(jj,8)) )
                 m = m + 1
              end do
           end do

           ! HR(i,j,k)-> xm(i),ym(j),zm(k) at t=n+1/2
           i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
           j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
           k = index_zm(n); wz2 = interp_zm(n); wz1 = 1.0_WP-wz2
           mval(31+nvars+nscalar+(n-1)*(ncols)) = &
                + wx1 * ( wy1 * (wz1*HR_gp(i  ,j  ,k)+wz2*HR_gp(i  ,j  ,k+1)) &
                        + wy2 * (wz1*HR_gp(i  ,j+1,k)+wz2*HR_gp(i  ,j+1,k+1)) )&
                + wx2 * ( wy1 * (wz1*HR_gp(i+1,j  ,k)+wz2*HR_gp(i+1,j  ,k+1)) &
                        + wy2 * (wz1*HR_gp(i+1,j+1,k)+wz2*HR_gp(i+1,j+1,k+1)) )
        end if
     end if
  end do
  
  ! Get the global value
  do n=1,(ncols)*nprobes
     call parallel_sum(mval(n),tmp)
     mval(n) = tmp
  end do

  ! Transfer to monitor
  nprev = 0
  do ipfile = 1,nprobefiles
     call monitor_select_file(probe_filename(ipfile))
     call monitor_set_array_values( mval(nprev+1 : nprev+(ncols)*lengths(ipfile) ) )
     nprev = nprev + (ncols)*lengths(ipfile)
  end do

  ! z-line (spanwise) probes
  !  -- Only interpolate in x and y
  !  -- Assumes no decomposition in z
  if (probe_line) then
     do n=1,nprobes
        mval_line = 0.0_WP
        if (probe_local(n)) then
           ! Velocity - U
           i = index_x (n); wx2 = interp_x (n); wx1 = 1.0_WP-wx2
           j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
           mval_line(1, kmin:kmax) = &
                + wx1 * ( wy1*U(i  ,j,kmin:kmax) + wy2*U(i  ,j+1,kmin:kmax) )&
                + wx2 * ( wy1*U(i+1,j,kmin:kmax) + wy2*U(i+1,j+1,kmin:kmax) )
           ! Velocity - V
           i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
           j = index_y (n); wy2 = interp_y (n); wy1 = 1.0_WP-wy2
           mval_line(2, kmin:kmax) = &
                + wx1 * ( wy1*V(i  ,j,kmin:kmax) + wy2*V(i  ,j+1,kmin:kmax) )&
                + wx2 * ( wy1*V(i+1,j,kmin:kmax) + wy2*V(i+1,j+1,kmin:kmax) )
           ! Velocity - W
           i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
           j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
           mval_line(3, kmin:kmax) = &
                + wx1 * ( wy1*W(i  ,j,kmin:kmax) + wy2*W(i  ,j+1,kmin:kmax) )&
                + wx2 * ( wy1*W(i+1,j,kmin:kmax) + wy2*W(i+1,j+1,kmin:kmax) )
           ! Pressure
           i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
           j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
           mval_line(4, kmin:kmax) = &
                + wx1 * ( wy1*P(i  ,j,kmin:kmax) + wy2*P(i  ,j+1,kmin:kmax) )&
                + wx2 * ( wy1*P(i+1,j,kmin:kmax) + wy2*P(i+1,j+1,kmin:kmax) )
           ! Scalars
           do isc=1,nscalar
              mval_line(4+isc, kmin:kmax) = &
                + wx1 * ( wy1*SC(i  ,j,kmin:kmax,isc) + wy2*SC(i  ,j+1,kmin:kmax,isc) )&
                + wx2 * ( wy1*SC(i+1,j,kmin:kmax,isc) + wy2*SC(i+1,j+1,kmin:kmax,isc) )
           end do
           ! Other variables
           do m=1,nvars
              mval_line(4+nscalar+m, kmin:kmax) = &
                + wx1 * ( wy1*varsdata(i  ,j,kmin:kmax,m) + wy2*varsdata(i  ,j+1,kmin:kmax,m) )&
                + wx2 * ( wy1*varsdata(i+1,j,kmin:kmax,m) + wy2*varsdata(i+1,j+1,kmin:kmax,m) )
           end do

           ! Information for spectra
           if (probe_spectra) then
              ! Density
              mval_line(5 +nvars+nscalar, kmin:kmax) = &
                   + wx1 * ( wy1*RHO(i  ,j,kmin:kmax) + wy2*RHO(i  ,j+1,kmin:kmax) )&
                   + wx2 * ( wy1*RHO(i+1,j,kmin:kmax) + wy2*RHO(i+1,j+1,kmin:kmax) )

              ! Viscosity
              mval_line(6 +nvars+nscalar, kmin:kmax) = &
                   + wx1 * ( wy1*VISC(i  ,j,kmin:kmax) + wy2*VISC(i  ,j+1,kmin:kmax) )&
                   + wx2 * ( wy1*VISC(i+1,j,kmin:kmax) + wy2*VISC(i+1,j+1,kmin:kmax) )

              ! Velocity gradients
              ! dUdx computed at cell centers using i,j,k indices at cell faces
              i = index_x(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
              j = index_y(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
              do k=kmin,kmax
                 call vel_grad_local(i  ,j  ,k, dUdx(:,:,1))
                 call vel_grad_local(i  ,j+1,k, dUdx(:,:,2))
                 call vel_grad_local(i+1,j  ,k, dUdx(:,:,3))
                 call vel_grad_local(i+1,j+1,k, dUdx(:,:,4))
                 m = 7
                 do ii=1,3
                    do jj=1,3
                       mval_line(m +nvars+nscalar, k) = &
                            + wx1 * ( wy1*dUdx(ii,jj,1) + wy2*dUdx(ii,jj,2) )&
                            + wx2 * ( wy1*dUdx(ii,jj,3) + wy2*dUdx(ii,jj,4) )
                       m = m + 1
                    end do
                 end do
              end do              
              
              ! Density gradient
              ! dRHOdx computed at cell faces using i,j,k indices at cell centers
              i = index_xm(n)
              j = index_ym(n)
              do k=kmin,kmax
                 dSCdx(1,1) = sum(grad_x(i  ,j  ,:)*RHO(i  -st2:i  +st1,j  ,k  ))
                 dSCdx(1,2) = sum(grad_x(i  ,j+1,:)*RHO(i  -st2:i  +st1,j+1,k  ))
                 dSCdx(1,3) = sum(grad_x(i+1,j  ,:)*RHO(i+1-st2:i+1+st1,j  ,k  ))
                 dSCdx(1,4) = sum(grad_x(i+1,j+1,:)*RHO(i+1-st2:i+1+st1,j+1,k  ))
                 dSCdx(2,1) = sum(grad_y(i  ,j  ,:)*RHO(i  ,j  -st2:j  +st1,k  ))
                 dSCdx(2,2) = sum(grad_y(i  ,j+1,:)*RHO(i  ,j+1-st2:j+1+st1,k  ))
                 dSCdx(2,3) = sum(grad_y(i+1,j  ,:)*RHO(i+1,j  -st2:j  +st1,k  ))
                 dSCdx(2,4) = sum(grad_y(i+1,j+1,:)*RHO(i+1,j+1-st2:j+1+st1,k  ))
                 dSCdx(3,1) = sum(grad_z(i  ,j  ,:)*RHO(i  ,j  ,k  -st2:k  +st1))
                 dSCdx(3,2) = sum(grad_z(i  ,j+1,:)*RHO(i  ,j+1,k  -st2:k  +st1))
                 dSCdx(3,3) = sum(grad_z(i+1,j  ,:)*RHO(i+1,j  ,k  -st2:k  +st1))
                 dSCdx(3,4) = sum(grad_z(i+1,j+1,:)*RHO(i+1,j+1,k  -st2:k  +st1))
                 m = 16
                 do ii=1,3
                    if (ii.eq.1) then
                       wx2 = interp_x(n);  wx1 = 1.0_WP-wx2
                       wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                    elseif (ii.eq.2) then
                       wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                       wy2 = interp_y(n);  wy1 = 1.0_WP-wy2
                    elseif (ii.eq.3) then
                       wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                       wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                    end if
                    mval_line(m +nvars+nscalar, k) = &
                         + wx1 * ( wy1*dSCdx(ii,1) + wy2*dSCdx(ii,2) )&
                         + wx2 * ( wy1*dSCdx(ii,3) + wy2*dSCdx(ii,4) )
                    m = m + 1
                 end do
              end do

              ! Pressure gradient
              ! dPdx computed at cell faces using i,j,k indices at cell centers
              i = index_xm(n)
              j = index_ym(n)
              do k=kmin,kmax
                 dSCdx(1,1) = sum(grad_x(i  ,j  ,:)*P(i  -st2:i  +st1,j  ,k  ))
                 dSCdx(1,2) = sum(grad_x(i  ,j+1,:)*P(i  -st2:i  +st1,j+1,k  ))
                 dSCdx(1,3) = sum(grad_x(i+1,j  ,:)*P(i+1-st2:i+1+st1,j  ,k  ))
                 dSCdx(1,4) = sum(grad_x(i+1,j+1,:)*P(i+1-st2:i+1+st1,j+1,k  ))
                 dSCdx(2,1) = sum(grad_y(i  ,j  ,:)*P(i  ,j  -st2:j  +st1,k  ))
                 dSCdx(2,2) = sum(grad_y(i  ,j+1,:)*P(i  ,j+1-st2:j+1+st1,k  ))
                 dSCdx(2,3) = sum(grad_y(i+1,j  ,:)*P(i+1,j  -st2:j  +st1,k  ))
                 dSCdx(2,4) = sum(grad_y(i+1,j+1,:)*P(i+1,j+1-st2:j+1+st1,k  ))
                 dSCdx(3,1) = sum(grad_z(i  ,j  ,:)*P(i  ,j  ,k  -st2:k  +st1))
                 dSCdx(3,2) = sum(grad_z(i  ,j+1,:)*P(i  ,j+1,k  -st2:k  +st1))
                 dSCdx(3,3) = sum(grad_z(i+1,j  ,:)*P(i+1,j  ,k  -st2:k  +st1))
                 dSCdx(3,4) = sum(grad_z(i+1,j+1,:)*P(i+1,j+1,k  -st2:k  +st1))
                 m = 19
                 do ii=1,3
                    if (ii.eq.1) then
                       wx2 = interp_x(n);  wx1 = 1.0_WP-wx2
                       wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                    elseif (ii.eq.2) then
                       wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                       wy2 = interp_y(n);  wy1 = 1.0_WP-wy2
                    elseif (ii.eq.3) then
                       wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                       wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                    end if
                    mval_line(m +nvars+nscalar, k) = &
                         + wx1 * ( wy1*dSCdx(ii,1) + wy2*dSCdx(ii,2) )&
                         + wx2 * ( wy1*dSCdx(ii,3) + wy2*dSCdx(ii,4) )
                    m = m + 1
                 end do
              end do

              ! Stress tensor gradients
              ! dUdx computed at cell centers using i,j,k indices at cell faces
              do k=kmin,kmax
                 i = index_x(n)
                 j = index_y(n)
                 do ii=-st2,1+st1
                    do jj=-st2,1+st1
                       do kk=-st2,st1
                          call vel_grad_local(i+ii,j+jj,k+kk, dUdxs(:,:,ii,jj,kk))
                       end do
                    end do
                 end do
                 i = index_xm(n)
                 j = index_ym(n)
                 m = 22
                 do ii=1,3
                    if (ii.eq.1) then
                       dSCdx(1,1) = sum(grad_x(i  ,j  ,:)*VISC(i  -st2:i  +st1,j  ,k  )*(dUdxs(ii,1, -st2: +st1,0,0) + dUdxs(1,ii, -st2: +st1,0,0) - 2.0_WP/3.0_WP*dUdxs(1,1, -st2: +st1,0,0)))
                       dSCdx(1,2) = sum(grad_x(i  ,j+1,:)*VISC(i  -st2:i  +st1,j+1,k  )*(dUdxs(ii,1, -st2: +st1,1,0) + dUdxs(1,ii, -st2: +st1,1,0) - 2.0_WP/3.0_WP*dUdxs(1,1, -st2: +st1,1,0)))
                       dSCdx(1,3) = sum(grad_x(i+1,j  ,:)*VISC(i+1-st2:i+1+st1,j  ,k  )*(dUdxs(ii,1,1-st2:1+st1,0,0) + dUdxs(1,ii,1-st2:1+st1,0,0) - 2.0_WP/3.0_WP*dUdxs(1,1,1-st2:1+st1,0,0)))
                       dSCdx(1,4) = sum(grad_x(i+1,j+1,:)*VISC(i+1-st2:i+1+st1,j+1,k  )*(dUdxs(ii,1,1-st2:1+st1,1,0) + dUdxs(1,ii,1-st2:1+st1,1,0) - 2.0_WP/3.0_WP*dUdxs(1,1,1-st2:1+st1,1,0)))
                    else
                       dSCdx(1,1) = sum(grad_x(i  ,j  ,:)*VISC(i  -st2:i  +st1,j  ,k  )*(dUdxs(ii,1, -st2: +st1,0,0) + dUdxs(1,ii, -st2: +st1,0,0)))
                       dSCdx(1,2) = sum(grad_x(i  ,j+1,:)*VISC(i  -st2:i  +st1,j+1,k  )*(dUdxs(ii,1, -st2: +st1,1,0) + dUdxs(1,ii, -st2: +st1,1,0)))
                       dSCdx(1,3) = sum(grad_x(i+1,j  ,:)*VISC(i+1-st2:i+1+st1,j  ,k  )*(dUdxs(ii,1,1-st2:1+st1,0,0) + dUdxs(1,ii,1-st2:1+st1,0,0)))
                       dSCdx(1,4) = sum(grad_x(i+1,j+1,:)*VISC(i+1-st2:i+1+st1,j+1,k  )*(dUdxs(ii,1,1-st2:1+st1,1,0) + dUdxs(1,ii,1-st2:1+st1,1,0)))
                    end if
                    if (ii.eq.2) then
                       dSCdx(2,1) = sum(grad_y(i  ,j  ,:)*VISC(i  ,j  -st2:j  +st1,k  )*(dUdxs(ii,2,0, -st2: +st1,0) + dUdxs(2,ii,0, -st2: +st1,0) - 2.0_WP/3.0_WP*dUdxs(2,2,0, -st2: +st1,0)))
                       dSCdx(2,2) = sum(grad_y(i  ,j+1,:)*VISC(i  ,j+1-st2:j+1+st1,k  )*(dUdxs(ii,2,0,1-st2:1+st1,0) + dUdxs(2,ii,0,1-st2:1+st1,0) - 2.0_WP/3.0_WP*dUdxs(2,2,0,1-st2:1+st1,0)))
                       dSCdx(2,3) = sum(grad_y(i+1,j  ,:)*VISC(i+1,j  -st2:j  +st1,k  )*(dUdxs(ii,2,1, -st2: +st1,0) + dUdxs(2,ii,1, -st2: +st1,0) - 2.0_WP/3.0_WP*dUdxs(2,2,1, -st2: +st1,0)))
                       dSCdx(2,4) = sum(grad_y(i+1,j+1,:)*VISC(i+1,j+1-st2:j+1+st1,k  )*(dUdxs(ii,2,1,1-st2:1+st1,0) + dUdxs(2,ii,1,1-st2:1+st1,0) - 2.0_WP/3.0_WP*dUdxs(2,2,1,1-st2:1+st1,0)))
                    else
                       dSCdx(2,1) = sum(grad_y(i  ,j  ,:)*VISC(i  ,j  -st2:j  +st1,k  )*(dUdxs(ii,2,0, -st2: +st1,0) + dUdxs(2,ii,0, -st2: +st1,0)))
                       dSCdx(2,2) = sum(grad_y(i  ,j+1,:)*VISC(i  ,j+1-st2:j+1+st1,k  )*(dUdxs(ii,2,0,1-st2:1+st1,0) + dUdxs(2,ii,0,1-st2:1+st1,0)))
                       dSCdx(2,3) = sum(grad_y(i+1,j  ,:)*VISC(i+1,j  -st2:j  +st1,k  )*(dUdxs(ii,2,1, -st2: +st1,0) + dUdxs(2,ii,1, -st2: +st1,0)))
                       dSCdx(2,4) = sum(grad_y(i+1,j+1,:)*VISC(i+1,j+1-st2:j+1+st1,k  )*(dUdxs(ii,2,1,1-st2:1+st1,0) + dUdxs(2,ii,1,1-st2:1+st1,0)))
                    end if
                    if (ii.eq.3) then
                       dSCdx(3,1) = sum(grad_z(i  ,j  ,:)*VISC(i  ,j  ,k  -st2:k  +st1)*(dUdxs(ii,3,0,0, -st2: +st1) + dUdxs(3,ii,0,0, -st2: +st1) - 2.0_WP/3.0_WP*dUdxs(3,3,0,0, -st2: +st1)))
                       dSCdx(3,2) = sum(grad_z(i  ,j+1,:)*VISC(i  ,j+1,k  -st2:k  +st1)*(dUdxs(ii,3,0,1, -st2: +st1) + dUdxs(3,ii,0,1, -st2: +st1) - 2.0_WP/3.0_WP*dUdxs(3,3,0,1, -st2: +st1)))
                       dSCdx(3,3) = sum(grad_z(i+1,j  ,:)*VISC(i+1,j  ,k  -st2:k  +st1)*(dUdxs(ii,3,1,0, -st2: +st1) + dUdxs(3,ii,1,0, -st2: +st1) - 2.0_WP/3.0_WP*dUdxs(3,3,1,0, -st2: +st1)))
                       dSCdx(3,4) = sum(grad_z(i+1,j+1,:)*VISC(i+1,j+1,k  -st2:k  +st1)*(dUdxs(ii,3,1,1, -st2: +st1) + dUdxs(3,ii,1,1, -st2: +st1) - 2.0_WP/3.0_WP*dUdxs(3,3,1,1, -st2: +st1)))
                    else
                       dSCdx(3,1) = sum(grad_z(i  ,j  ,:)*VISC(i  ,j  ,k  -st2:k  +st1)*(dUdxs(ii,3,0,0, -st2: +st1) + dUdxs(3,ii,0,0, -st2: +st1)))
                       dSCdx(3,2) = sum(grad_z(i  ,j+1,:)*VISC(i  ,j+1,k  -st2:k  +st1)*(dUdxs(ii,3,0,1, -st2: +st1) + dUdxs(3,ii,0,1, -st2: +st1)))
                       dSCdx(3,3) = sum(grad_z(i+1,j  ,:)*VISC(i+1,j  ,k  -st2:k  +st1)*(dUdxs(ii,3,1,0, -st2: +st1) + dUdxs(3,ii,1,0, -st2: +st1)))
                       dSCdx(3,4) = sum(grad_z(i+1,j+1,:)*VISC(i+1,j+1,k  -st2:k  +st1)*(dUdxs(ii,3,1,1, -st2: +st1) + dUdxs(3,ii,1,1, -st2: +st1)))
                    end if
                    do jj=1,3
                       if (jj.eq.1) then
                          wx2 = interp_x(n);  wx1 = 1.0_WP-wx2
                          wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                       elseif (jj.eq.2) then
                          wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                          wy2 = interp_y(n);  wy1 = 1.0_WP-wy2
                       elseif (jj.eq.3) then
                          wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
                          wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
                       end if
                       mval_line(m +nvars+nscalar, k) = &
                                + wx1 * ( wy1*dSCdx(jj,1) + wy2*dSCdx(jj,2) )&
                                + wx2 * ( wy1*dSCdx(jj,3) + wy2*dSCdx(jj,4) )
                       m = m + 1
                    end do
                 end do
              end do

              ! Heat release rate
              i = index_xm(n); wx2 = interp_xm(n); wx1 = 1.0_WP-wx2
              j = index_ym(n); wy2 = interp_ym(n); wy1 = 1.0_WP-wy2
              mval_line(31+nvars+nscalar, kmin:kmax) = &
                   + wx1 * ( wy1*HR_gp(i  ,j,kmin:kmax) + wy2*HR_gp(i  ,j+1,kmin:kmax) )&
                   + wx2 * ( wy1*HR_gp(i+1,j,kmin:kmax) + wy2*HR_gp(i+1,j+1,kmin:kmax) )
           end if
        end if

        ! Transfer this probe's data to monitor
        do m=1,ncols
           ! Get the global value
           do k=kmin,kmax
              call parallel_sum(mval_line(m,k), tmp)
              mval_line(m,k) = tmp
           end do
           call monitor_select_file( probe_line_filenames(indx_pfolder(n), m+(indx_pfile(n)-1)*ncols) )
           call monitor_set_array_values( mval_line(m, :) )
        end do
     end do
  end if

  return
end subroutine probe_monitor

