module dump_ensight_str
  use parallel
  use precision
  use geometry
  use partition
  use masks
  use fileio
  implicit none
  
  ! Time info
  integer :: nout_time
  real(WP), dimension(:), allocatable :: out_time
  
  ! SP buffers
  real(SP), dimension(:,:,:), pointer :: buffer1_SP
  real(SP), dimension(:,:,:,:), pointer :: buffer3_SP
  
  ! Fileviews
  integer :: fileview,datasize,gdatasize
  
  ! Monitoring
  integer  :: nfiles
  real(WP) :: time_open,time_close,time_write
  
end module dump_ensight_str


! ================================================= !
! Dump 3D binary ensight gold data - Initialization !
! ================================================= !
subroutine dump_ensight_str_3D_init
  use dump_ensight_str
  use parser
  use parallel
  use time_info
  implicit none
  integer :: ierr,i
  logical :: file_is_there
  integer, dimension(3) :: gsizes,lsizes,start
  
  ! Create & Start the timer
  call timing_create('ensight-str')
  call timing_start ('ensight-str')

  ! Allocate buffers
  allocate(buffer1_SP(imin_:imax_,jmin_:jmax_,kmin_:kmax_))
  allocate(buffer3_SP(imin_:imax_,jmin_:jmax_,kmin_:kmax_,3))
  
  ! Open the file
  inquire(file='ensight-3D/arts.case',exist=file_is_there)
  if (file_is_there) then
     ! Read the file
     call parser_parsefile('ensight-3D/arts.case')
     ! Get the time
     call parser_getsize('time values',nout_time)
     allocate(out_time(nout_time))
     call parser_read('time values',out_time)
     ! Remove future time
     future: do i=1,size(out_time)
        if (out_time(i).GE.time*0.99999_WP) then
           nout_time = i-1
           exit future
        end if
     end do future
  else
     ! Create directory
     !if (irank.eq.iroot) call CREATE_FOLDER("ensight-3D")
     if (irank.eq.iroot) call system("mkdir -p ensight-3D")
     ! Set the time
     nout_time = 0
     allocate(out_time(1))
  end if
  
  ! Write the geometry
  if (irank.eq.iroot) call dump_ensight_str_3D_geometry
  
  ! Create the view
  gsizes(1) = nx
  gsizes(2) = ny
  gsizes(3) = nz
  lsizes(1) = nx_
  lsizes(2) = ny_
  lsizes(3) = nz_
  start(1) = imin_-imin
  start(2) = jmin_-jmin
  start(3) = kmin_-kmin
  datasize = lsizes(1)*lsizes(2)*lsizes(3)
  gdatasize= gsizes(1)*gsizes(2)*gsizes(3)
  call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,start,MPI_ORDER_FORTRAN,MPI_REAL_SP,fileview,ierr)
  call MPI_TYPE_COMMIT(fileview,ierr)
  
  ! Stop the timer
  call timing_stop ('ensight-str')
  
  ! Create nw file to monitor
  call monitor_create_file_step('ensight-str',4)
  call monitor_set_header(1, 'nfiles','i')
  call monitor_set_header(2, 'open [s/file]','r')
  call monitor_set_header(3, 'write [s/file]','r')
  call monitor_set_header(4, 'close [s/file]','r')

  ! Save the first field
  call dump_ensight_str_3D
  
  return
end subroutine dump_ensight_str_3D_init


! ================================ !
! Dump 3D binary ensight gold data !
! ================================ !
subroutine dump_ensight_str_3D
  use dump_ensight_str
  use data
  use partition
  use string
  use interpolate
  use memory
  use combustion
  use sgsmodel
  use velocity
  use soot
  use strainrate
  use pressure
  use strainrate
  implicit none
  character(len=str_short) :: name
  integer :: i,j,k,isc
  
  ! Reset timing
  nfiles = 0
  time_open  = 0.0_WP
  time_write = 0.0_WP
  time_close = 0.0_WP
  
  ! Start a timer
  call timing_start('ensight-str')

  ! Update the case file
  call dump_ensight_str_3D_case
  
  ! Pressure
  name = "P"
  buffer1_SP = P(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
  call dump_ensight_str_3D_scalar(buffer1_SP,name)
     
  ! Viscosity
  name = "VISC"
  buffer1_SP = VISC(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
  call dump_ensight_str_3D_scalar(buffer1_SP,name)
  if (pope_crit) then
     name = "KSGS"
     buffer1_SP = k_sgs(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
     call dump_ensight_str_3D_scalar(buffer1_SP,name)
  end if
  
  ! Velocity
  name = "V"
  if (icyl.eq.0) then
     buffer3_SP(:,:,:,1) = Ui(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
     buffer3_SP(:,:,:,2) = Vi(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
     buffer3_SP(:,:,:,3) = Wi(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
     call dump_ensight_str_3D_vector(buffer3_SP,name)
  else
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              buffer3_SP(i,j,k,1) = Ui(i,j,k)
              buffer3_SP(i,j,k,2) = Vi(i,j,k) * cos(zm(k)) - Wi(i,j,k) * sin(zm(k))
              buffer3_SP(i,j,k,3) = Vi(i,j,k) * sin(zm(k)) + Wi(i,j,k) * cos(zm(k))
           end do
        end do
     end do
  end if
  call dump_ensight_str_3D_vector(buffer3_SP,name)
  
  ! Scalars
  do isc=1,nscalar
     name = SC_name(isc)
     buffer1_SP = SC(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc)
     call dump_ensight_str_3D_scalar(buffer1_SP,name)
  end do
  if (nscalar.ge.1) then
     name = 'DIFF'
     buffer1_SP = DIFF(imin_:imax_,jmin_:jmax_,kmin_:kmax_,isc_ZMIX)
     call dump_ensight_str_3D_scalar(buffer1_SP,name)
  end if
  
  ! Combustion
  if (combust) then
     if (isc_T.eq.0) then
        name = 'Temp'
        buffer1_SP = T(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
        call dump_ensight_str_3D_scalar(buffer1_SP,name)
     end if
     name = 'RHO'
     buffer1_SP = RHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
     call dump_ensight_str_3D_scalar(buffer1_SP,name)
     name = 'dRHO'
     buffer1_SP = dRHO(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
     call dump_ensight_str_3D_scalar(buffer1_SP,name)
     if (chemistry.eq.'finite chem') then
        name = "HR"
        buffer1_SP = HR_gp(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
        call dump_ensight_str_3D_scalar(buffer1_SP,name)
     end if
     if (isc_ZVAR.eq.0 .and. isc_ZMIX2.eq.0 .and. use_ZVAR) then
        name = 'ZVAR'
        buffer1_SP = ZVAR(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
        call dump_ensight_str_3D_scalar(buffer1_SP,name)
     end if
     if (isc_ZMIX.ne.0) then
        ! Output chi as CHI = DIFF*|nabla(Z)|^2
        name = 'CHI'
        buffer1_SP = CHI(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
        call dump_ensight_str_3D_scalar(buffer1_SP,name)
     end if
     if (trim(chemistry).eq.'chemtable') then
        call die("bugs with precision: do later...")
!!$        name = 'Y_F'
!!$        call chemtable_lookup('Y_F',buffer1_SP)
!!$        call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$        name = 'Y_O2'
!!$        call chemtable_lookup('Y_O2',buffer1_SP)
!!$        call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$        name = 'Y_CO'
!!$        call chemtable_lookup('Y_CO',buffer1_SP)
!!$        call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$        name = 'Y_CO2'
!!$        call chemtable_lookup('Y_CO2',buffer1_SP)
!!$        call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$        name = 'Y_H2'
!!$        call chemtable_lookup('Y_H2',buffer1_SP)
!!$        call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$        name = 'Y_H2O'
!!$        call chemtable_lookup('Y_H2O',buffer1_SP)
!!$        call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$        name = 'Y_OH'
!!$        call chemtable_lookup('Y_OH',buffer1_SP)
!!$        call dump_ensight_str_3D_scalar(buffer1_SP,name)
     end if
  end if
  
!!$  ! Soot
!!$  if (use_soot) then
!!$     name = "numdens"
!!$     buffer1_SP = numdens(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
!!$     call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$     name = "fv"
!!$     buffer1_SP = fv(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
!!$     call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$     name = "partdiam"
!!$     buffer1_SP = partdiam(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
!!$     call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$     name = "partaggr"
!!$     buffer1_SP = partaggr(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
!!$     call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$     name = "intermit"
!!$     buffer1_SP = intermit(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
!!$     call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$     name = "Y_PAH"
!!$     buffer1_SP = Y_PAH(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
!!$     call dump_ensight_str_3D_scalar(buffer1_SP,name)
!!$  end if

  ! Divergence
  name = "divg"
  buffer1_SP = divg
  call dump_ensight_str_3D_scalar(buffer1_SP,name)
  
  ! PQR criterion
  call strainrate_PQRcrit
  name = 'Pcrit'
  buffer1_SP = Pcrit(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
  call dump_ensight_str_3D_scalar(buffer1_SP,name)
  name = 'Qcrit'
  buffer1_SP = Qcrit(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
  call dump_ensight_str_3D_scalar(buffer1_SP,name)
  name = 'Rcrit'
  buffer1_SP = Rcrit(imin_:imax_,jmin_:jmax_,kmin_:kmax_)
  call dump_ensight_str_3D_scalar(buffer1_SP,name)
  
  ! Stop a timer
  call timing_stop('ensight-str')

  ! Transfer the values to monitor
  call monitor_select_file('ensight-str')
  call monitor_set_single_value(1,real(nfiles,WP))
  call monitor_set_single_value(2,time_open /real(nfiles,WP))
  call monitor_set_single_value(3,time_write/real(nfiles,WP))
  call monitor_set_single_value(4,time_close/real(nfiles,WP))

  return
end subroutine dump_ensight_str_3D


! =========================================== !
! Dump 3D binary ensight gold data - geometry !
! => one processor only - test before         !
! =========================================== !
subroutine dump_ensight_str_3D_geometry
  use dump_ensight_str
  implicit none
  integer :: ierr,ipart,iunit,i,j,k,ii,jj,kk
  integer, dimension(:,:,:), pointer :: iblank
  character(len=80) :: buffer
  real(SP), dimension(:,:,:), pointer :: xbuf,ybuf,zbuf
  real(SP) :: max_x,max_y,max_z
  real(SP) :: min_x,min_y,min_z
  
  ! Get single precision mesh
  allocate(xbuf(1:nx+1,1:ny+1,1:nz+1),ybuf(1:nx+1,1:ny+1,1:nz+1),zbuf(1:nx+1,1:ny+1,1:nz+1))
  if (icyl.eq.0) then
     do i=imin,imax+1
        do j=jmin,jmax+1
           do k=kmin,kmax+1
              ii=i-imin+1
              jj=j-jmin+1
              kk=k-kmin+1
              xbuf(ii,jj,kk)=x(i)
              ybuf(ii,jj,kk)=y(j)
              zbuf(ii,jj,kk)=z(k)
           end do
        end do
     end do
  else
     do i=imin,imax+1
        do j=jmin,jmax+1
           do k=kmin,kmax+1
              ii=i-imin+1
              jj=j-jmin+1
              kk=k-kmin+1
              xbuf(ii,jj,kk)=x(i)
              ybuf(ii,jj,kk)=y(j)*cos(z(k))
              zbuf(ii,jj,kk)=y(j)*sin(z(k))
           end do
        end do
     end do
  end if
  max_x = maxval(xbuf)
  max_y = maxval(ybuf)
  max_z = maxval(zbuf)
  min_x = minval(xbuf)
  min_y = minval(ybuf)
  min_z = minval(zbuf)
  
  ! Get the iblank in 3D
  allocate(iblank(nx+1,ny+1,nz+1))
  do k=1,nz+1
     iblank(1:nx+1,1:ny+1,k) = 1-mask_node(imin:imax+1,jmin:jmax+1)
  end do
  
  ! Open the file
  call BINARY_FILE_OPEN(iunit,"ensight-3D/geometry","w",ierr)
  
  ! Write the geometry
  buffer = 'C Binary'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  buffer = 'Ensight Gold Geometry File'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  buffer = 'Structured Geometry from ARTS'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  buffer = 'node id off'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  buffer = 'element id off'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  
  buffer = 'extents'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  call BINARY_FILE_WRITE(iunit,min_x,1,kind(min_x),ierr)
  call BINARY_FILE_WRITE(iunit,max_x,1,kind(max_x),ierr)
  call BINARY_FILE_WRITE(iunit,min_y,1,kind(min_y),ierr)
  call BINARY_FILE_WRITE(iunit,max_y,1,kind(max_y),ierr)
  call BINARY_FILE_WRITE(iunit,min_z,1,kind(min_z),ierr)
  call BINARY_FILE_WRITE(iunit,max_z,1,kind(max_z),ierr)
  
  buffer = 'part'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  ipart = 1
  call BINARY_FILE_WRITE(iunit,ipart,1,kind(ipart),ierr)
  
  buffer = 'Complete geometry'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  buffer = 'block curvilinear iblanked'
  call BINARY_FILE_WRITE(iunit,buffer,80,kind(buffer),ierr)
  
  call BINARY_FILE_WRITE(iunit,nx+1,1,kind(nx),ierr)
  call BINARY_FILE_WRITE(iunit,ny+1,1,kind(ny),ierr)
  call BINARY_FILE_WRITE(iunit,nz+1,1,kind(nz),ierr)
  
  call BINARY_FILE_WRITE(iunit,xbuf,(nx+1)*(ny+1)*(nz+1),kind(xbuf),ierr)
  call BINARY_FILE_WRITE(iunit,ybuf,(nx+1)*(ny+1)*(nz+1),kind(ybuf),ierr)
  call BINARY_FILE_WRITE(iunit,zbuf,(nx+1)*(ny+1)*(nz+1),kind(zbuf),ierr)
  call BINARY_FILE_WRITE(iunit,iblank,(nx+1)*(ny+1)*(nz+1),kind(iblank),ierr)
  
  ! Close the file
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  return
end subroutine dump_ensight_str_3D_geometry


! ============================================ !
! Dump 3D binary ensight gold data - case file !
! ============================================ !
subroutine dump_ensight_str_3D_case
  use dump_ensight_str
  use combustion
  use sgsmodel
  use time_info
  use string
  implicit none
  integer :: iunit,ierr,isc
  real(WP), dimension(:), allocatable :: buffer
  character(len=80) :: str
  
  ! Update the time info
  allocate(buffer(nout_time))
  buffer(1:nout_time) = out_time(1:nout_time)
  deallocate(out_time)
  nout_time = nout_time + 1
  allocate(out_time(nout_time))
  out_time(1:nout_time-1) = buffer(1:nout_time-1)
  out_time(nout_time) = time
  
  ! Write - Single proc & parallel => only root writes (in ASCII)
  if (irank==iroot) then
     ! Open the file
     iunit = iopen()
     open(iunit,file="ensight-3D/arts.case",form="formatted",iostat=ierr,status="REPLACE")
     ! Write the case
     str='FORMAT'
     write(iunit,'(a80)') str
     str='type: ensight gold'
     write(iunit,'(a80)') str
     str='GEOMETRY'
     write(iunit,'(a80)') str
     str='model: geometry'
     write(iunit,'(a80)') str
     str='VARIABLE'
     write(iunit,'(a80)') str
     ! Pressure
     str='scalar per element: 1 Pressure P/P.******'
     write(iunit,'(a80)') str
     ! Viscosity
     str='scalar per element: 1 Viscosity VISC/VISC.******'
     write(iunit,'(a80)') str
     if (pope_crit) then
        str='scalar per element: 1 KSGS KSGS/KSGS.******'
        write(iunit,'(a80)') str
     end if
     ! Diffusivity
     if (nscalar.ge.1) then
        str='scalar per element: 1 Diffusivity DIFF/DIFF.******'
        write(iunit,'(a80)') str
     end if
     ! Scalars
     do isc=1,nscalar
        str = 'scalar per element: 1 '//trim(SC_name(isc))//' '//trim(SC_name(isc))//'/'//trim(SC_name(isc))//'.******'
        write(iunit,'(a80)') str
     end do
     ! Chemistry
     if (combust) then
        if (isc_T.eq.0) then
           str = 'scalar per element: 1 Temp Temp/Temp.******'
           write(iunit,'(a80)') str
        end if
        str = 'scalar per element: 1 RHO RHO/RHO.******'
        write(iunit,'(a80)') str
        str = 'scalar per element: 1 dRHO dRHO/dRHO.******'
        write(iunit,'(a80)') str
        if (chemistry.eq.'finite chem') then
           str='scalar per element: 1 HR HR/HR.******'
           write(iunit,'(a80)') str
        end if
        str = 'scalar per element: 1 ZVAR ZVAR/ZVAR.******'
        if (isc_ZVAR.eq.0 .and. isc_ZMIX2.eq.0 .and. use_ZVAR) write(iunit,'(a80)') str
        str = 'scalar per element: 1 CHI CHI/CHI.******'
        if (isc_ZMIX.ne.0) write(iunit,'(a80)') str
        if (trim(chemistry).eq.'chemtable') then
           str = 'scalar per element: 1 Y_F Y_F/Y_F.******'
           write(iunit,'(a80)') str
           str = 'scalar per element: 1 Y_O2 Y_O2/Y_O2.******'
           write(iunit,'(a80)') str
           str = 'scalar per element: 1 Y_CO Y_CO/Y_CO.******'
           write(iunit,'(a80)') str
           str = 'scalar per element: 1 Y_CO2 Y_CO2/Y_CO2.******'
           write(iunit,'(a80)') str
           str = 'scalar per element: 1 Y_H2 Y_H2/Y_H2.******'
           write(iunit,'(a80)') str
           str = 'scalar per element: 1 Y_H2O Y_H2O/Y_H2O.******'
           write(iunit,'(a80)') str
           str = 'scalar per element: 1 Y_OH Y_OH/Y_OH.******'
           write(iunit,'(a80)') str
        end if
     end if
     ! Velocity
     str='vector per element: 1 V V/V.******'
     write(iunit,'(a80)') str
!!$     ! Soot stuff
!!$     if (use_soot) then
!!$        str='scalar per element: 1 numdens numdens/numdens.******'
!!$        write(iunit,'(a80)') str
!!$        str='scalar per element: 1 fv fv/fv.******'
!!$        write(iunit,'(a80)') str
!!$        str='scalar per element: 1 partdiam partdiam/partdiam.******'
!!$        write(iunit,'(a80)') str
!!$        str='scalar per element: 1 partaggr partaggr/partaggr.******'
!!$        write(iunit,'(a80)') str
!!$        str='scalar per element: 1 intermit intermit/intermit.******'
!!$        write(iunit,'(a80)') str
!!$        str='scalar per element: 1 Y_PAH Y_PAH/Y_PAH.******'
!!$        write(iunit,'(a80)') str
!!$     end if
     ! Divergence
     str='scalar per element: 1 divg divg/divg.******'
     write(iunit,'(a80)') str
     ! PQR criterion
     str='scalar per element: 1 Pcrit Pcrit/Pcrit.******'
     write(iunit,'(a80)') str
     str='scalar per element: 1 Qcrit Qcrit/Qcrit.******'
     write(iunit,'(a80)') str
     str='scalar per element: 1 Rcrit Rcrit/Rcrit.******'
     write(iunit,'(a80)') str
     ! Time section
     str='TIME'
     write(iunit,'(a80)') str
     str='time set: 1'
     write(iunit,'(a80)') str
     str='number of steps:'
     write(iunit,'(a16,x,i12)') str,nout_time
     str='filename start number: 1'
     write(iunit,'(a80)') str
     str='filename increment: 1'
     write(iunit,'(a80)') str
     str='time values:'
     write(iunit,'(a12,x,10000000(3(ES12.5,x),/))') str,out_time
     ! Close the file
     close(iclose(iunit))
  end if
  
  return
end subroutine dump_ensight_str_3D_case


! ========================================= !
! Dump 3D binary ensight gold data - scalar !
! ========================================= !
subroutine dump_ensight_str_3D_scalar(scalar,name)
  use dump_ensight_str
  use string
  use data
  use parallel
  implicit none
  integer :: iunit,ierr,size,ibuffer
  character(len=80) :: buffer
  character(len=str_short) :: name
  character(len=str_medium) :: file
  real(SP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_), intent(in) :: scalar
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: disp
  logical :: file_is_there
  
  ! Generate the file
  !call CREATE_FOLDER("ensight-3D/" // trim(adjustl(name)))
  call system("mkdir -p ensight-3D/" // trim(adjustl(name)))
  file = "ensight-3D/" // trim(adjustl(name)) // "/" // trim(adjustl(name)) // "."
  if (use_mpiiofs) file = trim(mpiiofs) // ":" // trim(file)
  write(file(len_trim(file)+1:len_trim(file)+6),'(i6.6)') nout_time
  
  ! Open the file
  inquire(file=file,exist=file_is_there)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(file,mpi_info,ierr)
  time_open = time_open-MPI_Wtime()
  call MPI_FILE_OPEN(comm,file,MPI_MODE_WRONLY+MPI_MODE_CREATE,mpi_info,iunit,ierr)
  time_open = time_open+MPI_Wtime()
  
  ! Write header (only root)
  time_write = time_write-MPI_Wtime()
  if (irank.eq.iroot) then
     buffer = trim(adjustl(name))
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'part'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     ibuffer = 1
     size = 1
     call MPI_FILE_WRITE(iunit,ibuffer,size,MPI_INTEGER,status,ierr)
     buffer = 'block'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
  end if
  
  ! Write the file
  disp = 3*80+4
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview,"native",mpi_info,ierr)
  call MPI_FILE_WRITE_ALL(iunit,scalar,datasize,MPI_REAL_SP,status,ierr)
  time_write = time_write+MPI_Wtime()
  
  ! Close the file
  time_close = time_close-MPI_Wtime()
  call MPI_FILE_CLOSE(iunit,ierr)
  time_close = time_close+MPI_Wtime()
  
  ! One more file written
  nfiles = nfiles+1
  
  return
end subroutine dump_ensight_str_3D_scalar


! ========================================= !
! Dump 3D binary ensight gold data - vector !
! ========================================= !
subroutine dump_ensight_str_3D_vector(vec,name)
  use dump_ensight_str
  use string
  use parallel
  implicit none
  integer :: iunit,ierr,size,ibuffer
  character(len=80) :: buffer
  character(len=str_short) :: name
  character(len=str_medium) :: file
  real(SP), dimension(imin_:imax_,jmin_:jmax_,kmin_:kmax_,1:3), intent(in) :: vec
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: disp
  logical :: file_is_there
  integer :: var
  
  ! Generate the file
  !call CREATE_FOLDER("ensight-3D/" // trim(adjustl(name)))
  call system("mkdir -p ensight-3D/" // trim(adjustl(name)))
  file = "ensight-3D/" // trim(adjustl(name)) // "/" // trim(adjustl(name)) // "."
  if (use_mpiiofs) file = trim(mpiiofs) // ":" // trim(file)
  write(file(len_trim(file)+1:len_trim(file)+6),'(i6.6)') nout_time
  
  ! Open the file
  inquire(file=file,exist=file_is_there)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(file,mpi_info,ierr)
  time_open = time_open-MPI_Wtime()
  call MPI_FILE_OPEN(comm,file,MPI_MODE_WRONLY+MPI_MODE_CREATE,mpi_info,iunit,ierr)
  time_open = time_open+MPI_Wtime()
  
  ! Write header (only root)
  time_write = time_write-MPI_Wtime()
  if (irank.eq.iroot) then
     buffer = trim(adjustl(name))
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'part'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     ibuffer = 1
     size = 1
     call MPI_FILE_WRITE(iunit,ibuffer,size,MPI_INTEGER,status,ierr)
     buffer = 'block'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
  end if
  
  ! Write the data
  do var=1,3
     disp = 3*80+4+(var-1)*gdatasize*4
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview,"native",mpi_info,ierr)
     call MPI_FILE_WRITE_ALL(iunit,vec(:,:,:,var),datasize,MPI_REAL_SP,status,ierr)
  end do
  time_write = time_write+MPI_Wtime()
  
  ! Close the file
  time_close = time_close-MPI_Wtime()
  call MPI_FILE_CLOSE(iunit,ierr)
  time_close = time_close+MPI_Wtime()
  
  ! One more file written
  nfiles = nfiles+1
  
  return
end subroutine dump_ensight_str_3D_vector

