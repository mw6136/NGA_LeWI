module dump_ensight
  use parallel
  use precision
  use geometry
  use partition
  use masks
  use fileio
  use unstruct
  implicit none
  
  ! Time info
  integer :: nout_time
  real(WP), dimension(:), pointer :: out_time
  
  ! SP buffers
  real(SP), dimension(:), pointer :: buffer_hexa
  real(SP), dimension(:), pointer :: buffer_wedge
  real(SP), dimension(:,:), pointer :: buffer3_hexa
  real(SP), dimension(:,:), pointer :: buffer3_wedge
  
  ! Fileviews
  integer :: fileview_hexa
  integer :: fileview_wedge
  integer :: fileview_node
  integer :: fileview_hexa_conn
  integer :: fileview_wedge_conn
  integer, dimension(:), pointer :: map
  
  ! Combustion model
  character(len=str_medium), dimension(:), pointer :: chemtable_ensight_vars
  integer :: n_ctable_ens_vars
  
end module dump_ensight


! ================================================= !
! Dump 3D binary ensight gold data - Initialization !
! ================================================= !
subroutine dump_ensight_3D_init
  use dump_ensight
  use parser
  use parallel
  use combustion
  use time_info
  implicit none
  integer :: ierr,i
  logical :: file_is_there, ens_vars_def
  integer, dimension(:), pointer :: blocklength
  
  ! Create & Start the timer
  call timing_create('ensight')
  call timing_start ('ensight')

  ! If chemtable combustion, get variables to dump
  if (combust) then
     if (trim(chemistry) .eq.'chemtable') then
        call parser_is_defined('Ensight variables', ens_vars_def)
        if (ens_vars_def) then
           call parser_getsize('Ensight variables', n_ctable_ens_vars)
           allocate(chemtable_ensight_vars(n_ctable_ens_vars))
           call parser_read('Ensight variables',chemtable_ensight_vars)
        else
           n_ctable_ens_vars=(7)
           allocate(chemtable_ensight_vars(n_ctable_ens_vars))
           chemtable_ensight_vars(1) = 'Y_F'
           chemtable_ensight_vars(2) = 'Y_O2'
           chemtable_ensight_vars(3) = 'Y_CO'
           chemtable_ensight_vars(4) = 'Y_CO2'
           chemtable_ensight_vars(5) = 'Y_H2'
           chemtable_ensight_vars(6) = 'Y_H2O'
           chemtable_ensight_vars(7) = 'Y_OH'
        end if
     end if
  end if
  
  ! Allocate buffers
  allocate(buffer_hexa(ncells_hexa_))
  allocate(buffer_wedge(ncells_wedge_))
  allocate(buffer3_hexa(ncells_hexa_,3))
  allocate(buffer3_wedge(ncells_wedge_,3))

  ! For MPI-1 derived datatypes
  allocate(blocklength(max(nnodes_,ncells_hexa_,ncells_wedge_)))
  
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
  
!!$  ! Create the views
!!$  allocate(map(1:ncells_hexa_))
!!$  map = hexa_-1
!!$  call MPI_TYPE_CREATE_INDEXED_BLOCK(ncells_hexa_,1,map,MPI_REAL_SP,fileview_hexa,ierr)
!!$  call MPI_TYPE_COMMIT(fileview_hexa,ierr)
!!$  map = map * 8
!!$  call MPI_TYPE_CREATE_INDEXED_BLOCK(ncells_hexa_,8,map,MPI_INTEGER,fileview_hexa_conn,ierr)
!!$  call MPI_TYPE_COMMIT(fileview_hexa_conn,ierr)
!!$  deallocate(map)
!!$  allocate(map(1:ncells_wedge_))
!!$  map = wedge_-1
!!$  call MPI_TYPE_CREATE_INDEXED_BLOCK(ncells_wedge_,1,map,MPI_REAL_SP,fileview_wedge,ierr)
!!$  call MPI_TYPE_COMMIT(fileview_wedge,ierr)
!!$  map = map * 6
!!$  call MPI_TYPE_CREATE_INDEXED_BLOCK(ncells_wedge_,6,map,MPI_INTEGER,fileview_wedge_conn,ierr)
!!$  call MPI_TYPE_COMMIT(fileview_wedge_conn,ierr)
!!$  deallocate(map)
!!$  allocate(map(1:nnodes_))
!!$  map = nodes_-1
!!$  !call MPI_TYPE_CREATE_INDEXED_BLOCK(nnodes_,1,map,MPI_REAL_SP,fileview_node,ierr)
!!$  blocklength = 1
!!$  call MPI_TYPE_INDEXED(nnodes_,blocklength,map,MPI_INTEGER,fileview_node,ierr)
!!$  call MPI_TYPE_COMMIT(fileview_node,ierr)
!!$  deallocate(blocklength)
!!$
  ! Create the views  - hex
  allocate(map(ncells_hexa_))
  map = hexa_-1
  blocklength = 1
  call MPI_TYPE_INDEXED(ncells_hexa_,blocklength,map,MPI_REAL_SP,fileview_hexa,ierr)
  call MPI_TYPE_COMMIT(fileview_hexa,ierr)
  map = map * 8
  blocklength = 8
  call MPI_TYPE_INDEXED(ncells_hexa_,blocklength,map,MPI_INTEGER,fileview_hexa_conn,ierr)
  call MPI_TYPE_COMMIT(fileview_hexa_conn,ierr)
  deallocate(map)
  ! Create the views  - wedge
  if (ncells_wedge_.gt.0) then
     allocate(map(ncells_wedge_))
     map = wedge_-1
     blocklength = 1
     call MPI_TYPE_INDEXED(ncells_wedge_,blocklength,map,MPI_REAL_SP,fileview_wedge,ierr)
     call MPI_TYPE_COMMIT(fileview_wedge,ierr)
     map = map * 6
     blocklength = 6
     call MPI_TYPE_INDEXED(ncells_wedge_,blocklength,map,MPI_INTEGER,fileview_wedge_conn,ierr)
     call MPI_TYPE_COMMIT(fileview_wedge_conn,ierr)
     deallocate(map)
  else
     allocate(map(1))
     map = 0
     blocklength = 1
     call MPI_TYPE_INDEXED(1,blocklength,map,MPI_REAL_SP,fileview_wedge,ierr)
     call MPI_TYPE_COMMIT(fileview_wedge,ierr)
     map = map * 6
     blocklength = 6
     call MPI_TYPE_INDEXED(1,blocklength,map,MPI_INTEGER,fileview_wedge_conn,ierr)
     call MPI_TYPE_COMMIT(fileview_wedge_conn,ierr)
     deallocate(map)
  end if
  ! Create the views  - nodes
  allocate(map(nnodes_))
  map = nodes_-1
  blocklength = 1
  call MPI_TYPE_INDEXED(nnodes_,blocklength,map,MPI_INTEGER,fileview_node,ierr)
  call MPI_TYPE_COMMIT(fileview_node,ierr)
  deallocate(blocklength)

  ! Stop the timer
  call timing_stop ('ensight')

  ! Write the geometry
  call dump_ensight_3D_geometry
  
  ! Save the first field
  call dump_ensight_3D
  
  return
end subroutine dump_ensight_3D_init


! ================================ !
! Dump 3D binary ensight gold data !
! ================================ !
subroutine dump_ensight_3D
  use dump_ensight
  use data
  use partition
  use string
  use combustion
  use sgsmodel
  use interpolate
  use memory
  use pressure
  use soot
  use pah
  use finitechem
  use strainrate
  implicit none
  character(len=str_short) :: name
  integer :: i,j,k,n
  ! Start the timer
  call timing_start('ensight')

  ! Update the case file
  call dump_ensight_3D_case
  
  ! Pressure
  name = "P"
  call dump_ensight_3D_scalar(P,name)
  
  ! Viscosity
  name = "VISC"
  call dump_ensight_3D_scalar(VISC,name)
  if (pope_crit) then
     name = "KSGS"
     call dump_ensight_3D_scalar(k_sgs,name)
  end if
  
  ! Velocity
  name = "V"
  if (icyl.eq.0) then
     call dump_ensight_3D_vector(Ui,Vi,Wi,name)
  else
     do k=kmin_,kmax_
        do j=jmin_,jmax_
           do i=imin_,imax_
              tmp1(i,j,k) = Ui(i,j,k)
              tmp2(i,j,k) = Vi(i,j,k) * cos(zm(k)) - Wi(i,j,k) * sin(zm(k))
              tmp3(i,j,k) = Vi(i,j,k) * sin(zm(k)) + Wi(i,j,k) * cos(zm(k))
           end do
        end do
     end do
     call dump_ensight_3D_vector(tmp1,tmp2,tmp3,name)
  end if
  
  ! Scalars
  do i=1,nscalar
     name = SC_name(i)
     if (name(1:2).ne.'S_') call dump_ensight_3D_scalar(SC(:,:,:,i),name)
  end do
!!$  if (nscalar.ge.1) then
!!$     name = 'DIFF'
!!$     call dump_ensight_3D_scalar(DIFF(:,:,:,1),name)
!!$  end if

  ! Compressible
  if (compressible .and. .not.combust) then
     name = 'RHO'
     call dump_ensight_3D_scalar(RHO,name)
  end if
  
  ! Combustion
  if (combust) then
     if (isc_T.eq.0) then
        name = 'Temp'
        call dump_ensight_3D_scalar(T,name)
     end if
     name = 'RHO'
     call dump_ensight_3D_scalar(RHO,name)
     name = 'dRHO'
     call dump_ensight_3D_scalar(dRHO,name)
     if (trim(chemistry).eq.'finite chem') then
        name = "HR"
        call dump_ensight_3D_scalar(HR_gp,name)
     end if
!!$     if (isc_ZVAR.eq.0 .and. isc_ZMIX.eq.0 .and. use_ZVAR) then
!!$        name = 'ZVAR'
!!$        call dump_ensight_3D_scalar(ZVAR,name)
!!$     end if
     if (isc_ZMIX.ne.0) then
        ! Output it as CHI = DIFF*|nabla(Z)|^2
        name = 'CHI'
!!$        CHI = DIFF(:,:,:,isc_ZMIX)*CHI
        call dump_ensight_3D_scalar(CHI,name)
!!$        where(DIFF(:,:,:,isc_ZMIX).ne.0.0_WP) CHI = CHI/DIFF(:,:,:,isc_ZMIX)
     end if
     if (trim(chemistry).eq.'chemtable') then
        do i=1,n_ctable_ens_vars
           name = trim(chemtable_ensight_vars(i))
           call chemtable_lookup(name,tmp1)
           call dump_ensight_3D_scalar(tmp1,name)
        end do
     end if
  end if
  
  ! Soot
  if (use_soot) then
     name = "numdens"
     call dump_ensight_3D_scalar(numdens,name)
     name = "volfrac"
     call dump_ensight_3D_scalar(volfrac,name)
     name = "partdiam"
     call dump_ensight_3D_scalar(partdiam,name)
     name = "partaggr"
     call dump_ensight_3D_scalar(partaggr,name)
     name = "intermit"
     call dump_ensight_3D_scalar(intermit,name)
     if (.not.use_pah) then
        name = "Y_PAH"
        call chemtable_lookup('Y_PAH',tmp1)
        call dump_ensight_3D_scalar(tmp1,name)
     end if
  end if

!!$  ! Divergence
!!$  name = "divg"
!!$  tmp1(imin_:imax_,jmin_:jmax_,kmin_:kmax_) = divg
!!$  call dump_ensight_3D_scalar(tmp1,name)
  
!!$  ! PQR criterion
!!$  call strainrate_PQRcrit
!!$  name = "Pcrit"
!!$  call dump_ensight_3D_scalar(Pcrit,name)
!!$  name = "Qcrit"
!!$  call dump_ensight_3D_scalar(Qcrit,name)
!!$  name = "Rcrit"
!!$  call dump_ensight_3D_scalar(Rcrit,name)
  
  ! Stop the timer
  call timing_stop('ensight')

  return
end subroutine dump_ensight_3D


! =========================================== !
! Dump 3D binary ensight gold data - geometry !
! => one processor only - test before         !
! =========================================== !
subroutine dump_ensight_3D_geometry
  use dump_ensight
  implicit none
  integer :: ierr,ibuffer,iunit,i
  character(len=80) :: buffer
  character(len=str_medium) :: file
  real(SP), dimension(:), pointer :: xbuf,ybuf,zbuf
  real(SP) :: max_x,max_y,max_z
  real(SP) :: min_x,min_y,min_z
  logical  :: file_is_there
  integer  :: size
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer, dimension(MPI_STATUS_SIZE) :: status
  
  ! Get single precision mesh
  allocate(xbuf(nnodes_),ybuf(nnodes_),zbuf(nnodes_))
  if (icyl.eq.0) then
     do i=1,nnodes_
        xbuf(i)=x(unstr2str(i,1))
        ybuf(i)=y(unstr2str(i,2))
        zbuf(i)=z(unstr2str(i,3))
     end do
  else
     do i=1,nnodes_
        xbuf(i)=x(unstr2str(i,1))
        ybuf(i)=y(unstr2str(i,2))*cos(z(unstr2str(i,3)))
        zbuf(i)=y(unstr2str(i,2))*sin(z(unstr2str(i,3)))
     end do
  end if
  min_x=minval(xbuf)
  min_y=minval(ybuf)
  min_z=minval(zbuf)
  max_x=maxval(xbuf)
  max_y=maxval(ybuf)
  max_z=maxval(zbuf)
  
  ! Generate the geometry file in parallel
  file = "ensight-3D/geometry"
  if (use_mpiiofs) file=trim(mpiiofs) // ":" // trim(file)
  inquire(file=file,exist=file_is_there)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(file,MPI_INFO_NULL,ierr)
  call MPI_FILE_OPEN(comm,file,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),mpi_info,iunit,ierr)
  
  ! Write header (only root)
  if (irank.eq.iroot) then
     ! Global header
     buffer = 'C Binary'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'Ensight Gold Geometry File'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'Unstructured Geometry from ARTS'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'node id given'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     buffer = 'element id given'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     ! Extents
     buffer = 'extents'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     size = 1
     call MPI_FILE_WRITE(iunit,min_x,size,MPI_REAL_SP,status,ierr)
     call MPI_FILE_WRITE(iunit,max_x,size,MPI_REAL_SP,status,ierr)
     call MPI_FILE_WRITE(iunit,min_y,size,MPI_REAL_SP,status,ierr)
     call MPI_FILE_WRITE(iunit,max_y,size,MPI_REAL_SP,status,ierr)
     call MPI_FILE_WRITE(iunit,min_z,size,MPI_REAL_SP,status,ierr)
     call MPI_FILE_WRITE(iunit,max_z,size,MPI_REAL_SP,status,ierr)
     ! Part header
     buffer = 'part'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     ibuffer = 1
     size = 1
     call MPI_FILE_WRITE(iunit,ibuffer,size,MPI_INTEGER,status,ierr)
     buffer = 'ARTS 3D geometry'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     ! Nodes list
     buffer = 'coordinates'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     size = 1
     call MPI_FILE_WRITE(iunit,nnodes,size,MPI_INTEGER,status,ierr)
  end if
  
  ! Write the node positions
  disp = 752
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_INTEGER,fileview_node,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,nodes_,nnodes_,MPI_INTEGER,status,ierr)
  disp = disp+4*nnodes
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_node,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,xbuf,nnodes_,MPI_REAL_SP,status,ierr)
  disp = disp+4*nnodes
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_node,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,ybuf,nnodes_,MPI_REAL_SP,status,ierr)
  disp = disp+4*nnodes
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_node,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,zbuf,nnodes_,MPI_REAL_SP,status,ierr)
  
  ! Write the hexa connectivity
  disp = disp+4*nnodes
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_CHARACTER,MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
  if (irank.eq.iroot) then
     buffer = 'hexa8'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     size = 1
     call MPI_FILE_WRITE(iunit,ncells_hexa,size,MPI_INTEGER,status,ierr)
     !do i=1,ncells_hexa
     !   call MPI_FILE_WRITE(iunit,i,size,MPI_INTEGER,status,ierr)
     !end do
  end if
  disp = disp+84
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_INTEGER,fileview_hexa_conn,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,hexa_,ncells_hexa_,MPI_INTEGER,status,ierr)
  disp = disp+4*ncells_hexa
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_INTEGER,fileview_hexa_conn,"native",MPI_INFO_NULL,ierr)
  size = 8*ncells_hexa_
  call MPI_FILE_WRITE_ALL(iunit,conn_hexa,size,MPI_INTEGER,status,ierr)
  
  ! Write the file - wedge
  if (ncells_wedge.gt.0) then
     disp = disp+8*4*ncells_hexa
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_CHARACTER,MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
     if (irank.eq.iroot) then
        buffer = 'penta6'
        size = 80
        call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
        size = 1
        call MPI_FILE_WRITE(iunit,ncells_wedge,size,MPI_INTEGER,status,ierr)
        !do i=1,ncells_wedge
        !   call MPI_FILE_WRITE(iunit,i,size,MPI_INTEGER,status,ierr)
        !end do
     end if
     disp = disp+84
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_INTEGER,fileview_wedge_conn,"native",MPI_INFO_NULL,ierr)
     if (ncells_wedge_.gt.0) call MPI_FILE_WRITE(iunit,wedge_,ncells_wedge_,MPI_INTEGER,status,ierr)
     disp = disp+4*ncells_wedge
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_INTEGER,fileview_wedge_conn,"native",MPI_INFO_NULL,ierr)
     size = 6*ncells_wedge_
     if (ncells_wedge_.gt.0) call MPI_FILE_WRITE(iunit,conn_wedge,size,MPI_INTEGER,status,ierr)
  end if
  
  ! Close the file
  call MPI_FILE_CLOSE(iunit,ierr)

  ! Deallocate arrays
  deallocate(xbuf,ybuf,zbuf)
  
  return
end subroutine dump_ensight_3D_geometry


! ============================================ !
! Dump 3D binary ensight gold data - case file !
! ============================================ !
subroutine dump_ensight_3D_case
  use dump_ensight
  use combustion
  use sgsmodel
  use time_info
  implicit none
  integer :: iunit,ierr,i,n
  real(WP), dimension(:), allocatable :: buffer
  character(len=80) :: str, str2
  character(len=str_short) :: temp, name
  
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
!!$     ! Diffusivity
!!$     str='scalar per element: 1 Diffusivity DIFF/DIFF.******'
!!$     if (nscalar.ge.1) write(iunit,'(a80)') str
     ! Scalars
     do i=1,nscalar
        name = SC_name(i)
        if (name(1:2).ne.'S_') then
           str = 'scalar per element: 1 '//trim(SC_name(i))//' '//trim(SC_name(i))//'/'//trim(SC_name(i))//'.******'
           write(iunit,'(a80)') str
        end if
     end do
     ! Compressible
     if (compressible .and. .not.combust) then
        str = 'scalar per element: 1 RHO RHO/RHO.******'
        write(iunit,'(a80)') str
     endif
     ! Combustion stuff
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
!!$        str = 'scalar per element: 1 ZVAR ZVAR/ZVAR.******'
!!$        if (isc_ZVAR.eq.0 .and. isc_ZMIX2.eq.0 .and. use_ZVAR) write(iunit,'(a80)') str
        str = 'scalar per element: 1 CHI CHI/CHI.******'
        if (isc_ZMIX.ne.0) write(iunit,'(a80)') str
        if (trim(chemistry).eq.'chemtable') then
           do i =1,n_ctable_ens_vars
              str2 = trim(chemtable_ensight_vars(i)) 
              str ='scalar per element: 1 '//trim(str2)//' '//trim(str2)//'/'//trim(str2)//'.******'
              write(iunit,'(a80)') str
           end do
        end if
     end if
     ! Velocity
     str='vector per element: 1 V V/V.******'
     write(iunit,'(a80)') str
     ! Soot stuff
     if (use_soot) then
        str='scalar per element: 1 numdens numdens/numdens.******'
        write(iunit,'(a80)') str
        str='scalar per element: 1 volfrac volfrac/volfrac.******'
        write(iunit,'(a80)') str
        str='scalar per element: 1 partdiam partdiam/partdiam.******'
        write(iunit,'(a80)') str
        str='scalar per element: 1 partaggr partaggr/partaggr.******'
        write(iunit,'(a80)') str
        str='scalar per element: 1 intermit intermit/intermit.******'
        write(iunit,'(a80)') str
        if (.not.use_pah) then
           str='scalar per element: 1 Y_PAH Y_PAH/Y_PAH.******'
           write(iunit,'(a80)') str
        end if
     end if
!!$     ! Divergence
!!$     str='scalar per element: 1 divg divg/divg.******'
!!$     write(iunit,'(a80)') str
!!$     ! PQR criterion
!!$     str='scalar per element: 1 Pcrit Pcrit/Pcrit.******'
!!$     write(iunit,'(a80)') str
!!$     str='scalar per element: 1 Qcrit Qcrit/Qcrit.******'
!!$     write(iunit,'(a80)') str
!!$     str='scalar per element: 1 Rcrit Rcrit/Rcrit.******'
!!$     write(iunit,'(a80)') str
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
end subroutine dump_ensight_3D_case


! ========================================= !
! Dump 3D binary ensight gold data - scalar !
! ========================================= !
subroutine dump_ensight_3D_scalar(scalar,name)
  use dump_ensight
  use string
  use data
  use parallel
  implicit none
  integer :: iunit,ierr,size,ibuffer,i,j,k,count_wedge_,count_hexa_
  character(len=80) :: buffer
  character(len=str_short) :: name
  character(len=str_medium) :: file
  real(WP),dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: scalar
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: disp
  logical :: file_is_there
  
  ! Extract the data
  count_wedge_ = 0
  count_hexa_ = 0
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              if (icyl.eq.1 .and. j.eq.jmin) then
                 count_wedge_ = count_wedge_ + 1
                 buffer_wedge(count_wedge_) = real(scalar(i,j,k),SP)
              else
                 count_hexa_ = count_hexa_ + 1
                 buffer_hexa(count_hexa_) = real(scalar(i,j,k),SP)
              end if
           end if
        end do
     end do
  end do
  
  ! Generate the file
  !call CREATE_FOLDER("ensight-3D/" // trim(adjustl(name)))
  call system("mkdir -p ensight-3D/" // trim(adjustl(name)))
  file = "ensight-3D/" // trim(adjustl(name)) // "/" // trim(adjustl(name)) // "."
  if (use_mpiiofs) file = trim(mpiiofs) // ":" // trim(file)
  write(file(len_trim(file)+1:len_trim(file)+6),'(i6.6)') nout_time
  
  ! Open the file
  inquire(file=file,exist=file_is_there)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(file,MPI_INFO_NULL,ierr)
  call MPI_FILE_OPEN(comm,file,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),mpi_info,iunit,ierr)
  ! Write header (only root)
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
     buffer = 'hexa8'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
  end if
  
  ! Write the file - hexa
  disp = 3*80+4
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_hexa,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,buffer_hexa,ncells_hexa_,MPI_REAL_SP,status,ierr)
  
  ! Write the file - wedge
  if (ncells_wedge.gt.0) then
     disp = 3*80+4+4*ncells_hexa
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_CHARACTER,MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
     if (irank.eq.iroot) then
        buffer = 'penta6'
        size = 80
        call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     end if
     disp = 3*80+4+80+4*ncells_hexa
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_wedge,"native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(iunit,buffer_wedge,ncells_wedge_,MPI_REAL_SP,status,ierr)
  end if
  
  ! Close the file
  call MPI_FILE_CLOSE(iunit,ierr)
  
  return
end subroutine dump_ensight_3D_scalar


! ========================================= !
! Dump 3D binary ensight gold data - vector !
! ========================================= !
subroutine dump_ensight_3D_vector(vec1,vec2,vec3,name)
  use dump_ensight
  use string
  use parallel
  implicit none
  integer :: iunit,ierr,size,ibuffer,i,j,k,count_wedge_,count_hexa_
  character(len=80) :: buffer
  character(len=str_short) :: name
  character(len=str_medium) :: file
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: vec1
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: vec2
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_) :: vec3
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer(kind=MPI_OFFSET_KIND) :: disp
  logical :: file_is_there
  
  ! Extract the data
  count_wedge_ = 0
  count_hexa_ = 0
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           if (mask(i,j).eq.0) then
              if (icyl.eq.1 .and. j.eq.jmin) then
                 count_wedge_ = count_wedge_ + 1
                 buffer3_wedge(count_wedge_,1) = real(vec1(i,j,k),SP)
                 buffer3_wedge(count_wedge_,2) = real(vec2(i,j,k),SP)
                 buffer3_wedge(count_wedge_,3) = real(vec3(i,j,k),SP)
              else
                 count_hexa_ = count_hexa_ + 1
                 buffer3_hexa(count_hexa_,1) = real(vec1(i,j,k),SP)
                 buffer3_hexa(count_hexa_,2) = real(vec2(i,j,k),SP)
                 buffer3_hexa(count_hexa_,3) = real(vec3(i,j,k),SP)
              end if
           end if
        end do
     end do
  end do
  
  ! Generate the file
  !call CREATE_FOLDER("ensight-3D/" // trim(adjustl(name)))
  call system("mkdir -p ensight-3D/" // trim(adjustl(name)))
  file = "ensight-3D/" // trim(adjustl(name)) // "/" // trim(adjustl(name)) // "."
  if (use_mpiiofs) file = trim(mpiiofs) // ":" // trim(file)
  write(file(len_trim(file)+1:len_trim(file)+6),'(i6.6)') nout_time
  
  ! Open the file
  inquire(file=file,exist=file_is_there)
  if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(file,MPI_INFO_NULL,ierr)
  call MPI_FILE_OPEN(comm,file,IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),mpi_info,iunit,ierr)
  
  ! Write header (only root)
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
     buffer = 'hexa8'
     size = 80
     call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
  end if
  
  ! Write the data
  disp = 3*80+4+0*ncells_hexa*4
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_hexa,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,buffer3_hexa(:,1),ncells_hexa_,MPI_REAL_SP,status,ierr)
  disp = 3*80+4+1*ncells_hexa*4
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_hexa,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,buffer3_hexa(:,2),ncells_hexa_,MPI_REAL_SP,status,ierr)
  disp = 3*80+4+2*ncells_hexa*4
  call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_hexa,"native",MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(iunit,buffer3_hexa(:,3),ncells_hexa_,MPI_REAL_SP,status,ierr)
  
  ! Write the file - wedge
  if (ncells_wedge.gt.0) then
     disp = 3*80+4+3*ncells_hexa*4+0*ncells_wedge*4
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_CHARACTER,MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
     if (irank.eq.iroot) then
        buffer = 'penta6'
        size = 80
        call MPI_FILE_WRITE(iunit,buffer,size,MPI_CHARACTER,status,ierr)
     end if
     disp = 3*80+4+3*ncells_hexa*4+80+0*ncells_wedge*4
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_wedge,"native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(iunit,buffer3_wedge(:,1),ncells_wedge_,MPI_REAL_SP,status,ierr)
     disp = 3*80+4+3*ncells_hexa*4+80+1*ncells_wedge*4
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_wedge,"native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(iunit,buffer3_wedge(:,2),ncells_wedge_,MPI_REAL_SP,status,ierr)
     disp = 3*80+4+3*ncells_hexa*4+80+2*ncells_wedge*4
     call MPI_FILE_SET_VIEW(iunit,disp,MPI_REAL_SP,fileview_wedge,"native",MPI_INFO_NULL,ierr)
     call MPI_FILE_WRITE_ALL(iunit,buffer3_wedge(:,3),ncells_wedge_,MPI_REAL_SP,status,ierr)
  end if
  
  ! Close the file
  call MPI_FILE_CLOSE(iunit,ierr)
  
  return
end subroutine dump_ensight_3D_vector
