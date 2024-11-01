module param
  use precision
  use string
  use parallel
  implicit none

  ! --- CONFIG FILE ---

  ! Cylindrical or cartesian
  integer :: icyl
  ! Number of grid points
  integer :: nx,ny,nz
  ! Periodicity
  integer :: xper,yper,zper

  ! Grid/mesh variables for data file
  real(WP), dimension(:), pointer :: x
  real(WP), dimension(:), pointer :: y
  real(WP), dimension(:), pointer :: z
  real(WP), dimension(:), pointer :: xm
  real(WP), dimension(:), pointer :: ym
  real(WP), dimension(:), pointer :: zm
  integer, dimension(:,:), pointer :: mask

  ! --- DATA FILE ---

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nvar
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:,:), pointer :: data
  
  ! --- OPTDATA FILE ---
  
  ! Data array with all the optional variables
  integer :: nod
  character(len=str_short), dimension(:), pointer :: OD_names
  real(WP), dimension(:,:,:,:), pointer :: OD
  
  ! --- INFLOW FILE ---

  ! Inflow grid related parameters
  integer :: ntime
  real(WP), dimension(:), pointer :: t_inflow
  
  ! Number of time steps
  real(WP) :: dt_inflow,time_inflow

  ! Inflow array with all variables
  integer :: nvar_inflow
  character(len=str_short), dimension(:), pointer :: names_inflow
  real(WP), dimension(:,:,:,:), pointer :: inflow
  
  ! --- CHEMTABLE ---

  ! Coordinates of the chemtable
  integer :: n1,n2,n3
  real(WP), dimension(:), pointer :: x1,x2,x3
  ! Mask of the chemtable
  integer, dimension(:,:,:), pointer :: chem_mask

  ! Names in the chemtable
  integer :: nvar_chem
  character(len=str_medium), dimension(:), pointer :: names_chem

  ! Chemtable model
  character(len=str_medium) :: combModel

  ! Table of variables
  real(WP), dimension(:,:,:,:), pointer :: table
  
  ! JFM added 2/22/2016
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
  
  ! Array with each variables
  integer :: MPI_IO_NVARS
  type(MPI_IO_VAR), dimension(:), pointer :: MPI_IO_DATA

  ! Local Cartesian coordinates
  integer :: imin_, imax_, nx_
  integer :: jmin_, jmax_, ny_
  integer :: kmin_, kmax_, nz_

  logical :: use_MPI = .false.

contains

  subroutine param_init_topology(xper,yper,zper)
    use parser
    implicit none
    integer :: ierr
    integer, dimension(3) :: dims
    logical, dimension(3) :: isper
    logical, dimension(3) :: dir
    logical :: reorder
    integer :: ndims
    integer, dimension(3) :: coords
    integer, intent(in) :: xper,yper,zper
    
    ! Save periodicity
    periodicity(1) = xper
    periodicity(2) = yper
    periodicity(3) = zper
    
    ! Read topology from input file
    call parser_read('Init processors along X',npx)
    call parser_read('Init processors along Y',npy)
    call parser_read('Init processors along Z',npz)
    
    ! Test if nproc is correct
    if (nproc .ne. npx*npy*npz) call parallel_kill('Wrong number of cpus specified in input file')
    
    ! Set MPI topology
    ndims = 3
    dims(1) = npx
    dims(2) = npy
    dims(3) = npz
    if (xper.EQ.1) then
       isper(1) = .true.
    else
       isper(1) = .false.
    end if
    if (yper.EQ.1) then
       isper(2) = .true.
    else
       isper(2) = .false.
    end if
    if (zper.EQ.1) then
       isper(3) = .true.
    else
       isper(3) = .false.
    end if
    reorder = .true.
    call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,isper,reorder,comm,ierr)
    call MPI_COMM_RANK(comm,irank,ierr)
    call MPI_CART_COORDS(comm,irank,ndims,coords,ierr)
    irank = irank + 1
    iproc = coords(1) + 1
    jproc = coords(2) + 1
    kproc = coords(3) + 1
    
    ! Define a root processor at coordinates (1,1,1)
    dims = 0
    call MPI_CART_RANK(comm,dims,iroot,ierr)
    iroot = iroot + 1
    
    ! Create line communicators
    ! Along x
    dir(1) = .true.
    dir(2) = .false.
    dir(3) = .false.
    call MPI_CART_SUB(comm,dir,comm_x,ierr)
    call MPI_COMM_RANK(comm_x,irank_x,ierr)
    ! Along y
    dir(1) = .false.
    dir(2) = .true.
    dir(3) = .false.
    call MPI_CART_SUB(comm,dir,comm_y,ierr)
    call MPI_COMM_RANK(comm_y,irank_y,ierr)
    ! Along z
    dir(1) = .false.
    dir(2) = .false.
    dir(3) = .true.
    call MPI_CART_SUB(comm,dir,comm_z,ierr)
    call MPI_COMM_RANK(comm_z,irank_z,ierr)
    
    ! Create planar communicators
    ! Along xy
    dir(1) = .true.
    dir(2) = .true.
    dir(3) = .false.
    call MPI_CART_SUB(comm,dir,comm_xy,ierr)
    ! Along yz
    dir(1) = .false.
    dir(2) = .true.
    dir(3) = .true.
    call MPI_CART_SUB(comm,dir,comm_yz,ierr)
    ! Along xz
    dir(1) = .true.
    dir(2) = .false.
    dir(3) = .true.
    call MPI_CART_SUB(comm,dir,comm_xz,ierr)
    
    return
  end subroutine param_init_topology


  ! =============================== !
  ! Partition the CPUs for MPI      !
  ! =============================== !
  subroutine param_partition_init()
    implicit none
    
    integer :: q,r

    ! Set initial partitions along x
    if (npx.gt.nx) stop "param_partition_init: nx has to be greater than or equal to npx"
    q = nx/npx
    r = mod(nx,npx)
    if (iproc<=r) then
       nx_   = q+1
       imin_ = 1 + (iproc-1)*(q+1)
    else
       nx_   = q
       imin_ = 1 + r*(q+1) + (iproc-r-1)*q
    end if
    imax_ = imin_ + nx_ - 1

    ! Set initial partitions along y
    if (npy.gt.ny) stop "param_partition_init: ny has to be greater than or equal to npy"
    q = ny/npy
    r = mod(ny,npy)
    if (jproc<=r) then
       ny_   = q+1
       jmin_ = 1 + (jproc-1)*(q+1)
    else
       ny_   = q
       jmin_ = 1 + r*(q+1) + (jproc-r-1)*q
    end if
    jmax_  = jmin_ + ny_ - 1

    ! Set initial partitions along z
    if (npz.gt.nz) stop "param_partition_init: nz has to be greater than or equal to npz"
    if (icyl.eq.1 .and. npz.ne.1) stop 'param_partition_init: npz has to be equal to 1 for cylindrical configuration'
    q = nz/npz
    r = mod(nz,npz)
    if (kproc<=r) then
       nz_   = q+1
       kmin_ = 1 + (kproc-1)*(q+1)
    else
       nz_   = q
       kmin_ = 1 + r*(q+1) + (kproc-r-1)*q
    end if
    kmax_  = kmin_ + nz_ - 1

  end subroutine param_partition_init


  ! =============================== !
  ! Initialize MPI data             !
  ! =============================== !
  subroutine param_MPI_data_init
    implicit none

    integer, dimension(3) :: gsizes, lsizes, start
    integer :: var, ierr

    MPI_IO_NVARS = nvar
    allocate(MPI_IO_DATA(MPI_IO_NVARS))

    ! Assign names and link the variables
    ! Requires 'names' and 'data' to be set in <configuration>_data subroutine
    do var=1,MPI_IO_NVARS
       MPI_IO_DATA(var)%name =  names(var)
       MPI_IO_DATA(var)%var  => data(imin_:imax_,jmin_:jmax_,kmin_:kmax_,var)
    end do

    ! Define global(g) and local(l) sizes
    gsizes(1) = nx
    gsizes(2) = ny
    gsizes(3) = nz
    
    lsizes(1) = nx_
    lsizes(2) = ny_
    lsizes(3) = nz_
    
    ! Define starting points
    start(1) = imin_-1
    start(2) = jmin_-1
    start(3) = kmin_-1
    
    ! Define the view for each variable
    do var=1,MPI_IO_NVARS
       
       call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,start,&
            MPI_ORDER_FORTRAN,MPI_REAL_WP,MPI_IO_DATA(var)%view,ierr)
       call MPI_TYPE_COMMIT(MPI_IO_DATA(var)%view,ierr)
       
    end do
    
    return
  end subroutine param_MPI_data_init


  ! =============================== !
  ! Write the mesh file to the disk !
  ! =============================== !
  subroutine param_write_config(simulation)
    use parser
    implicit none

    character(len=str_medium) :: simulation
    integer :: ierr,iunit
    character(len=str_medium) :: filename

    ! Open the mesh file
    call parser_read('Init config file',filename)
    call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)

    ! Write the mesh
    call BINARY_FILE_WRITE(iunit,simulation,str_medium,kind(simulation),ierr)
    call BINARY_FILE_WRITE(iunit,icyl,1,kind(icyl),ierr)
    call BINARY_FILE_WRITE(iunit,xper,1,kind(xper),ierr)
    call BINARY_FILE_WRITE(iunit,yper,1,kind(yper),ierr)
    call BINARY_FILE_WRITE(iunit,zper,1,kind(zper),ierr)
    call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
    call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
    call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
    call BINARY_FILE_WRITE(iunit,x,nx+1,kind(x),ierr)
    call BINARY_FILE_WRITE(iunit,y,ny+1,kind(y),ierr)
    call BINARY_FILE_WRITE(iunit,z,nz+1,kind(z),ierr)
    call BINARY_FILE_WRITE(iunit,mask,nx*ny,kind(mask),ierr)

    ! Close the file
    call BINARY_FILE_CLOSE(iunit,ierr)

    return
  end subroutine param_write_config


  ! =============================== !
  ! Write the data file to the disk !
  ! =============================== !
  subroutine param_write_data
    use parser
    implicit none

    integer :: ierr,iunit
    integer :: var
    character(len=str_medium) :: filename
    real(WP) :: dt,time

    ! Open the data file
    call parser_read('Init data file',filename)
    call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)

    ! Write sizes
    call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
    call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
    call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
    call BINARY_FILE_WRITE(iunit,nvar,1,kind(nvar),ierr)
    ! Write additional stuff
    dt = 0.0_WP
    time = 0.0_WP
    call BINARY_FILE_WRITE(iunit,dt,1,kind(dt),ierr)
    call BINARY_FILE_WRITE(iunit,time,1,kind(time),ierr)
    ! Write variable names
    do var=1,nvar
       call BINARY_FILE_WRITE(iunit,names(var),str_short,kind(names),ierr)
    end do
    ! Write data field
    do var=1,nvar
       call BINARY_FILE_WRITE_LONG(iunit,data(:,:,:,var),nx*ny*nz,kind(data),ierr)
    end do
    
    call BINARY_FILE_CLOSE(iunit,ierr)

    return
  end subroutine param_write_data


  ! =========================================== !
  ! Write the data file to the disk in parallel !
  ! =========================================== !
  subroutine param_MPI_write_data
    use parser
    implicit none

    integer :: ifile,ierr,var,data_size
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer, dimension(4) :: dims
    integer(kind=MPI_Offset_kind) :: disp
    integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK
    integer(kind=MPI_Offset_kind) :: WP_MOK, var_MOK, str_MOK
    integer(kind=MPI_Offset_kind) :: NVARS_MOK
    character(len=str_medium) :: filename,buffer
    integer :: overwrite
    logical :: file_is_there
    real(WP) :: dt, time
    
    call parser_read('Init data file',filename)
    if (use_mpiiofs) filename = trim(mpiiofs)//":" // trim(filename)

    ! Open the file to write
    inquire(file=filename,exist=file_is_there)
    if (file_is_there .and. irank.eq.iroot) call MPI_FILE_DELETE(filename,mpi_info,ierr)
    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,mpi_info,ifile,ierr)

    ! Write header
    dt = 0.0_WP
    time = 0.0_WP
    if (irank.eq.iroot) then
       ! Write dimensions
       dims(1) = nx
       dims(2) = ny
       dims(3) = nz
       dims(4) = MPI_IO_NVARS
       call MPI_FILE_WRITE(ifile,dims,4,MPI_INTEGER,status,ierr)
       ! Write additional stuff
       call MPI_FILE_WRITE(ifile,dt,1,MPI_REAL_WP,status,ierr)
       call MPI_FILE_WRITE(ifile,time,1,MPI_REAL_WP,status,ierr)
       ! Write variable names
       do var=1,MPI_IO_NVARS
          call MPI_FILE_WRITE(ifile,MPI_IO_DATA(var)%name,str_short,MPI_CHARACTER,status,ierr)
       end do
    end if

    ! Size of local arrays
    data_size = nx_*ny_*nz_
    
    ! Resize some integers so MPI can write even the biggest files
    nx_MOK    = int(nx,          MPI_Offset_kind)
    ny_MOK    = int(ny,          MPI_Offset_kind)
    nz_MOK    = int(nz,          MPI_Offset_kind)
    WP_MOK    = int(WP,          MPI_Offset_kind)
    str_MOK   = int(str_short,   MPI_Offset_kind)
    NVARS_MOK = int(MPI_IO_NVARS,MPI_Offset_kind)
    
    ! Write the data for each variable
    do var=1,MPI_IO_NVARS
       var_MOK = int(var,MPI_Offset_kind)
       disp = 4*4 + str_MOK*NVARS_MOK + 2*WP_MOK + & 
                    nx_MOK*ny_MOK*nz_MOK*WP_MOK*(var_MOK-1)
       call MPI_FILE_SET_VIEW(ifile,disp,MPI_REAL_WP,MPI_IO_DATA(var)%view, &
            "native",mpi_info,ierr)
       call MPI_FILE_WRITE_ALL(ifile,MPI_IO_DATA(var)%var,data_size, &
            MPI_REAL_WP,status,ierr)
    end do
  
    ! Close the file
    call MPI_FILE_CLOSE(ifile,ierr)
    
  end subroutine param_MPI_write_data
  
  
  ! =============================== !
  ! Write the data file to the disk !
  ! =============================== !
  subroutine param_write_optdata
    use parser
    implicit none
    
    integer :: ierr,iunit
    integer :: var
    character(len=str_medium) :: filename
    real(WP) :: dt,time
    
    ! If optdata then continue
    if (.not.associated(OD)) return
    
    ! Open the data file
    call parser_read('Init optional data',filename)
    call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
    
    ! Write sizes
    call BINARY_FILE_WRITE(iunit,nx,1,kind(nx),ierr)
    call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
    call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
    call BINARY_FILE_WRITE(iunit,nod,1,kind(nod),ierr)
    ! Write additional stuff
    dt = 0.0_WP
    time = 0.0_WP
    call BINARY_FILE_WRITE(iunit,dt,1,kind(dt),ierr)
    call BINARY_FILE_WRITE(iunit,time,1,kind(time),ierr)
    ! Write variable names
    do var=1,nod
       call BINARY_FILE_WRITE(iunit,OD_names(var),str_short,kind(OD_names),ierr)
    end do
    ! Write data field
    do var=1,nod
       call BINARY_FILE_WRITE(iunit,OD(:,:,:,var),nx*ny*nz,kind(OD),ierr)
    end do
    
    call BINARY_FILE_CLOSE(iunit,ierr)
    
    return
  end subroutine param_write_optdata
  
  
  ! ================================= !
  ! Write the inflow file to the disk !
  ! ================================= !
  subroutine param_write_inflow
     use parser
     implicit none
     
     integer :: ierr,iunit,var,n
     character(len=str_medium) :: filename
     
     ! If inflow then continue
     if (.not.associated(inflow)) return
          
     ! Read filename and frequency
     call parser_read('Init inflow file',filename)
     call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
     
     ! Write sizes
     call BINARY_FILE_WRITE(iunit,ntime,1,kind(ntime),ierr)
     call BINARY_FILE_WRITE(iunit,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit,nvar_inflow,1,kind(nvar_inflow),ierr)

     ! Write additional stuff
     call BINARY_FILE_WRITE(iunit,dt_inflow,1,kind(dt_inflow),ierr)
     call BINARY_FILE_WRITE(iunit,time_inflow,1,kind(time_inflow),ierr)
     ! Write variable names
     do var=1,nvar_inflow
        call BINARY_FILE_WRITE(iunit,names_inflow(var),str_short,kind(names_inflow),ierr)
     end do
     ! Write the grid
     call BINARY_FILE_WRITE(iunit,icyl,1,kind(icyl),ierr)
     call BINARY_FILE_WRITE(iunit,y,ny+1,kind(y),ierr)
     call BINARY_FILE_WRITE(iunit,z,nz+1,kind(z),ierr)
     ! Write data field
     do n=1,ntime
        do var=1,nvar_inflow
           call BINARY_FILE_WRITE(iunit,inflow(n,:,:,var),ny*nz,kind(inflow),ierr)
        end do
     end do
    
     call BINARY_FILE_CLOSE(iunit,ierr)

     return
  end subroutine param_write_inflow 


  ! ==================================== !
  ! Write the chemtable file to the disk !
  ! ==================================== !
  subroutine param_write_chemtable
    use parser
    implicit none

      integer :: ierr,var,iunit
      character(len=str_medium) :: filename
      
      ! If chemtable then continue
      if (.not.associated(table)) return

      ! Open the data file
      call parser_read('Init chemtable file', filename)
      call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
      
      ! Write sizes
      call BINARY_FILE_WRITE(iunit,n1,1,kind(n1),ierr)
      call BINARY_FILE_WRITE(iunit,n2,1,kind(n2),ierr)
      call BINARY_FILE_WRITE(iunit,n3,1,kind(n3),ierr)
      call BINARY_FILE_WRITE(iunit,nvar_chem,1,kind(nvar_chem),ierr)
      ! Write the axis coordinates
      call BINARY_FILE_WRITE(iunit,x1,n1,kind(x1),ierr)
      call BINARY_FILE_WRITE(iunit,x2,n2,kind(x2),ierr)
      call BINARY_FILE_WRITE(iunit,x3,n3,kind(x3),ierr)
      ! Write the mask of the chemtable
      call BINARY_FILE_WRITE(iunit,chem_mask,n1*n2*n3,kind(chem_mask),ierr)
      ! Write additional stuff
      call BINARY_FILE_WRITE(iunit,combModel,str_medium,kind(combModel),ierr)
      ! Write variable names
      do var=1,nvar_chem
         call BINARY_FILE_WRITE(iunit,names_chem(var),str_medium,kind(names_chem),ierr)
      end do
      ! Write data field
      do var=1,nvar_chem
         call BINARY_FILE_WRITE(iunit,table(:,:,:,var),n1*n2*n3,kind(table),ierr)
      end do
      
      call BINARY_FILE_CLOSE(iunit,ierr)
      
      return
    end subroutine param_write_chemtable
    
end module param


! ================================= !
! Start the parallel module         !
! ================================= !
subroutine param_parallel_init
  use param
  implicit none
  
  call parallel_init
  call parallel_init_filesystem

end subroutine param_parallel_init


! ================================= !
! Initialize everything for MPI     !
! ================================= !
subroutine param_MPI_init
  use param
  implicit none

  use_MPI = .true.
  call param_init_topology(xper,yper,zper)
  call param_partition_init

end subroutine param_MPI_init
