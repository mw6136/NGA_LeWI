module dnsbox_corr_mpi
  use dnsbox_corr
  implicit none

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
  
contains
  
  subroutine dnsbox_corr_init_topology(xper,yper,zper)
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
    call parser_read('Processors along X',npx)
    call parser_read('Processors along Y',npy)
    call parser_read('Processors along Z',npz)
    
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
  end subroutine dnsbox_corr_init_topology


  ! =============================== !
  ! Partition the CPUs for MPI      !
  ! =============================== !
  subroutine dnsbox_corr_partition_init()
    implicit none
    
    integer :: q,r

    ! Set initial partitions along x
    if (npx.gt.nx) stop "dnsbox_corr_partition_init: nx has to be greater than or equal to npx"
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
    if (npy.gt.ny) stop "dnsbox_corr_partition_init: ny has to be greater than or equal to npy"
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
    if (npz.gt.nz) stop "dnsbox_corr_partition_init: nz has to be greater than or equal to npz"
    if (icyl.eq.1 .and. npz.ne.1) stop 'dnsbox_corr_partition_init: npz has to be equal to 1 for cylindrical configuration'
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

  end subroutine dnsbox_corr_partition_init


  ! =============================== !
  ! Initialize MPI data             !
  ! =============================== !
  subroutine dnsbox_corr_MPI_data_init
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
  end subroutine dnsbox_corr_MPI_data_init

end module dnsbox_corr_mpi



! ================================= !
! Initialize everything for MPI     !
! ================================= !
subroutine dnsbox_corr_MPI_init
  use dnsbox_corr_mpi
  implicit none

  call dnsbox_corr_init_topology(xper,yper,zper)
  call dnsbox_corr_partition_init
  call dnsbox_corr_MPI_data_init

end subroutine dnsbox_corr_MPI_init
