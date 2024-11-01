module inflow
  use boundary
  use borders
  implicit none

  ! Type of inlets
  character(len=str_medium), dimension(:), allocatable :: inlet_type

  ! Plane of the variables at the inflow
  real(WP), dimension(:,:), pointer   :: Uin
  real(WP), dimension(:,:), pointer   :: Vin
  real(WP), dimension(:,:), pointer   :: Win
  real(WP), dimension(:,:,:), pointer :: SCin
  
  ! Bulk velocities
  real(WP), dimension(:), pointer :: ubulk,wbulk

  ! Additive velocity
  real(WP), dimension(:), pointer :: uadd
  
  ! Rescale inflow files
  logical :: rescale_inflow

  ! Add to inflow file
  logical :: additive_inflow

contains

  ! Read the velocity values from the input file
  subroutine inflow_read_velocity
    use parser
    implicit none

    integer  :: nflow,nu,nw,j
    real(WP) :: factor,ymin,ymax
    logical  :: read_bulk
    
    ! Rescaling of the inflow
    call parser_read('Inflow file rescale',rescale_inflow,.false.)
    
    read_bulk = rescale_inflow
    do nflow=1,ninlet
       if ( trim(inlet_type(nflow)).ne.'file' .and. &
            trim(inlet_type(nflow)).ne.'tanh' ) read_bulk = .true.
    end do
    
    if (read_bulk) then
       ! Read the values for 'bulk' and 'laminar' and rescaled 'file'
       call parser_getsize('Inlet u velocity',nu)
       call parser_getsize('Inlet w velocity',nw)
       if (nu.ne.ninlet .or. nw.ne.ninlet) &
            call die('inflow_read_velocity: not enough Inlet u/w velocities')
       allocate(ubulk(ninlet))
       allocate(wbulk(ninlet))
       call parser_read('Inlet u velocity',ubulk)
       call parser_read('Inlet w velocity',wbulk)
    end if

    ! Additive value to inflow
    call parser_read('Inflow file addition',additive_inflow,.false.)

    if (additive_inflow) then
       call parser_getsize('Inlet add velocity',nu)
       if (nu.ne.ninlet) &
          call die('inflow_read_velocity: not enough Inlet add velocities')
       allocate(uadd(ninlet))
       call parser_read('Inlet add velocity',uadd)
    end if
    
    ! Set the values
    do nflow=1,ninlet
       select case(trim(inlet_type(nflow)))

       ! Bulk velocity
       ! Constant velocity for each inflow
       case('bulk')
          do j=max(jmino_,inlet(nflow)%jmino),min(jmaxo_,inlet(nflow)%jmaxo)
             Uin(j,:) = ubulk(nflow)
             Vin(j,:) = 0.0_WP
             Win(j,:) = wbulk(nflow)
          end do

       ! Laminar velocity profile
       ! Quadratic profile for each inflow
       case('laminar')
          ymin = y(inlet(nflow)%jmin)
          ymax = y(inlet(nflow)%jmax+1)
          if (icyl.eq.1) then
             if (inlet(nflow)%jmin.eq.jmin) then
                ymin = -ymax
                factor = 8.0_WP
             else
                factor = 6.0_WP
             end if
          else
             factor = 6.0_WP
          end if
          
          do j=max(jmino_,inlet(nflow)%jmino),min(jmaxo_,inlet(nflow)%jmaxo)
             Uin(j,:) = factor*ubulk(nflow)*(ym(j)-ymin)*(ymax-ym(j))/(ymax-ymin)**2
             Vin(j,:) = 0.0_WP
             Win(j,:) = factor*wbulk(nflow)*(ym(j)-ymin)*(ymax-ym(j))/(ymax-ymin)**2
          end do
          
       ! Mean turbulent inflow 
       ! Power law function for each inflow
       ! TO BE CHANGED - Only for central pipe in cyl
       case('power')
          ymax = y(inlet(nflow)%jmax+1)
          factor = 60.0_WP/49.0_WP
          
          do j=max(jmino_,inlet(nflow)%jmino),min(jmaxo_,inlet(nflow)%jmaxo)
             Uin(j,:) = factor*ubulk(nflow)*(1.0_WP-ym(j)/ymax)**(1.0_WP/7.0_WP)
             Vin(j,:) = 0.0_WP
             Win(j,:) = 0.0_WP
          end do
          
       ! Tanh + noise profile
       ! Will be dealt with later
       case ('tanh')
          
       ! Profile from a file
       ! Will be dealt with later
       case('file')

       case default
          call die('inflow_read_velocity: unknown velocity inlet type')
       end select
    end do
    
    return
  end subroutine inflow_read_velocity
  
  ! Prepare a tanh velocity profile for injection
  subroutine inflow_tanh_velocity
    use data
    implicit none
    
    integer  :: j,k,nflow
    real(WP) :: Ujet,Ucof
    real(WP) :: radius,theta,r
    real(WP) :: rand,amp
    
    !$OMP PARALLEL PRIVATE(Ujet,Ucof,radius,theta,amp,r,rand)

    do nflow=1,ninlet
       
       if (trim(inlet_type(nflow)).eq.'tanh') then
          
          ! Set parameters
          Ujet=1.0_WP
          Ucof=0.069_WP
          radius=0.5_WP
          theta=0.025_WP
          amp=0.05_WP
          
          ! Base profile
          !$OMP DO
          do k=kmino_,kmaxo_
             do j=jmino_,jmaxo_
                r = sqrt(ym(j)**2+zm(k)**2)
                Uin(j,k) = 0.5_WP*(Ujet+Ucof)-0.5_WP*(Ujet-Ucof)*tanh(0.25_WP*radius/theta*(r/radius-radius/r))
                Vin(j,k) = 0.0_WP
                Win(j,k) = 0.0_WP
             end do
          end do
          !$OMP END DO
          
          ! Add noise
          !$OMP DO
          do k=kmino_,kmaxo_
             do j=jmino_,jmaxo_
                r = sqrt(ym(j)**2+zm(k)**2)
                if (r.lt.0.4_WP) then
                   call random_number(rand); rand=amp*(2.0_WP*rand-1.0_WP)
                   Uin(j,k) = Uin(j,k) + 0.2_WP*rand
                   call random_number(rand); rand=amp*(2.0_WP*rand-1.0_WP)
                   Vin(j,k) = Vin(j,k) + 0.2_WP*rand
                   call random_number(rand); rand=amp*(2.0_WP*rand-1.0_WP)
                   Win(j,k) = Win(j,k) + 0.2_WP*rand
                else if (r.ge.0.4_WP .and. r.le.0.6_WP) then
                   call random_number(rand); rand=amp*(2.0_WP*rand-1.0_WP)
                   Uin(j,k) = Uin(j,k) + 1.0_WP*rand
                   call random_number(rand); rand=amp*(2.0_WP*rand-1.0_WP)
                   Vin(j,k) = Vin(j,k) + 1.0_WP*rand
                   call random_number(rand); rand=amp*(2.0_WP*rand-1.0_WP)
                   Win(j,k) = Win(j,k) + 1.0_WP*rand
                end if
             end do
          end do
          !$OMP END DO
          
       end if
       
    end do

    !$OMP END PARALLEL
    
    return
  end subroutine inflow_tanh_velocity
  
  ! Read the scalar values from the input file
  subroutine inflow_read_scalar
    use parser
    use data
    implicit none

    integer :: isc,j,nflow,n,ninput
    character(len=str_medium), dimension(:), pointer :: list
    character(len=str_medium) :: name
    real(WP) :: val

    ! Return if no scalar
    if (nscalar.eq.0) return
    
    ! Read the input and check it
    allocate(list(nscalar*(ninlet+1)))
    call parser_getsize('Inlet scalar values',ninput)
    if (ninput.ne.(1+ninlet)*nscalar) then
       call die('inflow_init: Wrong number of Inlet scalar values')
    else
       call parser_read('Inlet scalar values',list)
    end if
    
    do n=1,nscalar
       ! Find the right scalar
       loop:do isc=1,nscalar
          name = list((n-1)*(ninlet+1)+1)
          if (trim(SC_name(isc)).eq.trim(name)) exit loop
       end do loop
       if (isc.eq.nscalar+1) &
            call die('inflow_read_scalar: unknown scalar name in input')

       ! Set the values
       do nflow=1,ninlet
          read(list((n-1)*(ninlet+1)+nflow+1),*) val
          do j=max(jmino_,inlet(nflow)%jmino),min(jmaxo_,inlet(nflow)%jmaxo)
             SCin(j,:,isc) = val
          end do
       end do
    end do

    return
  end subroutine inflow_read_scalar
  
end module inflow


! ============================ !
! Initilaize the inflow        !
! -> split the communicator    !
! -> detect the inflows        !
! ============================ !
subroutine inflow_init
  use inflow
  use data
  use parser
  implicit none

  integer :: ninput

  ! If periodic in x => no inflow
  if (xper.eq.1 .or. ninlet.eq.0) return

  ! Read the types of velocity inflow from the input file 
  call parser_getsize('Inlet velocity type',ninput)
  allocate(inlet_type(ninlet))
  if (ninput.eq.ninlet) then
     call parser_read('Inlet velocity type',inlet_type)
  else
     call die('inflow_init: Wrong number of Inlet velocity type defined.')
  end if
  
  ! Allocate the inflow variables
  allocate(Uin(jmino_:jmaxo_,kmino_:kmaxo_));Uin=0.0_WP
  allocate(Vin(jmino_:jmaxo_,kmino_:kmaxo_));Vin=0.0_WP
  allocate(Win(jmino_:jmaxo_,kmino_:kmaxo_));Win=0.0_WP
  allocate(SCin(jmino_:jmaxo_,kmino_:kmaxo_,nscalar));SCin=0.0_WP
  
  ! Read the velocity values from the input file or inflow file
  ! Read the scalar values from the input file
  call inflow_read_velocity
  call inflow_read_scalar
  
  ! Initialize the inflow profile and read the first time step
  call inflow_file_init

  return
end subroutine inflow_init
 

! ======================================== !
! Compute the Inflow for the new time step !
! -> velocity components                   !
! -> momentum components                   !
! ======================================== !
subroutine inflow_velocity
  use inflow
  use velocity
  use data
  implicit none

  integer :: i,j,k,nflow
  
  ! Nothing to do if:
  ! -> not the first proc in x
  ! -> periodic in x
  if (iproc.ne.1 .or. xper.eq.1) return
  
  ! Precompute the inflow profile for a tanh profile
  call inflow_tanh_velocity
  
  ! Precompute the inflow profile from a file
  call inflow_file_velocity
  
  ! Different types of inflows
  do nflow=1,ninlet
     do j=max(jmino_,inlet(nflow)%jmino),min(jmaxo_,inlet(nflow)%jmaxo)
        !$OMP PARALLEL DO
        do k=kmino_,kmaxo_
           ! All left ghost cells
           do i=imino,imin-1
              U(i,j,k)    = Uin(j,k)
              rhoU(i,j,k) = Uin(j,k) * RHOmid(imino,j,k)
              V(i,j,k)    = Vin(j,k)
              rhoV(i,j,k) = Vin(j,k) * RHOmid(imino,j,k)
              W(i,j,k)    = Win(j,k)
              rhoW(i,j,k) = Win(j,k) * RHOmid(imino,j,k)
           end do
           ! Last grid point only for U/rhoU
           U(imin,j,k) = Uin(j,k)
           rhoU(imin,j,k) = Uin(j,k) * RHOmid(imino,j,k)
        end do
        !$OMP END PARALLEL DO
     end do
  end do

  return
end subroutine inflow_velocity


! ======================================== !
! Compute the Inflow for the new time step !
! -> scalars                               !
! ======================================== !
subroutine inflow_scalar
  use inflow
  use velocity
  use data
  implicit none

  integer :: i,j,k,nflow,isc

  ! Nothing to do if:
  ! -> not the first proc in x
  ! -> periodic in x
  if (iproc.ne.1 .or. xper.eq.1) return

  ! Precompute the inflow profile from a file
  call inflow_file_scalar

  ! Different types of inflows
  do nflow=1,ninlet
     do isc=1,nscalar
        do j=max(jmino_,inlet(nflow)%jmino),min(jmaxo_,inlet(nflow)%jmaxo)
           !$OMP PARALLEL DO
           do k=kmino_,kmaxo_
              do i=imino,imin-1
                 SC(i,j,k,isc) = SCin(j,k,isc)
              end do
           end do
           !$OMP END PARALLEL DO
        end do
     end do
  end do

  ! NSCBC for density
  if (compressible) call nscbc_inflow_scalar

  return
end subroutine inflow_scalar
