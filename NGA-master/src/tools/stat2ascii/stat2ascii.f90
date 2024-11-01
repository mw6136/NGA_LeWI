module stat2ascii
  use precision
  use string
  use fileio
  use parser
  use cli_reader
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: nx,ny,nz,nvar
  character(len=str_medium), dimension(:), pointer :: names
  integer :: var
  character(len=str_medium) :: fdata,fconfig,config,input_name,output_name
  character(len=str_short) :: of_name
  real(WP) :: dt,time
  
  ! Sizes
  integer :: xper,yper,zper
  integer :: icyl

  ! Mesh
  real(WP), dimension(:), pointer :: x,xm
  real(WP), dimension(:), pointer :: y,ym
  real(WP), dimension(:), pointer :: z,zm
  integer, dimension(:,:), pointer :: mask

  ! Planes
  real(WP), dimension(:), pointer :: pnloc,pnlocm
  integer,  dimension(:), pointer :: pnindx,pnindxm
  real(WP), dimension(:,:), pointer :: interp_pnx,interp_pnxm

  ! Data
  integer :: isc_RHO,isc_PROG,isc_rhoU,isc_rhoU2,isc_H2O
  real(WP), dimension(:,:,:), pointer :: data
  real(WP), dimension(:,:), pointer :: y_condm,RHO,U,V,W,U2,V2,W2

contains

  ! ========================================================== !
  ! Get the filename header                                    !
  ! ========================================================== !
  subroutine get_name(iplane,of_name)
    implicit none
    integer, intent(in) :: iplane
    character(len=str_short) :: of_name

    if (iplane<10) then
       write(of_name,'(A1,I1)') '0', iplane
    else
       write(of_name,'(I2)') iplane
    end if

    return
  end subroutine get_name

end module stat2ascii


program main
  use stat2ascii
  implicit none

  integer :: i,j,k
  integer :: iunit,ierr
  integer :: pnmin,pnmax,iplane
  real(WP) :: prog_upper,prog_lower
  real(WP) :: tmp1

  ! Parse the input file
  call get_command_argument(1,input_name)
  call parser_init
  call parser_parsefile(input_name)

  ! Get file names
  call parser_read('Stat data file',fdata)
  call parser_read('Config file',fconfig)
  call parser_read('Output name',output_name)

  call system("mkdir -p output_ascii")
  
  ! Read the source config file
  call BINARY_FILE_OPEN(iunit,trim(fconfig),"r",ierr)
  call BINARY_FILE_READ(iunit,config,str_medium,kind(config),ierr)
  call BINARY_FILE_READ(iunit,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit,xper,1,kind(xper),ierr)
  call BINARY_FILE_READ(iunit,yper,1,kind(yper),ierr)
  call BINARY_FILE_READ(iunit,zper,1,kind(zper),ierr)
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)
  print*,'Source file'
  print*,'Config : ',config
  print*,'Grid :',nx,'x',ny,'x',nz
  allocate(x(nx+1),y(ny+1),z(nz+1))
  allocate(mask(nx,ny))
  call BINARY_FILE_READ(iunit,x,nx+1,kind(x),ierr)
  call BINARY_FILE_READ(iunit,y,ny+1,kind(y),ierr)
  call BINARY_FILE_READ(iunit,z,nz+1,kind(z),ierr)
  call BINARY_FILE_READ(iunit,mask,nx*ny,kind(mask),ierr)
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  ! Finish to create the meshes
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

  ! Planes
  ! Get plane locations
  pnmin = 1
  call parser_getsize('Plane x locations',pnmax)
  allocate(pnloc (pnmax))
  allocate(pnlocm(pnmax))
  call parser_read('Plane x locations',pnloc)
  allocate(pnindx (pnmax))
  allocate(pnindxm(pnmax))
  allocate(interp_pnx (pnmax,2))
  allocate(interp_pnxm(pnmax,2))

  ! Compute the interpolation to y-z planes
  do iplane=pnmin,pnmax
     if (pnloc(iplane).gt.xm(nx-1)) stop 'plane x-index not in domain'
     ! Interpolation from midpoints
     i = 1
     do while (xm(i).lt.pnloc(iplane))
        i = i+1
     end do
     interp_pnxm(iplane,1) = (xm(i)-pnloc(iplane))/(xm(i)-xm(i-1))
     interp_pnxm(iplane,2) = 1.0_WP-interp_pnxm(iplane,1)
     pnindxm(iplane) = i
     pnlocm (iplane) = interp_pnxm(iplane,1)*xm(i-1) + interp_pnxm(iplane,2)*xm(i)

     ! Interpolation from x-faces
     i = 1
     do while (x(i).lt.pnloc(iplane))
        i = i+1
     end do
     interp_pnx(iplane,1) = (x(i)-pnloc(iplane))/(x(i)-x(i-1))
     interp_pnx(iplane,2) = 1.0_WP-interp_pnx(iplane,1)
     pnindx(iplane) = i
     pnloc (iplane) = interp_pnx(iplane,1)*x(i-1) + interp_pnx(iplane,2)*x(i)
  end do

  print '(10i15)', pnindx(1:10)
  print '(10ES15.5)', pnloc(1:10)
  print '(10i15)', pnindxm(1:10)
  print '(10ES15.5)', pnloc(1:10)

  
  ! ** Open the stat file to read **
  call BINARY_FILE_OPEN(iunit,trim(fdata),"r",ierr)
  
  ! Read sizes
  call BINARY_FILE_READ(iunit,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit,nvar,1,kind(nvar),ierr)
  print*,'Grid :',nx,'x',ny,'x',nz
  print*,'nvar :',nvar
  if (nz.ne.1) stop "stat should be 2D"
  
  ! Read additional stuff
  call BINARY_FILE_READ(iunit,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit,time,1,kind(time),ierr)
  print*,'Data file at time :',time
  
  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit,names(var),str_medium,kind(names),ierr)
  end do
  !print*,'Variables : ',names
  
  ! Allocate arrays
  allocate(data(nx,ny,nvar))
  allocate(RHO(nx,ny))
  allocate(U(nx,ny))
  allocate(V(nx,ny))
  allocate(W(nx,ny))
  allocate(U2(nx,ny))
  allocate(V2(nx,ny))
  allocate(W2(nx,ny))
  
  do var=1,nvar
     ! Read data field
     call BINARY_FILE_READ(iunit,data(:,:,var),nx*ny,kind(data),ierr)
     !print*,maxval(abs(data8)), ' at ', maxloc(abs(data8))

     ! Find the progress variable
     select case(trim(names(var)))
     case ('rhoS-O2')
        isc_PROG = var
     case ('RHO')
        RHO = data(:,:,var)
        isc_RHO = var
     case('rhoU')
        isc_rhoU = var
        U = data(:,:,var)
     case('rhoV')
        V = data(:,:,var)
     case('rhoW')
        W = data(:,:,var)
     case('rhoU^2')
        isc_rhoU2 = var
        U2 = data(:,:,var)
     case('rhoV^2')
        V2 = data(:,:,var)
     case('rhoW^2')
        W2 = data(:,:,var)
     case('rhoS-H2O')
        isc_H2O = var
     case default
     end select
  end do
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Compute the density weighting
  U = U/RHO
  V = V/RHO
  W = W/RHO
  U2 = U2/RHO
  V2 = V2/RHO
  W2 = W2/RHO

  ! Compute the progress variable
  call parser_read('Prog upper', prog_upper)
  call parser_read('Prog lower', prog_lower)
  allocate(y_condm(nx,ny))
  do j=1,ny
     do i=1,nx
        ! Define the progress variable based on O2
        tmp1 = (prog_upper - data(i,j,isc_PROG)/RHO(i,j))/(prog_upper-prog_lower)
        y_condm(i,j) = tmp1
     end do
  end do
 

  ! Output data at centerline


  ! Output data at the planes
  do iplane=pnmin,pnmax
     call get_name(iplane,of_name)
     open(unit=iunit, file='output_ascii/'//trim(output_name)//'_RHO_'//trim(of_name), action='write')
     i = pnindxm(iplane)
     do j=1,ny
        write (iunit,'(ES22.13)',advance='no'), y(j) ! y-coordinate
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*y_condm(i-1:i,j))
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*data(i-1:i,j,isc_RHO))
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnx (iplane,:)*U(i-1:i,j))   ! U
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*V(i-1:i,j)) ! V
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*W(i-1:i,j)) ! W
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*data(i-1:i,j,isc_H2O)/data(i-1:i,j,isc_RHO)) ! H2O
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnx (iplane,:)*(U2(i-1:i,j)-U(i-1:i,j)*U(i-1:i,j)))   ! U_rms
        write (iunit,'(ES22.13)',advance='no'), sum(interp_pnxm(iplane,:)*(V2(i-1:i,j)-V(i-1:i,j)*V(i-1:i,j)))   ! V_rms
        write (iunit,'(ES22.13)'), sum(interp_pnxm(iplane,:)*(W2(i-1:i,j)-W(i-1:i,j)*W(i-1:i,j)))   ! W_rms
     end do
     close(iunit)
  end do

end program main
