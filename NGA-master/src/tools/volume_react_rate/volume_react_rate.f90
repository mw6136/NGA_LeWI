! ============================================================ !
!                volume_react_rate.f90                         !
!                                                              !
! Author: Temistocle Grenga                                    !
! Date:   June, 2017                                           !
! ============================================================ !

program volume_react_rate
  use parallel
  use string
  use precision
!  use geometry
  use partition
  use masks
  use fileio
!  use data
!  use masks

  implicit none
  character(len=str_medium) :: filename1, filename_base, folder
  integer :: i, j, k, l, ll, n
  real(WP) :: w_snap
  integer :: iunit1, ierr, iunit2
  integer :: n_snap, start_file, step_file, nvars
  integer, dimension(:), pointer :: var_list
  integer :: idum, nvars1, nx1, ny1, nz1
  real(WP), dimension(:), pointer :: x1, y1, z1
  real(WP), dimension(:,:,:,:), pointer :: data8, prod_rate
  real(SP), dimension(:,:,:),   pointer :: data4
  character(len=str_short), dimension(:), pointer :: names
  real(WP) :: ddum
  real(WP), dimension(:,:,:), pointer :: datadum
  real(WP), dimension(:), pointer :: W_sp, h_sp, Cp_sp, conc, conc_dot, sol, sol_t
  real(WP), dimension(:), pointer :: K_rxn, omega_rxn, M_rxn
  real(WP) :: P_thermo, W_mix, Cp_mix, rho_gp, tmp
  integer :: N_nons, Nreactions
  real(WP), parameter :: R_cst =8314.34_WP  ! [J/(kmol.K)]
  character(len=80) :: buffer

  ! Initialize MPI enironment
  call parallel_init()

  start_file = 8830   !6465 High !8830 Low  !15687 NR  !21525 React
  step_file = 1
  n_snap = 400
  nvars = 11
  N_nons = 9
  P_thermo = 1.01325e5

  allocate(var_list(nvars))
  ! 4=P, 6=N2, 7=H, 8=O2, 9=O, 10=OH, 11=H2, 12=H2O, 13=HO2, 14=H2O2, 15=T
  var_list = (/ 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 /)

  filename_base = 'vol_data.1'

  if (nproc /= n_snap+1) then
     if (irank == 1) print*, 'nproc .neq. nsnap+1'
        call parallel_final()
        GO TO 110
  end if

  ! 1 snapshot per core
  w_snap = MOD(irank,n_snap+1)

  l = start_file +step_file*w_snap
  write(filename1, '(''_'', I8.8)') l
  filename1 = trim(filename_base)//trim(filename1)
  filename1 = "volume_data/"//trim(filename1)
  print*, filename1

  !Read scalars
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

  call BINARY_FILE_READ(iunit1, idum, 1, kind(idum), ierr)
  call BINARY_FILE_READ(iunit1, nx1, 1, kind(idum), ierr)
  call BINARY_FILE_READ(iunit1, ny1, 1, kind(idum), ierr)
  call BINARY_FILE_READ(iunit1, nz1, 1, kind(idum), ierr)
  call BINARY_FILE_READ(iunit1, nvars1, 1, kind(nvars1), ierr)
  allocate (x1(nx1), y1(ny1), z1(nz1))
  if (irank == 1) print*, nx1, ny1, nz1
  call BINARY_FILE_READ(iunit1, x1, nx1, kind(x1), ierr)
  call BINARY_FILE_READ(iunit1, y1, ny1, kind(y1), ierr)
  call BINARY_FILE_READ(iunit1, z1, nz1, kind(z1), ierr)
  deallocate(x1,y1,z1)
  allocate(names(nvars1))
  do i = 1,nvars1
     call BINARY_FILE_READ(iunit1, names(i), str_short, kind(names), ierr)
  end do
  do i = 1, 2
     call BINARY_FILE_READ(iunit1, ddum, 1, kind(ddum), ierr)
  end do

  allocate(data8(nvars,nx1,ny1,nz1), datadum(nx1,ny1,nz1))
  do ll = 1, var_list(1)-1
     call BINARY_FILE_READ(iunit1, datadum, nx1*ny1*nz1, kind(datadum), ierr)
  end do

  data8(:,:,:,:) = 0.0_WP
  do l = 1,nvars
     print*, 'Reading', var_list(l)
     call BINARY_FILE_READ(iunit1, data8(l,:,:,:), nx1*ny1*nz1, kind(data8), ierr)
     if ( l == nvars)  GO TO 212
     ll = var_list(l+1) - var_list(l) -1
     do j = 1,ll
        call BINARY_FILE_READ(iunit1, datadum(:,:,:), nx1*ny1*nz1, kind(datadum), ierr)
     end do
  end do
 212 CONTINUE

  call BINARY_FILE_CLOSE(iunit1,ierr)
  deallocate(datadum)
  print*, 'closed', filename1


  ! Calculate Reaction Rate
  ! Variables order in data8: 1=P, 2=N2, 3=H, 4=O2, 5=O, 6=OH, 7=H2, 8=H2O, 9=HO2, 10=H2O2, 11=T
  allocate(prod_rate(N_nons+1,nx1,ny1,nz1))
  allocate(W_sp(1:N_nons), h_sp(1:N_nons), Cp_sp(1:N_nons))
  allocate(sol(1:N_nons+1), sol_t(1:N_nons), conc(1:N_nons), conc_dot(1:N_nons))
  call GETMOLARMASS(W_sp)
  call GETNREACTIONS(Nreactions)
  allocate(K_rxn(1:Nreactions), omega_rxn(1:Nreactions), M_rxn(1:Nreactions))

  print*,'Evaluating RR'

  do k=1,nz1
     do j=1,ny1
        do i=1,nx1
           ! Get species mass fractions
           sol(1:N_nons+1) = data8(2:11,i,j,k)

           ! Get mixture-averaged W and Cp
           ! from finitechem_W
           sol_t = min(max(sol(1:N_nons),0.0_WP),1.0_WP)
           W_mix = sum(sol_t)/sum(sol_t/W_sp)
           ! from finitechem_Cp
           call COMPTHERMODATA(h_sp, Cp_sp, sol(N_nons+1))
           Cp_mix = sum(sol(1:N_nons)*Cp_sp)
           rho_gp = P_thermo*W_mix/(R_cst*sol(N_nons+1))

           ! Compute species concentrations
           do n=1,N_nons
              conc(n) = rho_gp*sol(n)/W_sp(n)
           end do

           ! Obtain source terms
           call PRODRATES(conc_dot, omega_rxn, K_rxn, conc, M_rxn, sol(N_nons+1), P_thermo)

           ! Convert to mass fraction source terms
           do n=1,N_nons
              prod_rate(n,i,j,k) = conc_dot(n)*W_sp(n)
           end do

           ! Source term for temperature
           tmp = 0.0_WP
           do n=1,N_nons
              tmp = tmp - h_sp(n)*prod_rate(n,i,j,k)
           end do
           prod_rate(N_nons+1,i,j,k) = tmp/Cp_mix
        end do
     end do
  end do

  print*,'Printing RR'
  ! Print Reaction Rate in Ensight format
  allocate(data4(nx1,ny1,nz1))
  folder = 'react_rate_snap/'
  do l = 2,nvars
     filename1 = trim(adjustl(names(var_list(l)))) // "_"
     i = start_file+w_snap
     write(filename1(len_trim(filename1)+1:len_trim(filename1)+6),'(i6.6)') i
     print*, irank, filename1
     filename1 = trim(folder)//trim(filename1)
     call BINARY_FILE_OPEN(iunit2,trim(filename1),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
     buffer = 'part'
     call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
     idum = 1
     call BINARY_FILE_WRITE(iunit2,idum,1,kind(idum),ierr)
     buffer = 'block'
     call BINARY_FILE_WRITE(iunit2,buffer,80,kind(buffer),ierr)
     do k=1,nz1
        do j=1,ny1
           do i=1,nx1
              data4(i,j,k) = prod_rate(l-1,i,j,k)
           end do
        end do
     end do
     call BINARY_FILE_WRITE(iunit2,data4,nx1*ny1*nz1,kind(data4),ierr)
     call BINARY_FILE_CLOSE(iunit2,ierr)
  end do
  print*, irank, filename1, 'written'

  deallocate(data4)

  call parallel_final()

110 print*, irank, 'Complete'

end program volume_react_rate

