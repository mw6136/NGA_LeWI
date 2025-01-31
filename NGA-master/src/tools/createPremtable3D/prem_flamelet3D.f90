module prem_flamelet3D
  use precision
  use string
  implicit none
  
  ! Combustion model
  character(len=str_medium) :: combModel
  
  ! Number of points in the flamelet
  integer :: nPoints
  
  ! Coordinate in progress variable space
  real(WP), dimension(:), pointer :: C
  
  ! List of names/variables to get from FlameMaster files
  integer :: nvar_in
  character(len=str_medium), dimension(:), pointer :: input_name
  real(WP), dimension(:,:), pointer :: input_data
  logical, dimension(:), pointer :: found
  
  ! FPVA constituant variables
  integer :: nFPVA
  character(len=str_short), dimension(:), pointer :: FPVA_name
  
  ! Particular arrays for steady flamelets / FPVA
  !real(WP), dimension(:), pointer :: PROG
  real(WP), dimension(:), pointer :: SRC_PROG
  real(WP), dimension(:), pointer :: SRC_PRO2

  ! Particular variable for premixed flamelets
 ! real(WP) :: TOT_ENTHALPY

  ! Particular variables for PFPVA-3D
  real(WP) :: PHI
  real(WP) :: Z
  real(WP) :: Zst = 0.0551538

end module prem_flamelet3D


! ============================================== !
! FLAMELET INITIALIZATION                        !
! Convert the names to be those from FlameMaster !
! ============================================== !
subroutine prem_flamelet3D_init
  use prem_flamelet3D
  use parser
  implicit none
  
  ! Combustion model dependent parameters
  call parser_read('Combustion model',combModel)
  
  ! Get number of additional variables to store in the table
  call parser_getsize("FlameMaster variables", nvar_in)
  
  ! Treat the case the model is FPVA
  select case(trim(combModel))
  case ('PFPVA-3D')
     call parser_getsize('FPVA variables',nFPVA)
     allocate(FPVA_name(nFPVA))
     call parser_read('FPVA variables',FPVA_name)
     
     allocate(input_name(nvar_in+2))
     call parser_read("FlameMaster variables", input_name(1:nvar_in))

     nvar_in = nvar_in + 2
     input_name(nvar_in-1) = 'SRC_PRO2' ! SRC_PRO2 = PROG*SRC_PROG
     input_name(nvar_in) = 'SRC_PROG'
     
  case default
     stop "flamelet_init: Unknown combustion model"
  end select
  
  ! Allocate array to specify wether the variables have been found
  allocate(found(nvar_in))
  
  return
end subroutine prem_flamelet3D_init



! ================================================ !
! Read a flamelet file and store its value in data !
! ================================================ !
subroutine prem_flamelet3D_readfile(filename)
  use prem_flamelet3D
  use fileio
  use parser
  implicit none
  
  character(len=*), intent(in) :: filename
  integer :: iunit, ierr, var, n, nlines, index1
  character(len=str_long) :: buffer
  character(len=20*str_long) :: line
  character(len=str_medium) :: varname, tmpname
  real(WP), dimension(:), pointer :: tmp
  
  ! Open the file
  iunit = iopen()
  open(iunit,file=trim(filename),form='formatted',status='old',iostat=ierr)
  if (ierr.ne.0) then
     print*,"flamelet_readfile: Error opening the file : " // filename
     stop
  end if
  
  nlines = 0
  found = .false.
  ierr = 0
  buffer = ''
  do while(index(buffer,'body').eq.0)
     
     ! First get some parameters
     read(iunit,'(a)',iostat=ierr) buffer
     
     ! Get nPoints and allocate arrays
     if (index(buffer,'gridPoints').ne.0) then
        read(buffer(index(buffer,'=')+1:),*) nPoints
        nlines = ceiling(nPoints / 5.0_WP)
        
        allocate(input_data(nPoints,nvar_in))
        allocate(tmp(nPoints))
        allocate(C(nPoints))
        C = 0.0_WP
        
        if (trim(combModel).eq.'PFPVA-3D') then
           SRC_PROG => input_data(:,nvar_in)
           SRC_PRO2 => input_data(:,nvar_in-1)
           SRC_PROG = 0.0_WP
           SRC_PRO2 = 0.0_WP
        end if
        
     end if

     ! If prem flame, read in phi
     if (trim(combModel).eq.'PFPVA-3D') then
        if (index(buffer,'fuel-air-equivalence-ratio').ne.0) then
           read(buffer(index(buffer,'=')+1:),*) PHI
           call parser_read('Stoichiometric mixture fraction',Zst)
           Z = PHI/((1-Zst)/Zst + PHI)
        end if
     end if
     
  end do
     
  
  ! Test
  if (nlines.eq.0) stop "flamelet_readfile: missing gridPoints in flamemet file"
  
  ! Preset diffusivity to 1
  loop0:do var=1,nvar_in
     if (trim(input_name(var)).eq.'diffusivity') exit loop0
  end do loop0
  if (var.le.nvar_in) input_data(:,var) = 1.0_WP

  loop1:do while (ierr .eq. 0)
     
     ! Read name of variable
     read(iunit,'(a)',iostat=ierr) buffer
     if (trim(buffer).eq.'trailer') exit loop1
     
     ! Read name of variable
     read(buffer,'(a)',iostat=ierr) varname
     index1 = index(varname,' ')
     if (index1.ne.0) varname(index1:) = ''
     
     ! Read the array
     line = ''
     do n=1,nlines
        read(iunit,'(a)',iostat=ierr) buffer
        line = trim(line) // adjustl(trim(buffer))
     end do

     ! Is it the coordinate Z?
     if (trim(varname).eq.'Z') then
        read(line,*) Z
     end if
     
     ! Is it part of the diffusivity?
     if (trim(varname).eq.'lambda') then
        read(line,*) tmp
        loop4:do var=1,nvar_in
           if (trim(input_name(var)).eq.'diffusivity') exit loop4
        end do loop4
        if (var.le.nvar_in) then
           input_data(:,var) = input_data(:,var) * tmp
           found(var) = .true.
        end if
     end if
     if (trim(varname).eq.'cp') then
        read(line,*) tmp
        loop5:do var=1,nvar_in
           if (trim(input_name(var)).eq.'diffusivity') exit loop5
        end do loop5
        if (var.le.nvar_in) then
           input_data(:,var) = input_data(:,var) / tmp
           found(var) = .true.
        end if
     end if
     
     if (trim(varname).eq.'lambdaOverCp') then
        ! WARNING - this is a non-dimensionalized diffusivity
        read(line,*) tmp
        loop9:do var=1,nvar_in
           if (trim(input_name(var)).eq.'diffusivity') exit loop9
        end do loop9
        if (var.le.nvar_in) then
           input_data(:,var) = input_data(:,var) / tmp
           found(var) = .true.
        end if
     end if

     ! Is it a variable for FPVA ?
     if (trim(combModel).eq.'PFPVA-3D') then
        loop2:do var=1,nFPVA
           tmpname = 'massfraction-' // trim(FPVA_name(var))
           if (trim(varname).eq.trim(tmpname)) then
              read(line,*) tmp
              C = C + tmp
              !found(nvar_in) = .true.
              exit loop2
           end if
           tmpname = 'ProdRate' // trim(FPVA_name(var))
           if (trim(varname).eq.trim(tmpname)) then
              read(line,*) tmp
              SRC_PROG = SRC_PROG + tmp
              found(nvar_in) = .true.
              exit loop2
           end if
        end do loop2
     end if
     
     ! Do we want that variable?
     loop3:do var=1,nvar_in
        if (trim(input_name(var)).eq.varname) then
           read(line,*) input_data(:,var)
           found(var) = .true.
           exit loop3
        end if
     end do loop3
  end do loop1
  
  ! For FPVA, divide source term by density
  if (trim(combModel).eq.'PFPVA-3D') then
     loop6:do var=1,nvar_in
        if (trim(input_name(var)).eq.'density') exit loop6
     end do loop6
     if (var.le.nvar_in) SRC_PROG = SRC_PROG / input_data(:,var)
     ! Calculate source term for PROG^2 equation
     if (found(nvar_in)) then
        SRC_PRO2 = C * SRC_PROG
        found(nvar_in -1) = .true.
     end if
  end if
  
  if (trim(combModel).eq.'Enthalpy Flamelet') then
     do var=1,nvar_in-1
        if (.not.found(var)) then
           print*,"Variable",trim(input_name(var))," not found in flamelet file"
           stop
        end if
     end do
  else
     do var=1,nvar_in
        if (.not.found(var)) then
           print*,"Variable ",trim(input_name(var))," not found in flamelet file"
           stop
        end if
     end do
  end if
  
  ! Deallocate
  deallocate(tmp)
  nullify(tmp)
  close(iclose(iunit))
  
  return
end subroutine prem_flamelet3D_readfile


! ========================================= !
! Deallocate the data array for new reading !
! ========================================= !
subroutine prem_flamelet3D_cleanup
  use prem_flamelet3D
  implicit none
  
  deallocate(input_data)
  !deallocate(Z)
  nullify(input_data)
  !nullify(Z)
  
  return
end subroutine prem_flamelet3D_cleanup
