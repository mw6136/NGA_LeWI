module flamelet
  use precision
  use string
  implicit none
  
  ! Combustiopn model
  character(len=str_medium) :: combModel
  
  ! Number of points in the flamelet
  integer :: nPoints
  
  ! Coordinate in mixture fraction space
  real(WP), dimension(:), pointer :: Z

  !!! Parameters for multiple mixture fraction models
  real(WP) :: FMIX
  
  ! List of names/variables to get from FlameMaster files
  integer :: nvar_in
  character(len=str_medium), dimension(:), pointer :: input_name
  real(WP), dimension(:,:), pointer :: input_data
  logical, dimension(:), pointer :: found
  
  ! FPVA constituant variables
  logical :: useFPVA = .false. ! true for single or multiple mixture fraction FPVA
  integer :: nFPVA
  character(len=str_short), dimension(:), pointer :: FPVA_name

  ! Modified progress variable definition
  logical :: modprog

  ! Mixture fraction variance transport equation
  logical :: vartran
  
  ! Particular arrays for steady flamelets / FPVA
  real(WP), dimension(:), pointer :: PROG
  real(WP), dimension(:), pointer :: SRC_PROG

  ! Particular arrays for variance transport equation
  real(WP), dimension(:), pointer :: Z2RHODOT
  real(WP), dimension(:), pointer :: ZSRC_ZMIX
  real(WP), dimension(:), pointer :: ZS2RHODOT
  real(WP), dimension(:), pointer :: ZSRC_ZS

end module flamelet


! ============================================== !
! FLAMELET INITIALIZATION                        !
! Convert the names to be those from FlameMaster !
! ============================================== !
subroutine flamelet_init
  use flamelet
  use parser
  implicit none
  
  ! Combustion model dependent parameters
  call parser_read('Combustion model',combModel)
  
  ! Get number of additionnal variables to store in the table
  call parser_getsize("FlameMaster variables", nvar_in)
  
  ! Treat the case the model is FPVA
  select case(trim(combModel))
  case ('FPVA')
     call parser_getsize('FPVA variables',nFPVA)
     allocate(FPVA_name(nFPVA))
     call parser_read('FPVA variables',FPVA_name)

     call parser_read('Modified progress variable',modprog,.false.)

     call parser_read('Solve variance transport equation',vartran,.false.)

     if (vartran) then
        allocate(input_name(nvar_in+4))
     else
        allocate(input_name(nvar_in+2))
     end if

     call parser_read("FlameMaster variables", input_name(1:nvar_in))

     if (vartran) then
        nvar_in = nvar_in + 4
        input_name(nvar_in-3) = 'Z2RHODOT'
        input_name(nvar_in-2) = 'ZSRC_ZMIX'
     else
        nvar_in = nvar_in + 2
     end if     
     input_name(nvar_in-1) = 'SRC_PROG'
     input_name(nvar_in)   = 'PROG'
     useFPVA = .true.

  ! Treat either multiple mixture fraction case
  case ('MMFPVA-F', 'MMFPVA-Z*')
     call parser_getsize('FPVA variables',nFPVA)
     allocate(FPVA_name(nFPVA))
     call parser_read('FPVA variables',FPVA_name)

     call parser_read('Modified progress variable',modprog,.false.)

     call parser_read('Solve variance transport equation',vartran,.false.)

     if (vartran) then
        allocate(input_name(nvar_in+7))
     else
        allocate(input_name(nvar_in+3))
     end if
     
     call parser_read("FlameMaster variables", input_name(1:nvar_in))

     if (vartran) then
        nvar_in = nvar_in + 7
        input_name(nvar_in-6) = 'Z2RHODOT'
        input_name(nvar_in-5) = 'ZSRC_ZMIX'
        input_name(nvar_in-4) = 'ZS2RHODOT'
        input_name(nvar_in-3) = 'ZSRC_ZS'
     else
        nvar_in = nvar_in + 3
     end if
     input_name(nvar_in-2) = 'FMIX'
     input_name(nvar_in-1) = 'SRC_PROG'
     input_name(nvar_in)   = 'PROG'
     useFPVA = .true.
     
  case ('Steady Flamelet')
     allocate(input_name(nvar_in+1))
     call parser_read("FlameMaster variables", input_name(1:nvar_in))
     
     nvar_in = nvar_in + 1
     input_name(nvar_in)   = 'chi'
     
  case default
     stop "flamelet_init: Unknown combustion model"
  end select
  
  ! Allocate array to specify wether the variables have been found
  allocate(found(nvar_in))
  
  return
end subroutine flamelet_init



! ================================================ !
! Read a flamelet file and store its value in data !
! ================================================ !
subroutine flamelet_readfile(filename)
  use flamelet
  use fileio
  implicit none
  
  character(len=*), intent(in) :: filename
  integer :: iunit, ierr, var, n, nlines, index1, index_zmixsrc, index_dimer
  integer :: index_rhodot, index_pahsrc_pos, index_pahsrc_neg
  character(len=str_long) :: buffer
  character(len=4*str_long) :: line
  character(len=str_medium) :: varname, tmpname
  real(WP), dimension(:), pointer :: tmp
  character(len=str_medium) :: FuelName
  logical :: infs
  
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
  FMIX= -1.0_WP
  FuelName='NotInit'
  infs= .false.
  
  do while(index(buffer,'body').eq.0)
     
     ! First get some parameters
     read(iunit,'(a)',iostat=ierr) buffer
     
     ! Get nPoints and allocate arrays
     if (index(buffer,'gridPoints').ne.0) then
        read(buffer(index(buffer,'=')+1:),*) nPoints
        nlines = ceiling(nPoints / 5.0_WP)
        
        allocate(input_data(nPoints,nvar_in))
        allocate(tmp(nPoints))
        allocate(Z(nPoints))
        
        if (useFPVA) then
           if (vartran .and. trim(combModel).eq.'FPVA') then
              Z2RHODOT => input_data(:,nvar_in-3)
              ZSRC_ZMIX => input_data(:,nvar_in-2)
           else if (vartran .and. (trim(combModel).eq.'MMFPVA-F' .or.  trim(combModel).eq.'MMFPVA-Z*')) then
              Z2RHODOT => input_data(:,nvar_in-6)
              ZSRC_ZMIX => input_data(:,nvar_in-5)
              ZS2RHODOT => input_data(:,nvar_in-4)
              ZSRC_ZS => input_data(:,nvar_in-3)
           end if
           SRC_PROG => input_data(:,nvar_in-1)
           PROG => input_data(:,nvar_in)
           PROG = 0.0_WP
           SRC_PROG = 0.0_WP
        end if
     end if

     !!! If using multiple mixture fractions, get value of FMIX
     if ((trim(combModel).eq.'MMFPVA-F') .or. (trim(combModel).eq.'MMFPVA-Z*')) then
        if (index(buffer,'fuel = ').ne.0) then
           read(buffer(index(buffer,'=')+1:),*) FuelName
           FuelName = ('Massfraction-'//trim(FuelName)//' = ')
        else if (index(buffer,'FuelSide').ne.0) then
           infs = .true.
        end if

        if (infs .eq. .true.) then
           if (index(buffer,trim(FuelName)).ne.0) then
              read(buffer(index(buffer,'=')+1:),*) FMIX
              found(nvar_in-2) = .true.
           else if (index(buffer,'end').ne.0) then
              infs = .false.
           end if
        end if
     end if
  end do
  
  ! Test
  if (nlines.eq.0) stop "flamelet_readfile: missing gridPoints in flamelet file"
  if ( ((trim(combModel).eq.'MMFPVA-F') .or. (trim(combModel).eq.'MMFPVA-Z*')) .and. (FMIX.lt.0)) &
       stop "flamelet_readfile: Using multiple mixture fractions, could not determine valid FMIX in flamelet file"

  
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
     if (useFPVA) then
        if (.not.modprog) then
           loop2:do var=1,nFPVA
              tmpname = 'massfraction-' // trim(FPVA_name(var))
              if (trim(varname).eq.trim(tmpname)) then
                 read(line,*) tmp
                 PROG = PROG + tmp
                 found(nvar_in) = .true.
                 exit loop2
              end if
              tmpname = 'ProdRate' // trim(FPVA_name(var))
              if (trim(varname).eq.trim(tmpname)) then
                 read(line,*) tmp
                 SRC_PROG = SRC_PROG + tmp
                 found(nvar_in-1) = .true.
                 exit loop2
              end if
           end do loop2
        else
           if (trim(varname).eq.'ProgRat') then
              read(line,*) tmp
              PROG = tmp
              found(nvar_in) = .true.
           end if
           if (trim(varname).eq.'ProgSrc') then
              read(line,*) tmp
              SRC_PROG = tmp
              found(nvar_in-1) = .true.
           end if
        end if
     end if
     
     ! Do we want that variable?
     loop3:do var=1,nvar_in
        if (trim(input_name(var)).eq.varname) then
           read(line,*) input_data(:,var)
           found(var) = .true.

           if (trim(varname).eq.'ZBilgerSrc') then
              index_zmixsrc = var
           end if

           if (trim(varname).eq.'Dimer_ProdRate') then
              index_dimer = var
           end if

           if (trim(varname).eq.'RhoDot') then
              index_rhodot = var
           end if

           if (trim(varname).eq.'ProdRatePos-PAH') then
              index_pahsrc_pos = var
           end if

           if (trim(varname).eq.'ProdRateNeg-PAH') then
              index_pahsrc_neg = var
           end if

           exit loop3
        end if
     end do loop3

     if (vartran .and. trim(combModel).eq.'FPVA') then
        found(nvar_in-3) = .true.
        found(nvar_in-2) = .true.
     else if (vartran .and. (trim(combModel).eq.'MMFPVA-F' .or.  trim(combModel).eq.'MMFPVA-Z*')) then
        found(nvar_in-6) = .true.
        found(nvar_in-5) = .true.
        found(nvar_in-4) = .true.
        found(nvar_in-3) = .true.
     end if
     
  end do loop1

  ! Check if all variables were found
  do var=1,nvar_in
     if (.not.found(var)) then
        print*,"Variable ",trim(input_name(var))," not found in flamelet file"
        stop
     end if
  end do
  
  ! For FPVA, divide source term by density
  if (useFPVA) then
     loop6:do var=1,nvar_in
        if (trim(input_name(var)).eq.'density') exit loop6
     end do loop6
     if (var.le.nvar_in) SRC_PROG = SRC_PROG / input_data(:,var)
  end if

  ! If using multiple mixture fractions, copy second mixture fraction coordinate into input_data
  if ((trim(combModel).eq.'MMFPVA-F') .or. (trim(combModel).eq.'MMFPVA-Z*')) then
     input_data(:,nvar_in-2) = FMIX  
  end if

  ! For mixture fraction source term, divide by density
  if (vartran .and. found(index_zmixsrc)) then
     loop7:do var=1,nvar_in
        if (trim(input_name(var)).eq.'density') exit loop7
     end do loop7
     if (var.le.nvar_in) input_data(:,index_zmixsrc) = input_data(:,index_zmixsrc) / input_data(:,var)
  end if

  ! Mixture fraction source term for mixture fraction squared equation
  if (vartran .and. found(index_zmixsrc)) then
     ZSRC_ZMIX = input_data(:,index_zmixsrc) * Z
     ZSRC_ZS = input_data(:,index_zmixsrc) * Z * FMIX
  end if

  ! For density source term, divide by density, REMULTIPLY AFTER CONVOLUTION
  if (vartran .and. found(index_rhodot)) then
     loop8:do var=1,nvar_in
        if (trim(input_name(var)).eq.'density') exit loop8
     end do loop8
     if (var.le.nvar_in) input_data(:,index_rhodot) = input_data(:,index_rhodot) / input_data(:,var)
  end if

  ! Density source term for mixture fraction squared equation
  if (vartran .and. found(index_rhodot)) then
     Z2RHODOT = input_data(:,index_rhodot) * Z**2
     ZS2RHODOT = input_data(:,index_rhodot) * (Z * FMIX)**2
  end if

  ! For dimer production rate, divide by density, REMULTIPLY AFTER CONVOLUTION
  if (vartran .and. found(index_dimer)) then
     loop10:do var=1,nvar_in
        if (trim(input_name(var)).eq.'density') exit loop10
     end do loop10
     if (var.le.nvar_in) input_data(:,index_dimer) = input_data(:,index_dimer) / input_data(:,var)
  end if

  ! For PAH source terms, divide by density, REMULTIPLY AFTER CONVOLUTION
  if (vartran .and. found(index_pahsrc_pos) .and. found(index_pahsrc_neg)) then
     loop11:do var=1,nvar_in
        if (trim(input_name(var)).eq.'density') exit loop11
     end do loop11
     if (var.le.nvar_in) input_data(:,index_pahsrc_pos) = input_data(:,index_pahsrc_pos) / input_data(:,var)
     if (var.le.nvar_in) input_data(:,index_pahsrc_neg) = input_data(:,index_pahsrc_neg) / input_data(:,var)
  elseif (vartran .and. found(index_pahsrc_pos) .or. vartran .and. found(index_pahsrc_neg)) then
     write(*,*) 'Both positive and negative source terms needed for PAH transport equation'
  end if

  ! For dimer production rate, convert from kmol to mol
  if (vartran .and. found(index_dimer)) input_data(:,index_dimer) = input_data(:,index_dimer) * 1000.0_WP
  
  ! Force 0 at Z=1 for chi
  if (trim(combModel).eq.'Steady Flamelet') then
     input_data(nPoints,nvar_in) = 0.0_WP
  end if

  ! Deallocate
  deallocate(tmp)
  nullify(tmp)
  close(iclose(iunit))
  
  return
end subroutine flamelet_readfile


! ========================================= !
! Deallocate the data array for new reading !
! ========================================= !
subroutine flamelet_cleanup
  use flamelet
  implicit none
  
  deallocate(input_data)
  deallocate(Z)
  nullify(input_data)
  nullify(Z)
  
  return
end subroutine flamelet_cleanup
