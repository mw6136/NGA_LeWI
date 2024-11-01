module prem_table3D
  use prem_flamelet3D
  use precision
  use string
  implicit none
  
  ! Table parameters
  integer :: nC, nCVar, nZ, np_pdf
  real(WP), dimension(:), pointer :: Cdex, CVar, Zdex, Ztmp
  
  ! Beta pdf
  real(WP), dimension(:), pointer :: pdf, pdfaxis
  
  ! List of Flamelet files
  integer :: nfiles
  character(len=str_long), dimension(:), pointer :: files
  real(WP), dimension(:), pointer :: filesCmax
  
  ! Variables to be wrtitten into the table
  integer :: nvar_out
  character(len=str_medium), dimension(:), pointer :: output_name
  real(WP), dimension(:,:,:,:), pointer :: output_data
  integer, dimension(:,:,:), pointer :: mask
  
  ! Variable after convolution
  real(WP), dimension(:,:,:,:), pointer :: postconv
  
contains
  
  ! uniform mesh
  subroutine create_prog
    use parser
    implicit none
    
    integer :: cc, icv
    real(WP) :: Cmax

    call parser_read("Maximum Progress Variable", Cmax)
        
    do cc=1,nC
       Cdex(cc) = Cmax* real(cc-1,WP)/real(nC-1,WP)
    end do

    ! Quadratic scaling between 0 and max possible given Cmax
    ! max variance = delta functions at 0, Cmax: (Cmax/2)^2
    do icv=1,nCVar
       CVar(icv) = (Cmax/2.0_WP)**2 * (real(icv-1,WP) / real(nCVar-1,WP))**2
    end do

    return
  end subroutine create_prog
  
  
  
  ! Precompute the PDF for a beta distribution with
  ! -> mean mixture fraction : zm
  ! -> variance of mixture fraction : zv
  subroutine create_beta_pdf(Bmean, Bvar)
    use math
    implicit none
    
    real(WP), intent(in) :: Bmean, Bvar ! mean and variance of Beta PDF
    real(WP) :: a,b,factor,tmp,dx,mean
    integer :: index1,n
    
    pdf = 0.0_WP
    
    ! Zero mean : delta at C=0
    if (Bmean.le.1.0E-10_WP) then
       pdf(1) = 1.0_WP
       return
    end if
    
    ! Max mean : delta at C=1
    if (Bmean.ge.1.0_WP-1.0E-10_WP) then
       pdf(np_pdf) = 1.0_WP
       return
    end if
    
    ! Zero variance : delta at C=cm
    if (Bvar.le.1.0E-10_WP) then
       index1 = 1
       do while (pdfaxis(index1).lt.Bmean) 
          index1 = index1+1
       end do
       pdf(index1-1) = (pdfaxis(index1)-Bmean)  /(pdfaxis(index1)-pdfaxis(index1-1))
       pdf(index1)   = (Bmean-pdfaxis(index1-1))/(pdfaxis(index1)-pdfaxis(index1-1))
       return
    end if
        
    ! Impossible cases => two delta at 0 and 1
    if (Bvar.ge.Bmean*(1.0_WP-Bmean)) then
       pdf(1) = 1.0_WP-Bmean
       pdf(np_pdf) = Bmean
       return
    end if

    
    a = Bmean*(Bmean*(1.0_WP-Bmean)/Bvar - 1.0_WP)
    b = a/Bmean - a
    factor = gammaln(a+b) - gammaln(a) - gammaln(b)
    
    ! Left BC : explicit integration
    dx = 0.5_WP*(pdfaxis(2)-pdfaxis(1))
    tmp = a*log(dx) + factor
    pdf(1) = exp(tmp) / a
    ! Right BC : explicit integration
    dx = 0.5_WP*(pdfaxis(np_pdf)-pdfaxis(np_pdf-1))
    tmp = b*log(dx) + factor
    pdf(np_pdf) = exp(tmp) / b
    ! Other Points
    do n=2,np_pdf-1
       dx = 0.5_WP*(pdfaxis(n+1)-pdfaxis(n-1))
       tmp = (a-1.0_WP)*log(pdfaxis(n)) + (b-1.0_WP)*log(1.0_WP-pdfaxis(n))
       tmp = tmp + factor
       pdf(n) = exp(tmp) * dx
    end do

    ! Normalize the pdf
    pdf = pdf / sum(pdf)

    if (isnan(pdf(1))) print *, "nan in pdf"
    
    ! Check mean
    !mean = sum(pdf*Z)
    !pdf(nPoints) = pdf(nPoints) + (zm-mean)
    !pdf(1) = pdf(1) - (zm-mean)
    
    return
  end subroutine create_beta_pdf
  
end module prem_table3D


subroutine prem_table3D_init
  use prem_table3D
  use parser
  implicit none
  
  ! Read the dimension of the final table
  call parser_read('Number of points for C', nC)
  call parser_read('Number of points for Z', nZ)
  call parser_read('Number of points for Cvar', nCVar)
  allocate(Cdex(nC))
  allocate(Zdex(nZ))
  allocate(CVar(nCVar))
  allocate(Ztmp(nFiles))
  allocate(filesCmax(nFiles))
  
  ! Create the first two directions of the table
  call create_prog
  
  ! Allocate arrays
  allocate(postconv(nvar_in,nC,nCVar,nfiles))
  
  return
end subroutine prem_table3D_init
  

! ================================================ !
! Convolute the data with pdfs                     !
! ================================================ !
subroutine prem_table3D_convolute(file)
  use prem_table3D
  implicit none
  
  integer, intent(in) :: file
  integer :: icm, icv, k, var, noccur
  real(WP) :: meanVal, Cmax

  logical :: test !!!!!!!

  ! Fix if C does not monotonically increase and give a warning
  Cmax = C(nPoints)
  print*, 'Cmx =', Cmax
  test = .true.
  noccur = 0
  do icm = 2,nPoints-1
     if (C(icm).lt.C(icm-1) .or. C(icm).ge.Cmax ) then
        if (test) print *, 'Warning: C does not monotonically increase in flamelet. First occurence at C/CMax = ', C(icm)/Cmax, icm
        test = .false.
        noccur = noccur+1
        C(icm) = C(icm-1)
     end if
  end do
  if (noccur.ne.0) print *, 'Number of occurences:', noccur
  print *, ''
  
  ! Prepare the convolution
  np_pdf=Npoints
  allocate(pdf(np_pdf))
  allocate(pdfaxis(np_pdf))
  filesCmax(file) = Cmax
  pdfaxis = C / Cmax


  
  ! Convolutes
  do icv=1,NCvar
     do icm=1,nC
        call create_beta_pdf(Cdex(icm)/Cmax, CVar(icv)/Cmax/Cmax)
        
        do var=1,nvar_in
           meanVal  = 0.0_WP
           
           if (trim(input_name(var)).eq.'density') then
              do k=1,np_pdf
                 meanVal = meanVal + pdf(k)/input_data(k,var)
              end do
              postconv(var,icm,icv,file)  = 1.0_WP / meanVal
           else
              do k=1,np_pdf
                 meanVal = meanVal + pdf(k)*input_data(k,var)
              end do
              postconv(var,icm,icv,file)  = meanVal
           end if
        end do
        
        if (isnan(postconv(4,icm,icv,file))) print *, 'its a nan', isnan(pdf(1)), icm, Cmax, Cdex(icm)/Cmax, CVar(icv)/Cmax/Cmax !!!!
     end do

  end do

  
  ! Add to the 3rd dimension
  if (trim(combModel).eq.'PFPVA-3D') then
     Ztmp(file) = Z
  end if

  ! Finish the convolution
  deallocate(pdf)
  nullify(pdf)
  deallocate(pdfaxis)
  nullify(pdfaxis)
  
  return
end subroutine prem_table3D_convolute


! ========================================================== !
! Convert the names from FlameMaster to user specified names !
! ========================================================== !
subroutine prem_table3D_convert_names
  use prem_table3D
  use parser
  implicit none
  
  character(len=str_medium), dimension(:), pointer :: conversion
  character(len=str_medium) :: varname
  integer :: i,n,var
  
  ! Default name as in FlameMaster
  nvar_out = nvar_in
  allocate(output_name(nvar_out))
  output_name = input_name(1:nvar_in)
  
  ! Get the number of name conversions
  call parser_getsize('Name conversion',n)
  if (mod(n,3).ne.0) stop "table_convert_names: Problem in the definition of conversion names"
  
  ! Allocate array and read
  allocate(conversion(n))
  call parser_read('Name conversion',conversion)
  
  ! Convert the names
  n = n / 3
  do i=1,n
     varname = trim(conversion((i-1)*3+3))
     loop1:do var=1,nvar_in
        if (trim(input_name(var)).eq.trim(varname)) exit loop1
     end do loop1
     if (var.eq.nvar_in+1) then
        print*, "table_convert_names: Unknown variable name : " // varname
        stop
     end if
     output_name(var) = trim(conversion((i-1)*3+1))
  end do
  
  return
end subroutine prem_table3D_convert_names


! ============================================== !
! Setup the table by mapping the second direction !
! ============================================== !
subroutine prem_table3D_setup
  use prem_table3D
  use parser
  implicit none
  
  integer :: ic,iz,var,zcut,ic2,icv
  real(WP) :: minZ, maxZ
  real(WP) :: m11,m12,m21,m22,r1,r2,delta
  real(WP) :: a,b,d, dz
  
  integer :: file, file_up, file_down
  real(WP) :: err, err_up, err_down
  real(WP) :: alpha_up,alpha_down
  real(WP) :: alpha_up2,alpha_down2
  real(WP), dimension(:), pointer  :: interp_up, interp_down
  real(WP) :: Cproj
  
  !real(WP), dimension(:), pointer :: tmp
  character(str_short) :: scale
  
  ! Allocate final table
  allocate(output_data(nC,nCVar,nZ,nvar_out))
  allocate(mask(nC,nCVar,nZ))
  allocate(interp_up(nvar_out))
  allocate(interp_down(nvar_out))
  mask=0
  
  ! Find min and max
  select case (trim(combModel))
  case ('PFPVA-3D')
     maxZ = maxval(Ztmp)
     minZ = 0.0_WP
  end select
  
  ! Linear or log progression
  call parser_read('Scale for Z direction',scale)
  call parser_read('Max value for Z direction',maxZ,maxZ)
  select case(trim(scale))
  case ('lin')
     if (nZ.ne.1) then
        do iz=1,nZ
           Zdex(iz) = minZ + real(iz-1,WP)*(maxZ-minZ)/real(nZ-1,WP)
        end do
     else
        Zdex(nZ) = maxZ 
     end if
  case ('log')
     if (nZ.ne.1) then
        do iz=1,nZ
           Zdex(iz) = minZ * (maxZ/minZ) ** (real(iz-1,WP)/real(nZ-1,WP))
        end do
     else
        Zdex(iz) = maxZ
     end if
  case ('stretch') ! linear and fine on lean side, stretched on rich side
     if (nz.ne.1) then
        zcut = nZ/3+1
        dz = (Zst - minZ) / real(zcut-1,WP)
        ! Constant mesh spacing for [minZ,Zst]
        do iz=1,zcut 
           Zdex(iz) = minZ + real(iz-1,WP) * dz
        end do
        ! Mesh with linear growth to reach Z=maxZ
        m11 = real(nZ**2-zcut**2,WP)
        m12 = real(nZ-zcut,WP)
        m21 = real(2*zcut+1,WP)
        m22 = real(1,WP)
        r1 = maxZ - Zst
        r2 = dz
        delta = m11*m22-m12*m21
        a = (+ m22*r1 - m12*r2 )/delta
        b = (- m21*r1 + m11*r2 )/delta
        d = Zst - a*zcut**2-b*zcut
        do iz=zcut+1,nZ
           Zdex(iz) = a*real(iz,WP)**2 + b*real(iz,WP) + d
        end do
     else
        Zdex(iz) = maxZ
     end if
     print *, 'Zdex ', Zdex
  case default
     stop "table_setup: Unknown Scale for second direction"
  end select

  ! Loop over the two mapping directions
  do iz=1,nZ
     do icv=1,nCVar
        do ic=1,nC
           
           !tmp => postconv(nvar_in,ic,:)
           
           ! Find the two files right above and right below
           err_up = huge(1.0_WP)
           err_down = -huge(1.0_WP)
           
           file_up   = 0
           file_down = 0
           
           do file=1,nfiles
              err = Ztmp(file) - Zdex(iz)
              
              if ((err.ge.0.0_WP) .and. (err.le.err_up)) then
                 file_up = file
                 err_up = err
              end if
              if ((err.le.0.0_WP) .and. (err.ge.err_down)) then
                 file_down = file
                 err_down = err
              end if
           end do
           
           ! Interpolate
           if (file_up.eq.0 .or. file_down.eq.0) then
              if (file_up.eq.0) then
                 alpha_up   = 0.0_WP
                 alpha_down = 1.0_WP
                 file_up = 1
              end if
              if (file_down.eq.0) then
                 alpha_up   = 1.0_WP
                 alpha_down = 0.0_WP
                 file_down = 1
              end if
              ! Mask it
              mask(ic,icv,iz) = 1
           else
              if (file_up.eq.file_down) then
                 alpha_up   = 1.0_WP
                 alpha_down = 0.0_WP
              else
                 alpha_up   = (Zdex(iz)-Ztmp(file_down)) / (Ztmp(file_up)-Ztmp(file_down))
                 alpha_down = (Ztmp(file_up)-Zdex(iz))   / (Ztmp(file_up)-Ztmp(file_down))
              end if
           end if

           if (file_up.eq.file_down) then
              interp_up = postconv(:,ic,icv,file_up)
              interp_down = postconv(:,ic,icv,file_down)

           !    ! region with higher than equilibrium progress variable
           !    ! project upward from equilibrium
           ! else if ((Cdex(ic).gt. &
           !      filesCmax(file_down) + (filesCmax(file_up) - filesCmax(file_down)) * &
           !      (Zdex(iz) - Ztmp(file_down)) / (Ztmp(file_up) - Ztmp(file_down) ))   &
           !      .or. (file_up.eq.file_down) ) then
           !    interp_up = postconv(:,Nc,icv,file_up)
           !    interp_down = postconv(:,Nc,icv,file_down)
              
           !    ! deal with region between lowest burning flamelet and Z=0 flamelet
           !    ! diagonal interpolation
           ! else if (Cdex(ic).gt.filesCmax(file_down)) then
           !    Cproj = filesCmax(file_down) + (Cdex(ic) - filesCmax(file_down)) * &
           !         (Ztmp(file_up) - Ztmp(file_down)) / (Zdex(iz) - Ztmp(file_down))
           !    ic2 = ic
           !    do while (Cdex(ic2).lt.Cproj)
           !       ic2 = ic2+1
           !    end do
           !    alpha_up2   = (Cproj - Cdex(ic2)) / (Cdex(ic2+1) - Cdex(ic2))
           !    alpha_down2 = 1.0_WP - alpha_up2
           !    do var=1,nvar_out
           !       if (trim(input_name(var)).eq.'density') then
           !          interp_up(var) =  1.0_WP/( &
           !               alpha_up2   / postconv(var,ic2+1,icv,file_up) + &
           !               alpha_down2 / postconv(var,ic2,icv,file_up) )
           !       else
           !          interp_up(var) =  &
           !               alpha_up2  *postconv(var,ic2+1,icv,file_up) + &
           !               alpha_down2*postconv(var,ic2,icv,file_up)
           !       end if
           !       interp_down(var) = postconv(var,ic,icv,file_down)
           !    end do
              
           !    ! Region between highest burning flamelet and Z=1 flamelet
           !    ! diagonal interpolation
           ! else if (Cdex(ic).gt.filesCmax(file_up)) then
           !    Cproj = filesCmax(file_up) + (Cdex(ic) - filesCmax(file_up)) * &
           !         (Ztmp(file_down) - Ztmp(file_up)) / (Zdex(iz) - Ztmp(file_up))
           !    ic2 = ic
           !    do while (Cdex(ic2).lt.Cproj)
           !       ic2 = ic2+1
           !    end do
           !    alpha_up2   = (Cproj - Cdex(ic2)) / (Cdex(ic2+1) - Cdex(ic2))
           !    alpha_down2 = 1.0_WP - alpha_up2
           !    do var=1,nvar_out
           !       if (trim(input_name(var)).eq.'density') then
           !          interp_down(var) =  1.0_WP/( &
           !               alpha_up2   / postconv(var,ic2+1,icv,file_down) + &
           !               alpha_down2 / postconv(var,ic2,icv,file_down) )
           !       else
           !          interp_down(var) =  &
           !               alpha_up2  *postconv(var,ic2+1,icv,file_down) + &
           !               alpha_down2*postconv(var,ic2,icv,file_down)
           !       end if
           !       interp_up(var) = postconv(var,ic,icv,file_up)
           !    end do
              
              ! Region between burning flamelets
              ! horizontal interpolation
           else
              interp_up = postconv(:,ic,icv,file_up)
              interp_down = postconv(:,ic,icv,file_down)
           end if

           ! Interpolate
           do var=1,nvar_out
              if (trim(input_name(var)).eq.'density') then
                 output_data(ic,icv,iz,var) =  1.0_WP/( &
                      alpha_up   / interp_up(var) + &
                      alpha_down / interp_down(var) )
              else
                 output_data(ic,icv,iz,var) =  &
                      alpha_up   * interp_up(var)+ &
                      alpha_down * interp_down(var)
              end if
           end do
        end do
     end do
  end do
  
  return
end subroutine prem_table3D_setup



! ===================== !
! Print some statistics !
! ===================== !
subroutine prem_table3D_stats
  use prem_table3D
  implicit none
  
  integer :: var
  
  print*,''
  
  ! Min and Max of the coordinates
  print*, '** Coordinates of the table **'
  write(*,10) 'Variable    ','Min    ','Max    '
  write(*,10) '------------','------------','------------'
  write(*,11) 'C           ', minval(Cdex), maxval(Cdex)
  write(*,11) 'CVar        ', minval(CVar), maxval(CVar)
  write(*,11) 'Z           ', minval(Zdex), maxval(Zdex)
  print*,''
  
  ! Min and Max of all the mapped quantities
  print*, '** Mapped quantities **'
  write(*,10) 'Variable    ','Min    ','Max    '
  write(*,10) '------------','------------','------------'
  do var=1,nvar_out
     write(*,11) output_name(var),minval(output_data(:,:,:,var)),maxval(output_data(:,:,:,var))
  end do
  print*,''
  
10 format (A12,'  ',A12,'  ',A12)
11 format (A12,'  ',ES12.4,'  ',ES12.4)
  
  return
end subroutine prem_table3D_stats


! ===================================== !
! Write the table back to a binary file !
! ===================================== !
subroutine prem_table3D_write
  use prem_table3D
  use parser
  implicit none
  
  integer :: ierr,var,iunit
  character(len=str_medium) :: filename
  real(WP) :: tmp_val
  character(len=str_medium) :: tmp_str
  !integer :: i,j,k
  
  ! Open the data file
  call parser_read('Table filename', filename)
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)
  
!!$  ! Write VIDA header
!!$  tmp_val = 0.0_WP
!!$  call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
!!$  tmp_str = 'FPVA'
!!$  call BINARY_FILE_WRITE(iunit,tmp_str,str_medium,kind(tmp_str),ierr)
!!$  tmp_str = 'SPECIES'
!!$  call BINARY_FILE_WRITE(iunit,tmp_str,str_medium,kind(tmp_str),ierr)
!!$  tmp_val = 1.0133e5_WP
!!$  call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
!!$  tmp_val = 298.0_WP
!!$  call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
!!$  call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)

  ! Write sizes
  call BINARY_FILE_WRITE(iunit,nC,1,kind(nC),ierr)
  call BINARY_FILE_WRITE(iunit,nCVar,1,kind(nCVar),ierr)
  call BINARY_FILE_WRITE(iunit,nZ,1,kind(nZ),ierr)
  call BINARY_FILE_WRITE(iunit,nvar_out,1,kind(nvar_out),ierr)
  ! Write the axis coordinates
  call BINARY_FILE_WRITE(iunit,Cdex,nC,kind(Cdex),ierr)
  call BINARY_FILE_WRITE(iunit,CVar,nCVar,kind(CVar),ierr)
  call BINARY_FILE_WRITE(iunit,Zdex,nZ,kind(Zdex),ierr)
  ! Masks
!!$  call BINARY_FILE_WRITE(iunit,mask,nC*nZ,kind(mask),ierr)
  ! Write additional stuff
  call BINARY_FILE_WRITE(iunit,combModel,str_medium,kind(combModel),ierr)
  ! Write variable names
  do var=1,nvar_out
     call BINARY_FILE_WRITE(iunit,output_name(var),str_medium,kind(output_name),ierr)
  end do
  ! Write data field
  do var=1,nvar_out
     call BINARY_FILE_WRITE(iunit,output_data(:,:,:,var),nC*nCVar*nZ,kind(output_data),ierr)
  end do
!!$  ! Backwards order for VIDA (written in C++)
!!$  do i=1,nC
!!$     do j=1,nZ
!!$           do var=1,nvar_out
!!$              call BINARY_FILE_WRITE(iunit,output_data(i,j,var),1,kind(output_data),ierr)
!!$           end do
!!$     end do
!!$  end do
  
  call BINARY_FILE_CLOSE(iunit,ierr)
  
  return
end subroutine prem_table3D_write
