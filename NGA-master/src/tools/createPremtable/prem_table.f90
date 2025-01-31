module prem_table
  use prem_flamelet
  use precision
  use string
  implicit none
  
  ! Table parameters
  integer :: nC, nZ
  real(WP), dimension(:), pointer :: Cdex, Zdex, Ztmp
  
  ! Beta pdf
  real(WP), dimension(:), pointer :: pdf
  
  ! List of Flamelet files
  integer :: nfiles
  character(len=str_long), dimension(:), pointer :: files
  
  ! Variables to be wrtitten into the table
  integer :: nvar_out
  character(len=str_medium), dimension(:), pointer :: output_name
  real(WP), dimension(:,:,:), pointer :: output_data
  integer, dimension(:,:), pointer :: mask
  
  ! Variable after convolution
  real(WP), dimension(:,:,:), pointer :: postconv
  
contains
  
  ! Uniform mesh up to ...hmmm... 0.3 for now
  subroutine create_prog
    use parser
    implicit none
    
    integer :: cc
    real(WP) :: Cmax

    call parser_read("Maximum Progress Variable", Cmax)
        
    do cc=1,nC
       Cdex(cc) = Cmax* real(cc-1,WP)/real(nC-1,WP)
    end do
    
    return
  end subroutine create_prog
  
  
  ! Uniform mesh between 0 and 0.07
!!$  subroutine create_Z
!!$    implicit none
!!$    
!!$    integer :: zv
!!$    
!!$    do zv=1,nZ
!!$       Zdex(zv) = 0.07_WP* real(zv-1,WP) / real(nZ-1,WP)
!!$    end do
!!$    
!!$    return
!!$  end subroutine create_Z
  
  
  ! Precompute the PDF for a beta distribution with
  ! -> mean mixture fraction : zm
  ! -> variance of mixture fraction : zv
  subroutine create_beta_pdf(vc)
    use math
    implicit none
    
    real(WP), intent(in) :: vc
    !real(WP) :: a,b,factor,tmp,dz,mean
    integer :: index1
    
    pdf = 0.0_WP
    
!!$    ! Zero mean : delta at Z=0
!!$    if (zm.le.1.0E-10_WP) then
!!$       pdf(1) = 1.0_WP
!!$       return
!!$    end if
!!$    
!!$    ! Max mean : delta at Z=1
!!$    if (zm.ge.1.0_WP-1.0E-10_WP) then
!!$       pdf(nPoints) = 1.0_WP
!!$       return
!!$    end if
!!$    
    ! Zero variance : delta at C=cm
    if (C(nPoints).lt.vc) then
       pdf(nPoints) = 1.0_WP
       return
    elseif (vc.lt.C(1)) then
       pdf(1) = 1.0_WP
       return
    else
       index1 = 1
       do while (C(index1).lt.vc) 
          index1 = index1+1
       end do
       pdf(index1-1) = (C(index1)-vc)  /(C(index1)-C(index1-1))
       pdf(index1)   = (vc-C(index1-1))/(C(index1)-C(index1-1))
       return
    endif
        
!!$    ! Impossible cases => two delta at 0 and 1
!!$    if (zv.ge.zm*(1.0_WP-zm)) then
!!$       pdf(1) = 1.0_WP-zm
!!$       pdf(nPoints) = zm
!!$       return
!!$    end if
!!$    
!!$    a = zm*(zm*(1.0_WP-zm)/zv - 1.0_WP)
!!$    b = a/zm - a
!!$    factor = gammaln(a+b) - gammaln(a) - gammaln(b)
!!$    
!!$    ! Left BC : explicit integration
!!$    dz = 0.5_WP*(Z(2)-Z(1))
!!$    tmp = a*log(dz) + factor
!!$    pdf(1) = exp(tmp) / a
!!$    ! Right BC : explicit integration
!!$    dz = 0.5_WP*(Z(nPoints)-Z(nPoints-1))
!!$    tmp = b*log(dz) + factor
!!$    pdf(nPoints) = exp(tmp) / b
!!$    ! Other Points
!!$    do n=2,nPoints-1
!!$       dz = 0.5_WP*(Z(n+1)-Z(n-1))
!!$       tmp = (a-1.0_WP)*log(Z(n)) + (b-1.0_WP)*log(1.0_WP-Z(n))
!!$       tmp = tmp + factor
!!$       pdf(n) = exp(tmp) * dz
!!$    end do
    
!!$    ! Normalize the pdf
!!$    pdf = pdf / sum(pdf)
!!$    
!!$    ! Check mean
!!$    !mean = sum(pdf*Z)
!!$    !pdf(nPoints) = pdf(nPoints) + (zm-mean)
!!$    !pdf(1) = pdf(1) - (zm-mean)
!!$    
!!$    return
  end subroutine create_beta_pdf
  
end module prem_table


subroutine prem_table_init
  use prem_table
  use parser
  implicit none
  
  ! Read the dimension of the final table
  call parser_read('Number of points for C', nC)
  call parser_read('Number of points for Z', nZ)
  allocate(Cdex(nC))
  allocate(Zdex(nZ))
  allocate(Ztmp(nFiles))
  
  ! Create the first two directions of the table
  call create_prog
  !call create_Z
  
  ! Allocate arrays
  allocate(postconv(nvar_in,nC,nfiles))
  
  return
end subroutine prem_table_init
  

! ================================================ !
! Convolute the data with pdfs                     !
! ================================================ !
subroutine prem_table_convolute(file)
  use prem_table
  implicit none
  
  integer, intent(in) :: file
  integer :: ic, k, var
  real(WP) :: meanVal
  
  ! Prepare the convolution
  allocate(pdf(nPoints))
  print*, 'Cmx =', C(nPoints)
  
  ! Convolutes
     do ic=1,nC
        call create_beta_pdf(Cdex(ic))
        
        do var=1,nvar_in
           meanVal  = 0.0_WP
           
           if (trim(input_name(var)).eq.'density') then
              do k=1,nPoints
                 meanVal = meanVal + pdf(k)/input_data(k,var)
              end do
              postconv(var,ic,file)  = 1.0_WP / meanVal
           else
              do k=1,nPoints
                 meanVal = meanVal + pdf(k)*input_data(k,var)
              end do
              postconv(var,ic,file)  = meanVal
           end if
        end do
     end do
  
  ! Add to the 3rd dimension
  if (trim(combModel).eq.'PFPVA') then
     Ztmp(file) = Z
  end if

  ! Finish the convolution
  deallocate(pdf)
  nullify(pdf)
  
  return
end subroutine prem_table_convolute


! ========================================================== !
! Convert the names from FlameMaster to user specified names !
! ========================================================== !
subroutine prem_table_convert_names
  use prem_table
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
end subroutine prem_table_convert_names


! ============================================== !
! Setup the table by mapping the second direction !
! ============================================== !
subroutine prem_table_setup
  use prem_table
  use parser
  implicit none
  
  integer :: ic,iz,var
  real(WP) :: minZ, maxZ
  
  integer :: file, file_up, file_down
  real(WP) :: err, err_up, err_down
  real(WP) :: alpha_up,alpha_down
  
  !real(WP), dimension(:), pointer :: tmp
  character(str_short) :: scale
  
  ! Allocate final table
  allocate(output_data(nC,nZ,nvar_out))
  allocate(mask(nC,nZ))
  mask=0
  
  ! Find min and max
  select case (trim(combModel))
  case ('PFPVA')
     maxZ = maxval(Ztmp)
     minZ = 0.0_WP
  case default
     maxZ = maxval(postconv(nvar_in,:,:))
     minZ = minval(postconv(nvar_in,:,:))
  end select
  
  ! Linear or log progression
  call parser_read('Scale for second direction',scale)
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
  case default
     stop "table_setup: Unknown Scale for second direction"
  end select
  
  ! Loop over the two mapping directions
  do iz=1,nZ
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
              mask(ic,iz) = 1
           else
              if (file_up.eq.file_down) then
                 alpha_up   = 1.0_WP
                 alpha_down = 0.0_WP
              else
                 alpha_up   = (Zdex(iz)-Ztmp(file_down)) / (Ztmp(file_up)-Ztmp(file_down))
                 alpha_down = (Ztmp(file_up)-Zdex(iz))   / (Ztmp(file_up)-Ztmp(file_down))
              end if
           end if
           
           do var=1,nvar_out
              if (trim(input_name(var)).eq.'density') then
                 output_data(ic,iz,var) =  1.0_WP/( &
                      alpha_up  /postconv(var,ic,file_up) + &
                      alpha_down/postconv(var,ic,file_down) )
              else
                 output_data(ic,iz,var) =  &
                      alpha_up  *postconv(var,ic,file_up) + &
                      alpha_down*postconv(var,ic,file_down)
              end if
           end do
           
        end do
  end do
  
  return
end subroutine prem_table_setup



! ===================== !
! Print some statistics !
! ===================== !
subroutine prem_table_stats
  use prem_table
  implicit none
  
  integer :: var
  
  print*,''
  
  ! Min and Max of the coordinates
  print*, '** Coordinates of the table **'
  write(*,10) 'Variable    ','Min    ','Max    '
  write(*,10) '------------','------------','------------'
  write(*,11) 'C           ', minval(Cdex), maxval(Cdex)
  write(*,11) 'Z           ', minval(Zdex), maxval(Zdex)
  print*,''
  
  ! Min and Max of all the mapped quantities
  print*, '** Mapped quantities **'
  write(*,10) 'Variable    ','Min    ','Max    '
  write(*,10) '------------','------------','------------'
  do var=1,nvar_out
     write(*,11) output_name(var),minval(output_data(:,:,var)),maxval(output_data(:,:,var))
  end do
  print*,''
  
10 format (A12,'  ',A12,'  ',A12)
11 format (A12,'  ',ES12.4,'  ',ES12.4)
  
  return
end subroutine prem_table_stats


! ===================================== !
! Write the table back to a binary file !
! ===================================== !
subroutine prem_table_write
  use prem_table
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
  call BINARY_FILE_WRITE(iunit,nZ,1,kind(nZ),ierr)
 ! call BINARY_FILE_WRITE(iunit,n3,1,kind(n3),ierr)
  call BINARY_FILE_WRITE(iunit,nvar_out,1,kind(nvar_out),ierr)
  ! Write the axis coordinates
  call BINARY_FILE_WRITE(iunit,Cdex,nC,kind(Cdex),ierr)
  call BINARY_FILE_WRITE(iunit,Zdex,nZ,kind(Zdex),ierr)
!  call BINARY_FILE_WRITE(iunit,Z3,n3,kind(Z3),ierr)
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
     call BINARY_FILE_WRITE(iunit,output_data(:,:,var),nC*nZ,kind(output_data),ierr)
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
end subroutine prem_table_write
