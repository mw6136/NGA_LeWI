module unsteady_prem_table
  use unsteady_prem_flamelet
  use precision
  use string
  implicit none
  
  ! Table parameters
  integer :: nC, nZ, nH
  real(WP), dimension(:), pointer :: Cdex, Zdex, Hdex, Ztmp, Htmp
  
  ! Beta pdf
  real(WP), dimension(:), pointer :: pdf
  
  ! List of Flamelet files
  integer :: nfiles
  character(len=str_xl), dimension(:), pointer :: files


  
  ! Variables to be written into the table
  integer :: nvar_out
  character(len=str_medium), dimension(:), pointer :: output_name
  character(len=str_medium) :: mask_name
  real(WP), dimension(:,:,:,:), pointer :: output_data
  integer, dimension(:,:), pointer :: mask1
  integer, dimension(:,:,:), pointer :: mid_mask
  real(WP), dimension(:,:,:), pointer :: mask
  
  ! Variable after convolution
  real(WP), dimension(:,:,:), pointer :: postconv

  ! VIDA or NGA?
  logical :: vida_table, modenth, masking
  
contains

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

!    print*, C(nPoints)
    
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
    
    return
  end subroutine create_beta_pdf
  
end module unsteady_prem_table


subroutine unsteady_prem_table_init
  use unsteady_prem_table
  use parser
  implicit none
  
  ! Read the dimension of the final table
  call parser_read('Number of points for C', nC)
  !call parser_read('Number of points for variance of C', nCVar)
  allocate(Cdex(nC))
  !allocate(CVar(nCVar))
  call parser_read('Number of points for Z', nZ)
  call parser_read('Number of points for H', nH)
  call parser_read('Use masking',masking,.false.)
  allocate(Zdex(nZ))
  allocate(Hdex(nH))
  allocate(Ztmp(nFiles))
  allocate(Htmp(nFiles))
  
  ! Create the first two directions of the table
  call create_prog
  !call create_cvar
  
  ! Allocate arrays
  allocate(postconv(nvar_in,nC,nfiles))

  if (masking) then
     allocate(mask1(nC,nfiles))
     allocate(mask(nC,nZ,nH))
  end if

  ! Table for VIDA or NGA
  call parser_read('VIDA table',vida_table,.false.)
  
  return
end subroutine unsteady_prem_table_init
  

! ================================================ !
! Convolute the data with pdfs of Mixture Fraction !
! ================================================ !
subroutine unsteady_prem_table_convolute(file)
  use unsteady_prem_table
  implicit none
  
  integer, intent(in) :: file
  integer :: ic, k, var
  real(WP) :: meanVal
  
  ! Prepare the convolution
  allocate(pdf(nPoints))

!  print*, nPoints

!!$  if (file.eq.2) then
!!$     do var = 1,nvar_in
!!$        if (trim(input_name(var)).eq.'temperature') then
!!$           print*, 'inputdata', input_data(:,var)
!!$        end if
!!$     end do
!!$  end if
  
  ! Convolutes
     do ic=1,nC
        call create_beta_pdf(Cdex(ic))

!!$        if (file.eq.2 .and. ic.eq.60) then
!!$           print*, Cdex(ic)
!!$           do var = 1,nvar_in
!!$              if (trim(input_name(var)).eq.'temperature') then
!!$                 print*, 'pdf of 1', pdf(:)
!!$              end if
!!$           end do
!!$        end if
        
        do var=1,nvar_in
           meanVal  = 0.0_WP
           
           if (trim(input_name(var)).eq.'density') then
              do k=1,nPoints
                 meanVal = meanVal + pdf(k)/input_data(k,var)
              end do
              postconv(var,ic,file)  = 1.0_WP / meanVal
!!$              if (file.eq.39 .or. file.eq.73) then
!!$                 print*,'postconv1' , postconv(var,ic,file)
!!$              end if
           else
              do k=1,nPoints
                 meanVal = meanVal + pdf(k)*input_data(k,var)
              end do
              postconv(var,ic,file)  = meanVal
              if (masking) then
                 if (ic.eq.1 .or. ic.eq.nC) then
                    if (pdf(1).eq.1.0_WP .and. Cdex(ic).lt.C(1)) then
                       mask1(ic,file) = -1
                    elseif (pdf(nPoints).eq.1.0_WP .and. C(nPoints).lt.Cdex(ic)) then
                       mask1(ic,file) = -1
                    else
                       mask1(ic,file) = 1
                    end if
                 else
                    if (pdf(1).eq.1.0_WP .and. Cdex(ic).lt.C(1)) then ! .and. Cdex(ic+1).lt.C(1)
                       mask1(ic,file) = -1
                    elseif (pdf(nPoints).eq.1.0_WP .and. C(nPoints).lt.Cdex(ic)) then ! .and. C(nPoints).lt.Cdex(ic-1)
                       mask1(ic,file) = -1
                    else
                       mask1(ic,file) = 1
                    end if 
                 end if
              end if
           end if
        end do
     end do

  ! Add to the 3rd dimension
  if (trim(combModel).eq.'RPFPVA') then
     Ztmp(file) = Z
     if (prodconv) then
        Htmp(file) = ENTH
     end if
  end if

  ! Finish the convolution
  deallocate(pdf)
  nullify(pdf)
  
  return
end subroutine unsteady_prem_table_convolute


! ========================================================== !
! Convert the names from FlameMaster to user specified names !
! ========================================================== !
subroutine unsteady_prem_table_convert_names
  use unsteady_prem_table
  use parser
  implicit none
  
  character(len=str_medium), dimension(:), pointer :: conversion
  character(len=str_medium) :: varname
  integer :: i,n,var
  
  ! Default name as in FlameMaster
  nvar_out = nvar_in
  allocate(output_name(nvar_out))
  output_name = input_name(1:nvar_in)
  mask_name = 'mask'
  
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
end subroutine unsteady_prem_table_convert_names


! ========================================================== !
! Setup the table by mapping the third and fourth directions !
! ========================================================== !
subroutine unsteady_prem_table_setup
  use unsteady_prem_table
  use parser
  use fileio
  implicit none
  
  integer :: ic,iz,ih,var,i,j,ierr,ip
  integer :: var_density,var_temp
  real(WP) :: minZ, maxZ, minH, maxH, Hmin, ZTop, ZBot

  integer, dimension(:), pointer :: plane_num

  integer :: file, file_up, file_down, plane_up, plane_down, XLZ, SMZ
  real(WP) :: err, err_up, err_down
  real(WP) :: alpha_up,alpha_down
  
!!$  integer :: file, file1, file2, file3
!!$  real(WP) :: err, err1, err2, err3
!!$  real(WP) :: alpha_1,alpha_2,alpha_3,denom
  
  real(WP), dimension(:), pointer :: tmpZ,tmpH,Zn,Hn,Zplane,tmp,ZEnth
  character(str_short) :: scale

  ! Mid-interpolation structure
  real(WP), dimension(:,:,:,:), pointer :: mid_data

!!$  do var = 1,nvar_out
!!$     if (trim(input_name(var)).eq.'temperature') then
!!$        print*, 'postconv', postconv(var,:,2)
!!$        print*, 'postconv', postconv(var,:,5)
!!$     end if
!!$  end do

  ZTop = 0.0_WP
  ZBot = 0.01_WP

!  print*, 'postconv imm postconv', postconv(nvar_in,:,10)

  ! Associate files with c-h planes
  allocate(plane_num(nfiles))
  plane_num = 0
  do i = 1,nfiles
     if (i.eq.1) then
        plane_num(i) = 1
     else
        loop_check: do j=1,i-1
           if (abs(Ztmp(i)-Ztmp(j)).lt.0.00001_WP) then
              plane_num(i) = plane_num(j)
              exit loop_check
           endif
           if (j .eq. i-1) then
              plane_num(i) = maxval(plane_num(:)) + 1
           endif
        end do loop_check
     end if
     if (Ztmp(i).gt.ZTop) then
        ZTop = Ztmp(i)
        XLZ = plane_num(i)
     end if
     if (Ztmp(i).lt.ZBot) then
        ZBot = Ztmp(i)
        SMZ = plane_num(i)
     end if
  end do

  if (prodconv) then
     ! Find adiabatic flamelets
     allocate(ZEnth(maxval(plane_num(:))))
     do i = 1,nfiles
        if (postconv(nvar_in,1,i).eq.0.0_WP) then
           ZEnth(plane_num(i)) = Htmp(i)
        end if
     end do

     ! Get H
     do i = 1,nfiles
        do j = 1,nC
           postconv(nvar_in,j,i) =  Htmp(i) - ZEnth(plane_num(i))
        end do
     end do
  end if

!  print*,'postconv after try' , postconv(nvar_in,:,10)

  ! Allocate mid and final tables
  allocate(mid_data(nC,nH,maxval(plane_num(:)),nvar_out))
  allocate(output_data(nC,nZ,nH,nvar_out))
!!$  if (nZ.eq.1) then
!!$     allocate(output_data(nC,nZ+1,nH,nvar_out))
!!$  end if
  !allocate(mask(nC,nZ,nH))
  allocate(Zplane(maxval(plane_num)))
  if (masking) then
     allocate(mid_mask(nC,nH,maxval(plane_num(:))))
  end if
 ! mask=0
  
  ! Find min and max
  call parser_read('Hmin defined',modenth,.false.)
  if (modenth) then
     call parser_read("Maximum Enthalpy Loss", Hmin)
  end if

  select case (trim(combModel))
  case ('RPFPVA')
     maxZ = maxval(Ztmp)
     minZ = 0.0_WP
     maxH = maxval(postconv(nvar_in,:,:))
     if (modenth) then
        minH = Hmin
     else
        minH = minval(postconv(nvar_in,:,:))
     end if
  case default
     maxZ = maxval(postconv(nvar_in-2,:,:))
     minZ = minval(postconv(nvar_in-2,:,:))
     maxH = maxval(postconv(nvar_in,:,:))
     minH = minval(postconv(nvar_in,:,:))
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
  case ('quad')
     if (nZ.ne.1) then
        do iZ=1,nZ
           Zdex(iz) = minZ + (maxZ-minZ) * (real(iz-1,WP)/real(nZ-1,WP))**2
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
        Zdex(nZ) = maxZ 
     end if
  case default
     stop "table_setup: Unknown Scale for third direction"
  end select

  ! Linear or log progression
  call parser_read('Scale for third direction',scale)
  select case(trim(scale))
  case ('lin')
     do ih=1,nH
        Hdex(ih) = minH + real(ih-1,WP)*(maxH-minH)/real(nH-1,WP)
     end do
  case ('quad')
     do ih=1,nH
        !Z4(i4) = min4 + (max4-min4) * (real(i4-1,WP)/real(n4-1,WP))**2
        Hdex(ih) = maxH - (maxH-minH) * ((real(ih-1,WP)/real(nH-1,WP))-1.0_WP)**2
     end do
  case ('log')
     do ih=1,nH
        Hdex(ih) = minH * (maxH/minH) ** (real(ih-1,WP)/real(nH-1,WP))
     end do
  case default
     stop "table_setup: Unknown Scale for second direction"
  end select

!!$  print*, 'nPlanes', maxval(plane_num(:))

  ! Loop over the two mapping directions
  do ip=1,maxval(plane_num(:))
     do ih=1,nH
        do ic=1,nC
           
           tmp => postconv(nvar_in,ic,:)
           
           ! Find the two files right above and right below
           err_up = huge(1.0_WP)
           err_down = -huge(1.0_WP)
           
           file_up   = 0
           file_down = 0
           
           do file=1,nfiles
              if (ip.eq.plane_num(file)) then

                 Zplane(ip) = Ztmp(file)

                 err = tmp(file) - Hdex(ih)

                 if ((err.ge.0.0_WP) .and. (err.le.err_up)) then
                    file_up = file
                    err_up = err
                 end if
                 if ((err.le.0.0_WP) .and. (err.ge.err_down)) then
                    file_down = file
                    err_down = err
                 end if
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
              if (masking) then
                 mid_mask(ic,ih,ip) = -1
              end if
           else
              if (file_up.eq.file_down) then
                 alpha_up   = 1.0_WP
                 alpha_down = 0.0_WP
              else
                 alpha_up   = (Hdex(ih)-tmp(file_down)) / (tmp(file_up)-tmp(file_down))
                 alpha_down = (tmp(file_up)-Hdex(ih))   / (tmp(file_up)-tmp(file_down))
              end if
              
              if (masking) then
                 if (mask1(ic,file_up).eq.-1 .or. mask1(ic,file_down).eq.-1) then
                    mid_mask(ic,ih,ip) = -1
                 else
                    mid_mask(ic,ih,ip) = 1
                 end if
              end if
           end if
           
           do var=1,nvar_out
              if (trim(input_name(var)).eq.'density') then
                 mid_data(ic,ih,ip,var) =  1.0_WP/( &
                      alpha_up  /postconv(var,ic,file_up) + &
                      alpha_down/postconv(var,ic,file_down) )
              else
                 mid_data(ic,ih,ip,var) =  &
                      alpha_up  *postconv(var,ic,file_up) + &
                      alpha_down*postconv(var,ic,file_down)
              end if
           end do
           
        end do
     end do
  end do

!!$  do var = 1,nvar_out
!!$     if (trim(input_name(var)).eq.'temperature') then
!!$        var_density = var
!!$        print*, 'mid', mid_data(:,nH-1,1,var)
!!$        print*, 'mid', mid_data(:,nH/2,1,var)
!!$     end if
!!$  end do

!!$  print*, 'ZPlanes', Zplane(:)
  
!!$  do var=1,nvar_out
!!$     if (trim(input_name(var)).eq.'density') then
!!$        print*, 'midconv1', mid_data(:,nH,39,var), 'midconv2', mid_data(:,nH-50,39,var)
!!$     end if
!!$  end do

  !call interpolate_unsteady_prem_flamelets_delaunay

  do iz = 1,nZ
     ! Find the two files right above and right below
     err_up = huge(1.0_WP)
     err_down = -huge(1.0_WP)

     plane_up   = 0
     plane_down = 0

     do ip=1,maxval(plane_num(:))
        
        err = Zplane(ip) - Zdex(iz)

        if ((err.ge.0.0_WP) .and. (err.le.err_up)) then
           plane_up = ip
           err_up = err
        end if
        if ((err.le.0.0_WP) .and. (err.ge.err_down)) then
           plane_down = ip
           err_down = err
        end if

     end do

     ! Interpolate
     if (plane_up.eq.0 .or. plane_down.eq.0) then
        if (plane_up.eq.0) then
           alpha_up   = 0.0_WP
           alpha_down = 1.0_WP
           plane_up = XLZ
        end if
        if (plane_down.eq.0) then
           alpha_up   = 1.0_WP
           alpha_down = 0.0_WP
           plane_down = SMZ
        end if

        ! Mask it
        !mask(ic,ih,iz) = 1
     else
        if (plane_up.eq.plane_down) then
           alpha_up   = 1.0_WP
           alpha_down = 0.0_WP
        else
           alpha_up   = (Zdex(iz)-Zplane(plane_down)) / (Zplane(plane_up)-Zplane(plane_down))
           alpha_down = (Zplane(plane_up)-Zdex(iz))   / (Zplane(plane_up)-Zplane(plane_down))
        end if
     end if

!!$     if (iz.eq.nZ) then
!!$        print*, 'plane up', plane_up, 'plane down', plane_down, 'alpha_up', alpha_up, 'alpha_down', alpha_down
!!$     end if

     do var=1,nvar_out
        do ic=1,nC
           do ih=1,nH
              if (trim(input_name(var)).eq.'density') then
                 output_data(ic,iz,ih,var) =  1.0_WP/( &
                      alpha_up  /mid_data(ic,ih,plane_up,var) + &
                      alpha_down/mid_data(ic,ih,plane_down,var) )
              else
                 output_data(ic,iz,ih,var) =  &
                      alpha_up  *mid_data(ic,ih,plane_up,var) + &
                      alpha_down*mid_data(ic,ih,plane_down,var)
              end if

              if (masking) then
                 if (mid_mask(ic,ih,plane_up).eq.-1 .and. mid_mask(ic,ih,plane_down).eq.-1) then
                    mask(ic,iz,ih) = -1.0_WP
                 else if (plane_up.eq.0 .or. plane_down.eq.0) then
                    mask(ic,iz,ih) = -1.0_WP
                 else
                    mask(ic,iz,ih) = 1.0_WP
                 end if
              end if
           end do
        end do
     end do
  end do

!!$  do var = 1,nvar_out
!!$     if (trim(input_name(var)).eq.'temperature') then
!!$        var_density = var
!!$        print*, 'out', output_data(:,nZ,nH-1,var)
!!$        print*, 'out', output_data(:,nZ,nH/2,var)
!!$     end if
!!$  end do

  ! Remultiply by density for density and PAH source terms
  do var = 1,nvar_out
     if (trim(input_name(var)).eq.'RhoDot') output_data(:,:,:,var) = output_data(:,:,:,var) * output_data(:,:,:,var_density)
     if (trim(input_name(var)).eq.'Dimer_ProdRate') output_data(:,:,:,var) = output_data(:,:,:,var) * output_data(:,:,:,var_density)
     if (trim(input_name(var)).eq.'ProdRatePos-PAH') output_data(:,:,:,var) = output_data(:,:,:,var) * output_data(:,:,:,var_density)
     if (trim(input_name(var)).eq.'ProdRateNeg-PAH') output_data(:,:,:,var) = output_data(:,:,:,var) * output_data(:,:,:,var_density)
  end do

  ! Compute dT/dC
  do var = 1,nvar_out
     if (trim(input_name(var)).eq.'temperature') var_temp = var
  end do
  do var = 1,nvar_out
     if (trim(input_name(var)).eq.'dTdC') then
           do ih = 1,nH
              do iz = 1,nZ
                 output_data(1,iz,ih,var) = (output_data(2,iz,ih,var_temp) - output_data(2,iz,ih,var_temp)) / (Cdex(2)-Cdex(1))
                 do ic = 2,nC-1
                    output_data(ic,iz,ih,var) = (output_data(ic+1,iz,ih,var_temp) - output_data(ic-1,iz,ih,var_temp)) / (Cdex(ic+1)-Cdex(ic-1))
                 end do
                 output_data(nC,iz,ih,var) = (output_data(nC,iz,ih,var_temp) - output_data(nC-1,iz,ih,var_temp)) / (Cdex(nC)-Cdex(nC-1))
                 do ic = 2,nC
                    if (output_data(ic,iz,ih,var).eq.0.0_WP) then
                       do i = 1,ic-1
                          if (output_data(i,iz,ih,var).ne.0.0_WP) output_data(ic,iz,ih,var) = output_data(i,iz,ih,var)
                       end do
                    end if
                 end do
              end do
           end do
      end if
  end do

!!$  if (nZ.eq.1) then
!!$     do var = 1,nvar_out
!!$        do ih = 1,nH
!!$           do ic = 1,nC
!!$              output_data(ic,2,ih,var) = output_data(ic,1,ih,var)
!!$           end do
!!$        end do
!!$     end do
!!$     nZ = 2
!!$  end if

  return
end subroutine unsteady_prem_table_setup


! ============================================================= !
! Extent the value of the chemtable outside the physical bounds !
! ============================================================= !
subroutine unsteady_prem_table_extent
  use unsteady_prem_table
  implicit none
  
  integer  :: ic,i3,i4
  integer  :: i,j,imin,jmin
  real(WP) :: d,dmin
  
  ! Loop over the three mapping directions
!!$  do i4=1,n4
!!$     do i3=1,n3
!!$        do izv=1,nZVar
!!$           do izm=1,nZMean
!!$           
!!$              ! If masked recompute value from nearest neighboor
!!$              if (mask(ic,i3,i4).eq.1) then
!!$
!!$                 dmin = huge(WP)
!!$                 do i=1,nZMean
!!$                    do j=1,nZVar
!!$                       if (mask(i,j,i3).eq.0) then
!!$                          d = sqrt((ZMean(izm)-ZMean(i))**2+(ZVar(izv)-ZVar(j))**2)
!!$                          if (d.lt.dmin) then
!!$                             imin = i
!!$                             jmin = j
!!$                             dmin = d
!!$                          end if
!!$                       end if
!!$                    end do
!!$                 end do
!!$                 
!!$                 output_data(ic,i3,:) = output_data(imin,jmin,i3,:) 
!!$
!!$              end if
!!$           
!!$           end do
!!$        end do
!!$     end do
!!$  end do
  
  return
end subroutine unsteady_prem_table_extent



! ===================== !
! Print some statistics !
! ===================== !
subroutine unsteady_prem_table_stats
  use unsteady_prem_table
  implicit none
  
  integer :: var
  
  print*,''
  
  ! Min and Max of the coordinates
  print*, '** Coordinates of the table **'
  write(*,10) 'Variable    ','Min    ','Max    '
  write(*,10) '------------','------------','------------'
  write(*,11) 'C           ', minval(Cdex), maxval(Cdex)
  write(*,11) 'Z           ', minval(Zdex), maxval(Zdex)
  write(*,11) input_name(nvar_in), minval(Hdex), maxval(Hdex)
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
end subroutine unsteady_prem_table_stats


! ===================================== !
! Write the table back to a binary file !
! ===================================== !
subroutine unsteady_prem_table_write
  use unsteady_prem_table
  use parser
  use string
  use fileio
  implicit none
  
  logical :: CTprint
  integer :: ierr,var,iunit,i,j,k,l,iunit20
  character(len=str_medium) :: filename,tmp_str,filename4
  real(WP) :: tmp_val

  call parser_read('Print Chemtable', CTprint, .false.)
  if (CTprint) then
     call parser_read('Chemtable print name', filename4)
  end if
  
  ! Open the data file
  call parser_read('Table filename', filename)
  call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierr)

!!$  ! Header for VIDA table
!!$  if (vida_table) then
!!$     tmp_val = 0.0_WP
!!$     call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
!!$     call BINARY_FILE_WRITE(iunit,combModel,str_medium,kind(combModel),ierr)
!!$     tmp_str = 'SPECIES'
!!$     call BINARY_FILE_WRITE(iunit,tmp_str,str_medium,kind(tmp_str),ierr)
!!$     call parser_read('Reference pressure',tmp_val)
!!$     call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
!!$     call parser_read('Oxidizer temperature',tmp_val)
!!$     call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
!!$     call parser_read('Fuel temperature',tmp_val)
!!$     call BINARY_FILE_WRITE(iunit,tmp_val,1,kind(tmp_val),ierr)
!!$  end if
  
  ! Write sizes
  call BINARY_FILE_WRITE(iunit,nC,1,kind(nC),ierr)
  call BINARY_FILE_WRITE(iunit,nZ,1,kind(nZ),ierr)
  call BINARY_FILE_WRITE(iunit,nH,1,kind(nH),ierr)
  !call BINARY_FILE_WRITE(iunit,n4,1,kind(n4),ierr)
  if (masking) then
     call BINARY_FILE_WRITE(iunit,nvar_out+1,1,kind(nvar_out),ierr)
  else
     call BINARY_FILE_WRITE(iunit,nvar_out,1,kind(nvar_out),ierr)
  end if
  ! Write the axis coordinates
  call BINARY_FILE_WRITE(iunit,Cdex,nC,kind(Cdex),ierr)
  call BINARY_FILE_WRITE(iunit,Zdex,nZ,kind(Zdex),ierr)
  call BINARY_FILE_WRITE(iunit,Hdex,nH,kind(Hdex),ierr)
  !call BINARY_FILE_WRITE(iunit,Z4,n4,kind(Z4),ierr)
  if (.not.vida_table) then
     ! Write additional stuff
     call BINARY_FILE_WRITE(iunit,combModel,str_medium,kind(combModel),ierr)
  end if
  ! Write variable names
  do var=1,nvar_out
     call BINARY_FILE_WRITE(iunit,output_name(var),str_medium,kind(output_name),ierr)
  end do
  ! Write mask
  if (masking) then
     call BINARY_FILE_WRITE(iunit,mask_name,str_medium,kind(mask_name),ierr)
  end if
  if (.not.vida_table) then
     ! Write data field: Fortran -> column-major
     do var=1,nvar_out
        call BINARY_FILE_WRITE(iunit,output_data(:,:,:,var),nC*nZ*nH,kind(output_data),ierr)
     end do
     if (masking) then
        call BINARY_FILE_WRITE(iunit,mask(:,:,:),nC*nZ*nH,kind(mask),ierr)
     end if
  else
     ! Write data field: C++ -> row-major
     do i=1,nC
        do j=1,nZ
           do k=1,nH
              do var=1,nvar_out
                 call BINARY_FILE_WRITE(iunit,output_data(i,j,k,var),1,kind(output_data),ierr)
              end do
           end do
        end do
     end do
  end if
  
  call BINARY_FILE_CLOSE(iunit,ierr)

!!$  ! Write a 2D slice -- Zero heat loss
!!$  call BINARY_FILE_OPEN(iunit,'chemtable.3Dslice',"w",ierr)
!!$  call BINARY_FILE_WRITE(iunit,nC,1,kind(nC),ierr)
!!$  call BINARY_FILE_WRITE(iunit,nZ,1,kind(nZ),ierr)
!!$  !call BINARY_FILE_WRITE(iunit,n3,1,kind(n3),ierr)
!!$  call BINARY_FILE_WRITE(iunit,nvar_out,1,kind(nvar_out),ierr)
!!$  call BINARY_FILE_WRITE(iunit,Cdex,nC,kind(Cdex),ierr)
!!$  call BINARY_FILE_WRITE(iunit,Zdex,nZ,kind(Zdex),ierr)
!!$  !call BINARY_FILE_WRITE(iunit,Z3,n3,kind(Z3),ierr)
!!$  call BINARY_FILE_WRITE(iunit,combModel,str_medium,kind(combModel),ierr)
!!$  do var=1,nvar_out
!!$     call BINARY_FILE_WRITE(iunit,output_name(var),str_medium,kind(output_name),ierr)
!!$  end do
!!$  if (masking) then
!!$     call BINARY_FILE_WRITE(iunit,mask_name,str_medium,kind(mask_name),ierr)
!!$  end if
!!$  do var=1,nvar_out
!!$     call BINARY_FILE_WRITE(iunit,output_data(:,:,nH,var),nC*nZ,kind(output_data),ierr)
!!$  end do
!!$  if (masking) then
!!$     call BINARY_FILE_WRITE(iunit,mask(:,:,nH),nC*nZ,kind(mask),ierr)
!!$  end if
!!$  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Print the Chemtable for processing, if you check the option
  if (CTprint) then
     open (iunit20, file=filename4)
     do i=1,nC
        do j=1,nZ
           do k=1,nH
              write(iunit20,"(9999(ES15.8,:,','))") CDex(i), ZDex(j), HDex(k), (output_data(i,j,k,l) , l=1,nvar_out), mask(i,j,k)
           end do
        end do
     end do
     close(iclose(iunit20))
  end if
  
  return
end subroutine unsteady_prem_table_write

