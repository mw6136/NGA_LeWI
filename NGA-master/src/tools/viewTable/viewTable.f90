program viewTable
  use table_1D
  use table
  use mmm_table
  use prem_table
  use unsteady_table
  use unsteady_prem_table
  use parser
  implicit none

  ! =============================================================================== !
  ! program viewTable summarizes properties of a binary file containing tabulated   !
  ! thermochemical states. viewTable is useful when one does not have access to the !
  ! original input file or flamelet files used to create the table.                 !
  !                                                                                 ! 
  ! viewTable provides:                                                             !
  !    - the resolution of each dimension in the table                              !
  !    - the number of variables stored in the table                                !
  !    - an option to save the coordinates of each dimension in the table to a file !
  !    - the combustion model used to set up the table                              !
  !    - a list of all variables stored in the table                                ! 
  !    - the minimum and maximum values for each variable in the table              !
  !    - a check for whether any element is NaN                                     !
  !    - an option to create a 3D slice of the table (for 4D tables)                !
  ! =============================================================================== !

  character(str_long) :: filename
  character(str_medium) :: combustModel, combustModelslice
  character(str_short) :: coords,vida,chem3Dlastaxis,slice

  character(len=str_medium), dimension(:), pointer :: var_name

  integer :: tabtype,maxtypes,slicevar,iunit,ierr
  integer :: maxdim,nvars
  integer :: iter,iter1,iter2,iter3,iter4

  integer, dimension(4) :: dim, itervec, itervec2, reorder
  
  real(WP) :: val1,val2,val3,val4,sliceval,alpha,interpval
  real(WP) :: tol = 1e-7
  
  
  real(WP), dimension(:,:),       pointer :: axis
  real(WP), dimension(:,:),       pointer :: var1_data
  real(WP), dimension(:,:,:),     pointer :: var2_data
  real(WP), dimension(:,:,:,:),   pointer :: var3_data
  real(WP), dimension(:,:,:,:,:), pointer :: var4_data
  
  logical :: nanflag

  print *, " Please run viewTable in the same directory as the binary table file. "
  print *, ""
  print *, " If you input a path to the table file, there may be too many characters "
  print *, " and you may receive a segmentation fault with an error message stating "
  print *, " BINARY_FILE_OPEN : Error openning the file. "
  print *, ""

  ! Read table file name and type from standard input
  print "(a29,$)", " Name of binary table file : "
  read "(a)", filename
  print *, " Select table type : "
  print *, "    1) Chemtable, 1D with ZMEAN "
  print *, "    2) Chemtable, 3D with ZMEAN, ZVAR, PROG or CHI or ENTH "
  print *, "    3) Chemtable, 4D with ZMEAN, ZVAR, PROG, ENTH "
  print *, "    4) MMM Table, 4D with ZMEAN, ZVAR, FMIX, PROG "
  print *, "    5) Premtable, 2D with PROG, ZMEAN "
  print *, "    6) Premtable, 3D with PROG, ZMEAN, ENTH "
  print "(a18,$)", " Your selection : "
  read "(i1)", tabtype

  ! Number of table types currently supported
  maxtypes = 6
  if (tabtype .lt. 1 .or. tabtype .gt. maxtypes) then
     print *, "Error: Selection must be an integer in the range from 1 to", maxtypes
     stop
  end if

  ! Open the table file to read
  call BINARY_FILE_OPEN(iunit,trim(filename),"r",ierr)

  dim = [0,0,0,0]
  ! Read size of each dimension of table
  select case(tabtype)
     case (1)
        call BINARY_FILE_READ(iunit,dim(1),1,kind(dim(1)),ierr)
        print *, "Size of ZMEAN :", dim(1)
     case (2)
        call BINARY_FILE_READ(iunit,dim(1),1,kind(dim(1)),ierr)
        print *, "Size of ZMEAN :", dim(1)
        call BINARY_FILE_READ(iunit,dim(2),1,kind(dim(2)),ierr)
        print *, "Size of ZVAR :", dim(2)
        call BINARY_FILE_READ(iunit,dim(3),1,kind(dim(3)),ierr)
        print *, "Size of PROG/CHI/ENTH :", dim(3)
     case (3)
        call BINARY_FILE_READ(iunit,dim(1),1,kind(dim(1)),ierr)
        print *, "Size of ZMEAN :", dim(1)
        call BINARY_FILE_READ(iunit,dim(2),1,kind(dim(2)),ierr)
        print *, "Size of ZVAR :", dim(2)
        call BINARY_FILE_READ(iunit,dim(3),1,kind(dim(3)),ierr)
        print *, "Size of PROG :", dim(3)
        call BINARY_FILE_READ(iunit,dim(4),1,kind(dim(4)),ierr)
        print *, "Size of ENTH :", dim(4)
     case (4)
        call BINARY_FILE_READ(iunit,dim(1),1,kind(dim(1)),ierr)
        print *, "Size of ZMEAN :", dim(1)
        call BINARY_FILE_READ(iunit,dim(2),1,kind(dim(2)),ierr)
        print *, "Size of ZVAR :", dim(2)
        call BINARY_FILE_READ(iunit,dim(3),1,kind(dim(3)),ierr)
        print *, "Size of FMIX :", dim(3)
        call BINARY_FILE_READ(iunit,dim(4),1,kind(dim(4)),ierr)
        print *, "Size of PROG :", dim(4)
     case (5)
        call BINARY_FILE_READ(iunit,dim(1),1,kind(dim(1)),ierr)
        print *, "Size of PROG :", dim(1)
        call BINARY_FILE_READ(iunit,dim(2),1,kind(dim(2)),ierr)
        print *, "Size of ZMEAN :", dim(2)
     case (6)
        call BINARY_FILE_READ(iunit,dim(1),1,kind(dim(1)),ierr)
        print *, "Size of PROG :", dim(1)
        call BINARY_FILE_READ(iunit,dim(2),1,kind(dim(2)),ierr)
        print *, "Size of ZMEAN :", dim(2)
        call BINARY_FILE_READ(iunit,dim(3),1,kind(dim(3)),ierr)
        print *, "Size of ENTH :", dim(3)
     case default
        stop "Error: Unsupported table type selected"
     end select

     ! Read number of variables stored in table
     call BINARY_FILE_READ(iunit,nvars,1,kind(nvars),ierr)
     print *, "Number of variables in table :", nvars
     maxdim = max(dim(1),dim(2),dim(3),dim(4))
     allocate(axis(4,maxdim))
     
     ! Option to save axis coordinates to a file
     print "(a44,$)", " Save the axis coordinates to a file? (y/n) "
     read "(a)", coords
     if (trim(adjustl(coords)) .ne. 'y' .and. trim(adjustl(coords)) .ne. 'n') then
        print *, "Error: Please enter y or n"
        stop
     end if

     select case(tabtype)
     case (1)
        call BINARY_FILE_READ(iunit,axis(1,:),dim(1),kind(axis),ierr)
        if (coords .eq. 'y') then
           open(unit=2000,file="axiscoords.txt")
           write(2000,'(a22)') "ZMEAN"
           do iter1=1,dim(1)
              write(2000,'(es22.15)') axis(1,iter1) 
           end do
           close(2000)
        end if
     case (2)
        call BINARY_FILE_READ(iunit,axis(1,:),dim(1),kind(axis),ierr)
        call BINARY_FILE_READ(iunit,axis(2,:),dim(2),kind(axis),ierr)
        call BINARY_FILE_READ(iunit,axis(3,:),dim(3),kind(axis),ierr)
        ! Delay writing axis coordinates until combustion model determined
     case (3)
        call BINARY_FILE_READ(iunit,axis(1,:),dim(1),kind(axis),ierr)
        call BINARY_FILE_READ(iunit,axis(2,:),dim(2),kind(axis),ierr)
        call BINARY_FILE_READ(iunit,axis(3,:),dim(3),kind(axis),ierr)
        call BINARY_FILE_READ(iunit,axis(4,:),dim(4),kind(axis),ierr)
        if (coords .eq. 'y') then
           open(unit=2002,file="axiscoords.txt")
           write(2002,'(4(a22,1x))') "ZMEAN", "ZVAR", "PROG", "ENTH"
           do iter=1,maxdim
              ! Account for axes with differing sizes
              if (iter .le. dim(1)) then
                 val1 = axis(1,iter)
              else
                 ! Fill in space with 0.0
                 val1 = 0.0_WP
              end if
              if (iter .le. dim(2)) then
                 val2 = axis(2,iter)
              else
                 ! Fill in space with 0.0
                 val2 = 0.0_WP
              end if
              if (iter .le. dim(3)) then
                 val3 = axis(3,iter)
              else
                 ! Fill in space with 0.0
                 val3 = 0.0_WP
              end if
              if (iter .le. dim(4)) then
                 val4 = axis(4,iter)
              else
                 ! Fill in space with 0.0
                 val4 = 0.0_WP
              end if
              write(2002,'(4(es22.15,1x))') val1, val2, val3, val4
           end do
           close(2002)
        end if
     case (4)
        call BINARY_FILE_READ(iunit,axis(1,:),dim(1),kind(axis),ierr)
        call BINARY_FILE_READ(iunit,axis(2,:),dim(2),kind(axis),ierr)
        call BINARY_FILE_READ(iunit,axis(3,:),dim(3),kind(axis),ierr)
        call BINARY_FILE_READ(iunit,axis(4,:),dim(4),kind(axis),ierr)
        if (coords .eq. 'y') then
           open(unit=2003,file="axiscoords.txt")
           write(2003,'(4(a22,1x))') "ZMEAN", "ZVAR", "FMIX", "PROG"
           do iter=1,maxdim
              ! Account for axes with differing sizes
              if (iter .le. dim(1)) then
                 val1 = axis(1,iter)
              else
                 ! Fill in space with 0.0
                 val1 = 0.0_WP
              end if
              if (iter .le. dim(2)) then
                 val2 = axis(2,iter)
              else
                 ! Fill in space with 0.0
                 val2 = 0.0_WP
              end if
              if (iter .le. dim(3)) then
                 val3 = axis(3,iter)
              else
                 ! Fill in space with 0.0
                 val3 = 0.0_WP
              end if
              if (iter .le. dim(4)) then
                 val4 = axis(4,iter)
              else
                 ! Fill in space with 0.0
                 val4 = 0.0_WP
              end if
              write(2003,'(4(es22.15,1x))') val1, val2, val3, val4
           end do
           close(2003)
        end if
     case (5)
        call BINARY_FILE_READ(iunit,axis(1,:),dim(1),kind(axis),ierr)
        call BINARY_FILE_READ(iunit,axis(2,:),dim(2),kind(axis),ierr)
        if (coords .eq. 'y') then
           open(unit=2004,file="axiscoords.txt")
           write(2004,'(2(a22,1x))') "PROG","ZMEAN"
           do iter=1,maxdim
              ! Account for axes with differing sizes
              if (iter .le. dim(1)) then
                 val1 = axis(1,iter)
              else
                 ! Fill in space with 0.0
                 val1 = 0.0_WP
              end if
              if (iter .le. dim(2)) then
                 val2 = axis(2,iter)
              else
                 ! Fill in space with 0.0
                 val2 = 0.0_WP
              end if
              write(2004,'(2(es22.15,1x))') val1, val2
           end do
           close(2004)
        end if
     case (6)
        call BINARY_FILE_READ(iunit,axis(1,:),dim(1),kind(axis),ierr)
        call BINARY_FILE_READ(iunit,axis(2,:),dim(2),kind(axis),ierr)
        call BINARY_FILE_READ(iunit,axis(3,:),dim(3),kind(axis),ierr)
        if (coords .eq. 'y') then
           open(unit=2005,file="axiscoords.txt")
           write(2005,'(3(a22,1x))') "PROG", "ZMEAN", "ENTH"
           maxdim = max(dim(1),dim(2),dim(3))
           do iter=1,maxdim
              ! Account for axes with differing sizes
              if (iter .le. dim(1)) then
                 val1 = axis(1,iter)
              else
                 ! Fill in space with 0.0
                 val1 = 0.0_WP
              end if
              if (iter .le. dim(2)) then
                 val2 = axis(2,iter)
              else
                 ! Fill in space with 0.0
                 val2 = 0.0_WP
              end if
              if (iter .le. dim(3)) then
                 val3 = axis(3,iter)
              else
                 ! Fill in space with 0.0
                 val3 = 0.0_WP
              end if
              write(2005,'(3(es22.15,1x))') val1, val2, val3
           end do
           close(2005)
        end if
     case default
        stop "Error: Unsupported table type selected"
     end select

     ! Read combustion model
     print "(a29,$)", " Is this a VIDA table? (y/n) "
     read "(a)", vida
     if (trim(adjustl(vida)) .ne. 'y' .and. trim(adjustl(vida)) .ne. 'n') then
        print *, "Error: Please enter y or n"
        stop
     end if

     if (trim(adjustl(vida)) .eq. 'n') then
        call BINARY_FILE_READ(iunit,combustModel,str_medium,kind(combustModel),ierr)
        print *, "Combustion Model : ", trim(adjustl(combustModel))
     end if
     
     ! Special case for writing axis coordinates
     if (tabtype .eq. 2) then
        ! Determine name of last axis of table based on combustion model
        if (trim(adjustl(combustModel)) .eq. "FPVA") then
           chem3Dlastaxis = "PROG"
        elseif (trim(adjustl(combustModel)) .eq. "Steady Flamelet") then
           chem3Dlastaxis = "CHI"
        elseif (trim(adjustl(combustModel)) .eq. "Enthalpy Flamelet") then
           chem3Dlastaxis = "ENTH"
        else
           stop "Error: Unsupported combustion model"
        end if

        if (coords .eq. 'y') then
           open(unit=2001,file="axiscoords.txt")   
           write(2001,'(3(a22,1x))') "ZMEAN", "ZVAR", trim(adjustl(chem3Dlastaxis))
           do iter=1,maxdim
              ! Account for axes with differing sizes
              if (iter .le. dim(1)) then
                 val1 = axis(1,iter)
              else
                 ! Fill in space with 0.0
                 val1 = 0.0_WP
              end if
              if (iter .le. dim(2)) then
                 val2 = axis(2,iter)
              else
                 ! Fill in space with 0.0
                 val2 = 0.0_WP
              end if
              if (iter .le. dim(3)) then
                 val3 = axis(3,iter)
              else
                 ! Fill in space with 0.0
                 val3 = 0.0_WP
              end if
              write(2001,'(3(es22.15,1x))') val1, val2, val3
           end do
           close(2001)
        end if
     end if
     
     ! Read variable names
     allocate(var_name(nvars))
     do iter=1,nvars
        call BINARY_FILE_READ(iunit,var_name(iter),str_medium,kind(var_name),ierr)
     end do

     ! Initialize NaN indicator
     nanflag = .false.

     ! Loop through all dimensions of the table for each variable and find the 
     ! minimum and maximum values. Indicate if any element is NaN.
     print *, "Variables in table :"
     select case(tabtype)
     case (1)
        allocate(var1_data(dim(1),nvars))
        do iter=1,nvars
           call BINARY_FILE_READ(iunit,var1_data(:,iter),dim(1),kind(var1_data),ierr)
           print 1000, trim(adjustl(var_name(iter))), " -> Min:", &
                minval(var1_data(:,iter)), "  Max:", maxval(var1_data(:,iter))
           1000 format (a20,a8,g22.15,a6,g22.15)
           do iter1=1,dim(1)
              nanflag = nanflag .or. isnan(var1_data(iter1,iter))
              if (nanflag) then
                 print *, "An element in the chemtable is NaN. Its location is"
                 print *, "ZMEAN index :", iter1
                 stop
              end if
           end do
        end do
     case (2)
        allocate(var3_data(dim(1),dim(2),dim(3),nvars))
        do iter=1,nvars
           call BINARY_FILE_READ(iunit,var3_data(:,:,:,iter), &
                dim(1)*dim(2)*dim(3),kind(var3_data),ierr)
           print 1001, trim(adjustl(var_name(iter))), " -> Min:", &
                minval(var3_data(:,:,:,iter)), "  Max:", maxval(var3_data(:,:,:,iter))
           1001 format (a20,a8,g22.15,a6,g22.15)
           do iter1=1,dim(1)
              do iter2=1,dim(2)
                 do iter3=1,dim(3)
                    nanflag = nanflag .or. isnan(var3_data(iter1,iter2,iter3,iter))
                    if (nanflag) then
                       print *, "An element in the chemtable is NaN. Its location is"
                       print *, "ZMEAN index :", iter1
                       print *, "ZVAR index :", iter2
                       print *, trim(adjustl(chem3Dlastaxis))," index :", iter3
                       stop
                    end if
                 end do
              end do
           end do
        end do
     case (3)
        allocate(var4_data(dim(1),dim(2),dim(3),dim(4),nvars))
        do iter=1,nvars
           call BINARY_FILE_READ(iunit,var4_data(:,:,:,:,iter), &
                dim(1)*dim(2)*dim(3)*dim(4),kind(var4_data),ierr)
           print 1002, trim(adjustl(var_name(iter))), " -> Min:", &
                minval(var4_data(:,:,:,:,iter)), "  Max:", maxval(var4_data(:,:,:,:,iter))
           1002 format (a20,a8,g22.15,a6,g22.15)
           do iter1=1,dim(1)
              do iter2=1,dim(2)
                 do iter3=1,dim(3)
                    do iter4=1,dim(4)
                       nanflag = nanflag .or. isnan(var4_data(iter1,iter2,iter3,iter4,iter))
                       if (nanflag) then
                          print *, "An element in the chemtable is NaN. Its location is"
                          print *, "ZMEAN index :", iter1
                          print *, "ZVAR index :", iter2
                          print *, "PROG index :", iter3
                          print *, "ENTH index :", iter4
                          stop
                       end if
                    end do
                 end do
              end do
           end do
        end do
     case (4)
        allocate(var4_data(dim(1),dim(2),dim(3),dim(4),nvars))
        do iter=1,nvars
           call BINARY_FILE_READ(iunit,var4_data(:,:,:,:,iter), &
                dim(1)*dim(2)*dim(3)*dim(4),kind(var4_data),ierr)
           print 1003, trim(adjustl(var_name(iter))), " -> Min:", &
                minval(var4_data(:,:,:,:,iter)), "  Max:", maxval(var4_data(:,:,:,:,iter))
           1003 format (a20,a8,g22.15,a6,g22.15)
           do iter1=1,dim(1)
              do iter2=1,dim(2)
                 do iter3=1,dim(3)
                    do iter4=1,dim(4)
                       nanflag = nanflag .or. isnan(var4_data(iter1,iter2,iter3,iter4,iter))
                       if (nanflag) then
                          print *, "An element in the chemtable is NaN. Its location is"
                          print *, "ZMEAN index :", iter1
                          print *, "ZVAR index :", iter2
                          print *, "FMIX index :", iter3
                          print *, "PROG index :", iter4
                          stop
                       end if
                    end do
                 end do
              end do
           end do
        end do
     case (5)
        allocate(var2_data(dim(1),dim(2),nvars))
        do iter=1,nvars
           call BINARY_FILE_READ(iunit,var2_data(:,:,iter), &
                dim(1)*dim(2),kind(var2_data),ierr)
           print 1004, trim(adjustl(var_name(iter))), " -> Min:", &
                minval(var2_data(:,:,iter)), "  Max:", maxval(var2_data(:,:,iter))
           1004 format (a20,a8,g22.15,a6,g22.15)
           do iter1=1,dim(1)
              do iter2=1,dim(2)
                 nanflag = nanflag .or. isnan(var2_data(iter1,iter2,iter))
                 if (nanflag) then
                    print *, "An element in the chemtable is NaN. Its location is"
                    print *, "PROG index :", iter1
                    print *, "ZMEAN index :", iter2
                    stop
                 end if
              end do
           end do
        end do
     case (6)
        allocate(var3_data(dim(1),dim(2),dim(3),nvars))
        do iter=1,nvars
           call BINARY_FILE_READ(iunit,var3_data(:,:,:,iter), &
                dim(1)*dim(2)*dim(3),kind(var3_data),ierr)
           print 1005, trim(adjustl(var_name(iter))), " -> Min:", &
                minval(var3_data(:,:,:,iter)), "  Max:", maxval(var3_data(:,:,:,iter))
           1005 format (a20,a8,g22.15,a6,g22.15)
           do iter1=1,dim(1)
              do iter2=1,dim(2)
                 do iter3=1,dim(3)
                    nanflag = nanflag .or. isnan(var3_data(iter1,iter2,iter3,iter))
                    if (nanflag) then
                       print *, "An element in the chemtable is NaN. Its location is"
                       print *, "PROG index :", iter1
                       print *, "ZMEAN index :", iter2
                       print *, "ENTH index :", iter3
                       stop
                    end if
                 end do
              end do
           end do
        end do
     case default
        stop "Error: Unsupported table type selected"
     end select

     if (.not. nanflag) print *, "No elements in the table are NaN."  
     call BINARY_FILE_CLOSE(iunit,ierr)

     ! If table is 4D, give option to create a 3D slice
     if (tabtype .eq. 4) then
        print "(a26,$)", " Create a 3D slice? (y/n) "
        read "(a)", slice
        if (trim(adjustl(slice)) .ne. 'y' .and. trim(adjustl(slice)) .ne. 'n') then
           print *, "Error: Please enter y or n"
           stop
        end if
        
        print "(a40,$)", " Combustion Model for Slice (Optional) :"
        read "(a)", combustModelslice

        if (trim(adjustl(slice)) .eq. 'y') then

           ! Get axis to slice from user
           print "(a45,$)", "Slice at constant value of which axis? (1-4) "
           read "(i1)", slicevar
           if (slicevar.eq.1) then
              reorder=[2,3,4,1]
           else if (slicevar.eq.2) then
              reorder=[1,3,4,2]
           else if (slicevar.eq.3) then
              reorder=[1,2,4,3]
           else if (slicevar.eq.4) then
              reorder=[1,2,3,4]
           else
              print *, "Error: must chose axis 1, 2, 3, or 4"
              stop
           end if

           ! Get value to slice at and ensure it is within bounds
           print "(a21,$)", "Slice at what value? "
           read "(F)", sliceval
           itervec(slicevar)  = 0
           itervec2(slicevar) = 0
           if (sliceval.lt.(axis(slicevar,1)-tol) .or. sliceval.gt.(axis(slicevar,dim(slicevar))+tol)) then
              print *, "Error: slice value out of bounds for axis", slicevar
              print *, "min = ", axis(slicevar,1), " max = ", axis(slicevar,dim(slicevar))
              stop
           else if (axis(slicevar,1) .gt. axis(slicevar,dim(slicevar))) then
              print *, "Error: slice axis must be in increasing order"
           end if
              
           ! Find interpolation for slice value on slice axis
           do iter=1,dim(slicevar)
              if (sliceval.lt.(axis(slicevar,iter)+tol) .and. itervec2(slicevar).eq.0) itervec2(slicevar) = iter
              if (sliceval.ge.(axis(slicevar,iter)-tol)                              ) itervec(slicevar)  = iter
           end do
           if (itervec2(slicevar).eq.itervec(slicevar)) then
              alpha = 1.0_WP
           else
              alpha = (sliceval - axis(slicevar,itervec(slicevar)))/(axis(slicevar,itervec2(slicevar))-axis(slicevar,itervec(slicevar)))
           end if
           print *, "interp", sliceval, itervec(slicevar), itervec2(slicevar), alpha
                      
           ! Open file for chemtable slice
           call BINARY_FILE_OPEN(iunit,trim(filename)//'.3Dslice',"w",ierr)

           ! Write size of each dimension of slice
           do iter = 1,3
              call BINARY_FILE_WRITE(iunit,dim(reorder(iter)),1,kind(dim(reorder(iter))),ierr)
           end do

           ! Write number of variables stored in slice
           call BINARY_FILE_WRITE(iunit,nvars,1,kind(nvars),ierr)

           ! Write axis coordinates
           do iter = 1,3
              call BINARY_FILE_WRITE(iunit,axis(reorder(iter),1:dim(reorder(iter))),dim(reorder(iter)),kind(axis),ierr)
           end do

           ! Write combustion model
           if (trim(adjustl(vida)) .eq. 'n') then
              call BINARY_FILE_WRITE(iunit,combustModelslice,str_medium,kind(combustModelslice),ierr)
           end if

           ! Write variable names
           do iter=1,nvars
              call BINARY_FILE_WRITE(iunit,var_name(iter),str_medium,kind(var_name),ierr)
           end do

           ! Interpolate and write data to file
           do iter = 1,nvars
              do iter3 = 1,dim(reorder(3)) 
                 do iter2 = 1,dim(reorder(2))
                    do iter1 = 1,dim(reorder(1))
                       itervec2(reorder(1)) = iter1; itervec(reorder(1)) = iter1
                       itervec2(reorder(2)) = iter2; itervec(reorder(2)) = iter2
                       itervec2(reorder(3)) = iter3; itervec(reorder(3)) = iter3
                       interpval   =  alpha  * var4_data(itervec2(1), itervec2(2), itervec2(3), itervec2(4), iter) &
                            + (1.0_WP-alpha) * var4_data(itervec(1),  itervec(2),  itervec(3),  itervec(4),  iter)
                       call BINARY_FILE_WRITE(iunit,interpval,1,kind(interpval),ierr)
                    end do
                 end do
              end do
           end do
           
           call BINARY_FILE_CLOSE(iunit,ierr)
           
        end if
     end if

end program viewTable
