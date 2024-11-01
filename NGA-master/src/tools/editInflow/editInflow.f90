program editInflow
  use precision
  use string
  use fileio
  implicit none

  ! Data array with all the variables (velocity,pressure,...)
  integer :: ntime,ny,nz,nvar,nvar2,n_del,icyl
  character(len=str_short), dimension(:), pointer :: names,names2
  character(len=str_short), dimension(10) :: names_del
  real(WP), dimension(:,:), pointer :: data
  real(WP), dimension(:), pointer :: y,z
  integer :: iunit1,iunit2,ierr,var,choice,iunit3,var2,itime
  character(len=str_medium) :: filename1,filename2,filename3
  character(len=str_short) :: varname,varname2
  real(WP) :: dt,time,value

  ! Read file name from standard input
  print*,'========================'
  print*,'| ARTS - inflow Editor |'
  print*,'========================'
  print*
  print "(a30,$)", " inflow file before edition : "
  read "(a)", filename1
  print "(a29,$)", " inflow file after edition : "
  read "(a)", filename2

  ! ** Open the data file to read **
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit1,ntime,1,kind(ntime),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
  print*,'Grid :',ntime,'ntime',ny,'ny',nz,'nz'
 
  ! Read additional stuff
  call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  print*,'Inflow file at time :',time
  
  ! Read variable names
  allocate(names(nvar))
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
  end do
  print*,'Variables : ',names
  print*,'There are ',nvar,' variables.'

  ! Allocate arrays, read spatial coordinates
  allocate(data(ny,nz))
  allocate(y(ny+1))
  allocate(z(nz+1))
  call BINARY_FILE_READ(iunit1,icyl,1,kind(icyl),ierr)
  call BINARY_FILE_READ(iunit1,y,ny+1,kind(y),ierr)
  call BINARY_FILE_READ(iunit1,z,nz+1,kind(z),ierr)

  ! ** Ask what to do **
  print*
  print*, "1. Print Min/Max of variable"
  print*, "2. Delete variable"
  print*, "3. Rename variable" 
  print*, "4. Copy variable to new name"
  print*, "5. Reset simulation time"
  print "(a9,$)", "Choice : "
  read "(i1)", choice
  
  ! Case dependent operation
  select case(choice)
     
  case(1) ! Print min/max of all variables for first and last timesteps and every 1/100th of the way through the file
     do  itime=1,ntime
        if (itime.eq.1) then
           print *, "first timestep:"
        else if (itime.eq.ntime) then
           print *, "last timestep:"
        else if (modulo(itime,int(ntime/100_WP)).eq.0) then
           print *, "timestep", itime, "/", ntime, ":"
        end if
        do var=1,nvar
           call BINARY_FILE_READ(iunit1,data(:,:),ny*nz,kind(data),ierr)
           if (ISNAN(data(1,1))) then ! Check for NaNs as well
              print *, "WARNING: NAN at timestep", itime, "for variable", names(var)
           else if ((itime.eq.1).or.(itime.eq.ntime).or.(modulo(itime,int(ntime/100_WP)).eq.0)) then
              print*, names(var)
              print*,"min: ",minval(data)," - max: ",maxval(data)
              print*,"maxloc: ",maxloc(data)
           end if
        end do
        if(itime.eq.1) print*,''
     end do

  case (2) ! Delete variables 
     print*,'You can delete at most 10 variables, press q to exit'
     n_del = 0
     do while (n_del .le. 10)
        print "(a16,$)", "Variable name : "
        read "(a)", varname
        if (varname .ne. 'q') then
           n_del = n_del + 1
           names_del(n_del) = varname
        else
           exit   
        end if
     end do

     print*,'how many are deleted? ',n_del
     print*,'They are :',names_del(1:n_del)

     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,ntime,1,kind(ntime),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar-n_del,1,kind(nvar),ierr)
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)

     do var=1,nvar
        if (.not. any(names_del.eq.(trim(adjustl(names(var)))))) &
             call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
     end do

     call BINARY_FILE_WRITE(iunit2,icyl,1,kind(icyl),ierr)
     call BINARY_FILE_WRITE(iunit2,y,ny+1,kind(y),ierr)
     call BINARY_FILE_WRITE(iunit2,z,nz+1,kind(z),ierr)

     do itime = 1,ntime
        do var=1,nvar
           call BINARY_FILE_READ(iunit1,data,ny*nz,kind(data),ierr)
           if (.not. any(names_del.eq.(trim(adjustl(names(var)))))) &
                call BINARY_FILE_WRITE(iunit2,data,ny*nz,kind(data),ierr)
        end do
     end do
     call BINARY_FILE_CLOSE(iunit2,ierr)
          
  case (3)  ! Rename a variable 
     print "(a16,$)", "Variable name : "
     read "(a)", varname
     print "(a20,$)", "New variable name : "
     read "(a)", varname2

     call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
     call BINARY_FILE_WRITE(iunit2,ntime,1,kind(ntime),ierr)
     call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
     call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
     call BINARY_FILE_WRITE(iunit2,nvar,1,kind(nvar),ierr)
     call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
     call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
     do var=1,nvar
        if (trim(adjustl(names(var))).ne.trim(adjustl(varname))) then
           call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
        else
           call BINARY_FILE_WRITE(iunit2,varname2,str_short,kind(names),ierr)
        end if
     end do

     call BINARY_FILE_WRITE(iunit2,icyl,1,kind(icyl),ierr)
     call BINARY_FILE_WRITE(iunit2,y,ny+1,kind(y),ierr)
     call BINARY_FILE_WRITE(iunit2,z,nz+1,kind(z),ierr)
     do itime = 1,ntime
        do var=1,nvar
           call BINARY_FILE_READ(iunit1,data,ny*nz,kind(data),ierr)
           call BINARY_FILE_WRITE(iunit2,data,ny*nz,kind(data),ierr)
        end do
     end do
     call BINARY_FILE_CLOSE(iunit2,ierr)
      
 case (4)
    print "(a24,$)", "Variable name to copy : "
    read "(a)", varname
    print "(a20,$)", "New variable name : "
    read "(a)", varname2
    call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
    call BINARY_FILE_WRITE(iunit2,ntime,1,kind(ntime),ierr)
    call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
    call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
    call BINARY_FILE_WRITE(iunit2,nvar+1,1,kind(nvar),ierr)
    call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
    call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
    do var=1,nvar
       if (trim(adjustl(names(var))).ne.trim(adjustl(varname))) then
          call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
       else
          call BINARY_FILE_WRITE(iunit2,varname,str_short,kind(names),ierr)
          call BINARY_FILE_WRITE(iunit2,varname2,str_short,kind(names),ierr)
       end if
    end do

    call BINARY_FILE_WRITE(iunit2,icyl,1,kind(icyl),ierr)
    call BINARY_FILE_WRITE(iunit2,y,ny+1,kind(y),ierr)
    call BINARY_FILE_WRITE(iunit2,z,nz+1,kind(z),ierr)
    do itime = 1,ntime
       do var=1,nvar
          call BINARY_FILE_READ(iunit1,data,ny*nz,kind(data),ierr)
          call BINARY_FILE_WRITE(iunit2,data,ny*nz,kind(data),ierr)
          if (trim(adjustl(names(var))).eq.trim(adjustl(varname))) then
             call BINARY_FILE_WRITE(iunit2,data,ny*nz,kind(data),ierr)
          end if
       end do
    end do
    call BINARY_FILE_CLOSE(iunit2,ierr)

 case (5) ! Reset inflow time to zero
    time = 0.0_WP

    call BINARY_FILE_OPEN(iunit2,trim(filename2),"w",ierr)
    call BINARY_FILE_WRITE(iunit2,ntime,1,kind(ntime),ierr)
    call BINARY_FILE_WRITE(iunit2,ny,1,kind(ny),ierr)
    call BINARY_FILE_WRITE(iunit2,nz,1,kind(nz),ierr)
    call BINARY_FILE_WRITE(iunit2,nvar,1,kind(nvar),ierr)
    call BINARY_FILE_WRITE(iunit2,dt,1,kind(dt),ierr)
    call BINARY_FILE_WRITE(iunit2,time,1,kind(time),ierr)
    do var=1,nvar
       call BINARY_FILE_WRITE(iunit2,names(var),str_short,kind(names),ierr)
    end do

    call BINARY_FILE_WRITE(iunit2,icyl,1,kind(icyl),ierr)
    call BINARY_FILE_WRITE(iunit2,y,ny+1,kind(y),ierr)
    call BINARY_FILE_WRITE(iunit2,z,nz+1,kind(z),ierr)
    do itime = 1,ntime
       do var=1,nvar
          call BINARY_FILE_READ(iunit1,data,ny*nz,kind(data),ierr)
          call BINARY_FILE_WRITE(iunit2,data,ny*nz,kind(data),ierr)
       end do
    end do
    call BINARY_FILE_CLOSE(iunit2,ierr)
    
 case default
     stop "Unknown choice"
  end select
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit1,ierr)
  
end program editInflow
