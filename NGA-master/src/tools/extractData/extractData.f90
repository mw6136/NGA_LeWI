! ============================================================ !
!                   extractData.f90                            !
!                                                              !
! Author: Temistocle Grenga                                    !
! Date:   April 2017                                           !
! ============================================================ !

program extData
  use precision
  use string
  use fileio
  implicit none

  integer :: nx, ny, nz, nvar, i, j
  character(len=str_short), dimension(:), pointer :: names
  real(WP), dimension(:,:,:,:), pointer :: data
  integer :: iunit1,iunit2,ierr,var,iunit3,var2
  character(len=str_medium) :: filename1, filename2, fconfig1
  character(len=str_medium) :: config1
  character(len=str_short) :: varname,varname2
  real(WP) :: dt, time, value, minx, maxx, dum1, dum2, DZ, a1, a2
  real(WP) :: miny, maxy, selx, sely, selz
  real(WP), dimension(:), pointer :: x1, y1, z1, newdata
  integer :: minnx, maxnx, selnx, selny, selnz, minny, maxny, sel_case
  integer :: icyl1, xper1, yper1, zper1, nx1, ny1, nz1
  integer, dimension(:,:), pointer :: mask1
  character(15):: frmt

  ! Read file name from standard input
  print*,'======================='
  print*,'| ARTS - Extract Data |'
  print*,'======================='
  print*
  print "(a28,$)", " data file from NGA : "
  read "(a)", filename1
  print "(a27,$)", " data file to write : "
  read "(a)", filename2
  print "(a28,$)", " Enter source config file : "
  read "(a)", fconfig1

  ! Read the source config file
  call BINARY_FILE_OPEN(iunit3,trim(fconfig1),"r",ierr)
  call BINARY_FILE_READ(iunit3,config1,str_medium,kind(config1),ierr)
  call BINARY_FILE_READ(iunit3,icyl1,1,kind(icyl1),ierr)
  call BINARY_FILE_READ(iunit3,xper1,1,kind(xper1),ierr)
  call BINARY_FILE_READ(iunit3,yper1,1,kind(yper1),ierr)
  call BINARY_FILE_READ(iunit3,zper1,1,kind(zper1),ierr)
  call BINARY_FILE_READ(iunit3,nx1,1,kind(nx1),ierr)
  call BINARY_FILE_READ(iunit3,ny1,1,kind(ny1),ierr)
  call BINARY_FILE_READ(iunit3,nz1,1,kind(nz1),ierr)
  print*,'Source file'
  print*,'Config : ',config1
  print*,'Grid :',nx1,'x',ny1,'x',nz1
  allocate(x1(nx1+1),y1(ny1+1),z1(nz1+1))
  allocate(mask1(nx1,ny1))
  call BINARY_FILE_READ(iunit3,x1,nx1+1,kind(x1),ierr)
  call BINARY_FILE_READ(iunit3,y1,ny1+1,kind(y1),ierr)
  call BINARY_FILE_READ(iunit3,z1,nz1+1,kind(z1),ierr)
  call BINARY_FILE_READ(iunit3,mask1,nx1*ny1,kind(mask1),ierr)
  call BINARY_FILE_CLOSE(iunit3,ierr)
  print*,'x -> ', 'min: ',minval(x1(:)),'max : ',maxval(x1(:))
  print*,'y -> ', 'min: ',minval(y1(:)),'max : ',maxval(y1(:)) 
  print*,'z -> ', 'min: ',minval(z1(:)),'max : ',maxval(z1(:))
  call BINARY_FILE_CLOSE(iunit3,ierr)

  ! ** Open the data file to read **
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)

  ! Read sizes
  call BINARY_FILE_READ(iunit1,nx,1,kind(nx),ierr)
  call BINARY_FILE_READ(iunit1,ny,1,kind(ny),ierr)
  call BINARY_FILE_READ(iunit1,nz,1,kind(nz),ierr)
  call BINARY_FILE_READ(iunit1,nvar,1,kind(nvar),ierr)
  print*,'Grid :',nx,'x',ny,'x',nz
    
  ! Read additional stuff
  call BINARY_FILE_READ(iunit1,dt,1,kind(dt),ierr)
  call BINARY_FILE_READ(iunit1,time,1,kind(time),ierr)
  print*,'Timestep size     :',dt
  print*,'Data file at time :',time
  
  ! Read variable names
  allocate(names(nvar))
  print*,'There are ',nvar,' variables.'
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,names(var),str_short,kind(names),ierr)
     print*,'Variables : ', var, ' - ', names(var)
  end do

  ! Allocate arrays
  allocate(data(nx,ny,nz,nvar))

  !print "(a9,$)", "Choice : "
  !read "(i1)", choice
  ! Case dependent operation
  !select case(choice)
     
  do var=1,nvar
     call BINARY_FILE_READ(iunit1,data(:,:,:,var),nx*ny*nz,kind(data),ierr)
     print*,"min: ",minval(data(:,:,:,var))," - max: ",maxval(data(:,:,:,var))
     !print*,maxloc(data(:,:,:,var))
  end do
  call BINARY_FILE_CLOSE(iunit1,ierr)

  open(iunit2,file=trim(filename2),form="formatted",iostat=ierr,status="REPLACE")
  print*,'Open file : ',filename2, iunit2
  allocate(newdata(nvar+3))

  print*,'Select case:'
  print*,'1 - Along Z, from xmin to xmax'
  print*,'2 - Along Z, from ymin to ymax'
  print*,'3 - Along y, for given x,z'
  print*,'4 - Along x, for given y,z'
  print "(a16,$)", "Select case: "
  read(*,*) sel_case

  select case(sel_case)
  case(1)

     print "(a16,$)", "Minimum x : "
     read(*,*) minx  
     print "(a16,$)", "Maximum x : "
     read(*,*) maxx

     do i = 1, nx
        if (x1(i) .le. minx) then
           minnx = i
        else
           exit
        end if
     end do
     do i = nx, minnx, -1
        if (x1(i) .ge. maxx) then
           maxnx = i
        else
           exit
        end if
     end do
     print*,'nx -> min : ',minnx,'max : ',maxnx

     if (nz .gt. 1) then
        print "(a16,$)", " Select z : "
        read(*,*) selz
        do i = 1, nz
           if (z1(i) .gt. selz) exit
           selnz = i
        end do
     else
        selnz = nz
     end if
     print*,'nz selected :',selnz
  
     print "(a16,$)", " Select Z : "
     read(*,*) value

     !WRITE(frmt, '(1I5)' ) nvar+3
     !frmt = "("//frmt(4:5)//"E24.12E4)"

     do i = minnx, maxnx
        newdata(:) = 0.0_WP
        dum1 = 0.0_WP
        dum2 = 1.0_WP
        do j = 1, ny
           if ((data(i,j,selnz,nvar).ge.dum1).and.(data(i,j,selnz,nvar).le.value)) then
              minny = j
              dum1 = data(i,j,selnz,nvar)
           end if
           if ((data(i,j,selnz,nvar).le.dum2).and.(data(i,j,selnz,nvar).ge.value)) then
              maxny = j
              dum2 = data(i,j,selnz,nvar)
           end if
        end do
        print*,'ny(',i,')  min : ',minny,'max : ',maxny
        print*,'Z(',i,minny,selnz,')  min : ',dum1
        print*,'Z(',i,maxny,selnz,')  max : ',dum2
        if (minny .eq. maxny) then
           DZ = 0.0_WP
           a1 = 1.0_WP
           a2 = 0.0_WP
        else
           DZ = dabs(dum2-dum1)
           a1 = (value - dum1)/DZ
           a2 = (dum2 - value)/DZ
        end if
        !print*,'DZ',DZ
        !print*,'a1',a1
        !print*,'a2',a2
        newdata(1) = x1(i)
        newdata(2) = a1*y1(minny)+a2*y1(maxny)
        newdata(3) = z1(selnz)
        newdata(4:nvar+3) = a1*data(i,minny,selnz,1:nvar)+a2*data(i,maxny,selnz,1:nvar)
        !print*,'newdata',newdata(1:nvar+3)

        write(iunit2,'(50E24.12E4)') newdata(1:nvar+3)
        !write(iunit2,frmt) newdata(1:nvar+3)
     end do

  case(2)

     print "(a16,$)", "Minimum y : "
     read(*,*) miny
     print "(a16,$)", "Maximum y : "
     read(*,*) maxy

     do i = 1, ny
        if (y1(i) .le. miny) then
           minny = i
        else
           exit
        end if
     end do
     do i = ny, minny, -1
        if (y1(i) .ge. maxy) then
           maxny = i
        else
           exit
        end if
     end do
     print*,'ny -> min : ',minny,'max : ',maxny

     if (nz .gt. 1) then
        print "(a16,$)", " Select z : "
        read(*,*) selz
        do i = 1, nz
           if (z1(i) .gt. selz) exit
           selnz = i
        end do
     else
        selnz = nz
     end if
     print*,'nz selected :',selnz

     print "(a16,$)", " Select Z : "
     read(*,*) value

     do j = minny, maxny
        newdata(:) = 0.0_WP
        dum1 = 0.0_WP
        dum2 = 1.0_WP
        do i = 1, nx
           if ((data(i,j,selnz,nvar).ge.dum1).and.(data(i,j,selnz,nvar).le.value)) then
              minnx = i
              dum1 = data(i,j,selnz,nvar)
           end if
           if ((data(i,j,selnz,nvar).le.dum2).and.(data(i,j,selnz,nvar).ge.value)) then
              maxnx = i
              dum2 = data(i,j,selnz,nvar)
           end if
        end do
        print*,'nx(',j,')  min : ',minnx,'max : ',maxny
        print*,'Z(',minnx,j,selnz,')  min : ',dum1
        print*,'Z(',maxnx,j,selnz,')  max : ',dum2
        if (minnx .eq. maxnx) then
           DZ = 0.0_WP
           a1 = 1.0_WP
           a2 = 0.0_WP
        else
           DZ = dabs(dum2-dum1)
           a1 = (value - dum1)/DZ
           a2 = (dum2 - value)/DZ
        end if

        newdata(1) = a1*x1(minnx)+a2*x1(maxnx)
        newdata(2) = y1(j)
        newdata(3) = z1(selnz)
        newdata(4:nvar+3) = a1*data(minnx,j,selnz,1:nvar)+a2*data(maxnx,j,selnz,1:nvar)

        write(iunit2,'(50E24.12E4)') newdata(1:nvar+3)
     end do

  case(3)
     if (nx .gt. 1) then
        print "(a16,$)", " Select x : "
        read(*,*) selx
        do i = 1, nx
           if (x1(i) .gt. selx) exit
           selnx = i
        end do
     else
        selnx = nx
     end if
     print*,'nx selected :',selnx

     if (nz .gt. 1) then
        print "(a16,$)", " Select z : "
        read(*,*) selz
        do i = 1, nz
           if (z1(i) .gt. selz) exit
           selnz = i
        end do
     else
        selnz = nz
     end if
     print*,'nz selected :',selnz

     do j = 1, ny
        newdata(1) = x1(selnx)
        newdata(2) = y1(j)
        newdata(3) = z1(selnz)
        newdata(4:nvar+3) = data(selnx,j,selnz,1:nvar)

        write(iunit2,'(50E24.12E4)') newdata(1:nvar+3)
     end do

  case(4)
     if (nx .gt. 1) then
        print*, 'y min = ',y1(1),'; y max = ',y1(ny)
        print "(a16,$)", " Select y : "
        read(*,*) sely
        do i = 1, ny
           if (y1(i) .gt. sely) exit
           selny = i
        end do
     else
        selny = ny
     end if
     print*,'ny selected :',selny

     if (nz .gt. 1) then
        print*, 'z min = ',z1(1),'; z max = ',z1(nz)
        print "(a16,$)", " Select z : "
        read(*,*) selz
        do i = 1, nz
           if (z1(i) .gt. selz) exit
           selnz = i
        end do
     else
        selnz = nz
     end if
     print*,'nz selected :',selnz

     do i = 1, nx
        newdata(1) = x1(i)
        newdata(2) = y1(selny)
        newdata(3) = z1(selnz)
        newdata(4:nvar+3) = data(i,selny,selnz,1:nvar)

        write(iunit2,'(50E24.12E4)') newdata(1:nvar+3)
     end do

  end select

  close(iclose(iunit2))
  
end program extData
