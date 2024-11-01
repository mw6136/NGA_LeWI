! -------------------------------------------------------------------------- !
!                         COMPAREPLANEDATA.F90                               !
!     Purpose     : Compares old and new style binary plane files            !
!     Date        : July 20, 2016                                            !
!     Modified by : Jonathan F. MacArt                                       !
! -------------------------------------------------------------------------- !
program comparePlaneData
  use precision
  use string
  use fileio
!  use cli_reader
  implicit none

  character(len=str_long) :: filename1, filename2
  integer :: iunit1, iunit2, ierr
  integer :: ny1, nz1, nplanes1, nvars1, ntime1, planetype1
  integer :: ny2, nz2, nplanes2, nvars2, ntime2, planetype2
  integer :: it,n,i,j,jj,k,kk
  real(WP) :: dt1, time1
  real(WP) :: dt2, time2, tmp, tmp1, tmp2, max1, min1, max2, min2
  real(WP), dimension(:),     pointer :: x1, x2
  real(WP), dimension(:),     pointer :: y1, y2
  real(WP), dimension(:),     pointer :: z1, z2
  real(WP), dimension(:,:,:), pointer :: data1, data2
  character(len=str_short), dimension(:), pointer :: names1, names2
  integer, dimension(:), pointer :: jmatch

  print *
!!$  print "(a24,$)", " plane file of type 1 : "
!!$  read "(a)", filename1
!!$  print "(a24,$)", " plane file of type 2 : "
!!$  read "(a)", filename2
  !filename1 = "./plane_data_jet5k_PM_20pN2dil.43"
  !filename2 = "./plane_data_jet5k_PM_20pN2dil.44"
  filename1 = "./jet5k_NP_20pN2dil/plane_data_jet5k_NP_20pN2dil.19"
  filename2 = "./test2/plane_data_jet5k_NP_20pN2dil.13"
  print *, trim(filename1)
  print *, trim(filename2)

  ! Open the binary files
  call BINARY_FILE_OPEN(iunit1, trim(filename1), "r", ierr)
  call BINARY_FILE_OPEN(iunit2, trim(filename2), "r", ierr)

  ! Read data sizes
  ! File 1
  call BINARY_FILE_READ(iunit1, ntime1, 1,   kind(ntime1), ierr)
  call BINARY_FILE_READ(iunit1, ny1, 1,      kind(ny1), ierr)
  call BINARY_FILE_READ(iunit1, nz1, 1,      kind(nz1), ierr)
  call BINARY_FILE_READ(iunit1, nplanes1, 1, kind(nplanes1), ierr)
  call BINARY_FILE_READ(iunit1, planetype1, 1,kind(planetype1), ierr)
  call BINARY_FILE_READ(iunit1, nvars1, 1,   kind(nvars1), ierr)
  ! File 2
  call BINARY_FILE_READ(iunit2, ntime2, 1,   kind(ntime2), ierr)
  call BINARY_FILE_READ(iunit2, ny2, 1,      kind(ny2), ierr)
  call BINARY_FILE_READ(iunit2, nz2, 1,      kind(nz2), ierr)
  call BINARY_FILE_READ(iunit2, nplanes2, 1, kind(nplanes2), ierr)
  call BINARY_FILE_READ(iunit2, planetype2, 1,kind(planetype2), ierr)
  call BINARY_FILE_READ(iunit2, nvars2, 1,   kind(nvars2), ierr)

  print *, 'ntime :   ', ntime1, ntime2
  print *, 'ny :      ', ny1, ny2
  print *, 'nz :      ', nz1, nz2
  print *, 'nplanes : ', nplanes1, nplanes2
  print *, 'nvars :   ', nvars1, nvars2
  if (planetype2.eq.1) then
     print *, ' file 2 plane is in xy'
  elseif (planetype2.eq.2) then
     print *, ' file 2 plane is in xz'
  elseif (planetype2.eq.3) then
     print *, ' file 2 plane is in yz'
  else
     print *, ' file 2 plane type not recognized'
  end if

  ! Read y-locations
  allocate(y1(ny1))
  do j=1,ny1
     call BINARY_FILE_READ(iunit1, y1(j), 1, kind(y1), ierr)
  end do
  allocate(y2(ny2))
  do j=1,ny2
     call BINARY_FILE_READ(iunit2, y2(j), 1, kind(y2), ierr)
  end do
  print *
  !print *, 'plane 1 y = ', y1
  !print *, 'plane 2 y = ', y2
  tmp1 = sum(y1)/ny1
  tmp2 = sum(y2)/ny2
  min1 = minval(y1); max1 = maxval(y1)
  min2 = minval(y2); max2 = maxval(y2)
  tmp = 0.0_WP
  do j=1,minval([ny1,ny2])
     tmp = tmp + y1(j) - y2(j)
  end do
  print "(a11,a11,a11,a11,a11,a11,a11)", 'y-min-1', 'y-min-2', 'y-max-1', 'y-max-2', 'y-avg-1', 'y-avg-2', 'delta'
  print "(ES11.3,ES11.3,ES11.3,ES11.3,ES11.3,ES11.3,ES11.3)", min1, min2, max1, max2, tmp1, tmp2, abs(tmp)
  ! Verify y-locations
  allocate(jmatch(1:ny1))
  jmatch = 0
  do j=1,ny1
     do jj=1,ny2
        if (y2(jj).eq.y1(j)) jmatch(j) = jj
     end do
  end do

  ! Read z
  allocate(z1(nz1)) !! Changed the format of the old style file (in branch test_spectrum) on 8/18/16 to compare z-locations!
  do j=1,nz1
     call BINARY_FILE_READ(iunit1, z1(j), 1, kind(z1), ierr)
  end do
  !print *, 'plane 1 z = ', z1
  allocate(z2(nz2))
  do j=1,nz2
     call BINARY_FILE_READ(iunit2, z2(j), 1, kind(z2), ierr)
  end do
  !print *, 'plane 2 z = ', z2
  tmp1 = sum(z1)/nz1
  tmp2 = sum(z2)/nz2
  min1 = minval(z1); max1 = maxval(z1)
  min2 = minval(z2); max2 = maxval(z2)
  tmp = 0.0_WP
  do j=1,minval([nz1,nz2])
     tmp = tmp + z1(j) - z2(j)
  end do
  print *, ' '
  print "(a11,a11,a11,a11,a11,a11,a11)", 'z-min-1', 'z-min-2', 'z-max-1', 'z-max-2', 'z-avg-1', 'z-avg-2', 'delta'
  print "(ES11.3,ES11.3,ES11.3,ES11.3,ES11.3,ES11.3,ES11.3)", min1, min2, max1, max2, tmp1, tmp2, abs(tmp)

  ! Read x-locations of planes
  allocate(x1(nplanes1))
  do n=1,nplanes1
     call BINARY_FILE_READ(iunit1, x1(n), 1, kind(x1), ierr)
  end do
  allocate(x2(nplanes2))
  do n=1,nplanes2
     call BINARY_FILE_READ(iunit2, x2(n), 1, kind(x2), ierr)
  end do
  print *
  print *, 'planes in file 1 at x = ', x1
  print *, 'planes in file 2 at x = ', x2

  ! Read variable names
  allocate(names1(nvars1))
  do n=1,nvars1
     call BINARY_FILE_READ(iunit1, names1(n), str_short, kind(names1), ierr)
  end do
  allocate(names2(nvars2))
  do n=1,nvars2
     call BINARY_FILE_READ(iunit2, names2(n), str_short, kind(names2), ierr)
  end do
  print *, ' '
  print *, 'Variables in file 1 : ', names1
  print *, 'Variables in file 2 : ', names2

  ! Read the data
  allocate(data1(nplanes1,ny1,nz1))
  allocate(data2(nplanes2,ny2,nz2))
  do it=1,minval([ntime1,ntime2])
     ! Read time info
     call BINARY_FILE_READ(iunit1, dt1,   1, kind(dt1),   ierr)
     call BINARY_FILE_READ(iunit1, time1, 1, kind(time1), ierr)
     call BINARY_FILE_READ(iunit2, dt2,   1, kind(dt2),   ierr)
     call BINARY_FILE_READ(iunit2, time2, 1, kind(time2), ierr)
     print *, '****************************'
     print *, 'ntime : ', it
     print *, 'dt    : ', dt1, dt2
     print *, 'time  : ', time1, time2
     print "(a12,a11,a11,a11,a11,a11,a11,a11)", '            ', 'min-1', 'min-2', 'max-1', 'max-2', 'avg-1', 'avg-2', 'delta'
     do n=1,nvars1
        ! Read from file 1
        do k=1,nz1
           do j=1,ny1
              do i=1,nplanes1
                 call BINARY_FILE_READ(iunit1, data1(i,j,k), 1, kind(data1), ierr)
              end do
           end do
        end do

        ! Read from file 2
        do k=1,nz2
           do j=1,ny2
              do i=1,nplanes2
                 call BINARY_FILE_READ(iunit2, data2(i,j,k), 1, kind(data2), ierr)
              end do
           end do
        end do

        ! Compare the two
!!$        print *, ' var = ', names1(n)
!!$        do i=1,minval([nplanes1,nplanes2])
!!$           do j=1,ny1
!!$              do k=1,nz1
!!$                 do jj=1,ny2
!!$                    do kk=1,nz2
!!$                       if (data1(i,j,k).eq.data2(i,jj,kk)) then
!!$                          print *, 'match at j,jj,k,kk = ', j, jj, k, kk
!!$                       end if
!!$                    end do
!!$                 end do
!!$              end do
!!$           end do
!!$        end do

!!$        print *, ' var = ', names1(n)
!!$        do i=1,minval([nplanes1,nplanes2])
!!$           do j=1,ny1
!!$              do k=1,nz1
!!$                 print '(i4,i4,i4,ES21.7,ES21.7)', i, j, k, data1(i,j,k), data2(i,jmatch(j),k)
!!$              end do
!!$           end do
!!$        end do

        ! Compute normalized delta
        tmp  = 0.0_WP
        do i=1,1!minval([nplanes1,nplanes2])
           do j=1,ny1
              do k=1,nz1
                 tmp  = tmp + data1(i,j,k) - data2(i,j,k)
!!$                 if (abs(data1(i,j,k)).gt.100000_WP .and. names1(n).ne.'HR') then
!!$                    print "(i2,a12,ES13.5,ES13.5,ES13.5,ES13.5)", 1, names1(n), x1(i), y1(j), z1(k), data1(i,j,k)
!!$                 end if
!!$                 if (abs(data2(i,j,k)).gt.100000_WP .and. names2(n).ne.'HR') then
!!$                    print "(i2,a12,ES13.5,ES13.5,ES13.5,ES13.5)", 2, names2(n), x2(i), y2(j), z2(k), data2(i,j,k)
!!$                 end if
              end do
           end do
        end do
        tmp1 = sum(data1(1:1,1:ny1,1:nz1))/(1*ny1*nz1)
        min1 = minval(data1(1:1,1:ny1,1:nz1))
        max1 = maxval(data1(1:1,1:ny1,1:nz1))
        tmp2 = sum(data2(1:1,1:ny1,1:nz1))/(1*ny1*nz1)
        min2 = minval(data2(1:1,1:ny1,1:nz1))
        max2 = maxval(data2(1:1,1:ny1,1:nz1))
        print "(a12,ES11.3,ES11.3,ES11.3,ES11.3,ES11.3,ES11.3,ES11.3)", names1(n), min1, min2, max1, max2, tmp1, tmp2, abs(tmp)

     end do


  end do

  ! Close the file
  call BINARY_FILE_CLOSE(iunit1, ierr)
  call BINARY_FILE_CLOSE(iunit2, ierr)

end program comparePlaneData
