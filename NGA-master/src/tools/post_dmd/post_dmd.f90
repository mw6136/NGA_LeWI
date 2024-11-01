! ============================================================ !
!                         dmd.f90                              !
!                                                              !
! Author: Temistocle Grenga                                    !
! Date:   March 9, 2017                                        !
! ============================================================ !

program post_dmd
  use parallel
  use string
  use precision
  use geometry
  use partition
  use masks
  use fileio
  use data
  use masks

  implicit none
  character(len=str_medium) :: folder, filename1, filename_base, fname
  integer :: i, j, k, l, m, pr, o
  integer :: iunit0, iunit1, iunit2, iunit3, ierr, ifile
  integer :: n_snap, start_file, step_file, nvars, nvars1, ntime1, nx1, ny1, nz1
  real(WP), dimension(:,:), pointer :: Lambda
  real(WP), dimension(:), pointer :: x1, y1, z1
  character(len=str_short), dimension(:), pointer :: names
  character(len=80) :: cdum
  real(SP) :: ddum, IntV, IntV2, VolD, a, ampl, ampl2
  integer :: idum
  real(SP), dimension(:,:,:), pointer :: vv, vm, vs, res
  real(SP), dimension(:,:), pointer :: array_res, array2_res, array_ampl, array_ampl2
  integer, dimension(:), pointer :: var_list, array_idx

  !Initializa MPI enironment
  call parallel_init()

  start_file = 15286  !6465 !8830 !15686 !21500
  step_file = 1
  n_snap = 400
  nvars = 1

  allocate(var_list(nvars))
  var_list = (/ 1 /)
  !var_list = (/ 1, 2, 3, 4, 10, 12 /)

  write(folder, '(''U_'',I4.4,''_'', I2.2,''/'')') n_snap, step_file
  !n_snap = 60
  filename1 = 'Lambda'
  filename1 = trim(folder)//trim(filename1)
  open(iunit0,file=trim(filename1),form="formatted",iostat=ierr)
  allocate(Lambda(n_snap,2))
  do i = 1,n_snap
     read(iunit0,'(2E24.12E4)') Lambda(i,1), Lambda(i,2)
  end do
  close(iunit0)

  print*, 'End file 1'

  filename_base = 'vol_data.1'
  l = start_file +step_file*n_snap
  write(fname, '(''_'', I8.8)') l
  fname = trim(filename_base)//trim(fname)
  filename1 = "volume_data/"//trim(fname)
  print*, filename1
  call BINARY_FILE_OPEN(iunit1,trim(filename1),"r",ierr)
  print*, filename1
  call BINARY_FILE_READ(iunit1, ntime1, 1, kind(ntime1), ierr)
  call BINARY_FILE_READ(iunit1, nx1, 1,    kind(nx1), ierr)
  call BINARY_FILE_READ(iunit1, ny1, 1,    kind(ny1), ierr)
  call BINARY_FILE_READ(iunit1, nz1, 1,    kind(nz1), ierr)
  call BINARY_FILE_READ(iunit1, nvars1, 1, kind(nvars1), ierr)
  print*, nx1, ny1, nz1
  allocate(x1(0:nx1+1), y1(0:ny1+1), z1(0:nz1+1))
  do i=1,nx1
     call BINARY_FILE_READ(iunit1, x1(i), 1, kind(x1), ierr)
  end do
  do i=1,ny1
     call BINARY_FILE_READ(iunit1, y1(i), 1, kind(y1), ierr)
  end do
  do i=1,nz1
     call BINARY_FILE_READ(iunit1, z1(i), 1, kind(z1), ierr)
  end do
  x1(0) = x1(1) - (x1(2)-x1(1)); x1(nx1+1) = x1(nx1) + x1(nx1)-x1(nx1-1)
  y1(0) = y1(1) - (y1(2)-y1(1)); y1(ny1+1) = y1(ny1) + y1(ny1)-y1(ny1-1)
  z1(0) = z1(1) - (z1(2)-z1(1)); z1(nz1+1) = z1(nz1) + z1(nz1)-z1(nz1-1)

  allocate(names(nvars1))
  do i = 1,nvars1
     call BINARY_FILE_READ(iunit1, names(i), str_short, kind(names), ierr)
  end do
  call BINARY_FILE_CLOSE(iunit1,ierr)
  print*, 'End file 2'

  print*, 'Min x ', minval(x1)
  print*, 'Max x ', maxval(x1)
  print*, 'Min y ', minval(y1)
  print*, 'Max y ', maxval(y1)
  print*, 'Min z ', minval(z1)
  print*, 'Max z ', maxval(z1)
  VolD = REAL((maxval(x1)-minval(x1))*(maxval(y1)-minval(y1))*(maxval(z1)-minval(z1)),SP)
  print*, 'Domain volume ', VolD

!  filename1 = 'Residual'
!  filename1 = trim(folder)//trim(filename1)
!  open(iunit2,file=trim(filename1),form="formatted",iostat=ierr)

!  fname = 'Newdmd.case'
!  filename1 = trim(folder)//trim(fname)
!  open(iunit3,file=trim(filename1),form="formatted",iostat=ierr,status="REPLACE")
!  print*,filename1
!  write(iunit3,'(a)') 'FORMAT'
!  write(iunit3,'(a)') 'type: ensight gold'
!  write(iunit3,'(a)') 'GEOMETRY'
!  write(iunit3,'(a)') 'model: 1 1 geometry'
!  write(iunit3,'(a)') 'VARIABLE'

  allocate (array_res(nvars,n_snap), array2_res(nvars,n_snap), array_ampl(nvars,n_snap), array_ampl2(nvars,n_snap),array_idx(n_snap))
  array_idx = -1
  do m = 1,nvars
!     write(iunit3,'(10a)') 'scalar per node: ', trim(names(var_list(m))),' ',trim(names(var_list(m)))

     allocate(vv(nx1,ny1,nz1))
     fname = trim(folder)//trim(adjustl(names(var_list(m))))
     !fname = trim(folder)//"rr_"//trim(adjustl(names(var_list(m))))
     print*,fname
     call BINARY_FILE_OPEN(ifile,trim(fname),"r",ierr)
     call BINARY_FILE_READ(ifile, cdum, 80, kind(cdum), ierr)
     call BINARY_FILE_READ(ifile, cdum, 80, kind(cdum), ierr)
     call BINARY_FILE_READ(ifile, cdum, 1, kind(idum), ierr)
     call BINARY_FILE_READ(ifile, cdum, 80, kind(cdum), ierr)
     call BINARY_FILE_READ(ifile, vv, nx1*ny1*nz1, kind(ddum), ierr)
     call BINARY_FILE_CLOSE(ifile,ierr)

     allocate(vs(nx1,ny1,nz1),vm(nx1,ny1,nz1),res(nx1,ny1,nz1))
     vs = 0.0_SP
     l = 1; pr = 0
     do
        pr = pr+1
        o = l
        vm = 0.0_SP
        fname = trim(folder)//trim(adjustl(names(var_list(m))))// "_"
        !fname = trim(folder)//"rr_"//trim(adjustl(names(var_list(m))))! // "_"
        write(fname(len_trim(fname)+1:len_trim(fname)+3),'(i3.3)') l
        print*,fname
        call BINARY_FILE_OPEN(iunit1,trim(fname),"r",ierr)
        call BINARY_FILE_READ(iunit1, cdum, 80, kind(cdum), ierr)
        call BINARY_FILE_READ(iunit1, cdum, 80, kind(cdum), ierr)
        call BINARY_FILE_READ(iunit1, cdum, 1, kind(idum), ierr)
        call BINARY_FILE_READ(iunit1, cdum, 80, kind(cdum), ierr)
        call BINARY_FILE_READ(iunit1, vm, nx1*ny1*nz1, kind(vm),ierr)
        call BINARY_FILE_CLOSE(iunit1,ierr)
        !print*, 'vm',vm(1,1,1), vm(nx1,ny1,nz1)
        print*, 'close ', fname
        if (DABS(Lambda(l,2)) .gt. 1.0D-13) then
           l = l+1
           a = 2.0_SP
        else
           a = 1.0_SP
        end if
        vs(1:nx1,1:ny1,1:nz1) = vs(1:nx1,1:ny1,1:nz1) + a*vm(1:nx1,1:ny1,1:nz1)
        !print*, 'vs',vs(1,1,1), vs(nx1,ny1,nz1)

        if ((MOD(pr,5) == 0) .or. (l == n_snap)) then
!           cdum = trim(adjustl(names(var_list(m)))) // "_1_"
!           write(cdum(len_trim(cdum)+1:len_trim(cdum)+3),'(i3.3)') l
!           write(iunit3,'(10a)') 'scalar per node: ', trim(cdum),' ',trim(cdum)

           fname = trim(folder)//trim(adjustl(names(var_list(m)))) // "_1_"
           write(fname(len_trim(fname)+1:len_trim(fname)+3),'(i3.3)') l
           call BINARY_FILE_OPEN(ifile,trim(fname),"w",ierr)
           cdum = trim(adjustl(names(var_list(m)))) // "_1_"
           write(cdum(len_trim(cdum)+1:len_trim(cdum)+3),'(i3.3)') l
           call BINARY_FILE_WRITE(ifile,cdum,80,kind(cdum),ierr)
           cdum = 'part'
           call BINARY_FILE_WRITE(ifile,cdum,80,kind(cdum),ierr)
           idum = 1
           call BINARY_FILE_WRITE(ifile,idum,1,kind(idum),ierr)
           cdum = 'block'
           call BINARY_FILE_WRITE(ifile,cdum,80,kind(cdum),ierr)
           call BINARY_FILE_WRITE(ifile,vs,nx1*ny1*nz1,kind(vs),ierr)
           call BINARY_FILE_CLOSE(ifile,ierr)
        end if

        res = 0.0_SP
        res(1:nx1,1:ny1,1:nz1) = vs(1:nx1,1:ny1,1:nz1)-vv(1:nx1,1:ny1,1:nz1)
        !print*, 'res',res(1,1,1), res(nx1,ny1,nz1)
        IntV = 0.0_SP; IntV2 = 0.0_SP
        ddum = 0.0_SP
        ampl = 0.0_SP; ampl2 = 0.0_SP
        do k = 1,nz1
           do j = 1,ny1
              do i = 1,nx1
                   IntV = IntV + ABS(res(i,j,k)/vv(i,j,k))* &
                         REAL(0.5_WP*(x1(i+1)-x1(i-1))*0.5_WP*(y1(j+1)-y1(j-1))*0.5_WP*(z1(k+1)-z1(k-1))/VolD,SP)
                   IntV2 = IntV2 + ((res(i,j,k)/vv(i,j,k))*(res(i,j,k)/vv(i,j,k)))* &
                          REAL(0.5_WP*(x1(i+1)-x1(i-1))*0.5_WP*(y1(j+1)-y1(j-1))*0.5_WP*(z1(k+1)-z1(k-1))/VolD,SP)
                   !if ( ddum .le. ABS(res(i,j,k)) ) ddum = ABS(res(i,j,k))
                   ampl = ampl + ABS(vm(i,j,k))* &
                         REAL(0.5_WP*(x1(i+1)-x1(i-1))*0.5_WP*(y1(j+1)-y1(j-1))*0.5_WP*(z1(k+1)-z1(k-1))/VolD,SP)
                   ampl2 = ampl2 + (vm(i,j,k)**2.0_WP)* &
                          REAL(0.5_WP*(x1(i+1)-x1(i-1))*0.5_WP*(y1(j+1)-y1(j-1))*0.5_WP*(z1(k+1)-z1(k-1))/VolD,SP)
              end do
           end do
        end do
        array_res(m,l) = IntV
        array2_res(m,l) = IntV2
        array_ampl(m,l) = ampl
        array_ampl2(m,l) = ampl2
        array_idx(l) = o
        !write(iunit2,'(2I5,E24.12E4)') m, l, IntV

        l = l+1
        if( l > n_snap ) exit
     end do
        
     deallocate(vv, vm, vs, res)
  end do

!  write(iunit3,'(a)') 'TIME'
!  write(iunit3,'(a)') 'time set: 1'
!  write(iunit3,'(a)') 'number of steps: 1'
!  write(iunit3,'(a)') 'filename start number: 1'
!  write(iunit3,'(a)') 'filename increment: 1'
!  ddum = 2.0D-6
!  write(iunit3,'(a,ES12.5)') 'time values:',ddum
!  close(iunit3)

  print*, 'Complete evaluation'
  !print*, array_res(1:nvars,1)

  filename1 = 'T_Residual'
  filename1 = trim(folder)//trim(filename1)
  print*, filename1
  open(iunit2,file=trim(filename1),form="formatted",iostat=ierr)
  print*, filename1
  do i = 1, n_snap
     if (array_idx(i) > 0) write(iunit2,'(I5,6E24.12E4)') i, array_res(1:nvars,i)
  end do
  close(iunit2)

  filename1 = 'T_Residual2'
  filename1 = trim(folder)//trim(filename1)
  open(iunit2,file=trim(filename1),form="formatted",iostat=ierr)
  do i = 1, n_snap
     if (array_idx(i) > 0) write(iunit2,'(I5,6E24.12E4)') i, array2_res(1:nvars,i)
  end do
  close(iunit2)

  filename1 = 'T_Amplitude'
  filename1 = trim(folder)//trim(filename1)
  open(iunit2,file=trim(filename1),form="formatted",iostat=ierr)
  do i = 1, n_snap
     if (array_idx(i) > 0) write(iunit2,'(I5,6E24.12E4)') array_idx(i), array_ampl(1:nvars,i)
  end do
  close(iunit2)

  filename1 = 'T_Amplitude2'
  filename1 = trim(folder)//trim(filename1)
  open(iunit2,file=trim(filename1),form="formatted",iostat=ierr)
  do i = 1, n_snap
     if (array_idx(i) > 0) write(iunit2,'(I5,6E24.12E4)') array_idx(i), array_ampl2(1:nvars,i)
  end do
  close(iunit2)

  deallocate(Lambda)
  deallocate(x1, y1, z1, names)


  ! Finalize the parallel environment
  call parallel_final()
  print *, irank, 'Completed post_dmd'

end program post_dmd
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
