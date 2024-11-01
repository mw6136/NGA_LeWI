module dump_scalar_pdf
  use parallel
  use precision
  use geometry
  use partition
  use masks
  use fileio
  use unstruct
  implicit none
  
  integer :: njoint,nsingle,nbins
  integer, dimension(:), pointer :: Z0_isc,Z1_isc,Z2_isc

end module dump_scalar_pdf


! ================================================= !
! Dump 3D binary ensight gold data - Initialization !
! ================================================= !
subroutine dump_scalar_pdf_init
  use dump_scalar_pdf
  use parser
  use parallel
  use time_info
  use data
  implicit none
  integer :: ii,isc
  logical, dimension(:), pointer :: found0,found1,found2
  character(len=str_medium), dimension(:), pointer :: vars_in_single,vars_in_joint

  ! Create & Start the timer
  call timing_create('scalarpdf')
  call timing_start ('scalarpdf')

  ! Create directory
  if (irank.eq.iroot) call CREATE_FOLDER("scalarPDF")

  ! Find the variable indices for Joint case
  call parser_getsize('Joint scalar PDF variables',njoint)
  if (mod(njoint,2) .ne. 0) then
     call die('dump_scalar_pdf_init: must provide scalar pairs for joint pdfs')
  else
     allocate(vars_in_joint(njoint))
     call parser_read('Joint scalar PDF variables',vars_in_joint)
  end if
  njoint = njoint/2
  allocate(Z1_isc(njoint))
  allocate(Z2_isc(njoint))
  allocate(found1(njoint))
  allocate(found2(njoint))
  found1 = .false.
  found2 = .false.
  do ii = 1,njoint
     do isc=1,nscalar
        if (trim(SC_name(isc)).eq.vars_in_joint(ii*2-1)) then
           Z1_isc(ii) = isc
           found1(ii) = .true.
        end if
        if (trim(SC_name(isc)).eq.vars_in_joint(ii*2)  ) then
           Z2_isc(ii) = isc
           found2(ii) = .true.
        end if
     end do
  end do
  do ii=1,njoint
     if ((found1(ii) .eq. .false.) .or. (found2(ii) .eq. .false.)) then
        call die('dump_scalar_pdf_init: invalid scalar requested')
     end if
  end do
  deallocate(found1)
  deallocate(found2)
  
  ! Find the variable indices for single variable case
  call parser_getsize('Single scalar PDF variables',nsingle)
  allocate(vars_in_single(nsingle))
  call parser_read('Single scalar PDF variables',vars_in_single)
  allocate(Z0_isc(nsingle))
  allocate(found0(nsingle))
  found0=.false.
  do ii = 1,nsingle
     do isc=1,nscalar
        if (trim(SC_name(isc)).eq.vars_in_single(ii)) then
           Z0_isc(ii) = isc
           found0(ii) = .true.
        end if
     end do
  end do
  do ii=1,nsingle
     if (found0(ii) .eq. .false.) then
        call die('dump_scalar_pdf_init: invalid scalar requested')
     end if
  end do
  deallocate(found0)

  ! get number of bins to use
  call parser_read('Number of bins',nbins)

  ! Stop the timer
  call timing_stop ('scalarpdf')
  
  ! Save the first field
  call dump_scalar_pdf_dump

  return
end subroutine dump_scalar_pdf_init


! =========================================== !
! Dump count of cells that fall into each bin !
! =========================================== !
subroutine dump_scalar_pdf_dump
  use dump_scalar_pdf
  use data
  use partition
  use string
  use parallel
  use time_info
  implicit none

  integer :: i,j,k,isc,ibin1,ibin2
  real(WP) :: zvalue
  integer, dimension(:,:,:), pointer :: bins2D,bins2Dout
  integer, dimension(:,:), pointer :: bins1D,bins1Dout
  
  ! variables for printing to file
  character(len=str_medium) :: filename2,filename0,buffer
  integer :: iunit, ierr
  real(WP) :: rbin1, rbin2

  ! Start the timer
  call timing_start('scalarpdf')

  allocate(bins2D(nbins,nbins,njoint))
  allocate(bins1D(nbins,nsingle))
  allocate(bins2Dout(nbins,nbins,njoint))
  allocate(bins1Dout(nbins,nsingle))
  bins2D = 0
  bins1D = 0
  bins2Dout = 0
  bins1Dout = 0

  ! Find the appropriate bins, put values <0 and >1 in the first and last bins
  do k= kmin_,kmax_
     do j= jmin_,jmax_
        do i= imin_,imax_
           do isc=1,njoint
              zvalue = min(max(SC(i,j,k,Z1_isc(isc)), 1.0E-6_WP),(1.0_WP-1.0E-6_WP))
              ibin1 = int( zvalue * real(nbins)) + 1
              zvalue = min(max(SC(i,j,k,Z2_isc(isc)), 1.0E-6_WP),(1.0_WP-1.0E-6_WP))
              ibin2 = int( zvalue * real(nbins)) + 1
              bins2D(ibin1,ibin2,isc) = bins2D(ibin1,ibin2,isc) + 1
           end do
           do isc=1,nsingle
              zvalue = min(max(SC(i,j,k,Z0_isc(isc)), 1.0E-6_WP),(1.0_WP-1.0E-6_WP))
              ibin1 = int( zvalue * real(nbins)) + 1
              bins1D(ibin1,isc) = bins1D(ibin1,isc) + 1
           end do
        end do
     end do
  end do

  ! Sum over all processors and print
  call parallel_sum(bins2D,bins2Dout)
  call parallel_sum(bins1D,bins1Dout)
  if (irank.eq.iroot) then
     write(buffer,'(e13.5e2)'), time
     filename2 = 'scalarPDF/joint-'    // trim(adjustl(buffer))
     filename0 = 'scalarPDF/marginal-' // trim(adjustl(buffer))
     
     iunit = iopen()
     open(iunit, file=filename2,form="formatted",iostat=ierr,status="REPLACE")
     write(iunit,'(1000(a18))') "bins-v1","bins-v2", &      
          (trim(SC_name(Z1_isc(isc))) // '-' // trim(SC_name(Z2_isc(isc))), isc=1,njoint)
     do ibin1=1,nbins
        do ibin2=1,nbins
           rbin1 = (real(ibin1,WP) - 0.5_WP)/real(nbins,WP)
           rbin2 = (real(ibin2,WP) - 0.5_WP)/real(nbins,WP)
           write(iunit,'(F18.8,F18.8,1000(I18))') rbin1, rbin2, & 
                (bins2Dout(ibin1,ibin2,isc), isc=1,njoint)
        end do
     end do
     close(iclose(iunit))

     iunit = iopen()
     open(iunit, file=filename0,form="formatted",iostat=ierr,status="REPLACE")
     write(iunit,'(1000(a18))') "bins", (trim(SC_name(Z0_isc(isc))), isc=1,nsingle)
     do ibin1=1,nbins
        rbin1 = (real(ibin1,WP) - 0.5_WP)/real(nbins,WP)
        write(iunit,'(F18.8,1000(I18))') rbin1, (bins1Dout(ibin1,isc), isc=1,nsingle)  
     end do
     close(iclose(iunit))
  end if

  deallocate(bins2D)
  deallocate(bins1D)
  deallocate(bins2Dout)
  deallocate(bins1Dout)

  ! Stop the timer
  call timing_stop('scalarpdf')

  return
end subroutine dump_scalar_pdf_dump
