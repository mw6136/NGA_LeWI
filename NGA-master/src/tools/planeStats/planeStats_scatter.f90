module planeStats_scatter
  use planeStats

  integer :: Nvar_mpdf, Nvar_jpdf
  real(WP), dimension(:,:),   pointer :: mpdf
  real(WP), dimension(:,:,:), pointer :: jpdf

contains
  
  ! ========================================================== !
  ! Write the header for scatter data                          !
  ! ========================================================== !
  subroutine write_header(name)
    implicit none

    integer :: iplane, j, k, iunit
    character(len=str_medium) :: of_name
    character(len=*), intent(in) :: name

    ! Initialize the file and write headers
    do iplane=pnmin_,pnmax_
       if (iplane<10) then
          write(of_name,'(A1,I1)') '0', iplane
       else
          write(of_name,'(I2)') iplane
       end if
       of_name = trim(output_name)//'_'//trim(name)//'_ip'//trim(of_name)
       open(unit=iunit, file=trim(of_name), form='formatted', action='write')
       write(iunit,'(a15)',advance='no') 'time'
       do j=1,ny
          do k=1,nz
             write(iunit,'(a2,I5,a1,I5,a2)',advance='no') '  ',j,',',k,'  '
          end do
       end do
       close(iunit)
    end do

    return
  end subroutine write_header

  ! ========================================================== !
  ! Write a line of scatter data                               !
  ! ========================================================== !
  subroutine write_step(iplane,name,tmpxy)
    implicit none

    integer, intent(in) :: iplane
    integer :: j, k, iunit
    real(WP), dimension(1:nz,1:ny), intent(in) :: tmpxy
    character(len=str_medium) :: of_name
    character(len=*), intent(in) :: name

    if (iplane<10) then
       write(of_name,'(A1,I1)') '0', iplane
    else
       write(of_name,'(I2)') iplane
    end if
    of_name = trim(output_name)//'_'//trim(name)//'_ip'//trim(of_name)
    open(unit=iunit, file=trim(of_name), form='formatted', status='old', position='append', action='write')
    write(iunit,'(ES15.5)',advance='no') time(ntime_curr)
    do j=1,ny
       do k=1,nz
          write(iunit,'(ES15.5)',advance='no') tmpxy(k,j)
       end do
    end do
    close(iunit)

    return
  end subroutine write_step

end module planeStats_scatter


! ========================================================== !
! Initialize the module                                      !
! ========================================================== !
subroutine planeStats_scatter_init
  use planeStats_scatter
  implicit none

  integer :: isc

  if (output_scatter) then
     call write_header('RHO')
     call write_header('U1')
     call write_header('U2')
     call write_header('U3')
     do isc=1,nscalar
        call write_header(SC_name(isc))
     end do
  end if

  if (output_pdfs) then
     ! Do the pdfs need to be computed in parallel? Probably not. plane min/max in init
     ! ALSO:  need to read all data, THEN compute pdfs -- as in matlab. Need correct min/max
  end if

  return
end subroutine planeStats_scatter_init



! ========================================================== !
! Write the instantaneous statistics to text files           !
!   and/or update pdfs                                       !
! ========================================================== !
subroutine planeStats_scatter_write
  use planeStats_scatter
  implicit none

  integer :: iplane, j, k, isc, iunit
!!$  real(WP), dimension(1:nz,1:ny) :: tmpxy
!!$  character(len=str_medium) :: of_name

  if (output_scatter) then
     ! Write something
     do iplane=pnmin_,pnmax_
        call write_step(iplane,'RHO',RHO(:,:,iplane))
     end do
     do iplane=pnmin_,pnmax_
        call write_step(iplane,'U1',U(:,:,iplane,1))
     end do
     do iplane=pnmin_,pnmax_
        call write_step(iplane,'U2',U(:,:,iplane,2))
     end do
     do iplane=pnmin_,pnmax_
        call write_step(iplane,'U3',U(:,:,iplane,3))
     end do
     do isc=1,nscalar
        do iplane=pnmin_,pnmax_
           call write_step(iplane,SC_name(isc),SC(:,:,iplane,isc))
        end do
     end do
  end if

  if (output_pdfs) then
     call cond_map
     ! Instantaneous progress variable : y_cond(k,j,iplane)
  end if

  return
end subroutine planeStats_scatter_write
