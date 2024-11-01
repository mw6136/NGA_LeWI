! ========================================================== !
! Read input data from text files                            !
! ========================================================== !
subroutine spectrum_span_ascii_read
  use spectrum_span
  implicit none
  
  character(len=str_long) :: filename, filename2
  integer :: iunit, ifile, iprobe, istart, iend, irec
  integer :: tmp_int
  real(WP) :: tmp_real

  do ifile=1,n_monitors
     print *, 'Moving to ', trim(monitor_dirs(ifile))
     if (ifile.eq.1) then
        istart = 1
     else
        istart = sum(data_length(1:ifile-1)) + 1
     end if
     iend = sum(data_length(1:ifile))
     do iprobe=1,n_probes
        print *, 'Reading ', trim(probe_loc(iprobe))
        filename2 = trim(data_dir) // '/' // trim(monitor_dirs(ifile)) // '/' // trim(probe_loc(iprobe))
        ! RHO
        filename = trim(filename2)  // '-RHO'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, RHO(iprobe,irec,1:nz)
        end do
        close(iunit)
        ! dRHOdx
        filename = trim(filename2)  // '-dRHOdx1'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dRHOdx1(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dRHOdx2'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dRHOdx2(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dRHOdx3'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dRHOdx3(iprobe,irec,1:nz)
        end do
        close(iunit)
        ! VISC
        filename = trim(filename2)  // '-VISC'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, VISC(iprobe,irec,1:nz)
        end do
        close(iunit)
        ! P
        filename = trim(filename2)  // '-P'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, P(iprobe,irec,1:nz)
        end do
        close(iunit)
        ! U
        filename = trim(filename2)  // '-U'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, U(iprobe,irec,1:nz)
        end do
        close(iunit)
        ! V
        filename = trim(filename2)  // '-V'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, V(iprobe,irec,1:nz)
        end do
        close(iunit)
        ! W
        filename = trim(filename2)  // '-W'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, W(iprobe,irec,1:nz)
        end do
        close(iunit)
        ! dUdx
        filename = trim(filename2)  // '-dUdx11'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dUdx11(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dUdx12'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dUdx12(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dUdx13'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dUdx13(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dUdx21'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dUdx21(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dUdx22'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dUdx22(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dUdx23'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dUdx23(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dUdx31'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dUdx31(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dUdx32'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dUdx32(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dUdx33'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dUdx33(iprobe,irec,1:nz)
        end do
        close(iunit)
        ! dPdx
        filename = trim(filename2)  // '-dPdx1'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dPdx1(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dPdx2'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dPdx2(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dPdx3'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dPdx3(iprobe,irec,1:nz)
        end do
        close(iunit)
        ! dTAUdx
        filename = trim(filename2)  // '-dTAUdx11'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dTAUdx11(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dTAUdx12'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dTAUdx12(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dTAUdx13'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dTAUdx13(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dTAUdx21'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dTAUdx21(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dTAUdx22'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dTAUdx22(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dTAUdx23'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dTAUdx23(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dTAUdx31'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dTAUdx31(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dTAUdx32'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dTAUdx32(iprobe,irec,1:nz)
        end do
        close(iunit)
        filename = trim(filename2)  // '-dTAUdx33'
        open(unit=iunit, file=trim(filename), action='read')
        read(iunit,*); read(iunit,*); read(iunit,*)
        do irec=istart,iend
           read(iunit,*) tmp_int, tmp_real, dTAUdx33(iprobe,irec,1:nz)
        end do
        close(iunit)
     end do
  end do

end subroutine spectrum_span_ascii_read
