module planeStats_cond
  use planeStats

  ! Index for conditional variable
  integer :: isc_COND

  ! Gradient in conditional sample space
  real(WP), dimension(:,:), pointer :: grad1_C, grad2_C

  ! Progress variable pdf
  real(WP), dimension(:,:), pointer :: pdf_C

  ! Conditional scalar budgets
  integer :: nterms_SCc
  real(WP), dimension(:,:,:,:), pointer :: SCc

  ! Conditional velocity budgets
  integer :: nterms_Uc
  real(WP), dimension(:,:,:,:), pointer :: Uc

  ! For testing
!!$  real(WP), dimension(:,:,:), pointer :: DIFF_c, DIFF_fc_c
!!$  real(WP), dimension(:,:,:), pointer :: src_SC_c, src_SC_fc_c
!!$  real(WP), dimension(:,:,:), pointer :: dq_c

contains

  ! ========================================================== !
  ! Compute the conditional gradient (first deriv)             !
  ! ========================================================== !
  subroutine cgrad1(sol)
    implicit none
    integer :: bin
    real(WP), dimension(1:nbins_cond), intent(inout) :: sol
    real(WP), dimension(1:nbins_cond) :: data

    data = sol
    sol(1) = sum(grad1_C(1,2:3)*data(1:2))
    do bin=2,nbins_cond-1
       sol(bin) = sum(grad1_C(bin,:)*data(bin-1:bin+1))
    end do
    sol(nbins_cond) = sum(grad1_C(nbins_cond,1:2)*data(nbins_cond-1:nbins_cond))

  end subroutine cgrad1

  ! ========================================================== !
  ! Compute the conditional gradient (second deriv)            !
  ! ========================================================== !
  subroutine cgrad2(sol)
    implicit none
    integer :: bin
    real(WP), dimension(1:nbins_cond), intent(inout) :: sol
    real(WP), dimension(1:nbins_cond) :: data

    data = sol
    sol(1) = sum(grad2_C(1,:)*data(1:3))
    do bin=2,nbins_cond-1
       sol(bin) = sum(grad2_C(bin,:)*data(bin-1:bin+1))
    end do
    sol(nbins_cond) = sum(grad2_C(nbins_cond,:)*data(nbins_cond-2:nbins_cond))

  end subroutine cgrad2

end module planeStats_cond


! ========================================================== !
! Initialize the module                                      !
! ========================================================== !
subroutine planeStats_cond_init
  use planeStats_cond
  implicit none

  integer :: bin
  real(WP) :: DL

  ! Gradients in conditional sample space
  DL = 1.0_WP/real(nbins_cond,WP)
  allocate(grad1_C(1:nbins_cond,1:3))
  grad1_C = 0.0_WP
  grad1_C(1,2) = -1.0_WP/DL
  grad1_C(1,3) =  1.0_WP/DL
  do bin=2,nbins_cond-1
     grad1_C(bin,1) = -1.0_WP/(2.0_WP*DL)
     grad1_C(bin,3) =  1.0_WP/(2.0_WP*DL)
  end do
  grad1_C(nbins_cond,1) = -1.0_WP/DL
  grad1_C(nbins_cond,2) =  1.0_WP/DL

  allocate(grad2_C(1:nbins_cond,1:3))
  grad2_C = 0.0_WP
  grad2_C(:,1) =  1.0_WP/DL**2
  grad2_C(:,2) = -2.0_WP/DL**2
  grad2_C(:,3) =  1.0_WP/DL**2

  ! Progress variable pdf
  allocate(pdf_C(1:nbins_cond,pnmin_:pnmax_))

  ! Conditional scalar budgets
  nterms_SCc = 8
  allocate(SCc(1:nbins_cond,pnmin_:pnmax_,1:nterms_SCc,1:nscalar))
  SCc = 0.0_WP

  ! Conditional velocity budgets
  nterms_Uc = 8
  allocate(Uc(1:nbins_cond,pnmin_:pnmax_,1:nterms_Uc,1:2))
  Uc = 0.0_WP

  if (trim(cond_var).eq.'PROG') then
     isc_COND = isc_PROG
  else
     isc_COND = isc_ZMIX
  end if

  ! Initialize the finitechem module
  call planeStats_finitechem_init

  ! for testing
!!$  allocate(DIFF_c   (1:nbins_cond,pnmin_:pnmax_,1:nscalar))
!!$  allocate(DIFF_fc_c(1:nbins_cond,pnmin_:pnmax_,1:nscalar))
!!$  allocate(src_SC_c   (1:nbins_cond,pnmin_:pnmax_,1:nscalar-1))
!!$  allocate(src_SC_fc_c(1:nbins_cond,pnmin_:pnmax_,1:nscalar-1))
!!$  allocate(dq_c(1:nbins_cond,pnmin_:pnmax_,1:nscalar))

  return
end subroutine planeStats_cond_init


! ========================================================== !
! Compute the conditional budgets                            !
! ========================================================== !
subroutine planeStats_cond_compute
  use planeStats_cond
  implicit none

  integer :: iplane, j, k, ii, jj, isc, bin
  real(WP), dimension(1:nz,1:ny) :: tmpxy
  real(WP), dimension(1:nbins_cond) :: tmpc1, tmpc2

  ! Need scalar gradients - compute diffusivities below
!  if (.not.use_qSC) return

  ! Compute the conditional mapping
  call cond_map

  ! Evaluate diffusivity
  call planeStats_finitechem_diffusivity

  ! Evaluate source terms
  call planeStats_finitechem_source


!!$  if (use_qSC) then
!!$     do isc=1,nscalar
!!$        do iplane=pnmin_,pnmax_
!!$           call cmean(iplane, diff_flux_div(:,:,iplane,isc), dq_c(:,iplane,isc), icount_qSC)
!!$        end do
!!$     end do
!!$  end if
!!$  ! for testing
!!$  if (use_qSC) then
!!$     do isc=1,nscalar
!!$        do iplane=pnmin_,pnmax_
!!$           call cmean(iplane, DIFF(:,:,iplane,isc), DIFF_c(:,iplane,isc), icount_qSC)
!!$        end do
!!$     end do
!!$  end if
!!$  do isc=1,nscalar
!!$     do iplane=pnmin_,pnmax_
!!$        call cmean(iplane, DIFF_fc(:,:,iplane,isc), DIFF_fc_c(:,iplane,isc), icount_all)
!!$     end do
!!$  end do
!!$  if (use_qSC) then
!!$     do isc=1,nscalar-1
!!$        do iplane=pnmin_,pnmax_
!!$           call cmean(iplane, src_SC(:,:,iplane,isc), src_SC_c(:,iplane,isc), icount_qSC)
!!$        end do
!!$     end do
!!$  end if
!!$  do isc=1,nscalar-1
!!$     do iplane=pnmin_,pnmax_
!!$        call cmean(iplane, src_SC_fc(:,:,iplane,isc), src_SC_fc_c(:,iplane,isc), icount_all)
!!$     end do
!!$  end do

!!$  ! for testing - reduce number of data fields read in mpi - do same for DIFF
!!$  if (ntime_curr.eq.1) then
!!$     print *, 'saving SC_old'
!!$     SC_old = SC
!!$  else
!!$     print *, 'ntime=',ntime_curr,', computing source terms at n+1'
!!$     SC = 0.5_WP*(SC+SC_old)
!!$     call planeStats_finitechem_source


!!$     call planeStats_finitechem_diffusivity
!!$  do isc=1,nscalar
!!$     do iplane=pnmin_,pnmax_
!!$        call cmean(iplane, DIFF(:,:,iplane,isc), DIFF_c(:,iplane,isc), icount_dSC)
!!$        call cmean(iplane, DIFF_fc(:,:,iplane,isc), DIFF_fc_c(:,iplane,isc), icount_dSC)
!!$     end do
!!$  end do
!!$  end if

  ! Count the distribution of the conditional variable
  do iplane=pnmin_,pnmax_
     do j=1,ny
        do k=1,nz
           bin = bin_index(k,j,iplane)
           pdf_C(bin,iplane) = pdf_C(bin,iplane) + 1.0_WP
        end do
     end do
  end do


  ! Update the conditional means of terms in the CMC scalar equations
  do isc=1,nscalar
     ! TERM 1
     if (use_dSC) then
        do iplane=pnmin_,pnmax_
           tmpxy = 0.0_WP
           do jj=1,3
              tmpxy = tmpxy - RHO(:,:,iplane)*( &
                   + U(:,:,iplane,jj)*dSCdx(:,:,iplane,isc,jj) &
                   + SC(:,:,iplane,isc)*dUdx(:,:,iplane,jj,jj) ) &
                   - U(:,:,iplane,jj)*SC(:,:,iplane,isc)*dRHOdx(:,:,iplane,jj)
           end do
           call cmean(iplane, tmpxy, SCc(:,iplane,1,isc), icount_dSC)
        end do
        do iplane=pnmin_,pnmax_
           tmpxy = 0.0_WP
           do jj=1,3
              ! Negative from def. of C
              tmpxy = tmpxy - U(:,:,iplane,jj)*dSCdx(:,:,iplane,isc_COND,jj)
           end do
           tmpxy = tmpxy * RHO(:,:,iplane)*SC(:,:,iplane,isc)/(prog_upper-prog_lower)
           call cmean(iplane, tmpxy, SCc(:,iplane,8,isc), icount_dSC)
        end do
     end if
     
     ! TERM 2 - diffusive flux divergence
     if (use_qSC) then
        do iplane=pnmin_,pnmax_
          call cmean(iplane, diff_flux_div(:,:,iplane,isc), SCc(:,iplane,2,isc), icount_qSC)
       end do
    end if

    ! TERM 3 - scalar source term
    do iplane=pnmin_,pnmax_
       call cmean(iplane, src_SC_fc(:,:,iplane,isc), SCc(:,iplane,3,isc), icount_all)
    end do

     ! TERM 4
     if (use_dSC) then
        do iplane=pnmin_,pnmax_
           tmpxy = 0.0_WP
           do jj=1,3
              tmpxy = tmpxy - dSCdx(:,:,iplane,isc_COND,jj)**2
           end do
           tmpxy = tmpxy * SC(:,:,iplane,isc)*DIFF_fc(:,:,iplane,isc_COND)/(prog_upper-prog_lower)**2
           call cmean(iplane, tmpxy, SCc(:,iplane,4,isc), icount_dSC)
        end do
     end if

     ! TERM 5
     if (use_dSC) then
        do iplane=pnmin_,pnmax_
           tmpxy = 0.0_WP
           do jj=1,3
              tmpxy = tmpxy + dSCdx(:,:,iplane,isc,jj)*dSCdx(:,:,iplane,isc_COND,jj)
           end do
           tmpxy = tmpxy * DIFF_fc(:,:,iplane,isc_COND)/(prog_upper-prog_lower)
           call cmean(iplane, tmpxy, SCc(:,iplane,5,isc), icount_dSC)
        end do
     end if

     ! TERM 6

     ! TERM 7 - progress variable source 
     if (trim(cond_var).eq.'PROG') then
        do iplane=pnmin_,pnmax_
           tmpxy = -SC(:,:,iplane,isc)*src_SC_fc(:,:,iplane,isc_PROG)/(prog_upper-prog_lower)
           call cmean(iplane, tmpxy, SCc(:,iplane,7,isc), icount_all)
        end do
     end if
        
  end do
  

  ! Update the conditional means of terms in the CMC velocity equations
  do ii=1,2
     ! TERM 1
     if (use_dSC) then
        do iplane=pnmin_,pnmax_
           tmpxy = 0.0_WP
           do jj=1,3
              tmpxy = tmpxy - RHO(:,:,iplane)*( &
                   + U(:,:,iplane,jj)*dUdx(:,:,iplane,ii,jj) &
                   + U(:,:,iplane,ii)*dUdx(:,:,iplane,jj,jj) ) &
                   - U(:,:,iplane,ii)*U(:,:,iplane,jj)*dRHOdx(:,:,iplane,jj)
           end do
           call cmean(iplane, tmpxy, Uc(:,iplane,1,ii), icount_dSC)
        end do
        do iplane=pnmin_,pnmax_
           tmpxy = 0.0_WP
           do jj=1,3
              ! Negative from def. of C
              tmpxy = tmpxy - U(:,:,iplane,jj)*dSCdx(:,:,iplane,isc_COND,jj)
           end do
           tmpxy = tmpxy * RHO(:,:,iplane)*U(:,:,iplane,ii)/(prog_upper-prog_lower)
           call cmean(iplane, tmpxy, Uc(:,iplane,8,ii), icount_dSC)
        end do
     end if

     ! TERM 2 - PG
     do iplane=pnmin_,pnmax_
        call cmean(iplane, dPdx(:,:,iplane,ii), Uc(:,iplane,2,ii), icount_all)
     end do

     ! TERM 3 - stress tensor gradient
     do iplane=pnmin_,pnmax_
        tmpxy = 0.0_WP
        do jj=1,3
           tmpxy = tmpxy + dTAUdx(:,:,iplane,ii,jj)
        end do
        call cmean(iplane, tmpxy, Uc(:,iplane,3,ii), icount_all)
     end do

     ! TERM 4
     if (use_dSC) then
        do iplane=pnmin_,pnmax_
           tmpxy = 0.0_WP
           do jj=1,3
              tmpxy = tmpxy - dSCdx(:,:,iplane,isc_COND,jj)**2
           end do
           tmpxy = tmpxy * U(:,:,iplane,ii)*DIFF_fc(:,:,iplane,isc_COND)/(prog_upper-prog_lower)**2
           call cmean(iplane, tmpxy, Uc(:,iplane,4,ii), icount_dSC)
        end do
     end if

     !! DIFF is RHO*D

     ! TERM 5
     if (use_dSC) then
        do iplane=pnmin_,pnmax_
           tmpxy = 0.0_WP
           do jj=1,3
              tmpxy = tmpxy + dUdx(:,:,iplane,ii,jj)*dSCdx(:,:,iplane,isc_COND,jj)
           end do
           tmpxy = tmpxy * DIFF_fc(:,:,iplane,1)
           call cmean(iplane, tmpxy, Uc(:,iplane,5,ii), icount_dSC)
        end do
     end if

     ! TERM 6

     ! TERM 7 - source for reacting conditional variable
     if (trim(cond_var).eq.'PROG') then
        do iplane=pnmin_,pnmax_
           tmpxy = -RHO(:,:,iplane)*U(:,:,iplane,ii)*src_SC_fc(:,:,iplane,isc_PROG)/(prog_upper-prog_lower)
           call cmean(iplane, tmpxy, Uc(:,iplane,7,ii), icount_all)
        end do
     end if
  end do

  return
end subroutine planeStats_cond_compute


! ========================================================== !
! Finalize the conditional budgets                           !
! ========================================================== !
subroutine planeStats_cond_finalize
  use planeStats_cond
  implicit none

  integer :: iplane, ii, jj, isc

  ! Finalize the pdf
  pdf_C = pdf_C*real(nbins_cond,WP)/real(ny*nz*ntime_curr,WP)

  ! Scalars
!!$  do isc=1,nscalar
!!$     do ii=4,nterms_SCc
!!$        do iplane=pnmin_,pnmax_
!!$           SCc(:,iplane,ii,isc) = SCc(:,iplane,ii,isc)*pdf_C(:,iplane)
!!$        end do
!!$     end do
!!$     do iplane=pnmin_,pnmax_
!!$        call cgrad2(SCc(:,iplane,4,isc))
!!$        call cgrad1(SCc(:,iplane,5,isc))
!!$        call cgrad1(SCc(:,iplane,7,isc))
!!$        call cgrad1(SCc(:,iplane,8,isc))
!!$     end do
!!$     do ii=4,nterms_SCc
!!$        do iplane=pnmin_,pnmax_
!!$           SCc(:,iplane,ii,isc) = SCc(:,iplane,ii,isc)/pdf_C(:,iplane)
!!$        end do
!!$     end do
!!$  end do

  ! Momentum
  do jj=1,2
     do ii=4,nterms_Uc
        do iplane=pnmin_,pnmax_
           Uc(:,iplane,ii,jj) = Uc(:,iplane,ii,jj)*pdf_C(:,iplane)
        end do
     end do
     do iplane=pnmin_,pnmax_
        call cgrad2(Uc(:,iplane,4,jj))
        call cgrad1(Uc(:,iplane,5,jj))
        call cgrad1(Uc(:,iplane,7,jj))
        call cgrad1(Uc(:,iplane,8,jj))
     end do
     do ii=4,nterms_Uc
        do iplane=pnmin_,pnmax_
           Uc(:,iplane,ii,jj) = Uc(:,iplane,ii,jj)/pdf_C(:,iplane)
        end do
     end do
  end do

  return
end subroutine planeStats_cond_finalize


! ========================================================== !
! Output the conditional budgets                             !
! ========================================================== !
subroutine planeStats_cond_output
  use planeStats_cond
  implicit none

  integer :: iunit, iplane, ii, jj, isc, j, ierr
  
  ! Arrays to collect
  real(WP), dimension(:,:), pointer :: o_pdf_C
  real(WP), dimension(:,:,:,:), pointer :: o_SCc
  real(WP), dimension(:,:,:,:), pointer :: o_Uc

!!$  real(WP), dimension(:,:,:), pointer :: o_DIFF_c, o_DIFF_fc_c
!!$  real(WP), dimension(:,:,:), pointer :: o_src_SC_c, o_src_SC_fc_c
!!$  real(WP), dimension(:,:,:), pointer :: o_dq_c

  ! Allocate them
  allocate(o_pdf_C(1:nplanes,1:nbins_cond))
  allocate(o_SCc  (1:nplanes,1:nbins_cond,1:nterms_SCc,1:nscalar))
  allocate(o_Uc   (1:nplanes,1:nbins_cond,1:nterms_Uc,1:2))

!!$  allocate(o_DIFF_c  (1:nplanes,1:nbins_cond,1:nscalar))
!!$  allocate(o_DIFF_fc_c(1:nplanes,1:nbins_cond,1:nscalar))
!!$  allocate(o_src_SC_c  (1:nplanes,1:nbins_cond,1:nscalar-1))
!!$  allocate(o_src_SC_fc_c(1:nplanes,1:nbins_cond,1:nscalar-1))
!!$  allocate(o_dq_c  (1:nplanes,1:nbins_cond,1:nscalar))

  ! Collect from the processes
  do j=1,nbins_cond
     call MPI_GATHER(pdf_C(j,:), nplanes_, MPI_REAL_WP, o_pdf_C(:,j), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
  end do
  do isc=1,nscalar
     do ii=1,nterms_SCc
        do j=1,nbins_cond
           call MPI_GATHER(SCc(j,:,ii,isc), nplanes_, MPI_REAL_WP, o_SCc(:,j,ii,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do
  do jj=1,2
     do ii=1,nterms_Uc
        do j=1,nbins_cond
           call MPI_GATHER(Uc(j,:,ii,jj), nplanes_, MPI_REAL_WP, o_Uc(:,j,ii,jj), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
        end do
     end do
  end do

!!$  do isc=1,nscalar
!!$     do j=1,nbins_cond
!!$        call MPI_GATHER(DIFF_c(j,:,isc), nplanes_, MPI_REAL_WP, o_DIFF_c(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
!!$     end do
!!$  end do
!!$  do isc=1,nscalar
!!$     do j=1,nbins_cond
!!$        call MPI_GATHER(DIFF_fc_c(j,:,isc), nplanes_, MPI_REAL_WP, o_DIFF_fc_c(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
!!$     end do
!!$  end do
!!$  do isc=1,nscalar-1
!!$     do j=1,nbins_cond
!!$        call MPI_GATHER(src_SC_c(j,:,isc), nplanes_, MPI_REAL_WP, o_src_SC_c(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
!!$     end do
!!$  end do
!!$  do isc=1,nscalar-1
!!$     do j=1,nbins_cond
!!$        call MPI_GATHER(src_SC_fc_c(j,:,isc), nplanes_, MPI_REAL_WP, o_src_SC_fc_c(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
!!$     end do
!!$  end do
!!$  do isc=1,nscalar
!!$     do j=1,nbins_cond
!!$        call MPI_GATHER(dq_c(j,:,isc), nplanes_, MPI_REAL_WP, o_dq_c(:,j,isc), nplanes_, MPI_REAL_WP, iroot-1, MPI_COMM_WORLD, ierr)
!!$     end do
!!$  end do

  ! Only the root process writes output
  if (irank.ne.iroot) return
  
  ! Output the files

  ! Scalar equation budgets - Conditional space
  open(unit=iunit, file=trim(output_name)//'_cond_SCc', action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        write(iunit,'(ES22.13)',advance='no') o_pdf_C(iplane,j)
        do isc=1,nscalar
           do ii=1,nterms_SCc
              if (iplane.eq.nplanes.and.isc.eq.nscalar.and.ii.eq.nterms_SCc) then
                 write(iunit,'(ES22.13)') o_SCc(iplane,j,ii,isc)
              else
                 write(iunit,'(ES22.13)',advance='no') o_SCc(iplane,j,ii,isc)
              end if
           end do
        end do
     end do
  end do
  close(iunit)

  ! Momentum equation budgets - Conditional space
  open(unit=iunit, file=trim(output_name)//'_cond_Uc1', action='write')
  do j=1,nbins_cond
     do iplane=1,nplanes
        if (iplane.eq.1) then
           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
        end if
        write(iunit,'(ES22.13)',advance='no') o_pdf_C(iplane,j)
        do jj=1,2
           do ii=1,nterms_Uc
              if (iplane.eq.nplanes.and.jj.eq.2.and.ii.eq.nterms_Uc) then
                 write(iunit,'(ES22.13)') o_Uc(iplane,j,ii,jj)
              else
                 write(iunit,'(ES22.13)',advance='no') o_Uc(iplane,j,ii,jj)
              end if
           end do
        end do
     end do
  end do
  close(iunit)


!!$  open(unit=iunit, file=trim(output_name)//'_cond_dq', action='write')
!!$  do j=1,nbins_cond
!!$     do iplane=1,nplanes
!!$        if (iplane.eq.1) then
!!$           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
!!$        end if
!!$        do isc=1,nscalar
!!$           if (iplane.eq.nplanes.and.isc.eq.nscalar) then 
!!$              write(iunit,'(2ES22.13)') o_dq_c(iplane,j,isc)
!!$           else
!!$              write(iunit,'(2ES22.13)',advance='no') o_dq_c(iplane,j,isc)
!!$           end if
!!$        end do
!!$     end do
!!$  end do

!!$  ! Diffusivity -- NGA output vs planestats computed
!!$  open(unit=iunit, file=trim(output_name)//'_cond_DIFF', action='write')
!!$  do j=1,nbins_cond
!!$     do iplane=1,nplanes
!!$        if (iplane.eq.1) then
!!$           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
!!$        end if
!!$        do isc=1,nscalar
!!$           if (iplane.eq.nplanes.and.isc.eq.nscalar) then 
!!$              write(iunit,'(2ES22.13)') o_DIFF_c(iplane,j,isc), o_DIFF_fc_c(iplane,j,isc)
!!$           else
!!$              write(iunit,'(2ES22.13)',advance='no') o_DIFF_c(iplane,j,isc), o_DIFF_fc_c(iplane,j,isc)
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  ! Source term -- NGA output vs planestats computed
!!$  open(unit=iunit, file=trim(output_name)//'_cond_src_SC', action='write')
!!$  do j=1,nbins_cond
!!$     do iplane=1,nplanes
!!$        if (iplane.eq.1) then
!!$           write(iunit,'(ES22.13)',advance='no') bins_cond(j)
!!$        end if
!!$        do isc=1,nscalar-1
!!$           if (iplane.eq.nplanes.and.isc.eq.nscalar-1) then 
!!$              write(iunit,'(2ES22.13)') o_src_SC_c(iplane,j,isc), o_src_SC_fc_c(iplane,j,isc)
!!$           else
!!$              write(iunit,'(2ES22.13)',advance='no') o_src_SC_c(iplane,j,isc), o_src_SC_fc_c(iplane,j,isc)
!!$           end if
!!$        end do
!!$     end do
!!$  end do


  ! These are instantaneous values below -- for debugging
!!$  ! Source term -- NGA output vs planestats computed
!!$  open(unit=iunit, file=trim(output_name)//'_src_SC', action='write')
!!$  do j=1,ny
!!$     do iplane=1,nplanes
!!$        if (iplane.eq.1) then
!!$           write(iunit,'(ES22.13)',advance='no') y(j)
!!$        end if
!!$        do isc=1,nscalar-1
!!$           if (iplane.eq.nplanes.and.isc.eq.nscalar-1) then 
!!$              write(iunit,'(2ES22.13)') src_SC(1,j,iplane,isc), src_SC_fc(1,j,iplane,isc)
!!$           else
!!$              write(iunit,'(2ES22.13)',advance='no') src_SC(1,j,iplane,isc), src_SC_fc(1,j,iplane,isc)
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!!$  ! Source term -- NGA output vs planestats computed
!!$  open(unit=iunit, file=trim(output_name)//'_src_DIFF', action='write')
!!$  do j=1,ny
!!$     do iplane=1,nplanes
!!$        if (iplane.eq.1) then
!!$           write(iunit,'(ES22.13)',advance='no') y(j)
!!$        end if
!!$        do isc=1,nscalar
!!$           if (iplane.eq.nplanes.and.isc.eq.nscalar) then 
!!$              write(iunit,'(2ES22.13)') DIFF(1,j,iplane,isc), DIFF_fc(1,j,iplane,isc)
!!$           else
!!$              write(iunit,'(2ES22.13)',advance='no') DIFF(1,j,iplane,isc), DIFF_fc(1,j,iplane,isc)
!!$           end if
!!$        end do
!!$     end do
!!$  end do

  return
end subroutine planeStats_cond_output




!! OLD STUFF

! CMC SCALAR AND VELOCITY EQUATIONS FROM 9/26 NOTES
!!$
!!$  ! Update the conditional means of terms in the CMC scalar equations
!!$  do isc=1,nscalar
!!$     ! TERM 1
!!$     if (use_dSC) then
!!$        do iplane=pnmin_,pnmax_
!!$           tmpxy = 0.0_WP
!!$           do jj=1,3
!!$              tmpxy = tmpxy - RHO(:,:,iplane)*( &
!!$                   + U(:,:,iplane,jj)*dSCdx(:,:,iplane,isc,jj) &
!!$                   + SC(:,:,iplane,isc)*dUdx(:,:,iplane,jj,jj) ) &
!!$                   - U(:,:,iplane,jj)*SC(:,:,iplane,isc)*dRHOdx(:,:,iplane,jj)
!!$           end do
!!$           call cmean(iplane, tmpxy, SCc(:,iplane,1,isc), icount_dSC)
!!$        end do
!!$        do iplane=pnmin_,pnmax_
!!$           tmpxy = 0.0_WP
!!$           do jj=1,3
!!$              ! Negative from def. of C
!!$              tmpxy = tmpxy - U(:,:,iplane,jj)*dSCdx(:,:,iplane,isc_COND,jj)
!!$           end do
!!$           tmpxy = tmpxy * RHO(:,:,iplane)*SC(:,:,iplane,isc)/(prog_upper-prog_lower)
!!$           call cmean(iplane, tmpxy, SCc(:,iplane,8,isc), icount_dSC)
!!$        end do
!!$     end if
!!$     
!!$     ! TERM 2 - diffusive flux divergence
!!$     if (use_qSC) then
!!$        do iplane=pnmin_,pnmax_
!!$          call cmean(iplane, diff_flux_div(:,:,iplane,isc), SCc(:,iplane,2,isc), icount_qSC)
!!$       end do
!!$    end if
!!$
!!$    ! TERM 3 - scalar source term
!!$    do iplane=pnmin_,pnmax_
!!$       call cmean(iplane, src_SC_fc(:,:,iplane,isc), SCc(:,iplane,3,isc), icount_all)
!!$    end do
!!$
!!$     ! TERM 4
!!$     if (use_dSC) then
!!$        do iplane=pnmin_,pnmax_
!!$           tmpxy = 0.0_WP
!!$           do jj=1,3
!!$              tmpxy = tmpxy - dSCdx(:,:,iplane,isc_COND,jj)**2
!!$           end do
!!$           tmpxy = tmpxy * SC(:,:,iplane,isc)*DIFF_fc(:,:,iplane,isc_COND)/(prog_upper-prog_lower)**2
!!$           call cmean(iplane, tmpxy, SCc(:,iplane,4,isc), icount_dSC)
!!$        end do
!!$     end if
!!$
!!$     ! TERM 5
!!$     if (use_dSC) then
!!$        do iplane=pnmin_,pnmax_
!!$           tmpxy = 0.0_WP
!!$           do jj=1,3
!!$              tmpxy = tmpxy + dSCdx(:,:,iplane,isc,jj)*dSCdx(:,:,iplane,isc_COND,jj)
!!$           end do
!!$           tmpxy = tmpxy * DIFF_fc(:,:,iplane,isc_COND)/(prog_upper-prog_lower)
!!$           call cmean(iplane, tmpxy, SCc(:,iplane,5,isc), icount_dSC)
!!$        end do
!!$     end if
!!$
!!$     ! TERM 6
!!$
!!$     ! TERM 7 - progress variable source 
!!$     if (trim(cond_var).eq.'PROG') then
!!$        do iplane=pnmin_,pnmax_
!!$           tmpxy = -SC(:,:,iplane,isc)*src_SC_fc(:,:,iplane,isc_PROG)/(prog_upper-prog_lower)
!!$           call cmean(iplane, tmpxy, SCc(:,iplane,7,isc), icount_all)
!!$        end do
!!$     end if
!!$        
!!$  end do
!!$
!!$  do ii=1,2
!!$     ! TERM 1
!!$     if (use_dSC) then
!!$        do iplane=pnmin_,pnmax_
!!$           tmpxy = 0.0_WP
!!$           do jj=1,3
!!$              tmpxy = tmpxy - RHO(:,:,iplane)*( &
!!$                   + U(:,:,iplane,jj)*dUdx(:,:,iplane,ii,jj) &
!!$                   + U(:,:,iplane,ii)*dUdx(:,:,iplane,jj,jj) ) &
!!$                   - U(:,:,iplane,ii)*U(:,:,iplane,jj)*dRHOdx(:,:,iplane,jj)
!!$           end do
!!$           call cmean(iplane, tmpxy, Uc(:,iplane,1,ii), icount_dSC)
!!$        end do
!!$        do iplane=pnmin_,pnmax_
!!$           tmpxy = 0.0_WP
!!$           do jj=1,3
!!$              ! Negative from def. of C
!!$              tmpxy = tmpxy - U(:,:,iplane,jj)*dSCdx(:,:,iplane,isc_COND,jj)
!!$           end do
!!$           tmpxy = tmpxy * RHO(:,:,iplane)*U(:,:,iplane,ii)/(prog_upper-prog_lower)
!!$           call cmean(iplane, tmpxy, Uc(:,iplane,8,ii), icount_dSC)
!!$        end do
!!$     end if
!!$
!!$     ! TERM 2 - PG
!!$     do iplane=pnmin_,pnmax_
!!$        call cmean(iplane, dPdx(:,:,iplane,ii), Uc(:,iplane,2,ii), icount_all)
!!$     end do
!!$
!!$     ! TERM 3 - stress tensor gradient
!!$     do iplane=pnmin_,pnmax_
!!$        tmpxy = 0.0_WP
!!$        do jj=1,3
!!$           tmpxy = tmpxy + dTAUdx(:,:,iplane,ii,jj)
!!$        end do
!!$        call cmean(iplane, tmpxy, Uc(:,iplane,3,ii), icount_all)
!!$     end do
!!$
!!$     ! TERM 4
!!$     if (use_dSC) then
!!$        do iplane=pnmin_,pnmax_
!!$           tmpxy = 0.0_WP
!!$           do jj=1,3
!!$              tmpxy = tmpxy - dSCdx(:,:,iplane,isc_COND,jj)**2
!!$           end do
!!$           tmpxy = tmpxy * U(:,:,iplane,ii)*DIFF_fc(:,:,iplane,isc_COND)/(prog_upper-prog_lower)**2
!!$           call cmean(iplane, tmpxy, Uc(:,iplane,4,ii), icount_dSC)
!!$        end do
!!$     end if
!!$
!!$     !! DIFF is RHO*D
!!$
!!$     ! TERM 5
!!$     if (use_dSC) then
!!$        do iplane=pnmin_,pnmax_
!!$           tmpxy = 0.0_WP
!!$           do jj=1,3
!!$              tmpxy = tmpxy + dUdx(:,:,iplane,ii,jj)*dSCdx(:,:,iplane,isc_COND,jj)
!!$           end do
!!$           tmpxy = tmpxy * DIFF_fc(:,:,iplane,1)
!!$           call cmean(iplane, tmpxy, Uc(:,iplane,5,ii), icount_dSC)
!!$        end do
!!$     end if
!!$
!!$     ! TERM 6
!!$
!!$     ! TERM 7 - source for reacting conditional variable
!!$     if (trim(cond_var).eq.'PROG') then
!!$        do iplane=pnmin_,pnmax_
!!$           tmpxy = -RHO(:,:,iplane)*U(:,:,iplane,ii)*src_SC_fc(:,:,iplane,isc_PROG)/(prog_upper-prog_lower)
!!$           call cmean(iplane, tmpxy, Uc(:,iplane,7,ii), icount_all)
!!$        end do
!!$     end if
!!$  end do
