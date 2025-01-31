module soot
  use combustion
  use chemtable
  use math
  implicit none
  
  ! Size of problem
  integer, parameter :: dim = 2
  integer :: nMoments
  integer :: nEquations
  
  ! Indices
  integer, dimension(:), pointer :: isc_mom

  ! Array of derived mean quantities
  real(WP), dimension(:,:,:), pointer :: numdens
  real(WP), dimension(:,:,:), pointer :: volfrac
  real(WP), dimension(:,:,:), pointer :: partdiam
  real(WP), dimension(:,:,:), pointer :: partaggr
  real(WP), dimension(:,:,:), pointer :: intermit

  ! Arrays of source terms (for computing statistics)
  real(WP), dimension(:,:,:), pointer :: Nsrc_nucl
  real(WP), dimension(:,:,:), pointer :: Nsrc_coag
  real(WP), dimension(:,:,:), pointer :: Nsrc_ox
  real(WP), dimension(:,:,:), pointer :: Nsrc_frag
  real(WP), dimension(:,:,:), pointer :: FVsrc_nucl
  real(WP), dimension(:,:,:), pointer :: FVsrc_cond
  real(WP), dimension(:,:,:), pointer :: FVsrc_sg
  real(WP), dimension(:,:,:), pointer :: FVsrc_ox

  ! Soot intermittency threshold
  real(WP) :: soot_int

  ! Values to monitor
  real(WP) :: min_N, max_N
  real(WP) :: min_fv,max_fv
  real(WP) :: min_dp,max_dp
  real(WP) :: min_np,max_np
  real(WP) :: min_int,max_int

contains
  
  ! Get the indices of variables from their names
  ! ---------------------------------------------
  subroutine soot_get_indices
    implicit none
    character(len=str_short) :: name
    integer :: n,isc
    
    allocate(isc_mom(nMoments+1))
    
    isc_mom = -1
    
    do isc=1,nscalar
       read(SC_name(isc),*) name
       if (name(1:2).eq.'S_') then
          select case(name(3:5))
             case ('M00')
                isc_mom(1) = isc
             case ('M10')
                isc_mom(2) = isc
             case ('M01')
                isc_mom(3) = isc
             case ('M20')
                isc_mom(4) = isc
             case ('M11')
                isc_mom(5) = isc
             case ('M02')
                isc_mom(6) = isc
             case ('N00')
                isc_mom(nMoments) = isc
             case ('MSQ')
                isc_mom(nMoments+1) = isc
             case default
                print*, SC_name(isc)
                call die('soot_get_indices: unknown soot variable')
          end select
       end if
    end do
    
    do n=1,nEquations
       if(isc_mom(n).eq.-1) call die('soot_get_indices: missing soot variable')
    end do

    return
  end subroutine soot_get_indices

  ! Soot momentum
  ! -------------
  subroutine soot_momentum
    use metric_generic
    use masks
    implicit none

    integer :: i,j,k,n

    !$OMP PARALLEL

    ! Gas + thermophoresis
    !$OMP DO
    do k=kmino_,kmaxo_
       do j=jmino_,jmaxo_
          do i=imino_,imaxo_
             if (mask_u(i,j).eq.0) &
                  rhoUt(i,j,k,isc_mom(1)) = rhoUt(i,j,k,isc_mom(1)) - 0.55_WP*sum(interp_sc_x(i,j,:)*VISCmol(i-st2:i+st1,j,k)) * &
                  sum(grad_x(i,j,:)*T(i-st2:i+st1,j,k)) / sum(interp_sc_x(i,j,:)*T(i-st2:i+st1,j,k))
             if (mask_v(i,j).eq.0) &
                  rhoVt(i,j,k,isc_mom(1)) = rhoVt(i,j,k,isc_mom(1)) - 0.55_WP*sum(interp_sc_y(i,j,:)*VISCmol(i,j-st2:j+st1,k)) * &
                  sum(grad_y(i,j,:)*T(i,j-st2:j+st1,k)) / sum(interp_sc_y(i,j,:)*T(i,j-st2:j+st1,k))
             if (mask_w(i,j).eq.0) &
                  rhoWt(i,j,k,isc_mom(1)) = rhoWt(i,j,k,isc_mom(1)) - 0.55_WP*sum(interp_sc_z(i,j,:)*VISCmol(i,j,k-st2:k+st1)) * &
                  sum(grad_z(i,j,:)*T(i,j,k-st2:k+st1)) / sum(interp_sc_z(i,j,:)*T(i,j,k-st2:k+st1))
          end do
       end do
    end do
    !$OMP END DO
    
    ! Copy for the other moments
    do n=2,nEquations
       !$OMP DO
       do k=kmino_,kmaxo_
          do j=jmino_,jmaxo_
             do i=imino_,imaxo_
                rhoUt(i,j,k,isc_mom(n)) = rhoUt(i,j,k,isc_mom(1))
                rhoVt(i,j,k,isc_mom(n)) = rhoVt(i,j,k,isc_mom(1))
                rhoWt(i,j,k,isc_mom(n)) = rhoWt(i,j,k,isc_mom(1))
             end do
          end do
       end do
       !$OMP END DO
    end do

    !$OMP END PARALLEL
       
    return
  end subroutine soot_momentum

end module soot


! ========================== !
! Initialize the soot module !
! ========================== !
subroutine soot_init
  use soot
  implicit none

  ! If not using soot model, the return
  call parser_read('Use soot model',use_soot,.false.)

  if (.not.use_soot) return

  ! Check for proper chemistry model
  if (trim(chemistry).ne.'chemtable' .and. trim(chemistry).ne.'finite chem') then
     call die('soot_init: soot requires chemtable or finite rate chemistry')
  end if

  ! Create & Start the timer
  call timing_create('soot')
  call timing_start ('soot')

  ! Get the number of moments: only four moments supported now
  call parser_read('Number of moments',nMoments)
  nEquations = nMoments
  if (nMoments.ne.4 .and. nMoments.ne.7) then
     call die('soot_init: only four or seven moments supported for now')
  end if

  ! Soot intermittency lower cut-off
  call parser_read('Intermittency threshold',soot_int,1.0e-10_WP)

  ! Solve for number density squared for double delta soot subfilter PDF
  ! if (use_sgs) then
     call parser_read('Soot subfilter PDF',use_soot_sgs)
  ! else
  !    use_soot_sgs = .false.
  ! end if
  if (use_soot_sgs) nEquations = nEquations + 1

  ! Allocate general arrays for soot
  allocate(numdens (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(volfrac (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(partdiam(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(partaggr(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
  allocate(intermit(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))

  ! Allocate arrays for source terms for statistics
  allocate(Nsrc_nucl(imino_:imaxo_,jmino_:jmaxo_,kmino:kmaxo_))
  allocate(Nsrc_coag(imino_:imaxo_,jmino_:jmaxo_,kmino:kmaxo_))
  allocate(Nsrc_ox  (imino_:imaxo_,jmino_:jmaxo_,kmino:kmaxo_))
  allocate(Nsrc_frag(imino_:imaxo_,jmino_:jmaxo_,kmino:kmaxo_))
  allocate(FVsrc_nucl(imino_:imaxo_,jmino_:jmaxo_,kmino:kmaxo_))
  allocate(FVsrc_cond(imino_:imaxo_,jmino_:jmaxo_,kmino:kmaxo_))
  allocate(FVsrc_sg  (imino_:imaxo_,jmino_:jmaxo_,kmino:kmaxo_))
  allocate(FVsrc_ox  (imino_:imaxo_,jmino_:jmaxo_,kmino:kmaxo_))

  ! Detect the variable names
  call soot_get_indices
  
  ! Initialize the HMOM part of soot
  call soot_hmom_init

  ! Soot diffusivity
  call soot_diffusivity

  ! Compute derived soot quantities
  call soot_poststep

  ! Create a new file to monitor at each iteration
  call monitor_create_file_step('soot',10)
  call monitor_set_header ( 1,'min_N', 'r')
  call monitor_set_header ( 2,'max_N', 'r')
  call monitor_set_header ( 3,'min_fv','r')
  call monitor_set_header ( 4,'max_fv','r')
  call monitor_set_header ( 5,'min_dp','r')
  call monitor_set_header ( 6,'max_dp','r')
  call monitor_set_header ( 7,'min_np','r')
  call monitor_set_header ( 8,'max_np','r')
  call monitor_set_header ( 9,'min_int','r')
  call monitor_set_header (10,'max_int','r')

  ! Stop the timer
  call timing_stop('soot')

  return
end subroutine soot_init


! ============================= !
! Pre-timestep routine for soot !
! - Compute scalar momentum     !
! - Compute scalar diffusivity  !
! ============================= !
subroutine soot_prestep
  use soot
  implicit none

  ! If not soot, then exit
  if (.not.use_soot) return

  ! Start the timer
  call timing_start('soot')

!!$  ! Compute soot momentum
!!$  call soot_momentum

  ! Compute soot diffusivity
!!$  call soot_diffusivity

  ! Stop the timer
  call timing_stop('soot')

  return
end subroutine soot_prestep


! ============================================== !
! Compute the soot source terms for the momentum !
! ============================================== !
subroutine soot_source_momentum
  use soot
  use velocity
  use metric_generic
  implicit none

  integer :: i,j,k

  ! Momentum source due to mass transfer
  !$OMP PARALLEL DO
  do k=kmin_,kmax_
     do j=jmin_,jmax_
        do i=imin_,imax_
           srcU(i,j,k) = srcU(i,j,k) + dt_uvw*0.5_WP*(Uold(i,j,k)+U(i,j,k))*sum(interp_sc_x(i,j,:)*srcP(i-st2:i+st1,j,k))
           srcV(i,j,k) = srcV(i,j,k) + dt_uvw*0.5_WP*(Vold(i,j,k)+V(i,j,k))*sum(interp_sc_y(i,j,:)*srcP(i,j-st2:j+st1,k))
           srcW(i,j,k) = srcW(i,j,k) + dt_uvw*0.5_WP*(Wold(i,j,k)+W(i,j,k))*sum(interp_sc_z(i,j,:)*srcP(i,j,k-st2:k+st1))
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine soot_source_momentum


! ============================================= !
! Compute the soot source terms for the scalars !
! ============================================= !
subroutine soot_source_scalar(SC_,srcSC_,srcP_)
  use soot
  implicit none

  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(in) :: SC_
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,nScalar), intent(inout) :: srcSC_
  real(WP), dimension(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_), intent(inout) :: srcP_

  ! Soot source terms
  call soot_hmom_source_scalar(SC_,srcSC_,srcP_)

  return
end subroutine soot_source_scalar


! Soot diffusivity
! ----------------
subroutine soot_diffusivity
  use soot
  implicit none

  integer :: i,j,k,n

  do n=1,nEquations
     !$OMP PARALLEL DO
     do k=kmino_,kmaxo_
        do j=jmino_,jmaxo_
           do i=imino_,imaxo_
              DIFF(i,j,k,isc_mom(n)) = 0.0_WP
              DIFFmol(i,j,k,isc_mom(n)) = 0.0_WP
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end do

  return
end subroutine soot_diffusivity


! ============================== !
! Post-timestep routine for soot !
! - Compute derived quantities   !
! ============================== !
subroutine soot_poststep
  use soot
  implicit none

  ! Save derived quantities and source terms for post-processing
  call soot_hmom_poststep

  return
end subroutine soot_poststep


! ======================= !
! Monitor the soot module !
! ======================= !
subroutine soot_monitor
  use soot
  use time_info
  use parallel
  implicit none

  ! If not soot, then exit
  if (.not.use_soot) return

  call soot_poststep

  ! Compute min/max
  call parallel_max( numdens,max_N)
  call parallel_max(-numdens,min_N)
  min_N = -min_N

  call parallel_max( volfrac,max_fv)
  call parallel_max(-volfrac,min_fv)
  min_fv = -min_fv

  call parallel_max( partdiam,max_dp)
  call parallel_max(-partdiam,min_dp)
  min_dp = -min_dp

  call parallel_max( partaggr,max_np)
  call parallel_max(-partaggr,min_np)
  min_np = -min_np

  call parallel_max( intermit,max_int)
  call parallel_max(-intermit,min_int)
  min_int = -min_int

  ! Transfer values to monitor
  call monitor_select_file('soot')
  call monitor_set_single_value( 1,min_N)
  call monitor_set_single_value( 2,max_N)
  call monitor_set_single_value( 3,min_fv)
  call monitor_set_single_value( 4,max_fv)
  call monitor_set_single_value( 5,min_dp)
  call monitor_set_single_value( 6,max_dp)
  call monitor_set_single_value( 7,min_np)
  call monitor_set_single_value( 8,max_np)
  call monitor_set_single_value( 9,min_int)
  call monitor_set_single_value(10,max_int)

  return
end subroutine soot_monitor
