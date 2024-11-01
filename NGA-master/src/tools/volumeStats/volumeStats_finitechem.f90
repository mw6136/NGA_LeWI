module volumeStats_finitechem
  use volumeStats

  integer :: Nreactions
  real(WP) :: W_mix
  real(WP), parameter :: R_cst=8314.34_WP
  real(WP), parameter :: P_thermo=1.0132e5_WP
  real(WP), dimension(:), pointer :: Cp_sp, h_sp, Lewis
  real(WP),dimension(:),pointer :: kOverEps, mucoeff

contains
  
  ! ======================================================= !
  ! Compute the average molecular weight in a cell(kg/kmol) !
  ! ======================================================= !
  subroutine finitechem_W(scalar, Wmix)
    implicit none
    real(WP), intent(out) :: Wmix
    real(WP), dimension(1:N_tot), intent(in) :: scalar
    real(WP) :: tmp
    integer :: i

    Wmix = 0.0_WP
    tmp  = 0.0_WP
    do i=1,N_tot
       Wmix = Wmix + scalar(i)/W_sp(i)
       tmp  = tmp  + scalar(i)
    end do
    if (Wmix.gt.0.0_WP) then
       Wmix = tmp/Wmix
    else
       Wmix = 0.0_WP
    end if

    return
  end subroutine finitechem_W

end module volumeStats_finitechem


! ========================================================== !
! Initialize the module                                      !
! ========================================================== !
subroutine volumeStats_finitechem_init
  use volumeStats_finitechem
  implicit none

  allocate(Cp_sp(1:N_tot))
  allocate(h_sp (1:N_tot))
  call GETNREACTIONS(Nreactions)

  ! kOverOmega is calculated for all the species
  allocate(kOverEps(1:N_tot))
  ! mucoeff is calculated for all the species
  allocate(mucoeff(1:N_tot))
  call GETMUCOEFF(mucoeff)
  call GETKOVEREPS(kOverEps)

  allocate(Lewis(1:N_tot)) ! Can read these in from input file
  Lewis = 1.0_WP
  Lewis(1) = 1.28645  ! N2
  Lewis(2) = 0.190182 ! H
  Lewis(3) = 1.11965  ! O2
  Lewis(4) = 0.729858 ! O
  Lewis(5) = 0.743361 ! OH
  Lewis(6) = 0.315702 ! H2
  Lewis(7) = 0.858191 ! H2O
  Lewis(8) = 1.12405  ! HO2
  Lewis(9) = 1.13135  ! H2O2

  return
end subroutine volumeStats_finitechem_init


! ========================================================== !
! Compute chemical source terms                              !
! ========================================================== !
subroutine volumeStats_finitechem_source
  use volumeStats_finitechem
  implicit none

  integer :: i,j,k,isc
  real(WP), dimension(1:Nreactions) :: K_rxn, omega_rxn, M_rxn
  real(WP), dimension(1:N_tot) :: scalar, conc, conc_dot
  real(WP) :: Cp_mix, rho_gp, temp, tmp

  !$OMP THREADPRIVATE(Cp_sp,h_sp,W_mix)
  !$OMP PARALLEL DO PRIVATE(j,i,isc,K_rxn,omega_rxn,M_rxn,scalar,conc,conc_dot,Cp_mix,rho_gp,temp,tmp)
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           ! Get W and Cp of the mixture
           call finitechem_W(SC(i,j,k,1:N_tot), W_mix)
           temp = min(max(SC(i,j,k,isc_T),0.0_WP),5000.0_WP)
           call COMPTHERMODATA(h_sp, Cp_sp, temp)
           scalar = min(max(SC(i,j,k,1:N_tot),0.0_WP),1.0_WP)
           Cp_mix = sum(scalar*Cp_sp)
           rho_gp = P_thermo*W_mix/(R_cst*temp)
           
           ! Concentrations
           do isc=1,N_tot
              conc(isc) = rho_gp*scalar(isc)/W_sp(isc)
           end do
           
           ! Obtain source terms
           call PRODRATES(conc_dot, omega_rxn, K_rxn, conc, M_rxn, temp, P_thermo)
           
           ! Concentration to mass fraction
           do isc=1,N_tot
              src_SC(i,j,k,isc) = conc_dot(isc)*W_sp(isc)
           end do
           
           ! Temperature source for isobaric system
           tmp = 0.0_WP
           do isc=1,N_tot
              tmp = tmp - h_sp(isc)*src_SC(i,j,k,isc)
           end do
           src_SC(i,j,k,isc_T) = tmp/Cp_mix
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine volumeStats_finitechem_source


! ========================================================== !
! Compute the diffusivity using transport data from FM       !
! ========================================================== !
subroutine volumeStats_finitechem_diffusivity
  use volumeStats_finitechem
  implicit none

  integer :: i,j,k,isc
  real(WP), dimension(1:N_tot) :: scalar
  real(WP) :: temp,sum_1,sum_2,lambda, Cp_mix, omegamu, cond, eta

  !$OMP THREADPRIVATE(Cp_sp,h_sp,W_mix)
  !$OMP PARALLEL DO PRIVATE(j,i,isc,scalar,cond,eta,temp,sum_1,sum_2,lambda,Cp_mix,omegamu)
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           ! Get mixture molar mass
           scalar = min(max(SC(i,j,k,1:N_tot),0.0_WP),1.0_WP)
           call finitechem_W(SC(i,j,k,1:N_tot), W_mix)
           temp = min(max(SC(i,j,k,isc_T),0.0_WP),5000.0_WP)
           call COMPTHERMODATA(h_sp,Cp_sp,temp)
           sum_1 = 0.0_WP
           sum_2 = 0.0_WP
           do isc=1,N_tot
              call volumeStats_finitechem_omegamu(temp*kOverEps(isc),omegamu)
              eta  = mucoeff(isc)*sqrt(temp)/omegamu
              cond = eta * (Cp_sp(isc)+1.25_WP*R_cst/W_sp(isc))
              sum_1 = sum_1 + (SC(i,j,k,isc)*W_mix/(cond*W_sp(isc)))
              sum_2 = sum_2 + (SC(i,j,k,isc)*cond*W_mix/W_sp(isc))
           end do
           lambda = 0.5_WP*(sum_2+(1.0_WP/sum_1))
           ! Average Cp
           Cp_mix = sum(scalar*Cp_sp)
           do isc=1,N_tot
              DIFF(i,j,k,isc) = (lambda/Cp_mix)/Lewis(isc)
           end do
           DIFF(i,j,k,isc_T) = lambda/Cp_mix
           DIFF(i,j,k,isc_ZMIX) = lambda/Cp_mix
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine volumeStats_finitechem_diffusivity


! ========================================================= !
! Compute the viscosity by fitting using FlameMaster method !
! ========================================================= !
subroutine volumeStats_finitechem_viscosity
  use volumeStats_finitechem
  implicit none

  integer  :: i,j,k,nsc,j_sc,k_sc
  real(WP) :: temp,temp1,sum_1,sum_2
  real(WP),dimension(N_tot) :: eta
  real(WP),dimension(N_tot,N_tot) :: phi
  real(WP) :: omegamu

  !$OMP THREADPRIVATE(W_mix)
  !$OMP PARALLEL DO PRIVATE(j,i,nsc,j_sc,k_sc,eta,phi,omegamu,temp,temp1,sum_1,sum_2)
  do k=kmino_,kmaxo_
     do j=jmino_,jmaxo_
        do i=imino_,imaxo_
           temp=min(max(SC(i,j,k,isc_T),0.0_WP),5000.0_WP)
           do nsc=1,N_tot
              call volumeStats_finitechem_omegamu(temp*kOverEps(nsc),omegamu)
              eta(nsc)=mucoeff(nsc)*sqrt(temp)/omegamu
           end do
           call finitechem_W(SC(i,j,k,1:N_tot), W_mix)
           ! Use Wilke's method (same as Chemkin Manual 
           ! but obtained from Wilke's paper)
           do k_sc=1,N_tot
              do j_sc=1,N_tot
                 phi(k_sc,j_sc)=(1.0_WP/sqrt(8.0_WP))&
                      *((1.0_WP+(W_sp(k_sc)/W_sp(j_sc)))**-0.5_WP)
                 temp1=((eta(k_sc)/eta(j_sc))**0.5_WP)&
                      *((W_sp(j_sc)/W_sp(k_sc))**0.25_WP)
                 phi(k_sc,j_sc)=phi(k_sc,j_sc)*((1.0_WP+temp1)**2)
              end do
           end do
           sum_2=0.0_WP
           do k_sc=1,N_tot
              sum_1=0.0_WP
              do j_sc=1,N_tot
                 if (j_sc.ne.k_sc) then
                    sum_1=sum_1+(SC(i,j,k,j_sc)*W_mix*&
                         phi(k_sc,j_sc)/W_sp(j_sc))
                 end if
              end do
              sum_1=sum_1+SC(i,j,k,k_sc)*W_mix/W_sp(k_sc)
              if (SC(i,j,k,k_sc).gt.0.0_WP) then
                 sum_2=sum_2+(SC(i,j,k,k_sc)*W_mix*&
                      eta(k_sc)/(W_sp(k_sc)*sum_1))
              end if
           end do
           VISC(i,j,k)=sum_2
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine volumeStats_finitechem_viscosity


! ==================================================== !
! Compute the omegamu needed for transport calculation !
! Obtained from FlameMaster - Tspecies.c               !
! ==================================================== !
subroutine volumeStats_finitechem_omegamu(temp,omegamu)
  use volumeStats_finitechem
  implicit none

  real(WP) ::  m1 = 3.3530622607_WP
  real(WP) ::  m2 = 2.53272006_WP
  real(WP) ::  m3 = 2.9024238575_WP
  real(WP) ::  m4 = 0.11186138893_WP
  real(WP) ::  m5 = 0.8662326188_WP  ! = -0.1337673812 + 1.0
  real(WP) ::  m6 = 1.3913958626_WP
  real(WP) ::  m7 = 3.158490576_WP
  real(WP) ::  m8 = 0.18973411754_WP
  real(WP) ::  m9 = 0.00018682962894_WP
  real(WP) ::  num, den;
  real(WP) ::  temp, omegamu

  num = m1 + temp*(m2 + temp*(m3 + temp*m4))
  den = m5 + temp*(m6 + temp*(m7 + temp*(m8 + temp*m9)))
  omegamu = num/den

  return
end subroutine volumeStats_finitechem_omegamu
