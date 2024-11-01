!----------------------------------------------------------
! ======= HYDROGENF.f90 =======
!----------------------------------------------------------
subroutine PRODRATES( CDOT, W, K, C, M, TEMP, PRESSURE )
!----------------------------------------------------------
!     THIS SUBROUTINE COMPUTES RATES OF PRODUCTION CDOT
!     IN [KMOLE/(M^3S)]. THE PARAMETERS W ( REACTION RATE ),
!     K ( RATE COEFFICIENT ) AND M ( THIRD BODY CONCENTRATIONS ) ARE
!     JUST WORK SPACE FOR THIS FUNCTION.
!     C CONTAINS THE CONCENTRATIONS OF NON STEADY STATE SPECIES IN
!     [KMOLE/M^3] AND IS WORKSPACE FOR THE STEADY STATE 
!     CONCENTRATIONS, WHICH ARE COMPUTED IN THIS FUNCTION.
!     TEMP IS THE TEMPERATURE IN [K] AND
!     PRESSURE IS THE PRESSURE IN [PA].
!     CALLED FUNCTIONS ARE 'GETLINDRATECOEFF', 'COMPSTEADYSTATES',
!     'CTCHZERO'
!----------------------------------------------------------
      implicit none
      include 'HYDROGENF90.h'
      real(DP) :: CDOT(9), W(50), K(50), &
      C(9), M(5), TEMP, PRESSURE
      integer ::  I
      real(DP) :: GETLINDRATECOEFF, LT, RT
      real(DP), parameter ::  RGAS = 8314.34, CONCDEFAULT = -1.0 

      real(DP) ::  KINFTROE, K0TROE
      real(DP) ::  FCTROE

      LT = DLOG( TEMP )
      RT = RGAS * TEMP 


      M(MM1) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + C(SHO2) &
	    + C(SH2O2)
      M(MM2) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 6.3 * C(SH2O) + C(SHO2) &
	    + C(SH2O2)
      M(MM3) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2 * C(SH2) &
	    + 12 * C(SH2O) + C(SHO2) &
	    + C(SH2O2)
      M(MM4) = C(SN2) + C(SH) &
	    + C(SO2) + C(SO) &
	    + C(SOH) + 2.4 * C(SH2) &
	    + 15.4 * C(SH2O) + C(SHO2) &
	    + C(SH2O2)
      M(MM5) = C(SN2) + C(SH) &
	    + 0.85 * C(SO2) + C(SO) &
	    + C(SOH) + 0.75 * C(SH2) &
	    + 11.89 * C(SH2O) + C(SHO2) &
	    + C(SH2O2)


      K(R1F) = 2.6400000000D+13 &
	   * exp(-0.67 * LT - 71300000 / RT)
      K(R1B) = 5.2766713024D+10 &
	   * exp(-0.234291 * LT - 479604.2721 / RT)
      K(R2F) = 4.5900000000D+01 &
	   * exp(2.7 * LT - 26190000 / RT)
      K(R2B) = 2.7517349858D+01 &
	   * exp(2.66385 * LT - 20100582.44 / RT)
      K(R3F) = 1.7300000000D+05 &
	   * exp(1.51 * LT - 14350000 / RT)
      K(R3B) = 1.9032628245D+06 &
	   * exp(1.40625 * LT - 77394524.74 / RT)
      K(R4F) = 3.9700000000D+01 &
	   * exp(2.4 * LT + 8830000 / RT)
      K(R4B) = 7.2853303338D+02 &
	   * exp(2.33241 * LT - 60303942.3 / RT)
      K(R5F) = 1.7800000000D+12 * exp(-1 * LT)
      K(R5B) = 4.1604426150D+14 &
	   * exp(-0.649765 * LT - 433091487.8 / RT)
      K(R6F) = 9.0000000000D+10 * exp(-0.6 * LT)
      K(R6B) = 2.1035945806D+13 &
	   * exp(-0.249765 * LT - 433091487.8 / RT)
      K(R7F) = 5.6200000000D+13 * exp(-1.25 * LT)
      K(R7B) = 1.3135779492D+16 &
	   * exp(-0.899765 * LT - 433091487.8 / RT)
      K(R9F) = 4.4000000000D+16 * exp(-2 * LT)
      K(R9B) = 1.1314226588D+20 &
	   * exp(-1.75351 * LT - 496136012.5 / RT)
      K(R10F) = 9.4300000000D+12 * exp(-1 * LT)
      K(R10B) = 1.3213721422D+15 &
	   * exp(-0.685917 * LT - 427002070.2 / RT)
      K(R11F) = 1.2000000000D+11 * exp(-1 * LT)
      K(R11B) = 8.4127616577D+15 &
	   * exp(-1.12163 * LT - 497822466 / RT)
      K0TROE = 6.3300000000D+13 * exp(-1.4 * LT)
      KINFTROE = 5.1200000000D+09 * exp(0.44 * LT)
      FCTROE = 0.5 * EXP( -TEMP / 1e-30 ) &
	   + 0.5 * EXP( -TEMP / 1e+30 )
      K(R12F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM5) )
      K0TROE = 9.1184255084D+16 &
	   * exp(-1.41184 * LT - 204667815.3 / RT)
      KINFTROE = 7.3754089420D+12 &
	   * exp(0.428161 * LT - 204667815.3 / RT)
      K(R12B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM5) )
      K(R13F) = 5.9200000000D+02 &
	   * exp(2.43 * LT - 223850000 / RT)
      K(R13B) = 3.6485333254D+03 &
	   * exp(2.06793 * LT + 4573672.487 / RT)
      K0TROE = 2.0100000000D+11 &
	   * exp(-0.58 * LT + 9590000 / RT)
      KINFTROE = 1.1100000000D+11 * exp(-0.37 * LT)
      FCTROE = 0.2654 * EXP( -TEMP / 94 ) &
	   + 0.7346 * EXP( -TEMP / 1756 ) &
	   + 1 * EXP( -5182 / TEMP )
      K(R14F) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K0TROE = 7.1253840510D+19 &
	   * exp(-1.59836 * LT - 205148485.1 / RT)
      KINFTROE = 3.9349135804D+19 &
	   * exp(-1.38836 * LT - 214738485.1 / RT)
      K(R14B) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &
	   , FCTROE, M(mM3) )
      K(R15F) = 3.9700000000D+09 * exp(-2810000 / RT)
      K(R15B) = 1.4164546439D+07 &
	   * exp(0.694036 * LT - 223457801.5 / RT)
      K(R16F) = 7.4900000000D+10 * exp(-2660000 / RT)
      K(R16B) = 1.4562476643D+07 &
	   * exp(0.761631 * LT - 154173859.2 / RT)
      K(R17F) = 4.0000000000D+10
      K(R17B) = 3.8909647883D+09 &
	   * exp(0.325921 * LT - 222334254.9 / RT)
      K(R18F) = 2.3800000000D+10 * exp(2090000 / RT)
      K(R18B) = 4.2484744234D+10 &
	   * exp(0.258327 * LT - 289378197.2 / RT)
      K(R19F) = 1.0000000000D+13 * exp(-72510000 / RT)
      K(R19B) = 1.7850732872D+13 &
	   * exp(0.258327 * LT - 363978197.2 / RT)
      K(R20F) = 1.3000000000D+08 * exp(6820000 / RT)
      K(R20B) = 6.2200353343D+09 &
	   * exp(-0.244887 * LT - 154764529 / RT)
      K(R21F) = 3.6600000000D+11 * exp(-50210000 / RT)
      K(R21B) = 1.7511791788D+13 &
	   * exp(-0.244887 * LT - 211794529 / RT)
      K(R22F) = 6.0500000000D+03 &
	   * exp(2 * LT - 21760000 / RT)
      K(R22B) = 2.0516782971D+01 &
	   * exp(2.60696 * LT - 88599143.48 / RT)
      K(R23F) = 2.4100000000D+10 * exp(-16610000 / RT)
      K(R23B) = 1.7481432523D+05 &
	   * exp(1.26484 * LT - 298007527.4 / RT)
      K(R24F) = 9.6300000000D+03 &
	   * exp(2 * LT - 16610000 / RT)
      K(R24B) = 1.9578260237D+01 &
	   * exp(2.57081 * LT - 77359725.92 / RT)
      K(R25F) = 2.0000000000D+09 * exp(-1790000 / RT)
      K(R25B) = 7.4616787480D+07 &
	   * exp(0.503214 * LT - 131673668.2 / RT)
      K(R26F) = 2.6700000000D+38 &
	   * exp(-7 * LT - 157320000 / RT)
      K(R26B) = 9.9613411288D+36 &
	   * exp(-6.49679 * LT - 287203668.2 / RT)

      W(R1F) = K(R1F) * C(SH) * C(SO2)
      W(R1B) = K(R1B) * C(SOH) * C(SO)
      W(R2F) = K(R2F) * C(SO) * C(SH2)
      W(R2B) = K(R2B) * C(SOH) * C(SH)
      W(R3F) = K(R3F) * C(SOH) * C(SH2)
      W(R3B) = K(R3B) * C(SH2O) * C(SH)
      W(R4F) = K(R4F) * C(SOH) * C(SOH)
      W(R4B) = K(R4B) * C(SH2O) * C(SO)
      W(R5F) = K(R5F) * C(SH) * C(SH) * M(MM1)
      W(R5B) = K(R5B) * C(SH2) * M(MM1)
      W(R6F) = K(R6F) * C(SH) * C(SH) * C(SH2)
      W(R6B) = K(R6B) * C(SH2) * C(SH2)
      W(R7F) = K(R7F) * C(SH) * C(SH) * C(SH2O)
      W(R7B) = K(R7B) * C(SH2O) * C(SH2)
      W(R9F) = K(R9F) * C(SH) * C(SOH) * M(MM2)
      W(R9B) = K(R9B) * C(SH2O) * M(MM2)
      W(R10F) = K(R10F) * C(SH) * C(SO) * M(MM3)
      W(R10B) = K(R10B) * C(SOH) * M(MM3)
      W(R11F) = K(R11F) * C(SO) * C(SO) * M(MM4)
      W(R11B) = K(R11B) * C(SO2) * M(MM4)
      W(R12F) = K(R12F) * C(SH) * C(SO2)
      W(R12B) = K(R12B) * C(SHO2)
      W(R13F) = K(R13F) * C(SH2) * C(SO2)
      W(R13B) = K(R13B) * C(SH) * C(SHO2)
      W(R14F) = K(R14F) * C(SOH) * C(SOH)
      W(R14B) = K(R14B) * C(SH2O2)
      W(R15F) = K(R15F) * C(SHO2) * C(SH)
      W(R15B) = K(R15B) * C(SH2O) * C(SO)
      W(R16F) = K(R16F) * C(SHO2) * C(SH)
      W(R16B) = K(R16B) * C(SOH) * C(SOH)
      W(R17F) = K(R17F) * C(SHO2) * C(SO)
      W(R17B) = K(R17B) * C(SO2) * C(SOH)
      W(R18F) = K(R18F) * C(SHO2) * C(SOH)
      W(R18B) = K(R18B) * C(SO2) * C(SH2O)
      W(R19F) = K(R19F) * C(SHO2) * C(SOH)
      W(R19B) = K(R19B) * C(SO2) * C(SH2O)
      W(R20F) = K(R20F) * C(SHO2) * C(SHO2)
      W(R20B) = K(R20B) * C(SH2O2) * C(SO2)
      W(R21F) = K(R21F) * C(SHO2) * C(SHO2)
      W(R21B) = K(R21B) * C(SH2O2) * C(SO2)
      W(R22F) = K(R22F) * C(SH2O2) * C(SH)
      W(R22B) = K(R22B) * C(SH2) * C(SHO2)
      W(R23F) = K(R23F) * C(SH2O2) * C(SH)
      W(R23B) = K(R23B) * C(SOH) * C(SH2O)
      W(R24F) = K(R24F) * C(SH2O2) * C(SO)
      W(R24B) = K(R24B) * C(SOH) * C(SHO2)
      W(R25F) = K(R25F) * C(SH2O2) * C(SOH)
      W(R25B) = K(R25B) * C(SH2O) * C(SHO2)
      W(R26F) = K(R26F) * C(SH2O2) * C(SOH)
      W(R26B) = K(R26B) * C(SH2O) * C(SHO2)


      CDOT(SN2) = 0.0_DP
      CDOT(SH) = - W(R1F) + W(R1B) + W(R2F) & 
	    - W(R2B) + W(R3F) - W(R3B) & 
	    - 2 * W(R5F) + 2 * W(R5B) - 2 * W(R6F) & 
	    + 2 * W(R6B) - 2 * W(R7F) + 2 * W(R7B) & 
	    - W(R9F) + W(R9B) - W(R10F) & 
	    + W(R10B) - W(R12F) + W(R12B) & 
	    + W(R13F) - W(R13B) - W(R15F) & 
	    + W(R15B) - W(R16F) + W(R16B) & 
	    - W(R22F) + W(R22B) - W(R23F) & 
	    + W(R23B)
      CDOT(SO2) = - W(R1F) + W(R1B) + W(R11F) & 
	    - W(R11B) - W(R12F) + W(R12B) & 
	    - W(R13F) + W(R13B) + W(R17F) & 
	    - W(R17B) + W(R18F) - W(R18B) & 
	    + W(R19F) - W(R19B) + W(R20F) & 
	    - W(R20B) + W(R21F) - W(R21B)
      CDOT(SO) = W(R1F) - W(R1B) - W(R2F) & 
	    + W(R2B) + W(R4F) - W(R4B) & 
	    - W(R10F) + W(R10B) - 2 * W(R11F) & 
	    + 2 * W(R11B) + W(R15F) - W(R15B) & 
	    - W(R17F) + W(R17B) - W(R24F) & 
	    + W(R24B)
      CDOT(SOH) = W(R1F) - W(R1B) + W(R2F) & 
	    - W(R2B) - W(R3F) + W(R3B) & 
	    - 2 * W(R4F) + 2 * W(R4B) - W(R9F) & 
	    + W(R9B) + W(R10F) - W(R10B) & 
	    - 2 * W(R14F) + 2 * W(R14B) + 2 * W(R16F) & 
	    - 2 * W(R16B) + W(R17F) - W(R17B) & 
	    - W(R18F) + W(R18B) - W(R19F) & 
	    + W(R19B) + W(R23F) - W(R23B) & 
	    + W(R24F) - W(R24B) - W(R25F) & 
	    + W(R25B) - W(R26F) + W(R26B)
      CDOT(SH2) = - W(R2F) + W(R2B) - W(R3F) & 
	    + W(R3B) + W(R5F) - W(R5B) & 
	    - W(R6F) + 2 * W(R6F) - 2 * W(R6B) & 
	    + W(R6B) + W(R7F) - W(R7B) & 
	    - W(R13F) + W(R13B) + W(R22F) & 
	    - W(R22B)
      CDOT(SH2O) = W(R3F) - W(R3B) + W(R4F) & 
	    - W(R4B) - W(R7F) + W(R7F) & 
	    - W(R7B) + W(R7B) + W(R9F) & 
	    - W(R9B) + W(R15F) - W(R15B) & 
	    + W(R18F) - W(R18B) + W(R19F) & 
	    - W(R19B) + W(R23F) - W(R23B) & 
	    + W(R25F) - W(R25B) + W(R26F) & 
	    - W(R26B)
      CDOT(SHO2) = W(R12F) - W(R12B) + W(R13F) & 
	    - W(R13B) - W(R15F) + W(R15B) & 
	    - W(R16F) + W(R16B) - W(R17F) & 
	    + W(R17B) - W(R18F) + W(R18B) & 
	    - W(R19F) + W(R19B) - 2 * W(R20F) & 
	    + 2 * W(R20B) - 2 * W(R21F) + 2 * W(R21B) & 
	    + W(R22F) - W(R22B) + W(R24F) & 
	    - W(R24B) + W(R25F) - W(R25B) & 
	    + W(R26F) - W(R26B)
      CDOT(SH2O2) = W(R14F) - W(R14B) + W(R20F) & 
	    - W(R20B) + W(R21F) - W(R21B) & 
	    - W(R22F) + W(R22B) - W(R23F) & 
	    + W(R23B) - W(R24F) + W(R24B) & 
	    - W(R25F) + W(R25B) - W(R26F) & 
	    + W(R26B)
      END


subroutine CONSRATES( CDOT_m, C, TEMP, PRESSURE )
!----------------------------------------------------------
!     THIS SUBROUTINE COMPUTES RATES OF CONSUMPTION CDOT_m
!     IN [KMOLE/(M^3S)]. THE PARAMETERS W ( REACTION RATE ),
!     K ( RATE COEFFICIENT ) AND M ( THIRD BODY CONCENTRATIONS ) ARE
!     COMPUTED IN A CALL TO 'PRODRATES'.
!----------------------------------------------------------
      implicit none
      include 'HYDROGENF90.h'
      real(DP) :: CDOT(9), W(50), K(50), &
      C(9), M(5), TEMP, PRESSURE, CDOT_m(9)


      call PRODRATES( CDOT, W, K, C, M, TEMP, PRESSURE )


      CDOT_m(SN2) = 0.0_DP
      CDOT_m(SH) = - W(R1F) - W(R2B) - W(R3B) & 
	    - 2 * W(R5F) - 2 * W(R6F) - 2 * W(R7F) & 
	    - W(R9F) - W(R10F) - W(R12F) & 
	    - W(R13B) - W(R15F) - W(R16F) & 
	    - W(R22F) - W(R23F)
      CDOT_m(SO2) = - W(R1F) - W(R11B) - W(R12F) & 
	    - W(R13F) - W(R17B) - W(R18B) & 
	    - W(R19B) - W(R20B) - W(R21B)
      CDOT_m(SO) = - W(R1B) - W(R2F) - W(R4B) & 
	    - W(R10F) - 2 * W(R11F) - W(R15B) & 
	    - W(R17F) - W(R24F)
      CDOT_m(SOH) = - W(R1B) - W(R2B) - W(R3F) & 
	    - 2 * W(R4F) - W(R9F) - W(R10B) & 
	    - 2 * W(R14F) - 2 * W(R16B) - W(R17B) & 
	    - W(R18F) - W(R19F) - W(R23B) & 
	    - W(R24B) - W(R25F) - W(R26F)
      CDOT_m(SH2) = - W(R2F) - W(R3F) - W(R5B) & 
	    - W(R6F) - 2 * W(R6B) - W(R7B) & 
	    - W(R13F) - W(R22B)
      CDOT_m(SH2O) = - W(R3B) - W(R4B) - W(R7F) & 
	    - W(R7B) - W(R9B) - W(R15B) & 
	    - W(R18B) - W(R19B) - W(R23B) & 
	    - W(R25B) - W(R26B)
      CDOT_m(SHO2) = - W(R12B) - W(R13B) - W(R15F) & 
	    - W(R16F) - W(R17F) - W(R18F) & 
	    - W(R19F) - 2 * W(R20F) - 2 * W(R21F) & 
	    - W(R22B) - W(R24B) - W(R25B) & 
	    - W(R26B)
      CDOT_m(SH2O2) = - W(R14B) - W(R20B) - W(R21B) & 
	    - W(R22F) - W(R23F) - W(R24F) & 
	    - W(R25F) - W(R26F)
      END


      DOUBLE PRECISION FUNCTION GETLINDRATECOEFF( TEMP, PRESSURE, &
	   K0, KINF, FC, CONCIN )

      IMPLICIT NONE
      INTEGER ITROE
      DOUBLE PRECISION TEMP, PRESSURE, K0, KINF, FC
      DOUBLE PRECISION R
      PARAMETER (R = 8314.34)
!     [J / kmole K]
      DOUBLE PRECISION NTMP,CCOEFF,DCOEFF,LGKNULL
      DOUBLE PRECISION KL
      DOUBLE PRECISION F
      DOUBLE PRECISION CONC, CONCIN

      ITROE = 1
      IF (CONCIN.GT.0.0) THEN
        CONC = CONCIN
      ELSE
        CONC = PRESSURE / ( R * TEMP )
      END IF
      NTMP = 0.75 - 1.27 * DLOG10( FC )
      IF (ITROE.EQ.1) THEN
        CCOEFF = - 0.4 - 0.67 * DLOG10( FC )
        DCOEFF = 0.14
        K0 = K0 * CONC / MAX(KINF, 1.0D-60)
        LGKNULL = DLOG10( K0 )
        F=(LGKNULL+CCOEFF)/(NTMP-DCOEFF*(LGKNULL+CCOEFF))
        F = FC**(1.0 / ( F * F + 1.0 ))
        GETLINDRATECOEFF = KINF * F * K0 / ( 1.0 + K0 )
      ELSE
        K0 = K0 * CONC / KINF
        KL = K0 / ( 1.0 + K0 )
        F = DLOG10( K0 ) / NTMP
        F = FC ** ( 1.0 / ( F * F + 1.0 ) )
        GETLINDRATECOEFF = KINF * F * KL
      END IF
      END


      DOUBLE PRECISION FUNCTION CTCHZERO( X )

      IMPLICIT NONE
      DOUBLE PRECISION X
      IF ( X.EQ.0.0D0 ) THEN
         X = 1.0D-20
      END IF
      CTCHZERO = X
      END


      DOUBLE PRECISION FUNCTION SOLVQDRT( A, B, C )
!------------------------------------------------------------------
!	SOLVES THE QUADRATIC EQUATION A x^2 + B x - C = 0
!------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION RAD, A, B, C

      IF ( A.EQ.0.0D0 ) THEN
        IF ( B.EQ.0.0D0 ) THEN
          WRITE(*,*) '#WARNING: INVALID ARGUMENTS IN FUNCTION SOLVQDRT '
          SOLVQDRT = 0.0
        ELSE
          SOLVQDRT = C / B
        END IF
      ELSE

        B = B / A
        C = C / A

        RAD = 0.25D0 * B * B + C
        IF ( RAD.GE.0.0D0 ) THEN
          RAD = DSQRT( RAD ) + 0.5 * B
          IF ( RAD.NE.0.0D0 ) THEN
            SOLVQDRT = C / RAD
          ELSE
            SOLVQDRT = -0.5D0 * B + DSQRT( 0.25 * B * B + C )
          ENDIF
        ELSE
          SOLVQDRT = 0.0D0
        ENDIF
      END IF
      END


      SUBROUTINE GETMOLARMASS( MM )
!------------------------------------------------------------------
!	FILLS 'MM' WITH SPECIES MOLAR MASS IN KG/KMOLE
!------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION MM(9)
      INCLUDE 'HYDROGENF90.h'

      MM(SN2) =  2.80200000e+01
      MM(SH) =  1.00800000e+00
      MM(SO2) =  3.20000000e+01
      MM(SO) =  1.60000000e+01
      MM(SOH) =  1.70080000e+01
      MM(SH2) =  2.01600000e+00
      MM(SH2O) =  1.80160000e+01
      MM(SHO2) =  3.30080000e+01
      MM(SH2O2) =  3.40160000e+01

      END


      SUBROUTINE GETSPECIESNAMES( NAMES )
!------------------------------------------------------------------
!	FILLS 'NAMES' WITH SPECIES IDENTIFIER/KG
!------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER *20 NAMES(9)
      INCLUDE 'HYDROGENF90.h'

      NAMES(SN2)='N2                  '
      NAMES(SH)='H                   '
      NAMES(SO2)='O2                  '
      NAMES(SO)='O                   '
      NAMES(SOH)='OH                  '
      NAMES(SH2)='H2                  '
      NAMES(SH2O)='H2O                 '
      NAMES(SHO2)='HO2                 '
      NAMES(SH2O2)='H2O2                '

      END


      SUBROUTINE GETNSPECIES( NSPECIES )
!------------------------------------------------------------------
!	FILLS 'NSPECIES' WITH NUMBER OF SPECIES 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NSPECIES
      INCLUDE 'HYDROGENF90.h'

      NSPECIES = SEND - 1

      END


      SUBROUTINE GETNREACTIONS( NREACTIONS )
!------------------------------------------------------------------
!	FILLS 'NREACTIONS' WITH NUMBER OF REACTIONS 
!------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NREACTIONS

      NREACTIONS = 50

      END


      SUBROUTINE GETNSPECS(NSPECIES_NONS)
!------------------------------------------------------------------
!     RETURNS THE NUMBER OF THE NON STEADY STATE SPECIES
!------------------------------------------------------------------
      implicit none
      integer ::  NSPECIES_NONS
      include 'HYDROGENF90.h'

      NSPECIES_NONS = 9
      END


      SUBROUTINE GETMUCOEFF( MUCOEFF )
!------------------------------------------------------------------
!	FILLS 'MUCOEFF' WITH MUECOEFF IN KG/(M*S)
!------------------------------------------------------------------

      implicit none

      include 'HYDROGENF90.h'
      real(DP) :: MUCOEFF(9)

      MUCOEFF(SN2) =  1.07764173e-06
      MUCOEFF(SH) =  6.37705159e-07
      MUCOEFF(SO2) =  1.26276460e-06
      MUCOEFF(SO) =  1.41186116e-06
      MUCOEFF(SOH) =  1.45565556e-06
      MUCOEFF(SH2) =  4.44505304e-07
      MUCOEFF(SH2O) =  1.66959493e-06
      MUCOEFF(SHO2) =  1.28249894e-06
      MUCOEFF(SH2O2) =  1.30193418e-06

      END


      SUBROUTINE GETKOVEREPS( KOVEREPS )
!------------------------------------------------------------------
!	    FILLS 'KOVEREPS' WITH KOVEREPS
!------------------------------------------------------------------

      implicit none

      include 'HYDROGENF90.h'
      real(DP) :: KOVEREPS(9)

      KOVEREPS(SN2) =  1.02532554e-02
      KOVEREPS(SH) =  6.89655172e-03
      KOVEREPS(SO2) =  9.31098696e-03
      KOVEREPS(SO) =  1.25000000e-02
      KOVEREPS(SOH) =  1.25000000e-02
      KOVEREPS(SH2) =  2.63157895e-02
      KOVEREPS(SH2O) =  1.74703005e-03
      KOVEREPS(SHO2) =  9.31098696e-03
      KOVEREPS(SH2O2) =  9.31098696e-03

      END


      SUBROUTINE COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )
!------------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS FROM
!     THE CONCENTRATIONS OF COMPUTED SPECIES AND RATE COEFFICIENTS.
!     CONCENTRATIONS OF COMPUTED SPECIES MAY NOT BE ALTERED.
!------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'HYDROGENF90.h'
      integer  ::  NSPECIN, NSPEC, INOW
      real(DP) :: K(50), C(9), M(5)
      real(DP) :: TEMP, PRESSURE, CTOT
      real(DP), parameter ::  R= 8314.34
      END


      SUBROUTINE COMPTHERMODATA( H, CP, T )
!------------------------------------------------------------------
!     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS
!     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES
!     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.
!     THE PARAMETER H AND CP SHOULD PROVIDE WORKSPACE OF LENGTH 9
!------------------------------------------------------------------
      implicit none
      include 'HYDROGENF90.h'
      real(DP) :: H(9), CP(9), T

      IF (T.GT.1000.0_DP) THEN
      H(SN2) =  2.96728765e+02_DP * ( &
	   T * (  2.92664000e+00_DP + T * (  7.43988400e-04_DP &
	   + T * ( -1.89492000e-07_DP + T * (  2.52425950e-11_DP &
	   + T * ( -1.35067020e-15_DP ) ) ) ) ) -9.22797700e+02_DP )
      CP(SN2) =  2.96728765e+02_DP * ( &
	    2.92664000e+00_DP + T * (  1.48797680e-03_DP &
	   + T * ( -5.68476000e-07_DP + T * (  1.00970380e-10_DP &
	   + T * ( -6.75335100e-15_DP ) ) ) ) )
      H(SH) =  8.24835317e+03_DP * ( &
	   T * (  2.50000001e+00_DP + T * ( -1.15421486e-11_DP &
	   + T * (  5.38539827e-15_DP + T * ( -1.18378809e-18_DP &
	   + T * (  9.96394714e-23_DP ) ) ) ) ) +  2.54736599e+04_DP )
      CP(SH) =  8.24835317e+03_DP * ( &
	    2.50000001e+00_DP + T * ( -2.30842973e-11_DP &
	   + T * (  1.61561948e-14_DP + T * ( -4.73515235e-18_DP &
	   + T * (  4.98197357e-22_DP ) ) ) ) )
      H(SO2) =  2.59823125e+02_DP * ( &
	   T * (  3.28253784e+00_DP + T * (  7.41543770e-04_DP &
	   + T * ( -2.52655556e-07_DP + T * (  5.23676387e-11_DP &
	   + T * ( -4.33435588e-15_DP ) ) ) ) ) -1.08845772e+03_DP )
      CP(SO2) =  2.59823125e+02_DP * ( &
	    3.28253784e+00_DP + T * (  1.48308754e-03_DP &
	   + T * ( -7.57966669e-07_DP + T * (  2.09470555e-10_DP &
	   + T * ( -2.16717794e-14_DP ) ) ) ) )
      H(SO) =  5.19646250e+02_DP * ( &
	   T * (  2.56942078e+00_DP + T * ( -4.29870569e-05_DP &
	   + T * (  1.39828196e-08_DP + T * ( -2.50444497e-12_DP &
	   + T * (  2.45667382e-16_DP ) ) ) ) ) +  2.92175791e+04_DP )
      CP(SO) =  5.19646250e+02_DP * ( &
	    2.56942078e+00_DP + T * ( -8.59741137e-05_DP &
	   + T * (  4.19484589e-08_DP + T * ( -1.00177799e-11_DP &
	   + T * (  1.22833691e-15_DP ) ) ) ) )
      H(SOH) =  4.88848777e+02_DP * ( &
	   T * (  2.86472886e+00_DP + T * (  5.28252240e-04_DP &
	   + T * ( -8.63609193e-08_DP + T * (  7.63046685e-12_DP &
	   + T * ( -2.66391752e-16_DP ) ) ) ) ) +  3.71885774e+03_DP )
      CP(SOH) =  4.88848777e+02_DP * ( &
	    2.86472886e+00_DP + T * (  1.05650448e-03_DP &
	   + T * ( -2.59082758e-07_DP + T * (  3.05218674e-11_DP &
	   + T * ( -1.33195876e-15_DP ) ) ) ) )
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  3.33727920e+00_DP + T * ( -2.47012365e-05_DP &
	   + T * (  1.66485593e-07_DP + T * ( -4.48915985e-11_DP &
	   + T * (  4.00510752e-15_DP ) ) ) ) ) -9.50158922e+02_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    3.33727920e+00_DP + T * ( -4.94024731e-05_DP &
	   + T * (  4.99456778e-07_DP + T * ( -1.79566394e-10_DP &
	   + T * (  2.00255376e-14_DP ) ) ) ) )
      H(SH2O) =  4.61497558e+02_DP * ( &
	   T * (  3.03399249e+00_DP + T * (  1.08845902e-03_DP &
	   + T * ( -5.46908393e-08_DP + T * ( -2.42604967e-11_DP &
	   + T * (  3.36401984e-15_DP ) ) ) ) ) -3.00042971e+04_DP )
      CP(SH2O) =  4.61497558e+02_DP * ( &
	    3.03399249e+00_DP + T * (  2.17691804e-03_DP &
	   + T * ( -1.64072518e-07_DP + T * ( -9.70419870e-11_DP &
	   + T * (  1.68200992e-14_DP ) ) ) ) )
      H(SHO2) =  2.51888633e+02_DP * ( &
	   T * (  4.01721090e+00_DP + T * (  1.11991006e-03_DP &
	   + T * ( -2.11219383e-07_DP + T * (  2.85615925e-11_DP &
	   + T * ( -2.15817070e-15_DP ) ) ) ) ) +  1.11856713e+02_DP )
      CP(SHO2) =  2.51888633e+02_DP * ( &
	    4.01721090e+00_DP + T * (  2.23982013e-03_DP &
	   + T * ( -6.33658150e-07_DP + T * (  1.14246370e-10_DP &
	   + T * ( -1.07908535e-14_DP ) ) ) ) )
      H(SH2O2) =  2.44424389e+02_DP * ( &
	   T * (  4.16500285e+00_DP + T * (  2.45415847e-03_DP &
	   + T * ( -6.33797417e-07_DP + T * (  9.27964965e-11_DP &
	   + T * ( -5.75816610e-15_DP ) ) ) ) ) -1.78617877e+04_DP )
      CP(SH2O2) =  2.44424389e+02_DP * ( &
	    4.16500285e+00_DP + T * (  4.90831694e-03_DP &
	   + T * ( -1.90139225e-06_DP + T * (  3.71185986e-10_DP &
	   + T * ( -2.87908305e-14_DP ) ) ) ) )
      ELSE
      H(SN2) =  2.96728765e+02_DP * ( &
	   T * (  3.29867700e+00_DP + T * (  7.04120200e-04_DP &
	   + T * ( -1.32107400e-06_DP + T * (  1.41037875e-09_DP &
	   + T * ( -4.88970800e-13_DP ) ) ) ) ) -1.02089990e+03_DP )
      CP(SN2) =  2.96728765e+02_DP * ( &
	    3.29867700e+00_DP + T * (  1.40824040e-03_DP &
	   + T * ( -3.96322200e-06_DP + T * (  5.64151500e-09_DP &
	   + T * ( -2.44485400e-12_DP ) ) ) ) )
      H(SH) =  8.24835317e+03_DP * ( &
	   T * (  2.50000000e+00_DP + T * (  3.52666409e-13_DP &
	   + T * ( -6.65306547e-16_DP + T * (  5.75204080e-19_DP &
	   + T * ( -1.85546466e-22_DP ) ) ) ) ) +  2.54736599e+04_DP )
      CP(SH) =  8.24835317e+03_DP * ( &
	    2.50000000e+00_DP + T * (  7.05332819e-13_DP &
	   + T * ( -1.99591964e-15_DP + T * (  2.30081632e-18_DP &
	   + T * ( -9.27732332e-22_DP ) ) ) ) )
      H(SO2) =  2.59823125e+02_DP * ( &
	   T * (  3.78245636e+00_DP + T * ( -1.49836708e-03_DP &
	   + T * (  3.28243400e-06_DP + T * ( -2.42032377e-09_DP &
	   + T * (  6.48745674e-13_DP ) ) ) ) ) -1.06394356e+03_DP )
      CP(SO2) =  2.59823125e+02_DP * ( &
	    3.78245636e+00_DP + T * ( -2.99673416e-03_DP &
	   + T * (  9.84730201e-06_DP + T * ( -9.68129509e-09_DP &
	   + T * (  3.24372837e-12_DP ) ) ) ) )
      H(SO) =  5.19646250e+02_DP * ( &
	   T * (  3.16826710e+00_DP + T * ( -1.63965942e-03_DP &
	   + T * (  2.21435465e-06_DP + T * ( -1.53201656e-09_DP &
	   + T * (  4.22531942e-13_DP ) ) ) ) ) +  2.91222592e+04_DP )
      CP(SO) =  5.19646250e+02_DP * ( &
	    3.16826710e+00_DP + T * ( -3.27931884e-03_DP &
	   + T * (  6.64306396e-06_DP + T * ( -6.12806624e-09_DP &
	   + T * (  2.11265971e-12_DP ) ) ) ) )
      H(SOH) =  4.88848777e+02_DP * ( &
	   T * (  4.12530561e+00_DP + T * ( -1.61272470e-03_DP &
	   + T * (  2.17588230e-06_DP + T * ( -1.44963411e-09_DP &
	   + T * (  4.12474758e-13_DP ) ) ) ) ) +  3.38153812e+03_DP )
      CP(SOH) =  4.88848777e+02_DP * ( &
	    4.12530561e+00_DP + T * ( -3.22544939e-03_DP &
	   + T * (  6.52764691e-06_DP + T * ( -5.79853643e-09_DP &
	   + T * (  2.06237379e-12_DP ) ) ) ) )
      H(SH2) =  4.12417659e+03_DP * ( &
	   T * (  2.34433112e+00_DP + T * (  3.99026037e-03_DP &
	   + T * ( -6.49271700e-06_DP + T * (  5.03930235e-09_DP &
	   + T * ( -1.47522352e-12_DP ) ) ) ) ) -9.17935173e+02_DP )
      CP(SH2) =  4.12417659e+03_DP * ( &
	    2.34433112e+00_DP + T * (  7.98052075e-03_DP &
	   + T * ( -1.94781510e-05_DP + T * (  2.01572094e-08_DP &
	   + T * ( -7.37611761e-12_DP ) ) ) ) )
      H(SH2O) =  4.61497558e+02_DP * ( &
	   T * (  4.19864056e+00_DP + T * ( -1.01821705e-03_DP &
	   + T * (  2.17346737e-06_DP + T * ( -1.37199266e-09_DP &
	   + T * (  3.54395634e-13_DP ) ) ) ) ) -3.02937267e+04_DP )
      CP(SH2O) =  4.61497558e+02_DP * ( &
	    4.19864056e+00_DP + T * ( -2.03643410e-03_DP &
	   + T * (  6.52040211e-06_DP + T * ( -5.48797062e-09_DP &
	   + T * (  1.77197817e-12_DP ) ) ) ) )
      H(SHO2) =  2.51888633e+02_DP * ( &
	   T * (  4.30179801e+00_DP + T * ( -2.37456025e-03_DP &
	   + T * (  7.05276303e-06_DP + T * ( -6.06909735e-09_DP &
	   + T * (  1.85845025e-12_DP ) ) ) ) ) +  2.94808040e+02_DP )
      CP(SHO2) =  2.51888633e+02_DP * ( &
	    4.30179801e+00_DP + T * ( -4.74912051e-03_DP &
	   + T * (  2.11582891e-05_DP + T * ( -2.42763894e-08_DP &
	   + T * (  9.29225124e-12_DP ) ) ) ) )
      H(SH2O2) =  2.44424389e+02_DP * ( &
	   T * (  4.27611269e+00_DP + T * ( -2.71411208e-04_DP &
	   + T * (  5.57785670e-06_DP + T * ( -5.39427032e-09_DP &
	   + T * (  1.72490873e-12_DP ) ) ) ) ) -1.77025821e+04_DP )
      CP(SH2O2) =  2.44424389e+02_DP * ( &
	    4.27611269e+00_DP + T * ( -5.42822417e-04_DP &
	   + T * (  1.67335701e-05_DP + T * ( -2.15770813e-08_DP &
	   + T * (  8.62454363e-12_DP ) ) ) ) )
      END IF

      END


      SUBROUTINE COMPCHEMJACOBIAN( JAC, Y, T, P )
!------------------------------------------------------------------
!     THIS FUNCTION COMPUTES THE CHEMICAL JACOBIAN dM/dY AS A
!     FUNKTION OF CONCENTRATIONS 'C' OF ALL THE NON STEADY STATE
!     SPECIES IN [KMOLE/M^3], WITH DIMENSIONS (N_NONS+1,N_NONS+1).
!     THE PARAMETERS H AND CP SHOULD PROVIDE WORKSPACE OF LENGTH 9.
!------------------------------------------------------------------
      implicit none
      include 'HYDROGENF90.h'
      real(DP) :: JAC(10,10), Y(9), H(9), CP(9), T, T2
      real(DP) :: CDOT(9), W(50), K(50), K2(50), &
         C(9), M(5), P, MM(9), MMix, tmpv(9), rho_gp, Cp_mix
      real(DP), parameter :: RGAS = 8314.34
      integer :: i
      integer, parameter :: ST = 10

!     Compute  mixture molar mass and density
      call GETMOLARMASS( MM )
      tmpv = min(max(Y,0.0),1.0)
      MMix = sum(tmpv)/sum(tmpv/MM)
      rho_gp = P * MMix / ( RGAS * T )

!     Compute species concentrations from mass fractions
      do i=1,9
         C(i) = rho_gp * Y(i) / MM(i)
      end do

!     Compute mixture thermo data
      call COMPTHERMODATA( H, CP, T )
      Cp_mix = sum(Y*CP)

!     Get rate coefficients and reaction rates
      T2 = 1.00001 * T
      call PRODRATES( CDOT, W, K2, C, M, T2, P )
      call PRODRATES( CDOT, W, K, C, M, T, P )

      JAC(SN2,SN2) = MM(SN2) * ( 0.0 )
      JAC(SN2,SH) = MM(SN2) * ( 0.0 )
      JAC(SN2,SO2) = MM(SN2) * ( 0.0 )
      JAC(SN2,SO) = MM(SN2) * ( 0.0 )
      JAC(SN2,SOH) = MM(SN2) * ( 0.0 )
      JAC(SN2,SH2) = MM(SN2) * ( 0.0 )
      JAC(SN2,SH2O) = MM(SN2) * ( 0.0 )
      JAC(SN2,SHO2) = MM(SN2) * ( 0.0 )
      JAC(SN2,SH2O2) = MM(SN2) * ( 0.0 )
      JAC(SN2,ST) = MM(SN2) * ( 0.0 )

      JAC(SH,SN2) = MM(SH) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R3F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 2 * W(R5F) * ( - 3 * MMix/MM(SN2) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SN2) ) & 
	         + 2 * W(R5B) * ( - 2 * MMix/MM(SN2) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SN2) ) & 
	         - 2 * W(R6F) * ( - 3 * MMix/MM(SN2) ) & 
	         + 2 * W(R6B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 2 * W(R7F) * ( - 3 * MMix/MM(SN2) ) & 
	         + 2 * W(R7B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SN2) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SN2) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SN2) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SN2) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SN2) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SN2) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SN2) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SN2) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SN2) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SN2) ) )
      JAC(SH,SH) = MM(SH) * ( - 1 * W(R1F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R2B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R3F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R3B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 2 * W(R5F) * ( 2 / Y(SH) - 3 * MMix/MM(SH) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SH) ) & 
	         + 2 * W(R5B) * ( - 2 * MMix/MM(SH) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SH) ) & 
	         - 2 * W(R6F) * ( 2 / Y(SH) - 3 * MMix/MM(SH) ) & 
	         + 2 * W(R6B) * ( - 2 * MMix/MM(SH) ) & 
	         - 2 * W(R7F) * ( 2 / Y(SH) - 3 * MMix/MM(SH) ) & 
	         + 2 * W(R7B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R9F) * ( 1 / Y(SH) - 3 * MMix/MM(SH) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SH) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SH) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SH) ) & 
	         - 1 * W(R10F) * ( 1 / Y(SH) - 3 * MMix/MM(SH) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SH) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SH) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SH) ) & 
	         - 1 * W(R12F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SH) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R13B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R15F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R16F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R22F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R23F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SH) ) )
      JAC(SH,SO2) = MM(SH) * ( - 1 * W(R1F) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R3F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 2 * W(R5F) * ( - 3 * MMix/MM(SO2) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SO2) ) & 
	         + 2 * W(R5B) * ( - 2 * MMix/MM(SO2) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SO2) ) & 
	         - 2 * W(R6F) * ( - 3 * MMix/MM(SO2) ) & 
	         + 2 * W(R6B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 2 * W(R7F) * ( - 3 * MMix/MM(SO2) ) & 
	         + 2 * W(R7B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SO2) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SO2) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SO2) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SO2) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SO2) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SO2) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SO2) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SO2) ) & 
	         - 1 * W(R12F) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SO2) ) & 
	         + 1 * W(R13F) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SO2) ) )
      JAC(SH,SO) = MM(SH) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R1B) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R2F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R3F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SO) ) & 
	         - 2 * W(R5F) * ( - 3 * MMix/MM(SO) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SO) ) & 
	         + 2 * W(R5B) * ( - 2 * MMix/MM(SO) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SO) ) & 
	         - 2 * W(R6F) * ( - 3 * MMix/MM(SO) ) & 
	         + 2 * W(R6B) * ( - 2 * MMix/MM(SO) ) & 
	         - 2 * W(R7F) * ( - 3 * MMix/MM(SO) ) & 
	         + 2 * W(R7B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SO) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SO) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SO) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SO) ) & 
	         - 1 * W(R10F) * ( 1 / Y(SO) - 3 * MMix/MM(SO) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SO) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SO) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SO) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SO) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R15B) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SO) ) )
      JAC(SH,SOH) = MM(SH) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R1B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R2B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R3F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 2 * W(R5F) * ( - 3 * MMix/MM(SOH) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SOH) ) & 
	         + 2 * W(R5B) * ( - 2 * MMix/MM(SOH) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SOH) ) & 
	         - 2 * W(R6F) * ( - 3 * MMix/MM(SOH) ) & 
	         + 2 * W(R6B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 2 * W(R7F) * ( - 3 * MMix/MM(SOH) ) & 
	         + 2 * W(R7B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R9F) * ( 1 / Y(SOH) - 3 * MMix/MM(SOH) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SOH) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SOH) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SOH) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SOH) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SOH) ) & 
	         + 1 * W(R10B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SOH) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SOH) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R16B) * ( 2 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R23B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) )
      JAC(SH,SH2) = MM(SH) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R2F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R3F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 2 * W(R5F) * ( - 3 * MMix/MM(SH2) ) & 
	         + 2 * W(R5B) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 2 * W(R6F) * ( 1 / Y(SH2) - 3 * MMix/MM(SH2) ) & 
	         + 2 * W(R6B) * ( 2 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 2 * W(R7F) * ( - 3 * MMix/MM(SH2) ) & 
	         + 2 * W(R7B) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SH2) ) + CDOT(SH) * 2 * rho_gp / ( M(MM2) * MM(SH2) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SH2) ) + CDOT(SH) * 2 * rho_gp / ( M(MM2) * MM(SH2) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SH2) ) + CDOT(SH) * 2 * rho_gp / ( M(MM3) * MM(SH2) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SH2) ) + CDOT(SH) * 2 * rho_gp / ( M(MM3) * MM(SH2) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SH2) ) & 
	         + 1 * W(R13F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R22B) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SH2) ) )
      JAC(SH,SH2O) = MM(SH) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R3F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R3B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 2 * W(R5F) * ( - 3 * MMix/MM(SH2O) ) & 
	         + 2 * W(R5B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 2 * W(R6F) * ( - 3 * MMix/MM(SH2O) ) & 
	         + 2 * W(R6B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 2 * W(R7F) * ( 1 / Y(SH2O) - 3 * MMix/MM(SH2O) ) & 
	         + 2 * W(R7B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SH2O) ) + CDOT(SH) * 6.3 * rho_gp / ( M(MM2) * MM(SH2O) ) & 
	         + 1 * W(R9B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) + CDOT(SH) * 6.3 * rho_gp / ( M(MM2) * MM(SH2O) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SH2O) ) + CDOT(SH) * 12 * rho_gp / ( M(MM3) * MM(SH2O) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SH2O) ) + CDOT(SH) * 12 * rho_gp / ( M(MM3) * MM(SH2O) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SH2O) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R15B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R23B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) )
      JAC(SH,SHO2) = MM(SH) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R3F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 2 * W(R5F) * ( - 3 * MMix/MM(SHO2) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SHO2) ) & 
	         + 2 * W(R5B) * ( - 2 * MMix/MM(SHO2) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SHO2) ) & 
	         - 2 * W(R6F) * ( - 3 * MMix/MM(SHO2) ) & 
	         + 2 * W(R6B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 2 * W(R7F) * ( - 3 * MMix/MM(SHO2) ) & 
	         + 2 * W(R7B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SHO2) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SHO2) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SHO2) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SHO2) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SHO2) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SHO2) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SHO2) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SHO2) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R12B) * ( 1 / Y(SHO2) - 1 * MMix/MM(SHO2) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R13B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R15F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R16F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R22B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SHO2) ) )
      JAC(SH,SH2O2) = MM(SH) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R3F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 2 * W(R5F) * ( - 3 * MMix/MM(SH2O2) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SH2O2) ) & 
	         + 2 * W(R5B) * ( - 2 * MMix/MM(SH2O2) ) + CDOT(SH) * rho_gp / ( M(MM1) * MM(SH2O2) ) & 
	         - 2 * W(R6F) * ( - 3 * MMix/MM(SH2O2) ) & 
	         + 2 * W(R6B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 2 * W(R7F) * ( - 3 * MMix/MM(SH2O2) ) & 
	         + 2 * W(R7B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SH2O2) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SH2O2) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SH2O2) ) + CDOT(SH) * rho_gp / ( M(MM2) * MM(SH2O2) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SH2O2) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SH2O2) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SH2O2) ) + CDOT(SH) * rho_gp / ( M(MM3) * MM(SH2O2) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R22F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R23F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SH2O2) ) )
      JAC(SH,ST) = MM(SH) * ( - 1 * W(R1F) * ( -0.67 + -71300000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R1B) * ( -0.234291 + -479604.2721/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R2F) * ( 2.7 + -26190000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R2B) * ( 2.66385 + -20100582.44/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R3F) * ( 1.51 + -14350000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R3B) * ( 1.40625 + -77394524.74/(RGAS * T) - 2 ) / T  & 
	         - 2 * W(R5F) * ( -1 + 0/(RGAS * T) - 3 ) / T  & 
	         + 2 * W(R5B) * ( -0.649765 + -433091487.8/(RGAS * T) - 2 ) / T  & 
	         - 2 * W(R6F) * ( -0.6 + 0/(RGAS * T) - 3 ) / T  & 
	         + 2 * W(R6B) * ( -0.249765 + -433091487.8/(RGAS * T) - 2 ) / T  & 
	         - 2 * W(R7F) * ( -1.25 + 0/(RGAS * T) - 3 ) / T  & 
	         + 2 * W(R7B) * ( -0.899765 + -433091487.8/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R9F) * ( -2 + 0/(RGAS * T) - 3 ) / T  & 
	         + 1 * W(R9B) * ( -1.75351 + -496136012.5/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R10F) * ( -1 + 0/(RGAS * T) - 3 ) / T  & 
	         + 1 * W(R10B) * ( -0.685917 + -427002070.2/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R12F) * ((K2(SH)-K(SH))/(T2-T) / K(SH) - 2 / T ) & 
	         + 1 * W(R12B) * ((K2(SH)-K(SH))/(T2-T) / K(SH) - 1 / T ) & 
	         + 1 * W(R13F) * ( 2.43 + -223850000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R13B) * ( 2.06793 + 4573672.487/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R15F) * ( 0 + -2810000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R15B) * ( 0.694036 + -223457801.5/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R16F) * ( 0 + -2660000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R16B) * ( 0.761631 + -154173859.2/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R22F) * ( 2 + -21760000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R22B) * ( 2.60696 + -88599143.48/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R23F) * ( 0 + -16610000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R23B) * ( 1.26484 + -298007527.4/(RGAS * T) - 2 ) / T  )

      JAC(SO2,SN2) = MM(SO2) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R11F) * ( - 3 * MMix/MM(SN2) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SN2) ) & 
	         - 1 * W(R11B) * ( - 2 * MMix/MM(SN2) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SN2) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SN2) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SN2) ) )
      JAC(SO2,SH) = MM(SO2) * ( - 1 * W(R1F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R11F) * ( - 3 * MMix/MM(SH) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SH) ) & 
	         - 1 * W(R11B) * ( - 2 * MMix/MM(SH) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SH) ) & 
	         - 1 * W(R12F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SH) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R13B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SH) ) )
      JAC(SO2,SO2) = MM(SO2) * ( - 1 * W(R1F) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R11F) * ( - 3 * MMix/MM(SO2) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SO2) ) & 
	         - 1 * W(R11B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SO2) ) & 
	         - 1 * W(R12F) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SO2) ) & 
	         - 1 * W(R13F) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R17B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R18B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R19B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R20B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R21B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) )
      JAC(SO2,SO) = MM(SO2) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R1B) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R11F) * ( 2 / Y(SO) - 3 * MMix/MM(SO) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SO) ) & 
	         - 1 * W(R11B) * ( - 2 * MMix/MM(SO) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SO) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SO) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R17F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SO) ) )
      JAC(SO2,SOH) = MM(SO2) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R1B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R11F) * ( - 3 * MMix/MM(SOH) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SOH) ) & 
	         - 1 * W(R11B) * ( - 2 * MMix/MM(SOH) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SOH) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SOH) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R17B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R18F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R19F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SOH) ) )
      JAC(SO2,SH2) = MM(SO2) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R11F) * ( - 3 * MMix/MM(SH2) ) + CDOT(SO2) * 2.4 * rho_gp / ( M(MM4) * MM(SH2) ) & 
	         - 1 * W(R11B) * ( - 2 * MMix/MM(SH2) ) + CDOT(SO2) * 2.4 * rho_gp / ( M(MM4) * MM(SH2) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SH2) ) & 
	         - 1 * W(R13F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SH2) ) )
      JAC(SO2,SH2O) = MM(SO2) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R11F) * ( - 3 * MMix/MM(SH2O) ) + CDOT(SO2) * 15.4 * rho_gp / ( M(MM4) * MM(SH2O) ) & 
	         - 1 * W(R11B) * ( - 2 * MMix/MM(SH2O) ) + CDOT(SO2) * 15.4 * rho_gp / ( M(MM4) * MM(SH2O) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SH2O) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R18B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R19B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SH2O) ) )
      JAC(SO2,SHO2) = MM(SO2) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R11F) * ( - 3 * MMix/MM(SHO2) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SHO2) ) & 
	         - 1 * W(R11B) * ( - 2 * MMix/MM(SHO2) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SHO2) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R12B) * ( 1 / Y(SHO2) - 1 * MMix/MM(SHO2) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R13B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R17F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R18F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R19F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R20F) * ( 2 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R21F) * ( 2 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SHO2) ) )
      JAC(SO2,SH2O2) = MM(SO2) * ( - 1 * W(R1F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R1B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R11F) * ( - 3 * MMix/MM(SH2O2) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SH2O2) ) & 
	         - 1 * W(R11B) * ( - 2 * MMix/MM(SH2O2) ) + CDOT(SO2) * rho_gp / ( M(MM4) * MM(SH2O2) ) & 
	         - 1 * W(R12F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R12B) * ( - 1 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R20B) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R21B) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) )
      JAC(SO2,ST) = MM(SO2) * ( - 1 * W(R1F) * ( -0.67 + -71300000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R1B) * ( -0.234291 + -479604.2721/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R11F) * ( -1 + 0/(RGAS * T) - 3 ) / T  & 
	         - 1 * W(R11B) * ( -1.12163 + -497822466/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R12F) * ((K2(SO2)-K(SO2))/(T2-T) / K(SO2) - 2 / T ) & 
	         + 1 * W(R12B) * ((K2(SO2)-K(SO2))/(T2-T) / K(SO2) - 1 / T ) & 
	         - 1 * W(R13F) * ( 2.43 + -223850000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R13B) * ( 2.06793 + 4573672.487/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R17F) * ( 0 + 0/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R17B) * ( 0.325921 + -222334254.9/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R18F) * ( 0 + 2090000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R18B) * ( 0.258327 + -289378197.2/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R19F) * ( 0 + -72510000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R19B) * ( 0.258327 + -363978197.2/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R20F) * ( 0 + 6820000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R20B) * ( -0.244887 + -154764529/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R21F) * ( 0 + -50210000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R21B) * ( -0.244887 + -211794529/(RGAS * T) - 2 ) / T  )

      JAC(SO,SN2) = MM(SO) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R2F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SN2) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SN2) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SN2) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SN2) ) & 
	         - 2 * W(R11F) * ( - 3 * MMix/MM(SN2) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SN2) ) & 
	         + 2 * W(R11B) * ( - 2 * MMix/MM(SN2) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SN2) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SN2) ) )
      JAC(SO,SH) = MM(SO) * ( 1 * W(R1F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R2F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R2B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R10F) * ( 1 / Y(SH) - 3 * MMix/MM(SH) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SH) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SH) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SH) ) & 
	         - 2 * W(R11F) * ( - 3 * MMix/MM(SH) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SH) ) & 
	         + 2 * W(R11B) * ( - 2 * MMix/MM(SH) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SH) ) & 
	         + 1 * W(R15F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SH) ) )
      JAC(SO,SO2) = MM(SO) * ( 1 * W(R1F) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R2F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SO2) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SO2) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SO2) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SO2) ) & 
	         - 2 * W(R11F) * ( - 3 * MMix/MM(SO2) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SO2) ) & 
	         + 2 * W(R11B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SO2) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R17B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SO2) ) )
      JAC(SO,SO) = MM(SO) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R1B) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R2F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R4B) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R10F) * ( 1 / Y(SO) - 3 * MMix/MM(SO) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SO) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SO) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SO) ) & 
	         - 2 * W(R11F) * ( 2 / Y(SO) - 3 * MMix/MM(SO) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SO) ) & 
	         + 2 * W(R11B) * ( - 2 * MMix/MM(SO) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SO) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R15B) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R17F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R24F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SO) ) )
      JAC(SO,SOH) = MM(SO) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R1B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R2F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R2B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R4F) * ( 2 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SOH) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SOH) ) & 
	         + 1 * W(R10B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SOH) ) & 
	         - 2 * W(R11F) * ( - 3 * MMix/MM(SOH) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SOH) ) & 
	         + 2 * W(R11B) * ( - 2 * MMix/MM(SOH) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SOH) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R17B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R24B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) )
      JAC(SO,SH2) = MM(SO) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R2F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SH2) ) + CDOT(SO) * 2 * rho_gp / ( M(MM3) * MM(SH2) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SH2) ) + CDOT(SO) * 2 * rho_gp / ( M(MM3) * MM(SH2) ) & 
	         - 2 * W(R11F) * ( - 3 * MMix/MM(SH2) ) + CDOT(SO) * 2.4 * rho_gp / ( M(MM4) * MM(SH2) ) & 
	         + 2 * W(R11B) * ( - 2 * MMix/MM(SH2) ) + CDOT(SO) * 2.4 * rho_gp / ( M(MM4) * MM(SH2) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SH2) ) )
      JAC(SO,SH2O) = MM(SO) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R2F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R4B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SH2O) ) + CDOT(SO) * 12 * rho_gp / ( M(MM3) * MM(SH2O) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SH2O) ) + CDOT(SO) * 12 * rho_gp / ( M(MM3) * MM(SH2O) ) & 
	         - 2 * W(R11F) * ( - 3 * MMix/MM(SH2O) ) + CDOT(SO) * 15.4 * rho_gp / ( M(MM4) * MM(SH2O) ) & 
	         + 2 * W(R11B) * ( - 2 * MMix/MM(SH2O) ) + CDOT(SO) * 15.4 * rho_gp / ( M(MM4) * MM(SH2O) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R15B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SH2O) ) )
      JAC(SO,SHO2) = MM(SO) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R2F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SHO2) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SHO2) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SHO2) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SHO2) ) & 
	         - 2 * W(R11F) * ( - 3 * MMix/MM(SHO2) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SHO2) ) & 
	         + 2 * W(R11B) * ( - 2 * MMix/MM(SHO2) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SHO2) ) & 
	         + 1 * W(R15F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R17F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R24B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) )
      JAC(SO,SH2O2) = MM(SO) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R2F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R10F) * ( - 3 * MMix/MM(SH2O2) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SH2O2) ) & 
	         + 1 * W(R10B) * ( - 2 * MMix/MM(SH2O2) ) + CDOT(SO) * rho_gp / ( M(MM3) * MM(SH2O2) ) & 
	         - 2 * W(R11F) * ( - 3 * MMix/MM(SH2O2) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SH2O2) ) & 
	         + 2 * W(R11B) * ( - 2 * MMix/MM(SH2O2) ) + CDOT(SO) * rho_gp / ( M(MM4) * MM(SH2O2) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R24F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SH2O2) ) )
      JAC(SO,ST) = MM(SO) * ( 1 * W(R1F) * ( -0.67 + -71300000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R1B) * ( -0.234291 + -479604.2721/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R2F) * ( 2.7 + -26190000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R2B) * ( 2.66385 + -20100582.44/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R4F) * ( 2.4 + 8830000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R4B) * ( 2.33241 + -60303942.3/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R10F) * ( -1 + 0/(RGAS * T) - 3 ) / T  & 
	         + 1 * W(R10B) * ( -0.685917 + -427002070.2/(RGAS * T) - 2 ) / T  & 
	         - 2 * W(R11F) * ( -1 + 0/(RGAS * T) - 3 ) / T  & 
	         + 2 * W(R11B) * ( -1.12163 + -497822466/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R15F) * ( 0 + -2810000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R15B) * ( 0.694036 + -223457801.5/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R17F) * ( 0 + 0/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R17B) * ( 0.325921 + -222334254.9/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R24F) * ( 2 + -16610000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R24B) * ( 2.57081 + -77359725.92/(RGAS * T) - 2 ) / T  )

      JAC(SOH,SN2) = MM(SOH) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 2 * W(R4F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 2 * W(R4B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SN2) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SN2) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SN2) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SN2) ) & 
	         + 1 * W(R10F) * ( - 3 * MMix/MM(SN2) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SN2) ) & 
	         - 1 * W(R10B) * ( - 2 * MMix/MM(SN2) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SN2) ) & 
	         - 2 * W(R14F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 2 * W(R14B) * ( - 1 * MMix/MM(SN2) ) & 
	         + 2 * W(R16F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 2 * W(R16B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SN2) ) )
      JAC(SOH,SH) = MM(SOH) * ( 1 * W(R1F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R2B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R3B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 2 * W(R4F) * ( - 2 * MMix/MM(SH) ) & 
	         + 2 * W(R4B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R9F) * ( 1 / Y(SH) - 3 * MMix/MM(SH) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SH) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SH) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SH) ) & 
	         + 1 * W(R10F) * ( 1 / Y(SH) - 3 * MMix/MM(SH) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SH) ) & 
	         - 1 * W(R10B) * ( - 2 * MMix/MM(SH) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SH) ) & 
	         - 2 * W(R14F) * ( - 2 * MMix/MM(SH) ) & 
	         + 2 * W(R14B) * ( - 1 * MMix/MM(SH) ) & 
	         + 2 * W(R16F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 2 * W(R16B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R23F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SH) ) )
      JAC(SOH,SO2) = MM(SOH) * ( 1 * W(R1F) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 2 * W(R4F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 2 * W(R4B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SO2) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SO2) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SO2) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SO2) ) & 
	         + 1 * W(R10F) * ( - 3 * MMix/MM(SO2) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SO2) ) & 
	         - 1 * W(R10B) * ( - 2 * MMix/MM(SO2) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SO2) ) & 
	         - 2 * W(R14F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 2 * W(R14B) * ( - 1 * MMix/MM(SO2) ) & 
	         + 2 * W(R16F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 2 * W(R16B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R17B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R18B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R19B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SO2) ) )
      JAC(SOH,SO) = MM(SOH) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R1B) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R2F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SO) ) & 
	         - 2 * W(R4F) * ( - 2 * MMix/MM(SO) ) & 
	         + 2 * W(R4B) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SO) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SO) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SO) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SO) ) & 
	         + 1 * W(R10F) * ( 1 / Y(SO) - 3 * MMix/MM(SO) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SO) ) & 
	         - 1 * W(R10B) * ( - 2 * MMix/MM(SO) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SO) ) & 
	         - 2 * W(R14F) * ( - 2 * MMix/MM(SO) ) & 
	         + 2 * W(R14B) * ( - 1 * MMix/MM(SO) ) & 
	         + 2 * W(R16F) * ( - 2 * MMix/MM(SO) ) & 
	         - 2 * W(R16B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R17F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R24F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SO) ) )
      JAC(SOH,SOH) = MM(SOH) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R1B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R2B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R3F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 2 * W(R4F) * ( 2 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 2 * W(R4B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R9F) * ( 1 / Y(SOH) - 3 * MMix/MM(SOH) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SOH) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SOH) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SOH) ) & 
	         + 1 * W(R10F) * ( - 3 * MMix/MM(SOH) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SOH) ) & 
	         - 1 * W(R10B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SOH) ) & 
	         - 2 * W(R14F) * ( 2 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 2 * W(R14B) * ( - 1 * MMix/MM(SOH) ) & 
	         + 2 * W(R16F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 2 * W(R16B) * ( 2 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R17B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R18F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R19F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R23B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R24B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R25F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R26F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SOH) ) )
      JAC(SOH,SH2) = MM(SOH) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R2F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R3F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 2 * W(R4F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 2 * W(R4B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SH2) ) + CDOT(SOH) * 2 * rho_gp / ( M(MM2) * MM(SH2) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SH2) ) + CDOT(SOH) * 2 * rho_gp / ( M(MM2) * MM(SH2) ) & 
	         + 1 * W(R10F) * ( - 3 * MMix/MM(SH2) ) + CDOT(SOH) * 2 * rho_gp / ( M(MM3) * MM(SH2) ) & 
	         - 1 * W(R10B) * ( - 2 * MMix/MM(SH2) ) + CDOT(SOH) * 2 * rho_gp / ( M(MM3) * MM(SH2) ) & 
	         - 2 * W(R14F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 2 * W(R14B) * ( - 1 * MMix/MM(SH2) ) & 
	         + 2 * W(R16F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 2 * W(R16B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SH2) ) )
      JAC(SOH,SH2O) = MM(SOH) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R3B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 2 * W(R4F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 2 * W(R4B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SH2O) ) + CDOT(SOH) * 6.3 * rho_gp / ( M(MM2) * MM(SH2O) ) & 
	         + 1 * W(R9B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) + CDOT(SOH) * 6.3 * rho_gp / ( M(MM2) * MM(SH2O) ) & 
	         + 1 * W(R10F) * ( - 3 * MMix/MM(SH2O) ) + CDOT(SOH) * 12 * rho_gp / ( M(MM3) * MM(SH2O) ) & 
	         - 1 * W(R10B) * ( - 2 * MMix/MM(SH2O) ) + CDOT(SOH) * 12 * rho_gp / ( M(MM3) * MM(SH2O) ) & 
	         - 2 * W(R14F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 2 * W(R14B) * ( - 1 * MMix/MM(SH2O) ) & 
	         + 2 * W(R16F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 2 * W(R16B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R18B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R19B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R23B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R25B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R26B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) )
      JAC(SOH,SHO2) = MM(SOH) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 2 * W(R4F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 2 * W(R4B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SHO2) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SHO2) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SHO2) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SHO2) ) & 
	         + 1 * W(R10F) * ( - 3 * MMix/MM(SHO2) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SHO2) ) & 
	         - 1 * W(R10B) * ( - 2 * MMix/MM(SHO2) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SHO2) ) & 
	         - 2 * W(R14F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 2 * W(R14B) * ( - 1 * MMix/MM(SHO2) ) & 
	         + 2 * W(R16F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 2 * W(R16B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R17F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R18F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R19F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R24B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R25B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R26B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) )
      JAC(SOH,SH2O2) = MM(SOH) * ( 1 * W(R1F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R1B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R2F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R2B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 2 * W(R4F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 2 * W(R4B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R9F) * ( - 3 * MMix/MM(SH2O2) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SH2O2) ) & 
	         + 1 * W(R9B) * ( - 2 * MMix/MM(SH2O2) ) + CDOT(SOH) * rho_gp / ( M(MM2) * MM(SH2O2) ) & 
	         + 1 * W(R10F) * ( - 3 * MMix/MM(SH2O2) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SH2O2) ) & 
	         - 1 * W(R10B) * ( - 2 * MMix/MM(SH2O2) ) + CDOT(SOH) * rho_gp / ( M(MM3) * MM(SH2O2) ) & 
	         - 2 * W(R14F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 2 * W(R14B) * ( 1 / Y(SH2O2) - 1 * MMix/MM(SH2O2) ) & 
	         + 2 * W(R16F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 2 * W(R16B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R17F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R17B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R23F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R24F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R25F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R26F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SH2O2) ) )
      JAC(SOH,ST) = MM(SOH) * ( 1 * W(R1F) * ( -0.67 + -71300000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R1B) * ( -0.234291 + -479604.2721/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R2F) * ( 2.7 + -26190000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R2B) * ( 2.66385 + -20100582.44/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R3F) * ( 1.51 + -14350000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R3B) * ( 1.40625 + -77394524.74/(RGAS * T) - 2 ) / T  & 
	         - 2 * W(R4F) * ( 2.4 + 8830000/(RGAS * T) - 2 ) / T  & 
	         + 2 * W(R4B) * ( 2.33241 + -60303942.3/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R9F) * ( -2 + 0/(RGAS * T) - 3 ) / T  & 
	         + 1 * W(R9B) * ( -1.75351 + -496136012.5/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R10F) * ( -1 + 0/(RGAS * T) - 3 ) / T  & 
	         - 1 * W(R10B) * ( -0.685917 + -427002070.2/(RGAS * T) - 2 ) / T  & 
	         - 2 * W(R14F) * ((K2(SOH)-K(SOH))/(T2-T) / K(SOH) - 2 / T ) & 
	         + 2 * W(R14B) * ((K2(SOH)-K(SOH))/(T2-T) / K(SOH) - 1 / T ) & 
	         + 2 * W(R16F) * ( 0 + -2660000/(RGAS * T) - 2 ) / T  & 
	         - 2 * W(R16B) * ( 0.761631 + -154173859.2/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R17F) * ( 0 + 0/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R17B) * ( 0.325921 + -222334254.9/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R18F) * ( 0 + 2090000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R18B) * ( 0.258327 + -289378197.2/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R19F) * ( 0 + -72510000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R19B) * ( 0.258327 + -363978197.2/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R23F) * ( 0 + -16610000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R23B) * ( 1.26484 + -298007527.4/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R24F) * ( 2 + -16610000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R24B) * ( 2.57081 + -77359725.92/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R25F) * ( 0 + -1790000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R25B) * ( 0.503214 + -131673668.2/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R26F) * ( -7 + -157320000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R26B) * ( -6.49679 + -287203668.2/(RGAS * T) - 2 ) / T  )

      JAC(SH2,SN2) = MM(SH2) * ( - 1 * W(R2F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R5F) * ( - 3 * MMix/MM(SN2) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SN2) ) & 
	         - 1 * W(R5B) * ( - 2 * MMix/MM(SN2) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SN2) ) & 
	         - 1 * W(R6F) * ( - 3 * MMix/MM(SN2) ) & 
	         + 2 * W(R6F) * ( - 3 * MMix/MM(SN2) ) & 
	         - 2 * W(R6B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R6B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SN2) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SN2) ) )
      JAC(SH2,SH) = MM(SH2) * ( - 1 * W(R2F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R2B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R3B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R5F) * ( 2 / Y(SH) - 3 * MMix/MM(SH) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SH) ) & 
	         - 1 * W(R5B) * ( - 2 * MMix/MM(SH) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SH) ) & 
	         - 1 * W(R6F) * ( 2 / Y(SH) - 3 * MMix/MM(SH) ) & 
	         + 2 * W(R6F) * ( 2 / Y(SH) - 3 * MMix/MM(SH) ) & 
	         - 2 * W(R6B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R6B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R7F) * ( 2 / Y(SH) - 3 * MMix/MM(SH) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R13B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R22F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SH) ) )
      JAC(SH2,SO2) = MM(SH2) * ( - 1 * W(R2F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R5F) * ( - 3 * MMix/MM(SO2) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SO2) ) & 
	         - 1 * W(R5B) * ( - 2 * MMix/MM(SO2) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SO2) ) & 
	         - 1 * W(R6F) * ( - 3 * MMix/MM(SO2) ) & 
	         + 2 * W(R6F) * ( - 3 * MMix/MM(SO2) ) & 
	         - 2 * W(R6B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R6B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SO2) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R13F) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SO2) ) )
      JAC(SH2,SO) = MM(SH2) * ( - 1 * W(R2F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R5F) * ( - 3 * MMix/MM(SO) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SO) ) & 
	         - 1 * W(R5B) * ( - 2 * MMix/MM(SO) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SO) ) & 
	         - 1 * W(R6F) * ( - 3 * MMix/MM(SO) ) & 
	         + 2 * W(R6F) * ( - 3 * MMix/MM(SO) ) & 
	         - 2 * W(R6B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R6B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SO) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SO) ) )
      JAC(SH2,SOH) = MM(SH2) * ( - 1 * W(R2F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R2B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R3F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R5F) * ( - 3 * MMix/MM(SOH) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SOH) ) & 
	         - 1 * W(R5B) * ( - 2 * MMix/MM(SOH) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SOH) ) & 
	         - 1 * W(R6F) * ( - 3 * MMix/MM(SOH) ) & 
	         + 2 * W(R6F) * ( - 3 * MMix/MM(SOH) ) & 
	         - 2 * W(R6B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R6B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SOH) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SOH) ) )
      JAC(SH2,SH2) = MM(SH2) * ( - 1 * W(R2F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R3F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R5F) * ( - 3 * MMix/MM(SH2) ) & 
	         - 1 * W(R5B) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R6F) * ( 1 / Y(SH2) - 3 * MMix/MM(SH2) ) & 
	         + 2 * W(R6F) * ( 1 / Y(SH2) - 3 * MMix/MM(SH2) ) & 
	         - 2 * W(R6B) * ( 2 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R6B) * ( 2 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SH2) ) & 
	         - 1 * W(R7B) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R13F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R22B) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) )
      JAC(SH2,SH2O) = MM(SH2) * ( - 1 * W(R2F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R3B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R5F) * ( - 3 * MMix/MM(SH2O) ) & 
	         - 1 * W(R5B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R6F) * ( - 3 * MMix/MM(SH2O) ) & 
	         + 2 * W(R6F) * ( - 3 * MMix/MM(SH2O) ) & 
	         - 2 * W(R6B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R6B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R7F) * ( 1 / Y(SH2O) - 3 * MMix/MM(SH2O) ) & 
	         - 1 * W(R7B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SH2O) ) )
      JAC(SH2,SHO2) = MM(SH2) * ( - 1 * W(R2F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R5F) * ( - 3 * MMix/MM(SHO2) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SHO2) ) & 
	         - 1 * W(R5B) * ( - 2 * MMix/MM(SHO2) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SHO2) ) & 
	         - 1 * W(R6F) * ( - 3 * MMix/MM(SHO2) ) & 
	         + 2 * W(R6F) * ( - 3 * MMix/MM(SHO2) ) & 
	         - 2 * W(R6B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R6B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SHO2) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R13B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R22B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) )
      JAC(SH2,SH2O2) = MM(SH2) * ( - 1 * W(R2F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R2B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R3F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R3B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R5F) * ( - 3 * MMix/MM(SH2O2) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SH2O2) ) & 
	         - 1 * W(R5B) * ( - 2 * MMix/MM(SH2O2) ) + CDOT(SH2) * rho_gp / ( M(MM1) * MM(SH2O2) ) & 
	         - 1 * W(R6F) * ( - 3 * MMix/MM(SH2O2) ) & 
	         + 2 * W(R6F) * ( - 3 * MMix/MM(SH2O2) ) & 
	         - 2 * W(R6B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R6B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R13F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R13B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R22F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SH2O2) ) )
      JAC(SH2,ST) = MM(SH2) * ( - 1 * W(R2F) * ( 2.7 + -26190000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R2B) * ( 2.66385 + -20100582.44/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R3F) * ( 1.51 + -14350000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R3B) * ( 1.40625 + -77394524.74/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R5F) * ( -1 + 0/(RGAS * T) - 3 ) / T  & 
	         - 1 * W(R5B) * ( -0.649765 + -433091487.8/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R6F) * ( -0.6 + 0/(RGAS * T) - 3 ) / T  & 
	         + 2 * W(R6F) * ( -0.6 + 0/(RGAS * T) - 3 ) / T  & 
	         - 2 * W(R6B) * ( -0.249765 + -433091487.8/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R6B) * ( -0.249765 + -433091487.8/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R7F) * ( -1.25 + 0/(RGAS * T) - 3 ) / T  & 
	         - 1 * W(R7B) * ( -0.899765 + -433091487.8/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R13F) * ( 2.43 + -223850000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R13B) * ( 2.06793 + 4573672.487/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R22F) * ( 2 + -21760000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R22B) * ( 2.60696 + -88599143.48/(RGAS * T) - 2 ) / T  )

      JAC(SH2O,SN2) = MM(SH2O) * ( 1 * W(R3F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R7F) * ( - 3 * MMix/MM(SN2) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SN2) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R7B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R9F) * ( - 3 * MMix/MM(SN2) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SN2) ) & 
	         - 1 * W(R9B) * ( - 2 * MMix/MM(SN2) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SN2) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SN2) ) )
      JAC(SH2O,SH) = MM(SH2O) * ( 1 * W(R3F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R3B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R7F) * ( 2 / Y(SH) - 3 * MMix/MM(SH) ) & 
	         + 1 * W(R7F) * ( 2 / Y(SH) - 3 * MMix/MM(SH) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R7B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R9F) * ( 1 / Y(SH) - 3 * MMix/MM(SH) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SH) ) & 
	         - 1 * W(R9B) * ( - 2 * MMix/MM(SH) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SH) ) & 
	         + 1 * W(R15F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R23F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SH) ) )
      JAC(SH2O,SO2) = MM(SH2O) * ( 1 * W(R3F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R7F) * ( - 3 * MMix/MM(SO2) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SO2) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R7B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R9F) * ( - 3 * MMix/MM(SO2) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SO2) ) & 
	         - 1 * W(R9B) * ( - 2 * MMix/MM(SO2) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SO2) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R18B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R19B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SO2) ) )
      JAC(SH2O,SO) = MM(SH2O) * ( 1 * W(R3F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R4B) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R7F) * ( - 3 * MMix/MM(SO) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SO) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R7B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R9F) * ( - 3 * MMix/MM(SO) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SO) ) & 
	         - 1 * W(R9B) * ( - 2 * MMix/MM(SO) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SO) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R15B) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SO) ) )
      JAC(SH2O,SOH) = MM(SH2O) * ( 1 * W(R3F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R4F) * ( 2 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R7F) * ( - 3 * MMix/MM(SOH) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SOH) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R7B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R9F) * ( 1 / Y(SOH) - 3 * MMix/MM(SOH) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SOH) ) & 
	         - 1 * W(R9B) * ( - 2 * MMix/MM(SOH) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SOH) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R18F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R19F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R23B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R25F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R26F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SOH) ) )
      JAC(SH2O,SH2) = MM(SH2O) * ( 1 * W(R3F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R7F) * ( - 3 * MMix/MM(SH2) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SH2) ) & 
	         - 1 * W(R7B) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R7B) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R9F) * ( - 3 * MMix/MM(SH2) ) + CDOT(SH2O) * 2 * rho_gp / ( M(MM2) * MM(SH2) ) & 
	         - 1 * W(R9B) * ( - 2 * MMix/MM(SH2) ) + CDOT(SH2O) * 2 * rho_gp / ( M(MM2) * MM(SH2) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SH2) ) )
      JAC(SH2O,SH2O) = MM(SH2O) * ( 1 * W(R3F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R3B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R4B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R7F) * ( 1 / Y(SH2O) - 3 * MMix/MM(SH2O) ) & 
	         + 1 * W(R7F) * ( 1 / Y(SH2O) - 3 * MMix/MM(SH2O) ) & 
	         - 1 * W(R7B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R7B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R9F) * ( - 3 * MMix/MM(SH2O) ) + CDOT(SH2O) * 6.3 * rho_gp / ( M(MM2) * MM(SH2O) ) & 
	         - 1 * W(R9B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) + CDOT(SH2O) * 6.3 * rho_gp / ( M(MM2) * MM(SH2O) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R15B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R18B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R19B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R23B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R25B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R26B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) )
      JAC(SH2O,SHO2) = MM(SH2O) * ( 1 * W(R3F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R7F) * ( - 3 * MMix/MM(SHO2) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SHO2) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R7B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R9F) * ( - 3 * MMix/MM(SHO2) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SHO2) ) & 
	         - 1 * W(R9B) * ( - 2 * MMix/MM(SHO2) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SHO2) ) & 
	         + 1 * W(R15F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R18F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R19F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R23F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R25B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R26B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) )
      JAC(SH2O,SH2O2) = MM(SH2O) * ( 1 * W(R3F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R3B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R4F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R4B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R7F) * ( - 3 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R7F) * ( - 3 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R7B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R7B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R9F) * ( - 3 * MMix/MM(SH2O2) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SH2O2) ) & 
	         - 1 * W(R9B) * ( - 2 * MMix/MM(SH2O2) ) + CDOT(SH2O) * rho_gp / ( M(MM2) * MM(SH2O2) ) & 
	         + 1 * W(R15F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R15B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R18F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R18B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R19F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R19B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R23F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R23B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R25F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R26F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SH2O2) ) )
      JAC(SH2O,ST) = MM(SH2O) * ( 1 * W(R3F) * ( 1.51 + -14350000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R3B) * ( 1.40625 + -77394524.74/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R4F) * ( 2.4 + 8830000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R4B) * ( 2.33241 + -60303942.3/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R7F) * ( -1.25 + 0/(RGAS * T) - 3 ) / T  & 
	         + 1 * W(R7F) * ( -1.25 + 0/(RGAS * T) - 3 ) / T  & 
	         - 1 * W(R7B) * ( -0.899765 + -433091487.8/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R7B) * ( -0.899765 + -433091487.8/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R9F) * ( -2 + 0/(RGAS * T) - 3 ) / T  & 
	         - 1 * W(R9B) * ( -1.75351 + -496136012.5/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R15F) * ( 0 + -2810000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R15B) * ( 0.694036 + -223457801.5/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R18F) * ( 0 + 2090000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R18B) * ( 0.258327 + -289378197.2/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R19F) * ( 0 + -72510000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R19B) * ( 0.258327 + -363978197.2/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R23F) * ( 0 + -16610000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R23B) * ( 1.26484 + -298007527.4/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R25F) * ( 0 + -1790000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R25B) * ( 0.503214 + -131673668.2/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R26F) * ( -7 + -157320000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R26B) * ( -6.49679 + -287203668.2/(RGAS * T) - 2 ) / T  )

      JAC(SHO2,SN2) = MM(SHO2) * ( 1 * W(R12F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R12B) * ( - 1 * MMix/MM(SN2) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 2 * W(R20F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 2 * W(R20B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 2 * W(R21F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 2 * W(R21B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SN2) ) )
      JAC(SHO2,SH) = MM(SHO2) * ( 1 * W(R12F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R12B) * ( - 1 * MMix/MM(SH) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R13B) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R15F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R16F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SH) ) & 
	         - 2 * W(R20F) * ( - 2 * MMix/MM(SH) ) & 
	         + 2 * W(R20B) * ( - 2 * MMix/MM(SH) ) & 
	         - 2 * W(R21F) * ( - 2 * MMix/MM(SH) ) & 
	         + 2 * W(R21B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R22F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SH) ) )
      JAC(SHO2,SO2) = MM(SHO2) * ( 1 * W(R12F) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R12B) * ( - 1 * MMix/MM(SO2) ) & 
	         + 1 * W(R13F) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R17B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R18B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R19B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 2 * W(R20F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 2 * W(R20B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 2 * W(R21F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 2 * W(R21B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SO2) ) )
      JAC(SHO2,SO) = MM(SHO2) * ( 1 * W(R12F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R12B) * ( - 1 * MMix/MM(SO) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R15B) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R17F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SO) ) & 
	         - 2 * W(R20F) * ( - 2 * MMix/MM(SO) ) & 
	         + 2 * W(R20B) * ( - 2 * MMix/MM(SO) ) & 
	         - 2 * W(R21F) * ( - 2 * MMix/MM(SO) ) & 
	         + 2 * W(R21B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R24F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SO) ) )
      JAC(SHO2,SOH) = MM(SHO2) * ( 1 * W(R12F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R12B) * ( - 1 * MMix/MM(SOH) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R16B) * ( 2 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R17B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R18F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R19F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 2 * W(R20F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 2 * W(R20B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 2 * W(R21F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 2 * W(R21B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R24B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R25F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R26F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SOH) ) )
      JAC(SHO2,SH2) = MM(SHO2) * ( 1 * W(R12F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R12B) * ( - 1 * MMix/MM(SH2) ) & 
	         + 1 * W(R13F) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 2 * W(R20F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 2 * W(R20B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 2 * W(R21F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 2 * W(R21B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R22B) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SH2) ) )
      JAC(SHO2,SH2O) = MM(SHO2) * ( 1 * W(R12F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R12B) * ( - 1 * MMix/MM(SH2O) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R15B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R18B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R19B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 2 * W(R20F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 2 * W(R20B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 2 * W(R21F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 2 * W(R21B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R25B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R26B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) )
      JAC(SHO2,SHO2) = MM(SHO2) * ( 1 * W(R12F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R12B) * ( 1 / Y(SHO2) - 1 * MMix/MM(SHO2) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R13B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R15F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R16F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R17F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R18F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R19F) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 2 * W(R20F) * ( 2 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 2 * W(R20B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 2 * W(R21F) * ( 2 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 2 * W(R21B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R22F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R22B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R24F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R24B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R25F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R25B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R26F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R26B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) )
      JAC(SHO2,SH2O2) = MM(SHO2) * ( 1 * W(R12F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R12B) * ( - 1 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R13F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R13B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R15F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R15B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R16F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R16B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R17F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R17B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R18F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R18B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R19F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R19B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 2 * W(R20F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 2 * W(R20B) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 2 * W(R21F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 2 * W(R21B) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R22F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R22B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R24F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R24B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R25F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R25B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R26F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R26B) * ( - 2 * MMix/MM(SH2O2) ) )
      JAC(SHO2,ST) = MM(SHO2) * ( 1 * W(R12F) * ((K2(SHO2)-K(SHO2))/(T2-T) / K(SHO2) - 2 / T ) & 
	         - 1 * W(R12B) * ((K2(SHO2)-K(SHO2))/(T2-T) / K(SHO2) - 1 / T ) & 
	         + 1 * W(R13F) * ( 2.43 + -223850000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R13B) * ( 2.06793 + 4573672.487/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R15F) * ( 0 + -2810000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R15B) * ( 0.694036 + -223457801.5/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R16F) * ( 0 + -2660000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R16B) * ( 0.761631 + -154173859.2/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R17F) * ( 0 + 0/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R17B) * ( 0.325921 + -222334254.9/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R18F) * ( 0 + 2090000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R18B) * ( 0.258327 + -289378197.2/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R19F) * ( 0 + -72510000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R19B) * ( 0.258327 + -363978197.2/(RGAS * T) - 2 ) / T  & 
	         - 2 * W(R20F) * ( 0 + 6820000/(RGAS * T) - 2 ) / T  & 
	         + 2 * W(R20B) * ( -0.244887 + -154764529/(RGAS * T) - 2 ) / T  & 
	         - 2 * W(R21F) * ( 0 + -50210000/(RGAS * T) - 2 ) / T  & 
	         + 2 * W(R21B) * ( -0.244887 + -211794529/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R22F) * ( 2 + -21760000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R22B) * ( 2.60696 + -88599143.48/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R24F) * ( 2 + -16610000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R24B) * ( 2.57081 + -77359725.92/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R25F) * ( 0 + -1790000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R25B) * ( 0.503214 + -131673668.2/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R26F) * ( -7 + -157320000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R26B) * ( -6.49679 + -287203668.2/(RGAS * T) - 2 ) / T  )

      JAC(SH2O2,SN2) = MM(SH2O2) * ( 1 * W(R14F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R14B) * ( - 1 * MMix/MM(SN2) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SN2) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SN2) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SN2) ) )
      JAC(SH2O2,SH) = MM(SH2O2) * ( 1 * W(R14F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R14B) * ( - 1 * MMix/MM(SH) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R22F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R23F) * ( 1 / Y(SH) - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SH) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SH) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SH) ) )
      JAC(SH2O2,SO2) = MM(SH2O2) * ( 1 * W(R14F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R14B) * ( - 1 * MMix/MM(SO2) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R20B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R21B) * ( 1 / Y(SO2) - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SO2) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SO2) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SO2) ) )
      JAC(SH2O2,SO) = MM(SH2O2) * ( 1 * W(R14F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R14B) * ( - 1 * MMix/MM(SO) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R24F) * ( 1 / Y(SO) - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SO) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SO) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SO) ) )
      JAC(SH2O2,SOH) = MM(SH2O2) * ( 1 * W(R14F) * ( 2 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R14B) * ( - 1 * MMix/MM(SOH) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R23B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R24B) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R25F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SOH) ) & 
	         - 1 * W(R26F) * ( 1 / Y(SOH) - 2 * MMix/MM(SOH) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SOH) ) )
      JAC(SH2O2,SH2) = MM(SH2O2) * ( 1 * W(R14F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R14B) * ( - 1 * MMix/MM(SH2) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R22B) * ( 1 / Y(SH2) - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SH2) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SH2) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SH2) ) )
      JAC(SH2O2,SH2O) = MM(SH2O2) * ( 1 * W(R14F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R14B) * ( - 1 * MMix/MM(SH2O) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R23B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R25B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SH2O) ) & 
	         + 1 * W(R26B) * ( 1 / Y(SH2O) - 2 * MMix/MM(SH2O) ) )
      JAC(SH2O2,SHO2) = MM(SH2O2) * ( 1 * W(R14F) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R14B) * ( - 1 * MMix/MM(SHO2) ) & 
	         + 1 * W(R20F) * ( 2 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R20B) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R21F) * ( 2 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R21B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R22F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R22B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R23F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R24F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R24B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R25F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R25B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) & 
	         - 1 * W(R26F) * ( - 2 * MMix/MM(SHO2) ) & 
	         + 1 * W(R26B) * ( 1 / Y(SHO2) - 2 * MMix/MM(SHO2) ) )
      JAC(SH2O2,SH2O2) = MM(SH2O2) * ( 1 * W(R14F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R14B) * ( 1 / Y(SH2O2) - 1 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R20F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R20B) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R21F) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R21B) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R22F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R22B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R23F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R23B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R24F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R24B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R25F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R25B) * ( - 2 * MMix/MM(SH2O2) ) & 
	         - 1 * W(R26F) * ( 1 / Y(SH2O2) - 2 * MMix/MM(SH2O2) ) & 
	         + 1 * W(R26B) * ( - 2 * MMix/MM(SH2O2) ) )
      JAC(SH2O2,ST) = MM(SH2O2) * ( 1 * W(R14F) * ((K2(SH2O2)-K(SH2O2))/(T2-T) / K(SH2O2) - 2 / T ) & 
	         - 1 * W(R14B) * ((K2(SH2O2)-K(SH2O2))/(T2-T) / K(SH2O2) - 1 / T ) & 
	         + 1 * W(R20F) * ( 0 + 6820000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R20B) * ( -0.244887 + -154764529/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R21F) * ( 0 + -50210000/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R21B) * ( -0.244887 + -211794529/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R22F) * ( 2 + -21760000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R22B) * ( 2.60696 + -88599143.48/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R23F) * ( 0 + -16610000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R23B) * ( 1.26484 + -298007527.4/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R24F) * ( 2 + -16610000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R24B) * ( 2.57081 + -77359725.92/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R25F) * ( 0 + -1790000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R25B) * ( 0.503214 + -131673668.2/(RGAS * T) - 2 ) / T  & 
	         - 1 * W(R26F) * ( -7 + -157320000/(RGAS * T) - 2 ) / T  & 
	         + 1 * W(R26B) * ( -6.49679 + -287203668.2/(RGAS * T) - 2 ) / T  )

      JAC(ST,SN2) = -( H(SN2) * JAC(SN2,SN2) - CP(SN2) * MM(SN2) * CDOT(SN2) / Cp_mix & 
	             + H(SH) * JAC(SH,SN2) - CP(SN2) * MM(SH) * CDOT(SH) / Cp_mix & 
	             + H(SO2) * JAC(SO2,SN2) - CP(SN2) * MM(SO2) * CDOT(SO2) / Cp_mix & 
	             + H(SO) * JAC(SO,SN2) - CP(SN2) * MM(SO) * CDOT(SO) / Cp_mix & 
	             + H(SOH) * JAC(SOH,SN2) - CP(SN2) * MM(SOH) * CDOT(SOH) / Cp_mix & 
	             + H(SH2) * JAC(SH2,SN2) - CP(SN2) * MM(SH2) * CDOT(SH2) / Cp_mix & 
	             + H(SH2O) * JAC(SH2O,SN2) - CP(SN2) * MM(SH2O) * CDOT(SH2O) / Cp_mix & 
	             + H(SHO2) * JAC(SHO2,SN2) - CP(SN2) * MM(SHO2) * CDOT(SHO2) / Cp_mix & 
	             + H(SH2O2) * JAC(SH2O2,SN2) - CP(SN2) * MM(SH2O2) * CDOT(SH2O2) / Cp_mix) / Cp_mix
      JAC(ST,SH) = -( H(SN2) * JAC(SN2,SH) - CP(SH) * MM(SN2) * CDOT(SN2) / Cp_mix & 
	             + H(SH) * JAC(SH,SH) - CP(SH) * MM(SH) * CDOT(SH) / Cp_mix & 
	             + H(SO2) * JAC(SO2,SH) - CP(SH) * MM(SO2) * CDOT(SO2) / Cp_mix & 
	             + H(SO) * JAC(SO,SH) - CP(SH) * MM(SO) * CDOT(SO) / Cp_mix & 
	             + H(SOH) * JAC(SOH,SH) - CP(SH) * MM(SOH) * CDOT(SOH) / Cp_mix & 
	             + H(SH2) * JAC(SH2,SH) - CP(SH) * MM(SH2) * CDOT(SH2) / Cp_mix & 
	             + H(SH2O) * JAC(SH2O,SH) - CP(SH) * MM(SH2O) * CDOT(SH2O) / Cp_mix & 
	             + H(SHO2) * JAC(SHO2,SH) - CP(SH) * MM(SHO2) * CDOT(SHO2) / Cp_mix & 
	             + H(SH2O2) * JAC(SH2O2,SH) - CP(SH) * MM(SH2O2) * CDOT(SH2O2) / Cp_mix) / Cp_mix
      JAC(ST,SO2) = -( H(SN2) * JAC(SN2,SO2) - CP(SO2) * MM(SN2) * CDOT(SN2) / Cp_mix & 
	             + H(SH) * JAC(SH,SO2) - CP(SO2) * MM(SH) * CDOT(SH) / Cp_mix & 
	             + H(SO2) * JAC(SO2,SO2) - CP(SO2) * MM(SO2) * CDOT(SO2) / Cp_mix & 
	             + H(SO) * JAC(SO,SO2) - CP(SO2) * MM(SO) * CDOT(SO) / Cp_mix & 
	             + H(SOH) * JAC(SOH,SO2) - CP(SO2) * MM(SOH) * CDOT(SOH) / Cp_mix & 
	             + H(SH2) * JAC(SH2,SO2) - CP(SO2) * MM(SH2) * CDOT(SH2) / Cp_mix & 
	             + H(SH2O) * JAC(SH2O,SO2) - CP(SO2) * MM(SH2O) * CDOT(SH2O) / Cp_mix & 
	             + H(SHO2) * JAC(SHO2,SO2) - CP(SO2) * MM(SHO2) * CDOT(SHO2) / Cp_mix & 
	             + H(SH2O2) * JAC(SH2O2,SO2) - CP(SO2) * MM(SH2O2) * CDOT(SH2O2) / Cp_mix) / Cp_mix
      JAC(ST,SO) = -( H(SN2) * JAC(SN2,SO) - CP(SO) * MM(SN2) * CDOT(SN2) / Cp_mix & 
	             + H(SH) * JAC(SH,SO) - CP(SO) * MM(SH) * CDOT(SH) / Cp_mix & 
	             + H(SO2) * JAC(SO2,SO) - CP(SO) * MM(SO2) * CDOT(SO2) / Cp_mix & 
	             + H(SO) * JAC(SO,SO) - CP(SO) * MM(SO) * CDOT(SO) / Cp_mix & 
	             + H(SOH) * JAC(SOH,SO) - CP(SO) * MM(SOH) * CDOT(SOH) / Cp_mix & 
	             + H(SH2) * JAC(SH2,SO) - CP(SO) * MM(SH2) * CDOT(SH2) / Cp_mix & 
	             + H(SH2O) * JAC(SH2O,SO) - CP(SO) * MM(SH2O) * CDOT(SH2O) / Cp_mix & 
	             + H(SHO2) * JAC(SHO2,SO) - CP(SO) * MM(SHO2) * CDOT(SHO2) / Cp_mix & 
	             + H(SH2O2) * JAC(SH2O2,SO) - CP(SO) * MM(SH2O2) * CDOT(SH2O2) / Cp_mix) / Cp_mix
      JAC(ST,SOH) = -( H(SN2) * JAC(SN2,SOH) - CP(SOH) * MM(SN2) * CDOT(SN2) / Cp_mix & 
	             + H(SH) * JAC(SH,SOH) - CP(SOH) * MM(SH) * CDOT(SH) / Cp_mix & 
	             + H(SO2) * JAC(SO2,SOH) - CP(SOH) * MM(SO2) * CDOT(SO2) / Cp_mix & 
	             + H(SO) * JAC(SO,SOH) - CP(SOH) * MM(SO) * CDOT(SO) / Cp_mix & 
	             + H(SOH) * JAC(SOH,SOH) - CP(SOH) * MM(SOH) * CDOT(SOH) / Cp_mix & 
	             + H(SH2) * JAC(SH2,SOH) - CP(SOH) * MM(SH2) * CDOT(SH2) / Cp_mix & 
	             + H(SH2O) * JAC(SH2O,SOH) - CP(SOH) * MM(SH2O) * CDOT(SH2O) / Cp_mix & 
	             + H(SHO2) * JAC(SHO2,SOH) - CP(SOH) * MM(SHO2) * CDOT(SHO2) / Cp_mix & 
	             + H(SH2O2) * JAC(SH2O2,SOH) - CP(SOH) * MM(SH2O2) * CDOT(SH2O2) / Cp_mix) / Cp_mix
      JAC(ST,SH2) = -( H(SN2) * JAC(SN2,SH2) - CP(SH2) * MM(SN2) * CDOT(SN2) / Cp_mix & 
	             + H(SH) * JAC(SH,SH2) - CP(SH2) * MM(SH) * CDOT(SH) / Cp_mix & 
	             + H(SO2) * JAC(SO2,SH2) - CP(SH2) * MM(SO2) * CDOT(SO2) / Cp_mix & 
	             + H(SO) * JAC(SO,SH2) - CP(SH2) * MM(SO) * CDOT(SO) / Cp_mix & 
	             + H(SOH) * JAC(SOH,SH2) - CP(SH2) * MM(SOH) * CDOT(SOH) / Cp_mix & 
	             + H(SH2) * JAC(SH2,SH2) - CP(SH2) * MM(SH2) * CDOT(SH2) / Cp_mix & 
	             + H(SH2O) * JAC(SH2O,SH2) - CP(SH2) * MM(SH2O) * CDOT(SH2O) / Cp_mix & 
	             + H(SHO2) * JAC(SHO2,SH2) - CP(SH2) * MM(SHO2) * CDOT(SHO2) / Cp_mix & 
	             + H(SH2O2) * JAC(SH2O2,SH2) - CP(SH2) * MM(SH2O2) * CDOT(SH2O2) / Cp_mix) / Cp_mix
      JAC(ST,SH2O) = -( H(SN2) * JAC(SN2,SH2O) - CP(SH2O) * MM(SN2) * CDOT(SN2) / Cp_mix & 
	             + H(SH) * JAC(SH,SH2O) - CP(SH2O) * MM(SH) * CDOT(SH) / Cp_mix & 
	             + H(SO2) * JAC(SO2,SH2O) - CP(SH2O) * MM(SO2) * CDOT(SO2) / Cp_mix & 
	             + H(SO) * JAC(SO,SH2O) - CP(SH2O) * MM(SO) * CDOT(SO) / Cp_mix & 
	             + H(SOH) * JAC(SOH,SH2O) - CP(SH2O) * MM(SOH) * CDOT(SOH) / Cp_mix & 
	             + H(SH2) * JAC(SH2,SH2O) - CP(SH2O) * MM(SH2) * CDOT(SH2) / Cp_mix & 
	             + H(SH2O) * JAC(SH2O,SH2O) - CP(SH2O) * MM(SH2O) * CDOT(SH2O) / Cp_mix & 
	             + H(SHO2) * JAC(SHO2,SH2O) - CP(SH2O) * MM(SHO2) * CDOT(SHO2) / Cp_mix & 
	             + H(SH2O2) * JAC(SH2O2,SH2O) - CP(SH2O) * MM(SH2O2) * CDOT(SH2O2) / Cp_mix) / Cp_mix
      JAC(ST,SHO2) = -( H(SN2) * JAC(SN2,SHO2) - CP(SHO2) * MM(SN2) * CDOT(SN2) / Cp_mix & 
	             + H(SH) * JAC(SH,SHO2) - CP(SHO2) * MM(SH) * CDOT(SH) / Cp_mix & 
	             + H(SO2) * JAC(SO2,SHO2) - CP(SHO2) * MM(SO2) * CDOT(SO2) / Cp_mix & 
	             + H(SO) * JAC(SO,SHO2) - CP(SHO2) * MM(SO) * CDOT(SO) / Cp_mix & 
	             + H(SOH) * JAC(SOH,SHO2) - CP(SHO2) * MM(SOH) * CDOT(SOH) / Cp_mix & 
	             + H(SH2) * JAC(SH2,SHO2) - CP(SHO2) * MM(SH2) * CDOT(SH2) / Cp_mix & 
	             + H(SH2O) * JAC(SH2O,SHO2) - CP(SHO2) * MM(SH2O) * CDOT(SH2O) / Cp_mix & 
	             + H(SHO2) * JAC(SHO2,SHO2) - CP(SHO2) * MM(SHO2) * CDOT(SHO2) / Cp_mix & 
	             + H(SH2O2) * JAC(SH2O2,SHO2) - CP(SHO2) * MM(SH2O2) * CDOT(SH2O2) / Cp_mix) / Cp_mix
      JAC(ST,SH2O2) = -( H(SN2) * JAC(SN2,SH2O2) - CP(SH2O2) * MM(SN2) * CDOT(SN2) / Cp_mix & 
	             + H(SH) * JAC(SH,SH2O2) - CP(SH2O2) * MM(SH) * CDOT(SH) / Cp_mix & 
	             + H(SO2) * JAC(SO2,SH2O2) - CP(SH2O2) * MM(SO2) * CDOT(SO2) / Cp_mix & 
	             + H(SO) * JAC(SO,SH2O2) - CP(SH2O2) * MM(SO) * CDOT(SO) / Cp_mix & 
	             + H(SOH) * JAC(SOH,SH2O2) - CP(SH2O2) * MM(SOH) * CDOT(SOH) / Cp_mix & 
	             + H(SH2) * JAC(SH2,SH2O2) - CP(SH2O2) * MM(SH2) * CDOT(SH2) / Cp_mix & 
	             + H(SH2O) * JAC(SH2O,SH2O2) - CP(SH2O2) * MM(SH2O) * CDOT(SH2O) / Cp_mix & 
	             + H(SHO2) * JAC(SHO2,SH2O2) - CP(SH2O2) * MM(SHO2) * CDOT(SHO2) / Cp_mix & 
	             + H(SH2O2) * JAC(SH2O2,SH2O2) - CP(SH2O2) * MM(SH2O2) * CDOT(SH2O2) / Cp_mix) / Cp_mix
      JAC(ST,ST) = -( CP(SN2) * MM(SN2) * ( CDOT(SN2) ) + H(SN2)*JAC(SN2,ST)  & 
	             + CP(SH) * MM(SH) * ( CDOT(SH) ) + H(SH)*JAC(SH,ST)  & 
	             + CP(SO2) * MM(SO2) * ( CDOT(SO2) ) + H(SO2)*JAC(SO2,ST)  & 
	             + CP(SO) * MM(SO) * ( CDOT(SO) ) + H(SO)*JAC(SO,ST)  & 
	             + CP(SOH) * MM(SOH) * ( CDOT(SOH) ) + H(SOH)*JAC(SOH,ST)  & 
	             + CP(SH2) * MM(SH2) * ( CDOT(SH2) ) + H(SH2)*JAC(SH2,ST)  & 
	             + CP(SH2O) * MM(SH2O) * ( CDOT(SH2O) ) + H(SH2O)*JAC(SH2O,ST)  & 
	             + CP(SHO2) * MM(SHO2) * ( CDOT(SHO2) ) + H(SHO2)*JAC(SHO2,ST)  & 
	             + CP(SH2O2) * MM(SH2O2) * ( CDOT(SH2O2) ) + H(SH2O2)*JAC(SH2O2,ST) ) / Cp_mix

      END
