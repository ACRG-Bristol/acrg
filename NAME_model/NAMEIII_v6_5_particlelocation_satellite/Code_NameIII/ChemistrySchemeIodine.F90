! Module: NAME Chemistry Scheme Module
! This chemistry scheme converts gaseous elemental iodine to particulate iodine 

Module NameChemistrySchemeModule

! SpeciesListForChemistryScheme
!
! ChemistryScheme---(EQMCON
!                   (CHEMCO
!                   (PH3----PHCALC

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: SpeciesListForChemistryScheme
Public :: ChemistryScheme
Public :: MassReassignment
Public :: DistributeUniformly

!-------------------------------------------------------------------------------------------------------------

! local parameters (MOVED FROM CHEMISTRY SUBROUTINES)
! Note that the molecular masses of dry air and water are specified in the collection
! of Physics parameters (in g/mole). We use these values in the chemistry scheme but
! new parameters are set (in kg/mole) for consistency with the NAME chemistry code.
! We also declare a copy of the Avogadro constant here - for back-compatibility with
! the NAME chemistry code. Note that all these new parameters are REAL(8).
! $$ Change "AVOGAD" in NAME code to "Avogadro"? $$ DBLE ???
! $$ Checking needed for Chemistry Timestep - also would real be better?

Real(8), Parameter :: AVOGAD     = (Avogadro)            ! Copy of Avogadro number.
Real(8), Parameter :: RMM_AIR    = (MoleMassAir)/1.0E3   ! Molecular mass of dry air (kg/mole).
Real(8), Parameter :: RMM_W      = (MoleMassWater)/1.0E3 ! Molecular mass of water (kg/mole).
Real(8), Parameter :: RUGC       = 8.314                 ! Universal gas constant (J K_1 MOL_1)
Real(8), Parameter :: P_REF      = 101325.0              ! Reference pressure (Pa). $$ Also in physics module.
Real(8), Parameter :: MINCONC    = 1.7E-21               ! Minimum concentration value equivalent to
                                                         ! 1 molecule/cm3 (moles/litre).
                                                         ! $$ Calculate direct from Avogadro?
Integer, Parameter :: ChemistryTimestep        = 100     ! Chemistry timestep (secs).
Integer, Parameter :: NumIterationsOfChemistry = 8       ! Number of iterations of the chemistry code.

!species list

!IODINE_E_131 
!IODINE_E_132
!IODINE_E_133
!IODINE_E_135
!IODINE_P_131
!IODINE_P_132
!IODINE_P_133
!IODINE_P_135

!where IODINE_E is gaseous elemental iodine and IODINE_P is particulate iodine

!each species goes back to itself


Real(8), Parameter, Dimension(8,8) :: &                ! Array containing distribution rules for mass reassignment
  MassReassignment =                                                                           &
  Reshape(                                                                                     &
    (/                                                                                         &
      1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !1
      0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0, & !2
      0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0, & !3
      0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0, & !4
      1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !5
      0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0, & !6
      0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0, & !7
      0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0 & !8
    /),                                                                                        &
    (/ 8, 8 /)                                                                               &
  )


Logical, Parameter :: DistributeUniformly(8) = (/         &
                                                  .false., & !IODINE_E_131 
                                                  .false., & !IODINE_E_132      
                                                  .false., & !IODINE_E_133           
                                                  .false., & !IODINE_E_135        
                                                  .false., & !IODINE_P_131      
                                                  .false., & !IODINE_P_132            
                                                  .false., & !IODINE_P_133            
                                                  .false.  & !IODINE_P_135             
                                                /)

! DistributeUniformly :: Indicates that this species, if carried on particles, should be distributed uniformly 
!                        over particles.

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine SpeciesListForChemistryScheme(nChemistrySpecieses, ChemistrySpeciesNames, OutputOnly)

  Implicit None
  ! Argument list:
  Integer,                  Intent(Out) :: nChemistrySpecieses
  Character(MaxCharLength), Intent(Out) :: ChemistrySpeciesNames(MaxSpecieses) 
  Logical,                  Intent(Out) :: OutputOnly(MaxSpecieses)
  ! nChemistrySpecieses   ::
  ! ChemistrySpeciesNames ::
  ! OutputOnly            :: Indicates which species are output from the chemistry scheme but are not input. 
                                                                 
  !$$ all species are currently set as output false, as they are both input and output. However the user
  !$$ may wish to output H2O, or an intermediate chemistry species as a diagnostic, in which case that
  !$$ species could be included in the ChemistrySpeciesNames and set as OutputOnly .true.

  nChemistrySpecieses = 8

  If (nChemistrySpecieses > MaxSpecieses) Then
    Call Message('UNEXPECTED FATAL ERROR in SpeciesListForChemistryScheme.', 4)
  End If

  !$$ I haven't added all the "output" species. Also some species which are actually output species may
  !   not be so labelled.

  ChemistrySpeciesNames( 1) = 'IODINE_E_131          '
  ChemistrySpeciesNames( 2) = 'IODINE_E_132          '
  ChemistrySpeciesNames( 3) = 'IODINE_E_133          '
  ChemistrySpeciesNames( 4) = 'IODINE_E_135          '
  ChemistrySpeciesNames( 5) = 'IODINE_P_131          '
  ChemistrySpeciesNames( 6) = 'IODINE_P_132          '
  ChemistrySpeciesNames( 7) = 'IODINE_P_133          '
  ChemistrySpeciesNames( 8) = 'IODINE_P_135          '

  OutputOnly(1:8) = .false.


End Subroutine SpeciesListForChemistryScheme

!-------------------------------------------------------------------------------------------------------------

Subroutine ChemistryScheme(Conc, IDELT,                                         &
                           ZQFE_R, T_R, AIR_DEN_R, SPECH_R, CLDWAT_R, CLDICE_R, &
                           CF_R, CLD_COL_R, PAMBIENT_R, ZI_R, ZENITH_R)


! This is a simple scheme to convert gaseous Iodine to particulate
! Met input are not needed but left in for consistancy with other chemistry schemes

! All chemistry is done as REAL*8.
! Species and chemistry fields are now stored and read in as double precision.
! Met data are read in as reals (denoted by _R) and converted to double precision.
! Species concentrations are passed as molecules/cm3.

!$$ chem timestep is set to 100 sec, and needs to be divisible into idelt


  Implicit None
  ! Argument list:
  Real(8),   Intent(InOut) :: Conc(MaxSpecieses) ! Gridbox concs of various species.
  Integer,   Intent(In)    :: IDELT       ! Advection timestep.
                                        ! Met data (evaluated at gridbox centre) ...
  Real(Std), Intent(In)    :: ZQFE_R      !] ... height
  Real(Std), Intent(In)    :: T_R         !] ... ambient air temperature
  Real(Std), Intent(In)    :: AIR_DEN_R   !] ... molecular density of air (in molec/cm3)
  Real(Std), Intent(In)    :: SPECH_R     !] ... specific humidity
  Real(Std), Intent(In)    :: CLDWAT_R    !] ... cloud liquid water content
  Real(Std), Intent(In)    :: CLDICE_R    !] ... cloud ice content
  Real(Std), Intent(In)    :: CF_R        !] ... total cloud fraction (dyn+conv)
  Real(Std), Intent(In)    :: CLD_COL_R   !] ... max cloud fraction over cloud column
  Real(Std), Intent(In)    :: PAMBIENT_R  !] ... ambient air pressure
  Real(Std), Intent(In)    :: ZI_R        !] ... boundary layer depth
  Real(Std), Intent(In)    :: ZENITH_R    !] ... solar zenith angle


! met input

      REAL(8) :: ZQFE
      REAL(8) :: T
      REAL(8) :: AIR_DEN
      REAL(8) :: ZENITH
      REAL(8) :: SPECH
      REAL(8) :: CLDWAT
      REAL(8) :: CLDICE
      REAL(8) :: CF
      REAL(8) :: CLD_COL
      REAL(8) :: PAMBIENT
      REAL(8) :: ZI


! Chemistry species input

      REAL(8) :: IODINE_E_131,IODINE_E_132,IODINE_E_133,IODINE_E_135
      REAL(8) :: IODINE_P_131,IODINE_P_132,IODINE_P_133,IODINE_P_135 

! concentrations are stored in _P values between chem timesteps

      REAL(8) :: IODINE_E_131_P,IODINE_E_132_P,IODINE_E_133_P,IODINE_E_135_P
      REAL(8) :: IODINE_P_131_P,IODINE_P_132_P,IODINE_P_133_P,IODINE_P_135_P 


!Reaction rates: notation - generally the name gives the two species
!that are reacting - see CHEMCO for details

      REAL(8) :: RC_IODINE_E_131,RC_IODINE_E_132,RC_IODINE_E_133,RC_IODINE_E_135


! Misc

      REAL(8) ::  P          !Product
      REAL(8) ::  L          !Loss


      INTEGER ::  NIT,I,J    !iteration loop variable
      INTEGER ::  ICHEMT     !chemistry timestep
      INTEGER ::  ILOOP      !chem loop variable
      INTEGER ::  BLFLAG     !boundary layer flag



! initialise species concentrations

!IODINE_E_131 
!IODINE_E_132
!IODINE_E_133
!IODINE_E_135
!IODINE_P_131
!IODINE_P_132
!IODINE_P_133
!IODINE_P_135

      IODINE_E_131       = CONC(1)
      IODINE_E_132       = CONC(2)
      IODINE_E_133       = CONC(3)
      IODINE_E_135       = CONC(4)
      IODINE_P_131       = CONC(5)
      IODINE_P_132       = CONC(6)
      IODINE_P_133       = CONC(7)
      IODINE_P_135       = CONC(8)


      CALL CHEMCO(RC_IODINE_E_131,RC_IODINE_E_132,RC_IODINE_E_133,RC_IODINE_E_135)



! loop through the idelt timestep in chem timesteps

      ICHEMT = ChemistryTimestep
      ILOOP  = IDELT/ICHEMT

      ! $$ This test may not be in the best place. Also perhaps should only insist on imposing an upper limit
      ! on the chemistry time step - then this test wouldn't be needed.
      If (ILoop * IChemT /= IDelT) Then
        Call Message(                                                                  &
               'FATAL ERROR: The main time step (sync time) must be a multiple of ' // &
               Trim(Int2Char(ChemistryTimestep))                                    // &
               ' for runs with chemistry',                                             &
               3                                                                       &
             )
      End If

      DO J=1,ILOOP

! Fix the initial value (as species_P) for use later

        IODINE_E_131_P=IODINE_E_131 
        IODINE_E_132_P=IODINE_E_132
        IODINE_E_133_P=IODINE_E_133
        IODINE_E_135_P=IODINE_E_135
        IODINE_P_131_P=IODINE_P_131
        IODINE_P_132_P=IODINE_P_132
        IODINE_P_133_P=IODINE_P_133
        IODINE_P_135_P=IODINE_P_135


! loop over iterations of the chem code for each chem timestep

        NIT = NumIterationsOfChemistry
        DO I=1,NIT


! DO GASEOUS CHEMISTRY (nb: still inside the iteration loop)
! ********************

! loss of gaseous iodine species to particlulate iodine

! IODINE_E_131

          P=0.0
          L=RC_IODINE_E_131
          IODINE_E_131=(IODINE_E_131_P + (ICHEMT*P))/(1.0 + (ICHEMT*L)) 

! IODINE_E_132

          P=0.0
          L=RC_IODINE_E_132 
          IODINE_E_132=(IODINE_E_132_P + (ICHEMT*P))/(1.0 + (ICHEMT*L)) 

! IODINE_E_133

          P=0.0
          L=RC_IODINE_E_133 
          IODINE_E_133=(IODINE_E_133_P + (ICHEMT*P))/(1.0 + (ICHEMT*L)) 

! IODINE_E_135

          P=0.0
          L=RC_IODINE_E_135 
          IODINE_E_135=(IODINE_E_135_P + (ICHEMT*P))/(1.0 + (ICHEMT*L)) 

! production of particulate iodine species from gaseous iodine

! IODINE_P_131

          P=RC_IODINE_E_131 * IODINE_E_131
          L=0.0
          IODINE_P_131=(IODINE_P_131_P + (ICHEMT*P))/(1.0 + (ICHEMT*L)) 

! IODINE_P_132

          P=RC_IODINE_E_132 * IODINE_E_132
          L=0.0
          IODINE_P_132=(IODINE_P_132_P + (ICHEMT*P))/(1.0 + (ICHEMT*L)) 
! IODINE_P_133

          P=RC_IODINE_E_133 * IODINE_E_133
          L=0.0
          IODINE_P_133=(IODINE_P_133_P + (ICHEMT*P))/(1.0 + (ICHEMT*L)) 
! IODINE_P_135

          P=RC_IODINE_E_135 * IODINE_E_135
          L=0.0
          IODINE_P_135=(IODINE_P_135_P + (ICHEMT*P))/(1.0 + (ICHEMT*L)) 


        ENDDO      ! end of iteration loop

      ENDDO    ! END OF CHEM TIMESTEP LOOP

! Chem species 
      CONC(1)  = IODINE_E_131
      CONC(2)  = IODINE_E_132 
      CONC(3)  = IODINE_E_133
      CONC(4)  = IODINE_E_135  
      CONC(5)  = IODINE_P_131
      CONC(6)  = IODINE_P_132
      CONC(7)  = IODINE_P_133  
      CONC(8)  = IODINE_P_135

End Subroutine ChemistryScheme

!-------------------------------------------------------------------------------------------------------------


Subroutine CHEMCO(RC_IODINE_E_131,RC_IODINE_E_132,RC_IODINE_E_133,RC_IODINE_E_135)

  Implicit None

  ! Argument list:
  REAL(8),    Intent(Out) :: RC_IODINE_E_131        !OUT: rate of conversion to particulate
  REAL(8),    Intent(Out) :: RC_IODINE_E_132        !OUT: rate of conversion to particulate
  REAL(8),    Intent(Out) :: RC_IODINE_E_133        !OUT: rate of conversion to particulate
  REAL(8),    Intent(Out) :: RC_IODINE_E_135        !OUT: rate of conversion to particulate


! GAS PHASE REACTIONS:

! Gas to particulate conversion rate: assuming 17 days (Uematsu et al 1988) as an approximation
! for the complete conversion from gas to particle, gives a rate of 1/146880 seconds = 6.8083 x 10-7

! (Uematsu, M et al Aerosol residence times and iodine gas/particle conversion over the North Pacific
! as determined from Chernobyl radioactivity. Geochemical Journal, Vol. 22, pp. 157 to 163, 1988.)
! IODINE_E_131 converts to IODINE_P_131

      RC_IODINE_E_131=6.8083E-7


! IODINE_E_132 converts to IODINE_P_132

      RC_IODINE_E_132=6.8083E-7


! IODINE_E_133 converts to IODINE_P_133

      RC_IODINE_E_133=6.8083E-7


! IODINE_E_135 converts to IODINE_P_135

      RC_IODINE_E_135=6.8083E-7



End Subroutine CHEMCO

!-------------------------------------------------------------------------------------------------------------


End Module NameChemistrySchemeModule
