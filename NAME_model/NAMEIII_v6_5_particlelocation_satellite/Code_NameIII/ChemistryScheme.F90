! Module: NAME Chemistry Scheme Module

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
Real(8), Parameter :: RUGC       = 8.314                 ! Universal gas constant (J K-1 MOL-1)
Real(8), Parameter :: P_REF      = 101325.0              ! Reference pressure (Pa). $$ Also in physics module.
Real(8), Parameter :: MINCONC    = 1.7E-21               ! Minimum concentration value equivalent to
                                                         ! 1 molecule/cm3 (moles/litre).
                                                         ! $$ Calculate direct from Avogadro?
Integer, Parameter :: ChemistryTimestep        = 100     ! Chemistry timestep (secs).
Integer, Parameter :: NumIterationsOfChemistry = 8       ! Number of iterations of the chemistry code.
Real(8), Parameter, Dimension(42,42) :: &                ! Array containing distribution rules for mass reassignment
  MassReassignment =                                                                           &
  Reshape(                                                                                     &
    (/                                                                                         &
      1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !1
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !2
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !3
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !4
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !5
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !6
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !7
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !8
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !9
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !10
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !11
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !12
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !13
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !14
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !15
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0, & !16
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0, & !17
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0, & !18
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0, & !19
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0, & !20
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0,0.0, & !21
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !22
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !23
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !24
      0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !25
      0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !26
      0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !27
      0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !28
      0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !29
      0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !30
      0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !31
      0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !32
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !33
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !34
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !35
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !36
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !37
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !38
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !39
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !40
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !41
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,         &
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, & !42
      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0          &
    /),                                                                                        &
    (/ 42, 42 /)                                                                               &
  )

  ! mass reassignment rules order
  ! SO2,NH3,NO,(NH4)2SO4,SO4,NO2,NO3,N2O5,HNO3,NAER,NH4NO3,O3,
  ! CO,HCHO,C2H4,C3H6,C5H8,OXYL,TOLUEN,BD,CH3CHO,PAN,HONO,PM10,H2S,
  ! C5H8BIO,C10H16BIO,SOA,SOAP,AROPREC,NAMEARO,PM25,H2O2
  ! OH,HO2,CH3OOH,MVK,ISOPOOH,MVKOOH,MGLYOX,GLYOX,MEMALD
  ! Of these species, the last line above OH - MEMALD are never on particles
  ! and so the massreassignment number will never be used. ($$ this is now guidance only - code is more flexible)
  ! O3 and H2O2 can sometimes be on particles, sometimes not.

Logical, Parameter :: DistributeUniformly(42) = (/         &
                                                  .false., & ! SULPHUR-DIOXIDE
                                                  .false., & ! AMMONIA        
                                                  .false., & ! NO             
                                                  .false., & ! NH42SO4        
                                                  .false., & ! SULPHATE       
                                                  .false., & ! NO2            
                                                  .false., & ! NO3            
                                                  .false., & ! N2O5           
                                                  .false., & ! HNO3           
                                                  .false., & ! NAER           
                                                  .false., & ! NH4NO3         
                                                  .true. , & ! O3             
                                                  .false., & ! CO             
                                                  .false., & ! HCHO           
                                                  .false., & ! C2H4           
                                                  .false., & ! C3H6           
                                                  .false., & ! C5H8           
                                                  .false., & ! OXYL           
                                                  .false., & ! TOLUEN         
                                                  .false., & ! BD             
                                                  .false., & ! CH3CHO         
                                                  .false., & ! PAN            
                                                  .false., & ! HONO           
                                                  .false., & ! PM10           
                                                  .false., & ! H2S            
                                                  .false., & ! C5H8BIO        
                                                  .false., & ! C10H16BIO      
                                                  .false., & ! SOA            
                                                  .false., & ! SOAP           
                                                  .false., & ! AROPREC        
                                                  .false., & ! NAMEARO        
                                                  .false., & ! PM25           
                                                  .true. , & ! H2O2           
                                                  .false., & ! OH             
                                                  .false., & ! HO2            
                                                  .false., & ! CH3OOH         
                                                  .false., & ! MVK            
                                                  .false., & ! ISOPOOH        
                                                  .false., & ! MVKOOH         
                                                  .false., & ! MGLYOX         
                                                  .false., & ! GLYOX          
                                                  .false.  & ! MEMALD         
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

  nChemistrySpecieses = 42

  If (nChemistrySpecieses > MaxSpecieses) Then
    Call Message('UNEXPECTED FATAL ERROR in SpeciesListForChemistryScheme.', 4)
  End If

  !$$ I haven't added all the "output" species. Also some species which are actually output species may
  !   not be so labelled.

  ChemistrySpeciesNames( 1) = 'SULPHUR-DIOXIDE'
  ChemistrySpeciesNames( 2) = 'AMMONIA        '
  ChemistrySpeciesNames( 3) = 'NO             '
  ChemistrySpeciesNames( 4) = 'NH42SO4        '
  ChemistrySpeciesNames( 5) = 'SULPHATE       '
  ChemistrySpeciesNames( 6) = 'NO2            '
  ChemistrySpeciesNames( 7) = 'NO3            '
  ChemistrySpeciesNames( 8) = 'N2O5           '
  ChemistrySpeciesNames( 9) = 'HNO3           '
  ChemistrySpeciesNames(10) = 'NAER           '
  ChemistrySpeciesNames(11) = 'NH4NO3         '
  ChemistrySpeciesNames(12) = 'O3             '
  ChemistrySpeciesNames(13) = 'CO             '
  ChemistrySpeciesNames(14) = 'HCHO           '
  ChemistrySpeciesNames(15) = 'C2H4           '
  ChemistrySpeciesNames(16) = 'C3H6           '
  ChemistrySpeciesNames(17) = 'C5H8           '
  ChemistrySpeciesNames(18) = 'OXYL           '
  ChemistrySpeciesNames(19) = 'TOLUEN         '
  ChemistrySpeciesNames(20) = 'BD             '
  ChemistrySpeciesNames(21) = 'CH3CHO         '
  ChemistrySpeciesNames(22) = 'PAN            '
  ChemistrySpeciesNames(23) = 'HONO           '
  ChemistrySpeciesNames(24) = 'PM10           '
  ChemistrySpeciesNames(25) = 'H2S            '
  ChemistrySpeciesNames(26) = 'C5H8BIO        '
  ChemistrySpeciesNames(27) = 'C10H16BIO      '
  ChemistrySpeciesNames(28) = 'SOA            '
  ChemistrySpeciesNames(29) = 'SOAP           '
  ChemistrySpeciesNames(30) = 'AROPREC        '
  ChemistrySpeciesNames(31) = 'NAMEARO        '
  ChemistrySpeciesNames(32) = 'PM25           '
  ChemistrySpeciesNames(33) = 'H2O2           '
  ChemistrySpeciesNames(34) = 'OH             '
  ChemistrySpeciesNames(35) = 'HO2            '
  ChemistrySpeciesNames(36) = 'CH3OOH         '
  ChemistrySpeciesNames(37) = 'MVK            '
  ChemistrySpeciesNames(38) = 'ISOPOOH        '
  ChemistrySpeciesNames(39) = 'MVKOOH         '
  ChemistrySpeciesNames(40) = 'MGLYOX         '
  ChemistrySpeciesNames(41) = 'GLYOX          '
  ChemistrySpeciesNames(42) = 'MEMALD         '

  OutputOnly(1:42) = .false.

End Subroutine SpeciesListForChemistryScheme

!-------------------------------------------------------------------------------------------------------------

Subroutine ChemistryScheme(Conc, IDELT,                                         &
                           ZQFE_R, T_R, AIR_DEN_R, SPECH_R, CLDWAT_R, CLDICE_R, &
                           CF_R, CLD_COL_R, PAMBIENT_R, ZI_R, ZENITH_R)


! Performs all the gaseous and aqueous phase chemistry for a chemistry grid box
! (reactions based on STOCHEM).

! This is a Fortran 90 version of the subroutine CHEM in NAME V8.17.

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


! Arrays - local

      REAL(8) :: LAQ(10)     !LOSSES ARRAY
      REAL(8) :: A(17)       !photolysis parameter
      REAL(8) :: B(17)       !photolysis parameter
      REAL(8) :: DJ(17)      !phololysis rate

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

      REAL(8) ::  SO2,NH3,NO,NH42SO4,SO4,NO2,NO3,N2O5
      REAL(8) ::  HNO3,NAER,NH4NO3,CO,HCHO,C2H4
      REAL(8) ::  C3H6,C5H8,OXYL,TOLUEN,BD,CH3CHO,PAN,HONO,H2S
      REAL(8) ::  C5H8BIO,C10H16BIO,SOA,SOAP,AROPREC,NAMEARO
      REAL(8) ::  O3,H2O2
      Real(8) :: OH      ! Hydroxide.
      Real(8) :: HO2     ! Hydroperoxy.
      Real(8) :: CH3OOH  ! Methylhydroperoxide.
      Real(8) :: MVK     !
      Real(8) :: ISOPOOH !
      Real(8) :: MVKOOH  !
      Real(8) :: MGLYOX  !
      Real(8) :: GLYOX   !
      Real(8) :: MEMALD  !
      Real(8) :: H2O     !

! Chem species calculated only in CHEM

      REAL(8) ::  OP,OD,CH3COO2,CH2OO,CH3CHX,OXYL1,MEMALD1
      REAL(8) ::  RO2IP1,RO2IP2,TOLP1,BDPEROXY,CH3O2
      REAL(8) ::  CH2O2C,HO2NO2


! Chem species with values set/calced at the start of chem

      REAL(8) ::  CH4,H2,O2

! concentrations are stored in _P values between chem timesteps

      REAL(8) ::  SO2_P,SO4_P,NH3_P,NH42SO4_P,NAER_P,HNO3_P
      REAL(8) ::  N2O5_P,NO3_P,NO2_P,NO_P,H2O2_P,O3_P
      REAL(8) ::  NH4NO3_P,CO_P,HCHO_P,C2H4_P
      REAL(8) ::  CH4_P,CH3O2_P,CH2O2C_P,HO2NO2_P,CH3OOH_P
      REAL(8) ::  H2_P,CH3CHO_P,CH3COO2_P,PAN_P,C3H6_P
      REAL(8) ::  OXYL_P,CH2OO_P,MGLYOX_P,CH3CHX_P,GLYOX_P
      REAL(8) ::  OXYL1_P,MEMALD_P,MEMALD1_P,C5H8_P,RO2IP1_P
      REAL(8) ::  MVK_P,RO2IP2_P,ISOPOOH_P,MVKOOH_P,TOLUEN_P
      REAL(8) ::  TOLP1_P,BD_P,BDPEROXY_P,HONO_P,H2S_P,HO2_P
      REAL(8) ::  C5H8BIO_P,C10H16BIO_P,SOA_P,SOAP_P,AROPREC_P,NAMEARO_P


!used to calc the eqm values at the end of chem, after loop

      REAL(8) ::  NH4NO3_EQM,NH3_EQM,HNO3_EQM
!      REAL(8) ::  EQM_NH3HNO3               !NH3 + HNO3 EQUILIBRIUM


!cloud total of species in moles/litre of air

      REAL(8) ::  CT_SO2,CT_H2O2,CT_HNO3,CT_NH42SO4
      REAL(8) ::  CT_NH3

!These values are set equal to ct's within each timestep and used
!to calculate the new value (nb moles/litre of air)

      REAL(8) ::  PREV_SO2,PREV_H2O2,PREV_NH42SO4
      REAL(8) ::  PREV_NH3,PREV_SO4,PREV_HNO3

!These values are set equal to ct's at the start and used at the
!end to mix the aq and gas phases together(nb moles/litre of air)

      REAL(8) ::  INIT_SO2,INIT_H2O2,INIT_NH42SO4
      REAL(8) ::  INIT_NH3,INIT_SO4,INIT_HNO3

!Aqueous concentrations

      REAL(8) ::  AQ_O3,AQ_HNO3,AQ_H2O2,AQ_SO2,AQ_NH3
      REAL(8) ::  AQ_CO2,AQ_NO3,AQ_HSO3,AQ_SO3,AQ_NH42SO4
      REAL(8) ::  AQ_HCO3,AQ_SO4


!Equilibrium coefficients

      REAL(8) ::  KE_HNO3,KE_SO2,KE_HSO3,KE_NH3,KE_CO2,KE_H2O


!Henry's law coefficients

      REAL(8) ::  KH_O3,KH_HNO3,KH_H2O2,KH_SO2,KH_NH3,KH_CO2


!Reaction rates: notation - generally the name gives the two species
!that are reacting - see CHEMCO for details

      REAL(8) :: RC_HO2_HO2,RC_SO2_OH,RC_HSO3_H2O2,RC_HSO3_O3
      REAL(8) :: RC_SO3_O3,RC_NO_O3,RC_NO2_O3,RC_NO2_NO3
      REAL(8) :: RC_NO2_OH,RC_NAER1,RC_NAER2,RC_OH_HNO3,RC_NO3_NO3
      REAL(8) :: RC_NO3_HO2,RC2_NO3_HO2,RC_NO_NO3,RC2_NO2_NO3
      REAL(8) :: RC_N2O5,RC_OP_O2,RC_OP_NO,RC_OD,RC_OD_H2O
      REAL(8) :: RC_O3_C2H4,RC_OH_HCHO,RC_NO3_HCHO,RC_OH_CO
      REAL(8) :: RC_OH_CH4,RC_CH3O2_CH2O2C,RC_CH2O2C_NO
      REAL(8) :: RC_CH3O2_CH3O2,RC_NO_CH3O2,RC_CH3O2_2
      REAL(8) :: RC_SO2_CH3O2,RC_OH_CH3OOH,RC2_OH_CH3OOH
      REAL(8) :: RC_CH3O2_HO2,RC_OH_C2H4,RC_NO3_C2H4,RC_NO2_HO2
      REAL(8) :: RC_HO2NO2,RC_OH_HO2NO2,RC_OH_H2O2,RC_OH_H2
      REAL(8) :: RC_OH_O3,RC_HO2_O3,RC_NO_HO2,RC_HO2_OH
      REAL(8) :: RC_CH3O2_CH3CHX,RC_O3_C3H6,RC_CH3CHX_NO
      REAL(8) :: RC_OH_CH3CHO,RC_NO3_CH3CHO,RC_MGLYOX_OH
      REAL(8) :: RC_PAN,RC_CH3COO2_CH3COO2,RC_CH3O2_CH3COO2
      REAL(8) :: RC_CH3COO2_NO2,RC_CH3COO2_NO,RCA_CH3O2_CH3COO2
      REAL(8) :: RC_OH_PAN,RC_NO3_C3H6,RC2_O3_C3H6,RC_OH_C3H6
      REAL(8) :: RC_OH_OXYL,RC_O3_MVK,RC_O3_C5H8,RC_CH2OO_NO
      REAL(8) :: RC_CH2OO_NO2,RC_CH2OO_SO2,RC_CH2OO_H2O
      REAL(8) :: RC_RO2IP2_NO,RC_MVKOOH_OH,RC_MEMALD1_CH3O2
      REAL(8) :: RC_RO2IP2_CH3O2,RC_MEMALD1_NO,RC_OXYL1_NO
      REAL(8) :: RC_NO_TOLP1,RC_GLYOX_OH,RC_MEMALD_OH
      REAL(8) :: RC_NO3_C5H8,RC_OH_C5H8,RC_RO2IP1_CH3O2
      REAL(8) :: RC_RO2IP1_HO2,RC_RO2IP1_NO,RC_ISOPOOH_OH
      REAL(8) :: RC_OH_MVK,RC_RO2IP2_HO2,RC_OH_TOLUEN
      REAL(8) :: RC_NO2_OP,RC_OH_BD,RC_NO_BDPEROXY,RC_OH_HONO
      REAL(8) :: RC_OH_NO,RC_NO2,RC_OH_H2S,RC_N2O5_2
      REAL(8) :: RC_O3_C10H16BIO,RC_NO3_C10H16BIO,RC_OH_C10H16BIO
      REAL(8) :: RC_OH_SOAP,RC_SCAV,RC2_OH_TOLUEN,RC_OH_AROPREC,KIN,KOUT


! Misc

      REAL(8) ::  F
      REAL(8) ::  FSO2       !extra factor for aq so2
      REAL(8) ::  PH         !PH
      REAL(8) ::  HP         !first guess ph
      REAL(8) ::  REUS       !aq conversion factor
      REAL(8) ::  P          !Product
      REAL(8) ::  L          !Loss
      REAL(8) ::  TV         !phot calc variable
      REAL(8) ::  SEC        !1/COS(ZEN) to be consistant with stochem
      REAL(8) ::  R1,R2,L1,L2,L3 !GAS PHASE CALC COEFFS
      REAL(8) ::  RH         !relative humidity
      REAL(8) ::  BYPROD     !controls amount of anthropogenic SOA precursor created
      REAL(8) ::  VNO        !volume mixing ratio of NO
      REAL(8) ::  REACFRAC   !calculated to represent nox dependancy of anth SOA

! Required for new nitrate aerosol scheme

      REAL(8) ::  TNH3
      REAL(8) ::  THNO3
      REAL(8) ::  Kp_RT2
      REAL(8) ::  CONST


      INTEGER ::  NIT,I,J    !iteration loop variable
      INTEGER ::  ICHEMT     !chemistry timestep
      INTEGER ::  ILOOP      !chem loop variable
      INTEGER ::  BLFLAG     !boundary layer flag

      ZQFE     = DBLE(ZQFE_R)
      T        = DBLE(T_R)
      AIR_DEN  = DBLE(AIR_DEN_R)
      ZENITH   = DBLE(ZENITH_R)
      SPECH    = DBLE(SPECH_R)
      CLDWAT   = DBLE(CLDWAT_R)
      CLDICE   = DBLE(CLDICE_R)
      CF       = DBLE(CF_R)
      CLD_COL  = DBLE(CLD_COL_R)
      PAMBIENT = DBLE(PAMBIENT_R)
      ZI       = DBLE(ZI_R)


! set up sec to equal 1/cos(zen) so be in line with stochem
!$$ Slightly dangerous since cos(zenith) is zero when zenith is pi/2 (i.e. sun is
!$$ on the horizon) - then division by zero could occur.
      SEC=1.0/DCOS(ZENITH)


! set up photolysis parameters

!     DISSOCIATION RATE COEFFICIENTS CALCULATED FOR CLEAR SKY, 340 DU
!     TAKEN FROM STANDARD PHOTOCHEMICAL TRAJECTORY MODEL VERSION:
!     CHMIV (1994).

!      NO2 + hv = NO + OP ie reaction DJ(3) in Dick's code
      A(1)=1.3320E-02
      B(1)=0.4020
!      NO3 + hv = NO            DJ(13)
      A(2)=2.8590E-02
      B(2)=0.1980
!      NO3 + hv = NO2 + OP      DJ(14)
      A(3)=1.9870E-01
      B(3)=0.2080
!      N2O5 + hv = NO2 + NO3    DJ(15)
      A(4)=3.3240E-05
      B(4)=0.5664
!      O3 + hv = OP + O2        DJ(1)
      A(5)=5.8940E-04
      B(5)=0.2450
!      O3 + hv = OD + O2        DJ(2)
      A(6)=1.5100E-04
      B(6)=1.6550
!      H2O2 + hv = OH + OH      DJ(4)
      A(7)=1.4880E-05
      B(7)=0.6740
!      HNO3 + hv = NO2 + OH     DJ(5)
      A(8)=1.6950E-06
      B(8)=0.9660
!      HCHO + hv = HO2 + HO2 + CO   DJ(6)
      A(9)=6.7360E-05
      B(9)=0.7650
!      HCHO + hv = H2 + CO          DJ(7)
      A(10)=8.7440E-05
      B(10)=0.5860
!      CH3OOH + hv = HCHO + HO2 + OH  DJ(16)
      A(11)=2.2700E-05
      B(11)=0.6208
!      HO2NO2 + hv = HO2 + NO2    DJ(10)
      A(12)=1.3470E-04
      B(12)=0.6020
!      CH3CHO + hv = CH3O2 + HO2 + CO  DJ(8)
      A(13)=1.3480E-05
      B(13)=1.0530
!      PAN + hv = CH3COO2 + NO2  DJ(17)
      A(14)=2.2700E-05
      B(14)=0.6208
!      MGLYOX + hv = CH3CHO + HO2 + CO DJ(11)
      A(15)=1.7140E-04
      B(15)=0.3010
!      GLYOX + hv = HCHO + CO  DJ(12)
      A(16)=1.1340E-05
      B(16)=0.2730
!      HONO = OH + NO
! reaction not in stochem - A and B values estimated using
! data in 4th PORG Report 1997
      A(17)=2.135E-03
      B(17)=0.294

! Calculate photolysis rate
! nb this is adjusted for cloud later in the code

      DO I=1,17
      IF(SEC.GT.0.0)THEN

        TV=-B(I)*SEC
        IF(TV.LT.-90.0)THEN
            DJ(I)=1.0E-20
        ELSE
            DJ(I)=A(I)*DEXP(-B(I)*SEC)
        ENDIF
      ELSE
        DJ(I)=1.0E-20
      ENDIF
      ENDDO

!set DJ(17) to equal DJ(1)*0.18 as Dick has

      DJ(17)=DJ(1)*0.18

! set water vapour in molecules/cm3

      H2O=SPECH*AIR_DEN*RMM_AIR/RMM_W

! set molecular oxygen

      O2=0.2095*AIR_DEN

! initialise species concentrations
! note that CONC(24) is PM10 and CONC(32) is PM25 and so not needed in the chemistry

      SO2       = CONC(1)
      NH3       = CONC(2)
      NO        = CONC(3)
      NH42SO4   = CONC(4)
      SO4       = CONC(5)
      NO2       = CONC(6)
      NO3       = CONC(7)
      N2O5      = CONC(8)
      HNO3      = CONC(9)
      NAER      = CONC(10)
      NH4NO3    = CONC(11)
      O3        = CONC(12)
      CO        = CONC(13)
      HCHO      = CONC(14)
      C2H4      = CONC(15)
      C3H6      = CONC(16)
      C5H8      = CONC(17)
      OXYL      = CONC(18)
      TOLUEN    = CONC(19)
      BD        = CONC(20)
      CH3CHO    = CONC(21)
      PAN       = CONC(22)
      HONO      = CONC(23)
      H2S       = CONC(25)
      C5H8BIO   = CONC(26)
      C10H16BIO = CONC(27)
      SOA       = CONC(28)
      SOAP      = CONC(29)
      AROPREC   = CONC(30)
      NAMEARO   = CONC(31)
      H2O2      = CONC(33)
      OH        = Conc(34)
      HO2       = Conc(35)
      CH3OOH    = Conc(36)
      MVK       = Conc(37)
      ISOPOOH   = Conc(38)
      MVKOOH    = Conc(39)
      MGLYOX    = Conc(40)
      GLYOX     = Conc(41)
      MEMALD    = Conc(42)

! combine the two isoprene values into c5h8

      C5H8=C5H8+C5H8BIO

! add in a regional background component for CO and HCHO
! (this will be subtracted again at the end of the subroutine)

      CO=CO+2.9529E12
      HCHO=HCHO+3.6594E10

! nb. OD, OP, OH, HO2 etc are calculated afresh each timestep

      HO2      = 0.0
      OH       = 0.0
      OP       = 0.0
      OD       = 0.0
      TOLP1    = 0.0
      MEMALD1  = 0.0
      OXYL1    = 0.0
      CH2OO    = 0.0
      RO2IP2   = 0.0
      RO2IP1   = 0.0
      CH3COO2  = 0.0
      CH3CHX   = 0.0
      CH3O2    = 0.0
      CH2O2C   = 0.0
      HO2NO2   = 0.0
      BDPEROXY = 0.0

      CH4 = 1.8E-06*AIR_DEN
      H2  = 0.5E-06*AIR_DEN

      NH4NO3_EQM = 0.0
      HNO3_EQM   = 0.0
      NH3_EQM    = 0.0

! calculate the relative humidity
! **using the NAME III physics function CalcRH**

      RH=DBLE(CalcRH(SPECH_R, T_R, PAMBIENT_R))
      RH=MIN(RH,100.0d0)

! Calculate equilibrium and Henry's law coefficients

      CALL EQMCON(T,KE_SO2,KE_HSO3,KE_NH3,KE_CO2,KE_H2O,    &
        KE_HNO3,KH_HNO3,KH_O3,KH_H2O2,KH_SO2,KH_NH3,KH_CO2)

      IF(ZQFE.GT.ZI)THEN
        BLFLAG=0
      ELSE
        BLFLAG=1
      ENDIF

      CALL CHEMCO(T,AIR_DEN,H2O,O2,RC_HO2_HO2,RC_SO2_OH,            &
        RC_HSO3_H2O2,RC_HSO3_O3,RC_SO3_O3,RC_NO_O3,RC_NO2_O3,       &
        RC_NO2_OH,RC_OH_HNO3,RC_NO3_NO3,RC_NO3_HO2,                 &
        RC_NO_NO3,RC2_NO2_NO3,RC_N2O5,RC_OP_O2,RC_OP_NO,RC_OD,      &
!        EQM_NH3HNO3,                                                &
        RC_O3_C2H4,RC_OH_HCHO,RC_NO3_HCHO,RC_OH_CO,                 &
        RC_OH_CH4,RC_CH3O2_CH2O2C,RC_CH2O2C_NO,RC_CH3O2_HO2,        &
        RC_CH3O2_CH3O2,RC_NO_CH3O2,RC_CH3O2_2,RC_SO2_CH3O2,         &
        RC_OH_CH3OOH,RC2_OH_CH3OOH,RC_OH_C2H4,RC_NO3_C2H4,          &
        RC_NO2_HO2,RC_HO2NO2,RC_OH_HO2NO2,RC_OH_H2O2,RC_OH_H2,      &
        RC_HO2_O3,RC_NO_HO2,RC_HO2_OH,RC_NO2_NO3,RC2_NO3_HO2,       &
        RC_OD_H2O,RC_OH_O3,RC_CH3O2_CH3CHX,RC_O3_C3H6,              &
        RC_CH3CHX_NO,RC_OH_CH3CHO,RC_NO3_CH3CHO,RC_MGLYOX_OH,       &
        RC_PAN,RC_CH3COO2_CH3COO2,RC_CH3O2_CH3COO2,RC_CH3COO2_NO2,  &
        RC_CH3COO2_NO,RCA_CH3O2_CH3COO2,RC_OH_PAN,RC_NO3_C3H6,      &
        RC2_O3_C3H6,RC_OH_C3H6,RC_OH_OXYL,RC_O3_MVK,                &
        RC_O3_C5H8,RC_CH2OO_NO,RC_CH2OO_NO2,RC_CH2OO_SO2,           &
        RC_CH2OO_H2O,RC_RO2IP2_NO,RC_MVKOOH_OH,RC_MEMALD1_CH3O2,    &
        RC_RO2IP2_CH3O2,RC_MEMALD1_NO,RC_OXYL1_NO,                  &
        RC_NO_TOLP1,RC_GLYOX_OH,RC_MEMALD_OH,RC_NO3_C5H8,           &
        RC_OH_C5H8,RC_RO2IP1_CH3O2,RC_RO2IP1_HO2,RC_RO2IP1_NO,      &
        RC_ISOPOOH_OH,RC_OH_MVK,RC_RO2IP2_HO2,RC_OH_TOLUEN,         &
        RC_NO2_OP,RC_OH_BD,RC_NO_BDPEROXY,RC_OH_HONO,RC_OH_NO,      &
        RC_NO2,ZI,RH,Kp_RT2,RC_NAER1,RC_NAER2,RC_OH_H2S,BLFLAG,     &
        RC_N2O5_2,RC_O3_C10H16BIO,RC_NO3_C10H16BIO,RC_OH_C10H16BIO, &
        RC_OH_SOAP,RC_SCAV,RC2_OH_TOLUEN,RC_OH_AROPREC,KIN,KOUT)



! adjust the J values if its cloudy - nb the adjustments are
! an average of the reductions for the individual reactions taken
! from the table provided by Dick Derwent

!$$ this was found to reduce the chemistry too much so currently assuming clear skies

!      IF(CLD_COL.GT.0.0)THEN
!        IF(CLD_COL.GT.0.8)THEN
!          DO I=1,17
!            DJ(I)=DJ(I)*0.181
!          ENDDO
!        ELSEIF((CLD_COL.GT.0.6).AND.(CLD_COL.LE.0.8))THEN
!          DO I=1,17
!            DJ(I)=DJ(I)*0.359
!          ENDDO
!        ELSEIF((CLD_COL.GT.0.4).AND.(CLD_COL.LE.0.6))THEN
!          DO I=1,17
!            DJ(I)=DJ(I)*0.552
!          ENDDO
!        ELSEIF((CLD_COL.GT.0.2).AND.(CLD_COL.LE.0.4))THEN
!          DO I=1,17
!            DJ(I)=DJ(I)*0.734
!          ENDDO
!        ELSEIF((CLD_COL.GT.0.0).AND.(CLD_COL.LE.0.2))THEN
!          DO I=1,17
!            DJ(I)=DJ(I)*0.895
!          ENDDO
!        ENDIF
!      ENDIF

! Convert concentrations in molecules/cm3 to moles/litre of air
! to use in calculation of aqueous chemistry - Cloud Total values

      CT_SO2=SO2*1.0E3/AVOGAD
      CT_H2O2=H2O2*1.0E3/AVOGAD
      CT_HNO3=HNO3*1.0E3/AVOGAD
      CT_NH42SO4=NH42SO4*1.0E3/AVOGAD
      CT_NH3=NH3*1.0E3/AVOGAD
      IF(CLDWAT.GT.0.0)THEN
        AQ_SO4=SO4*1.0E3/(AVOGAD*CLDWAT)
      ELSE
        AQ_SO4=0.0
      ENDIF

! now check that concentrations arn't zero: if they are set them
! to a tiny value equivalent to 1 molecule per cm3 so that the
! aqueous phase code will still work.
! 1 molec/cm3 = 1.7E-21 moles/litre

      CT_SO2=MAX(CT_SO2,MINCONC)
      CT_H2O2=MAX(CT_H2O2,MINCONC)
      CT_NH42SO4=MAX(CT_NH42SO4,MINCONC)
      CT_NH3=MAX(CT_NH3,MINCONC)
      AQ_SO4=MAX(AQ_SO4,MINCONC)
      CT_HNO3=MAX(CT_HNO3,MINCONC)

      INIT_SO2=CT_SO2
      INIT_H2O2=CT_H2O2
      INIT_HNO3=CT_HNO3
      INIT_NH42SO4=CT_NH42SO4
      INIT_NH3=CT_NH3
      INIT_SO4=AQ_SO4

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

        SO2_P=SO2
        NH3_P=NH3
        NO_P=NO
        NH42SO4_P=NH42SO4
        SO4_P=SO4
        NO2_P=NO2
        NO3_P=NO3
        N2O5_P=N2O5
        HNO3_P=HNO3
        NAER_P=NAER
        NH4NO3_P=NH4NO3
        CO_P=CO
        HCHO_P=HCHO
        C2H4_P=C2H4
        C3H6_P=C3H6
        C5H8_P=C5H8
        OXYL_P=OXYL
        TOLUEN_P=TOLUEN
        BD_P=BD
        CH3CHO_P=CH3CHO
        PAN_P=PAN
        HONO_P=HONO
        H2O2_P=H2O2
        O3_P=O3
        CH3OOH_P=CH3OOH
        MVK_P=MVK
        ISOPOOH_P=ISOPOOH
        MVKOOH_P=MVKOOH
        MGLYOX_P=MGLYOX
        GLYOX_P=GLYOX
        MEMALD_P=MEMALD
        TOLP1_P=TOLP1
        MEMALD1_P=MEMALD1
        OXYL1_P=OXYL1
        CH2OO_P=CH2OO
        RO2IP2_P=RO2IP2
        RO2IP1_P=RO2IP1
        CH3COO2_P=CH3COO2
        CH3CHX_P=CH3CHX
        CH3O2_P=CH3O2
        CH2O2C_P=CH2O2C
        HO2NO2_P=HO2NO2
        BDPEROXY_P=BDPEROXY
        CH4_P=CH4
        H2_P=H2
        H2S_P=H2S
        HO2_P=HO2
        C10H16BIO_P=C10H16BIO
        SOA_P=SOA
        SOAP_P=SOAP
        AROPREC_P=AROPREC
        NAMEARO_P=NAMEARO

! if we have cloud and cldwat then can do aq chem:

        IF(CF.GT.0.0.AND.CLDWAT.GT.0.0)THEN


! also set prev_ variables equal to the initial values so that they
! can be used later

          PREV_SO2=CT_SO2
          PREV_H2O2=CT_H2O2
          PREV_NH42SO4=CT_NH42SO4
          PREV_NH3=CT_NH3
          PREV_SO4=AQ_SO4
          PREV_HNO3=CT_HNO3

! Calculate ph

          CALL PH3(KE_SO2,KE_HSO3,KE_NH3,KE_CO2,KE_H2O,KE_HNO3, &
            KH_O3,KH_H2O2,KH_SO2,KH_NH3,KH_CO2,KH_HNO3,CLDWAT,  &
            AIR_DEN,T,CT_SO2,CT_NH3,AQ_SO4,CT_HNO3,HP)

!nb do not put in HP** as when you complile code for PC it comments
!out lines containing HP**

!02/06/09 Allowing pH to be calculated if it is more acidic than 5.8,
!but not allowing it to go more basic as over-production of sulphate
!results. Poor modelling of basic values is thought to be due to poorly
!resolved ammonia emissions.

          IF (HP.LT.1.58E-6) THEN
            HP=1.58E-6
          ENDIF

          FSO2=(1.0+KE_SO2/HP+(KE_SO2*KE_HSO3/(HP*HP)))
          F=1.0
          PH=-DLOG10(HP)

        ENDIF  !setting up the aq chem


! loop over iterations of the chem code (gas & aq) for each chem timestep

        NIT = NumIterationsOfChemistry
        DO I=1,NIT

! if we have cloud and cldwat then can do aq chem:

          IF(CF.GT.0.0.AND.CLDWAT.GT.0.0)THEN

! obtain dissolved species concentrations in moles/litre

            REUS=RUGC*T*1.0E6/(AVOGAD*P_REF)

            AQ_O3=KH_O3*O3*REUS

!nb nothing happens to hno3 in the aq phase... just used for ph

            AQ_HNO3=CT_HNO3/CLDWAT                  !all dissolved

            AQ_H2O2=CT_H2O2/(CLDWAT+P_REF/(KH_H2O2*RUGC*T*1.0E3))

            AQ_SO2=CT_SO2/(CLDWAT+P_REF/(KH_SO2*RUGC*T*FSO2))

            AQ_NH3=CT_NH3/CLDWAT

            AQ_CO2=KH_CO2*360.0E-6*AIR_DEN*REUS

            AQ_NO3=AQ_HNO3

            AQ_HSO3=(AQ_SO2/FSO2)*KE_SO2/HP

            AQ_SO3=AQ_HSO3*KE_HSO3/HP

            AQ_NH42SO4=AQ_NH3/(1.0+KE_H2O/(HP*KE_NH3))

            AQ_HCO3=AQ_CO2/(1.0+HP/KE_CO2)

! Calculate the production of so4 and the loss of h2o2

! Loss rates:

!  NOT USED (ie in STOCHEM)    LAQ(1) - loss of gaseous SO2,
!   LAQ(2) - gain of SA,
!   LAQ(3) - loss of in-cloud total H2O2,
!   LAQ(4) - loss of gaseous H2O2,
!   LAQ(5) - loss of in-cloud total SO2.
!   LAQ(6) - loss of in-cloud total NH3,
!   LAQ(7) - loss of in-cloud SO4,
!   LAQ(8) - gain of in-cloud total (NH4)2SO4,
!   LAQ(9) - loss of SA.(in molecules/cm3 to use in the gas phase)
!   LAQ(10) - gain of Ammonium Sulphate.

            IF(PREV_NH3.GT.1.0E-20.AND.PREV_SO4.GT.1.0E-20) THEN

! LAQ(6) - loss of in-cloud total NH3
! this is a fractional loss per second
! CLDWAT needs to be in litres of water / litres of air

              LAQ(6) = 2.0*MIN((AQ_NH42SO4/2.0),AQ_SO4)*CLDWAT/(ICHEMT*CT_NH3)

! LAQ(7) - loss of in-cloud SO4
! this is a fractional loss per second

              LAQ(7) = MIN((AQ_NH42SO4/2.0),AQ_SO4)/(ICHEMT*AQ_SO4)

! LAQ(8) - gain of in-cloud total (NH4)2SO4
! this is used for a production and is in moles/litre air/second

              LAQ(8) = LAQ(6)*CT_NH3/2.0

! LAQ(9) - loss of SA.!
! this is fractional loss of SA per second * cloud fraction
! THIS IS BECAUSE ITS USED IN THE GAS PHASE...

              LAQ(9) = LAQ(7)*AQ_SO4*CLDWAT*AVOGAD*1.0E-3*CF/SO4

! LAQ(10) - gain of Ammonium Sulphate.
! this is gain in molecules/cm3/second in cloud

              LAQ(10)= LAQ(9)*SO4


            ELSE
              LAQ(6) = 0.0
              LAQ(7) = 0.0
              LAQ(8) = 0.0
              LAQ(9) = 0.0
              LAQ(10) = 0.0
            ENDIF

! SO4(Aq)  (AQ_SO4)
! now work out production of SO4 occuring

            P= AQ_HSO3*AQ_H2O2*RC_HSO3_H2O2*HP/(HP+0.1)             &
               + AQ_HSO3*AQ_O3*RC_HSO3_O3 + AQ_SO3*AQ_O3*RC_SO3_O3

            L = LAQ(7)

! new value for aqueous so4

            AQ_SO4 = (PREV_SO4+ICHEMT*P*CF)/(1.0+ICHEMT*L*CF)

! LAQ(2) - gain of SA,
            LAQ(2) = P*AVOGAD*CLDWAT*1.0E-3*CF

! LAQ(3) - loss of in-cloud total H2O2,

            LAQ(3) = CLDWAT*F*(AQ_HSO3*AQ_H2O2*RC_HSO3_H2O2*HP/(HP+0.1))/CT_H2O2

! LAQ(4) - loss of gaseous H2O2,

            LAQ(4) = LAQ(3)*CT_H2O2*AVOGAD*1.0E-3*CF/H2O2

! LAQ(5) - loss of in-cloud total SO2

            LAQ(5) = CLDWAT*P/CT_SO2

!NB we don't model loss of O3 in the aq phase (neither does STOCHEM)

! Total cloud peroxide  CT_H2O2
            P=0.0
            L=LAQ(3)
            CT_H2O2 = (PREV_H2O2+ICHEMT*P)/(1.0+ICHEMT*L*CF)

! Total cloud SO2  CT_SO2
            P=0.0
            L=LAQ(5)
            CT_SO2 = (PREV_SO2+ICHEMT*P)/(1.0+ICHEMT*L*CF)

! Total cloud NH3  CT_NH3
            P=0.0
            L=LAQ(6)
            CT_NH3 = (PREV_NH3+ICHEMT*P)/(1.0+ICHEMT*L*CF)

! Total cloud NH42SO4  CT_NH42SO4
            P=LAQ(8)
            L=0.0
            CT_NH42SO4 = (PREV_NH42SO4+ICHEMT*P*CF)/(1.0+ICHEMT*L)


          ELSE               ! No cloud water, set rates to zero
            LAQ(1:10)=0.0
            PH=0.11111111           ! Missing data value.

          ENDIF        ! end of if cloud loop

! DO GASEOUS CHEMISTRY (nb: still inside the iteration loop)
! ********************

! OD           Y(JB, 1)
          P = DJ(6)*O3
          L = 0.0 + RC_OD + (RC_OD_H2O * H2O)
          IF(L.GT.1.0)THEN
            OD=P/L
          ELSE
            OD=P
          ENDIF

! OP           Y(JB, 2)
          P = (DJ(1)*NO2) + (DJ(3)*NO3) + (RC_OD*OD) + (DJ(5)*O3)
          L = 0.0 + RC_OP_O2 + (RC_OP_NO*NO) + RC_NO2_OP*NO2
          IF(L.GT.1.0)THEN
            OP = P/L
          ELSE
            OP=P
          ENDIF

! OH          Y(3)
          P = (DJ(7)*H2O2*2.0) + (DJ(8)*HNO3)                       &
              + 2.0*RC_OD_H2O*OD*H2O + DJ(11)*CH3OOH                &
              + RC_OH_CH3OOH*OH*CH3OOH + RC_NO_HO2*NO*HO2           &
              + RC2_NO3_HO2*NO3*HO2 + RC_HO2_O3*HO2*O3              &
              + DJ(11)*ISOPOOH + DJ(11)*MVKOOH                      &
              + RC_MVKOOH_OH*MVKOOH*OH + RC_O3_MVK*O3*MVK*0.36      &
              + RC_ISOPOOH_OH*ISOPOOH*OH + RC2_O3_C3H6*O3*C3H6*0.28 &
              + RC_O3_C5H8*O3*C5H8*0.27 + DJ(17)*HONO
          L = 0.0 + RC_SO2_OH*SO2 + RC_NO2_OH*NO2 + RC_OH_HNO3*HNO3 &
              + RC_OH_H2O2*H2O2 + RC_OH_O3*O3                       &
              + RC_OH_CO*CO + RC_OH_CH4*CH4 + RC_OH_HCHO*HCHO       &
              + RC_OH_H2*H2 + RC2_OH_CH3OOH*CH3OOH                  &
              + RC_OH_CH3OOH*CH3OOH + RC_OH_C2H4*C2H4               &
              + RC_OH_HO2NO2*HO2NO2 + RC_HO2_OH*HO2                 &
              + RC_MGLYOX_OH*MGLYOX + RC_GLYOX_OH*GLYOX             &
              + RC_OH_C5H8*C5H8 + RC_OH_MVK*MVK                     &
              + RC_ISOPOOH_OH*ISOPOOH + RC_MVKOOH_OH*MVKOOH         &
              + RC_OH_OXYL*OXYL + RC_MEMALD_OH*MEMALD               &
              + RC_OH_TOLUEN*TOLUEN + RC_OH_C3H6*C3H6               &
              + RC_OH_PAN*PAN + RC_OH_CH3CHO                        &
              + RC_OH_BD*BD + RC_OH_HONO*HONO                       &
              + RC_OH_NO*NO + RC_OH_H2S*H2S
          IF(L.GT.1.0)THEN
            OH=P/L
          ELSE
            OH=P
          ENDIF

! NO  Y(4) (9nb dj is e-20 if not daytime)

          P=DJ(1)*NO2 +DJ(2)*NO3 + RC2_NO2_NO3*NO2*NO3        &
            + RC_NO2_OP*OP*NO2 + DJ(17)*HONO
          L=RC_NO_O3*O3 + RC_NO_NO3*NO3 + RC_OP_NO*OP         &
            + RC_NO_CH3O2*CH3O2 + RC_CH2O2C_NO*CH2O2C         &
            + RC_NO_HO2*HO2 + RC_RO2IP1_NO*RO2IP1             &
            + RC_RO2IP2_NO*RO2IP2 + RC_NO_TOLP1*TOLP1         &
            + RC_OXYL1_NO*OXYL1 + RC_CH3CHX_NO*CH3CHX         &
            + RC_CH2OO_NO*CH2OO + RC_MEMALD1_NO*MEMALD1       &
            + RC_CH3COO2_NO*CH3COO2 + RC_NO_BDPEROXY*BDPEROXY &
            + RC_OH_NO*OH

          NO=(NO_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

! NO2  Y(5)

          P=DJ(3)*NO3 + DJ(4)*N2O5 + RC_NO_O3*O3*NO             &
            + RC2_NO2_NO3*NO2*NO3 + RC_NO_NO3*NO*NO3            &
            + RC_NO_NO3*NO*NO3 + RC2_NO3_HO2*NO3*HO2            &
            + RC_NO3_NO3*NO3*NO3+ RC_NO3_NO3*NO3*NO3            &
            + RC_N2O5*N2O5 + RC_OP_NO*OP*NO + DJ(8)*HNO3        &
            + RC_NO_CH3O2*NO*CH3O2 + RC_CH2O2C_NO*CH2O2C*NO     &
            + DJ(12)*HO2NO2 + RC_OH_HO2NO2*OH*HO2NO2            &
            + RC_HO2NO2*HO2NO2 + RC_NO_HO2*NO*HO2               &
            + DJ(14)*PAN + RC_OXYL1_NO*NO*OXYL1                 &
            + RC_NO_TOLP1*NO*TOLP1 + RC_RO2IP1_NO*RO2IP1*NO     &
            + RC_RO2IP2_NO*RO2IP2*NO + RC_MEMALD1_NO*MEMALD1*NO &
            + RC_CH2OO_NO*CH2OO*NO + RC_CH3CHX_NO*CH3CHX*NO     &
            + RC_PAN*PAN + RC_CH3COO2_NO*CH3COO2*NO             &
            + RC_NO_BDPEROXY*NO*BDPEROXY                        &
            + RC_OH_HONO*OH*HONO
          L=DJ(1) + (RC_NO2_NO3*NO3) + (RC_NO2_O3*O3)           &
            + (RC_NO2_OH*OH) + RC2_NO2_NO3*NO3                  &
            +  RC_NO2_HO2*HO2 + RC_CH3COO2_NO2*CH3COO2          &
            + RC_CH2OO_NO2*CH2OO + RC_NO2_OP*OP

          NO2=(NO2_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

! NO3 Y(6)

          P=RC_CH2OO_NO2*CH2OO*NO2 + DJ(4)*N2O5       &
            + RC_OH_HNO3*OH*HNO3 + RC_OH_PAN*OH*PAN   &
            + RC_NO2_O3*NO2*O3 + RC_N2O5*N2O5
          L=DJ(2) + DJ(3)                             &
            + RC_NO3_CH3CHO*CH3CHO + RC_NO3_C5H8*C5H8 &
            + RC_NO3_C2H4*C2H4 + RC_NO3_C3H6*C3H6     &
            + RC_NO3_NO3*NO3 + RC_NO3_NO3*NO3         &
            + RC_NO3_HO2*HO2 + RC2_NO3_HO2*HO2        &
            + RC_NO3_HCHO*HCHO + RC_NO_NO3*NO         &
            + RC2_NO2_NO3*NO2 + RC_NO2_NO3*NO2

          NO3=(NO3_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

! N2O5 Y(7)

          P=RC_NO2_NO3*NO2*NO3
          L=RC_N2O5 + RC_NAER2 + DJ(4) + RC_N2O5_2

          N2O5=(N2O5_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

! CO - Y(8)

          P=DJ(9)*HCHO + DJ(10)*HCHO                              &
            + RC_O3_C2H4*O3*C2H4*0.31 + RC_OH_HCHO*OH*HCHO        &
            + RC_NO3_HCHO*NO3*HCHO + DJ(16)*GLYOX + DJ(13)*CH3CHO &
            + DJ(15)*MGLYOX + RC_MGLYOX_OH*MGLYOX*OH              &
            + RC_GLYOX_OH*GLYOX*OH*2.0 + RC_O3_MVK*O3*MVK*0.76    &
            + RC_O3_C3H6*O3*C3H6*0.58 + RC_O3_C5H8*O3*C5H8*0.78   &
            + RC2_O3_C3H6*O3*C3H6*0.4
          L=RC_OH_CO*OH
          CO=(CO_P+ICHEMT*P)/(1.0+ICHEMT*L)

! CH4 - Y(9)

          P=RC2_O3_C3H6*O3*C3H6*0.3
          L=RC_OH_CH4*OH
          CH4=CH4_P/(1.0+ICHEMT*L)

! HCHO - formaldehyde y(10)

          P=DJ(11)*CH3OOH + RC_CH3O2_CH2O2C*CH3O2*CH2O2C*3.0          &
            + RC_O3_C2H4*O3*C2H4 + RC_CH2O2C_NO*CH2O2C*NO*2.0         &
            + RC_CH3O2_CH3O2*CH3O2*CH3O2 + RC_NO_CH3O2*NO*CH3O2       &
            + RC_CH3O2_2*CH3O2*CH3O2*2.0 + RC_SO2_CH3O2*SO2*CH3O2     &
            + RC_OH_CH3OOH*OH*CH3OOH + DJ(11)*ISOPOOH                 &
            + DJ(11)*MVKOOH + DJ(16)*GLYOX + RC_RO2IP2_NO*RO2IP2*NO   &
            + RC_MVKOOH_OH*MVKOOH*OH + RC_RO2IP1_NO*RO2IP1*NO         &
            + RC_RO2IP2_CH3O2*RO2IP2*CH3O2 + RC_ISOPOOH_OH*ISOPOOH*OH &
            + RC_MEMALD1_CH3O2*MEMALD1*CH3O2                          &
            + RC_RO2IP1_CH3O2*RO2IP1*CH3O2 + RC_CH2OO_SO2*CH2OO*SO2   &
            + RC_CH2OO_NO2*CH2OO*NO2 + RC_CH2OO_H2O*CH2OO             &
            + RC_CH3O2_CH3CHX*CH3O2*CH3CHX*2.0                        &
            + RC_CH2OO_NO*CH2OO*NO + RC2_O3_C3H6*O3*C3H6              &
            + RC_CH3CHX_NO*CH3CHX*NO + RC_OH_PAN*OH*PAN               &
            + RCA_CH3O2_CH3COO2*CH3O2*CH3COO2
          L=DJ(10) + RC_OH_HCHO*OH + RC_NO3_HCHO*NO3 + DJ(9)
          HCHO=(HCHO_P+ICHEMT*P)/(1.0+ICHEMT*L)

! HONO - nitrous acid

          P=RC_OH_NO*OH*NO + RC_NO2*NO2
          L=RC_OH_HONO*OH + DJ(17)

          HONO=(HONO_P+ICHEMT*P)/(1.0+ICHEMT*L)

! O3           Y(11)

          P = RC_OP_O2*OP
          L = DJ(6) + DJ(5) + (RC_NO_O3*NO) + (RC_NO2_O3*NO2)    &
             + RC_O3_C2H4*C2H4 + RC_HO2_O3*HO2 + RC_OH_O3*OH     &
             + RC_O3_MVK*MVK + RC_O3_C3H6*C3H6 + RC_O3_C5H8*C5H8 &
             + RC2_O3_C3H6*C3H6
          O3 = (O3_P+ICHEMT*P)/(1.0+ICHEMT*L)

! H2   Y(12)

          P=DJ(10)*HCHO + RC_O3_C2H4*O3*C2H4*0.13 + RC_O3_C3H6*O3*C3H6
          L=RC_OH_H2*OH
          H2=(H2_P+ICHEMT*P)/(1.0+ICHEMT*L)

! HNO3 -nitric acid - Y(13)

          P=RC_NO2_OH*NO2*OH + RC_NO3_HO2*NO3*HO2             &
            + RC_NO3_HCHO*NO3*HCHO + RC_NO3_CH3CHO*NO3*CH3CHO &
            + RC_N2O5_2*N2O5 + RC_N2O5_2*N2O5

          L=RC_NAER1 + RC_OH_HNO3*OH + DJ(8)
          HNO3=(HNO3_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

! H2O2 - hydrogen peroxide Y(14)

          P=RC_HO2_HO2*HO2*HO2
          L=LAQ(4) + DJ(7) + RC_OH_H2O2*OH
          H2O2=(H2O2_P+(ICHEMT*P))/(1.0 + ICHEMT*L)

! CH3O2 - methylperoxy Y(15)

          P=RC2_OH_CH3OOH*OH*CH3OOH + RC_OH_CH4*OH*CH4               &
            + DJ(13)*CH3CHO + RC_CH3COO2_CH3COO2*CH3COO2*CH3COO2*2.0 &
            + RC2_O3_C3H6*O3*C3H6 + RC_CH3COO2_NO*CH3COO2*NO         &
            + RCA_CH3O2_CH3COO2*CH3O2*CH3COO2
          L=RC_CH3O2_CH3O2*CH3O2 + RC_CH3O2_CH3O2*CH3O2              &
            + RC_CH3O2_HO2*HO2 + RC_SO2_CH3O2*SO2                    &
            + RC_NO_CH3O2*NO + RC_CH3O2_2*CH3O2                      &
            + RC_CH3O2_2*CH3O2 + RC_CH3O2_CH2O2C*CH2O2C              &
            + RC_MEMALD1_CH3O2*MEMALD1 + RC_RO2IP1_CH3O2*RO2IP1      &
            + RC_RO2IP2_CH3O2*RO2IP2 + RC_CH3O2_CH3CHX*CH3CHX
          CH3O2=(CH3O2_P+ICHEMT*P)/(1.0 + (ICHEMT*L))

! HO2 hydroperoxy Y(16)

          P=DJ(11)*CH3OOH + DJ(12)*HO2NO2                                 &
            + DJ(9)*HCHO + DJ(9)*HCHO                                     &
            + RC_NO3_C2H4*NO3*C2H4 + RC_CH2O2C_NO*CH2O2C*NO               &
            + RC_O3_C2H4*O3*C2H4*0.2 + RC_NO3_HCHO*NO3*HCHO               &
            + RC_OH_CO*OH*CO + RC_OH_HCHO*OH*HCHO                         &
            + RC_NO_CH3O2*NO*CH3O2 + RC_CH3O2_2*CH3O2*CH3O2*2.0           &
            + RC_SO2_OH*OH*SO2 + RC_SO2_CH3O2*SO2*CH3O2                   &
            + RC_OH_H2O2*OH*H2O2 + RC_OH_H2*OH*H2                         &
            + RC_OH_O3*OH*O3 + RC_HO2NO2*HO2NO2                           &
            + DJ(11)*ISOPOOH + DJ(11)*MVKOOH + DJ(15)*MGLYOX              &
            + DJ(13)*CH3CHO + RC_RO2IP2_NO*RO2IP2*NO                      &
            + RC_RO2IP1_NO*RO2IP1*NO + RC_RO2IP2_CH3O2*RO2IP2*CH3O2*2.0   &
            + RC_GLYOX_OH*GLYOX*OH + RC_MEMALD1_CH3O2*MEMALD1*CH3O2*2.0   &
            + RC_RO2IP1_CH3O2*RO2IP1*CH3O2*2.0                            &
            + RC_NO_TOLP1*NO*TOLP1 + RC_OXYL1_NO*OXYL1*NO                 &
            + RC_NO3_C5H8*NO3*C5H8 + RC_NO3_C3H6*NO3*C3H6                 &
            + RC_O3_C5H8*O3*C5H8*0.27 + RC_MEMALD1_NO*MEMALD1*NO          &
            + RC_O3_MVK*O3*MVK*0.36 + RC_CH3CHX_NO*CH3CHX*NO              &
            + RC_CH3O2_CH3CHX*CH3O2*CH3CHX*2.0                            &
            + RC2_O3_C3H6*O3*C3H6*0.3 + RC_O3_C3H6*O3*C3H6*0.18           &
            + RCA_CH3O2_CH3COO2*CH3O2*CH3COO2                             &
            + RC_NO_BDPEROXY*NO*BDPEROXY + RC_OH_H2S*H2S*OH

          L=RC_CH3O2_HO2*CH3O2 + RC_NO3_HO2*NO3                           &
            + RC2_NO3_HO2*NO3 + RC_HO2_HO2*HO2 + RC_HO2_HO2*HO2           &
            + RC_HO2_O3*O3 + RC_NO_HO2*NO + RC_NO2_HO2*NO2                &
            + RC_HO2_OH*OH + RC_RO2IP2_HO2*RO2IP2                         &
            + RC_RO2IP1_HO2*RO2IP1

! Changed to calculate amount using value from previous timestep
!          IF(L.GT.1.0)THEN
!            HO2=P/L
!          ELSE
!            HO2=P
!          ENDIF

           HO2=(HO2_P+ICHEMT*P)/(1.0 + (ICHEMT*L))

! CH3CHO Y(19)

          P=RC_CH3O2_CH3CHX*CH3O2*CH3CHX+RC_O3_C3H6*O3*C3H6+RC_CH3CHX_NO*CH3CHX*NO
          L=RC_OH_CH3CHO*OH + RC_NO3_CH3CHO*NO3 + DJ(13)
          CH3CHO=(CH3CHO_P+ICHEMT*P)/(1.0+(ICHEMT*L))

! CH3COO2 Y(20)

          P=DJ(14)*PAN + DJ(15)*MGLYOX                              &
            + RC_NO3_CH3CHO*NO3*CH3CHO + RC_MGLYOX_OH*MGLYOX*OH     &
            + RC_OH_CH3CHO*OH*CH3CHO + RC_PAN*PAN
          L=RC_CH3COO2_CH3COO2*CH3COO2 + RC_CH3COO2_CH3COO2*CH3COO2 &
            + RC_CH3O2_CH3COO2*CH3O2 + RC_CH3COO2_NO2*NO2           &
            + RC_CH3COO2_NO*NO + RCA_CH3O2_CH3COO2*CH3O2
          CH3COO2=(CH3COO2_P+ICHEMT*P)/(1.0+ICHEMT*L)

! PAN - CH3COO2NO2 Y(21)

          P=RC_CH3COO2_NO2*CH3COO2*NO2
          L=RC_PAN + RC_OH_PAN*OH + DJ(14)
          PAN=(PAN_P+ICHEMT*P)/(1.0+ICHEMT*L)

! CH3OOH - methyl hydroperoxide Y(22)

          P=RC_CH3O2_HO2*CH3O2*HO2
          L=RC2_OH_CH3OOH*OH + RC_OH_CH3OOH*OH + DJ(11)
          CH3OOH=(CH3OOH_P+ICHEMT*P)/(1.0 + (ICHEMT*L))

! SO2 - sulphur dioxide Y(26)

          P=RC_OH_H2S*H2S*OH
          L=RC_SO2_OH*OH + RC_SO2_CH3O2*CH3O2 + RC_CH2OO_SO2*CH2OO
          SO2=SO2_P/(1.0+ICHEMT*L)

! C2H4 - ethylene Y(27)

          P=0.0
          L=RC_OH_C2H4*OH + RC_O3_C2H4*O3 + RC_NO3_C2H4*NO3
          C2H4=C2H4_P/(1.0+ICHEMT*L)

! C3H6 - propylene Y(28)

          P=0.0
          L=RC_NO3_C3H6*NO3 + RC2_O3_C3H6*O3 + RC_O3_C3H6*O3 + RC_OH_C3H6*OH
          C3H6=C3H6/(1.0+ICHEMT*L)

! OXYL - o-xylene Y(28) (C6H4)(CH3)2

          P=0.0
          L=RC_OH_OXYL*OH
          OXYL=OXYL_P/(1.0+ICHEMT*L)

! SO4 -  sulphate aerosol Y(30)

          P=RC_SO2_OH*OH*SO2 + LAQ(2) + RC_SO2_CH3O2*SO2*CH3O2 + RC_CH2OO_SO2*CH2OO*SO2
          L=LAQ(9)
          SO4=(SO4_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! CH2OO - Y(39)

          P=RC_O3_MVK*O3*MVK*0.24 + RC_O3_C2H4*O3*C2H4*0.47 + RC_O3_C5H8*O3*C5H8*0.22
          L=RC_CH2OO_NO*NO + RC_CH2OO_NO2*NO2                              &
                  + RC_CH2OO_SO2*SO2 + RC_CH2OO_H2O*H2O + RC_CH2OO_SO2*SO2
          CH2OO=(CH2OO_P+(ICHEMT*P))/(1.0+ICHEMT*L)

! CH2O2C - ethylene peroxy Y(41)

          P=RC_OH_C2H4*OH*C2H4
          L=RC_CH2O2C_NO*NO + RC_CH3O2_CH2O2C*CH3O2
          CH2O2C=(CH2O2C_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! MGLYOX - Y(42)-CH3COCHO

          P=RC_RO2IP2_NO*RO2IP2*NO + RC_MVKOOH_OH*MVKOOH*OH           &
            + DJ(11)*MVKOOH + RC_MEMALD1_CH3O2*MEMALD1*CH3O2          &
            + RC_RO2IP2_CH3O2*RO2IP2*CH3O2 + RC_MEMALD1_NO*MEMALD1*NO &
            + RC_O3_MVK*O3*MVK + RC_OXYL1_NO*NO*MEMALD1
          L=RC_MGLYOX_OH*OH + DJ(15)
          MGLYOX=(MGLYOX_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! CH3CHX - Y(43)

          P=RC_OH_C3H6*OH*C3H6
          L=RC_CH3CHX_NO*NO + RC_CH3O2_CH3CHX*CH3O2
          CH3CHX=(CH3CHX_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! GLYOX - Y(44)

          P=RC_MEMALD1_CH3O2*MEMALD1*CH3O2                    &
            + RC_NO_TOLP1*NO*TOLP1 + RC_MEMALD1_NO*MEMALD1*NO &
            + RC_NO_BDPEROXY*NO*BDPEROXY
          L=RC_GLYOX_OH*OH + DJ(16)
          GLYOX=(GLYOX_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! OXYL1 - Y(45)

          P=RC_OH_OXYL*OH*OXYL
          L=RC_OXYL1_NO*NO
          OXYL1=(OXYL1_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! MEMALD - Y(46)

          P=RC_OXYL1_NO*NO*MEMALD1 + RC_NO_TOLP1*NO*TOLP1
          L=RC_MEMALD_OH*OH
          MEMALD=(MEMALD_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! MEMALD1 - Y(47)

          P=RC_MEMALD_OH*OH*MEMALD
          L=RC_MEMALD1_NO*NO + RC_MEMALD1_CH3O2
          MEMALD1=(MEMALD1_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! ISOPRENE C5H8 - Y(48)

          P=0.0
          L=RC_O3_C5H8*O3 + RC_NO3_C5H8*NO3 + RC_OH_C5H8*OH
          C5H8=(C5H8_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! RO2IP1 - Y(49)

          P=RC_OH_C5H8*OH*C5H8
          L=RC_RO2IP1_CH3O2*CH3O2 + RC_RO2IP1_HO2*HO2 + RC_RO2IP1_NO*NO
          RO2IP1=(RO2IP1_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! MVK - Y(50)

          P=RC_RO2IP1_NO*RO2IP1*NO + DJ(11)*ISOPOOH                   &
            + RC_RO2IP1_CH3O2*RO2IP1*CH3O2 + RC_ISOPOOH_OH*ISOPOOH*OH &
            + RC_O3_C5H8*O3*C5H8
          L=RC_O3_MVK*O3 + RC_OH_MVK*OH
          MVK=(MVK_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! RO2IP2 - Y(51)

          P=RC_OH_MVK*OH*MVK
          L=RC_RO2IP2_CH3O2*CH3O2 + RC_RO2IP2_HO2*HO2 + RC_RO2IP2_NO*NO
          RO2IP2=(RO2IP2_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! ISOPOOH - Y(52)

          P=RC_RO2IP1_HO2*RO2IP1*HO2
          L=RC_ISOPOOH_OH*OH + DJ(11)
          ISOPOOH=(ISOPOOH_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! MVKOOH - Y(53)

          P=RC_RO2IP2_HO2*RO2IP2*HO2
          L=RC_MVKOOH_OH*OH + DJ(11)
          MVKOOH=(MVKOOH_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! TOLUEN - Y(54)

          P=0.0
          L=RC_OH_TOLUEN*OH + RC2_OH_TOLUEN*OH + RC_OH_AROPREC*OH
          TOLUEN=(TOLUEN_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! AROPREC - Y(65)
!          BYPROD=0.35
!          BYPROD=1.84
          VNO=NO/2.55E+10
          REACFRAC=1.84*(0.165+(0.886*VNO)-(0.289*VNO*VNO) &
             +(0.023*VNO*VNO*VNO))
          BYPROD=MIN(1.84_P64,REACFRAC)         
          P=BYPROD*RC2_OH_TOLUEN*OH*TOLUEN + KOUT*NAMEARO
          L=RC_OH_AROPREC*OH + KIN
          AROPREC=(AROPREC_P + (ICHEMT*P))/(1.0+ICHEMT*L)

! NAMEARO - Y(66) anthropogenic SOA
          P=KIN*AROPREC
          L=KOUT
          NAMEARO=(NAMEARO_P + (ICHEMT*P))/(1.0+ICHEMT*L)


! NAER - NITRATE AEROSOL Y(58)
! lost only by wet and dry dep
! N2O5 is in p twice as 1 n2o5 gives 2 no3-'s

          P=(RC_NAER2*N2O5) + (RC_NAER2*N2O5) + (RC_NAER1*HNO3)
          L=0.0
          NAER=(NAER_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

! TOLP1 - Y(59)

          P=RC_OH_TOLUEN*OH*TOLUEN
          L=RC_NO_TOLP1*NO
          TOLP1=(TOLP1_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

! HO2NO2 - peroxynitric acid Y(65)

          P=RC_NO2_HO2*NO2*HO2
          L=RC_HO2NO2 + RC_OH_HO2NO2*OH + DJ(12)
          HO2NO2=(HO2NO2_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

! BD - 1,3 butadiene (not in STOCHEM)

          P=0.0
          L=RC_OH_BD*OH
          BD=(BD_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

! BDPEROXY - 1,3-butadiene-peroxy (not in STOCHEM)

          P=RC_OH_BD*OH*BD
          L=RC_NO_BDPEROXY*NO
          BDPEROXY=(BDPEROXY_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

! H2S - Hydrogen Sulfide (not in STOCHEM, added by Alex Archibald)

          P=0.0
          L=RC_OH_H2S*OH
          H2S=(H2S_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

! APINENE C10H16BIO

          P=0.0
          L=RC_O3_C10H16BIO*O3 + RC_NO3_C10H16BIO*NO3 &
           +RC_OH_C10H16BIO*OH
          C10H16BIO=(C10H16BIO_P + (ICHEMT*P))/(1.0+ICHEMT*L)


!  SOA PRECURSOR

          P=RC_O3_C10H16BIO*O3*C10H16BIO &
           +RC_NO3_C10H16BIO*NO3*C10H16BIO    &
           +RC_OH_C10H16BIO*OH*C10H16BIO
          L= RC_OH_SOAP*OH &
           +RC_SCAV
          SOAP=(SOAP_P+ICHEMT*P)/(1.0 + (ICHEMT*L))

! SOA

          P=RC_SCAV*SOAP
          L=0.0
          SOA=(SOA_P + (ICHEMT*P))/(1.0 + (ICHEMT*L))

        ENDDO      ! end of iteration loop

      ENDDO    ! END OF CHEM TIMESTEP LOOP

! now have new values of ct's for aq phase and conc's for gas
! phase - need to mix these in for SO2
! NH42SO4 doesn't do anything in the gas phase
! SO4, H2O2  are mixed into the gas phase in the iteration loop
! so don't have to do anything to these (not sure why so2
! isn't too - but thats just how it is in STOCHEM

! convert ct's (moles/litre) back to molecules/cm3 for output:
! allow the gas phase to occur everywhere regardless of
! the cloud fraction (changed 22/11/00) like STOCHEM

      IF(CF.GT.0.0.AND.CLDWAT.GT.0.0)THEN
        IF(SO2_P.GT.0.0)THEN
          SO2 = SO2 - ((INIT_SO2-CT_SO2)*(AVOGAD/1.0E3)*CF)
        ELSE
          SO2=0.0
        ENDIF
        IF(SO2.LT.0.0)THEN
          SO2=0.0
!         WRITE(86,*),'Cloud mixing, SO2 -ve; now reset',D,E
        ENDIF
        IF(NH3_P.GT.0.0)THEN
          NH3 = NH3 - ((INIT_NH3-CT_NH3)*(AVOGAD/1.0E3)*CF)
        ELSE
          NH3=0.0
        ENDIF
        IF(NH3.LT.0.0)THEN
          NH3=0.0
!         WRITE(86,*),'Cloud mixing, NH3 -ve; now reset',D,E
        ENDIF
        NH42SO4 = CT_NH42SO4 * AVOGAD/1.0E3
      ENDIF

! now calc the eqm between nh3 hno3 and nh4no3:
! new nitrate scheme is from Ackermann et al (1995)
! NAME currently has molecules/cm3
! We are doing the conversion by modifying the rate and not the concentration units

         TNH3=(NH3+NH4NO3)      ! Total available NH3
         THNO3=(HNO3+NH4NO3)    ! Total available HNO3


         IF(TNH3*THNO3 > Kp_RT2) THEN

           CONST=((TNH3+THNO3)**2)-4.0*((TNH3*THNO3)-KP_RT2)

           IF ( CONST > (((TNH3+THNO3)*(TNH3+THNO3))*10E-8)) THEN

! If b squared is significant then we can solve quadratic in traditional manner


             IF ( CONST <0) THEN
               print*,'complex roots!'
             ENDIF

             NH4NO3_EQM=0.5*((TNH3+THNO3)-SQRT(CONST))

           ELSE

! Otherwise neglect the b2-4ac and just solve the rest of the equation
!x=-b/2a

             NH4NO3_EQM=(TNH3+THNO3)/2.0
           ENDIF


           NH3_EQM=TNH3-NH4NO3_EQM
           HNO3_EQM=THNO3-NH4NO3_EQM

         ELSE

           NH4NO3_EQM=0.0
           NH3_EQM=TNH3
           HNO3_EQM=THNO3

         ENDIF

! Subtract the regional background component for CO and HCHO

      IF(CO.GT.2.9529E12)THEN
        CO=CO-2.9529E12
      ELSE
        CO=0.0
      ENDIF
      IF(HCHO.GT.3.6595E10)THEN
        HCHO=HCHO-3.6594E10
      ELSE
        HCHO=0.0
      ENDIF

! Chem species 
      CONC(1)  = SO2
      CONC(2)  = NH3_EQM
      CONC(3)  = NO
      CONC(4)  = NH42SO4
      CONC(5)  = SO4
      CONC(6)  = NO2
      CONC(7)  = NO3
      CONC(8)  = N2O5
      CONC(9)  = HNO3_EQM
      CONC(10) = NAER
      CONC(11) = NH4NO3_EQM
      CONC(12) = O3
      CONC(13) = CO
      CONC(14) = HCHO
      CONC(15) = C2H4
      CONC(16) = C3H6
      CONC(17) = 0.0
      CONC(18) = OXYL
      CONC(19) = TOLUEN
      CONC(20) = BD
      CONC(21) = CH3CHO
      CONC(22) = PAN
      CONC(23) = HONO
      CONC(25) = H2S
      CONC(26) = C5H8
      CONC(27) = C10H16BIO
      CONC(28) = SOA
      CONC(29) = SOAP
      CONC(30) = AROPREC
      CONC(31) = NAMEARO
      CONC(33) = H2O2
      Conc(34) = OH
      Conc(35) = HO2
      Conc(36) = CH3OOH
      Conc(37) = MVK
      Conc(38) = ISOPOOH
      Conc(39) = MVKOOH
      Conc(40) = MGLYOX
      Conc(41) = GLYOX
      Conc(42) = MEMALD

!nb the isoprene is now a combination of anthropogenic contribution
!and a biogenic contribution & hence is only carried on one species now - it
!doesn't matter which.

End Subroutine ChemistryScheme

!-------------------------------------------------------------------------------------------------------------

Subroutine EQMCON(T, KE_SO2, KE_HSO3, KE_NH3, KE_CO2, KE_H2O, KE_HNO3, &
                  KH_HNO3, KH_O3, KH_H2O2, KH_SO2, KH_NH3, KH_CO2)
! Calculates Henry's law and equilibrium constants.

! This is a Fortran 90 version of the subroutine EQMCON in NAME V8.08.

  Implicit None
  ! Argument list:
  REAL(8),    Intent(In)  :: T            ! Absolute temperature.
  REAL(8),    Intent(Out) :: KE_HNO3      !} Equilibrium coefficients.
  REAL(8),    Intent(Out) :: KE_SO2       !}
  REAL(8),    Intent(Out) :: KE_HSO3      !}
  REAL(8),    Intent(Out) :: KE_NH3       !}
  REAL(8),    Intent(Out) :: KE_CO2       !}
  REAL(8),    Intent(Out) :: KE_H2O       !}
  REAL(8),    Intent(Out) :: KH_O3        !] Henry's law coefficients.
  REAL(8),    Intent(Out) :: KH_HNO3      !]
  REAL(8),    Intent(Out) :: KH_H2O2      !]
  REAL(8),    Intent(Out) :: KH_SO2       !]
  REAL(8),    Intent(Out) :: KH_NH3       !]
  REAL(8),    Intent(Out) :: KH_CO2       !]
  ! Locals:


! local scalar

      REAL(8) ::  REUS1             !conversion factor


! Calculate equilibrium coefficients:

      REUS1= ((1.0/T)-(1.0/298.0))

! HNO3 = NO3(-) + H(+)
      KE_HNO3 = 1.8E-5*DEXP(-450.0*REUS1)

! SO2(aq) = H(+) + HSO3(-)
      KE_SO2 = 1.7E-02*DEXP(2090.0*REUS1)

! HSO3(-) = H(+) + SO3(2-)
      KE_HSO3 = 6.0E-08*DEXP(1120.0*REUS1)

! NH3 + H2O = NH4(+) + OH(-)
      KE_NH3 = 1.8E-5*DEXP(-450.0*REUS1)

! CO2 = H(+) + HCO3(-)
      KE_CO2 = 4.3E-7*DEXP(-913.0*REUS1)

! H2O = H(+) + OH(-)
      KE_H2O = (1000.0/18.0)*1.8E-16*DEXP(-6716.0*REUS1)

! Henrys law constants (mol/(l.atm)).
! Aqueous phase equilibria.

! O3
      KH_O3=1.1E-2*DEXP(2300.0*REUS1)
! HNO3
      KH_HNO3=3.3E+6*DEXP(8700.0*REUS1)
! H2O2
      KH_H2O2=7.36E+4*DEXP(6621.0*REUS1)
! SO2
      KH_SO2=1.23E+0*DEXP(3120.0*REUS1)
! NH3
      KH_NH3=7.5E+1*DEXP(3400.0*REUS1)
! CO2
      KH_CO2=3.4E-2*DEXP(2420.0*REUS1)

End Subroutine EQMCON

!-------------------------------------------------------------------------------------------------------------

Subroutine CHEMCO(T, AIR_DEN, H2O, O2, RC_HO2_HO2, RC_SO2_OH,          &
        RC_HSO3_H2O2, RC_HSO3_O3, RC_SO3_O3, RC_NO_O3, RC_NO2_O3,      &
        RC_NO2_OH, RC_OH_HNO3, RC_NO3_NO3, RC_NO3_HO2,                 &
        RC_NO_NO3, RC2_NO2_NO3, RC_N2O5, RC_OP_O2, RC_OP_NO, RC_OD,    &
!        EQM_NH3HNO3,                                                   &
        RC_O3_C2H4, RC_OH_HCHO, RC_NO3_HCHO, RC_OH_CO,                 &
        RC_OH_CH4, RC_CH3O2_CH2O2C, RC_CH2O2C_NO, RC_CH3O2_HO2,        &
        RC_CH3O2_CH3O2, RC_NO_CH3O2, RC_CH3O2_2, RC_SO2_CH3O2,         &
        RC_OH_CH3OOH, RC2_OH_CH3OOH, RC_OH_C2H4, RC_NO3_C2H4,          &
        RC_NO2_HO2, RC_HO2NO2, RC_OH_HO2NO2, RC_OH_H2O2, RC_OH_H2,     &
        RC_HO2_O3, RC_NO_HO2, RC_HO2_OH, RC_NO2_NO3, RC2_NO3_HO2,      &
        RC_OD_H2O, RC_OH_O3, RC_CH3O2_CH3CHX, RC_O3_C3H6,              &
        RC_CH3CHX_NO, RC_OH_CH3CHO, RC_NO3_CH3CHO, RC_MGLYOX_OH,       &
        RC_PAN, RC_CH3COO2_CH3COO2, RC_CH3O2_CH3COO2, RC_CH3COO2_NO2,  &
        RC_CH3COO2_NO, RCA_CH3O2_CH3COO2, RC_OH_PAN, RC_NO3_C3H6,      &
        RC2_O3_C3H6, RC_OH_C3H6, RC_OH_OXYL, RC_O3_MVK,                &
        RC_O3_C5H8, RC_CH2OO_NO, RC_CH2OO_NO2, RC_CH2OO_SO2,           &
        RC_CH2OO_H2O, RC_RO2IP2_NO, RC_MVKOOH_OH, RC_MEMALD1_CH3O2,    &
        RC_RO2IP2_CH3O2, RC_MEMALD1_NO, RC_OXYL1_NO,                   &
        RC_NO_TOLP1, RC_GLYOX_OH, RC_MEMALD_OH, RC_NO3_C5H8,           &
        RC_OH_C5H8, RC_RO2IP1_CH3O2, RC_RO2IP1_HO2, RC_RO2IP1_NO,      &
        RC_ISOPOOH_OH, RC_OH_MVK, RC_RO2IP2_HO2, RC_OH_TOLUEN,         &
        RC_NO2_OP, RC_OH_BD, RC_NO_BDPEROXY, RC_OH_HONO, RC_OH_NO,     &
        RC_NO2, ZI, RH, Kp_RT2, RC_NAER1, RC_NAER2, RC_OH_H2S,         &
        BLFLAG, RC_N2O5_2,RC_O3_C10H16BIO,RC_NO3_C10H16BIO,            &
        RC_OH_C10H16BIO,RC_OH_SOAP, RC_SCAV, RC2_OH_TOLUEN,            &
        RC_OH_AROPREC, KIN, KOUT)

! Calculates all the rate coefficients for gas and aqueous phase reactions.

! This is a Fortran 90 version of the subroutine CHEMCO in NAME V8.08.

!====================================================================
! LIST OF CHANGES TO CHEMCO 
!(1) Included an Intent attribute for each argument.
!(2) 30/1/07 New nitrate calculation ref: Ackermann et al 1995
!(3) 30/1/07 Reactions 5,12,22,45,131,245,247 re-instated
!(4) 6/11/08 H2S reaction added by Alex Archibald, ref: IUPAC subcommittee 2004
!====================================================================

  Implicit None
  ! Argument list:
  REAL(8),    Intent(In)  :: T                 !IN: absolute temp
  REAL(8),    Intent(In)  :: AIR_DEN           !IN: air density
  REAL(8),    Intent(In)  :: H2O               !IN: water concentration
  REAL(8),    Intent(In)  :: O2                !IN: oxygen density
  REAL(8),    Intent(In)  :: ZI                !IN: boundary layer depth in m
  REAL(8),    Intent(In)  :: RH                !IN: RELATIVE HUMIDITY
  REAL(8),    Intent(Out) :: RC_HO2_HO2        !OUT: ho2+ho2 rate
  REAL(8),    Intent(Out) :: RC_SO2_OH          !OUT: SO2+ OH rate
  REAL(8),    Intent(Out) :: RC_HSO3_H2O2       !OUT: HSO3+H2O2 rate
  REAL(8),    Intent(Out) :: RC_HSO3_O3         !OUT: HSO3+O3 rate
  REAL(8),    Intent(Out) :: RC_SO3_O3          !OUT: SO3+O3 rate
  REAL(8),    Intent(Out) :: RC_NO_O3           !OUT: NO+O3 rate
  REAL(8),    Intent(Out) :: RC_NO2_O3          !OUT: NO2+O3 rate
  REAL(8),    Intent(Out) :: RC_NO2_NO3         !OUT: NO2+NO3 rate
  REAL(8),    Intent(Out) :: RC_NO2_OH          !OUT: NO2+OH rate
  REAL(8),    Intent(Out) :: RC_NAER1           !OUT: NAER1 rate
  REAL(8),    Intent(Out) :: RC_NAER2           !OUT: NAER2 rate
  REAL(8),    Intent(Out) :: RC_OH_HNO3         !OUT: OH+HNO3 rate
  REAL(8),    Intent(Out) :: RC_NO3_NO3         !OUT: NO3+NO3
  REAL(8),    Intent(Out) :: RC_NO3_HO2         !OUT: NO3+HO2
  REAL(8),    Intent(Out) :: RC2_NO3_HO2        !OUT: NO3+HO2
  REAL(8),    Intent(Out) :: RC_NO_NO3          !OUT: NO + NO3
  REAL(8),    Intent(Out) :: RC2_NO2_NO3        !OUT: NO2+NO3
  REAL(8),    Intent(Out) :: RC_N2O5           !OUT: N2O5
  REAL(8),    Intent(Out) :: RC_N2O5_2         !OUT: N2O5=2HNO3
  REAL(8),    Intent(Out) :: RC_OP_O2           !OUT: OP+O2 RATE
  REAL(8),    Intent(Out) :: RC_OP_NO           !OUT: OP+NO RATE
  REAL(8),    Intent(Out) :: RC_OD             !OUT: OD RATE
  REAL(8),    Intent(Out) :: RC_OD_H2O          !OUT: OD + H2O
  REAL(8),    Intent(Out) :: RC_O3_C2H4        !OUT: O3+C2H4
  REAL(8),    Intent(Out) :: RC_OH_HCHO        !OUT: OH+HCHO
  REAL(8),    Intent(Out) :: RC_NO3_HCHO       !OUT: NO3+HCHO
  REAL(8),    Intent(Out) :: RC_OH_CO          !OUT: OH+CO
  REAL(8),    Intent(Out) :: RC_OH_CH4         !OUT: OH+CH4
  REAL(8),    Intent(Out) :: RC_CH3O2_CH2O2C
  REAL(8),    Intent(Out) :: RC_CH2O2C_NO
  REAL(8),    Intent(Out) :: RC_CH3O2_CH3O2
  REAL(8),    Intent(Out) :: RC_NO_CH3O2
  REAL(8),    Intent(Out) :: RC_CH3O2_2
  REAL(8),    Intent(Out) :: RC_SO2_CH3O2
  REAL(8),    Intent(Out) :: RC_OH_CH3OOH
  REAL(8),    Intent(Out) :: RC2_OH_CH3OOH
  REAL(8),    Intent(Out) :: RC_CH3O2_HO2
  REAL(8),    Intent(Out) :: RC_OH_C2H4
  REAL(8),    Intent(Out) :: RC_NO3_C2H4
  REAL(8),    Intent(Out) :: RC_NO2_HO2
  REAL(8),    Intent(Out) :: RC_HO2NO2
  REAL(8),    Intent(Out) :: RC_OH_HO2NO2
  REAL(8),    Intent(Out) :: RC_OH_H2O2
  REAL(8),    Intent(Out) :: RC_OH_H2
  REAL(8),    Intent(Out) :: RC_OH_O3
  REAL(8),    Intent(Out) :: RC_HO2_O3
  REAL(8),    Intent(Out) :: RC_NO_HO2
  REAL(8),    Intent(Out) :: RC_HO2_OH
  REAL(8),    Intent(Out) :: RC_CH3O2_CH3CHX
  REAL(8),    Intent(Out) :: RC_O3_C3H6
  REAL(8),    Intent(Out) :: RC_CH3CHX_NO
  REAL(8),    Intent(Out) :: RC_OH_CH3CHO
  REAL(8),    Intent(Out) :: RC_NO3_CH3CHO
  REAL(8),    Intent(Out) :: RC_MGLYOX_OH
  REAL(8),    Intent(Out) :: RC_PAN
  REAL(8),    Intent(Out) :: RC_CH3COO2_CH3COO2
  REAL(8),    Intent(Out) :: RC_CH3O2_CH3COO2
  REAL(8),    Intent(Out) :: RC_CH3COO2_NO2
  REAL(8),    Intent(Out) :: RC_CH3COO2_NO
  REAL(8),    Intent(Out) :: RCA_CH3O2_CH3COO2
  REAL(8),    Intent(Out) :: RC_OH_PAN
  REAL(8),    Intent(Out) :: RC_NO3_C3H6
  REAL(8),    Intent(Out) :: RC2_O3_C3H6
  REAL(8),    Intent(Out) :: RC_OH_C3H6
  REAL(8),    Intent(Out) :: RC_OH_OXYL
  REAL(8),    Intent(Out) :: RC_O3_MVK
  REAL(8),    Intent(Out) :: RC_O3_C5H8
  REAL(8),    Intent(Out) :: RC_CH2OO_NO
  REAL(8),    Intent(Out) :: RC_CH2OO_NO2
  REAL(8),    Intent(Out) :: RC_CH2OO_SO2
  REAL(8),    Intent(Out) :: RC_CH2OO_H2O
  REAL(8),    Intent(Out) :: RC_RO2IP2_NO
  REAL(8),    Intent(Out) :: RC_MVKOOH_OH
  REAL(8),    Intent(Out) :: RC_MEMALD1_CH3O2
  REAL(8),    Intent(Out) :: RC_RO2IP2_CH3O2
  REAL(8),    Intent(Out) :: RC_MEMALD1_NO
  REAL(8),    Intent(Out) :: RC_NO_TOLP1
  REAL(8),    Intent(Out) :: RC_GLYOX_OH
  REAL(8),    Intent(Out) :: RC_MEMALD_OH
  REAL(8),    Intent(Out) :: RC_NO3_C5H8
  REAL(8),    Intent(Out) :: RC_OH_C5H8
  REAL(8),    Intent(Out) :: RC_RO2IP1_CH3O2
  REAL(8),    Intent(Out) :: RC_RO2IP1_HO2
  REAL(8),    Intent(Out) :: RC_RO2IP1_NO
  REAL(8),    Intent(Out) :: RC_OXYL1_NO
  REAL(8),    Intent(Out) :: RC_ISOPOOH_OH
  REAL(8),    Intent(Out) :: RC_OH_MVK
  REAL(8),    Intent(Out) :: RC_RO2IP2_HO2
  REAL(8),    Intent(Out) :: RC_OH_TOLUEN
  REAL(8),    Intent(Out) :: RC_NO2_OP
  REAL(8),    Intent(Out) :: RC_OH_BD
  REAL(8),    Intent(Out) :: RC_NO_BDPEROXY
  REAL(8),    Intent(Out) :: RC_OH_HONO
  REAL(8),    Intent(Out) :: RC_OH_NO
  REAL(8),    Intent(Out) :: RC_NO2
  REAL(8),    Intent(Out) :: RC_OH_H2S      !New reaction for H2S
!  REAL(8),    Intent(Out) :: EQM_NH3HNO3   !OUT: NH3 + HNO3
  REAL(8),    Intent(Out) :: Kp_RT2
  REAL(8),    Intent(Out) :: RC_O3_C10H16BIO
  REAL(8),    Intent(Out) :: RC_NO3_C10H16BIO
  REAL(8),    Intent(Out) :: RC_OH_C10H16BIO
  REAL(8),    Intent(Out) :: RC_OH_SOAP
  REAL(8),    Intent(Out) :: RC_SCAV
  REAL(8),    Intent(Out) :: RC2_OH_TOLUEN
  REAL(8),    Intent(Out) :: RC_OH_AROPREC
  REAL(8),    Intent(Out) :: KIN
  REAL(8),    Intent(Out) :: KOUT

  ! Locals:

! local scalar

      REAL(8) :: FC             !interpolation parameter between
                                !low and high pressure limits
      REAL(8) :: RKLOW          !low pressure reaction rate limit
                                !molecules/cm3/s
      REAL(8) :: RKHIGH         !high pressure limit
      REAL(8) :: BRN            !dry rate calc variable
      REAL(8) :: FAC1,FAC2,FAC3 !dry rate calc variables
      REAL(8) :: REUS1          !conversion factor
      REAL(8) :: RK0,RK2,RK3
      REAL(8) :: RH_nitrate     !RELATIVE HUMIDITY as a fraction
      REAL(8) :: RHd            !deliquescence point
      REAL(8) :: P1             !following scalars used in
      REAL(8) :: P2             !ammonium nitrate eqm calc
      REAL(8) :: P3
      REAL(8) :: KP1
      REAL(8) :: KP(1)
      INTEGER :: BLFLAG


! GAS PHASE REACTIONS:

! These first two reactions are not numbered as they are not in
! the STOCHEM code

! OH + 1,3-BUTADIENE = 1,3-BUTADIENE-PEROXY
! 1,3-BUTADIENE=C4H6

      RC_OH_BD=1.48E-11*DEXP(484/T)

! 1,3-BUTADIENE-PEROXY + NO = NO2 + GLYOX + HO2

      RC_NO_BDPEROXY=4.1E-12*DEXP(180/T)

! Reaction   1 OP + O2 + M = O3 + M
      RC_OP_O2  = 6.0E-34*(T/300)**(-2.3)
      if (rc_op_o2.lt.0.0)then
        print*,'rc_op_o2=',rc_op_o2,'air_den=',air_den,'o2=',o2
      endif
      RC_OP_O2 = RC_OP_O2*AIR_DEN*O2

! Reaction   5 O + NO + M = NO2 + M
      RKLOW  = 9.0E-32*((T/300)**(-1.5))*AIR_DEN
      RKHIGH = 3.0E-11
      FC     = 0.6
      BRN    = 0.75-1.27*DLOG10(FC)
      FAC1   = RKLOW/RKHIGH
      FAC2   = RKLOW/(1.0+FAC1)
      FAC3   = 1.0+((DLOG10(FAC1)/BRN)**2)
      RC_OP_NO = FAC2*FC**(1/FAC3)

! Reaction   7 O(1D) + M = O(3P) + M
      RC_OD  = 0.20948*3.2E-11*DEXP(70/T)+0.78084*1.8E-11*DEXP(110/T)
      RC_OD = RC_OD*AIR_DEN

! Reaction   8 O(1D) + H2O = OH + OH
      RC_OD_H2O  = 2.2E-10

! Reaction  11 NO + O3 = NO2 + O2

      RC_NO_O3 = 2.0E-12*DEXP(-1400/T)

! Reaction  12 NO2 + O3 = NO3 + O2

      RC_NO2_O3 = 1.2E-13*DEXP(-2450/T)

! Reaction  13 : OH + O3 = HO2 + O2

      RC_OH_O3 = 1.6E-12*DEXP(-940/T)

! Reaction  14 : HO2 + O3 = OH + O2 + O2

      RC_HO2_O3 = 1.1E-14*DEXP(-500/T)

! Reaction  15 NO + NO3 = NO2 + NO2

      RC_NO_NO3= 1.5E-11*DEXP(170/T)

! Reaction  16 : NO2 + O = NO + O2    - reaction added with JPL 1994

      RC_NO2_OP = 6.5E-12*DEXP(120/T)

! Reaction  17 : NO + HO2 = OH + NO2

      RC_NO_HO2 = 3.7E-12*DEXP(250/T)

! Reaction  19 NO2 + NO3 = NO + NO2 + O2

      RC2_NO2_NO3=4.5E-14*DEXP(-1260/T)

! Reaction  20 NO2 + NO3 + M = N2O5 + M

      RKLOW  = 2.2E-30*((T/300)**(-3.9))*AIR_DEN
      RKHIGH = 1.5E-12*((T/300)**(-0.7))
      FC     = 0.6
      BRN    = 0.75-1.27*DLOG10(FC)
      FAC1   = RKLOW/RKHIGH
      FAC2   = RKLOW/(1.0+FAC1)
      FAC3   = 1.0+((DLOG10(FAC1)/BRN)**2)
      RC_NO2_NO3 = FAC2*FC**(1/FAC3)

! Reaction  21 NO2 + OH + M = HNO3 + M

      RKLOW  = 2.6E-30*((T/300)**(-3.2))*AIR_DEN
      RKHIGH = 2.4E-11*((T/300)**(-1.3))
      FC     = 0.6
      BRN    = 0.75-1.27*DLOG10(FC)
      FAC1   = RKLOW/RKHIGH
      FAC2   = RKLOW/(1.0+FAC1)
      FAC3   = 1.0+((DLOG10(FAC1)/BRN)**2)
      RC_NO2_OH = FAC2*FC**(1/FAC3)

! Reaction 22 : NO2 + HO2 + M = HO2NO2 + M

      RKLOW  = 1.8E-31*((T/300)**(-3.2))*AIR_DEN
      RKHIGH = 4.7E-12*((T/300)**(-1.4))
      FC     = 0.6
      BRN    = 0.75-1.27*DLOG10(FC)
      FAC1   = RKLOW/RKHIGH
      FAC2   = RKLOW/(1.0+FAC1)
      FAC3   = 1.0+((DLOG10(FAC1)/BRN)**2)
      RC_NO2_HO2 = FAC2*FC**(1/FAC3)

! Reaction 23 : HO2NO2 + M = HO2 + NO2 + M

      RC_HO2NO2 = RC_NO2_HO2/(1.8E-27*DEXP(10900/T))

! Reaction 24 : OH + HO2NO2 = H2O + NO2 + O2

      RC_OH_HO2NO2 = 1.3E-12*DEXP(380/T)

! Reaction 27 NO3 + NO3 = NO2 + NO2 + O2

      RC_NO3_NO3=8.5E-13*DEXP(-2450/T)

! Reaction 29 N2O5 + M = NO2 + NO3 + M

      RKLOW  = 2.2E-03*((T/300)**(-4.4))*DEXP(-11080/T)*AIR_DEN
      RKHIGH = 9.7E+14*((T/300)**(+0.1))*DEXP(-11080/T)
      FC     = DEXP(-T/250)+DEXP(-1050/T)
      BRN    = 0.75-1.27*DLOG10(FC)
      FAC1   = RKLOW/RKHIGH
      FAC2   = RKLOW/(1.0+FAC1)
      FAC3   = 1.0+((DLOG10(FAC1)/BRN)**2)
      RC_N2O5 = FAC2*FC**(1/FAC3)

! Reaction 30 : HO2 + OH = H2O + O2

      RC_HO2_OH = 4.8E-11*DEXP(250/T)

! Reaction 31 : OH + H2O2 = H2O + HO2

      RC_OH_H2O2 = 2.9E-12*DEXP(-160/T)

! Reaction 32 NO3 + HO2 = HNO3 + O2

      RC_NO3_HO2 = 9.2E-13

! Reaction 33 : OH + H2 (+O2) = HO2 + H2O

      RC_OH_H2 = 5.5E-12*DEXP(-2000/T)

! Reaction 34 NO3 + HO2 = OH + NO2 + O2

      RC2_NO3_HO2 = 3.6E-12

! Reaction 35 OH + HNO3 = NO3 + H2O

      RK0  = 7.2E-15*DEXP(785/T)
      RK2  = 4.1E-16*DEXP(1440/T)
      RK3  = 1.9E-33*DEXP(725/T)
      FAC1 = RK3*AIR_DEN/RK2
      FAC2 = RK3*AIR_DEN/(1.0+FAC1)
      RC_OH_HNO3= RK0 + FAC2

! Reaction 36 HO2 + HO2 (+M,H2O) = H2O2 + O2

      RC_HO2_HO2=(2.3E-13*DEXP(600/T) + AIR_DEN*1.9E-33*DEXP(890/T) )  &
                   *( 1.0 + H2O*1.4E-21*DEXP(2200/T) )

! Reaction 39  OH + SO2 + M = HOSO2 + M
!(nb: RC_SO2_OH * OH concentration gives a dry oxd rate per second)

      RKLOW  = 3.0E-31*((T/300.0)**(-3.3))*AIR_DEN
      RKHIGH = 1.5E-12
      FC     = 0.6
      BRN    = 0.75-1.27*DLOG10(FC)
      FAC1   = RKLOW/RKHIGH
      FAC2   = RKLOW/(1.0+FAC1)
      FAC3   = 1.0+((DLOG10(FAC1)/BRN)**2)
      RC_SO2_OH=FAC2*FC**(1.0/FAC3)

! Reaction 40 : SO2 + CH3O2 = PRODUCTS

      RC_SO2_CH3O2 = 4.0E-17

! Reaction 43 NITRATE AEROSOL NAER - uptake of NOy
!(ie NO3-) by aerosol
!originally used   RC_NAER= 5.0E-06 for hno3 and n2o5 = NA
!as STOCHEM does
!7/09/2006 now trying 6.0e-6 for hno3 and 4.0e-4 for n2o5 as
!Salah used.
!       RC_NAER1=6.0E-6
!       RC_NAER2=4.0E-4

!new method from Dick 24/1/08
!keeping hno3=naer at 6.0e-6
!but n2o5=2naer to 6.0e-6 and a new reaction to do n2o5=2hno3
!which gives rate 1.0e-4 at 0 degrees and v small e-11 at 25 degrees

      RC_NAER1=6.0E-6
      RC_NAER2=6.0E-6

      RC_N2O5_2=EXP(49057.0_P64/T + LOG(9.07E-83_P64))

! Reaction 44 : OH + CH3OOH = CH3O2

      RC2_OH_CH3OOH = 1.90E-12*DEXP(190/T)

! Reaction 45 : OH + CH3OOH = OH + HCHO

      RC_OH_CH3OOH = 1.0E-12*DEXP(190/T)

! Reaction 59 : OH + CH4 = CH3O2 + H2O - Atkinson 1994 update

      RC_OH_CH4 = 7.44E-18*T*T*DEXP(-1361/T)

! Reaction 60 : NO + CH3O2 = CH3O + NO2
!                               = HCHO + HO2 + NO2

      RC_NO_CH3O2 = 4.2E-12*DEXP(180/T)

! Reaction 61 : CH3O2 + CH3O2 = CH3O + CH3O + O2 - change JPL
!                                  = HCHO + HCHO + HO2 + HO2
      FAC1  = 25*DEXP(-1165/T)
      FAC2  = FAC1/(1.0+FAC1)
      RC_CH3O2_2 = 9.1E-14*DEXP(416/T)*FAC2

! Reaction 62 : CH3O2 + CH3O2 = CH3OH + HCHO + O2

      RC_CH3O2_CH3O2 = 9.1E-14*DEXP(416/T)*(1.0-FAC2)

! Reaction 65 : CH3O2 + HO2 = CH3OOH + O2  - change JPL 1994

      RC_CH3O2_HO2 = 3.8E-13*DEXP(800/T)

! Reaction 66 : OH + HCHO = HCO + H2O
!                              = HO2 + CO + H2O

      RC_OH_HCHO = 1.0E-11

! Reaction 67 : NO3 + HCHO = HO2 + CO + HNO3

      RC_NO3_HCHO = 5.8E-16

! Reaction 70 : OH + CO = H + CO2
!                            = HO2 + CO2

      RC_OH_CO = 1.5E-13*(1.0+0.6*(AIR_DEN/2.55E+19))


! Reaction  75 : OH + CH3CHO = CH3CO + H2O   - JPL 1994
!                                = CH3COO2 + H2O

      RC_OH_CH3CHO = 6.0E-12*DEXP(250/T)

! Reaction  77 : CH3COO2 + NO2 + M = CH3COO2NO2 (PAN) + M
      RKLOW  = 2.7E-28*((T/300)**(-7.1))*AIR_DEN
      RKHIGH = 1.2E-11*((T/300)**(-0.9))
      FC     = 0.3
      BRN    = 0.75-1.27*DLOG10(FC)
      FAC1   = RKLOW/RKHIGH
      FAC2   = RKLOW/(1.0+FAC1)
      FAC3   = 1.0+((DLOG10(FAC1)/BRN)**2)
      RC_CH3COO2_NO2 = FAC2*FC**(1/FAC3)

! Reaction  78 : CH3COO2NO2 + M = CH3COO2 + NO2 + M
!   nb CH3COO2NO2=PAN

      RKLOW  = 5.5E-03*DEXP(-12064/T)*AIR_DEN
      RKHIGH = 3.9E+16*DEXP(-13628/T)
      FC     = 0.3
      BRN    = 0.75-1.27*DLOG10(FC)
      FAC1   = RKLOW/RKHIGH
      FAC2   = RKLOW/(1.0+FAC1)
      FAC3   = 1.0+((DLOG10(FAC1)/BRN)**2)
      RC_PAN = FAC2*FC**(1/FAC3)

! Reaction  79 : CH3COO2 + NO = CH3COO + NO2
!                                 = CH3 + CO2 + NO2
!                                 = CH3O2 + CO2 + NO2
      RC_CH3COO2_NO = 2.0E-11

! Reaction 80A : CH3O2 + CH3COO2 = CH3O + CH3COO + O2
!                                    = HCHO + HO2 + CH3 + CO2 +O2
!                                    = HCHO + HO2 + CH3O2 + CO2 +O2
! THIS IS RC(80)IN BACKIT
      FAC1  = 4.4E+05*DEXP(-3910/T)
      FAC2  = 1.0/(1.0+FAC1)
      RCA_CH3O2_CH3COO2 = 5.1E-12*DEXP(272/T)*FAC2

! Reaction 80B : CH3O2 + CH3COO2 = HCHO + CH3COOH + O2
! ACETIC ACID CONVERTED TO FORMALDEHYDE IN ONE STEP
! this is RC(74) in backit despite being called reaction 80B

      RC_CH3O2_CH3COO2 = 5.1E-12*DEXP(272/T)*(1.0-FAC2)

! Reaction  91 : CH3COO2 + CH3COO2 = CH3COO + CH3COO + O2
!                                      = CH3 + CH3 + CO2 + CO2 + O2
!                                      = CH3O2 + CH3O2 + CO2 + CO2 + O2

      RC_CH3COO2_CH3COO2 = 2.8E-12*DEXP(530/T)

! Reaction  98 : OH + PAN = CH3COO2 + HNO3

      RC_OH_PAN = 9.5E-13*DEXP(-650/T)

! Reaction 109 : OH + C2H4 (+O2) = HOC2H4O2 (ch2o2c)

      RKLOW  = 7.0E-29*((T/300)**(-3.1))*AIR_DEN
      RKHIGH = 9.0E-12
      FC    = 0.7
      BRN     = 0.75-1.27*DLOG10(FC)
      FAC1  = RKLOW/RKHIGH
      FAC2  = RKLOW/(1.0+FAC1)
      FAC3  = 1.0+((DLOG10(FAC1)/BRN)**2)
      RC_OH_C2H4 = FAC2*FC**(1/FAC3)

! Reaction 110 : HOC2H4O2 + NO = HOC2H4O + NO2
! (CH2O2C+NO)                    = HCHO + HO2 + HCHO + NO2

      RC_CH2O2C_NO = 9.0E-12

! Reaction 111 : CH3O2 + HOC2H4O2 = CH3O + HOC2H4O + O2 - JPL
!                            = HCHO + HO2 + HCHO + HCHO + HO2 + O2
! (diff way of putting ch3o2+ch2o2c really)

      RC_CH3O2_CH2O2C = 1.0E-12

! Reaction 112 : O3 + C2H4 = HCHO + 0.47 CH2O2 + 0.31 CO
!            + 0.22 CO2 + 0.31 H2O + 0.13 H2 + 0.20 HO2

      RC_O3_C2H4 = 1.2E-14*DEXP(-2630/T)

! Reaction 123 : O3 + C3H6 = HCHO + 0.30 CH4 + 0.40 CO
!        + 0.60 CO2 + 0.28 OH + 0.12 CH3OH + 0.30 HO2 + 0.58 CH3O2

      RC2_O3_C3H6 = 4.0E-15*DEXP(-1900/T)


! Reaction 124 : O3 + C3H6 = CH3CHO + 0.24 H2 + 0.58 CO
!        + 0.42 CO2 + 0.58 H2O +0.18 HO2

      RC_O3_C3H6 = 2.6E-15*DEXP(-1900/T)

! Reaction 125 : OH + C3H6 (+O2) = HOC3H6O2   - JPL

      RKLOW  = 8.0E-27*((T/298)**(-3.5))*AIR_DEN
      RKHIGH = 3.0E-11
      FC     = 0.5
      BRN    = 0.75-1.27*DLOG10(FC)
      FAC1   = RKLOW/RKHIGH
      FAC2   = RKLOW/(1.0+FAC1)
      FAC3   = 1.0+((DLOG10(FAC1)/BRN)**2)
      RC_OH_C3H6 = FAC2*FC**(1/FAC3)

! Reaction 126 : HOC3H6O2 + NO = HOC3H6O + NO2
!                                  = CH2OH + CH3CHO + NO2
!                                  = HCHO + HO2 + CH3CHO + NO2

      RC_CH3CHX_NO = 9.0E-12


! Reaction 127 : CH3O2 + HOC3H6O2 = CH3O + HOC3H6O + O2  - JPL
!                          = HCHO + HO2 + CH2OH + CH3CHO + O2
!                          = HCHO + HO2 + HCHO + HO2 + CH3CHO + O2

      RC_CH3O2_CH3CHX = 1.0E-14

! Reaction 128 : O3+C5H8=MVK+0.78*CO+0.22*CH2OO+0.27*HO2+0.27*OH

      RC_O3_C5H8 = 7.86E-15*DEXP(-1913/T)


! Reaction 129 : O3+MVK=MGLYOX+0.76*CO+0.24*CH2OO+0.36*HO2+0.36*OH

      RC_O3_MVK = 7.56E-16*DEXP(-1521/T)

! Reaction 130 : CH2OO + NO = NO2 + HCHO

      RC_CH2OO_NO = 1.0E-14

! Reaction 131 : CH2OO + NO2 = NO3 + HCHO

      RC_CH2OO_NO2 = 1.0E-15

! Reaction 132 : CH2OO + H2O = HCOOH + H2O

      RC_CH2OO_H2O = 5.8E-17

! Reaction 133 : CH2OO + SO2 = SA + HCHO

      RC_CH2OO_SO2 = 7.0E-14


! Reaction 203 : NO3 + C2H4 = C2H4NO3           - JPL
!                               = CH2(NO3)CHO + HO2

      RC_NO3_C2H4 = 4.88E-18*T*T*DEXP(-2282/T)

! Reaction 205 : NO3 + C3H6 = C3H6NO3         - JPL
!                               = CH3CH(NO3)CHO + HO2

      RC_NO3_C3H6 = 4.59E-13*DEXP(-1156/T)


! Reaction 208 : NO3 + CH3CHO = CH3CO + HNO3         JPL
!                                 = CH3COO2 + HNO3

      RC_NO3_CH3CHO = 1.4E-12*DEXP(-1900/T)

! Reaction 209 : NO3 + C5H8 = (NO3)C4H6CHO + HO2

      RC_NO3_C5H8 = 3.03E-12*DEXP(-446/T)


! Reaction 213 : NO + TOLP1 = MEMALD + GLYOX + HO2 + NO2 - JPL

      RC_NO_TOLP1 = 4.0E-12

! Reaction 230 : OH + (C6H4)(CH3)2 = HO(C6H4)(CH3)2
! {OXYL=(C6H4)(CH3)2}                     = HO(C6H4)(CH3)2O2
! {OXYL1=HO(C6H4)(CH3)2O2}

      RC_OH_OXYL = 1.37E-11

! Reaction 231 : HO(C6H4)(CH3)2(O2) + NO
!                  = HO(C6H4)(CH3)2(O) + NO2
!                  = CH3COCHO + HO2 + CH3COCH=CHCHO + NO2

      RC_OXYL1_NO = 4.0E-12

! Reaction 232 : OH + CH3COCH=CHCHO = CH3COCH(OH)CHCHO
! MEMALD+OH=MEMALD1                      = CH3COCH(OH)CH(O2)CHO

      RC_MEMALD_OH = 5.6E-11

! Reaction 233 : CH3COCH(OH)CH(O2)CHO + NO
!                           = CH3COCH(OH)CH(O)CHO + NO2
!                           = HO2 + CH3COCHO + CHOCHO + NO2

      RC_MEMALD1_NO = 9.0E-12

! Reaction 234 : OH + TOLUENE = TOLP1
! TOLUENE=C7H8

      RC_OH_TOLUEN = 5.96E-12

! Additional reactions of Toluen with OH to represent formation of SOA
! SOURCE: DICK DERWENT reaction(78) toleune+oh=RA16O2
! reaction(79) AROPREC + OH where aroprec is formed as a by-product from OH + toluene and 
! removed by reversible take up into the aerosol and by oxidation by OH

      RC2_OH_TOLUEN = 1.81E-12*DEXP(338/T)

      RC_OH_AROPREC = 1.0E-11

      KIN=0.0558
      KOUT=2.065*DEXP(-10282/T)

! Reaction 241 : MEMALDIAL1 + CH3O2 = 2HO2 + HCHO + MGLYOX + GLYOX

      RC_MEMALD1_CH3O2 = 1.0E-13

! Reaction 242 : RO2IP1 + CH3O2 = 2HO2 + HCHO + MVK - MCM 1996

      RC_RO2IP1_CH3O2 = 5.0E-13

! Reaction 243 : RO2IP2 + CH3O2 = 2HO2 + HCHO + MGLYOX  - MCM 1996

      RC_RO2IP2_CH3O2 = 2.0E-12

! Reaction 244 : RO2IP1 + HO2 = ISOPOOH + O2

      RC_RO2IP1_HO2 = 2.45E-13*DEXP(1250/T)

! Reaction 245 : ISOPOOH + OH = MVK + HCHO + OH

      RC_ISOPOOH_OH = 4.20E-11

! Reaction 246 : RO2IP2 + HO2 = MVKOOH + O2

      RC_RO2IP2_HO2 = 2.23E-13*DEXP(1250/T)

! Reaction 247 : MKVOOH + OH = MGLYOX + HCHO + OH

      RC_MVKOOH_OH = 5.77E-11

! Reaction 248 : MGLYOX + OH = CH3COO2 + CO

      RC_MGLYOX_OH = 1.72E-11

! Reaction 249 : GLYOX + OH = HO2 + CO + CO

      RC_GLYOX_OH = 1.14E-11

! Reaction 251 : OH + C5H8 = (HO)C5H8O2 + H2O

      RC_OH_C5H8 = 2.54E-11*DEXP(410/T)

! Reaction 252 : (HO)C5H8O2 + NO = MVK + HO2 + HCHO + NO2 - MCM 1996

      RC_RO2IP1_NO = 2.08E-12*DEXP(180/T)

! Reaction 253 : OH + MVK = (HO)MVKO2 + H2O
! (HO)MVKO2=RO2IP2

      RC_OH_MVK = 4.13E-12*DEXP(452/T)

! Reaction 254 : (HO)MVKO2 + NO = CH3COCHO + CH2O + HO2 + NO2
! (HO)MVKO2=RO2IP2

      RC_RO2IP2_NO = 2.46E-12*DEXP(180/T)

! Reactions for formation of SOA (no numbers as not from the original STOCHEM
! scheme to which the earlier reaction numbers refer):

! Reaction soa1 : OH + C10H16BIO = precursor

      RC_OH_C10H16BIO = 1.21E-11*DEXP(444/T)

! Reaction soa2 : O3 + C10H16BIO = precursor

      RC_O3_C10H16BIO = 1.01E-15*DEXP(-732/T)

! Reaction soa3 : NO3 + C10H16BIO = precursor

      RC_NO3_C10H16BIO = 1.19E-12*DEXP(490/T)

! Reaction soa4 : OH + precursor = products

      RC_OH_SOAP = 1.0E-11

! Reaction soa5 : Scavenging of the SOA precursor

      RC_SCAV = 4.0E-4

! NH3 + HNO3 TO NH4NO3 or NH4NO3 TO NH3 + HNO3
! this rate is dependent on humidity
! this is not a rate coefficient but an equilibrium one

      RHd = exp((618.3/T)-2.551)

! This is setting the RH as a fraction as opposed to a %
      RH_nitrate=RH/100.0

      RH_nitrate=min(RH_nitrate,0.9999d0)

      IF (RH_nitrate < RHd) THEN       ! NH4NO3 is below deliquescence point
        KP(1)=(T**(-6.025))*EXP(118.87-24084/T)
      ELSE                     ! NH4NO3 is above deliquescence point
        P1=exp(-135.94+(8763/T))*T**19.12
        P2=exp(-122.65+(9969/T))*T**16.22
        P3=exp(-182.61+(13875/T))*T**24.46
        KP1=(T**(-6.025))*EXP(118.87-24084/T)
        KP(1)=(P1-(P2*(1.0-RH_nitrate))+P3*((1.0-RH_nitrate)**2))*((1.0-RH_nitrate)**1.75)*KP1
      ENDIF

      Kp_RT2=5.39E25*KP(1)/(T*T)

! HONO reactions (not in STOCHEM)
! ORIGINAL SET TRIED  - NOW REPLACED - SEE BELOW
! OH + HONO = H2O + NO2 rate from DeMore et al 1997 via Salisbury
! et al 2001 JGR vol 106 no.D12

!      RC_OH_HONO=1.8E-11*DEXP(-390.0/T)
!      RC_OH_HONO=0.0

! OH + NO (+M) = HONO (+M) rate origin as above

!      RC_OH_NO=4.87E-12
!      RC_OH_NO=0.0

! Thermal heterogenous reaction with water vapour
! ref Mike Jenkin & Ozone in the United Kingdom 4th PORG report
! 1997 page 22
! 2NO2(g) + H2O(ads) = HNO3(ads) + HONO(g)
! this rate is however only a fudge saying NO2 = HONO

!      RC_NO2=5.6E-6*100.0/ZI
!      RC_NO2=0.0

!NB we are not using this method of generating HONO now - we are
!directly emitting 1% NO as HONO instead
!FROM 9/02 NOW NOT DOING ANY HONO CHEM OR EMISSIONS

!15/01/2008 trying hono chemistry again to try to resolve HNO3
! and nitrate aerosol problems

! OH + HONO = H2O + NO2 rate from Dick but similar to old rate above
! so may revert to that after initial testing

      RC_OH_HONO=6.0E-12

! OH + NO (+M) = HONO (+M) rate from Dick

      RC_OH_NO=9.8E-12

! heterogeneous source of HONO - if in the boundary layer
! BLFLAG=1 IN BL 0 NOT IN BL
! reverting back to Harrison & Kitto 0.5e-3 rather then the
! 0.5e-2 as Dick suggested we might need (nb this is approx what
! we used to use above)

      RC_NO2=BLFLAG*0.5E-3/ZI

! H2S Reactions (not in STOCHEM)

! H2S + OH = HS + H2O
! Atkinson et al, ACP, (2004)
! HS + O2 = HO2 + SO2
! Chernysheva et al., Academy of Sciences USSR (1989)

  RC_OH_H2S = 6.1E-12*DEXP(-80.0/T)


! AQUEOUS PHASE REACTIONS:

      REUS1 = ((1.0/T)-(1.0/298.0))

! HSO3(-) + H2O2(aq) = H(+) + SO4(2-) + H2O
! nb:  The H+ dependence is in backward Euler code
! Martin & Damschen (1981), Bower et al., (1991)

      RC_HSO3_H2O2= 5.2E+6*DEXP(-3650*REUS1)

! HSO3(-) + O3 = H(+) + SO4(2-) + O2
! Maahs (1983)

      RC_HSO3_O3= 4.2E+5*DEXP(-4131*REUS1)

! SO3(2-) + O3 = SO4(2-) + O2
! Maahs (1983)

      RC_SO3_O3= 1.5E+9*DEXP(-996*REUS1)

End Subroutine CHEMCO

!-------------------------------------------------------------------------------------------------------------

Subroutine PH3(KE_SO2, KE_HSO3, KE_NH3, KE_CO2, KE_H2O, KE_HNO3,        &
               KH_O3, KH_H2O2, KH_SO2, KH_NH3, KH_CO2, KH_HNO3,         &
               CLDWAT, AIR_DEN, T, CT_SO2, CT_NH3, AQ_SO4, CT_HNO3, HP)
! Solves H+ by a bissection routine.

! This is a Fortran 90 version of the subroutine PH3 in NAME V8.08.

  Implicit None
  ! Argument list:
  REAL(8),    Intent(In)  :: KE_SO2     !} Equilibrium coefficients.
  REAL(8),    Intent(In)  :: KE_HSO3    !}
  REAL(8),    Intent(In)  :: KE_NH3     !}
  REAL(8),    Intent(In)  :: KE_CO2     !}
  REAL(8),    Intent(In)  :: KE_H2O     !}
  REAL(8),    Intent(In)  :: KE_HNO3    !}
  REAL(8),    Intent(In)  :: KH_O3      !] Henry's law coefficients.
  REAL(8),    Intent(In)  :: KH_H2O2    !]
  REAL(8),    Intent(In)  :: KH_SO2     !]
  REAL(8),    Intent(In)  :: KH_NH3     !]
  REAL(8),    Intent(In)  :: KH_CO2     !]
  REAL(8),    Intent(In)  :: KH_HNO3    !]
  REAL(8),    Intent(In)  :: CLDWAT     ! Cloud liquid water content.
  REAL(8),    Intent(In)  :: AIR_DEN    ! Air density (in molec/cm3).
  REAL(8),    Intent(In)  :: T          ! Absolute temperature.
  REAL(8),    Intent(In)  :: CT_SO2     !} Initial concentrations (in moles/litre).
  REAL(8),    Intent(In)  :: CT_NH3     !}
  REAL(8),    Intent(In)  :: CT_HNO3    !}
  REAL(8),    Intent(In)  :: AQ_SO4     !}
  REAL(8),    Intent(Out) :: HP         ! Calculated ph (H+).
  ! Locals:

! The method is one of
! trial and error search to find a minima in the modulus of the
! imposed minus the returned value of [H+], then use bissection to
! obtain the value to the supplied tolerance.

! local scalars

      INTEGER :: I,N          !N=counter for failures
      INTEGER :: COUNT
      INTEGER :: MAXCOUNT     !maximum iteration length
      REAL(8) :: XTOL         !fractional tolerance
      REAL(8) :: HP1          !first guess

! local arrays

      REAL(8) :: X(20)
      REAL(8) :: Y(20)
      REAL(8) :: X2(5)
      REAL(8) :: Y2(5)

      XTOL=0.001
      MAXCOUNT=20
      N=0

! set at largest expected [H+] *25:

      HP=1.0E-2
      COUNT=0

! find first minima

      IF(CLDWAT.GT.1.0E-10)THEN
        DO I=1,20
          HP=HP/5.0
          X(I)=HP


          CALL PHCALC(X(I),HP1,CT_SO2,CT_NH3,CT_HNO3,AQ_SO4,CLDWAT,      &
            AIR_DEN,T,KE_SO2,KE_HSO3,KE_NH3,KE_CO2,KE_H2O,KH_SO2,KH_CO2)

          Y(I)=DABS(X(I)-HP1)
          IF(I.GT.2)THEN
            IF(Y(I).GT.Y(I-1)) GOTO 20
          ENDIF
        ENDDO
        GOTO 30          ! shouldn't reach here

   20   CONTINUE
        X2(1)=X(I-2)
        X2(2)=X(I-1)
        X2(3)=X(I)
        Y2(1)=Y(I-2)
        Y2(2)=Y(I-1)
        Y2(3)=Y(I)

        IF(Y2(2).GT.Y2(1).OR.Y2(2).GT.Y2(3)) GOTO 30

! Iterate bisection until the interval is within the tolerance

   25   CONTINUE
        COUNT=COUNT+1
        X2(4)=(X2(1)+X2(2))/2.0


        CALL PHCALC(X2(4),HP1,CT_SO2,CT_NH3,CT_HNO3,AQ_SO4,CLDWAT,     &
          AIR_DEN,T,KE_SO2,KE_HSO3,KE_NH3,KE_CO2,KE_H2O,KH_SO2,KH_CO2)

        Y2(4)=DABS(X2(4)-HP1)
        X2(5)=(X2(3)+X2(2))/2.0


        CALL PHCALC(X2(5),HP1,CT_SO2,CT_NH3,CT_HNO3,AQ_SO4,CLDWAT,     &
          AIR_DEN,T,KE_SO2,KE_HSO3,KE_NH3,KE_CO2,KE_H2O,KH_SO2,KH_CO2)

        Y2(5)=DABS(X2(5)-HP1)

        IF(Y2(2).LT.Y2(4).AND.Y2(2).LT.Y2(5)) THEN
          X2(1)=X2(4)
          X2(3)=X2(5)
          Y2(1)=Y2(4)
          Y2(3)=Y2(5)
        ELSE IF(Y2(4).LT.Y2(5)) THEN
          X2(3)=X2(2)
          X2(2)=X2(4)
          Y2(3)=Y2(2)
          Y2(2)=Y2(4)
          IF((X2(1)-X2(2))/X2(2).LT.XTOL) GOTO 27
        ELSE
          X2(1)=X2(2)
          X2(2)=X2(5)
          Y2(1)=Y2(2)
          Y2(2)=Y2(5)
          IF((X2(1)-X2(2))/X2(2).LT.XTOL) GOTO 27
        ENDIF
        IF(COUNT.GT.MAXCOUNT) GOTO 27
        GOTO 25

   27   CONTINUE
        HP=X2(2)
        GOTO 40

   30   CONTINUE
        N=N+1
        HP=1.0E-6                          ! In case of trouble.
      ENDIF

   40 CONTINUE

End Subroutine PH3

!-------------------------------------------------------------------------------------------------------------

Subroutine PHCALC(HP, HP1, CT_SO2, CT_NH3, CT_HNO3, AQ_SO4,   &
                  CLDWAT, AIR_DEN, T,                         &
                  KE_SO2, KE_HSO3, KE_NH3, KE_CO2, KE_H2O,    &
                  KH_SO2, KH_CO2)
! Returns [H+].

! This is a Fortran 90 version of the subroutine PHCALC in NAME V8.08.

  Implicit None
  ! Argument list:
  REAL(8),    Intent(In)  :: HP       ! Trial pH in.
  REAL(8),    Intent(Out) :: HP1      ! pH out.
  REAL(8),    Intent(In)  :: CT_SO2   !} Initial concs in moles/litre of air.
  REAL(8),    Intent(In)  :: CT_NH3   !}
  REAL(8),    Intent(In)  :: CT_HNO3  !}
  REAL(8),    Intent(In)  :: AQ_SO4   ! Initial so4 conc in moles/litre aqueous.
  REAL(8),    Intent(In)  :: CLDWAT   ! Cloud liquid water content.
  REAL(8),    Intent(In)  :: AIR_DEN  ! Air density (molec/cm3).
  REAL(8),    Intent(In)  :: T        ! Absolute temperature.
  REAL(8),    Intent(In)  :: KE_SO2   !} Equilibrium coefficients.
  REAL(8),    Intent(In)  :: KE_HSO3  !}
  REAL(8),    Intent(In)  :: KE_NH3   !}
  REAL(8),    Intent(In)  :: KE_CO2   !}
  REAL(8),    Intent(In)  :: KE_H2O   !}
  REAL(8),    Intent(In)  :: KH_SO2   !] Henry's law coefficients.
  REAL(8),    Intent(In)  :: KH_CO2   !]
  ! Locals:

! local scalars

      REAL(8) :: A1
      REAL(8) :: REUS     !conversion factor
      REAL(8) :: FSO2     !extra factor for aq so2
      REAL(8) :: AQ_SO2   !dissolved concentrations in moles/litre
      REAL(8) :: AQ_NH3
      REAL(8) :: AQ_HNO3
      REAL(8) :: AQ_CO2
      REAL(8) :: AQ_NO3
      REAL(8) :: AQ_HSO3
      REAL(8) :: AQ_SO3
      REAL(8) :: AQ_NH4
      REAL(8) :: AQ_HCO3


! FSO2 is an extra factor to increase the solubility of so2 -
! so2 is more soluble than its Henry's law coeff suggests due to
! dissociation of the dissolved species

      A1=KE_SO2/HP
      FSO2=1.0+A1+A1*KE_HSO3/HP

! Obtain the dissolved species concentrations in (moles/litre)
! HNO3, SO2, NH3 from eqm equations

      REUS=RUGC*T*1.0E6/(AVOGAD*P_REF)

! HNO3(aq)+NO3

      AQ_HNO3=CT_HNO3/CLDWAT

! AQ_SO2 is [SO2(aq)] + [HSO3-] + [SO3--]

      AQ_SO2=CT_SO2/(CLDWAT+P_REF/(KH_SO2*RUGC*T*1.0E3*FSO2))

! NH3+NH4 (Aq)

      AQ_NH3=CT_NH3/CLDWAT

! CO2 (AQ) fixed conc of co2 is used here - check it

      AQ_CO2=KH_CO2*360.0E-6*AIR_DEN*REUS

! NO3  - HNO3 is assumed completly dissociated.

      AQ_NO3=AQ_HNO3

! HSO3

      AQ_HSO3=(AQ_SO2/FSO2)*KE_SO2/HP

! SO3

      AQ_SO3=AQ_HSO3*KE_HSO3/HP

! NH4

      AQ_NH4=AQ_NH3/(1+KE_H2O/(HP*KE_NH3))

! HCO3

      AQ_HCO3=AQ_CO2/(1+HP/KE_CO2)

! Calculate HP using:
!        [H+]=[NO3-]+[HSO3-]+2[SO3--]+2[SO4--]+[HCO3-]-[NH4+]


      HP1 = AQ_HSO3 + 2.0*AQ_SO3 + 2.0*AQ_SO4 + AQ_NO3 + AQ_HCO3 - AQ_NH4

End Subroutine PHCALC

!-------------------------------------------------------------------------------------------------------------

End Module NameChemistrySchemeModule
