! Module: Chemistry Module

Module ChemistryModule

! This module contains chemical constants and routines for the chemistry calculations.

! Note that some changes to the format of the NAME V8.08 code have been necessary
! to run under Fortran 90 (comment lines, continuation lines, End statements, etc.).
! Locally declared parameters (e.g. Avogadro's number) have been removed from each
! module procedure and are declared in the module or as global parameters instead.
! Any material changes to the NAME V8.08 code are noted at the start of each routine.

! Module call tree
! ----------------

! PreInitChemOpts
!
! InitChemOpts
!
! SetUpChemistryDefns---InitChemistryDefn
!
! PreInitChemistryState
!
! InitChemistryState
!
! WriteRestartFileChemistry
!
! ReadRestartFileChemistry
!
! SplitMassiveParticles
!
! UpdateChemistry---(SetBackgroundFieldsFromSTOCHEM---(FindSTOCHEMGridbox
!                   (InitChemistryGridboxes           (GetMetForChemistry
!                   (GetMetForChemistry
!                   (ConvertMassToConc
!                   (ConvertConcToMass
!                   (ConvertChemFieldToConc
!                   (ConvertConcToChemField
!                   (UpdateParticleMasses
!                   (Individual chemistry schemes !$$ 

!-------------------------------------------------------------------------------------------------------------

Use TimerModule
Use ServiceModule
Use FlowsModule, Only: Flows_, FlowMemory_,                    &
                       ConvertToZ, GetAttrib, ResetFlowMemory, &
                       A_Flow, A_Cloud,                        &
                       Mets_,                                  &
                       Flow_, Cloud_, Rain_,                   &
                       Surface_, Soil_, Plant_
Use SpeciesModule
Use ParticleModule
Use PuffModule
Use EulerianInterfaceModule
Use NameChemistrySchemeModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private

! Items from this module which need to be made available outside the module:
Public :: ChemOpts_                 ! Options for use with the chemistry scheme.
Public :: ChemistryDefn_            ! Information defining a chemistry scheme.
Public :: GridboxState_             ! Identifies particles in a chemistry grid box.
Public :: ChemistryState_           ! State of a chemistry calculation.
Public :: PreInitChemOpts           ! Pre-initialises the chemistry options
Public :: InitChemOpts              ! Initialises the chemistry options.
Public :: SetUpChemistryDefns       ! Sets up a definition of a chemistry scheme.
Public :: PreInitChemistryState     ! Pre-initialises the state of a chemistry
                                    ! calculation by setting the pointers to null.
Public :: InitChemistryState        ! Initialises the state of a chemistry calculation.
Public :: WriteRestartFileChemistry ! Writes chemistry fields to a restart file.
Public :: ReadRestartFileChemistry  ! Reads chemistry fields from a restart file.
Public :: SplitMassiveParticles     ! Particle splitting routine for massive particles.
Public :: UpdateChemistry           ! Call to the chemistry scheme (at each sync time).
Public :: ChemistryModuleTimerSummary
Public :: ChemistryModuleTimerInitialise

!-------------------------------------------------------------------------------------------------------------

! local parameters (MOVED FROM CHEMISTRY SUBROUTINES)
! Note that the molecular masses of dry air and water are specified in the collection
! of Physics parameters (in g/mole). We use these values in the chemistry scheme but
! new parameters are set (in kg/mole) for consistency with the NAME chemistry code.
! We also declare a copy of the Avogadro constant here - for back-compatibility with
! the NAME chemistry code. Note that all these new parameters are REAL(8).
! $$ Change "AVOGAD" in NAME code to "Avogadro"? $$ DBLE ???

Real(8), Parameter :: AVOGAD  = (Avogadro)          ! Copy of Avogadro number.
Real(8), Parameter :: RMM_AIR = (MoleMassAir)/1.0E3 ! Molecular mass of dry air (kg/mole).

!-------------------------------------------------------------------------------------------------------------

Type :: ChemOpts_ ! Chemistry options.
  Logical                      :: Initialised  ! Indicates that the chemistry options have been set.
  Character(MaxFileNameLength) :: ChemFolder   ! Name of folder storing the STOCHEM background fields.

End Type ChemOpts_

!-------------------------------------------------------------------------------------------------------------

Type :: ChemistryDefn_ ! Information defining a chemistry scheme (grids, species, etc.).
  Logical                  :: Initialised
  Logical                  :: ChangeUnitsParticle (MaxSpecieses)
  Logical                  :: ChangeUnitsField (MaxSpecieses)
  Character(MaxCharLength) :: Name
  Integer                  :: nParticleSpecies
  Integer                  :: nFieldSpecies
  Integer                  :: nBackFields
  Real(8)                  :: ParticleSpeciesMolWeights (MaxSpecieses)
  Real(8)                  :: ParticleSpeciesInvHalfLife (MaxSpecieses)
  Character(MaxCharLength) :: FieldSpeciesNames         (MaxSpecieses)
  Character(MaxCharLength) :: FieldSpeciesCategories    (MaxSpecieses)
  Real(8)                  :: FieldSpeciesMolWeights    (MaxSpecieses)
  Real(8)                  :: FieldSpeciesInvHalfLife   (MaxSpecieses)
  Character(MaxCharLength) :: BackFieldSpeciesNames     (MaxSpecieses)
  Real(8)                  :: BackFieldSpeciesMolWeights(MaxSpecieses)
  Integer                  :: nChemSpecieses
  Integer                  :: iChem2Particle(MaxSpecieses) 
  Integer                  :: iChem2Field   (MaxSpecieses)
  Integer                  :: iZCoordMagl
                              ! $$ might want to reconsider where/how nChemSpecieses, iChem2Particle and
                              ! iChem2Field are stored (e.g. multiple chemistry schemes, other uses of
                              ! indexing between species lists).
  ! Initialised                :: Indicates that the chemistry definition has been initialised.
  ! ChangeUnitsParticle        :: Indicates whether the input emission on particles needs changing from Bq to mass enable chemistry
  ! ChangeUnitsField           :: Indicates whether the input emission on field needs changing from Bq to mass enable chemistry
  ! Name                       :: Name of the chemistry definition.
  ! nParticleSpecies           :: Number of species held on particles.
  ! nFieldSpecies              :: Number of species held on fields.
  ! ParticleSpeciesMolWeights  :: Molecular weights of species held on particles.
  ! ParticleSpeciesInvHalfLife  :: Reciprocal of the half life of species held on particles.
  ! FieldSpeciesNames          :: Names of species held on fields.
  ! FieldSpeciesCategories     :: Categories of species held on fields,
  !                               used to indicate if a species needs initialising from a file
  ! FieldSpeciesMolWeights     :: Molecular weights of species held on fields.
  ! FieldSpeciesInvHalfLife    :: Reciprocal of the half life of species held on fields.
  ! BackFieldSpeciesNames      :: Names of species read in as a background field on chemistry grid.
  !                            :: Currently these are from STOCHEM but in future could be from other sources
  ! BackFieldSpeciesMolWeights :: Molecular weights of background species held on chemistry grid.
  ! nChemSpecieses             :} Total number of chemistry species and their indices in specieses
  ! iChem2Particle             :} held on particles and in specieses held on fields.
  ! iChem2Field                :} 0 indicates species not held on particles / fields.
  ! iZCoordMagl                :: Index of magl coord system. $$ Not nec best place to store,
  !                               esp if ChemistryDefn initialised/setup earlier than at present.
End Type ChemistryDefn_

!-------------------------------------------------------------------------------------------------------------

Type :: GridboxState_ ! Identifies the particles in a chemistry grid box.
  Integer :: nOld     ! Old number of particles in the chemistry grid box (at the previous sync time).
  Integer :: n        ! Number of particles in the chemistry grid box (at the current sync time).
  Integer :: k        ! Position in the ParticleIndex array of the latest particle
                      ! to be assigned for this grid box (ranges from N+1 to N+n,
                      ! where N is the sum of the n over all the previous grid boxes).
  Real(8) :: DeltaO3  ! Increment to relax O3 field in the chemistry grid box (in molecules/cm3).
End Type GridboxState_

!-------------------------------------------------------------------------------------------------------------

Type :: ChemistryState_ ! State of a chemistry calculation.
  ! If changing this type, remember restart read and write routines.
  Logical                      :: Initialised
  Logical                      :: InitialiseBackFields ! initialise these in ChemistryUpdate
  Integer,             Pointer :: ParticleIndex(:)
  Type(GridboxState_), Pointer :: GridboxState(:,:,:)
  Real(8),             Pointer :: Volume(:,:,:)
  Real(8),             Pointer :: Mass(:,:,:,:)
  Real(8),             Pointer :: OldMass(:,:,:,:)
  Real(8),             Pointer :: ChemBackField(:,:,:,:)
  ! Initialised   :: Indicates that the chemistry state has been initialised.
  ! ParticleIndex :: Indices of the particles (grouped by gridbox).
  ! GridboxState  :: Used with ParticleIndex to identify particles in each grid box.
  ! Volume        :: Volume of each grid box in cm3 (constant throughout model run).
  ! Mass          :: Mass (in g) in each grid box (for species held on particles).
  ! OldMass       :: Initial mass (in g) in each grid box (for species held on particles).
  ! ChemBackField :: Chemistry field concentrations (in g/m3) for background species read in.
End Type ChemistryState_

!-------------------------------------------------------------------------------------------------------------

! Timers

Type(Timer_), Save :: ChemistryLoopTimer

Logical :: TimersInitialised = .false.

!-------------------------------------------------------------------------------------------------------------

Contains

Subroutine ChemistryModuleTimerInitialise(TimerOpts)

  Implicit None

  Type(TimerOpts_) :: TimerOpts

  ! Initialise timers

  If (.not. TimersInitialised) Then
    Call TimerCreate(ChemistryLoopTimer,"WorkerThread_ChemistryLoop", TimerOpts)
    TimersInitialised = .true.
  End If

End Subroutine ChemistryModuleTimerInitialise

!-------------------------------------------------------------------------------------------------------------

Subroutine ChemistryModuleTimerSummary()

! Summarise timer information in Chemistry module

  Implicit None

  If (TimersInitialised) Then
    Call TimerWriteSummary(ChemistryLoopTimer)
  End If

End Subroutine ChemistryModuleTimerSummary

!-------------------------------------------------------------------------------------------------------------

Function PreInitChemOpts() Result(ChemOpts)
! Pre-initialises the chemistry options.

  Implicit None
  ! Function result:
  Type(ChemOpts_) :: ChemOpts ! Pre-initialised chemistry options.

  ChemOpts%Initialised = .false.

End Function PreInitChemOpts

!-------------------------------------------------------------------------------------------------------------

Subroutine InitChemOpts(ChemFolder, ChemOpts)
! Initialises the chemistry options: specifies the folder storing the STOCHEM background fields.

  Implicit None
  ! Argument list:
  Character(*),    Intent(In)    :: ChemFolder
  Type(ChemOpts_), Intent(InOut) :: ChemOpts
  ! ChemFolder       :: Full path to the STOCHEM folder or path relative to the current
  !                     folder. The later should start with '.'. The path can optionally
  !                     end in '\' or '/' and there is no need to distinguish between '\'
  !                     and '/' for different platforms. If the folder is not specified
  !                     explicitly (i.e. a blank is given in the input file), the default
  !                     folder '../Resources/Stochem/' is used.
  ! ChemOpts         :: Initialised chemistry options.

  If (ChemOpts%Initialised) Then
    Call Message(                                                        &
           'FATAL ERROR: The chemistry options have been given more ' // &
           'than once in the input file(s)',                             &
           3                                                             &
         )
  End If

  Call TokenLengthTest(ChemFolder, MaxFileNameLength, .false., 'Chemistry Options', ' ', 'ChemFolder')

  If (ChemFolder == ' ') Then  ! $$ Retained but not currently used.
    ChemOpts%ChemFolder = '.\'
  Else If (                                                           &
    ChemFolder(Len_Trim(ChemFolder):Len_Trim(ChemFolder)) == '\' .or. &
    ChemFolder(Len_Trim(ChemFolder):Len_Trim(ChemFolder)) == '/'      &
  ) Then
    ChemOpts%ChemFolder = Trim(ChemFolder)
  Else
    Call TokenLengthTest(                                                                     &
           ChemFolder, MaxFileNameLength - 1, .false., 'Chemistry Options', ' ', 'ChemFolder' &
         )
    ChemOpts%ChemFolder = Trim(ChemFolder) // '\'
  End If

  ChemOpts%Initialised = .true.

End Subroutine InitChemOpts

!-------------------------------------------------------------------------------------------------------------

Function InitChemistryDefn(Name) Result(ChemistryDefn)

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name        ! Name of chemistry definition.
  ! Function result:
  Type(ChemistryDefn_) :: ChemistryDefn ! Initialised chemistry definition.
  ! Locals:
  Integer :: i ! Loop index.

  If (Len_Trim(Name) == 0) Then
    Call Message('ERROR in InitChemistryDefn: Name is blank', 3)
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message(                                             &
           'ERROR in InitChemistryDefn: Name is given as ' // &
           Trim(Name)                                      // &
           ' and is too long',                                &
           3                                                  &
         )
  End If

  ChemistryDefn%Name        = Name

  ChemistryDefn%Initialised = .true.

End Function InitChemistryDefn

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpChemistryDefns(Coords, Grids, Specieses, ChemOpts, ChemistryDefn)
! Sets up a definition of a chemistry scheme.

  Implicit None
  ! Argument list:
  Type(Coords_),         Intent(In)    :: Coords
  Type(Grids_),          Intent(In)    :: Grids
  Type(Specieses_),      Intent(In)    :: Specieses
  Type(ChemOpts_),       Intent(In)    :: ChemOpts
  Type(ChemistryDefn_),  Intent(Out)   :: ChemistryDefn
  ! Grids         :: Collection of grids.
  ! Specieses     :: Collection of species.
  ! ChemOpts      :: Chemistry options.
  ! ChemistryDefn :: Initialised definition of a chemistry scheme.
  ! Locals:
  Character(MaxCharLength) :: ChemSpeciesNames(MaxSpecieses)
  Logical                  :: OutputOnly(MaxSpecieses)
  Integer                  :: i, j
  Integer                  :: iParticle
  Integer                  :: iField

  ! Initialise the chemistry definition. !$$ note hard wired grid names etc.
  ChemistryDefn = InitChemistryDefn(                                         &
                    Name        = 'NAME Chemistry Scheme'                    &
                  )

  ! Note m agl always available cos added by flows.
  ChemistryDefn%iZCoordMagl = FindZCoordIndex('m agl', Coords)
  
  ChemistryDefn%nParticleSpecies = Specieses%nParticleSpecieses

  ChemistryDefn%nFieldSpecies = Specieses%nFields
  Do i = 1, ChemistryDefn%nFieldSpecies
    ChemistryDefn%FieldSpeciesNames(i)      = Specieses%Specieses( Specieses%iField2Species(i) )%Name
    ChemistryDefn%FieldSpeciesCategories(i) = Specieses%Specieses( Specieses%iField2Species(i) )%Category
  End Do

  ! Set up molecular weights and half lifes of species in the chemistry definition.
  Do i = 1, Specieses%nParticleSpecieses
    ChemistryDefn%ParticleSpeciesMolWeights(i) = &
      Dble(Specieses%Specieses( Specieses%iParticle2Species(i) )%MolecularWeight)
    ChemistryDefn%ParticleSpeciesInvHalfLife(i) = &
      Dble(Specieses%Specieses( Specieses%iParticle2Species(i) )%InvHalfLife)
  End Do

  ! Set up molecular weights and half lifes of chemistry fields in the chemistry definition.
  ChemistryDefn%nBackFields = 0
  Do i = 1, Specieses%nFields
    ChemistryDefn%FieldSpeciesMolWeights(i) = &
      Dble(Specieses%Specieses( Specieses%iField2Species(i) )%MolecularWeight)
    ChemistryDefn%FieldSpeciesInvHalfLife(i) = &
      Dble(Specieses%Specieses( Specieses%iField2Species(i) )%InvHalfLife)
    If (ChemistryDefn%FieldSpeciesCategories(i) .CIEq. 'CHEMISTRY-FIELD-INIT') Then
      ChemistryDefn%nBackFields = ChemistryDefn%nBackFields + 1
      ChemistryDefn%BackFieldSpeciesNames(ChemistryDefn%nBackFields)      = &
        Specieses%Specieses( Specieses%iField2Species(i) )%Name
      ChemistryDefn%BackFieldSpeciesMolWeights(ChemistryDefn%nBackFields) = &
        Dble(Specieses%Specieses( Specieses%iField2Species(i) )%MolecularWeight)
    End If
  End Do

  ! Set up indices of the chemistry species in species on particles and in species on static fields.
  Call SpeciesListForChemistryScheme(ChemistryDefn%nChemSpecieses, ChemSpeciesNames, OutputOnly)
  ChemistryDefn%iChem2Particle(:) = 0
  ChemistryDefn%iChem2Field   (:) = 0
  Do i = 1, ChemistryDefn%nChemSpecieses
    Do j = 1, Specieses%nParticleSpecieses
      If (ChemSpeciesNames(i) .CIEq. Specieses%Specieses( Specieses%iParticle2Species(j) )%Name) Then
        ChemistryDefn%iChem2Particle(i) = j
      End If
    End Do
    Do j = 1, Specieses%nFields
      If (ChemSpeciesNames(i) .CIEq. Specieses%Specieses( Specieses%iField2Species(j) )%Name) Then
         ChemistryDefn%iChem2Field(i) = j
      End If
    End Do
    If (ChemistryDefn%iChem2Particle(i) /= 0) Then
      If (OutputOnly(i)) Then
        Call Message(                                                                                     &
               'WARNING: Species '                                                                     // &
               Trim(ChemSpeciesNames(i))                                                               // &
               ', which is output by but not input to the chemistry scheme, is a Lagrangian species. ' // &
               'Because it is only output there is no point in storing it on particles.',                 &
               1                                                                                          &
             ) ! $$ In due course add similar comment for advected Eulerian species.
      End If
    End If
    If (ChemistryDefn%iChem2Particle(i) == 0 .and. ChemistryDefn%iChem2Field(i) == 0) Then
      If (OutputOnly(i)) Then
        Call Message(                                                                                   &
               'Note species '                                                                       // &
               Trim(ChemSpeciesNames(i))                                                             // &
               ', which is output by but not input to the chemistry scheme, is not a Lagrangian or ' // &
               'Eulerian species and so is not available for output.',                                  &
               1                                                                                        &
             ) ! $$ This message is probably not needed (?) but for the moment helps to make structure clearer.
               ! Could be useful somewhere to write out ChemSpeciesNames so user can see what names chemistry
               ! is expecting.
      Else
        Call Message(                                                                                      &
               'WARNING: Species '                                                                      // &
               Trim(ChemSpeciesNames(i))                                                                // &
               ', which is used in the chemistry, is not a Lagrangian or Eulerian species. '            // &
               'This species will be passed to the chemistry with zero concentration and any chemical ' // &
               'product will be discarded.',                                                               &
               1                                                                                           &
             )
      End If
    End If
  End Do

  ! Check that the units for all the species used in the chemistry are in grams or Bq,
  ! otherwise fail and write error message.

  ! First initialise ChangeUnits... arrays to false.

  ChemistryDefn%ChangeUnitsParticle(:) = .false.
  ChemistryDefn%ChangeUnitsField(:)    = .false.

  ! If units in Bq then set ChangeUnits to true.

  Do i = 1, ChemistryDefn%nChemSpecieses
    If (ChemistryDefn%iChem2Particle(i) /= 0) Then
      iParticle = ChemistryDefn%iChem2Particle(i)
      If (Specieses%Specieses(Specieses%iParticle2Species(iParticle))%MaterialUnitName .CIEq. 'g') Then
      Else If (Specieses%Specieses(Specieses%iParticle2Species(iParticle))%MaterialUnitName .CIEq. 'Bq') Then
        ChemistryDefn%ChangeUnitsParticle(i)=.true.
      Else
        Call Message(                                                                           &
               'ERROR in SetUpChemistryDefn: a species material unit used by the chemistry ' // &
               'is not in grams or Bq. Check "Species:" block.',                                &
               3                                                                                &
             )
      End If 
    End If
    If (ChemistryDefn%iChem2Field(i) /= 0) Then
      iField = ChemistryDefn%iChem2Field(i)
      If (Specieses%Specieses(Specieses%iField2Species(iField))%MaterialUnitName .CIEq. 'g') Then
      Else If (Specieses%Specieses(Specieses%iField2Species(iField))%MaterialUnitName .CIEq. 'Bq') Then
        ChemistryDefn%ChangeUnitsField(i)=.true.
      Else
        Call Message(                                                                           &
               'ERROR in SetUpChemistryDefn: a species material unit used by the chemistry ' // &
               'is not in grams or Bq. Check "Species:" block.',                                &
               3                                                                                &
             )
      End If
    End If
  End Do

End Subroutine SetUpChemistryDefns

!-------------------------------------------------------------------------------------------------------------

Subroutine PreInitChemistryState(ChemistryState)
! Pre-initialises the state of a chemistry calculation by setting the pointers to
! null.

  Implicit None
  ! Argument list:
  Type(ChemistryState_), Intent(InOut) :: ChemistryState
  ! ChemistryState :: State of a chemistry calculation.

  ! Set pointers to null.
  ChemistryState%ParticleIndex => Null()
  ChemistryState%GridboxState  => Null()
  ChemistryState%Volume        => Null()
  ChemistryState%Mass          => Null()
  ChemistryState%OldMass       => Null()
  ChemistryState%ChemBackField => Null()

  ! Set status to uninitialised.
  ChemistryState%Initialised = .false.

End Subroutine PreInitChemistryState

!-------------------------------------------------------------------------------------------------------------

Subroutine InitChemistryState(             &
             MaxParticles,                 &
             Coords, Grids, Domains, Time, &
             ChemOpts, ChemistryDefn,      &
             Units, Mets, Flows,           &
             EulerianField,                &
             ChemistryState                &
           )

! Initialises the state of a chemistry calculation (at the start of each case running chemistry), including
! the reading of background O3 and H2O2 fields from STOCHEM, if required and the calculation of gridbox volumes.
! This subroutine requires flow information to be available (e.g. for converting volume mixing ratios to
! mass concentrations), and is therefore called immediately after the first met/flow update cycle.

  Implicit None
  ! Argument list:
  Integer,               Intent(In)    :: MaxParticles
  Type(Coords_),         Intent(In)    :: Coords
  Type(Grids_),          Intent(In)    :: Grids
  Type(Domains_),        Intent(In)    :: Domains
  Type(Time_),           Intent(In)    :: Time
  Type(ChemOpts_),       Intent(In)    :: ChemOpts
  Type(ChemistryDefn_),  Intent(In)    :: ChemistryDefn
  Type(Units_),          Intent(InOut) :: Units
  Type(Mets_),           Intent(InOut) :: Mets
  Type(Flows_),          Intent(InOut) :: Flows
  Type(EulerianField_),  Intent(In)    :: EulerianField
  Type(ChemistryState_), Intent(InOut) :: ChemistryState
  ! Coords         :: Collection of coord systems.
  ! Grids          :: Collection of grids.
  ! Domains        :: Collection of domains.
  ! Time           :: Current time.
  ! ChemOpts       :: Chemistry options.
  ! ChemistryDefn  :: Information defining a chemistry scheme.
  ! Units          :: Collection of information on input/output unit numbers.
  ! Mets           :: Collection of met module instance states.
  ! Flows          :: Collection of flows.
  ! ChemistryState :: Initialised state of a chemistry calculation.
  ! Locals:
  Type(HGrid_)    :: HGrid            !} Horizontal and vertical grids and interface
  Type(ZGrid_)    :: ZGrid            !} levels used by the chemistry scheme.
  Type(ZGrid_)    :: ZLevels          !}
  Integer         :: nX               !] Size of horizontal grid.
  Integer         :: nY               !]
  Integer         :: nZ               !  Size of vertical grid.
  Integer         :: nParticleSpecies !  Number of species held on particles.
  Integer         :: nFieldSpecies    !  Number of species held on chemistry grid.
  Integer         :: nBackFields      !  Number of species held on chemistry grid.
  Integer         :: Err              !  Error code from allocate statements.
  Integer         :: iX               !} Indices of a chemistry gridbox.
  Integer         :: iY               !}
  Integer         :: iZ               !}
  Integer         :: i                !
  Integer         :: j                !
  Real(Std)       :: HMax             ! Max metric coefficient.
  Real(Std)       :: H1               !] Metric coefficients (evaluated at each
  Real(Std)       :: H2               !] gridbox location).
  Real(Std)       :: DeltaX           !} Increments in X, Y, Z across each gridbox.
  Real(Std)       :: DeltaY           !}
  Real(Std)       :: DeltaZ           !}
  Integer         :: iHGridSTOCHEM    !] Indices of the horizontal and vertical grids
  Integer         :: iZGridSTOCHEM    !] and interface levels used by STOCHEM.
  Integer         :: iZLevelsSTOCHEM  !]
  Integer         :: nXSTOCHEM        !} Size of STOCHEM horizontal grid.
  Integer         :: nYSTOCHEM        !}
  Integer         :: nZSTOCHEM        !  Size of STOCHEM vertical grid.
  Type(Position_) :: Position         !  Position of gridbox centre.
  Integer         :: iXSTOCHEM        !} Indices of the coincident STOCHEM gridbox.
  Integer         :: iYSTOCHEM        !}
  Integer         :: iZSTOCHEM        !}
  Real(Std)       :: X(3)             ! Location of grid point.
  Type(FlowMemory_) :: FlowMemory ! Flow memory.

  ! $$ note the allocate statements will leak memory in multicase simulations

  ! Allocate the particle index array.
  Allocate(ChemistryState%ParticleIndex(MaxParticles), Stat = Err)
  If (Err /= 0) Then
    Call Message('ERROR in InitChemistryState: problem allocating particle index array', 3)
  End If

  ! Identify horizontal and vertical grids of the chemistry scheme.
  HGrid = Grids%HGrids(EulerianField%iHGrid)
  ZGrid = Grids%ZGrids(EulerianField%iZGrid)

  ! Identify vertical interface levels of the chemistry scheme.
  ZLevels = Grids%ZGrids(EulerianField%iZGridBoundary)

  ! Set up grid sizes for the chemistry grid.
  nX = HGrid%nX
  nY = HGrid%nY
  nZ = ZGrid%nZ
  
  ! Set up number of species and chemistry fields to be stored.
  nParticleSpecies = ChemistryDefn%nParticleSpecies
  nFieldSpecies    = ChemistryDefn%nFieldSpecies
  nBackFields      = ChemistryDefn%nBackFields

  ! Allocate the gridbox-state array.
  Allocate(ChemistryState%GridboxState(nX, nY, nZ), Stat = Err)
  If (Err /= 0) Then
    Call Message('ERROR in InitChemistryState: problem allocating state array', 3)
  End If

  ! Initialise the state of each individual gridbox: n, nOld, k set to -999; DeltaO3 set to 0.0.
  ChemistryState%GridboxState(:,:,:)%n       = -999
  ChemistryState%GridboxState(:,:,:)%nOld    = -999
  ChemistryState%GridboxState(:,:,:)%k       = -999
  ChemistryState%GridboxState(:,:,:)%DeltaO3 = 0.0

  ! Allocate the gridbox-volume array and calculate the volume of every gridbox.
  ! Note that the approach used here supports a variable grid spacing (since each
  ! gridbox is treated individually) - $$ currently not set up in the horizontal
  ! (DeltaX/Y set to fixed grid spacing at present).
  ! $$ Also doesn't aupport non-z-based vertical coordinates where the volume may change with time.
  Allocate(ChemistryState%Volume(nX, nY, nZ), Stat = Err)
  If (Err /= 0) Then
    Call Message('ERROR in InitChemistryState: problem allocating volume array',3)
  End If
  Do iX = 1, nX
  Do iY = 1, nY
  Do iZ = 1, nZ
    Call MetricCoeffs(Coords%HCoords(HGrid%iHCoord),  &
                        (/                            &
                          HGrid%X(iX),                & !$$ Std or Pos?
                          HGrid%Y(iY)                 &
                        /),                           &
                      HMax, H1, H2)
    DeltaX = HGrid%dX
    DeltaY = HGrid%dY

    If (ZLevels%iZCoord == ChemistryDefn%iZCoordMagl) Then
!    DeltaZ = ZLevels%Z(iZ+1) - ZLevels%Z(iZ)
      DeltaZ = ZLevels%AvZ(iZ)
    Else
      Call ResetFlowMemory(Flows, FlowMemory)
      Position = X2Position(                                                              &
                   Coords,                                                                &
                   (/ HGrid%X(iX), HGrid%Y(iY), ZLevels%Z(iZ) + 0.5 * ZLevels%AvZ(iZ) /), &
                   HGrid%iHCoord,                                                         &
                   ZLevels%iZCoord                                                        &
                 )
      Call ConvertToZ(                                                &
             Coords, Grids, Domains,                                  &
             ChemistryDefn%iZCoordMagl,                               &
             Time2ShortTime(Time), .true., ZeroShortTime(), Position, &
             Units, Mets, Flows,                                      &
             FlowMemory,                                              &
             Err                                                      &
           )
      If (Err /= 0) Then
        Call Message('ERROR in InitChemistryState: height conversion error', 3)
      End If
      X(:) = Position2X(Coords, Position, HGrid%iHCoord, ChemistryDefn%iZCoordMagl)
      DeltaZ = X(3)

      Call ResetFlowMemory(Flows, FlowMemory)
      Position = X2Position(                                                              &
                   Coords,                                                                &
                   (/ HGrid%X(iX), HGrid%Y(iY), ZLevels%Z(iZ) - 0.5 * ZLevels%AvZ(iZ) /), &
                   HGrid%iHCoord,                                                         &
                   ZLevels%iZCoord                                                        &
                 )
      Call ConvertToZ(                                                &
             Coords, Grids, Domains,                                  &
             ChemistryDefn%iZCoordMagl,                               &
             Time2ShortTime(Time), .true., ZeroShortTime(), Position, &
             Units, Mets, Flows,                                      &
             FlowMemory,                                              &
             Err                                                      &
           )
      If (Err /= 0) Then
        Call Message('ERROR in InitChemistryState: height conversion error', 3)
      End If
      X(:) = Position2X(Coords, Position, HGrid%iHCoord, ChemistryDefn%iZCoordMagl)
      DeltaZ = DeltaZ - X(3)
    End If

    ChemistryState%Volume(iX, iY, iZ) = Dble((H1 * DeltaX) * (H2 * DeltaY)  &
                                         * DeltaZ * 1.0E6)
  End Do
  End Do
  End Do

  ! Allocate the gridbox-mass arrays and initialise their values to zero.
  Allocate(ChemistryState%Mass(nParticleSpecies, nZ, nY, nX), Stat = Err)
  If (Err /= 0) Then
    Call Message('ERROR in InitChemistryState: problem allocating mass array', 3)
  End If
  ChemistryState%Mass(:,:,:,:) = 0.0
  Allocate(ChemistryState%OldMass(nParticleSpecies, nZ, nY, nX), Stat = Err)
  If (Err /= 0) Then
    Call Message('ERROR in InitChemistryState: problem allocating mass array', 3)
  End If
  ChemistryState%OldMass(:,:,:,:) = 0.0

  ! If using background fields, then allocate arrays for storing background fields.
  If (nBackFields >= 1) Then
!    Allocate(ChemistryState%ChemBackField(nBackFields, nZ, nY, nX), Stat = Err)
!    Allocate(ChemistryState%ChemBackField(nX, nY, nZ + 1, nBackFields), Stat = Err)
    Allocate(ChemistryState%ChemBackField(nX, nY, nZ, nBackFields), Stat = Err)
    If (Err /= 0) Then
      Call Message('ERROR in InitChemistryState: problem allocating background arrays', 3)
    End If
  End If
  ! $$ Should think about treating other fields similarly with "if >= 1" test.
  ! In principle could have nFieldSpecies = 0 if all species on particles, or could have no particles.
  ! Treatment of arrays in restart read/write would then also need to be similar.

  ! Check - only O3 and H2O2 currently possible on background fields $$.
  Do i = 1, ChemistryDefn%nFieldSpecies ! $$ could simplify by using nBackFields and BackFieldSpeciesNames.
    If(ChemistryDefn%FieldSpeciesCategories(i) .CIEq. 'CHEMISTRY-FIELD-INIT') Then
      If (ChemistryDefn%FieldSpeciesNames(i) .CIEq. 'O3') Then 
      Else If (ChemistryDefn%FieldSpeciesNames(i) .CIEq. 'H2O2') Then 
      Else
        Call Message('ERROR in InitChemistryState: currently only possible to initialise O3 or H2O2', 3)
      End If
    End If
  End Do

  Do i = 1, ChemistryDefn%nFieldSpecies
    If (ChemistryDefn%FieldSpeciesNames(i) .CIEq. 'O3') Then  
      If (.not. (ChemistryDefn%FieldSpeciesCategories(i) .CIEq. 'CHEMISTRY-FIELD-INIT')) Then
        Call Message('WARNING in InitChemistryState: O3 field has not been initialised', 1)
      End If
    End If
    If (ChemistryDefn%FieldSpeciesNames(i) .CIEq. 'H2O2') Then  
      If (.not. (ChemistryDefn%FieldSpeciesCategories(i) .CIEq. 'CHEMISTRY-FIELD-INIT')) Then
        Call Message('WARNING in InitChemistryState: H2O2 field has not been initialised', 1)
      End If
    End If
  End Do

  ChemistryState%InitialiseBackFields = .true.

  ChemistryState%Initialised = .true. ! $$ needed?

End Subroutine InitChemistryState

!-------------------------------------------------------------------------------------------------------------

Subroutine SetBackgroundFieldsFromSTOCHEM( &
             Coords, Grids, Domains,       &
             Time, PrevTime,               &
             ChemOpts, ChemistryDefn,      &
             Units, Mets, Flows,           &
             EulerianField,                &
             ChemistryState                &
           )
! Sets the background O3 and H2O2 fields from STOCHEM data (adapted from subroutine READFFCHEM in NAME V8.17).
!$$ Reading of STOCHEM file is hard-coded here. In particular, the grid sizes, etc. are fixed.

  Implicit None
  ! Argument list:
  Type(Coords_),         Intent(In)    :: Coords
  Type(Grids_),          Intent(In)    :: Grids
  Type(Domains_),        Intent(In)    :: Domains
  Type(Time_),           Intent(In)    :: Time
  Type(Time_),           Intent(In)    :: PrevTime
  Type(ChemOpts_),       Intent(In)    :: ChemOpts
  Type(ChemistryDefn_),  Intent(In)    :: ChemistryDefn
  Type(Units_),          Intent(InOut) :: Units
  Type(Mets_),           Intent(InOut) :: Mets
  Type(Flows_),          Intent(InOut) :: Flows
  Type(EulerianField_),  Intent(In)    :: EulerianField
  Type(ChemistryState_), Intent(InOut) :: ChemistryState
  ! Coords         :: Collection of coord systems.
  ! Grids          :: Collection of grids.
  ! Domains        :: Collection of domains.
  ! Time           :: Current time.
  ! PrevTime       :: Previous time (infinite past or previous synchronisation time).
  ! ChemOpts       :: Chemistry options.
  ! ChemistryDefn  :: Information defining chemistry scheme.
  ! Units          :: Collection of information on input/output unit numbers.
  ! Mets           :: Collection of met module instance states.
  ! Flows          :: Collection of flows.
  ! ChemistryState :: State of chemistry calculation.
  ! Locals:
  Character(MaxCharLength)     :: CharTime        ! Current time as a character string
  Character(MaxCharLength)     :: CharPrevTime    ! Previous time as a character string
  Character(4)                 :: CharYear        ! Current year as a character string
  Character(2)                 :: CharMonth       ! Current month as a character string
  Type(HGrid_)                 :: HGrid           !] Main chemistry grid - horizontal and vertical grids
  Type(ZGrid_)                 :: ZGrid           !]
  Integer                      :: iHGridSTOCHEM   !} STOCHEM grid - grid indices and sizes
  Integer                      :: iZGridSTOCHEM   !}
  Integer                      :: iZLevelsSTOCHEM !}
  Integer                      :: nXSTOCHEM       !}
  Integer                      :: nYSTOCHEM       !}
  Integer                      :: nZSTOCHEM       !}
  Character(MaxFileNameLength) :: Filename        ! Filename (including path) of STOCHEM monthly data files
  Integer                      :: I, J, N, L      ! Loop variables used in reading of STOCHEM files
  Real(8)                      :: TEMP(12)        ! Temporary storage array used in reading data
  Integer                      :: iX, iY, iZ      ! Indices for looping over main chemistry grid
  Type(Position_)              :: Position        ! Position of chemistry grid point
  Integer                      :: iXSTOCHEM       !] Index of coincident STOCHEM gridbox
  Integer                      :: iYSTOCHEM       !]
  Integer                      :: iZSTOCHEM       !]
  Real(8)                      :: ConvFactor      ! Conversion factor (volume mixing ratio -> moles/m3).
                                              ! MET DATA FOR CHEMISTRY GRIDPOINT:
  Real(Std) :: T                              ! Ambient air temperature.
  Real(Std) :: Q                              ! Specific humidity.
  Real(Std) :: P                              ! Ambient pressure.
  Real(Std) :: MolecularAirDensity            ! Air density (in molecules/cm3).
  Real(Std) :: H                              ! Boundary layer height.
  Real(Std) :: Topog                          ! Topographic height above sea level.
  Real(Std) :: TotalOrDynCloudWater           ! Total or dynamic cloud liquid water.
  Real(Std) :: TotalOrDynCloudIce             ! Total or dynamic cloud ice.
  Real(Std) :: CloudFraction                  ! Total (dyn+conv) cloud fraction.
  Real(Std) :: TotalColumnCloudFraction       ! Max cloud fraction in cloud column above gridpoint.
                                                                 ! STOCHEM DATA:
  Real(8)   :: O3LIM  (MaxSTOCHEMnZ, MaxSTOCHEMnY, MaxSTOCHEMnX) ! O3   values
  Real(8)   :: H2O2LIM(MaxSTOCHEMnZ, MaxSTOCHEMnY, MaxSTOCHEMnX) ! H2O2 values
  Logical   :: ReadO3, ReadH2O2


  ! Only update background fields at start of each month.
  CharTime     = FileNameTime(Time)
  CharPrevTime = FileNameTime(PrevTime)
  If (CharTime(5:6) == CharPrevTime(5:6)) Return

  ! Extract year and month from the given time.
  If (.not. IsCalendar()) Then
    Call Message('ERROR in SetBackgroundFieldsFromSTOCHEM: a calendar time frame is needed here', 3)
  End If
  CharYear  = CharTime(1:4)
  CharMonth = CharTime(5:6)

  !
  ! Properties of main chemistry grid.
  !

  ! Horizontal and vertical grids used.
  HGrid = Grids%HGrids(EulerianField%iHGrid)
  ZGrid = Grids%ZGrids(EulerianField%iZGrid)

  !
  ! Properties of STOCHEM grid.
  !

  ! Find the indices of the STOCHEM horizontal and vertical grids. $$ note grid names hardwired at present
  iHGridSTOCHEM = FindHGridIndex('STOCHEM - horizontal grid', Grids)
  iZGridSTOCHEM = FindZGridIndex('STOCHEM - vertical grid',   Grids)

  ! Find the index of the STOCHEM vertical interface levels.
  iZLevelsSTOCHEM = FindZGridIndex('STOCHEM - interface levels', Grids)

  ! Set up grid sizes of the STOCHEM grid.
  nXSTOCHEM = Grids%HGrids(iHGridSTOCHEM)%nX
  nYSTOCHEM = Grids%HGrids(iHGridSTOCHEM)%nY
  nZSTOCHEM = Grids%ZGrids(iZGridSTOCHEM)%nZ

  ! Check that STOCHEM horizontal grid is not too big.
  If (nXSTOCHEM > MaxSTOCHEMnX .or. nYSTOCHEM > MaxSTOCHEMnY) Then
    Call Message('ERROR in SetBackgroundFieldsFromSTOCHEM: STOCHEM horizontal grid is too big', 3)
  End If

  ! Check that there are not too many STOCHEM vertical levels.
  If (nZSTOCHEM > MaxSTOCHEMnZ) Then
    Call Message('ERROR in SetBackgroundFieldsFromSTOCHEM: too many STOCHEM vertical levels', 3)
  End If

  ! Check that STOCHEM horizontal grid agrees with hard coded values below. $$
  If (nXSTOCHEM /= 72 .or. nYSTOCHEM /= 36) Then
    Call Message('ERROR in SetBackgroundFieldsFromSTOCHEM: STOCHEM grid must be 72 by 36', 3)
  End If

  ! Identify fields to read.
  ReadO3   = .false.
  ReadH2O2 = .false.
  Do i = 1, ChemistryDefn%nBackFields 
    If (ChemistryDefn%BackFieldSpeciesNames(i) .CIeq. 'O3') Then 
      ReadO3 = .true. 
    ElseIf (ChemistryDefn%BackFieldSpeciesNames(i) .CIeq. 'H2O2') Then 
      ReadH2O2 = .true. 
    Else
      Call Message(                                                                                      &
             'UNEXPECTED ERROR in InitChemistryState: currently only possible to initialise O3 or H2O2', &
             4                                                                                           &
           )
      ! Error is "unexpected" because this is trapped earlier. But error message left in here as a reminder
      ! of existing restrictions for when we make this routine more flexible.
    Endif
  End Do

  !
  ! Get STOCHEM data.
  !

  ! 1. O3 data
  If (ReadO3) Then

    ! Generate filename for STOCHEM O3 data.
    ! $$ Directory structures on the various platforms need to be clarified.
#   ifdef CompaqPCCompiler
      Filename = Trim(ChemOpts%ChemFolder) // 'O3' // CharMonth // '.GRD'
#   else
      Filename = Trim(ChemOpts%ChemFolder) // 'O3' // CharMonth // '.GRD'
!     Filename = Trim(ChemOpts%ChemFolder) // CharYear // '/' // 'O3' // CharMonth // '.GRD'
#   endif

    ! Open O3 data file. $$ The unit number is hard-wired here at present.
    ! $$ use openfile in due course; trap errors

    OPEN(71, File = Trim(ConvertFileName(Filename)), Status = 'Old', Action = 'Read')

    ! Read O3

    READ(71,*)
    DO L=1,nZSTOCHEM
      READ(71,*)
      READ(71,*)
      DO N=1,6
        READ(71,*)
        READ(71,*)
        DO I=1,36
          READ(71,*) TEMP
          DO J=1,12
            O3LIM(L,I,(N-1)*12+J) = TEMP(J)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CLOSE(71)

  End If

  ! 2. H2O2 data
  If (ReadH2O2) Then

    ! Generate filename for STOCHEM H2O2LIM data.
    ! $$ Directory structures on the various platforms need to be clarified.
#   ifdef CompaqPCCompiler
      Filename = Trim(ChemOpts%ChemFolder) // 'H2O2LIM' // CharMonth // '.GRD'
#   else
      Filename = Trim(ChemOpts%ChemFolder) // 'H2O2LIM' // CharMonth // '.GRD'
!     Filename = Trim(ChemOpts%ChemFolder) // CharYear // '/' // 'H2O2LIM' // CharMonth // '.GRD'
#   endif

    ! Open H2O2LIM data file. $$ The unit number is hard-wired here at present.
    ! $$ use openfile in due course; trap errors

    OPEN(71, File = Trim(ConvertFileName(Filename)), Status = 'Old', Action = 'Read')

    ! Read H2O2LIM

    READ(71,*)
    DO L=1,nZSTOCHEM
      READ(71,*)
      READ(71,*)
      DO N=1,6
        READ(71,*)
        READ(71,*)
        DO I=1,36
          READ(71,*) TEMP
          DO J=1,12
            H2O2LIM(L,I,(N-1)*12+J) = TEMP(J)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CLOSE(71)

  End If

  !
  ! Interpolate STOCHEM data onto main chemistry grid - use nearest STOCHEM point for interpolation.
  !

  ! Loop over main chemistry grid.
  Do iX = 1, HGrid%nX
  Do iY = 1, HGrid%nY
  Do iZ = 1, ZGrid%nZ

    ! Calculate the coincident STOCHEM gridbox for setting background fields.
    Position = X2Position(                    &
                 Coords,                      &
                 (/                           &
                   HGrid%X(iX),               &
                   HGrid%Y(iY),               &
                   ZGrid%Z(iZ)                &
                  /),                         &
                 HGrid%iHCoord, ZGrid%iZCoord &
               )
    Call FindSTOCHEMGridbox(               &
           Coords, Grids, Domains,         &
           iHGridSTOCHEM, iZLevelsSTOCHEM, &
           Time2ShortTime(Time), Position, &
           Units, Mets, Flows,             &
           iXSTOCHEM, iYSTOCHEM, iZSTOCHEM &
         )

    ! Retrieve met data for chemistry gridbox - but only actually need air density here!
    Call GetMetForChemistry(                                       &
           Coords, Grids, Domains, Time2ShortTime(Time), Position, &
           Units, Mets, Flows,                                     &
           T, Q, P, MolecularAirDensity, H, Topog,                 &
           TotalOrDynCloudWater, TotalOrDynCloudIce,               &
           CloudFraction, TotalColumnCloudFraction                 &
         )

    ! STOCHEM concentrations are volume mixing ratios - need to convert to g/m3 by
    !  a) multiplying by conversion factor (volume mixing ratio -> moles/m3), then
    !  b) multiplying by molecular weight of species.
    ConvFactor = 1.0E6 * MolecularAirDensity / Dble(Avogadro)

    Do i = 1, ChemistryDefn%nBackFields 
      If (ChemistryDefn%BackFieldSpeciesNames(i) .CIEq. 'O3') Then 
        ChemistryState%ChemBackField(iX, iY, iZ, i) = O3LIM(iZSTOCHEM, iYSTOCHEM, iXSTOCHEM)   * ConvFactor &
                                                     * ChemistryDefn%BackFieldSpeciesMolWeights(i)
      Else If (ChemistryDefn%BackFieldSpeciesNames(i) .CIEq. 'H2O2') Then 
        ChemistryState%ChemBackfield(iX, iY, iZ, i) = H2O2LIM(iZSTOCHEM, iYSTOCHEM, iXSTOCHEM) * ConvFactor &
                                                     * ChemistryDefn%BackFieldSpeciesMolWeights(i)
      End If
    End Do

  End Do
  End Do
  End Do
  
End Subroutine SetBackgroundFieldsFromSTOCHEM

!-------------------------------------------------------------------------------------------------------------

Subroutine FindSTOCHEMGridbox(                         &
             Coords, Grids, Domains, iHGrid, iZLevels, &
             Time, Position,                           &
             Units, Mets, Flows,                       &
             iXSTOCHEM, iYSTOCHEM, iZSTOCHEM           &
           )
! Calculates the coincident STOCHEM gridbox at a given position.

  Implicit None
  ! Argument list:
  Type(Coords_),   Intent(In)    :: Coords     ! Collection of coord systems.
  Type(Grids_),    Intent(In)    :: Grids      ! Collection of grids.
  Type(Domains_),  Intent(In)    :: Domains    ! Collection of domains.
  Integer,         Intent(In)    :: iHGrid     ! Index of horizontal grid.
  Integer,         Intent(In)    :: iZLevels   ! Index of interface levels.
  Type(ShortTime_),Intent(In)    :: Time       ! Current time.
  Type(Position_), Intent(InOut) :: Position   ! Position (updated in routine).
  Type(Units_),    Intent(InOut) :: Units      ! Collection of information on input/output unit numbers.
  Type(Mets_),     Intent(InOut) :: Mets       ! Collection of met module instance states.
  Type(Flows_),    Intent(InOut) :: Flows      ! Collection of flows.
  Integer,         Intent(Out)   :: iXSTOCHEM  !} Indices of the nearest
  Integer,         Intent(Out)   :: iYSTOCHEM  !} STOCHEM gridpoint.
  Integer,         Intent(Out)   :: iZSTOCHEM  !}
  ! Locals:
  Type(HGrid_)      :: HGrid      !} Horizontal grid and vertical interface levels
  Type(ZGrid_)      :: ZLevels    !} of the STOCHEM fields.
  Integer           :: iHCoord    !} Indices of STOCHEM coord systems.
  Integer           :: iZCoord    !}
  Integer           :: Err        ! Error code from vertical coord conversion.
  Type(FlowMemory_) :: FlowMemory ! Flow memory.
  Real(Std)         :: X(3)       ! Coordinates in STOCHEM coord systems.
  Type(HCoeffs_)    :: HCoeffs    ! Horizontal interpolation coefficients.
  Type(ZCoeffs_)    :: ZCoeffs    ! Vertical interpolation coefficients.
  Type(ShortTime_)  :: TravelTime ! Dummy variable for TravelTime argument in
                                  ! ConvertToZ.

  ! Identify horizontal grid and vertical interface levels of the STOCHEM fields.
  HGrid   = Grids%HGrids(iHGrid)
  ZLevels = Grids%ZGrids(iZLevels)

  ! Set up indices of the STOCHEM coord systems.
  iHCoord = HGrid%iHCoord
  iZCoord = ZLevels%iZCoord

  ! Calculate position in STOCHEM coordinates (if not already available).
  If (.not.Position%XYValid(iHCoord)       .or. &
      .not.Position%ZValid (iZCoord)) Then
    Call ResetFlowMemory(Flows, FlowMemory)
    Call ConvertToH(Coords, iHCoord, Position)
    Call ConvertToZ(                           &
           Coords, Grids, Domains,             &
           iZCoord,                            &
           Time, .true., TravelTime, Position, &
           Units, Mets, Flows,                 &
           FlowMemory,                         &
           Err                                 &
         )
    If (Err /= 0) Then
      Call Message('ERROR in FindSTOCHEMGridbox: height conversion error', 3)
    End If
  End If

  ! Retrieve coordinates in STOCHEM coord systems.
  X = Position2X(Coords, Position, iHCoord, iZCoord)

  ! Calculate interpolation coeffs.
  Call GetHCoeffs(X(1), X(2), HGrid, HCoeffs)
  Call GetZCoeffs(X(3), ZLevels, ZCoeffs)

  iXSTOCHEM = HCoeffs%iXNear
  iYSTOCHEM = HCoeffs%iYNear
  iZSTOCHEM = Min(ZCoeffs%iZ1, ZCoeffs%iZ2)

End Subroutine FindSTOCHEMGridbox

!-------------------------------------------------------------------------------------------------------------

Subroutine WriteRestartFileChemistry(Unit, ChemistryState)
! Write information on the chemistry state to a restart file.

! $$ would be better to code the read and write routines so as to not make any assumptions about how the
! chemistry scheme works. A start's been made with ChemBackField (following approach in output.f90).

  Implicit None
  ! Argument list:
  Integer,               Intent(In) :: Unit           ! Input/output unit number.
  Type(ChemistryState_), Intent(In) :: ChemistryState ! State of chemistry calculation.

  ! Write control information and array bounds for the chemistry arrays.
  Write (Unit) ChemistryState%InitialiseBackFields
  Write (Unit) Size(ChemistryState%ParticleIndex)
  Write (Unit) Size(ChemistryState%Mass,      Dim = 4)
  Write (Unit) Size(ChemistryState%Mass,      Dim = 3)
  Write (Unit) Size(ChemistryState%Mass,      Dim = 2)
  Write (Unit) Size(ChemistryState%Mass,      Dim = 1)

  ! Gridbox-state array.
  Write (Unit) ChemistryState%GridboxState(:, :, :)

  ! Gridbox-volume array.
  Write (Unit) ChemistryState%Volume(:, :, :)

  ! Gridbox-mass arrays.
  Write (Unit) ChemistryState%Mass(:, :, :, :)
  Write (Unit) ChemistryState%OldMass(:, :, :, :)

  ! Background arrays.
  Write (Unit) Associated(ChemistryState%ChemBackField)
  If (Associated(ChemistryState%ChemBackField)) Then
    Write (Unit) Shape(ChemistryState%ChemBackField)
    Write (Unit) ChemistryState%ChemBackField(:, :, :, :)
  End If

End Subroutine WriteRestartFileChemistry

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadRestartFileChemistry(Unit, ChemistryState)
! Read information on the chemistry state from a restart file.

  Implicit None
  ! Argument list:
  Integer,               Intent(In)  :: Unit           ! Input/output unit number.
  Type(ChemistryState_), Intent(Out) :: ChemistryState ! State of chemistry calculation.
  ! Locals:
  Integer :: MaxParticles     ! Size of particle index array.
  Integer :: nX               !} Size of chemistry grid.
  Integer :: nY               !}
  Integer :: nZ               !}
  Integer :: nParticleSpecies ! Number of species held on particles.
  Integer :: Error            ! Error code.
  Logical :: IsAssociated     ! Indicates the allocatable array should be allocated.
  Integer :: ArrayShape(4)    ! Shape of array to be allocated.
  Integer :: ErrorCode        ! Error code.

  Read (Unit) ChemistryState%InitialiseBackFields
  ! Read array bounds for the chemistry arrays.
  Read (Unit) MaxParticles
  Read (Unit) nX
  Read (Unit) nY
  Read (Unit) nZ
  Read (Unit) nParticleSpecies

  ! Particle index array (only needs to be allocated here).
  Allocate(ChemistryState%ParticleIndex(MaxParticles), Stat = Error)
  If (Error /= 0) Then
    Call Message('FATAL ERROR: Unable to allocate particle index array for chemistry', 3)
  End If

  ! Gridbox-state array.
  Allocate(ChemistryState%GridboxState(nX, nY, nZ), Stat = Error)
  If (Error /= 0) Then
    Call Message('FATAL ERROR: Unable to allocate gridbox-state array for chemistry', 3)
  End If
  Read (Unit) ChemistryState%GridboxState(:, :, :)

  ! Gridbox-volume array.
  Allocate(ChemistryState%Volume(nX, nY, nZ), Stat = Error)
  If (Error /= 0) Then
    Call Message('FATAL ERROR: Unable to allocate gridbox-volume array for chemistry', 3)
  End If
  Read (Unit) ChemistryState%Volume(:, :, :)

  ! Gridbox-mass arrays.
  Allocate(ChemistryState%Mass(nParticleSpecies, nZ, nY, nX), Stat = Error)
  If (Error /= 0) Then
    Call Message('FATAL ERROR: Unable to allocate gridbox-mass array for chemistry', 3)
  End If
  Read (Unit) ChemistryState%Mass(:, :, :, :)
  Allocate(ChemistryState%OldMass(nParticleSpecies, nZ, nY, nX), Stat = Error)
  If (Error /= 0) Then
    Call Message('FATAL ERROR: Unable to allocate gridbox-mass array for chemistry', 3)
  End If
  Read (Unit) ChemistryState%OldMass(:, :, :, :)

  ! Chemistry background fields array.
  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1:4)
    Allocate(                       &
      ChemistryState%ChemBackField( &
        ArrayShape(1),              &
        ArrayShape(2),              &
        ArrayShape(3),              &
        ArrayShape(4)               &
      ),                            &
      Stat = ErrorCode              &
    )
    If (ErrorCode /= 0) Call Message('FATAL ERROR: Unable to allocate arrays for chemistry fields', 3)
    Read (Unit, End = 1, Err = 2) ChemistryState%ChemBackField(:, :, :, :)
  End If

  ChemistryState%Initialised = .true.

  Return

1 Continue
  Call Message('FATAL ERROR: The restart file is shorter than expected', 3)

2 Continue
  Call Message('FATAL ERROR: An error occurred when reading the restart file', 3)

End Subroutine ReadRestartFileChemistry

!-------------------------------------------------------------------------------------------------------------

Subroutine InitChemistryGridboxes(                     &
             Coords, Grids, Domains, Time, nParticles, &
             Particles, Masses,                        &
             ChemistryDefn,                            &
             Units, Mets, Flows,                       &
             EulerianField,                            &
             ChemistryState                            &
           )
! Initialises the chemistry gridbox information (at the beginning of each call
! to the chemistry scheme).

  Implicit None
  ! Argument list:
  Type(Coords_),         Intent(In)           :: Coords
  Type(Grids_),          Intent(In)           :: Grids
  Type(Domains_),        Intent(In)           :: Domains
  Type(Time_),           Intent(In)           :: Time
  Integer,               Intent(In)           :: nParticles
  Type(Particle_),       Intent(In),   Target :: Particles(:)
  Real(Std),             Intent(In)           :: Masses(:, :)
  Type(ChemistryDefn_),  Intent(In)           :: ChemistryDefn
  Type(Units_),          Intent(InOut)        :: Units
  Type(Mets_),           Intent(InOut)        :: Mets
  Type(Flows_),          Intent(InOut)        :: Flows
  Type(EulerianField_),  Intent(In)           :: EulerianField
  Type(ChemistryState_), Intent(InOut)        :: ChemistryState
  ! Coords         :: Collection of coord systems.
  ! Grids          :: Collection of grids.
  ! Domains        :: Collection of domains.
  ! Time           :: Current time.
  ! nParticles     :: Number of active particles currently defined in DispState.
  ! Particles      :: Collection of particles currently defined in DispState.
  ! Masses         :: Masses of the various species carried by the particles.
  ! ChemistryDefn  :: A chemistry definition.
  ! Units          :: Collection of information on input/output unit numbers.
  ! Mets           :: Collection of met module instance states.
  ! Flows          :: Collection of flows.
  ! ChemistryState :: State of a chemistry calculation.
  ! Locals:
  Type(HGrid_)             :: HGrid     !} Horizontal and vertical grids and interface
  Type(ZGrid_)             :: ZGrid     !} levels used by the chemistry grid.
  Type(ZGrid_)             :: ZLevels   !}
  Real(Std)                :: Xl        !] Values of X and Y at the edges of the
  Real(Std)                :: Xu        !] chemistry grid.
  Real(Std)                :: Yl        !]
  Real(Std)                :: Yu        !]
  Real(Std)                :: XMin      !} Minimum and maximum values of X, Y and Z
  Real(Std)                :: XMax      !} for the chemistry grid.
  Real(Std)                :: YMin      !}
  Real(Std)                :: YMax      !}
  Real(Std)                :: ZMin      !}
  Real(Std)                :: ZMax      !}
  Integer                  :: iX        !] Indices of chemistry gridbox.
  Integer                  :: iY        !]
  Integer                  :: iZ        !]
  Integer                  :: iP        ! Index of a particle.
  Integer                  :: j         ! Species loop index.
  Integer                  :: r         ! Temporary variable.
  Integer                  :: n         ! Number of particles in a gridbox.
  Integer                  :: k         ! Position of particle in ParticleIndex array.
  Integer                  :: Err       ! Error code from memory allocation statement
                                        ! or vertical coord conversion routine.
  Type(Particle_), Pointer :: Particle  ! A particle.
  Type(Position_)          :: Position  ! Position information for a particle.
  Type(FlowMemory_)        :: FlowMemory !
  Real(Std)                :: X(3)      ! Particle location in chemistry grid coords.
  Type(HCoeffs_)           :: HCoeffs   ! Horizontal interpolation coefficients.
  Type(ZCoeffs_)           :: ZCoeffs   ! Vertical interpolation coefficients.
  Integer                  :: GridboxID(nParticles) ! Identifies gridbox occupied
                                                    ! by each particle. A value of -1
                                                    ! indicates that the particle is
                                                    ! outside the chemistry grid.
  !$$ Note that the use of a single integer to identify a gridbox is more efficient
  !$$ on memory but does require extra calculation steps - which is more significant?
  Type(ShortTime_)  :: DummyTravelTime ! Dummy variable for TravelTime argument in
                                       ! ConvertToZ.
  Integer :: kMax, kMin, kk
  Real(Std) :: ZBoundary

  ! Identify horizontal and vertical grids of the chemistry scheme.
  HGrid = Grids%HGrids(EulerianField%iHGrid)
  ZGrid = Grids%ZGrids(EulerianField%iZGrid)

  ! Identify vertical interface levels of the chemistry scheme.
  ZLevels = Grids%ZGrids(EulerianField%iZGridBoundary)

  ! Calculate horizontal bounds of the chemistry grid.
  ! $$ This assumes a uniform horizontal grid - always the case for chemistry?
  Xl = HGrid%X(1)        - 0.5 * HGrid%dX
  Xu = HGrid%X(HGrid%nX) + 0.5 * HGrid%dX
  Yl = HGrid%Y(1)        - 0.5 * HGrid%dY
  Yu = HGrid%Y(HGrid%nY) + 0.5 * HGrid%dY
  XMin = Min(Xl, Xu)
  XMax = Max(Xl, Xu)
  YMin = Min(Yl, Yu)
  YMax = Max(Yl, Yu)

  ! Vertical bounds of the chemistry grid.
  ZMin = ZLevels%Z(1) - 0.5 * ZLevels%AvZ(1)
  ZMax = ZLevels%Z(ZLevels%nZ) + 0.5 * ZLevels%AvZ(ZLevels%nZ)

  ! Initialise gridbox mass fields to zero.
  ChemistryState%Mass(:, :, :, :) = 0.0

  ! Store previous particle count in the nOld component. Then re-initialise particle count to zero
  ! in each gridbox. Note that the array index k is also initialised to zero here.
  Do iX = 1, HGrid%nX
  Do iY = 1, HGrid%nY
  Do iZ = 1, ZGrid%nZ
    ChemistryState%GridboxState(iX, iY, iZ)%nOld = ChemistryState%GridboxState(iX, iY, iZ)%n
    ChemistryState%GridboxState(iX, iY, iZ)%n    = 0
    ChemistryState%GridboxState(iX, iY, iZ)%k    = 0
  End Do
  End Do
  End Do

  ! First loop over all particles - to determine their gridbox positions and
  ! count the total number in each gridbox.
  Do iP = 1, nParticles
    Particle => Particles(iP)
    If (.not.ParticleReallyActive(Particle)) Cycle
    ! $$ Note this test used to be If (.not.ParticleActive(Particle)). ParticleReallyActive is necessary to
    ! exclude particles outside the outer computational domain for which the conversions will probably fail.
    ! It also excludes particles which have failed an attempt to get flow info and have been marked.
    ! The exact criteria should however be reviewed in due course. If test further down in this routine
    ! needs to match this one. Is there a danger that due to precision issues, a particle could be
    ! really active but fail to convert? Could mark particle in this case and note particles being lost
    ! due to absence of a valid flow module.
    X = Particle%X
    If (Particle%iHCoord /= HGrid%iHCoord       .or. &
        Particle%iZCoord /= ZGrid%iZCoord) Then
      Call ResetFlowMemory(Flows, FlowMemory)
      Position = X2Position(Coords, X, Particle%iHCoord, Particle%iZCoord)
      Call ConvertToH(Coords, HGrid%iHCoord, Position)
      Call ConvertToZ(                                                &
             Coords, Grids, Domains,                                  &
             ZGrid%iZCoord,                                           &
             Time2ShortTime(Time), .true., DummyTravelTime, Position, &
             Units, Mets, Flows,                                      &
             FlowMemory,                                              &
             Err                                                      &
           ) !$$ Use MetTime or Time ??
      If (Err /= 0) Then
        GridboxID(iP) = -1
        Cycle
        ! Call Message('ERROR in InitChemistryGridboxes: height conversion error', 3)
      End If
      X = Position2X(Coords, Position, HGrid%iHCoord, ZGrid%iZCoord)
    End If
    If (X(1) < XMin .or. X(1) > XMax       .or. &
        X(2) < YMin .or. X(2) > YMax       .or. &
        X(3) < ZMin .or. X(3) > ZMax) Then
      GridboxID(iP) = -1
    Else
      Call GetHCoeffs(X(1), X(2), HGrid, HCoeffs)
  !    Call GetZCoeffs(X(3), ZLevels, ZCoeffs)
      iX = HCoeffs%iXNear
      iY = HCoeffs%iYNear
  !    iZ = Min(ZCoeffs%iZ1, ZCoeffs%iZ2)

      ! This next bit is like GetZCoeffs, but using boundary values to determine nearest point.
      ! See also similar bit in CalcEulerianResults. $$
      ! Find levels surrounding point.
      kMax = ZLevels%nZ + 2 ! Note there are ZLevels%nZ + 1 grid cell boundaries.
      kMin = 0
      Do
        kk = kMax - kMin
        If (kk <= 1) Exit
        kk = kMin + kk/2
        If (kk == ZLevels%nZ + 1) Then
          ZBoundary = ZLevels%Z(kk-1) + 0.5 * ZLevels%AvZ(kk-1) 
        Else
          ZBoundary = ZLevels%Z(kk) - 0.5 * ZLevels%AvZ(kk) 
        End If
        If ((X(3) - ZBoundary) * ZLevels%IncreasingValue > 0.0) Then
          kMin = kk
        Else
          kMax = kk
        End If
      End Do
      ! Set iZ.
      iZ = Max(kMin, 1)

      ChemistryState%GridboxState(iX, iY, iZ)%n =     &
        ChemistryState%GridboxState(iX, iY, iZ)%n + 1
      GridboxID(iP) = (iX - 1) + (iY - 1)*HGrid%nX + (iZ - 1)*HGrid%nX*HGrid%nY
    End If
  End Do

  ! Next assign the starting indices k for those gridboxes containing particles.
  k = 1
  Do iX = 1, HGrid%nX
  Do iY = 1, HGrid%nY
  Do iZ = 1, ZGrid%nZ
    n = ChemistryState%GridboxState(iX, iY, iZ)%n
    If (n == 0) Cycle
    ChemistryState%GridboxState(iX, iY, iZ)%k = k
    k = k + n
  End Do
  End Do
  End Do

  ! Second loop over all particles - to generate the gridbox mass fields and
  ! insert each particle into the appropriate position of the indexing array.
  Do iP = 1, nParticles
    Particle => Particles(iP)
    If (.not.ParticleReallyActive(Particle)) Cycle
    If (GridboxID(iP) == -1) Cycle
    r = GridboxID(iP)
    iZ = Int(Float(r)/Float(HGrid%nX*HGrid%nY))
    r = r - iZ*HGrid%nX*HGrid%nY
    iY = Int(Float(r)/Float(HGrid%nX))
    r = r - iY*HGrid%nX
    iX = r+1; iY = iY+1; iZ = iZ+1
    Do j = 1, ChemistryDefn%nParticleSpecies
      ChemistryState%Mass(j, iZ, iY, iX) = ChemistryState%Mass(j, iZ, iY, iX) + Masses(j, iP)
    End Do
    k = ChemistryState%GridboxState(iX, iY, iZ)%k
    ChemistryState%ParticleIndex(k) = iP
    ChemistryState%GridboxState(iX, iY, iZ)%k = k+1
  End Do

  ! Create a copy of the initial gridbox mass fields
  ! - to allow us to calculate mass changes during chemistry.
  ChemistryState%OldMass = ChemistryState%Mass

End Subroutine InitChemistryGridboxes

!-------------------------------------------------------------------------------------------------------------

Subroutine GetMetForChemistry(                         &
             Coords, Grids, Domains, Time, Position,   &
             Units, Mets, Flows,                       &
             T, Q, P, MolecularAirDensity, H, Topog,   &
             TotalOrDynCloudWater, TotalOrDynCloudIce, &
             CloudFraction, TotalColumnCloudFraction   &
           )
! Retrieves the meteorological data needed in the chemistry scheme. Met data (flow
! and cloud information) is retrieved at the centre of each chemistry gridbox.
! $$ Returns the same information as the subroutine METCHEM in NAME V8.08 but now
! $$ uses the NAME III met extraction routines (different interpolations???).

  Implicit None
  ! Argument list:
  Type(Coords_),   Intent(In)    :: Coords              ! Collection of coord systems
  Type(Grids_),    Intent(In)    :: Grids               ! Collection of grids
  Type(Domains_),  Intent(In)    :: Domains             ! Collection of domains
  Type(ShortTime_),Intent(In)    :: Time                ! Current time
  Type(Position_), Intent(InOut) :: Position            ! Position of chemistry gridbox
  Type(Units_),    Intent(InOut) :: Units               ! Collection of information on input/output
                                                        ! unit numbers.
  Type(Mets_),     Intent(InOut) :: Mets                ! Collection of met module instance states.
  Type(Flows_),    Intent(InOut) :: Flows               ! Collection of flows
  Real(Std),       Intent(Out)   :: T                   ! Temperature
  Real(Std),       Intent(Out)   :: Q                   ! Specific humidity
  Real(Std),       Intent(Out)   :: P                   ! Ambient pressure
  Real(Std),       Intent(Out)   :: MolecularAirDensity ! Air density (in molecules/cm3)
  Real(Std),       Intent(Out)   :: H                   ! Boundary layer height
  Real(Std),       Intent(Out)   :: Topog               ! Topographic height above sea level.
  Real(Std),       Intent(Out)   :: TotalOrDynCloudWater       ! Total or dynamic cloud liquid water
  Real(Std),       Intent(Out)   :: TotalOrDynCloudIce         ! Total or dynamic cloud ice
  Real(Std),       Intent(Out)   :: CloudFraction              ! Total (dyn+conv) cloudfraction
  Real(Std),       Intent(Out)   :: TotalColumnCloudFraction   ! Total cloud fraction of
                                                               ! vertical column above
                                                               ! the chemistry gridbox
  ! Locals:
  Integer           :: Err        ! Error code from met retrievals.
  Type(FlowMemory_) :: FlowMemory ! Flow memory.
  Type(Flow_)       :: Flow       ! Flow information.
  Type(Cloud_)      :: Cloud      ! Cloud information.
  Type(Rain_)       :: Rain       ! Rain information (dummy variable for call to GetAttrib).
  Type(Surface_)    :: Surface    ! Surface information (dummy variable for call to GetAttrib).
  Type(Soil_)       :: Soil       ! Soil information (dummy variable for call to GetAttrib).
  Type(Plant_)      :: Plant      ! Plant information (dummy variable for call to GetAttrib).
  Real(Std)         :: ConvFactor ! Conversion factor for air density
                                  ! (kg/m3 -> molecules/cm3).
  Type(ShortTime_)  :: TravelTime ! Dummy variable for TravelTime argument in GetAttrib.

  ! Calculate factor for converting air density from kg/m3 to molecules/cm3.
  ConvFactor = Avogadro / (Real(RMM_AIR)*1.0E6)

  ! Reset the flow memory prior to the met retrievals.
  Call ResetFlowMemory(Flows, FlowMemory)

  ! Get flow information.
  !$$ Check flag settings here.
  Call GetAttrib(                       &
         A_Flow,                        &
         Coords, Grids, Domains,        &
         Moisture      = .true.,        &
         Inhomog       = .true.,        &
         Homog         = .false.,       &
         Time          = Time,          &
         AnyTravelTime = .true.,        &
         TravelTime    = TravelTime,    &
         Position      = Position,      &
         Units         = Units,         &
         Mets          = Mets,          &
         Flows         = Flows,         &
         FlowMemory    = FlowMemory,    &
         Flow          = Flow,          &
         Cloud         = Cloud,         &
         Rain          = Rain,          &
         Surface       = Surface,       &
         Soil          = Soil,          &
         Plant         = Plant,         &
         ErrorCode     = Err            &
       )
  If (Err /= 0) Then
    Call Message('ERROR in GetMetForChemistry: problem while retrieving flow', 3)
  End If

  ! Get cloud information.
  Call GetAttrib(                            &
         A_Cloud,                            &
         Coords, Grids, Domains,             &
         .false., .false., .false.,          &
         Time, .true., TravelTime, Position, &
         Units, Mets, Flows,                 &
         FlowMemory,                         &
         Flow, Cloud, Rain,                  &
         Surface, Soil, Plant,               &
         Err                                 &
       )
  If (Err /= 0) Then
    Call Message('ERROR in GetMetForChemistry: problem while retrieving cloud', 3)
  End If

  ! Assign met variables.
  T                    = Flow%T
  Q                    = Flow%Q
  P                    = Flow%P
  MolecularAirDensity  = Flow%Rho * ConvFactor
  H                    = Flow%H
  Topog                = Flow%Topog           ! $$ Use land-use data when available?
  TotalOrDynCloudWater = Cloud%TotalOrDynCloudWater
  TotalOrDynCloudIce   = Cloud%TotalOrDynCloudIce
  CloudFraction        = Cloud%Cloud3d

  !$$ Total column cloud fraction currently set to gridbox value. We need to calculate
  !$$ this but it's probably best to do it in the GetAttrib routine - check with Dave?
  !$$ Also consider: 1) overlap assumptions? (current chemistry scheme uses Max val)
  !$$                2) are total cloud columns available directly from the UM?
  !
  ! $$ I think the total column means total above the particle (DJT) - see metchem.f
  !    If so need to think how to calculate this (but this should be responsibility
  !    of NWPMet/Flow).
  TotalColumnCloudFraction = CloudFraction

End Subroutine GetMetForChemistry

!-------------------------------------------------------------------------------------------------------------

Subroutine ConvertMassToConc(Conc, iX, iY, iZ, ChemOpts, ChemistryDefn, ChemistryState)
! Converts total gridbox masses (in g) to gridbox concentrations (in molecules/cm3)
! for species held on particles (provides the same functionality as the subroutine
! CALC_CONC in NAME V8.08). The routine is also able to convert any radiological species
! given in Bq units.

  Implicit None
  ! Argument list:
  Real(8),               Intent(Out) :: Conc(MaxSpecieses)
  Integer,               Intent(In)  :: iX
  Integer,               Intent(In)  :: iY
  Integer,               Intent(In)  :: iZ
  Type(ChemOpts_),       Intent(In)  :: ChemOpts
  Type(ChemistryDefn_),  Intent(In)  :: ChemistryDefn
  Type(ChemistryState_), Intent(In)  :: ChemistryState
  ! Conc           :: Species concentrations in chemistry gridbox.
  ! iX             ::} Indices of chemistry gridbox.
  ! iY             ::}
  ! iZ             ::}
  ! ChemistryDefn  :: A chemistry definition.
  ! ChemistryState :: State of a chemistry calculation.
  ! Locals:
  Integer :: iSpecies           ! Loop index.
  Real(8) :: ConvFactor         ! Conversion factor (moles -> molecules/cm3).
  Real(8) :: ConvFactorBq       ! Conversion factor (Bq -> g).

  ! Note: calculation of gridbox volume used to be done in this routine in
  ! NAME V8.08 - this is now performed in the InitChemistryState subroutine.

  ConvFactor   = Dble(Avogadro) / ChemistryState%Volume(iX, iY, iZ)

 ! Convert from mass in g (or Bq) to concentrations in molecules/cm3 

  Do iSpecies = 1, ChemistryDefn%nParticleSpecies
    If (ChemistryState%Mass(iSpecies, iZ, iY, iX) > 0.0) Then

      If (ChemistryDefn%ChangeUnitsParticle(iSpecies)) Then

        ! Convert from Bq to g
        ConvFactorBq = ChemistryDefn%ParticleSpeciesMolWeights(iSpecies) /                              &
                       (Log(2.0) * ChemistryDefn%ParticleSpeciesInvHalfLife(iSpecies) * Dble(Avogadro))

        Conc(iSpecies) = (ChemistryState%Mass(iSpecies, iZ, iY, iX) * ConvFactorBq /      &
                         ChemistryDefn%ParticleSpeciesMolWeights(iSpecies)) * ConvFactor 

      Else

        Conc(iSpecies) = (ChemistryState%Mass(iSpecies, iZ, iY, iX) /                     &
                         ChemistryDefn%ParticleSpeciesMolWeights(iSpecies)) * ConvFactor

      Endif

    Else
      Conc(iSpecies) = 0.0
    End If
  End Do

End Subroutine ConvertMassToConc

!-------------------------------------------------------------------------------------------------------------

Subroutine ConvertConcToMass(Conc, iX, iY, iZ, ChemOpts, ChemistryDefn, ChemistryState)
! Converts gridbox concentrations (in molecules/cm3) to total gridbox masses (in g)
! for species held on particles (provides the same functionality as the subroutine
! CALC_MASS in NAME V8.08). The routine is also able to convert back any radiological species
! given in Bq units.

  Implicit None
  ! Argument list:
  Real(8),               Intent(In)    :: Conc(MaxSpecieses)
  Integer,               Intent(In)    :: iX
  Integer,               Intent(In)    :: iY
  Integer,               Intent(In)    :: iZ
  Type(ChemOpts_),       Intent(In)    :: ChemOpts
  Type(ChemistryDefn_),  Intent(In)    :: ChemistryDefn
  Type(ChemistryState_), Intent(InOut) :: ChemistryState
  ! Conc           :: Species concentrations in chemistry gridbox.
  ! iX             ::} Indices of chemistry gridbox.
  ! iY             ::}
  ! iZ             ::}
  ! ChemistryDefn  :: A chemistry definition.
  ! ChemistryState :: State of a chemistry calculation.
  ! Locals:
  Integer :: iSpecies           ! Loop index.
  Real(8) :: ConvFactor         ! Conversion factor (molecules/cm3 -> moles).
  Real(8) :: ConvFactorBq       ! Conversion factor (g -> Bq).

  ConvFactor         = ChemistryState%Volume(iX, iY, iZ) / Dble(Avogadro)

  ! Concentrations are in molecules/cm3 - need to convert back to g (or Bq).

  Do iSpecies = 1, ChemistryDefn%nParticleSpecies

    If (ChemistryDefn%ChangeUnitsParticle(iSpecies)) Then

      ! If the original units were Bq, need to convert from mass back to Bq
      ConvFactorBq = (Log(2.0) * ChemistryDefn%ParticleSpeciesInvHalfLife(iSpecies) * Dble(Avogadro)) / &
                     ChemistryDefn%ParticleSpeciesMolWeights(iSpecies)

      ChemistryState%Mass(iSpecies, iZ, iY, iX) =                                                       &
        Conc(iSpecies) * ConvFactor * ChemistryDefn%ParticleSpeciesMolWeights(iSpecies) * ConvFactorBq

    Else

      ChemistryState%Mass(iSpecies, iZ, iY, iX) =                                       &
        Conc(iSpecies) * ConvFactor * ChemistryDefn%ParticleSpeciesMolWeights(iSpecies)

    End If

  End Do

End Subroutine ConvertConcToMass

!-------------------------------------------------------------------------------------------------------------

Subroutine ConvertChemFieldToConc(Conc, ConcBackFields, iX, iY, iZ, ChemOpts, ChemistryDefn, &
                                  ChemistryState, EulerianField)
! Converts gridbox concentrations from g/m3 to molecules/cm3 for species stored on chemistry fields
! (provides the same functionality as the subroutine GETCHEM in NAME V8.17). The routine is also able
! to convert any radiological species given in Bq units.

  Implicit None
  ! Argument list:
  Real(8),               Intent(Out) :: Conc(MaxSpecieses)
  Real(8),               Intent(Out) :: ConcBackFields(MaxSpecieses)
  Integer,               Intent(In)  :: iX
  Integer,               Intent(In)  :: iY
  Integer,               Intent(In)  :: iZ
  Type(ChemOpts_),       Intent(In)  :: ChemOpts
  Type(ChemistryDefn_),  Intent(In)  :: ChemistryDefn
  Type(ChemistryState_), Intent(In)  :: ChemistryState
  Type(EulerianField_),  Intent(In)  :: EulerianField
  ! Conc                :: Species concentrations in chemistry gridbox: elements 1 to
  !                        ChemistryDefn%nFieldSpecies store species held on fields.
  ! iX                  ::} Indices of chemistry gridbox.
  ! iY                  ::}
  ! iZ                  ::}
  ! ChemOpts            :: Chemistry options.
  ! ChemistryDefn       :: A chemistry definition.
  ! ChemistryState      :: State of a chemistry calculation.
  ! Locals:
  Integer :: iChemField   ! Loop index.
  Integer :: i            ! Loop index.
  Real(8) :: ConvFactor   ! Conversion factor (moles/m3 -> molecules/cm3).
  Real(8) :: ConvFactorBq ! Conversion factor (Bq/m3 -> g/m3).

  ! Concentrations needed in molecules/cm3 - convert from g/m3 by
  !  a) dividing by molecular weight of species, then
  !  b) multiplying by conversion factor (moles/m3 -> molecules/cm3).
  ConvFactor = AVOGAD / 1.0E6

  Do iChemField = 1, ChemistryDefn%nFieldSpecies

    If (ChemistryDefn%ChangeUnitsField(iChemField)) Then

      ! Convert from Bq/m3 to g/m3
      ConvFactorBq = ChemistryDefn%FieldSpeciesMolWeights(iChemField) /                              &
                     (Log(2.0) * ChemistryDefn%FieldSpeciesInvHalfLife(iChemField) * Dble(Avogadro))

      Conc(iChemField) = (EulerianField%Concentration(iX, iY, iZ, iChemField, EulerianField%iNew) * ConvFactorBq /  &
                         ChemistryDefn%FieldSpeciesMolWeights(iChemField)) * ConvFactor

    Else

      Conc(iChemField) = (EulerianField%Concentration(iX, iY, iZ, iChemField, EulerianField%iNew) / &
                         ChemistryDefn%FieldSpeciesMolWeights(iChemField)) * ConvFactor
    End If

  End Do

  ! Background fields. $$ Only needed for fields relaxed to. Not needed if only used for initialisation.
  Do i = 1, ChemistryDefn%nBackFields    
    ConcBackFields(i) = ChemistryState%ChemBackField(iX, iY, iZ, i) * ConvFactor   &
                              / ChemistryDefn%BackFieldSpeciesMolWeights(i)
  End Do


End Subroutine ConvertChemFieldToConc

!-------------------------------------------------------------------------------------------------------------

Subroutine ConvertConcToChemField(Conc, iX, iY, iZ, ChemOpts, ChemistryDefn, ChemistryState, EulerianField)
! Converts gridbox concentrations from molecules/cm3 to g/m3 for species stored on chemistry fields
! and updates chemistry fields (provides the same functionality as the subroutine NEW_CHEM in NAME V8.17).
! The routine is also able to convert back any radiological species given in Bq units.

  Implicit None
  ! Argument list:
  Real(8),               Intent(In)    :: Conc(MaxSpecieses)
  Integer,               Intent(In)    :: iX
  Integer,               Intent(In)    :: iY
  Integer,               Intent(In)    :: iZ
  Type(ChemOpts_),       Intent(In)    :: ChemOpts
  Type(ChemistryDefn_),  Intent(In)    :: ChemistryDefn
  Type(ChemistryState_), Intent(In)    :: ChemistryState
  Type(EulerianField_),  Intent(InOut) :: EulerianField
  ! Conc                :: Species concentrations in chemistry gridbox: elements 1 to
  !                        ChemistryDefn%nFieldSpecies store species held on fields.
  ! iX                  ::} Indices of chemistry gridbox.
  ! iY                  ::}
  ! iZ                  ::}
  ! ChemOpts            :: Chemistry options.
  ! ChemistryDefn       :: A chemistry definition.
  ! ChemistryState      :: State of a chemistry calculation.
  ! Locals:
  Integer  :: iChemField   ! Loop index.
  Real(8)  :: ConvFactor   ! Conversion factor (molecules/cm3 -> moles/m3).
  Real(8)  :: ConvFactorBq ! Conversion factor (g/m3 -> Bq/m3).

  ! Concentrations are in molecules/cm3 - need to convert to g/m3 by
  !  a) multiplying by conversion factor (molecules/cm3 -> moles/m3), then
  !  b) multiplying by molecular weight of species.
  ConvFactor = 1.0E6 / AVOGAD

  Do iChemField = 1, ChemistryDefn%nFieldSpecies

    If (ChemistryDefn%ChangeUnitsField(iChemField)) Then

      ! Convert from g/m3 back to Bq/m3
      ConvFactorBq = (Log(2.0) * ChemistryDefn%FieldSpeciesInvHalfLife(iChemField) * Dble(Avogadro)) / &
                     ChemistryDefn%FieldSpeciesMolWeights(iChemField)

      EulerianField%Concentration(iX, iY, iZ, iChemField, EulerianField%iNew) =          &
        Conc(iChemField) * ConvFactor * ChemistryDefn%FieldSpeciesMolWeights(iChemField) * ConvFactorBq

    Else

      EulerianField%Concentration(iX, iY, iZ, iChemField, EulerianField%iNew) =          &
        Conc(iChemField) * ConvFactor * ChemistryDefn%FieldSpeciesMolWeights(iChemField)

    End If

  End Do

End Subroutine ConvertConcToChemField

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateParticleMasses(                      &
             ChemOpts, ChemistryDefn, ChemistryState, &
             iX, iY, iZ,                              &
             Particles, Masses                        &
           )
! Reassigns the total grid box masses back to the particles in a grid box.
! (adapted from the subroutine ASSIGNMASS in NAME V8.08 - but now the mass
! reassignment is performed on each grid box individually).
! H2S added by Alex Archibald, University of Bristol

  Implicit None
  ! Argument list:
  Type(ChemOpts_),       Intent(In)            :: ChemOpts
  Type(ChemistryDefn_),  Intent(In)            :: ChemistryDefn
  Type(ChemistryState_), Intent(In)            :: ChemistryState
  Integer,               Intent(In)            :: iX
  Integer,               Intent(In)            :: iY
  Integer,               Intent(In)            :: iZ
  Type(Particle_),       Intent(InOut), Target :: Particles(:)
  Real(Std),             Intent(InOut)         :: Masses(:, :)
  ! ChemistryDefn  :: A chemistry definition.
  ! ChemistryState :: State of a chemistry calculation.
  ! iX             :: } Indices of chemistry gridbox.
  ! iY             :: }
  ! iZ             :: }
  ! Particles      :: Collection of particles currently defined in DispState.
  ! Masses         :: Masses of the various species carried by the particles.
  ! Locals:
  Logical   :: UpdateSpecies(MaxSpecieses)
  Real(8)   :: DeltaMass(MaxSpecieses)
  Real(8)   :: CurrentMass(MaxSpecieses)
  Real(8)   :: NegativeMass(MaxSpecieses)
  Real(8)   :: PrevParticleMass(MaxSpecieses)
  Real(8)   :: ProportionToParticle
  Real(8)   :: MassFraction
  Real(8)   :: ParticleMass
  Integer   :: i
  Integer   :: k
  Integer   :: iP
  Integer   :: iSpecies
  Integer   :: preSpecies
  Integer   :: iChemSpecies
  Integer   :: preChemSpecies
  Real(Std) :: ExcessNegativeMass
  Type(Particle_), Pointer :: Particle
  ! UpdateSpecies        :: Indicates that species needs updating.
  ! DeltaMass            :: Change in gridbox mass during chemistry.
  ! CurrentMass          :: Current mass reassigned to particles in gridbox.
  ! NegativeMass         :: Negative mass carried in calculation.
  ! PrevParticleMass     :: Initial particle mass.
  ! ProportionToParticle :: Proportion of gridbox mass assigned to a particle.
  ! MassFraction         :: Mass in a particle (relative to gridbox total).
  ! ParticleMass         :: Interim new mass of particle.
  ! i                    :: Index variable.
  ! k                    :: Specifies location of gridbox particles in ParticleIndex array.
  ! iP                   :: Particle index.
  ! iSpecies             ::} Species indices.
  ! preSpecies           ::}
  ! iChemSpecies         :] Species indices in the chemistry schemes list of species.
  ! preChemSpecies       :]
  ! ExcessNegativeMass   :: Excess negative mass in a gridbox.
  ! Particle             :: A particle.
  !
  ! MassReassignment     :: Mass reassignment rules: these are hard wired in ChemistryScheme.F90 to suit the
  !                         chemistry scheme being used
  ! $$ should check somewhere that array sums are 1.0.

  ! Check whether there are any particles in the chemistry grid box.
  ! Nothing to do if there are no particles in the grid box.
  If (ChemistryState%GridboxState(iX, iY, iZ)%n == 0) Return

  ! Set species update flags (avoids treating species that do not change).
  UpdateSpecies = .true.

  !Calculate changes in gridbox masses during chemistry.
  Do iSpecies = 1, ChemistryDefn%nParticleSpecies
    DeltaMass(iSpecies) = ChemistryState%Mass(iSpecies, iZ, iY, iX) -   &
                          ChemistryState%OldMass(iSpecies, iZ, iY, iX)
!   If (Abs(DeltaMass(iSpecies)) <= 1.0E-4) DeltaMass(iSpecies) = 0.0
    NegativeMass(iSpecies) = 0.0
    CurrentMass(iSpecies) = 0.0
    If (DeltaMass(iSpecies) == 0.0) UpdateSpecies(iSpecies) = .false.
  End Do

  ! Loop through all particles in the chemistry grid box. These are identified
  ! as elements k-1 to k-n of the ParticleIndex array.
  k = ChemistryState%GridboxState(iX, iY, iZ)%k
  Do i = 1, ChemistryState%GridboxState(iX, iY, iZ)%n
    ! Identify particle
    iP = ChemistryState%ParticleIndex(k-i)
    Particle => Particles(iP)

    ! Store previous masses of the particle
    Do iSpecies = 1, ChemistryDefn%nParticleSpecies
      PrevParticleMass(iSpecies) = Dble(Masses(iSpecies, iP))
    End Do

    ! Loop through all species and update particle masses
    Do iChemSpecies = 1, ChemistryDefn%nChemSpecieses

      If (ChemistryDefn%iChem2Particle(iChemSpecies) == 0) Cycle
      iSpecies = ChemistryDefn%iChem2Particle(iChemSpecies)

      If (.not.UpdateSpecies(iSpecies)) Cycle

      If (DistributeUniformly(iChemSpecies)) Then
        ! Distribute evenly between all particles in grid box.
        Masses(iSpecies, iP) = ChemistryState%Mass(iSpecies, iZ, iY, iX) /  &
                               ChemistryState%GridboxState(iX, iY, iZ)%n
        Cycle
      End If

      ! Calculate the proportion of the gridbox species to assign to this particular
      ! particle - this factor is based on the relative mass of its precursor species.
      ProportionToParticle = 0.0
      Do preChemSpecies = 1, ChemistryDefn%nChemSpecieses

        If (MassReassignment(preChemSpecies, iChemSpecies) == 0.0) Cycle

        ! Precursor not on particles - distribute product evenly over particles.
        If (ChemistryDefn%iChem2Particle(preChemSpecies) == 0) Then

          ProportionToParticle = ProportionToParticle +                           &
                                 MassReassignment(preChemSpecies, iChemSpecies) / &
                                 ChemistryState%GridboxState(iX, iY, iZ)%n

        ! Precursor on particles - distribute product in proportion to precursor mass.
        Else

          preSpecies = ChemistryDefn%iChem2Particle(preChemSpecies)

          If (PrevParticleMass(preSpecies) > 0.0) Then
            MassFraction         = PrevParticleMass(preSpecies) /                     &
                                   ChemistryState%OldMass(preSpecies, iZ, iY, iX)
            ProportionToParticle = ProportionToParticle +                             &
                                   MassReassignment(preChemSpecies, iChemSpecies) * MassFraction
          End If

        End If

      End Do

      ! $$ I'd like to delete this but just commented out for now (DJT). While it might trap a badly
      ! set up MassReassignment array (values adding to more than 1) its not the best way to do that.
      ! Check that the proportionality factor is <= 1
      ! (trapping out rounding errors first).
      ! If (ProportionToParticle > 0.999999 .and. ProportionToParticle < 1.000001) &
      !   ProportionToParticle = 1.0
      ! If (ProportionToParticle > 1.0) Then
      !   Call Message(                                                            &
      !      'WARNING in UpdateParticleMasses: Proportion of mass assigned to ' // &
      !      'particle exceeds unity - proportionality factor is '              // &
      !      Std2Char(Real(ProportionToParticle), 12)                           // &
      !      ' for species '                                                    // &
      !      Int2Char(iSpecies),                                                   &
      !      1                                                                     &
      !   )
      ! End If

      ! Calculate the interim particle mass (before any negative mass is treated).
      ParticleMass = PrevParticleMass(iSpecies) +                 &
                       ProportionToParticle * DeltaMass(iSpecies)
      ! If particle mass becomes negative, set to zero and carry negative component.
      If (ParticleMass < 0.0) Then
        NegativeMass(iSpecies)  = NegativeMass(iSpecies) + ParticleMass
        ParticleMass            = 0.0
      End If
      ! Update particle mass.
      Masses(iSpecies, iP) = Real(ParticleMass)
      CurrentMass(iSpecies) = CurrentMass(iSpecies) + ParticleMass
    End Do
  End Do

  ! Redistribute any negative mass carried forward from the mass update procedure.
  ! $$ Currently eliminates mass on a preferential basis (rather than redistributing
  ! $$ mass equally between particles).
  ! $$ Better would be to calculate a second DeltaMass and assign that based on Mass(species) 
  ! not Mass(prespecies). 
  ! $$ Note I don't think CurrentMass is used at all.
  Do iSpecies = 1, ChemistryDefn%nParticleSpecies
    i = 1
    Do While (NegativeMass(iSpecies) < 0.0)
      If (i <= ChemistryState%GridboxState(iX, iY, iZ)%n) Then
        ! Identify particle
        iP = ChemistryState%ParticleIndex(k-i)
        Particle => Particles(iP)
        ! If particle carries species then allocate negative mass to this particle
        If (Masses(iSpecies, iP) > 0.0) Then
          ParticleMass = Dble(Masses(iSpecies, iP))
          If (ParticleMass + NegativeMass(iSpecies) < 0.0) Then
            NegativeMass(iSpecies) = NegativeMass(iSpecies) + ParticleMass
            CurrentMass(iSpecies)  = CurrentMass(iSpecies)  - ParticleMass
            ParticleMass           = 0.0
          Else
            ParticleMass           = ParticleMass          + NegativeMass(iSpecies)
            CurrentMass(iSpecies)  = CurrentMass(iSpecies) + NegativeMass(iSpecies)
            NegativeMass(iSpecies) = 0.0
          End If
          Masses(iSpecies, iP) = Real(ParticleMass)
        End If
        i = i + 1
      Else
        ExcessNegativeMass = Real(NegativeMass(iSpecies))
        ! $$ I'd like to delete this bit (DJT). I think can only occur if all particles have
        ! mass set to zero for the species. One can't possibly then be overestimating the 
        ! concentration!
        If (ExcessNegativeMass < -0.05) Then
          Call Message(                                                     &
             'WARNING in UpdateParticleMasses: Excess negative mass of ' // &
             Std2Char(ExcessNegativeMass)                                // &
             ' g cannot be assigned to particles, for particle species ' // & ! $$ user might have trouble
             Int2Char(iSpecies),                                            & !    interpreting index
             1                                                              &
          )
        End If
        Exit
      End If
    End Do
  End Do

End Subroutine UpdateParticleMasses

!-------------------------------------------------------------------------------------------------------------
!$$ We use a separate ParticleMass within this routine for particle splitting
!$$ - since it's different to the 'release' ParticleMass. Here we only consider PM10
!$$ components (PM10, SO4, NH42SO4, NAER, NH4NO3) with all ParticleMass values set
!$$ to the input PM10 value. Either hardwire the new array here (as done) or specify
!$$ a second ParticleMass array in the Species file (more flexible)??? ARJ, 17/11/04.

Subroutine SplitMassiveParticles(                                                                        &
             Specieses, LastParticle, nParticles,                                                        &
             ParticleCeiling, ParticleFactor, MaxParticles,                                              &
             FreeParticleStack, Particles, Masses,                                                       &
             ParticleFactorType,                                                                         &
             nParticleExtras, LastParticleExtra, FreeParticleExtraStack, ParticleExtras, iParticleExtras &
           )
! Splits massive particles into multiple copies of lighter particles.

  Implicit None
  ! Argument list:
  Type(Specieses_), Intent(In)            :: Specieses
  Integer,          Intent(InOut)         :: LastParticle
  Integer,          Intent(InOut)         :: nParticles
  Integer,          Intent(In)            :: ParticleCeiling
  Real(Std),        Intent(In)            :: ParticleFactor
  Integer,          Intent(In)            :: MaxParticles
  Integer,          Intent(InOut)         :: FreeParticleStack(:)
  Type(Particle_),  Intent(InOut), Target :: Particles(:)
  Real(Std),        Intent(InOut)         :: Masses(:, :)
  Integer,          Intent(InOut)         :: ParticleFactorType

  Integer,          Intent(InOut)         :: nParticleExtras    ! $$ This bit needs tidying up.
  Integer,          Intent(InOut)         :: LastParticleExtra
  Integer,          Intent(InOut)         :: FreeParticleExtraStack(:)
  Integer,          Intent(InOut)         :: iParticleExtras(:) ! Indices of the Extra corresponding to
                                                                ! a given Particle.
  Type(Extra_),     Intent(InOut)         :: ParticleExtras(:)  ! Extra(0) is a default extra for
                                                                ! cheap particles, rest are
                                                                ! useable for expensive particles.

  ! Specieses          :: Collection of species.
  ! LastParticle       :: Index of last particle.
  ! nParticles         :: Number of active particles.
  ! ParticleCeiling    :: Particle number ceiling above which the number of particles
  !                       created in the splitting and their masses are adjusted
  !                       to try to prevent the model running out of particles.
  ! ParticleFactor     :: Factor by which the mass per particle is increased above
  !                       ParticleCeiling.
  ! MaxParticles       :: Maximum number of particles allowed.
  ! FreeParticleStack  :: Indices of free particles.
  ! Particles          :: Collection of all particles.
  ! Masses             :: Masses of the various species carried by the particles.
  ! ParticleFactorType :: Indicates type of adjustment to the splitting of particles
  !                       to try to prevent the model running out of particles:
  !                           0 - no adjustment
  !                           1 - ParticleCeiling reached, numbers reduced and masses
  !                               increased,
  !                           2 - MaxParticles reached and particle splitting stopped.
  ! Locals:
  Integer   :: i                            ! Particle loop index.
  Integer   :: j                            ! Particle loop index (for new particles).
  Integer   :: iSpecies                     ! Species index.
  Integer   :: iE                           ! Index of 'particle extra'.
  Integer   :: iLimSpecies                  ! Mass-limited species loop index.
  Integer   :: nLimitedSpecies              ! Number of mass-limited species.
  Integer   :: LimitedSpecies(MaxSpecieses) ! Indices of the mass-limited species.
  Integer   :: nNew                         ! Number of new particles for each particle.
  Integer   :: iNew                         ! Index of a new particle.
  Integer   :: nTotal                       ! Total number of new particles created.
  Integer   :: iPM10InSpecieses             ! Index of PM10 species in Specieses.
  Integer   :: iPM10                        ! Index of PM10 in species on particles.
  Integer   :: iSulphateAerosol             ! Index of sulphate aerosol in species on particles.
  Integer   :: iAmmoniumSulphate            ! Index of ammonium sulphate in species on particles.
  Integer   :: iNitrateAerosol              ! Index of nitrate aerosol in species on particles.
  Integer   :: iAmmoniumNitrate             ! Index of ammonium nitrate in species on particles.
  Real(Std) :: ParticleMass                 ! Limiting particle mass (= PM10 partmass).
  Real(Std) :: MassRatio                    ! Mass carried by particle relative to the
                                            ! desired mass limit.
  Integer   :: P1                           !} Integers used in adjusting the number
  Integer   :: P2                           !} of new particles when the particle
                                            !} ceiling is first exceeded.
  Integer   :: LastParticle1                ! Local copy of LastParticle.
  Integer   :: PFType                       ! Local copy of ParticleFactorType.
  Real(Std) :: MaxSulphateAerosol           !} Local diagnostics.
  Real(Std) :: MaxAmmoniumSulphate          !}
  Real(Std) :: MaxNitrateAerosol            !}
  Real(Std) :: MaxAmmoniumNitrate           !}
  Real(Std) :: MaxPM10                      !}
  Type(Particle_), Pointer :: Particle      ! A particle.
  Logical   :: Error                        !

  ! **** Code to set up mass limits via Species information (not currently used) ****
  !
  ! ! Identify those species which are mass-limited.
  ! nLimitedSpecies = 0
  ! Do iSpecies = 1, Specieses%nSpecieses
  !   If (Specieses%Specieses(iSpecies)%UseMassLimit) Then
  !     nLimitedSpecies                 = nLimitedSpecies + 1
  !     LimitedSpecies(nLimitedSpecies) = iSpecies
  !   End If
  ! End Do
  ! If (nLimitedSpecies == 0) Return
  !
  ! **** End of alternative code fragment ****

  ! Set up those species which are mass-limited (hard-wired here).
  ! Note FindSpeciesIndex returns zero and sets Error = .true. if the species isn't found.

  iPM10InSpecieses    = FindSpeciesIndex('PM10', Specieses, Error)
  If (iPM10InSpecieses == 0) Return ! $$ Doesn't work without PM10 present.
  
  nLimitedSpecies     = 5
  iPM10               = Specieses%iSpecies2Particle( FindSpeciesIndex('PM10',     Specieses, Error) )
  iSulphateAerosol    = Specieses%iSpecies2Particle( FindSpeciesIndex('SULPHATE', Specieses, Error) )
  iAmmoniumSulphate   = Specieses%iSpecies2Particle( FindSpeciesIndex('NH42SO4',  Specieses, Error) )
  iNitrateAerosol     = Specieses%iSpecies2Particle( FindSpeciesIndex('NAER',     Specieses, Error) )
  iAmmoniumNitrate    = Specieses%iSpecies2Particle( FindSpeciesIndex('NH4NO3',   Specieses, Error) )
  LimitedSpecies(1:5) = (/                    &
                           iPM10,             &
                           iSulphateAerosol,  &
                           iAmmoniumSulphate, &
                           iNitrateAerosol,   &
                           iAmmoniumNitrate   &
                        /)

  If (iPM10 == 0) Return ! $$ Doesn't work without PM10 present.

  If (.not.Specieses%Specieses(iPM10InSpecieses)%UseMassLimit) Return

  ParticleMass = Specieses%Specieses(iPM10InSpecieses)%MassLimit

  ! Initialise local diagnostics.
  nTotal = 0
  MaxSulphateAerosol  = 0.0
  MaxAmmoniumSulphate = 0.0
  MaxNitrateAerosol   = 0.0
  MaxAmmoniumNitrate  = 0.0
  MaxPM10             = 0.0

  ! New variable to keep track of revised value of LastParticle.
  LastParticle1 = LastParticle

  ! Set PFType to zero.
  PFType = 0

  ! Create new particles if the masses carried by any existing particle are too large.
  Do i = 1, LastParticle
    Particle => Particles(i)
    If (.not.ParticleActive(Particle)) Cycle
    ! local diagnostics
    If (iSulphateAerosol  /= 0) MaxSulphateAerosol  = Max(MaxSulphateAerosol,  Masses(iSulphateAerosol,  i))
    If (iAmmoniumSulphate /= 0) MaxAmmoniumSulphate = Max(MaxAmmoniumSulphate, Masses(iAmmoniumSulphate, i))
    If (iNitrateAerosol   /= 0) MaxNitrateAerosol   = Max(MaxNitrateAerosol,   Masses(iNitrateAerosol,   i))
    If (iAmmoniumNitrate  /= 0) MaxAmmoniumNitrate  = Max(MaxAmmoniumNitrate,  Masses(iAmmoniumNitrate,  i))
    If (iPM10             /= 0) MaxPM10             = Max(MaxPM10,             Masses(iPM10,             i))

    ! determine desired number of new particles (assuming no further adjustments are needed)
    nNew = 0
    Do iLimSpecies = 1, nLimitedSpecies
      iSpecies = LimitedSpecies(iLimSpecies)
      If (iSpecies == 0) Cycle
      MassRatio = Masses(iSpecies, i) / ParticleMass
      nNew = Max(nNew, Ceiling(MassRatio) - 1)
    End Do
    If (nNew == 0) Cycle
    ! check if the new particles will exceed the particle ceiling and adjust their
    ! number accordingly. Also update PFType.
    P1 = nParticles        - ParticleCeiling
    P2 = nParticles + nNew - ParticleCeiling
    If (P1 >= 0) Then
      nNew = Ceiling(Real(nNew, Std) / ParticleFactor)
      PFType = 1
    Else If (P2 > 0) Then
      nNew = Ceiling(Real(P2, Std) / ParticleFactor - Real(P1, Std))
      PFType = 1
    End If
    ! check if adding the new particles will exceed MaxParticles and adjust their
    ! number accordingly. Also update PFType.
    If (nParticles + nNew >= MaxParticles) Then
      nNew = MaxParticles - nParticles
      PFType = 2
    End If
    ! reset particle mass array to fraction of its current value
    Masses(:, i) = Masses(:, i) / Real(nNew+1)
    ! create new particles as duplicates of the mass-reduced particle
    Do j = 1, nNew
      nParticles = nParticles + 1
      iNew = FreeParticleStack(nParticles)
      Particles(iNew) = Particles(i)
      Masses(:, iNew) = Masses(:, i)
      If (iNew > LastParticle1) LastParticle1 = iNew
      ! create new ParticleExtra if required
      If (iParticleExtras(i) /= 0) Then
        If (nParticleExtras + 1 > Size(FreeParticleExtraStack)) Then
          Call Message('FATAL ERROR: too many "full" particles', 3)
        End If
        nParticleExtras = nParticleExtras + 1
        iE = FreeParticleExtraStack(nParticleExtras)
        If (iE > LastParticleExtra) LastParticleExtra = iE
        iParticleExtras(iNew) = iE
        ParticleExtras(iParticleExtras(iNew)) = ParticleExtras(iParticleExtras(i))
      Else
        iParticleExtras(iNew) = 0
      End If
    End Do
    ! update total number of new particles added
    nTotal = nTotal + nNew
    ! exit the loop if the MaxParticles particle limit has been reached
    If (PFType == 2) Exit
  End Do
  ! Warning messages for reductions in number of particles created.
  If (PFType /= ParticleFactorType) Then
    If (PFType == 0) Then
      Call Message('Particle splitting returned to normal (this will apply to particle releases too)')
    Else If (PFType == 1 .and. ParticleFactorType == 0) Then
      Call Message(                                                              &
             'WARNING: Particle splitting reduced as nearing particle limit ' // &
             '(this will apply to particle releases too)',                       &
             1                                                                   &
           )
    Else If (PFType == 1 .and. ParticleFactorType == 2) Then
      Call Message(                                                                             &
             'Particle splitting restarted but at a reduced rate as nearing particle limit ' // &
             '(this will apply to particle releases too)'                                       &
           )
    Else If (PFType == 2) Then
      Call Message(                                                              &
             'WARNING: Particle splitting stopped as particle limit reached ' // &
             '(this will apply to particle releases too)',                       &
             1                                                                   &
           )
    End If
  End If
  ! Update the values of LastParticle and ParticleFactorType.
  LastParticle       = LastParticle1
  ParticleFactorType = PFType
  ! Output local diagnostics.
  If (nTotal > 0) Then
    Call Message(                                                           &
       'PARTICLE SPLITTING DIAGNOSTICS: particle splitting has created ' // &
       Int2Char(nTotal)                                                  // &
       ' new particles; '                                                // &
       'max masses initially carried on a single particle were '         // &
       Std2Char(MaxNitrateAerosol)                                       // &
       ' (NAER), '                                                       // &
       Std2Char(MaxAmmoniumNitrate)                                      // &
       ' (NH4NO3), '                                                     // &
       Std2Char(MaxSulphateAerosol)                                      // &
       ' (SO4), '                                                        // &
       Std2Char(MaxPM10)                                                 // &
       ' (PM10), '                                                       // &
       Std2Char(MaxAmmoniumSulphate)                                     // &
       ' ((NH4)2SO4)',                                                      &
       1                                                                    &
    )
  End If

End Subroutine SplitMassiveParticles

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateChemistry(                                                         &
             Coords, Grids, Domains, Units, Mets, Flows, Specieses, PrevTime, Time, &
             SyncInterval, nParticles, Particles, Masses,                           &
             ChemOpts, ChemistryDefn, ChemistryState, EulerianField,                &
             OpenMPOpts                                                             &
           )
! Calls the chemistry scheme to update particle species and chemistry fields.

  Use OpenMPModule, Only : OpenMPOpts_

  Implicit None
  ! Argument list:
  Type(Coords_),         Intent(In)    :: Coords
  Type(Grids_),          Intent(In)    :: Grids
  Type(Domains_),        Intent(In)    :: Domains
  Type(Units_),          Intent(InOut) :: Units
  Type(Mets_),           Intent(InOut) :: Mets
  Type(Flows_),          Intent(InOut) :: Flows
  Type(Specieses_),      Intent(In)    :: Specieses
  Type(Time_),           Intent(In)    :: PrevTime
  Type(Time_),           Intent(In)    :: Time
  Real(Std),             Intent(In)    :: SyncInterval
  Integer,               Intent(In)    :: nParticles
  Type(Particle_),       Intent(InOut) :: Particles(:)
  Real(Std),             Intent(InOut) :: Masses(:, :)
  Type(ChemOpts_),       Intent(In)    :: ChemOpts
  Type(ChemistryDefn_),  Intent(In)    :: ChemistryDefn
  Type(ChemistryState_), Intent(InOut) :: ChemistryState
  Type(EulerianField_),  Intent(InOut) :: EulerianField
  Type(OpenMPOpts_),     Intent(In)    :: OpenMPOpts
  ! Coords         :: Collection of coord systems.
  ! Grids          :: Collection of grids.
  ! Domains        :: Collection of domains.
  ! Units          :: Collection of information on input/output unit numbers.
  ! Mets           :: Collection of met module instance states.
  ! Flows          :: Collection of flows.
  ! Specieses      :: Collection of specieses.
  ! PrevTime       :: Previous time (start of latest synchronisation interval).
  ! Time           :: Current time (end of latest synchronisation interval).
  ! SyncInterval   :: Duration of synchronisation interval (in seconds).
  ! nParticles     :: Number of active particles currently defined in DispState.
  ! Particles      :: Collection of particles currently defined in DispState.
  ! Masses         :: Masses of the various species carried by the particles.
  ! ChemOpts       :: Chemistry options.
  ! ChemistryDefn  :: A chemistry definition.
  ! ChemistryState :: State of a chemistry calculation.
  ! OpenMPOpts     :: OpenMP options
  ! Locals:
  Integer          :: i
  Integer          :: j
  Integer          :: iLatLong
  Type(ShortTime_) :: MetTime
  Type(HGrid_)     :: HGrid
  Type(ZGrid_)     :: ZGrid
  Integer          :: iX
  Integer          :: iY
  Integer          :: iZ
  Type(Position_)  :: Position
  Real(Std)        :: X(3)
  Real(Std)        :: ZenithAngle
  Real(Std)        :: T
  Real(Std)        :: Q
  Real(Std)        :: P
  Real(Std)        :: MolecularAirDensity
  Real(Std)        :: H
  Real(Std)        :: Topog
  Real(Std)        :: TotalOrDynCloudWater
  Real(Std)        :: TotalOrDynCloudIce
  Real(Std)        :: CloudFraction
  Real(Std)        :: TotalColumnCloudFraction
  Real(8)          :: ConcSpeciesOnParticles(MaxSpecieses)
  Real(8)          :: ConcSpeciesOnFields(MaxSpecieses)
  Real(8)          :: ConcSpeciesOnBackFields(MaxSpecieses)
  Real(8)          :: ConcChemistrySpecies(MaxSpecieses)
  Integer          :: Err
  Integer          :: iBackO3
  ! i               :: Loop index.
  ! j               :: Loop index.
  ! iLatLong        :: Index (in Coords) of the standard lat-long coord system.
  ! MetTime         :: Time to be used for the met data retrieval.
  ! HGrid           ::} Horizontal and vertical grids used by the chemistry scheme.
  ! ZGrid           ::}
  ! iX              ::] Indices of a chemistry gridbox.
  ! iY              ::]
  ! iZ              ::]
  ! Position        :: Position of gridbox centre.
  ! X               :: Lat-long-height coords of gridbox centre.
  ! ZenithAngle     :: Solar zenith angle from gridbox centre.
  ! MET DATA FOR CHEMISTRY GRIDBOX:-
  !   T                        :: Ambient air temperature.
  !   Q                        :: Specific humidity.
  !   P                        :: Ambient pressure.
  !   MolecularAirDensity      :: Air density (in molecules/cm3).
  !   H                        :: Boundary layer height.
  !   Topog                    :: Topographic height above sea level.
  !   TotalOrDynCloudWater     :: Total or dynamic cloud liquid water.
  !   TotalOrDynCloudIce       :: Total or dynamic cloud ice.
  !   CloudFraction            :: Total (dyn+conv) cloud fraction.
  !   TotalColumnCloudFraction :: Max cloud fraction in cloud column above gridpoint.
  ! SPECIES CONCENTRATIONS IN CHEMISTRY GRIDBOX:-
  !   ConcChemistrySpecies     :: ... of species needed by chemistry.
  ! SPECIES CONCENTRATIONS IN CHEMISTRY COLUMN:-
  !   ConcSpeciesOnParticles   :: ... of species held on particles.
  !   ConcSpeciesOnFields      :: ... of species held on chemistry fields.
  !   ConcSpeciesOnBackFields  :: ... of background species read in and held on chemistry fields.
  ! Err        :: Error code from deallocate statements.
  ! iBackO3    :: Index of background ozone field.

  ! Find index of the standard lat-long coord system (position of each gridbox in
  ! standard lat-long coordinates is required for computing solar zenith angle).
  iLatLong = FindHCoordIndex('Lat-Long', Coords)

  ! Set the time for the met data retrieval.
  !$$ NAME uses the end of the timestep - would the mid-point be better here?
  !$$ When chemistry update is called, the "current" time is at the end
  !$$ of the sync time step.
  MetTime = Time2ShortTime(Time)

  ! If any chemistry background field species are being initialised, update these species from
  ! background fields in first call after the start of a new month.

  ! Identify horizontal and vertical grids of the chemistry scheme.
  HGrid = Grids%HGrids(EulerianField%iHGrid)
  ZGrid = Grids%ZGrids(EulerianField%iZGrid)

  If (ChemistryState%InitialiseBackFields) Then
    ! Initialise the background fields from nearest STOCHEM grid point if required.
    If (ChemistryDefn%nBackFields >= 1) Then
      Call SetBackgroundFieldsFromSTOCHEM( &
             Coords, Grids, Domains,       &
             Time, InfPastTime(),          &
             ChemOpts, ChemistryDefn,      &
             Units, Mets, Flows,           &
             EulerianField,                &
             ChemistryState                &
           )
    End If
    ! Initialise chemistry fields from background fields.
    Do i = 1, ChemistryDefn%nBackFields
      Do j = 1, ChemistryDefn%nFieldSpecies
        If (ChemistryDefn%FieldSpeciesNames(j) .CIEq. ChemistryDefn%BackFieldSpeciesNames(i)) Then
          EulerianField%Concentration(1:HGrid%nX, 1:HGrid%nY, 1:ZGrid%nZ, j, EulerianField%iNew) =  &
            ChemistryState%ChemBackField(:, :, :, i)
          Exit
        End If
      End Do
    End Do
    ChemistryState%InitialiseBackFields = .false.
  Else
    If (ChemistryDefn%nBackFields >= 1) Then
      Call SetBackgroundFieldsFromSTOCHEM( &
             Coords, Grids, Domains,       &
             Time, PrevTime,               &
             ChemOpts, ChemistryDefn,      &
             Units, Mets, Flows,           &
             EulerianField,                &
             ChemistryState                &
           )
    End If
  End If

  ! Initialise the chemistry gridbox information (count number of particles in each
  ! gridbox, label these particles and generate gridbox masses).
  Call InitChemistryGridboxes(                     &
         Coords, Grids, Domains, Time, nParticles, &
         Particles, Masses,                        &
         ChemistryDefn,                            &
         Units, Mets, Flows,                       &
         EulerianField,                            &
         ChemistryState                            &
       )

  Call TimerOn(ChemistryLoopTimer)
  ! Main chemistry loop (loop over all chemistry gridboxes column-by-column).

  !$OMP PARALLEL DO                                           &
  !$OMP NUM_THREADS(OpenMPOpts%nChemistryThreads)             &
  !$OMP SCHEDULE(RUNTIME)                                     &
  !$OMP DEFAULT (NONE)                                        &
  !$OMP SHARED  (HGrid,ZGrid,Coords,iLatLong,Time,Flows,Grids,&
  !$OMP          Mets,Units,                                  &
  !$OMP          Domains,MetTime,                             &
  !$OMP          ChemOpts,ChemistryDefn,ChemistryState,       &
  !$OMP          EulerianField,OpenMPOpts,                    &
  !$OMP          SyncInterval,Particles,Masses,Specieses)     &
  !$OMP PRIVATE (iX,iY,iZ,Position,X,ZenithAngle,T,Q,P,       &
  !$OMP          MolecularAirDensity,H,Topog,                 &
  !$OMP          TotalOrDynCloudWater,TotalOrDynCloudIce,     &
  !$OMP          CloudFraction,TotalColumnCloudFraction,      &
  !$OMP          i,ConcChemistrySpecies,                      &
  !$OMP          ConcSpeciesOnParticles,                      &
  !$OMP          ConcSpeciesOnFields,                         &
  !$OMP          ConcSpeciesOnBackFields,iBackO3)

  Do iX = 1, HGrid%nX
  Do iY = 1, HGrid%nY
  Do iZ = 1, ZGrid%nZ
   
    ! Set position of centre of gridbox.
    Position = X2Position(                    &
                 Coords,                      &
                 (/                           &
                   HGrid%X(iX),               &
                   HGrid%Y(iY),               &
                   ZGrid%Z(iZ)                &
                 /),                          &
                 HGrid%iHCoord, ZGrid%iZCoord &
               )

    ! Calculate lat-long of gridbox (if not already available).
    ! $$ Also need to convert height to metres above ground -
    ! $$ used in the main chemistry code. This is done later on. But lat-long is
    ! $$ no longer required since the zenith angle is calculated externally now.
    ! $$ Remove the following four lines here??
    If (.not.Position%XYValid(iLatLong)) Then
      Call ConvertToH(Coords, iLatLong, Position)
    End If
    X = Position2X(Coords, Position, iLatLong, ZGrid%iZCoord)

    ! Calculate the solar zenith angle.
    ZenithAngle = CalcZenithAngle(Coords, Time, Position)

    ! Only get met data if there are particles in the chemistry gridbox.
    If (ChemistryState%GridboxState(iX, iY, iZ)%n /= 0) Then
      ! Retrieve met data for chemistry gridbox.
      Call GetMetForChemistry(                           &
              Coords, Grids, Domains, MetTime, Position, &
              Units, Mets, Flows,                        &
              T, Q, P, MolecularAirDensity, H, Topog,    &
              TotalOrDynCloudWater, TotalOrDynCloudIce,  &
              CloudFraction, TotalColumnCloudFraction    &
              )
      ! LatLong will be available from ZenithAngle call and magl from GetMetForChemistry.
      X = Position2X(Coords, Position, iLatLong, ChemistryDefn%iZCoordMagl)
    End If

    ! Determine gridbox concentrations (in molecules/cm3) for species held on
    ! the particles (converts from the total gridbox masses in g).
    Call ConvertMassToConc(ConcSpeciesOnParticles(:), iX, iY, iZ, ChemOpts, &
                           ChemistryDefn, ChemistryState)

    ! Extract gridbox concentrations (in molecules/cm3) for species stored on
    ! the chemistry fields (note: chem fields stored as g/m3).
    Call ConvertChemFieldToConc(ConcSpeciesOnFields(:), ConcSpeciesOnBackFields(:), iX, iY, iZ, &
                                ChemOpts, ChemistryDefn, ChemistryState, EulerianField)

    ! If there are no particles in a chemistry gridbox, then ChemistryScheme is not called and ...
    !  - H2O2 (if stored on fields) takes the same value as previous timestep,
    !  - O3 (if stored on fields and there is no O3 background field) takes same value as previous timestep,
    !  - O3 (if stored on fields and there is an O3 background field) is relaxed back to the background O3
    !        field with relaxation time of 3 hours,
    !  - all other species stored on fields are reset to zero.

    If (ChemistryState%GridboxState(iX, iY, iZ)%n == 0) Then

      iBackO3 = 0
      Do i = 1, ChemistryDefn%nBackFields
        If (ChemistryDefn%BackFieldSpeciesNames(i) .CIEq. 'O3') iBackO3 = i
      End Do

      If (iBackO3 /= 0) Then

        ! If this is the first step when there are no particles in the gridbox, then set the
        ! O3 increment to relax back to background O3 in 3 hours
        ! $$ Use a logical flag here instead? Or use an exponential relaxation
        ! (i.e. recalc DeltaO3 each time
        ! step) which would be simpler.

        ! concspeciesonfields has the species in the same order as the names in chemistrydefn,
        ! so this can be used to find out which species is ozone

        Do i = 1, Specieses%nFields

          If (ChemistryDefn%FieldSpeciesNames(i) .CIEq. 'O3') Then 
       
            If (ChemistryState%GridboxState(iX, iY, iZ)%DeltaO3 == 0.0) Then
              ChemistryState%GridboxState(iX, iY, iZ)%DeltaO3 =             &
                (ConcSpeciesOnBackFields(iBackO3) - ConcSpeciesOnFields(i)) &
                * DBLE(SyncInterval / (3.0 * 3600.0))
            End If

            ! Increment O3 field.
            ConcSpeciesOnFields(i) = ConcSpeciesOnFields(i) +                     &
                                     ChemistryState%GridboxState(iX, iY, iZ)%DeltaO3

            ! Avoid O3 over-shooting the background ozone value.
            If (ChemistryState%GridboxState(iX, iY, iZ)%DeltaO3 > 0.0) Then
              If (ConcSpeciesOnFields(i) > ConcSpeciesOnBackFields(iBackO3)) Then
                ConcSpeciesOnFields(i) = ConcSpeciesOnBackFields(iBackO3)
              End If
            Else
              If (ConcSpeciesOnFields(i) < ConcSpeciesOnBackFields(iBackO3)) Then
                ConcSpeciesOnFields(i) = ConcSpeciesOnBackFields(iBackO3)
              End If
            End If
          End If

        End Do

      End If

      ! Species carried on fields other than H2O2 and ozone are reset to zero.

      Do i = 1, Specieses%nFields
        If (ChemistryDefn%FieldSpeciesNames(i) .CIEq. 'O3') Then 
        Else If (ChemistryDefn%FieldSpeciesNames(i) .CIEq. 'H2O2') Then 
        Else
          ConcSpeciesOnFields(i) = 0.0
        End If
      End Do

    Else

      ! Call the chemistry scheme routine (this performs all the gaseous and aqueous
      ! phase chemistry for a chemistry grid box). 

      !$$ Need to consider chemistry timestep more carefully - check that SyncDt is an
      !$$ exact integer and that it is exactly divisible by the chemistry timestep.

      ! Set O3 increment to zero, to ensure that any increment calculation will be repeated on
      ! the next occasion when there are no particles in the gridbox.

      ChemistryState%GridboxState(iX, iY, iZ)%DeltaO3 = 0.0

      ! Sort species for chemistry, do chemistry, and store results. 

      ! Note the chemistry species which are not either Lagrangian or Eulerian species will be passed to the
      ! chemistry with zero concentration and any chemical product will be discarded.

      Do i = 1, ChemistryDefn%nChemSpecieses
        ConcChemistrySpecies(i) = 0.0
        If (ChemistryDefn%iChem2Particle(i) /= 0) Then
          ConcChemistrySpecies(i) = ConcChemistrySpecies(i) +                               &
                                    ConcSpeciesOnParticles(ChemistryDefn%iChem2Particle(i))
        End If
        If (ChemistryDefn%iChem2Field   (i) /= 0) Then
          ConcChemistrySpecies(i) = ConcChemistrySpecies(i) +                               &
                                    ConcSpeciesOnFields   (ChemistryDefn%iChem2Field   (i))
        End If
      End Do

      Call ChemistryScheme(ConcChemistrySpecies,                                        &
                           NInt(SyncInterval), X(3), T, MolecularAirDensity,            &
                           Q, TotalOrDynCloudWater, TotalOrDynCloudIce,                 &
                           CloudFraction, TotalColumnCloudFraction, P, H, ZenithAngle)

      ! Assign mass back to particles or fields. For a species held on both particles and fields, the assignment is to fields.
      Do i = 1, ChemistryDefn%nChemSpecieses
        If (ChemistryDefn%iChem2Particle(i) /= 0 .and. ChemistryDefn%iChem2Field (i) /= 0) Then
          ConcSpeciesOnParticles(ChemistryDefn%iChem2Particle(i)) = 0.0
          ConcSpeciesOnFields   (ChemistryDefn%iChem2Field   (i)) = ConcChemistrySpecies(i)
        Else If (ChemistryDefn%iChem2Particle(i) /= 0) Then
          ConcSpeciesOnParticles(ChemistryDefn%iChem2Particle(i)) = ConcChemistrySpecies(i)
        Else If (ChemistryDefn%iChem2Field   (i) /= 0) Then
          ConcSpeciesOnFields   (ChemistryDefn%iChem2Field   (i)) = ConcChemistrySpecies(i)
        End If
      End Do

    End If

    ! Convert the updated gridbox concentrations (in molecules/cm3) of species
    ! held on chemistry fields back to g/m3.
    Call ConvertConcToChemField(ConcSpeciesOnFields(:), iX, iY, iZ,                     &
                                ChemOpts, ChemistryDefn, ChemistryState, EulerianField)

    ! Convert the updated gridbox concentrations (in molecules/cm3) of species
    ! held on particles back to total gridbox masses (in g).
    Call ConvertConcToMass(ConcSpeciesOnParticles(:), iX, iY, iZ, ChemOpts, &
                           ChemistryDefn, ChemistryState)

    ! Reassign the total gridbox masses back to the particles in the grid box.
    Call UpdateParticleMasses(                      &
           ChemOpts, ChemistryDefn, ChemistryState, &
           iX, iY, iZ,                              &
           Particles, Masses                        &
         )

  End Do
  End Do
  End Do

  !$OMP END PARALLEL DO
  Call TimerOff(ChemistryLoopTimer)
  Call TimerWriteLast(ChemistryLoopTimer)

End Subroutine UpdateChemistry

!-------------------------------------------------------------------------------------------------------------

End Module ChemistryModule
