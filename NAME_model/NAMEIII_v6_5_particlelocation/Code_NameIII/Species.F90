! Module:  Species Module

Module SpeciesModule

! This module provides code for handling different species of material.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use FlowsModule, Only: Flow_, Cloud_, Rain_, Surface_, Plant_
Use SizeDistModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: Species_                         ! Information defining one species.
Public  :: SpeciesUses_                     ! A set of species uses.
Public  :: Specieses_                       ! A collection of species.
Public  :: CloudGammaParams_                ! Information defining parameters for the calculation of
                                            ! cloud gamma dose.
Public  :: CloudGammaParamses_              ! A collection of photon energies and asscoiated information
                                            ! defining
                                            ! parameters for the calculation of cloud gamma dose.
Public  :: MassLimit_                       ! A particle mass limit.
Public  :: MassLimits_                      ! A collection of particle mass limits.
Public  :: InitSpecieses                    ! Initialises a collection of specieses.
Public  :: InitCloudGammaParamses           ! Initialises a collection of photon energies and associated
                                            ! cloud gamma parameters.                  ! $$ reorder list
Public  :: InitSpecies                      ! Initialises a species.
Public  :: InitCloudGammaParams             ! Initialises a photon energy and associated cloud gamma
                                            ! parameters.
Public  :: AddSpecies                       ! Adds a species to a collection of specieses.
Public  :: AddCloudGammaParams              ! Adds a photon energy and associated cloud gamma parameters to
                                            ! a collection of cloud
                                            ! gamma parameterses.
Public  :: FindSpeciesIndex                 ! Finds the index of a species.
Public  :: FindCloudGammaParamsIndex        ! Finds the index of a photon energy and associated cloud
                                            ! gamma parameters.
Public  :: InitSpeciesUses                  ! Initialises a set of species uses.
Public  :: AddSpeciesUses                   ! Adds a set of species uses to a collection of specieses.
Public  :: InitMassLimits                   ! Initialises a collection of particle mass limits.
Public  :: InitMassLimit                    ! Initialises a particle mass limit.
Public  :: AddMassLimit                     ! Adds a particle mass limit to a collection of particle mass
                                            ! limits.
Public  :: SetUpSpecieses                   ! Sets up Specieses.
Public  :: SetUpSpecieses_MassLimits        ! Sets up Specieses using information from MassLimits.
Public  :: SetUpSpecieses_DecayChains       ! Sets up full decay chains in Specieses.
Public  :: SetUpiSpecieses                  ! Sets up indices in Specieses.
Public  :: SetUpSpecieses_iCloudGammaParams ! Sets up indices in Species for referring to sets of cloud
                                            ! gamma parameters.
Public  :: SetUpSpecieses_iMaterialUnits    ! Sets up indices in Species for referring to material units.
Public  :: CalcDecayFactor                  ! Calculates decay factor for decay of species.
Public  :: CalcVd                           ! Calculates dry deposition velocity (including sedimentation
                                            ! contribution) for the species
Public  :: CalcWetScavCoeff                 ! Calculates wet deposition scavenging coefficient for species

Public :: Con
Integer, Parameter :: Con = 8 ! Precision of field concentration. 8 might be useful for chemistry, but perhaps remove if not useful $$
                              ! It's used in EulerianInterface.F90 and SemiLagrangian.F90

!-------------------------------------------------------------------------------------------------------------

Type :: Species_ ! A species.
  Character(MaxCharLength) :: Name
  Character(MaxCharLength) :: Category
  Character(MaxCharLength) :: MaterialUnitName
  Integer                  :: iMaterialUnit
  Real(Std)                :: MolecularWeight
  Real(Std)                :: InvHalfLife
  Logical                  :: HasDaughter
  Character(MaxCharLength) :: Daughter
  Real(Std)                :: BranchingRatio
  Integer                  :: DecayChain(MaxDecayChainLength)
  Integer                  :: DecayChainLength
  Character(MaxCharLength) :: CloudGammaParameters
  Integer                  :: iCloudGammaParams
  Real(Std)                :: UVLossRate
  Logical                  :: UsePowerLawDecay
  Real(Std)                :: PowerLawDecayExponent
  Real(Std)                :: PowerLawDecayDelay
  Logical                  :: UseVd
  Real(Std)                :: Rs
  Real(Std)                :: Vd
  Real(Std)                :: ABelowRain
  Real(Std)                :: BBelowRain
  Real(Std)                :: ABelowSnow
  Real(Std)                :: BBelowSnow
  Real(Std)                :: AInRain
  Real(Std)                :: BInRain
  Real(Std)                :: AInSnow
  Real(Std)                :: BInSnow
  Logical                  :: UseMassLimit
  Real(Std)                :: MassLimit
  Real(Std)                :: AerosolDiameter
  Logical                  :: LandUseDryDep
  ! Name                  :: Name of species.
  ! Category              :: Species category. This is a user defined quantity which is read
  !                          in with the species information and is written with the output
  !                          (e.g. the user may wish to clearly identify radioactive species
  !                          by giving them the category 'Radionuclide'). It has no other
  !                          significance. It can be blank.
  ! MaterialUnitName      :: Name of unit in which the quantity of material is measured.
  ! iMaterialUnit         :: Index of unit in which the quantity of material is measured.
  ! MolecularWeight       :: Molecular weight.
  ! InvHalfLife           :: Reciprocal of the half life.
  ! HasDaughter           :: Determines whether or not the radionuclide concerned decays to
  !                          a daughter product.
  ! Daughter              :: The name of the 1st daughter product resulting from the respective
  !                          parent.
  ! BranchingRatio        :: Fraction of the activity of the parent which is available for decay to
  !                          the daughter product.
  ! DecayChain            :: Details the decay chains (comprising of only the radionuclides
  !                          considered).
  ! DecayChainLength      :: Details the number of radionuclides, including the parent, in each
  !                          of the considered decay chains.
  ! CloudGammaParameters  :: Cloud gamma parameters of species.
  ! iCloudGammaParams     :: Index of cloud gamma parameters in collection of all
  !                          cloud gamma parameters.
  ! UVLossRate            :: Loss rate (per hour) from UV decay.
  ! UsePowerLawDecay      :: Indicates power law decay to be used.
  ! PowerLawDecayExponent :: Power law decay exponent.
  ! PowerLawDecayDelay    :: Travel time before power law decay starts.
  ! UseVd                 :: Indicates Vd rather than Rs is to be used.
  ! Rs                    :: Surface resistance. Defined only if UseVd is false.
  ! Vd                    :: Deposition velocity. Defined only if UseVd is true.
  ! ABelowRain            ::] Wet deposition coefficients. Scavenging coefficient is given by lambda = A * r^B
  ! BBelowRain            ::] where r is the rainfall rate (mm/hr). A and B take different values depending on
  ! ABelowSnow            ::] whether the particle is in cloud or below cloud and whether the precipitation
  ! BBelowSnow            ::] is snow (T<270K at particle location) or rain (T>270K at particle location).
  ! AInRain               ::] If specified, all A, B values should be >= 0. Leave all values blank for no
  ! BInRain               ::] wet deposition. Alternatively, setting any A value to zero will deactivate that
  ! AInSnow               ::] particular component of the wet deposition scheme (e.g. below-cloud washout
  ! BInSnow               ::] by rain is switched off by setting ABelowRain to zero).
  ! UseMassLimit          :: Indicates there is an upper limit on particle  } Note no upper
  !                          mass for this species.                         } limit applies
  ! MassLimit             :: Upper limit on particle mass for this species. } for fixed met.
  ! AerosolDiameter       :: Aerosol mean diameter (um) for dry deposition parameterisation
  ! LandUseDryDep         :: Indicates land use dependent dry deposition scheme to be used

  ! Note not all information on the species is given here. E.g. for complex chemicals
  ! or biological agents, various properties are hard wired in the code in association
  ! with the species name.

End Type Species_

!-------------------------------------------------------------------------------------------------------------

Type :: SpeciesUses_ ! A set of species uses. This gives information on the uses for a single species.
  Character(MaxCharLength) :: SpeciesName  ! Name of species.
  Logical                  :: OnParticles  ! Indicates species is carried on particles.
  Logical                  :: OnFields     ! Indicates species is carried on fields.
  Logical                  :: AdvectField  ! Indicates the species is advected on fields.
  Character(MaxCharLength) :: SizeDistName ! Name of particle size distribution. Note the distribution itself
                                           ! is not used but the bins associated with the distribution are
                                           ! used to bin the material when carried on fields.
  Integer                  :: iSizeDist    ! Index of particle size distribution in collection of all particle
                                           ! size distributions.
End Type SpeciesUses_

!-------------------------------------------------------------------------------------------------------------

Type :: Specieses_ ! A collection of species.
  Integer            :: nSpecieses                                                  ! Number of species.
  Type(Species_)     :: Specieses(MaxSpecieses)                                     ! Collection of species.
  Integer            :: nSpeciesUseses                                              ! Number of sets of species uses.
  Type(SpeciesUses_) :: SpeciesUseses(MaxSpecieses)                                 ! Collection of sets of species uses.
  Integer            :: nParticleSpecieses                                          ! Number of particle species.
  Integer            :: nFields                                                     ! Number of fields.
  Integer            :: iSpecies2SpeciesUses     (MaxSpecieses                    ) !} Cross referencing index
  Integer            :: iSpeciesUses2Species     (MaxSpecieses                    ) !} arrays - see below.
  Integer            :: iSpecies2Particle        (MaxSpecieses                    ) !}
  Integer            :: iParticle2Species        (MaxSpecieses                    ) !}
  Integer            :: iSpeciesUses2Particle    (MaxSpecieses                    ) !}
  Integer            :: iParticle2SpeciesUses    (MaxSpecieses                    ) !}
  Integer            :: iSpeciesAndSize2Field    (MaxSpecieses,  MaxDiameterRanges) !}
  Integer            :: iSpeciesUsesAndSize2Field(MaxSpecieses,  MaxDiameterRanges) !}
  Integer            :: iField2Species           (MaxSpecieses * MaxDiameterRanges) !}
  Integer            :: iField2SpeciesUses       (MaxSpecieses * MaxDiameterRanges) !}
  Integer            :: iField2Size              (MaxSpecieses * MaxDiameterRanges) !}
  Integer            :: iSpecies2SizeDist        (MaxSpecieses                    ) !}
  Integer            :: iSpeciesUses2SizeDist    (MaxSpecieses                    ) !}
  Integer            :: iField2SizeDist          (MaxSpecieses * MaxDiameterRanges) !}

  Integer            :: iFirstAdvectedNonSedimentingField !} used for eulerian evolution. 
  Integer            :: iLastAdvectedNonSedimentingField  !} Might be better stored there?? $$
  Logical            :: AdvectedFields                    !}
  Logical            :: AdvectedSedimentingFields         !}

  ! $$ currently assumes all nonsedimenting advected species are together in Field. See SetUpiSpecieses too.
  !    Could remove this requirement, even if actually coding this way for efficiency.
  
  ! The cross referencing index arrays give the index of an object in one set corresponding to a given object 
  ! in another set. For example iParticle2Species(5) = 7 indicates that the 5th species held on particles is 
  ! the 7th species in the collection of species.
  !
  ! Here is a diagram of indices relating species, species uses, species held on particles, particulate size 
  ! ranges, size distributions, and fields:
  !
  !                                                Field
  !                                               /  |  \
  !                                              /   |   \
  !                                             /    |    \
  !                                            /     |     \
  !                                           /      |      \
  !                                          /       |       \
  !                                         /    SizeDist     \
  !                                        /    /        \     \
  !                                       /    /          \     \
  !                                      /    /            \     \
  !                                     /    /              \     \
  !                                    /    /                \     \
  !                                   /    /                  \     \
  !                        SizeRange & Species --------- SpeciesUses & SizeRange
  !                                        \                  /
  !                                         \                /
  !                                          \              /
  !                                           \            /
  !                                            \          /
  !                                             \        /
  !                                              \      /  
  !                                           ParticleSpecies
  !
  ! Links to SizeDist are all one way. Other links are all two way. Values given by the index arrays are zero
  ! if the target object doesn't exist. However, for the domain of iSpeciesAndSize2Field and of 
  ! iSpeciesUsesAndSize2Field, the particulate size range index is 1 if there are no particulate size ranges 
  ! (and no particulate size distribution).

End Type Specieses_

!-------------------------------------------------------------------------------------------------------------

Type :: CloudGammaParams_ ! Defines all parameters, primarily photon energy dependent
                          ! parameters, required for the cloud gamma dose calculations
  Character(MaxCharLength) :: Name
  Integer                  :: nEnergies
  Real(Std)                :: PhotonEnergy(MaxEnergies)
  Real(Std)                :: PhotonIntensity(MaxEnergies)
  Real(Std)                :: LinearAttenCoeff(MaxEnergies)
  Real(Std)                :: BBuildUpFactora(MaxEnergies)
  Real(Std)                :: BBuildUpFactorb(MaxEnergies)
  Real(Std)                :: AirKermapuFluence(MaxEnergies)
  Real(Std)                :: AdEffDosepuAirKerma(MaxEnergies)
  Real(Std)                :: AdThyDosepuAirKerma(MaxEnergies)
  Real(Std)                :: AdLunDosepuAirKerma(MaxEnergies)
  Real(Std)                :: AdBoSDosepuAirKerma(MaxEnergies)
  ! Name                :: Name of the cloud gamma parameters block.
  ! nEnergies           :: Number of photon energies associated with a single radionculide.
  ! PhotonEnergy        :: The energy of a photon (or gamma-ray emission) from a species.
  !                        Note that a single species may emit multiple photons with varying energy's.
  !                        Note that photon energy is not used directly in the cloud gamma dose
  !                        calculations but is included to make it clear of the photon energy from which
  !                        all energy dependent parameters have been derived. Units = MeV
  ! PhotonIntensity     :: The frequency with which the photons are released. Note the frequency of the
  !                        photon emissions does not necessarily add up to 1. Only the frequnecy of all
  !                        emissions will add up to 1. No units.
  ! LinearAttenCoeff    :: The interaction probability per unit differential path length of a photon
  !                        i.e. how the radiation interacts with the matter through which it is passing.
  !                        Units = m^-1.
  ! BBuildUpFactor      :: Berger build-up factor. A factor describing the contribution to the photon flux at
  !                        the receptor from the scattered photons. Derived using coefficients: a and b.
  !                        No units.
  ! AirKermapuFluence   :: Air kerma per unit fluence is effectively the dose in air incident on a unit
  !                        area. Units = Gy m^2.
  ! AdEffDosepuAirKerma :: Adult effective dose per unit air kerma is the whole body dose to an adult
  !                        per unit dose in air. Units = Sv/Gy.
  ! AdThyDosepuAirKerma :: Adult thyroid dose per unit air kerma. Units = Gy/Gy.
  ! AdLunDosepuAirKerma :: Adult lung dose per unit air kerma. Units = Gy/Gy.
  ! AdBoSDosepuAirKerma :: Adult bone surface dose per unit air kerma. Units = Gy/Gy.

End Type CloudGammaParams_

!-------------------------------------------------------------------------------------------------------------

Type :: CloudGammaParamses_ ! A collection of sets of cloud gamma parameters.
  Integer                 :: nCloudGammaParamses               ! Number of sets of cloud gamma parameters.
  Type(CloudGammaParams_) :: CloudGammaParamses(MaxCGSpecies)  ! Sets of cloud gamma parameters.
End Type CloudGammaParamses_

!-------------------------------------------------------------------------------------------------------------

Type :: MassLimit_ ! A particle mass limit.
  Character(MaxCharLength) :: SpeciesName
  Real(Std)                :: Limit
  ! SpeciesName :: Name of species.
  ! Limit       :: Upper limit on particle mass for this species. Note no upper limit
  !                applies for fixed met.
End Type MassLimit_

!-------------------------------------------------------------------------------------------------------------

Type :: MassLimits_ ! A collection of particle mass limits.
  Private
  Integer          :: nMassLimits              ! Number of particle mass limits.
  Type(MassLimit_) :: MassLimits(MaxSpecieses) ! Collection of particle mass limits.
End Type MassLimits_

!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of specieses, sets of species uses and cloud gamma parameters.
  Module Procedure SpeciesEq
  Module Procedure SpeciesUsesEq
  Module Procedure CloudGammaParamsEq
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitSpecieses() Result (Specieses)
! Initialises a collection of specieses.

  Implicit None
  ! Function result:
  Type(Specieses_) :: Specieses ! Initialised collection of specieses.

  Specieses%nSpecieses     = 0
  Specieses%nSpeciesUseses = 0

End Function InitSpecieses

!-------------------------------------------------------------------------------------------------------------

Function InitSpecies(                                                   &
           Name, Category, MaterialUnitName,                            &
           MaterialUnits,                                               &
           MolecularWeight,                                             &
           InvHalfLife, HasDaughter,                                    &
           Daughter, BranchingRatio,                                    &
           CloudGammaParams,                                            &
           UVLossRate,                                                  &
           UsePowerLawDecay, PowerLawDecayExponent, PowerLawDecayDelay, &
           UseVd, Rs,  Vd,                                              &
           ABelowRain, BBelowRain,                                      &
           ABelowSnow, BBelowSnow, AInRain,                             &
           BInRain, AInSnow, BInSnow,                                   &
           LandUseDryDep, AerosolDiameter                               &
         )                                                              &
Result (Species)
! Initialises a species.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name                     ! Name of species.
  Character(*), Intent(In) :: Category                 ! Species category.
  Character(*), Intent(In) :: MaterialUnitName         ! Unit in which the quantity of material is measured
  Type(MaterialUnits_), Intent(InOut) :: MaterialUnits ! Collection of material units
  Real(Std),    Intent(In) :: MolecularWeight          ! Molecular weight.
  Real(Std),    Intent(In) :: InvHalfLife              ! Reciprocal of the half life.
  Logical,      Intent(In) :: HasDaughter              ! Determines whether or not the radionuclide
                                                       ! concerned decays to a daughter product.
  Character(*), Intent(In) :: Daughter                 ! The name of the 1st daughter product
                                                       ! resulting from the respective parent.
  Real(Std),    Intent(In) :: BranchingRatio           ! Fraction of the activity of the parent which
                                                       ! is available for decay to the daughter product.
  Character(*), Intent(In) :: CloudGammaParams         ! Name of set of cloud gamma parameters.
  Real(Std),    Intent(In) :: UVLossRate               ! Loss rate (per hour) from UV decay.
  Logical,      Intent(In) :: UsePowerLawDecay         ! Indicates power law decay to be used.
  Real(Std),    Intent(In) :: PowerLawDecayExponent    ! Power law decay exponent.
  Real(Std),    Intent(In) :: PowerLawDecayDelay       ! Travel time before power law decay starts.
  Logical,      Intent(In) :: UseVd                    ! Indicates Vd rather than Rs is to be
                                                       ! used.
  Real(Std),    Intent(In) :: Rs                       ! Surface resistance. Defined only if
                                                       ! UseVd is false.
  Real(Std),    Intent(In) :: Vd                       ! Deposition velocity. Defined only if
                                                       ! UseVd is true.
  Real(Std),    Intent(In) :: AerosolDiameter          ! Aerosol mean diameter (um) for dry deposition 
                                                       ! parametrisation
  Real(Std),    Intent(In) :: ABelowRain               !] Wet deposition coefficients
  Real(Std),    Intent(In) :: BBelowRain               !]
  Real(Std),    Intent(In) :: ABelowSnow               !] 
  Real(Std),    Intent(In) :: BBelowSnow               !]
  Real(Std),    Intent(In) :: AInRain                  !]
  Real(Std),    Intent(In) :: BInRain                  !]
  Real(Std),    Intent(In) :: AInSnow                  !]
  Real(Std),    Intent(In) :: BInSnow                  !]
  Logical,      Intent(In) :: LandUseDryDep            !} Indicates land use dependent 
                                                       !} dry deposition scheme is to be used
  Integer                  :: SpeciesMaterialUnitIndex ! Material unit of species
  Integer                  :: SpeciesMaterialUnitType  ! Type of material unit of species
  ! Function result:
  Type (Species_) :: Species ! Initialised species.

  If (Name == ' ') Then
    Call Message(                                                  &
           'FATAL ERROR in InitSpecies: name of species is blank', &
           3                                                       &
         )
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message(                             &
           'FATAL ERROR in InitSpecies: '  // &
           'name of species is given as "' // &
           Trim(Name)                      // &
           '" and is too long',               &
           3                                  &
         )
  End If
  Species%Name = Name

  If (Len_Trim(Category) > MaxCharLength) Then
    Call Message(                              &
           'FATAL ERROR in InitSpecies: '   // &
           'species category is given as "' // &
           Trim(Category)                   // &
           '" and is too long',                &
           3                                   &
         )
  End If
  Species%Category = Category
 
  If (Len_Trim(MaterialUnitName) > MaxCharLength) Then
    Call Message(                                      &
           'FATAL ERROR in InitSpecies: '           // &
           'material unit of species is given as "' // &
           Trim(MaterialUnitName)                   // &
           '" and is too long',                        &
           3                                           &
         )
  End If
  If (FindMaterialUnitIndex(MaterialUnitName, MaterialUnits) == -1) Then
    Call Message(                                                       &
           'Material unit "' // Trim(MaterialUnitName)  //              &
           '" in species "' // Trim(Species%Name)       //              &
           '" not recognized. Adding to allowed units.',                &
           1                                                            &
         )
    Call AddMaterialUnit(MaterialUnitName,'unknown unit',UnknownUnitType,1.0,MaterialUnits)
  End If
  SpeciesMaterialUnitIndex = FindMaterialUnitIndex(MaterialUnitName, MaterialUnits)
  SpeciesMaterialUnitType  = MaterialUnits%MaterialUnits(SpeciesMaterialUnitIndex)%UnitType
  If ( ( SpeciesMaterialUnitType .ne. UnknownUnitType  ) .and. &
       ( SpeciesMaterialUnitType .ne. MassUnitType     ) .and. &
       ( SpeciesMaterialUnitType .ne. ActivityUnitType )       &
     ) Then
    Call Message(                                                 &
           'ERROR: Material unit "' // Trim(MaterialUnitName)  // &
           '" in species "' // Trim(Species%Name)              // &
           '" not allowed.',                                      &
           3                                                      &
         )
  End If
  Species%MaterialUnitName = MaterialUnitName

  !$$ Should we check that mol. weight etc. are all positive here?? In fact these
  ! range checking issues apply throughout the code.
  Species%MolecularWeight = MolecularWeight
  Species%InvHalfLife     = InvHalfLife

  If (UVLossRate < 0.0 .or. UVLossRate > 100.0) Then
    Call Message(                                                 &
           'FATAL ERROR in InitSpecies: '                      // &
           'UV loss rate is given as '                         // &
           Std2Char(UVLossRate)                                // &
           ' and is outside the permitted range [0% -> 100%]',    &
           3                                                      &
         )
  End If
  Species%UVLossRate = UVLossRate

  Species%UsePowerLawDecay = UsePowerLawDecay
  If (UsePowerLawDecay) Then
    If (PowerLawDecayExponent <= 0.0) Then
      Call Message(                                     &
             'FATAL ERROR in InitSpecies: '          // &
             'Power Law Decay Exponent is given as ' // &
             Std2Char(PowerLawDecayExponent)         // &
             ' and is <= 0.',                           &
             3                                          &
           )
    End If
    If (PowerLawDecayDelay <= 0.0) Then
      Call Message(                                  &
             'FATAL ERROR in InitSpecies: '       // &
             'Power Law Decay Delay is given as ' // &
             Std2Char(PowerLawDecayDelay)         // &
             ' and is <= 0.',                        &
             3                                       &
           )
    End If
  End If
  Species%PowerLawDecayExponent = PowerLawDecayExponent
  Species%PowerLawDecayDelay    = PowerLawDecayDelay

  ! $$ probably should exclude power law decay and radioactive decay. 
  ! Not quite sure what decay options can be used together.

  Species%UseVd = UseVd

  Species%LandUseDryDep = LandUseDryDep

  Species%AerosolDiameter = AerosolDiameter


  If (.not. Species%LandUseDryDep .and. Species%AerosolDiameter == 0.0) Then
    If (UseVd) Then
      Species%Vd = Vd
    Else
      Species%Rs = Rs
    End If
  End If

  Species%ABelowRain = ABelowRain
  Species%BBelowRain = BBelowRain
  Species%ABelowSnow = ABelowSnow
  Species%BBelowSnow = BBelowSnow
  Species%AInRain    = AInRain
  Species%BInRain    = BInRain
  Species%AInSnow    = AInSnow
  Species%BInSnow    = BInSnow

  Species%HasDaughter    = HasDaughter
  Species%Daughter       = ''
  Species%BranchingRatio = -1.0
  If (HasDaughter) Then
    If (Len_Trim(Daughter) > MaxCharLength) Then
      Call Message(                                      &
             'FATAL ERROR in InitSpecies: '           // &
             'name of daughter species is given as "' // &
             Trim(Daughter)                           // &
             '" and is too long',                        &
             3                                           &
           )
    End If
    Species%Daughter = Daughter
    If (BranchingRatio < 0.0 .or. BranchingRatio > 1.0) Then
      Call Message(                                                 &
             'FATAL ERROR in InitSpecies: '                      // &
             'Branching Ratio is given as '                      // &
             Std2Char(BranchingRatio)                            // &
             ' and is outside the permitted range [0.0 -> 1.0]',    &
             3                                                      &
           )
    End If
    Species%BranchingRatio = BranchingRatio
  End If

  If (Len_Trim(CloudGammaParams) > MaxCharLength) Then
    Call Message(                                           &
           'FATAL ERROR in InitSpecies: '                // &
           'set of cloud gamma parameters is given as "' // &
           Trim(CloudGammaParams)                        // &
           '" and is too long',                             &
           3                                                &
         )
  End If
  Species%CloudGammaParameters = CloudGammaParams

End Function InitSpecies

!-------------------------------------------------------------------------------------------------------------

Subroutine AddSpecies(Species, Specieses)
! Adds a species to a collection of specieses.

  Implicit None
  ! Argument list:
  Type(Species_),   Intent(In)    :: Species   ! Species to be added.
  Type(Specieses_), Intent(InOut) :: Specieses ! Collection of specieses.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Specieses%nSpecieses
    If (Species%Name .CIEq. Specieses%Specieses(i)%Name) Then
      If (Species == Specieses%Specieses(i)) Then
        Return
      Else
        Call Message(                                                         &
               'FATAL ERROR in adding the species "'                       // &
               Trim(Species%Name)                                          // &
               '": a different species with the same name already exists',    &
               3                                                              &
             )
      End If
    Endif
  End Do

  If (Specieses%nSpecieses >= MaxSpecieses) Then
    Call Message(                                   &
           'FATAL ERROR in adding the species "' // &
           Trim(Species%Name)                    // &
           '": there are too many species',         &
           3                                        &
         )
  End If

  Specieses%nSpecieses                      = Specieses%nSpecieses + 1
  Specieses%Specieses(Specieses%nSpecieses) = Species

End Subroutine AddSpecies

!-------------------------------------------------------------------------------------------------------------

Function FindSpeciesIndex(Name, Specieses, Error)
! Finds the index of a species.

  Implicit None
  ! Argument list:
  Character(*),     Intent(In)            :: Name      ! Name of species.
  Type(Specieses_), Intent(In)            :: Specieses ! Collection of specieses.
  Logical,          Intent(Out), Optional :: Error     ! Error flag for species not found.
  ! Function result:
  Integer :: FindSpeciesIndex ! Index of the species.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Specieses%nSpecieses
    If (Name .CIEq. Specieses%Specieses(i)%Name) Then
      FindSpeciesIndex = i
      If (Present(Error)) Error = .false.
      Return
    End If
  End Do

  If (Present(Error)) Then
    FindSpeciesIndex = 0
    Error = .true.
  Else
    Call Message(                      &
           'FATAL ERROR: species "' // &
           Trim(Name)               // &
           '" not found',              &
           3                           &
         )
  End If

End Function FindSpeciesIndex

!-------------------------------------------------------------------------------------------------------------

Function SpeciesEq(Species1, Species2)
! Tests for equality of species.

  Implicit None
  ! Argument list:
  Type(Species_), Intent(In) :: Species1 !} The two species.
  Type(Species_), Intent(In) :: Species2 !}
  ! Function result:
  Logical :: SpeciesEq ! Indicates if species are equal.

  SpeciesEq = (Species1%Name                    .CIEq. Species2%Name)                   .and. &
              (Species1%Category                .CIEq. Species2%Category)               .and. &
              (Trim(Species1%MaterialUnitName)    ==   Trim(Species2%MaterialUnitName)) .and. & !$$ trim needed?
               Species1%MolecularWeight           ==   Species2%MolecularWeight         .and. &
               Species1%InvHalfLife               ==   Species2%InvHalfLife             .and. &
              (Species1%Daughter                .CIEq. Species2%Daughter)               .and. &
               Species1%BranchingRatio            ==   Species2%BranchingRatio          .and. &
              (Species1%CloudGammaParameters    .CIEq. Species2%CloudGammaParameters)   .and. &
               Species1%UVLossRate                ==   Species2%UVLossRate              .and. &
               Species1%UseVd                   .eqv.  Species2%UseVd                   .and. &
               Species1%ABelowRain                ==   Species2%ABelowRain              .and. &
               Species1%BBelowRain                ==   Species2%BBelowRain              .and. &
               Species1%ABelowSnow                ==   Species2%ABelowSnow              .and. &
               Species1%BBelowSnow                ==   Species2%BBelowSnow              .and. &
               Species1%AInRain                   ==   Species2%AInRain                 .and. &
               Species1%BInRain                   ==   Species2%BInRain                 .and. &
               Species1%AInSnow                   ==   Species2%AInSnow                 .and. &
               Species1%BInSnow                   ==   Species2%BInSnow                 .and. &
               Species1%LandUseDryDep           .eqv.  Species2%LandUseDryDep           .and. &
               Species1%AerosolDiameter           ==   Species2%AerosolDiameter

  ! $$ Mass limit? HasDaughter? Move derived quantities (decay chain) to end of structure?

  If (SpeciesEq .and. (.not. Species1%LandUseDryDep) .and. (Species1%AerosolDiameter == 0.0)) Then
    If (Species1%UseVd) Then
      SpeciesEq = SpeciesEq .and. Species1%Vd == Species2%Vd
    Else
      SpeciesEq = SpeciesEq .and. Species1%Rs == Species2%Rs
    End If
  End If

End Function SpeciesEq

!-------------------------------------------------------------------------------------------------------------

Function InitSpeciesUses(                                                &
           SpeciesName, OnParticles, OnFields, AdvectField, SizeDistName &
         )                                                               &
Result (SpeciesUses)
! Initialises a set of species uses.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: SpeciesName  ! Name of species.
  Logical,      Intent(In) :: OnParticles  ! Indicates species is carried on particles.
  Logical,      Intent(In) :: OnFields     ! Indicates species is carried on fields. 
  Logical,      Intent(In) :: AdvectField  ! Indicates species is advected. 
  Character(*), Intent(In) :: SizeDistName ! Name of particle size distribution. Note the distribution itself
                                           ! is not used but the bins associated with the  distribution are
                                           ! used to bin the material when carried on fields.
  ! Function result:
  Type (SpeciesUses_) :: SpeciesUses ! Initialised set of species uses.  
  ! Locals:
  Integer :: iParticleSize ! Index.
  
  Call TokenLengthTest(              &
         C         = SpeciesName,    & 
         Length    = MaxCharLength,  &
         Zero      = .true.,         &
         BlockKey  = 'Species Uses', &
         Item      = ' ',            &
         ColumnKey = 'Species'       &
       )
  Call TokenLengthTest(                           &
         C         = SizeDistName,                & 
         Length    = MaxCharLength,               &
         Zero      = .false.,                     &
         BlockKey  = 'Species Uses',              &
         Item      = ' ',                         &
         ColumnKey = 'Particle Size Distribution' &
       )
  
  If ((.not. OnFields) .and. SizeDistName /= ' ') Then
    Call Message(                                                                                    &
           'FATAL ERROR in inputting "Species Uses" information for species "'                    // & 
           Trim(SpeciesName)                                                                      // &
           '". "Particle Size Distribution" should not be given for species not held on fields.',    &
           3                                                                                         &
         )
  End If

  ! $$ currently non advected fields have no dep, sedimentation etc (but could have separate flags).
  !    Hence don't currently allow size-dist or on-particles.
  
  If (OnParticles .and. (OnFields .and. .not. AdvectField)) Then
    Call Message(                                                                  &
           'FATAL ERROR in inputting "Species Uses" information for species "'  // & 
           Trim(SpeciesName)                                                    // &
           '". Species held on static fields should not be held on particles.',    &
           3                                                                       &
         )
  End If

  If ((.not. (OnFields .and. AdvectField)) .and. SizeDistName /= ' ') Then
    Call Message(                                                                                               &
           'FATAL ERROR in inputting "Species Uses" information for species "'                               // & 
           Trim(SpeciesName)                                                                                 // &
           '". "Particle Size Distribution" should not be given for species not held on (advected) fields.',    &
           3                                                                                                    &
         )
  End If
  
  SpeciesUses%SpeciesName  = SpeciesName
  SpeciesUses%OnParticles  = OnParticles
  SpeciesUses%OnFields     = OnFields
  SpeciesUses%AdvectField  = AdvectField
  SpeciesUses%SizeDistName = SizeDistName

End Function InitSpeciesUses

!-------------------------------------------------------------------------------------------------------------

Subroutine AddSpeciesUses(SpeciesUses, Specieses)
! Adds a set of species uses to a collection of species.

  Implicit None
  ! Argument list:
  Type(SpeciesUses_), Intent(In)    :: SpeciesUses ! Set of species uses to be added.
  Type(Specieses_),   Intent(InOut) :: Specieses   ! Collection of species.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Specieses%nSpeciesUseses
    If (SpeciesUses%SpeciesName .CIEq. Specieses%SpeciesUseses(i)%SpeciesName) Then
      If (SpeciesUses == Specieses%SpeciesUseses(i)) Then
        Return
      Else
        Call Message(                                                                        &
               'FATAL ERROR in adding a set of species uses for species "'                // &
               Trim(SpeciesUses%SpeciesName)                                              // &
               '": a different set of species uses for the same species already exists.',    &
               3                                                                             &
             )
      End If
    Endif
  End Do

  If (Specieses%nSpeciesUseses >= Size(Specieses%SpeciesUseses)) Then
    Call Message(                                                         &
           'FATAL ERROR in adding a set of species uses for species "' // &
           Trim(SpeciesUses%SpeciesName)                               // &
           '": there are too many sets of species uses.',                 &
           3                                                              &
         )
  End If

  Specieses%nSpeciesUseses                          = Specieses%nSpeciesUseses + 1
  Specieses%SpeciesUseses(Specieses%nSpeciesUseses) = SpeciesUses

End Subroutine AddSpeciesUses

!-------------------------------------------------------------------------------------------------------------

Function FindSpeciesUsesIndex(Name, Specieses, Error)
! Finds the index of a set of species uses.

  Implicit None
  ! Argument list:
  Character(*),     Intent(In)            :: Name      ! Name of set of species uses (named after the species 
                                                       ! it corresponds to).
  Type(Specieses_), Intent(In)            :: Specieses ! Collection of specieses.
  Logical,          Intent(Out), Optional :: Error     ! Error flag for set of species uses not found.
  ! Function result:
  Integer :: FindSpeciesUsesIndex ! Index of the set of species uses.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Specieses%nSpeciesUseses
    If (Name .CIEq. Specieses%SpeciesUseses(i)%SpeciesName) Then
      FindSpeciesUsesIndex = i
      If (Present(Error)) Error = .false.
      Return
    End If
  End Do

  If (Present(Error)) Then
    FindSpeciesUsesIndex = 0
    Error = .true.
  Else
    Call Message(                                        &
           'FATAL ERROR: species "'                   // &
           Trim(Name)                                 // &
           '" not found in species uses input block',    &
           3                                             &
         )
  End If

End Function FindSpeciesUsesIndex

!-------------------------------------------------------------------------------------------------------------

Function SpeciesUsesEq(SpeciesUses1, SpeciesUses2)
! Tests for equality of sets of species uses.

  Implicit None
  ! Argument list:
  Type(SpeciesUses_), Intent(In) :: SpeciesUses1 !} The two sets of species uses.
  Type(SpeciesUses_), Intent(In) :: SpeciesUses2 !}
  ! Function result:
  Logical :: SpeciesUsesEq ! Indicates if sets of species uses are equal.

  SpeciesUsesEq = (SpeciesUses1%SpeciesName  .CIEq. SpeciesUses2%SpeciesName ) .and. &
                  (SpeciesUses1%OnParticles  .eqv.  SpeciesUses2%OnParticles ) .and. &
                  (SpeciesUses1%OnFields     .eqv.  SpeciesUses2%OnFields    ) .and. &
                  (SpeciesUses1%AdvectField  .eqv.  SpeciesUses2%AdvectField ) .and. &
                  (SpeciesUses1%SizeDistName .CIEq. SpeciesUses2%SizeDistName)

End Function SpeciesUsesEq

!-------------------------------------------------------------------------------------------------------------

Function InitCloudGammaParamses() Result (CloudGammaParamses)
! Initialises a collection of CloudGammaParamses.

  Implicit None
  ! Function result:
  Type(CloudGammaParamses_) :: CloudGammaParamses ! Initialised collection of CloudGammaParamses.

  CloudGammaParamses%nCloudGammaParamses = 0

End Function InitCloudGammaParamses

!-------------------------------------------------------------------------------------------------------------

Function InitCloudGammaParams(                        &
           Name, nEnergies,                           &
           PhotonEnergy, PhotonIntensity,             &
           LinearAttenCoeff, BBuildUpFactora,         &
           BBuildUpFactorb,                           &
           AirKermapuFluence, AdEffDosepuAirKerma,    &
           AdThyDosepuAirKerma, AdLunDosepuAirKerma,  &
           AdBoSDosepuAirKerma                        &
         )                                            &
Result (CloudGammaParams)
! Initialises a the cloud gamma parameters for a single photon energy.

  Implicit None
  ! Argument list:
  Character(MaxCharLength), Intent(In) :: Name         ! Name of cloud gamma parameters block.
  Integer,      Intent(In) :: nEnergies                ! } Number of photon energies associated with each
                                                       ! } radionuclide.
  Real(Std),    Intent(In) :: PhotonEnergy(:)          ! Energy of photon.
  Real(Std),    Intent(In) :: PhotonIntensity(:)       ! The frequency with which the photons are released.
  Real(Std),    Intent(In) :: LinearAttenCoeff(:)      ! Coefficient describing how the photon interacts with
                                                       ! the matter through which it is passing.
  Real(Std),    Intent(In) :: BBuildUpFactora(:)       ! } Berger build-up factor. A factor describing the
                                                       ! } contribution to the photon flux at
  Real(Std),    Intent(In) :: BBuildUpFactorb(:)       ! } the receptor from the scattered photons. Derived
                                                       ! } using coefficients: a and b. No units.
  Real(Std),    Intent(In) :: AirKermapuFluence(:)     ! Dose in air incident on a unit area.
  Real(Std),    Intent(In) :: AdEffDosepuAirKerma(:)   ! Whole body dose to an adult per unit dose in air
  Real(Std),    Intent(In) :: AdThyDosepuAirKerma(:)   ! Adult thyroid dose per unit air kerma.
  Real(Std),    Intent(In) :: AdLunDosepuAirKerma(:)   ! Adult lung dose per unit air kerma.
  Real(Std),    Intent(In) :: AdBoSDosepuAirKerma(:)   ! Adult bone surface dose per unit air kerma.

  ! Function result:
  Type (CloudGammaParams_) :: CloudGammaParams ! Initialised cloud gamma parameters.

  ! Locals:
  Integer :: i  ! Variable for looping over nEnergies.

  If (Name == ' ') Then
    Call Message(                                                                &
           'FATAL ERROR in InitCloudGammaParams: name of radionuclide is blank', &
           3                                                                     &
         )
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message(                                      &
           'FATAL ERROR in InitCloudGammaParams: '  // &
           'name of radionuclide is given as "'     // &
           Trim(Name)                               // &
           '" and is too long',                        &
           3                                           &
         )
  End If
  CloudGammaParams%Name = Name

  If (nEnergies > MaxEnergies) Then
    Call Message(                                                                                 &
           'FATAL ERROR in InitCloudGammaParams: too many energies in cloud gamma parameters ' // &
           'for radionuclide "'                                                                // &
           Trim(Name)                                                                          // &
           '"',                                                                                   &
           3                                                                                      &
         )
  End If
  CloudGammaParams%nEnergies = nEnergies

  Do i = 1, nEnergies

    If (PhotonEnergy(i) <= 0.0) Then
      Call Message(                                                                           &
             'FATAL ERROR in InitCloudGammaParams: photon energy of species is not positive', &
             3                                                                                &
           )
    Else If (PhotonEnergy(i) < 0.01 .or. PhotonEnergy(i) > 10.0) Then
      Call Message(                                                      &
             'WARNING in InitCloudGammaParams: '                      // &
             'Photon Energy is given as '                             // &
             Std2Char(PhotonEnergy(i))                                // &
             ' which is outside the expected range [0.01 -> 10] MeV',    &
             1                                                           &
           )
    End If

    CloudGammaParams%PhotonEnergy(i)        = PhotonEnergy(i)
    CloudGammaParams%PhotonIntensity(i)     = PhotonIntensity(i)
    CloudGammaParams%LinearAttenCoeff(i)    = LinearAttenCoeff(i)
    CloudGammaParams%BBuildUpFactora(i)     = BBuildUpFactora(i)
    CloudGammaParams%BBuildUpFactorb(i)     = BBuildUpFactorb(i)
    CloudGammaParams%AirKermapuFluence(i)   = AirKermapuFluence(i)
    CloudGammaParams%AdEffDosepuAirKerma(i) = AdEffDosepuAirKerma(i)
    CloudGammaParams%AdThyDosepuAirKerma(i) = AdThyDosepuAirKerma(i)
    CloudGammaParams%AdLunDosepuAirKerma(i) = AdLunDosepuAirKerma(i)
    CloudGammaParams%AdBoSDosepuAirKerma(i) = AdBoSDosepuAirKerma(i)

  End Do

End Function InitCloudGammaParams

!-------------------------------------------------------------------------------------------------------------

Subroutine AddCloudGammaParams(CloudGammaParams, CloudGammaParamses)
! Adds a set of cloud gamma parameters to a collection of sets of cloud gamma parameters.

  Implicit None
  ! Argument list:
  Type(CloudGammaParams_),   Intent(In)    :: CloudGammaParams   ! Set of cloud gamma parameters to be added.
  Type(CloudGammaParamses_), Intent(InOut) :: CloudGammaParamses ! Collection of sets of cloud gamma
                                                                 ! parameters.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, CloudGammaParamses%nCloudGammaParamses
    If (CloudGammaParams%Name .CIEq. CloudGammaParamses%CloudGammaParamses(i)%Name) Then
      If (CloudGammaParams == CloudGammaParamses%CloudGammaParamses(i)) Then
        Return
      Else
        Call Message(                                                           &
               'FATAL ERROR in adding the set of cloud gamma parameters "'   // &
               Trim(CloudGammaParams%Name)                                   // &
               '": a different set of cloud gamma parameters with the same ' // &
               'name already exists',                                           &
               3                                                                &
             )
      End If
    End If
  End Do

  If (CloudGammaParamses%nCloudGammaParamses >= MaxCGSpecies) Then
    Call Message(                                                         &
           'FATAL ERROR in adding the set of cloud gamma parameters "' // &
           Trim(CloudGammaParams%Name)                                 // &
           '": there are too many sets of cloud gamma parameters',        &
           3                                                              &
         )
  End If

  CloudGammaParamses%nCloudGammaParamses = CloudGammaParamses%nCloudGammaParamses + 1
  CloudGammaParamses%CloudGammaParamses(CloudGammaParamses%nCloudGammaParamses) = CloudGammaParams

End Subroutine AddCloudGammaParams

!-------------------------------------------------------------------------------------------------------------

Function FindCloudGammaParamsIndex(Name, CloudGammaParamses)
! Finds the index of a set of cloud gamma parameters in the collection of all such sets of parameters.

  Implicit None
  ! Argument list:
  Character(*),              Intent(In) :: Name               ! Name of set of cloud gamma parameters.
  Type(CloudGammaParamses_), Intent(In) :: CloudGammaParamses ! Collection of 'cloud gamma parameterses'.
  ! Function result:
  Integer :: FindCloudGammaParamsIndex ! Index of the set of cloud gamma parameters.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, CloudGammaParamses%nCloudGammaParamses
    If (Name .CIEq. CloudGammaParamses%CloudGammaParamses(i)%Name) Then
      FindCloudGammaParamsIndex = i
      Return
    End If
  End Do

  Call Message(                                            &
         'FATAL ERROR: set of cloud gamma parameters "' // &
         Trim(Name)                                     // &
         '" not found',                                    &
         3                                                 &
       )

End Function FindCloudGammaParamsIndex

!-------------------------------------------------------------------------------------------------------------

Function CloudGammaParamsEq(CloudGammaParams1, CloudGammaParams2)
! Tests for equality of cloud gamma parameters

  Implicit None
  ! Argument list:
  Type(CloudGammaParams_), Intent(In) :: CloudGammaParams1 !} Two sets of cloud gamma parameters to be checked
  Type(CloudGammaParams_), Intent(In) :: CloudGammaParams2 !}
  ! Function result:
  Logical :: CloudGammaParamsEq ! Indicates if two sets of cloud gamma parameters are equal.
  ! Locals:
  Logical :: ParamsEq  ! Local copy of function result
  Integer :: i         ! Loop index.

  If (CloudGammaParams1%nEnergies /= CloudGammaParams2%nEnergies) Then
    ParamsEq = .false.
  Else
    ParamsEq = .true.
    Do i = 1, CloudGammaParams1%nEnergies

      ParamsEq =  CloudGammaParams1%PhotonEnergy(i)    ==   CloudGammaParams2%PhotonEnergy(i)         .and. &
            CloudGammaParams1%LinearAttenCoeff(i)      ==   CloudGammaParams2%LinearAttenCoeff(i)     .and. &
            CloudGammaParams1%BBuildUpFactora(i)       ==   CloudGammaParams2%BBuildUpFactora(i)      .and. &
            CloudGammaParams1%BBuildUpFactorb(i)       ==   CloudGammaParams2%BBuildUpFactorb(i)      .and. &
            CloudGammaParams1%AirKermapuFluence(i)     ==   CloudGammaParams2%AirKermapuFluence(i)    .and. &
            CloudGammaParams1%AdEffDosepuAirKerma(i)   ==   CloudGammaParams2%AdEffDosepuAirKerma(i)  .and. &
            CloudGammaParams1%AdThyDosepuAirKerma(i)   ==   CloudGammaParams2%AdThyDosepuAirKerma(i)  .and. &
            CloudGammaParams1%AdLunDosepuAirKerma(i)   ==   CloudGammaParams2%AdLunDosepuAirKerma(i)  .and. &
            CloudGammaParams1%AdBoSDosepuAirKerma(i)   ==   CloudGammaParams2%AdBoSDosepuAirKerma(i)  .and. &
            ParamsEq
    End Do
  End If

  CloudGammaParamsEq = ParamsEq

End Function CloudGammaParamsEq

!-------------------------------------------------------------------------------------------------------------

Function InitMassLimits() Result (MassLimits)
! Initialises a collection of particle mass limits.

  Implicit None
  ! Function result:
  Type(MassLimits_) :: MassLimits ! Initialised collection of particle mass limits.

  MassLimits%nMassLimits = 0

End Function InitMassLimits

!-------------------------------------------------------------------------------------------------------------

Function InitMassLimit(SpeciesName, Limit) Result (MassLimit)
! Initialises a particle mass limit.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: SpeciesName ! Name of species.
  Real(Std),    Intent(In) :: Limit       ! Upper limit on particle mass for this
                                          ! species.
  ! Function result:
  Type (MassLimit_) :: MassLimit ! Initialised particle mass limit.

  If (SpeciesName == ' ') Then
    Call Message(                                                    &
           'FATAL ERROR in InitMassLimit: name of species is blank', &
           3                                                         &
         )
  End If
  If (Len_Trim(SpeciesName) > MaxCharLength) Then
    Call Message(                               &
           'FATAL ERROR in InitMassLimit: '  // &
           'name of species is given as "'   // &
           Trim(SpeciesName)                 // &
           '" and is too long',                 &
           3                                    &
         )
  End If
  MassLimit%SpeciesName = SpeciesName

  ! $$ Fixed met error? Particle masses are ignored for fixed met, but best to give error.
  If (Limit <= 0.0) Then
    Call Message(                                                                 &
           'FATAL ERROR: "Particle Mass Limit" is given as <= 0 for species "' // &
           Trim(SpeciesName)                                                   // &
           '"',                                                                   &
           3                                                                      &
         )
  End If
  MassLimit%Limit = Limit

End Function InitMassLimit

!-------------------------------------------------------------------------------------------------------------

Subroutine AddMassLimit(MassLimit, MassLimits)
! Adds a particle mass limit to a collection of particle mass limits.

  Implicit None
  ! Argument list:
  Type(MassLimit_),  Intent(In)    :: MassLimit  ! Particle mass limit to be added.
  Type(MassLimits_), Intent(InOut) :: MassLimits ! Collection of particle mass limits.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, MassLimits%nMassLimits
    If (MassLimit%SpeciesName .CIEq. MassLimits%MassLimits(i)%SpeciesName) Then
      If (MassLimit%Limit == MassLimits%MassLimits(i)%Limit) Then
        Return
      Else
        Call Message(                                                          &
               'FATAL ERROR in adding the particle mass limit for species"' // &
               Trim(MassLimit%SpeciesName)                                  // &
               '": a different particle mass limit already exists',            &
               3                                                               &
             )
      End If
    Endif
  End Do

  If (MassLimits%nMassLimits >= MaxSpecieses) Then
    Call Message(                                                           &
           'FATAL ERROR in adding the particle mass limit for species "' // &
           Trim(MassLimit%SpeciesName)                                   // &
           '": there are too many particle mass limits',                    &
           3                                                                &
         )
  End If

  MassLimits%nMassLimits                        = MassLimits%nMassLimits + 1
  MassLimits%MassLimits(MassLimits%nMassLimits) = MassLimit

End Subroutine AddMassLimit

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpSpecieses(SizeDists, Specieses)
! Sets up Specieses.

  Implicit None
  ! Argument list:
  Type(SizeDists_), Intent(In),    Target :: SizeDists ! Collection of particle size distribution.
  Type(Specieses_), Intent(InOut), Target :: Specieses ! Collection of species.
  ! Locals:
  Type(SizeDist_),    Pointer :: SizeDist     !} Abbreviations.
  Type(SpeciesUses_), Pointer :: SpeciesUses  !}
  Integer                     :: iSpeciesUses ! Index of set of species uses.

  Specieses%nParticleSpecieses = 0
  Specieses%nFields            = 0
  Do iSpeciesUses = 1, Specieses%nSpeciesUseses
    SpeciesUses => Specieses%SpeciesUseses(iSpeciesUses)
    If (SpeciesUses%OnParticles) Specieses%nParticleSpecieses = Specieses%nParticleSpecieses + 1
    If (SpeciesUses%OnFields) Then
      If (SpeciesUses%SizeDistName /= ' ') Then
        SizeDist => SizeDists%SizeDists(FindSizeDistIndex(SpeciesUses%SizeDistName, SizeDists))
        If (.not. SizeDist%DensityPresent) Then
          Call Message (                                                          &
                 'FATAL ERROR: The particle size distribution "'               // &
                 Trim(SpeciesUses%SizeDistName)                                // &
                 '" referred to in the Species Uses input block for Species "' // &
                 Trim(SpeciesUses%SpeciesName)                                 // &
                 '" is not suitable for use in the Species Uses input block.',    &
                 3                                                                &
               )
        End If
        Specieses%nFields = Specieses%nFields + SizeDist%nSizeRanges
      Else
        Specieses%nFields = Specieses%nFields + 1
      End If
    End If    
  End Do
  
End Subroutine SetUpSpecieses

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpSpecieses_MassLimits(MassLimits, Specieses)
! Sets up Specieses using information from MassLimits.

  Implicit None
  ! Argument list:
  Type(MassLimits_), Intent(In),    Target :: MassLimits ! Collection of particle mass
                                                         ! limits.
  Type(Specieses_),  Intent(InOut), Target :: Specieses  ! Collection of specieses.
  ! Locals:
  Type(Species_),   Pointer :: Species   !} Abbreviations for species and particle
  Type(MassLimit_), Pointer :: MassLimit !} mass limits.
  Integer                   :: i         !] Loop indices.
  Integer                   :: j         !]

  Do i = 1, Specieses%nSpecieses

    Species => Specieses%Specieses(i)

    Species%UseMassLimit = .false.

    Do j = 1, MassLimits%nMassLimits

      MassLimit => MassLimits%MassLimits(j)

      If (MassLimit%SpeciesName .CIEq. Species%Name) Then
        Species%UseMassLimit = .true.
        Species%MassLimit    = MassLimit%Limit
        Exit
      End If

    End Do

  End Do

End Subroutine SetUpSpecieses_MassLimits

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpSpecieses_DecayChains(Specieses)
! Sets up full decay chains in Specieses.

  Implicit None
  ! Argument list:
  Type(Specieses_), Intent(InOut), Target :: Specieses  ! Collection of specieses.
  ! Locals:
  Type(Species_), Pointer :: Species       ! Abbreviation for species.
  Integer                 :: i             ! Loop variable over species.
  Integer                 :: k             ! Generation of a daughter product.
  Integer                 :: iSpecies      ! Species index of a daughter product.
  Integer                 :: iSpeciesUses  ! SpeciesUses index of a daughter product.
  Integer                 :: j             ! Index of SpeciesUses associated with species 
  Integer                 :: iSizeDist0    ! Size distribution associated with parent species

  ! Loop over all species, and construct decay chain for those species having a daughter product.
  Do i = 1, Specieses%nSpecieses

    Species => Specieses%Specieses(i)

    Species%DecayChain(1) = i

    Species%DecayChainLength = 1

    j = Specieses%iSpecies2SpeciesUses(i) !$$ need to use findspeciesusesindex etc and rename/define some local vars
    If (j == 0) Cycle        !$$ Also call to this routine is now after SetUp_iSpecieses which works but breaks convention
    iSizeDist0 = Specieses%SpeciesUseses(j)%iSizeDist

    If (Species%HasDaughter) Then

      ! Need all decay chain species to be carried on particles.
      iSpeciesUses = FindSpeciesUsesIndex(Species%Name, Specieses)
      If ( .not. Specieses%SpeciesUseses(iSpeciesUses)%OnParticles ) Then
        Call Message(                                                                    &
               'FATAL ERROR: Species involved in radioactive decay chains, such as "' // &
               Trim(Species%Name)                                                     // &
               '", must be carried on particles',                                        &
               3                                                                         &
             )
      End If

      k = 2
      iSpecies = i

      Do While (Specieses%Specieses(iSpecies)%HasDaughter)

        If (k > MaxDecayChainLength) Then
          Call Message(                                                  &
                 'FATAL ERROR: decay chain is too long for species "' // &
                 Trim(Species%Name)                                   // &
                 '"',                                                    &
                 3                                                       &
               )
        End If

        ! $$ Needs all daughter species to exist in the input file. Update the FindSpeciesIndex function
        ! to have option of returning error code if species is missing rather than stop with a fatal error.
        ! Need all decay chain species to be carried on particles.
        iSpeciesUses = FindSpeciesUsesIndex(Specieses%Specieses(iSpecies)%Daughter, Specieses)
        If ( .not. Specieses%SpeciesUseses(iSpeciesUses)%OnParticles ) Then
          Call Message(                                                                    &
                 'FATAL ERROR: Species involved in radioactive decay chains, such as "' // &
                 Trim(Specieses%Specieses(iSpecies)%Daughter)                           // &
                 '", must be carried on particles',                                        &
                 3                                                                         &
               )
        End If
        iSpecies = FindSpeciesIndex(Specieses%Specieses(iSpecies)%Daughter, Specieses)

        ! Test that radioactive daughter species on fields have the same particle size distribution as parent species
        j = Specieses%iSpecies2SpeciesUses(iSpecies)
        If ( Specieses%SpeciesUseses(j)%iSizeDist /= iSizeDist0 ) Then 
          Call Message(                                                                            &
                 'FATAL ERROR: Species involved in radioactive decay chains, such as "'         // &
                 Trim(Specieses%Specieses(iSpecies)%Daughter)                                   // &
                 '", that are carried on fields with an associated particle size distribution ' // &
                 'must have the same size distribution as the parent species',                     &
                 3                                                                                 &
               )
        End If

        Species%DecayChain(k)    = iSpecies
        Species%DecayChainLength = k
        k = k + 1

      End Do

    End If

  End Do

End Subroutine SetUpSpecieses_DecayChains

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpiSpecieses(SizeDists, Specieses)
! Sets up indices in Specieses.

  Implicit None
  ! Argument list:
  Type(SizeDists_), Intent(In),    Target :: SizeDists ! Collection of particle size distributions.
  Type(Specieses_), Intent(InOut), Target :: Specieses ! Collection of species.
  ! Locals:
  Integer                     :: iSpecies     !} Indices of species, set of species uses, species stored on
  Integer                     :: iSpeciesUses !} particles, field and particulate size range.
  Integer                     :: iParticle    !}
  Integer                     :: iField       !}
  Integer                     :: iSizeRange   !}
  Type(SizeDist_),    Pointer :: SizeDist     !] Abbreviations.
  Type(SpeciesUses_), Pointer :: SpeciesUses  !]
  
  ! Size distributions.
  Do iSpeciesUses = 1, Specieses%nSpeciesUseses
    SpeciesUses => Specieses%SpeciesUseses(iSpeciesUses)
    If (SpeciesUses%SizeDistName == ' ') Then
      SpeciesUses%iSizeDist = 0
    Else
      SpeciesUses%iSizeDist = FindSizeDistIndex(SpeciesUses%SizeDistName, SizeDists)
    End If
  End Do
  
  ! Initialise indices.
  Specieses%iSpecies2SpeciesUses     (:)    = 0
  Specieses%iSpeciesUses2Species     (:)    = 0
  Specieses%iSpecies2Particle        (:)    = 0
  Specieses%iParticle2Species        (:)    = 0
  Specieses%iSpeciesUses2Particle    (:)    = 0
  Specieses%iParticle2SpeciesUses    (:)    = 0
  Specieses%iSpeciesAndSize2Field    (:, :) = 0
  Specieses%iSpeciesUsesAndSize2Field(:, :) = 0
  Specieses%iField2Species           (:)    = 0
  Specieses%iField2SpeciesUses       (:)    = 0
  Specieses%iField2Size              (:)    = 0
  Specieses%iSpecies2SizeDist        (:)    = 0
  Specieses%iSpeciesUses2SizeDist    (:)    = 0
  Specieses%iField2SizeDist          (:)    = 0 
 
  ! Indices associated with species, species uses and species on particles.

  iParticle = 0
  Do iSpeciesUses = 1, Specieses%nSpeciesUseses
    SpeciesUses => Specieses%SpeciesUseses(iSpeciesUses)
    iSpecies = FindSpeciesIndex(SpeciesUses%SpeciesName, Specieses)
    Specieses%iSpecies2SpeciesUses(iSpecies)     = iSpeciesUses
    Specieses%iSpeciesUses2Species(iSpeciesUses) = iSpecies
    If (SpeciesUses%OnParticles) Then
      iParticle                                     = iParticle + 1
      Specieses%iSpecies2Particle(iSpecies)         = iParticle
      Specieses%iParticle2Species(iParticle)        = iSpecies
      Specieses%iSpeciesUses2Particle(iSpeciesUses) = iParticle
      Specieses%iParticle2SpeciesUses(iParticle)    = iSpeciesUses
    End If    
  End Do

  ! Indices associated with species, species uses and fields.
  ! $$ Non-sedimenting advected fields need to be grouped together currently.
  ! $$ Blocks (i), (ii) and (iii) would be simpler (unified) if don't need to do this.
  ! $$ Might simplify code elsewhere too and allow more flexibility.
  ! $$ But keeping them together might be more efficient (?).

  iField = 0
    
  ! (i) Non-sedimenting advected fields.
  Do iSpeciesUses = 1, Specieses%nSpeciesUseses
    SpeciesUses => Specieses%SpeciesUseses(iSpeciesUses)
    iSpecies = Specieses%iSpeciesUses2Species(iSpeciesUses)
    If (SpeciesUses%OnFields .and. SpeciesUses%AdvectField) Then
      If (SpeciesUses%iSizeDist == 0) Then
        iField = iField + 1
        Specieses%iSpeciesAndSize2Field    (iSpecies,     1) = iField
        Specieses%iSpeciesUsesAndSize2Field(iSpeciesUses, 1) = iField
        Specieses%iField2Species           (iField         ) = iSpecies
        Specieses%iField2SpeciesUses       (iField         ) = iSpeciesUses
        Specieses%iField2Size              (iField         ) = 1
        Specieses%iSpecies2SizeDist        (iSpecies       ) = 0
        Specieses%iSpeciesUses2SizeDist    (iSpeciesUses   ) = 0
        Specieses%iField2SizeDist          (iField         ) = 0
      End If
    End If
  End Do

  Specieses%iFirstAdvectedNonSedimentingField = 1  
  Specieses%iLastAdvectedNonSedimentingField  = iField  
  
  ! (ii) Sedimenting advected fields.
  Do iSpeciesUses = 1, Specieses%nSpeciesUseses
    SpeciesUses => Specieses%SpeciesUseses(iSpeciesUses)
    iSpecies = Specieses%iSpeciesUses2Species(iSpeciesUses)
    If (SpeciesUses%OnFields .and. SpeciesUses%AdvectField) Then
      If (SpeciesUses%iSizeDist > 0) Then
        SizeDist => SizeDists%SizeDists(SpeciesUses%iSizeDist)
        Do iSizeRange = 1, SizeDist%nSizeRanges
          iField = iField + 1
          Specieses%iSpeciesAndSize2Field    (iSpecies,     iSizeRange) = iField
          Specieses%iSpeciesUsesAndSize2Field(iSpeciesUses, iSizeRange) = iField
          Specieses%iField2Species           (iField                  ) = iSpecies
          Specieses%iField2SpeciesUses       (iField                  ) = iSpeciesUses
          Specieses%iField2Size              (iField                  ) = iSizeRange
          Specieses%iSpecies2SizeDist        (iSpecies                ) = SpeciesUses%iSizeDist
          Specieses%iSpeciesUses2SizeDist    (iSpeciesUses            ) = SpeciesUses%iSizeDist
          Specieses%iField2SizeDist          (iField                  ) = SpeciesUses%iSizeDist
        End Do
      End If
    End If
  End Do

  Specieses%AdvectedFields            = iField > 0
  Specieses%AdvectedSedimentingFields = iField > Specieses%iLastAdvectedNonSedimentingField

  ! (iii) Non-sedimenting non-advected fields. 
  Do iSpeciesUses = 1, Specieses%nSpeciesUseses
    SpeciesUses => Specieses%SpeciesUseses(iSpeciesUses)
    iSpecies = Specieses%iSpeciesUses2Species(iSpeciesUses)
    If (SpeciesUses%OnFields .and. (.not. SpeciesUses%AdvectField)) Then
      If (SpeciesUses%iSizeDist == 0) Then
        iField = iField + 1
        Specieses%iSpeciesAndSize2Field    (iSpecies,     1) = iField
        Specieses%iSpeciesUsesAndSize2Field(iSpeciesUses, 1) = iField
        Specieses%iField2Species           (iField         ) = iSpecies
        Specieses%iField2SpeciesUses       (iField         ) = iSpeciesUses
        Specieses%iField2Size              (iField         ) = 1
        Specieses%iSpecies2SizeDist        (iSpecies       ) = 0
        Specieses%iSpeciesUses2SizeDist    (iSpeciesUses   ) = 0
        Specieses%iField2SizeDist          (iField         ) = 0
      Else
        SizeDist => SizeDists%SizeDists(SpeciesUses%iSizeDist)
        Do iSizeRange = 1, SizeDist%nSizeRanges
          iField = iField + 1
          Specieses%iSpeciesAndSize2Field    (iSpecies,     iSizeRange) = iField
          Specieses%iSpeciesUsesAndSize2Field(iSpeciesUses, iSizeRange) = iField
          Specieses%iField2Species           (iField                  ) = iSpecies
          Specieses%iField2SpeciesUses       (iField                  ) = iSpeciesUses
          Specieses%iField2Size              (iField                  ) = iSizeRange
          Specieses%iSpecies2SizeDist        (iSpecies                ) = SpeciesUses%iSizeDist
          Specieses%iSpeciesUses2SizeDist    (iSpeciesUses            ) = SpeciesUses%iSizeDist
          Specieses%iField2SizeDist          (iField                  ) = SpeciesUses%iSizeDist
        End Do
      End If
    End If
  End Do

  If (Specieses%nParticleSpecieses /= iParticle .or. Specieses%nFields /= iField) Then
    Call Message('UNEXPECTED FATAL ERROR in SetUpiSpecieses', 4)
  End If

End Subroutine SetUpiSpecieses

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpSpecieses_iCloudGammaParams(Specieses, CloudGammaParamses)
! Sets up indices in Species for referring to sets of cloud gamma parameters.

  Implicit None
  ! Argument list:
  Type(Specieses_),          Intent(InOut) :: Specieses           ! Collection of species.
  Type(CloudGammaParamses_), Intent(In)    :: CloudGammaParamses  ! Collection of sets of cloud gamma
                                                                  ! parameters.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Specieses%nSpecieses

    If (Specieses%Specieses(i)%CloudGammaParameters == ' ') Then
      Specieses%Specieses(i)%iCloudGammaParams = 0
    Else
      Specieses%Specieses(i)%iCloudGammaParams = FindCloudGammaParamsIndex(                     &
                                                   Specieses%Specieses(i)%CloudGammaParameters, &
                                                   CloudGammaParamses                           &
                                                 )
    End If

  End Do

End Subroutine SetUpSpecieses_iCloudGammaParams

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpSpecieses_iMaterialUnits(Specieses, MaterialUnits)
! Sets up indices in Species for referring to material units.

  Implicit None
  ! Argument list:
  Type(Specieses_),     Intent(InOut) :: Specieses     ! Collection of species.
  Type(MaterialUnits_), Intent(In)    :: MaterialUnits ! Collection of material units.

  ! Locals:
  Integer                      :: i                 ! Loop index.
  Integer                      :: iMaterialUnit     ! Index of material unit
  Character(len=MaxCharLength) :: MaterialUnitName  ! Name of material unit

  Do i = 1, Specieses%nSpecieses
    MaterialUnitName = Specieses%Specieses(i)%MaterialUnitName
    iMaterialUnit = FindMaterialUnitIndex(MaterialUnitName, MaterialUnits)
    If ( iMaterialUnit == -1 ) Then
      Call Message('ERROR: Unknown material unit '   // &
                   Trim(MaterialUnitName)            // &
                   ' in species '                    // &
                   Trim(Specieses%Specieses(i)%Name) // &
                   '.',                                 &
                   3                                    &
           )
    Else
      Specieses%Specieses(i)%iMaterialUnit = iMaterialUnit
    End If
  End Do

End Subroutine SetUpSpecieses_iMaterialUnits

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcDecayFactor(Species, Age, dT, DecayFactor)
! Calculates decay factor for decay of species.

! $$ currently used in output for decay of deposition - use in particle.f90 too (or have a 'time_'
!    version and a 'real' version)?
! $$ other sorts of decay?

  Implicit None
  ! Argument list:
  Type(Species_),   Intent(In)  :: Species     ! Species.
  Type(ShortTime_), Intent(In)  :: Age         !} Defines time interval for decay to be [Age, Age + dT]. Age
  Type(ShortTime_), Intent(In)  :: dT          !} and dT must be >= 0 (and can be infinite).
  Real(Std),        Intent(Out) :: DecayFactor ! Decay factor.
  ! Locals:
  Real(Std) :: RdT  ! Real version of dT.
  Real(Std) :: RAge ! Real version of Age.

  ! Power law decay.
  If (Species%UsePowerLawDecay) Then

    If (IsInfFuture(dT) .or. IsInfFuture(Age)) Then
      DecayFactor = 0.0
    Else If (dT == ZeroShortTime()) Then
      DecayFactor = 1.0
    Else
      RdT         = ShortTime2RealTime(dT)
      RAge        = ShortTime2RealTime(Age)
      DecayFactor = (                                                                                     &
                      Max(Species%PowerLawDecayDelay, RAge) / Max(Species%PowerLawDecayDelay, RAge + RDt) &
                    ) ** Species%PowerLawDecayExponent
    End If

  ! Exponential decay.
  Else If (Species%InvHalfLife > 0.0) Then

    If (IsInfFuture(dT)) Then
      DecayFactor = 0.0
    Else If (dT == ZeroShortTime()) Then
      DecayFactor = 1.0
    Else
      RdT         = ShortTime2RealTime(dT)
      DecayFactor = Exp( - RdT * Species%InvHalfLife * LOG(2.0))
    End If

  ! No decay.
  Else

    DecayFactor = 1.0

  End If

End Subroutine CalcDecayFactor

!-------------------------------------------------------------------------------------------------------------

Function CalcVd(                                           &
           Species, Flow, Surface, Plant, Time, Cloud3d,   &
           WSed, Diameter, Z, Zs, YLatLong                 &
         )                                                 &
Result(Vd)
! Calculates dry deposition velocity (including sedimentation contribution)

  Implicit None
  ! Argument list:
  Type(Species_), Intent(In) :: Species  ! A species.
  Type(Flow_),    Intent(In) :: Flow     ! Flow information.
  Type(Surface_), Intent(In) :: Surface  ! Surface information.
  Type(Plant_),   Intent(In) :: Plant    ! Plant information.
  Type(Time_),    Intent(In) :: Time     ! Time
  Real(Std),      Intent(In) :: Cloud3d  ! 3-d cloud amount (fraction).
  Real(Std),      Intent(In) :: WSed     ! Gravitational settling velocity.
  Real(Std),      Intent(In) :: Diameter ! Diameter of sedimenting particles.
  Real(Std),      Intent(In) :: Z        ! z position.
  Real(Std),      Intent(In) :: Zs       ! Deposition height
  Real(Std),      Intent(In) :: YLatLong ! Latitude position in lat-long coordinate system

  ! Function result:
  Real(Std) :: Vd ! Dry deposition velocity.
  ! Locals:
  Real(Std)            :: Psi
  Real(Std)            :: Psi0
  Real(Std)            :: Phi
  Real(Std)            :: Z0
  Real(Std)            :: ZRef
  Real(Std)            :: Ra
  Real(Std)            :: Rb
  Real(Std)            :: Rs
  Real(Std)            :: Rc(9)        ! Surface resistance for each land use type
  Real(Std)            :: VdLandUse(9) ! Vd for each land use type
  Integer              :: iLandUseType ! Loop variable over land use types
  Real(Std), Parameter :: VdHydrogen(12) =                                                  &
       (/ 2.7e-4, 1.8e-4, 2.0e-4, 3.0e-4, 4.7e-4, 5.5e-4, 6.4e-4, 4.5e-4, 2.7e-4, 3.4e-4, 3.2e-4, 2.3e-4 /)
!      Estimated monthly hydrogen deposition velocity m/s for N Hemisphere derived by EuroHydros project
  Real(Std), Parameter :: VdMax = 0.1

  ! Z0. $$ Here we limit Z0 to avoid too large a Z0 being used if Z0 has orographic
  ! enhancement. This should be unnecessary once we decide exactly what the met data Z0
  ! represents or if we introduce two z0's.
  Z0 = Min(Flow%Z0, 1.0)

  ZRef = Zs/2.0

! Calculate stability correction function

! Unstable

  If ( Flow%RecipLMO < 0.0 ) Then
    Psi  = ( 1.0 - 16.0 * (ZRef + Z0) * Flow%RecipLMO )**0.25
    Psi0 = ( 1.0 - 16.0 * Z0 * Flow%RecipLMO )**0.25
    Phi  = 2.0 * Log( (1.0 + Psi * Psi)/(1.0 + Psi0 * Psi0) )

! Stable

  Else If ( (ZRef + Z0) * Flow%RecipLMO < 1.0 ) Then
    Phi = -5.0 * ZRef * Flow%RecipLMO
  Else If ( Z0 * Flow%RecipLMO > 1.0 ) Then
    Phi = -5.0 * Log( (ZRef + Z0)/Z0 )
  Else
    Phi = -5.0 * (1.0/Flow%RecipLMO - Z0) * Flow%RecipLMO - &
          5.0 * Log( (ZRef + Z0) * Flow%RecipLMO )
  End If


! Calculate aerodynamic resistance
  Ra = Max( 1.0, ( Log( (ZRef + Z0)/Z0 ) - Phi )/( VK * Flow%UStar ) )
  
! Calculate default laminar resistance
  Rb = 8.0/Flow%UStar


  ! If test here requires either 
  ! (a) a non-zero deposition velocity to be specified for the species, or
  ! (b) for non-sedimenting sources only (diameter = 0), either 
  !     (i) a deposition velocity (>=0) to be specified for the species, or
  !     (ii) all dry deposition options (deposition velocity, surface resistance, mean 
  !          aerosol diameter or land use dry dep) to be left blank / set to false for the species
  If (Species%UseVd .and. ((Diameter == 0.0) .or. Species%Vd > 0.0)) Then

    ! Non-sedimenting only
    If (Species%Vd == 0.0) Then

      Vd = 0.0  ! WSed = 0 

    ! Non-zero Vd specified for species.
    ! Applies to sedimenting and non-sedimenting alike
    Else

      If (WSed / Species%Vd > 0.01) Then
        Vd = WSed / (1.0 - Exp(- WSed / Species%Vd))
      Else
        Vd = Species%Vd
      End If

      ! Specific Dry Deposition for Hydrogen
      ! Vd set from monthly EuroHydros values.
      ! $$ note this works only when Vd specified and non-zero. Otherwise code won't get to this point.
      ! Do not deposit over sea (defined as Absolute Topography less than 1.e-6)
      ! If over land and surface temperature is less than 0 degC but greater than -15 degC then halve Vd
      ! If over land and surface temperature is less than -15 degC then quarter Vd
      ! These figures come from Price et al 2007 JGR Vol 112
      !
      ! Requires Vd (>0) to be specified for HYDROGEN
      If (Species%Name(1:8) .CIEq. 'HYDROGEN') Then

        If (ABS(Flow%Topog) < 1.0e-6) Then
          Vd = 0.0
        Else
          If (Flow%T0 >= TKAtTCEq0) Then
            Vd = VdHydrogen(Time%Month)
          Else If (Flow%T0 > (TKAtTCEq0 - 15.0)) Then
            Vd = VdHydrogen(Time%Month) / 2.0
          Else
            Vd = VdHydrogen(Time%Month) / 4.0
          End If
        End If

      Else

        ! Specify dry deposition rate over sea for ozone      ! $$
        !
        ! Requires Vd (>0) to be specified for O3
        If (Species%Name(1:2) .CIEq. 'O3') Then
          If (ABS(Flow%Topog) < 1.0e-6) Then
            Vd = 7.0E-4
          End If
        End If

      End If

    End If

  ! This section invoked if either
  ! (a) the land use dry dep flag is set to true, a surface resistance is given or 
  !     a mean aerosol diameter is given for the species, or
  ! (b) for sedimenting sources (diameter > 0) only, either
  !     (i) a zero deposition velocity is given for the species (Vd = 0), or
  !     (ii) all dry deposition options for the species (deposition velocity, surface resistance, 
  !           mean aerosol diameter or land use dry dep) are left blank / set to false
  Else

    ! Occult deposition.
    If (Cloud3d > 0.5 .and. Z < 50.0) Then

      Vd = 0.8/Ra
      If (WSed / Vd > 0.01) Then
        Vd = WSed / (1.0 - Exp(- WSed / Vd))
      End If

    ! Normal case.
    Else

      ! New surface resistance parameterisation from STOCHEM
      ! Invoked if the land use dry dep flag is set to true.
      If (Species%LandUseDryDep) Then
        Rc = CalcSurfResistance(Species, Flow, Surface, Plant, YLatLong)

        Vd = 0.0

        Do iLandUseType = 1, 9
          If (Surface%LandUseFracs(iLandUseType) > 0.0) Then
            If (Rc(iLandUseType) < Huge(Zs) / 10.0) Then
              VdLandUse(iLandUseType) = 1.0 / (Ra + Rb + Rc(iLandUseType))
              !$$ - should limit on vd be applied? $$ Particularly for new surface resistance scheme
              If (WSed / VdLandUse(iLandUseType) > 0.01) Then
                VdLandUse(iLandUseType) = WSed / (1.0 - Exp(-WSed / VdLandUse(iLandUseType)))
              End If
            Else
              VdLandUse(iLandUseType) = WSed
            End If
            Vd = Vd + Surface%LandUseFracs(iLandUseType) * VdLandUse(iLandUseType)
          End If
        End Do



      Else

        ! Aerosol parametrisation: invoked if either
        ! (a) mean aerosol diameter set (>0) for the species, or
        ! (b) sedimenting source (diameter > 0)
        ! Parameterisation uses particle diameter, if available, or aerosol mean diameter otherwise
        If ((Species%AerosolDiameter > 0.0) .or. (Diameter > 0.0)) Then

          Rb = 300.0 / Flow%Ustar

          If (Diameter > 0.0) Then
            If (Diameter > 1.0) Rb = Rb / (Diameter * Diameter)
          Else If (Species%AerosolDiameter > 1.0) Then
            Rb = Rb / (Species%AerosolDiameter * Species%AerosolDiameter)
          End If

          Vd = Min(1.0 / (Ra + Rb), VdMax) !$$ - should limit on vd be removed?
                                           !$$ - should it be applied for new surface resistance scheme above?

        Else

          Rs = Species%Rs

          Vd = Min(1.0 / (Ra + Rb + Rs), VdMax) !$$ - should limit on vd be removed?
                                                !$$ - should it be applied for new surface resistance scheme above?
        End If
        If (WSed / Vd > 0.01) Then
          Vd = WSed / (1.0 - Exp(- WSed / Vd))
        End If
      End If

    End If

  End If


End Function CalcVd

!-------------------------------------------------------------------------------------------------------------

Function CalcSurfResistance(Species, Flow, Surface, Plant, YLatLong) Result(Rc)


! New surface resistance parameterisation from STOCHEM

  Implicit None
  ! Argument list:
  Type(Species_),   Intent(In)    :: Species
  Type(Flow_),      Intent(In)    :: Flow
  Type(Surface_),   Intent(In)    :: Surface
  Type(Plant_),     Intent(In)    :: Plant
  Real(Std),        Intent(In)    :: YLatLong
  ! Function result:
  Real(Std) :: Rc(9)    ! Surface resistance for each land use type

  ! Locals:
  Real(Std)             :: RelHumidity        ! Relative humidity
  Real(Std)             :: SO4_vd             ! Aerosol deposition velocity
  Real(Std)             :: RStom(5)           ! Stomatal resistance for 5 plant types
  Real(Std)             :: RCut(5)            ! Cuticular resistance for 5 plant types
  Real(Std)             :: RcTemp             ! Temporary variable
  Integer               :: iLandUseType       ! Loop over Land Use type
  Real(Std)             :: NO2Standard(9)     ! Standard surface resistance for NO2
  Real(Std)             :: COStandard(9)      ! Standard surface resistance for CO
  Real(Std)             :: CH4Standard(9)     ! Standard surface resistance for CH4
  Real(Std)             :: O3Standard(9)      ! Standard surface resistance for O3
  Real(Std)             :: PANStandard(9)     ! Standard surface resistance for PAN (peroxyacetyl nitrate
                                              ! (CH3COO2NO2))
!  Real(Std)             :: CH3OOHStandard(9)  ! Standard surface resistance for CH3OOH
  Real(Std)             :: SO2Standard(9)     ! Standard surface resistance for SO2
  Real(Std)             :: SAStandard(9)      ! Standard surface resistance for SA
  Real(Std)             :: HNO3Standard(9)    ! Standard surface resistance for HNO3 (nitric acid vapour)
!  Real(Std)             :: ROOHStandard(9)    ! Standard surface resistance for ROOH (organic hydroperoxides)

  Real(Std),  Parameter :: CH4UptakeFlux(9)  = (/ 39.5, 50.0, 30.0, 37.0, 27.5, 0.0, 0.0,   &
                                                  27.5, 0.0 /)
                                                  ! Methane uptake fluxes (ug/m2/hr)
  Real(Std),  Parameter :: H2dd_c(5)         = (/ 0.00197, 0.00197, 0.00177, 1.2346, 0.0001 &
                                               /)
                                                  ! Hydrogen coefficients
  Real(Std),  Parameter :: H2dd_m(5)         = (/ -0.00419, -0.00419, -0.00414, -0.472,     &
                                                  0.0 /)
                                                  ! Hydrogen coefficients
  Real(Std),  Parameter :: H2dd_q            = 0.27
                                                  ! Hydrogen coefficient
  Real(Std),  Parameter :: glmin             = 1.0E-6
                                                  ! Minimum stomatal conductance required for
                                                  ! stomatal deposition
  Real(Std),  Parameter :: mml               = 1.008E5
                                                  ! Factor to convert methane flux to dry depositiom
  Real(Std),  Parameter :: RUGC              = 8.314
                                                  ! Universal gas constant (J K-1 MOL-1)
  Real(Std),  Parameter :: CH4dd_tun(4)      = (/ -4.757E-6, 4.0288E-3, -1.13592, 106.636 /)
                                                  ! CH4 coefficients for tundra
  Real(Std),  Parameter :: hno3dd_ice(3)     = (/ -13.57, 6841.9, -857410.6 /)
                                                  ! HNO3 coefficients for ice
  Real(Std),  Parameter :: so2dd_ice(3)      = (/ 0.0001, 0.003308, 0.1637 /)
                                                  ! SO2 coefficients for ice


! Note that if additional parameterisations are added for other species then the list of species in
! InitSpecies
! needs to be updated so that if LandUseDryDep is set to be true in the NAMEIII input file and the
! species is not
! a species for which a parameterisation is given in this subroutine then LandUseDryDep is changed to false.


! Set default values for Rc(9), RStom(5) and RCut(5)
  Rc = (/ Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0,    &
          Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0 /)
  RStom = (/ Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0 /)
  RCut = (/ Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0 /)

! Set default values for standard surface resistances for NO2, CO, CH4, O3, PAN, CH3OOH, SO2, SA, HNO3

  NO2Standard    = (/ 225.0, 225.0, 400.0, 400.0, 600.0, 1200.0, 2600.0, 1200.0, 3500.0 /)
  COStandard     = (/ 3700.0, 7300.0, 4550.0, 1960.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0,    &
                      Huge(Rc(1))/3.0, Huge(Rc(1))/3.0 /)
  CH4Standard    = (/ Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0,  &
                      Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0 /)
  O3Standard     = (/ 200.0, 200.0, 200.0, 200.0, 400.0, 800.0, 2200.0, 800.0, 2500.0 /)
  PANStandard    = (/ 500.0, 500.0, 500.0, 500.0, 500.0, Huge(Rc(1))/3.0, 12500.0, 500.0, 12500.0 /)
!  CH3OOHStandard = (/ 30.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0 /)
  SO2Standard    = (/ 100.0, 100.0, 150.0, 350.0, 400.0, 400.0, 10.0, 700.0, Huge(Rc(1))/3.0 /)
  SAStandard     = (/ Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0, Huge(Rc(1))/3.0,  &
                      Huge(Rc(1))/3.0, 1000.0, Huge(Rc(1))/3.0, 20000.0 /)
  HNO3Standard   = (/ 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0 /)
!  ROOHStandard   = (/ 30.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0 /)

! Note that the standard surface resistances for some other species are set equal to
! those standard ones above as follows:-
!
! CH4 = H2
! CH3OOH = C3H7OOH = C2H5OOH = C4H9OOH = ISOPOOH = MVKOOH
! SA = AMMSUL = NAER = MSA = ORGNIT
! HNO3 = NH3 = H2O2 (hydrogen peroxide)
!
! NO and Formaldehyde (HCHO) resistance is set to be high (i.e. default)
! HNO3 - bulk surface resistance of 10 s m-1 to avoid unrealistically high estimates of deposition
! velocities over unusually rough surfaces  (Wesely 1989)
!


! Calculate relative humidity
  RelHumidity = CalcRH(Flow%Q, Flow%T, Flow%P)

  Do iLandUseType = 1, 5
    If (Surface%LandUseFracs(iLandUseType) > 1.0) Then ! $$
      Print*,'Land Use Fraction for Land Surface type ',iLandUseType,' > 1.0'
      Call Message(                                                                          &
        'FATAL ERROR: Land Use Fraction > 1.0 ', 3)
    End If
    If ((Surface%LandUseFracs(iLandUseType) > 0.0) .and. (Plant%StomataConduct(iLandUseType) > 1.0)) Then
      ! $$ This shouldn't happen -
      ! stomatal conductance is set to 1E6 but for ice
      Print*,'Stomatal conductance too large for non-zero land use fraction'
      Call Message(                                                                          &
        'FATAL ERROR: Stomatal conductance > 1.0 ', 3)
    End If
  End Do

! Ozone
  If (Species%Name(1:2) .CIEq. 'O3') Then

    Do iLandUseType = 1, 5
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then

        ! Stomatal Resistance
        If (Plant%StomataConduct(iLandUseType) > glmin) Then
          RStom(iLandUseType) = 1.6 / Plant%StomataConduct(iLandUseType)
        End If

        ! Cuticular Resistance
        If (Plant%LAI(iLanduseType) > 0.0) Then
          RCut(iLandUseType) = 5000.0 / Plant%LAI(iLandUseType)
        End If

        ! Surface Resistance
        ! wet surfaces
        ! Soil moisture is in kg / m2 = mm. To convert to fraction divide by depth of top layer (100mm)
        ! So soil moisture content is %
        If (Surface%SoilMoisture > 30.0) Then
          O3Standard(iLandUseType) = 500.0
        ! dry surfaces
        ! Adjust standard resistance for tundra
        Else If ((iLandUseType == 5) .and. (YLatLong > 60.0)) Then
          O3Standard(iLandUseType) = 800.0
        End If
        Rc(iLandUseType) = 1.0 / (1.0 / RStom(iLandUseType) + 1.0 / RCut(iLandUseType) &
                         + 1.0 / O3Standard(iLandUseType))

      End If
    End Do

    Do iLandUseType = 6, 9
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then
        ! Adjust standard resistance for tundra
        If ((iLandUseType == 8) .and. (Surface%SoilMoisture > 30.0)) Then
          O3Standard(iLandUseType) = 500.0
        End If
        Rc(iLandUseType) = O3Standard(iLandUseType)
      End If
    End Do



! SO2
  Else If (Species%Name(1:15) .CIEq. 'SULPHUR-DIOXIDE') Then

    Do iLandUseType = 1, 5
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then

        ! Stomatal Resistance
        If (Plant%StomataConduct(iLandUseType) > glmin) Then
          RStom(iLandUseType) = 1.9 / Plant%StomataConduct(iLandUseType)
        End If

        ! Cuticular Resistance
        If (Flow%T0 < 268.0) Then
          RCut(iLandUseType) = 500.0
        Else If (Flow%T0 < 272.0) Then
          RCut(iLandUseType) = 200.0
        Else
          ! If Canopy water > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
          If (Plant%CanopyWater(iLandUseType) > 0.25) Then
            RCut(iLandUseType) = 1.0
          Else
            If (RelHumidity < 81.3) Then
              RCut(iLandUseType) = 2.5E4 * Exp(-6.93 * RelHumidity / 100.0)
            Else
              RCut(iLandUseType) = 5.8E11 * Exp(-27.8 * RelHumidity / 100.0)
            End If
          End If
        End If

        ! Surface Resistance
        Rc(iLandUseType) = 1.0 / (1.0 / RStom(iLandUseType) + 1.0 / RCut(iLandUseType) &
                         + 1.0 / SO2Standard(iLandUseType))
      End If
    End Do

    Do iLandUseType = 6, 9
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then
        If (iLandUseType == 9) Then
          Rc(iLandUseType) = 1.0 / (so2dd_ice(1) + so2dd_ice(2) * Exp(so2dd_ice(3) * (Flow%T0 - 273.15)))
        Else
          Rc(iLandUseType) = SO2Standard(iLandUseType)
        End If
      End If
    End Do


! PAN (Peroxyacetyl nitrate CH3COO2NO2)
  Else If (Species%Name(1:3) .CIEq. 'PAN') Then

    Do iLandUseType = 1, 5
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then

        ! Stomatal Resistance
        If (Plant%StomataConduct(iLandUseType) > glmin) Then
          RStom(iLandUseType) = 2.6 / Plant%StomataConduct(iLandUseType)
        End If

        ! Adjust standard resistance for tundra
        If ((iLandUseType == 5) .and. (YLatLong > 60.0)) Then
          PANStandard(iLandUseType) = 1100.0
        End If

        ! Surface Resistance
        Rc(iLandUseType) = 1.0 / (1.0 / RStom(iLandUseType) + 1.0 / PANStandard(iLandUseType))
      End If
    End Do

    Do iLandUseType = 6, 9
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then
        ! Adjust standard resistance for tundra
        If ((iLandUseType == 8) .and. (YLatLong > 60.0)) Then
          PANStandard(iLandUseType) = 1100.0
        End If
        Rc(iLandUseType) = PANStandard(iLandUseType)
      End If
    End Do

! CH4
  Else If (Species%Name(1:7) .CIEq. 'CH4') Then

    ! microbes active
    If ((RelHumidity > 40.0) .and. (Flow%T0 > 278.0)) Then
      Do iLandUseType = 1, 5
        If (Surface%LandUseFracs(iLandUseType) > 0.0) Then
          If ((iLandUseType == 5) .and. (YLatLong > 60.0)) Then
            RcTemp = Max(3600.0 * (CH4dd_tun(4) + Flow%T0 * (CH4dd_tun(3) + Flow%T0 *(CH4dd_tun(2) &
                      + CH4dd_tun(1) * Flow%T0))),0.0)
          Else
            If (Surface%SoilMoisture < 16.0) Then
              RcTemp = CH4UptakeFlux(iLandUseType) * Surface%SoilMoisture / 16.0
            Else If (Surface%SoilMoisture > 30.0) Then
              RcTemp = CH4UptakeFlux(iLandUseType) * Max(0.0 ,50.0 - Surface%SoilMoisture) / 20.0
            Else
              RcTemp = CH4UptakeFlux(iLandUseType)
            End If
          End If

          ! Convert CH4 uptake fluxes (ug m-2 h-1) to resistance (s m-1)
          If (RcTemp > 0.0) Then
            Rc(iLandUseType) = Flow%PS * mml / (RUGC * Flow%T0 * RcTemp)
          End If
        End If
      End Do

      ! bare soil
      iLandUseType = 8
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then
        If (YLatLong > 60.0) Then
          RcTemp = Max(3600.0 * (CH4dd_tun(4) + Flow%T0 * (CH4dd_tun(3) + Flow%T0 *(CH4dd_tun(2) &
                    + CH4dd_tun(1) * Flow%T0))),0.0)
        Else
          If (Surface%SoilMoisture < 16.0) Then
            RcTemp = CH4UptakeFlux(8) * Surface%SoilMoisture / 16.0
          Else If (Surface%SoilMoisture > 30.0) Then
            RcTemp = CH4UptakeFlux(8) * Max(0.0 ,50.0 - Surface%SoilMoisture) / 20.0
          Else
            RcTemp = CH4UptakeFlux(8)
          End If
        End If

        ! Convert CH4 uptake fluxes (ug m-2 h-1) to resistance (s m-1)
        If (RcTemp > 0.0) Then
          Rc(iLandUseType) = Flow%PS * mml / (RUGC * Flow%T0 * RcTemp)
        End If
      End If
    End If


! NH3
  Else If (Species%Name(1:7) .CIEq. 'AMMONIA') Then

    Do iLandUseType = 1, 5
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then

        ! Stomatal Resistance
        If (Plant%StomataConduct(iLandUseType) > glmin) Then
          RStom(iLandUseType) = 0.97 / Plant%StomataConduct(iLandUseType)
        End If

        ! Cuticular Resistance
        If (Flow%T0 < 268.0) Then
          RCut(iLandUseType) = 1000.0
        Else If (Flow%T0 < 272.0) Then
          RCut(iLandUseType) = 200.0
        Else
          ! If Canopy water > 0.25 mm assume cuticle is wet (Wyers/Erisman Atm Env 1998)
          If (Plant%CanopyWater(iLandUseType) > 0.25) Then
            RCut(iLandUseType) = 10.0
          Else
            RCut(iLandUseType) = 5.0 * Log10(Flow%T0 -271.15) * Exp ((100.0 - RelHumidity) / 12.0)
          End If
        End If

        ! Surface Resistance
        Rc(iLandUseType) = 1.0 / (1.0 / RStom(iLandUseType) + 1.0 / RCut(iLandUseType) &
                         + 1.0 / HNO3Standard(iLandUseType))

      End If
    End Do

    Do iLandUseType = 6, 9
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then
        Rc(iLandUseType) = HNO3Standard(iLandUseType)
      End If
    End Do


! Hydrogen
  Else If (Species%Name(1:8) .CIEq. 'HYDROGEN') Then

    ! microbes active
    If ((RelHumidity > 40.0) .and. (Flow%T0 > 278.0)) Then
      Do iLandUseType = 1, 5
        If (Surface%LandUseFracs(iLandUseType) > 0.0) Then

          ! Limit soil moisture to avoid excessively high deposition velocities
          !quadratic-log dependence on soil moisture for C4 grasslands
          If (iLandUseType == 4) Then
            RcTemp = (H2dd_c(4) + Log(Max(0.1, Surface%SoilMoisture / 100.0)) * (H2dd_m(4) + H2dd_q * &
                      Log(Max(0.1, Surface%SoilMoisture / 100.0)))) * 1.0E-4
            Rc(iLandUseType) = 1.0 / Min(0.00131, RcTemp)

          ! Adjust surface resistance for tundra
          Else If ((iLandUseType == 5) .and. (YLatLong > 60.0)) Then
            Rc(iLandUseType) = 3850.0

          ! Linear dependence on soil moisture
          Else
            Rc(iLandUseType) = 1.0 / (H2dd_m(iLandUseType) * Max(0.1, Surface%SoilMoisture / 100.0) &
                             + H2dd_c(iLandUseType))
          End If

        End If

      End Do

      ! Bare soil
      iLandUseType = 8
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then
        If (YLatLong > 60.0) Then
          Rc(iLandUseType) = 3850.0
        Else
          Rc(iLandUseType) = 1.0 / (H2dd_m(3) * Max(0.1, Surface%SoilMoisture / 100.0) + H2dd_c(3))
        End If
      End If
    End If

! NO2
  Else If (Species%Name(1:3) .CIEq. 'NO2') Then
    Do iLandUseType = 1, 5
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then

        ! Stomatal Resistance
        If (Plant%StomataConduct(iLandUseType) > glmin) Then
          RStom(iLandUseType) = 1.5 * 1.6 / Plant%StomataConduct(iLandUseType)
        End If

        ! Adjust standard resistance for tundra
        If ((iLandUseType == 5) .and. (YLatLong > 60.0)) Then
          NO2Standard(iLandUseType) = 1200.0
        End If

        ! Surface resistance
        Rc(iLandUseType) = 1.0 / (1.0 / RStom(iLandUseType) + 1.0 / NO2Standard(iLandUseType))
      End If

    End Do

    Do iLandUseType = 6, 9
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then
        ! Adjust standard resistance for tundra
        If ((iLandUseType == 8) .and. (YLatLong > 60.0)) Then
          NO2Standard(iLandUseType) = 1200.0
        End If
        Rc(iLandUseType) = NO2Standard(iLandUseType)
      End If
    End Do

! CO
  Else If (Species%Name(1:2) .CIEq. 'CO') Then
    Do iLandUseType = 1, 5
      If ((Surface%LandUseFracs(iLandUseType) > 0.0) .and.(RelHumidity > 40.0) .and. (Flow%T0 > 278.0)) Then
        ! Adjust standard resistance for tundra
        If ((YLatLong > 60.0) .and. (iLandUseType == 5)) Then
          COStandard(iLandUseType) = 25000.0
        End If
        Rc(iLandUseType) = COStandard(iLandUseType)
      End If
    End Do

    iLandUseType = 8
    If ((Surface%LandUseFracs(iLandUseType) > 0.0) .and.(RelHumidity > 40.0) .and. (Flow%T0 > 278.0)) Then
      ! Adjust standard resistance for tundra
      If ((YLatLong > 60.0) .and. (iLandUseType == 8)) Then
        COStandard(iLandUseType) = 25000.0
      End If
      Rc(iLandUseType) = COStandard(iLandUseType)
    End If


! HNO3
  Else If (Species%Name(1:4) .CIEq. 'HNO3') Then
    Do iLandUseType = 1, 9
      If (Surface%LandUseFracs(iLandUseType) > 0.0) Then
        If (iLandUseType == 9) Then
          Rc(iLandUseType) = Max(10.0, hno3dd_ice(3) + Flow%T0 * (hno3dd_ice(2) +     &
                          hno3dd_ice(1) * Flow%T0))
        Else
          Rc(iLandUseType) = HNO3Standard(iLandUseType)
        End If
      End If
    End Do
  End If

End Function CalcSurfResistance

!-------------------------------------------------------------------------------------------------------------

Function CalcWetScavCoeff(Species, Cloud, Rain, Z, T) Result(Lambda)


! Calculates wet deposition scavenging coefficient

  Implicit None
  ! Argument list:
  Type(Species_),   Intent(In)           :: Species
  Type(Cloud_),     Intent(In)           :: Cloud
  Type(Rain_),      Intent(In)           :: Rain
  Real(Std),        Intent(In)           :: Z         ! z position.
  Real(Std),        Intent(In)           :: T         ! Temperature (K).
    
  ! Function result:
  Real(Std) :: Lambda ! Wet deposition scavenging coefficient.
  ! Locals:
  Real(Std) :: AMixedPhase
  Real(Std) :: BMixedPhase
  Real(Std) :: DynPPt              ! zero if dynamic precip < PptCrit
  Real(Std) :: TotalOrDynCloudBase
  Real(Std) :: TotalOrDynCloudTop
  Real(Std) :: ConCloud            ! zero if no convective precipitation, min of 0.05 otherwise
  Real(Std) :: CloudBase           ! total cloud base (only defined if convective precip > PptCrit)
  Real(Std) :: CloudTop            ! total cloud top (only defined if convective precip > PptCrit)
  
  ! Temporary locals $$:
  Real(Std), Parameter :: PptCrit = 0.03 ! $$ Units? In NAME 0.03 is the default and
                                         ! can be overwritten via a namelist.

  ! Note that Cloud and Rain may be inconsistent
  ! $$ note cloud bases/tops assumed in same units as Z and increasing upwards. 
  ! Z assumed to be height above ground.

! Check total or dynamic / convective cloud base / top are defined 
! and convective cloud amount > 0 if there is rain
! If not fix up: base = 2000.0 m and top = 6000.0 m   $$ arbitrary values
! ConCloud = Max(Cloud%ConCloud, 0.05)   $$ arbitrary values
! Note that under PC2 dynamic precipitation and total cloud are out of sync
  If (Rain%DynPpt > PptCrit) Then
    DynPpt=Rain%DynPpt
    If (Cloud%TotalOrDynCloudBase < 0.0) Then
      TotalOrDynCloudBase = 2000.0
    Else
      TotalOrDynCloudBase = Cloud%TotalOrDynCloudBase
    End If
    If (Cloud%TotalOrDynCloudTop < 0.0) Then
      TotalOrDynCloudTop = 6000.0
    Else
      TotalOrDynCloudTop = Cloud%TotalOrDynCloudTop
    End If
    If (Rain%ConPpt > PptCrit) Then
      If (Cloud%ConCloudBase < 0.0) Then
        CloudBase = Min(2000.0, TotalOrDynCloudBase)
      Else
        CloudBase = Min(Cloud%ConCloudBase, TotalOrDynCloudBase)
      End If
      If (Cloud%ConCloudTop < 0.0) Then
        CloudTop = Max(6000.0, TotalOrDynCloudTop)
      Else
        CloudTop = Max(Cloud%ConCloudTop, TotalOrDynCloudTop)
      End If
    End If
  Else
    DynPpt=0.0
    If (Rain%ConPpt > PptCrit) Then
      If (Cloud%ConCloudBase < 0.0) Then
        CloudBase = 2000.0
      Else
        CloudBase = Cloud%ConCloudBase
      End If
      If (Cloud%ConCloudTop < 0.0) Then
        CloudTop = 6000.0
      Else
        CloudTop = Cloud%ConCloudTop
      End If
    End If
  End If
      
  If (Rain%ConPpt > PptCrit) Then        ! $$ should we make this ConPpt > 0 
    ConCloud = Max(Cloud%ConCloud, 0.05) 
  Else
    ConCloud = 0.0
  End If

! Initialise scavenging coefficient

  Lambda = 0.0
                 
! Significant convective ppt within regions of convective ppt
  If (ConCloud > 0.0) Then 

! Below cloud
    If (Z < CloudBase) Then
                                    
! Rain
      If (T > 270.0) Then

        If (Species%ABelowRain > 0.0) Then
          Lambda = ConCloud * Species%ABelowRain * (DynPpt + Rain%ConPpt/ConCloud)**Species%BBelowRain
        End If
! Snow
      Else

        If (Species%ABelowSnow > 0.0) Then
          Lambda = ConCloud * Species%ABelowSnow * (DynPpt + Rain%ConPpt/ConCloud)**Species%BBelowSnow
        End If
              
      End If
            
! In cloud
    Else If (Z < CloudTop) Then

! Rain
      If (T > 273.15) Then

        If (Species%AInRain > 0.0) Then
          Lambda = ConCloud * Species%AInRain * (DynPpt + Rain%ConPpt/ConCloud)**Species%BInRain
        End If
! Snow
      Else If (T < 238.15) Then

        If (Species%AInSnow > 0.0) Then
          Lambda = ConCloud * Species%AInSnow * (DynPpt + Rain%ConPpt/ConCloud)**Species%BInSnow
        End If
              
! Mixed Phase
      Else

        If (Species%AInSnow > 0.0) Then
          AMixedPhase = (Species%AInSnow * (273.15 - T) +       &
            Species%AInRain * (T - 238.15)) / (273.15 - 238.15) 
          BMixedPhase = (Species%BInSnow * (273.15 - T) +       &
            Species%BInRain * (T - 238.15)) / (273.15 - 238.15) 
          Lambda = ConCloud * AMixedPhase * (DynPpt + Rain%ConPpt/ConCloud)**BMixedPhase
        End If
        
      End If
    End If
  End If
        
! Dynamic ppt outside of regions of convective ppt
  If (DynPpt > 0.0) Then
! Below cloud
    If (Z < TotalOrDynCloudBase) Then
            
! Rain
      If (T > 270.0) Then

        If (Species%ABelowRain > 0.0) Then
          Lambda = Lambda + (1.0 - ConCloud) * Species%ABelowRain * DynPpt**Species%BBelowRain
        End If
! Snow
      Else

        If (Species%ABelowSnow > 0.0) Then
          Lambda = Lambda + (1.0 - ConCloud) * Species%ABelowSnow * DynPpt**Species%BBelowSnow
        End If
              
      End If

! In cloud
    Else If (Z < TotalOrDynCloudTop) Then
            
! Rain
      If (T > 273.15) Then

        If (Species%AInRain > 0.0) Then
          Lambda = Lambda + (1.0 - ConCloud) * Species%AInRain * DynPpt**Species%BInRain
        End If
! Snow
      Else If (T < 238.15) Then

        If (Species%AInSnow > 0.0) Then
          Lambda = Lambda + (1.0 - ConCloud) * Species%AInSnow * DynPpt**Species%BInSnow
        End If
              
! Mixed Phase
      Else

        If (Species%AInSnow > 0.0) Then
          AMixedPhase = (Species%AInSnow * (273.15 - T) +       &
            Species%AInRain * (T - 238.15)) / (273.15 - 238.15) 
          BMixedPhase = (Species%BInSnow * (273.15 - T) +       &
            Species%BInRain * (T - 238.15)) / (273.15 - 238.15) 
          Lambda = Lambda + (1.0 - ConCloud) * AMixedPhase * DynPpt**BMixedPhase
        End If

      End If
    End If
  End If


End Function CalcWetScavCoeff

!-------------------------------------------------------------------------------------------------------------

End Module SpeciesModule
