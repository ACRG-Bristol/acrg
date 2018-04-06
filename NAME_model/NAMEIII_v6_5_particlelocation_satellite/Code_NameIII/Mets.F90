! Module: Mets Module

Module MetsModule

! This module provides code to make all met modules look identical.

! Module overview
! ---------------

! $$

! Certain pieces of code need to be repeated for some or all met modules. This code repetition is achieved by
! prepreocessing Mets.P90 with Preprocessor.exe to produce Mets.F90. Whenever a new met module is added,
! appropriate $UseBlock directives should be added for use by Preprocessor.exe.

! Module use
! ----------

! $$

! Module call tree
! ----------------

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use CommonMetModule

Use PrototypeMetModule
Use SingleSiteMetModule
Use NWPMetModule
Use RadarMetModule
Use AncillaryMetModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: CommonMetP_          ! A pointer to an instance of the type CommonMet_.
Public  :: Mets_                ! A collection of met module instance states.
Public  :: InitMets             ! Initialises a collection of met module instance
                                ! states.

Public  :: AddPrototypeMet           ! Adds the state of a Prototype met module instance to a
                                ! collection of met module instance states.
Public  :: AddSingleSiteMet           ! Adds the state of a SingleSite met module instance to a
                                ! collection of met module instance states.
Public  :: AddNWPMet           ! Adds the state of a NWP met module instance to a
                                ! collection of met module instance states.
Public  :: AddRadarMet           ! Adds the state of a Radar met module instance to a
                                ! collection of met module instance states.
Public  :: AddAncillaryMet           ! Adds the state of a Ancillary met module instance to a
                                ! collection of met module instance states.

Public  :: FindMetIndex         ! Finds the indices of a met module and a met module
                                ! instance.
Public  :: SetUpCoordsEtc_Mets  ! Sets up Coords and Grids by adding any extra coords
                                ! and grids which Mets wants to define.
Public  :: SetUpMets_CoordsEtc  ! Sets up Mets using information from EtaDefns, Coords
                                ! and Grids.
Public  :: SetUpMets_iCoordsEtc ! Sets up indices in the met module instances for
                                ! referring coord systems and grids.
Public  :: MetList              ! Returns a list of the met module instances.
Public  :: MetValid             ! Returns information on the validity of a met module
                                ! instance.
Public  :: UpdateMets           ! Updates the state of the met module instances.
Public  :: UpdateMet            ! Updates the state of a met module instance.
Public  :: ResetMets            ! Resets the state of the met module instances for a
                                ! new realisation.

!-------------------------------------------------------------------------------------------------------------

! The following items are items from the various met modules which need to be made
! available by the mets module.

Public  :: PrototypeMet_    ! Information describing the state of an instance of the Prototype
                       ! met module.
Public  :: SingleSiteMet_    ! Information describing the state of an instance of the SingleSite
                       ! met module.
Public  :: NWPMet_    ! Information describing the state of an instance of the NWP
                       ! met module.
Public  :: RadarMet_    ! Information describing the state of an instance of the Radar
                       ! met module.
Public  :: AncillaryMet_    ! Information describing the state of an instance of the Ancillary
                       ! met module.

Public  :: InitPrototypeMet ! Initialises an instance of PrototypeMet_.
Public  :: InitSingleSiteMet ! Initialises an instance of SingleSiteMet_.
Public  :: InitNWPMet ! Initialises an instance of NWPMet_.
Public  :: InitRadarMet ! Initialises an instance of RadarMet_.
Public  :: InitAncillaryMet ! Initialises an instance of AncillaryMet_.

!-------------------------------------------------------------------------------------------------------------

Type :: CommonMetP_ ! A pointer to an instance of the type CommonMet_.
  Type(CommonMet_), Pointer :: P ! Pointer to an instance of the type CommonMet_.
End Type CommonMetP_

!-------------------------------------------------------------------------------------------------------------

Type :: Mets_ ! A collection of met module instance states.
  Integer           :: nMetMods
  Integer           :: nMets(MaxMetMods)
  Type(CommonMetP_) :: C(MaxMetMods, MaxMetsPerMod)

  Integer        :: nPrototypeMets
  Type(PrototypeMet_) :: PrototypeMets(MaxPrototypeMets)
  Integer        :: nSingleSiteMets
  Type(SingleSiteMet_) :: SingleSiteMets(MaxSingleSiteMets)
  Integer        :: nNWPMets
  Type(NWPMet_) :: NWPMets(MaxNWPMets)
  Integer        :: nRadarMets
  Type(RadarMet_) :: RadarMets(MaxRadarMets)
  Integer        :: nAncillaryMets
  Type(AncillaryMet_) :: AncillaryMets(MaxAncillaryMets)

  ! nMetMods  :: Number of met modules.
  ! nMets     :: Number of instances of each met module.
  ! C         :: Pointers to the common parts of the met module instance states.
  ! n????Mets :: Number of instances of the ???? met module.
  ! ????Mets  :: States of all instances of the ???? met module.
End Type Mets_

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitMets() Result(Mets)
! Initialises a collection of met module instance states.

  Implicit None
  ! Function result:
  Type(Mets_) :: Mets ! The initialised collection of met module instance states.

  Mets%nMetMods = 0

  Mets%nPrototypeMets = 0
  Mets%nSingleSiteMets = 0
  Mets%nNWPMets = 0
  Mets%nRadarMets = 0
  Mets%nAncillaryMets = 0

End Function InitMets

!-------------------------------------------------------------------------------------------------------------


Subroutine AddPrototypeMet(PrototypeMet, Mets)
! Adds the state of a Prototype met module instance to a collection of met module instance
! states.

  Implicit None
  ! Argument list:
  Type(PrototypeMet_), Intent(In)            :: PrototypeMet ! The state of the Prototype met module
                                                   ! instance.
  Type(Mets_),    Intent(InOut), Target :: Mets    ! The collection of met module
                                                   ! instance states.
  ! Locals:
  Integer :: i       ! Loop index.
  Integer :: j       ! Loop index.
  Integer :: iMetMod ! Index of met module to be added.
  Integer :: iMet    ! Index of met module instance to be added.

  ! Calculate iMetMod and, if necessary, check for too many met modules, update
  ! Mets%nMetMods and initialise Mets%nMets(iMetMod).
  If (Mets%nPrototypeMets == 0) Then
    If (Mets%nMetMods == MaxMetMods) Then
      Call Message('FATAL ERROR in AddPrototypeMet: Too many met modules', 3)
    End If
    iMetMod             = Mets%nMetMods + 1
    Mets%nMetMods       = iMetMod
    Mets%nMets(iMetMod) = 0
  Else
    iMetMod = Mets%PrototypeMets(1)%C%iMetMod
  End If

  ! Check for too many met module instances of this type and for duplicate names.
  If (Mets%nMets(iMetMod) == MaxPrototypeMets) Then
    Call Message(                                          &
           'FATAL ERROR in AddPrototypeMet: '                // &
           'Too many instances of the Prototype met module',    &
           3                                               &
         )
  End If
  Do i = 1, Mets%nMetMods
  Do j = 1, Mets%nMets(i)
    If (Mets%C(i ,j)%P%MetName .CIEq. PrototypeMet%C%MetName) Then
      Call Message(                                                               &
             'FATAL ERROR in AddPrototypeMet: a met module instance with the same ' // &
             'name already exists',                                               &
             3                                                                    &
           )
    End If
  End Do
  End Do

  ! Calculate iMet.
  iMet = Mets%nMets(iMetMod) + 1

  ! Add met module instance.
  Mets%nMets(iMetMod)           =  iMet
  Mets%nPrototypeMets                =  iMet
  Mets%PrototypeMets(iMet)           =  PrototypeMet
  Mets%PrototypeMets(iMet)%C%iMetMod =  iMetMod
  Mets%PrototypeMets(iMet)%C%iMet    =  iMet
  Mets%C(iMetMod, iMet)%P       => Mets%PrototypeMets(iMet)%C

  ! Check TValid is the right type of time.
  If (IsTimeInterval(Mets%C(iMetMod, iMet)%P%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in AddPrototypeMet', 4)
  End If

  ! Check consistency of FixedMet between met module instances.
  If (Mets%C(1, 1)%P%FixedMet .neqv. Mets%C(iMetMod, iMet)%P%FixedMet) Then
    Call Message('UNEXPECTED FATAL ERROR in AddPrototypeFlow', 4)
  End If

End Subroutine AddPrototypeMet


Subroutine AddSingleSiteMet(SingleSiteMet, Mets)
! Adds the state of a SingleSite met module instance to a collection of met module instance
! states.

  Implicit None
  ! Argument list:
  Type(SingleSiteMet_), Intent(In)            :: SingleSiteMet ! The state of the SingleSite met module
                                                   ! instance.
  Type(Mets_),    Intent(InOut), Target :: Mets    ! The collection of met module
                                                   ! instance states.
  ! Locals:
  Integer :: i       ! Loop index.
  Integer :: j       ! Loop index.
  Integer :: iMetMod ! Index of met module to be added.
  Integer :: iMet    ! Index of met module instance to be added.

  ! Calculate iMetMod and, if necessary, check for too many met modules, update
  ! Mets%nMetMods and initialise Mets%nMets(iMetMod).
  If (Mets%nSingleSiteMets == 0) Then
    If (Mets%nMetMods == MaxMetMods) Then
      Call Message('FATAL ERROR in AddSingleSiteMet: Too many met modules', 3)
    End If
    iMetMod             = Mets%nMetMods + 1
    Mets%nMetMods       = iMetMod
    Mets%nMets(iMetMod) = 0
  Else
    iMetMod = Mets%SingleSiteMets(1)%C%iMetMod
  End If

  ! Check for too many met module instances of this type and for duplicate names.
  If (Mets%nMets(iMetMod) == MaxSingleSiteMets) Then
    Call Message(                                          &
           'FATAL ERROR in AddSingleSiteMet: '                // &
           'Too many instances of the SingleSite met module',    &
           3                                               &
         )
  End If
  Do i = 1, Mets%nMetMods
  Do j = 1, Mets%nMets(i)
    If (Mets%C(i ,j)%P%MetName .CIEq. SingleSiteMet%C%MetName) Then
      Call Message(                                                               &
             'FATAL ERROR in AddSingleSiteMet: a met module instance with the same ' // &
             'name already exists',                                               &
             3                                                                    &
           )
    End If
  End Do
  End Do

  ! Calculate iMet.
  iMet = Mets%nMets(iMetMod) + 1

  ! Add met module instance.
  Mets%nMets(iMetMod)           =  iMet
  Mets%nSingleSiteMets                =  iMet
  Mets%SingleSiteMets(iMet)           =  SingleSiteMet
  Mets%SingleSiteMets(iMet)%C%iMetMod =  iMetMod
  Mets%SingleSiteMets(iMet)%C%iMet    =  iMet
  Mets%C(iMetMod, iMet)%P       => Mets%SingleSiteMets(iMet)%C

  ! Check TValid is the right type of time.
  If (IsTimeInterval(Mets%C(iMetMod, iMet)%P%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in AddSingleSiteMet', 4)
  End If

  ! Check consistency of FixedMet between met module instances.
  If (Mets%C(1, 1)%P%FixedMet .neqv. Mets%C(iMetMod, iMet)%P%FixedMet) Then
    Call Message('UNEXPECTED FATAL ERROR in AddSingleSiteFlow', 4)
  End If

End Subroutine AddSingleSiteMet


Subroutine AddNWPMet(NWPMet, Mets)
! Adds the state of a NWP met module instance to a collection of met module instance
! states.

  Implicit None
  ! Argument list:
  Type(NWPMet_), Intent(In)            :: NWPMet ! The state of the NWP met module
                                                   ! instance.
  Type(Mets_),    Intent(InOut), Target :: Mets    ! The collection of met module
                                                   ! instance states.
  ! Locals:
  Integer :: i       ! Loop index.
  Integer :: j       ! Loop index.
  Integer :: iMetMod ! Index of met module to be added.
  Integer :: iMet    ! Index of met module instance to be added.

  ! Calculate iMetMod and, if necessary, check for too many met modules, update
  ! Mets%nMetMods and initialise Mets%nMets(iMetMod).
  If (Mets%nNWPMets == 0) Then
    If (Mets%nMetMods == MaxMetMods) Then
      Call Message('FATAL ERROR in AddNWPMet: Too many met modules', 3)
    End If
    iMetMod             = Mets%nMetMods + 1
    Mets%nMetMods       = iMetMod
    Mets%nMets(iMetMod) = 0
  Else
    iMetMod = Mets%NWPMets(1)%C%iMetMod
  End If

  ! Check for too many met module instances of this type and for duplicate names.
  If (Mets%nMets(iMetMod) == MaxNWPMets) Then
    Call Message(                                          &
           'FATAL ERROR in AddNWPMet: '                // &
           'Too many instances of the NWP met module',    &
           3                                               &
         )
  End If
  Do i = 1, Mets%nMetMods
  Do j = 1, Mets%nMets(i)
    If (Mets%C(i ,j)%P%MetName .CIEq. NWPMet%C%MetName) Then
      Call Message(                                                               &
             'FATAL ERROR in AddNWPMet: a met module instance with the same ' // &
             'name already exists',                                               &
             3                                                                    &
           )
    End If
  End Do
  End Do

  ! Calculate iMet.
  iMet = Mets%nMets(iMetMod) + 1

  ! Add met module instance.
  Mets%nMets(iMetMod)           =  iMet
  Mets%nNWPMets                =  iMet
  Mets%NWPMets(iMet)           =  NWPMet
  Mets%NWPMets(iMet)%C%iMetMod =  iMetMod
  Mets%NWPMets(iMet)%C%iMet    =  iMet
  Mets%C(iMetMod, iMet)%P       => Mets%NWPMets(iMet)%C

  ! Check TValid is the right type of time.
  If (IsTimeInterval(Mets%C(iMetMod, iMet)%P%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in AddNWPMet', 4)
  End If

  ! Check consistency of FixedMet between met module instances.
  If (Mets%C(1, 1)%P%FixedMet .neqv. Mets%C(iMetMod, iMet)%P%FixedMet) Then
    Call Message('UNEXPECTED FATAL ERROR in AddNWPFlow', 4)
  End If

End Subroutine AddNWPMet


Subroutine AddRadarMet(RadarMet, Mets)
! Adds the state of a Radar met module instance to a collection of met module instance
! states.

  Implicit None
  ! Argument list:
  Type(RadarMet_), Intent(In)            :: RadarMet ! The state of the Radar met module
                                                   ! instance.
  Type(Mets_),    Intent(InOut), Target :: Mets    ! The collection of met module
                                                   ! instance states.
  ! Locals:
  Integer :: i       ! Loop index.
  Integer :: j       ! Loop index.
  Integer :: iMetMod ! Index of met module to be added.
  Integer :: iMet    ! Index of met module instance to be added.

  ! Calculate iMetMod and, if necessary, check for too many met modules, update
  ! Mets%nMetMods and initialise Mets%nMets(iMetMod).
  If (Mets%nRadarMets == 0) Then
    If (Mets%nMetMods == MaxMetMods) Then
      Call Message('FATAL ERROR in AddRadarMet: Too many met modules', 3)
    End If
    iMetMod             = Mets%nMetMods + 1
    Mets%nMetMods       = iMetMod
    Mets%nMets(iMetMod) = 0
  Else
    iMetMod = Mets%RadarMets(1)%C%iMetMod
  End If

  ! Check for too many met module instances of this type and for duplicate names.
  If (Mets%nMets(iMetMod) == MaxRadarMets) Then
    Call Message(                                          &
           'FATAL ERROR in AddRadarMet: '                // &
           'Too many instances of the Radar met module',    &
           3                                               &
         )
  End If
  Do i = 1, Mets%nMetMods
  Do j = 1, Mets%nMets(i)
    If (Mets%C(i ,j)%P%MetName .CIEq. RadarMet%C%MetName) Then
      Call Message(                                                               &
             'FATAL ERROR in AddRadarMet: a met module instance with the same ' // &
             'name already exists',                                               &
             3                                                                    &
           )
    End If
  End Do
  End Do

  ! Calculate iMet.
  iMet = Mets%nMets(iMetMod) + 1

  ! Add met module instance.
  Mets%nMets(iMetMod)           =  iMet
  Mets%nRadarMets                =  iMet
  Mets%RadarMets(iMet)           =  RadarMet
  Mets%RadarMets(iMet)%C%iMetMod =  iMetMod
  Mets%RadarMets(iMet)%C%iMet    =  iMet
  Mets%C(iMetMod, iMet)%P       => Mets%RadarMets(iMet)%C

  ! Check TValid is the right type of time.
  If (IsTimeInterval(Mets%C(iMetMod, iMet)%P%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in AddRadarMet', 4)
  End If

  ! Check consistency of FixedMet between met module instances.
  If (Mets%C(1, 1)%P%FixedMet .neqv. Mets%C(iMetMod, iMet)%P%FixedMet) Then
    Call Message('UNEXPECTED FATAL ERROR in AddRadarFlow', 4)
  End If

End Subroutine AddRadarMet


Subroutine AddAncillaryMet(AncillaryMet, Mets)
! Adds the state of a Ancillary met module instance to a collection of met module instance
! states.

  Implicit None
  ! Argument list:
  Type(AncillaryMet_), Intent(In)            :: AncillaryMet ! The state of the Ancillary met module
                                                   ! instance.
  Type(Mets_),    Intent(InOut), Target :: Mets    ! The collection of met module
                                                   ! instance states.
  ! Locals:
  Integer :: i       ! Loop index.
  Integer :: j       ! Loop index.
  Integer :: iMetMod ! Index of met module to be added.
  Integer :: iMet    ! Index of met module instance to be added.

  ! Calculate iMetMod and, if necessary, check for too many met modules, update
  ! Mets%nMetMods and initialise Mets%nMets(iMetMod).
  If (Mets%nAncillaryMets == 0) Then
    If (Mets%nMetMods == MaxMetMods) Then
      Call Message('FATAL ERROR in AddAncillaryMet: Too many met modules', 3)
    End If
    iMetMod             = Mets%nMetMods + 1
    Mets%nMetMods       = iMetMod
    Mets%nMets(iMetMod) = 0
  Else
    iMetMod = Mets%AncillaryMets(1)%C%iMetMod
  End If

  ! Check for too many met module instances of this type and for duplicate names.
  If (Mets%nMets(iMetMod) == MaxAncillaryMets) Then
    Call Message(                                          &
           'FATAL ERROR in AddAncillaryMet: '                // &
           'Too many instances of the Ancillary met module',    &
           3                                               &
         )
  End If
  Do i = 1, Mets%nMetMods
  Do j = 1, Mets%nMets(i)
    If (Mets%C(i ,j)%P%MetName .CIEq. AncillaryMet%C%MetName) Then
      Call Message(                                                               &
             'FATAL ERROR in AddAncillaryMet: a met module instance with the same ' // &
             'name already exists',                                               &
             3                                                                    &
           )
    End If
  End Do
  End Do

  ! Calculate iMet.
  iMet = Mets%nMets(iMetMod) + 1

  ! Add met module instance.
  Mets%nMets(iMetMod)           =  iMet
  Mets%nAncillaryMets                =  iMet
  Mets%AncillaryMets(iMet)           =  AncillaryMet
  Mets%AncillaryMets(iMet)%C%iMetMod =  iMetMod
  Mets%AncillaryMets(iMet)%C%iMet    =  iMet
  Mets%C(iMetMod, iMet)%P       => Mets%AncillaryMets(iMet)%C

  ! Check TValid is the right type of time.
  If (IsTimeInterval(Mets%C(iMetMod, iMet)%P%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in AddAncillaryMet', 4)
  End If

  ! Check consistency of FixedMet between met module instances.
  If (Mets%C(1, 1)%P%FixedMet .neqv. Mets%C(iMetMod, iMet)%P%FixedMet) Then
    Call Message('UNEXPECTED FATAL ERROR in AddAncillaryFlow', 4)
  End If

End Subroutine AddAncillaryMet


!-------------------------------------------------------------------------------------------------------------

Subroutine FindMetIndex(MetModName, MetName, Mets, iMetMod, iMet)
! Finds the indices of a met module and a met module instance.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)  :: MetModName ! Met module name.
  Character(*), Intent(In)  :: MetName    ! Met module instance name.
  Type(Mets_),  Intent(In)  :: Mets       ! Collection of met module instance states.
  Integer,      Intent(Out) :: iMetMod    ! Met module index.
  Integer,      Intent(Out) :: iMet       ! Met module instance index.
  ! Locals:
  Integer :: i ! Loop index.

  iMetMod = 0
  Do i = 1, Mets%nMetMods
    If (Mets%nMets(i) > 0) Then
      If (MetModName .CIEq. Mets%C(i, 1)%P%MetModName) Then
        iMetMod = i
        Exit
      End If
    End If
  End Do
  If (iMetMod == 0) Then
    Call Message('FATAL ERROR in FindMetIndex: met module not found', 3)
  End If

  iMet = 0
  Do i = 1, Mets%nMets(iMetMod)
    If (MetName .CIEq. Mets%C(iMetMod, i)%P%MetName) Then
      iMet = i
      Exit
    End If
  End Do
  If (iMet == 0) Then
    Call Message('FATAL ERROR in FindMetIndex: met module instance not found', 3)
  End If

End Subroutine FindMetIndex

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpCoordsEtc_Mets(Mets, Coords, Grids)
! Sets up Coords and Grids by adding any extra coords and grids which Mets wants to
! define.

  Implicit None
  ! Argument list:
  Type(Mets_),   Intent(In)    :: Mets   ! Collection of met module instance states.
  Type(Coords_), Intent(InOut) :: Coords ! Collection of coord systems.
  Type(Grids_),  Intent(InOut) :: Grids  ! Collection of grids.
  ! Locals:
  Integer :: iMetMod ! Met module index.
  Integer :: iMet    ! Met module instance index.

  Do iMetMod = 1, Mets%nMetMods
  Do iMet    = 1, Mets%nMets(iMetMod)

    If (.false.) Then

    Else If (                               &
      Mets%nSingleSiteMets /= 0 .and.             &
      iMetMod == Mets%SingleSiteMets(1)%C%iMetMod &
    ) Then

      Call SetUpCoordsEtc_SingleSiteMet(Mets%SingleSiteMets(iMet), Coords, Grids)

    Else If (                               &
      Mets%nNWPMets /= 0 .and.             &
      iMetMod == Mets%NWPMets(1)%C%iMetMod &
    ) Then

      Call SetUpCoordsEtc_NWPMet(Mets%NWPMets(iMet), Coords, Grids)

    Else If (                               &
      Mets%nAncillaryMets /= 0 .and.             &
      iMetMod == Mets%AncillaryMets(1)%C%iMetMod &
    ) Then

      Call SetUpCoordsEtc_AncillaryMet(Mets%AncillaryMets(iMet), Coords, Grids)


    ! Check met module present.
    Else If (.not. MetPresent(iMetMod, iMet, Mets)) Then
      Call Message('UNEXPECTED FATAL ERROR in SetUpCoordsEtc_Mets: met module not found', 4)
    End If

  End Do
  End Do

End Subroutine SetUpCoordsEtc_Mets

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpMets_CoordsEtc(EtaDefns, Coords, Grids, MetEnsembleSize, Mets)
! Sets up Mets using information from EtaDefns, Coords and Grids.

  Implicit None
  ! Argument list:
  Type(EtaDefns_), Intent(In)    :: EtaDefns        ! Collection of eta definitions.
  Type(Coords_),   Intent(In)    :: Coords          ! Collection of coord systems.
  Type(Grids_),    Intent(In)    :: Grids           ! Collection of grids.
  Integer,         Intent(In)    :: MetEnsembleSize ! Size of the met ensemble (i.e. number of met
                                                    ! realisations).
  Type(Mets_),     Intent(InOut) :: Mets            ! Collection of met module instance states.
  ! Locals:
  Integer :: iMetMod ! Met module index.
  Integer :: iMet    ! Met module instance index.

  Do iMetMod = 1, Mets%nMetMods
  Do iMet    = 1, Mets%nMets(iMetMod)

    If (.false.) Then

    Else If (                               &
      Mets%nSingleSiteMets /= 0 .and.             &
      iMetMod == Mets%SingleSiteMets(1)%C%iMetMod &
    ) Then

      Call SetUpSingleSiteMet_CoordsEtc(EtaDefns, Coords, Grids, MetEnsembleSize, Mets%SingleSiteMets(iMet))

    Else If (                               &
      Mets%nNWPMets /= 0 .and.             &
      iMetMod == Mets%NWPMets(1)%C%iMetMod &
    ) Then

      Call SetUpNWPMet_CoordsEtc(EtaDefns, Coords, Grids, MetEnsembleSize, Mets%NWPMets(iMet))

    Else If (                               &
      Mets%nRadarMets /= 0 .and.             &
      iMetMod == Mets%RadarMets(1)%C%iMetMod &
    ) Then

      Call SetUpRadarMet_CoordsEtc(EtaDefns, Coords, Grids, MetEnsembleSize, Mets%RadarMets(iMet))

    Else If (                               &
      Mets%nAncillaryMets /= 0 .and.             &
      iMetMod == Mets%AncillaryMets(1)%C%iMetMod &
    ) Then

      Call SetUpAncillaryMet_CoordsEtc(EtaDefns, Coords, Grids, MetEnsembleSize, Mets%AncillaryMets(iMet))


    ! Check met module present.
    Else If (.not. MetPresent(iMetMod, iMet, Mets)) Then
      Call Message('UNEXPECTED FATAL ERROR in SetUpMets_CoordsEtc: met module not found', 4)
    End If

  End Do
  End Do

End Subroutine SetUpMets_CoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpMets_iCoordsEtc(Coords, Grids, Mets)
! Sets up indices in the met module instances for referring coord systems and grids.

  Implicit None
  ! Argument list:
  Type(Coords_), Intent(In)    :: Coords ! Collection of coord systems.
  Type(Grids_),  Intent(In)    :: Grids  ! Collection of grids.
  Type(Mets_),   Intent(InOut) :: Mets   ! Collection of met module instance states.
  ! Locals:
  Type(CommonMet_), Pointer :: CommonMet ! Abreviation for the part of the met state common to all met
                                         ! modules.
  Integer                   :: iMetMod   ! Met module index.
  Integer                   :: iMet      ! Met module instance index.
  Integer                   :: i         ! Loop Index.

  Do iMetMod = 1, Mets%nMetMods
  Do iMet    = 1, Mets%nMets(iMetMod)

    CommonMet => Mets%C(iMetMod, iMet)%P

    Do i = 1, CommonMet%nHCoords
      CommonMet%iHCoords(i) = FindHCoordIndex(CommonMet%HCoordNames(i), Coords)
    End Do

    Do i = 1, CommonMet%nZCoords
      CommonMet%iZCoords(i) = FindZCoordIndex(CommonMet%ZCoordNames(i), Coords)
    End Do

    Do i = 1, CommonMet%nHGrids
      CommonMet%iHGrids(i) = FindHGridIndex(CommonMet%HGridNames(i), Grids)
    End Do

    Do i = 1, CommonMet%nZGrids
      CommonMet%iZGrids(i) = FindZGridIndex(CommonMet%ZGridNames(i), Grids)
    End Do

  End Do
  End Do

End Subroutine SetUpMets_iCoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine MetList(Mets, nMets, MetNames)
! Returns a list of the met module instances.

  Implicit None
  ! Argument list:
  Type(Mets_),              Intent(In)  :: Mets
  Integer,                  Intent(Out) :: nMets
  Character(MaxCharLength), Intent(Out) :: MetNames(MaxMets)
  ! Mets     :: Collection of met module instance states.
  ! nMets    :: Number of the met module instance states.
  ! MetNames :: Names of the met module instance states.
  ! Locals:
  Integer :: iMetMod ! Met module index.
  Integer :: iMet    ! Met module instance index.

  nMets = 0
  Do iMetMod = 1, Mets%nMetMods
  Do iMet    = 1, Mets%nMets(iMetMod)
    nMets = nMets + 1
    MetNames(nMets) = Trim(Mets%C(iMetMod, iMet)%P%MetModName) // &
                      '.'                                      // &
                      Trim(Mets%C(iMetMod, iMet)%P%MetName)
  End Do
  End Do


End Subroutine MetList

!-------------------------------------------------------------------------------------------------------------

Function MetPresent(iMetMod, iMet, Mets)
! Checks whether a met module instance with given indices is present.

  Implicit None
  ! Argument list:
  Integer,     Intent(In) :: iMetMod ! Met module index.
  Integer,     Intent(In) :: iMet    ! Met module instance index.
  Type(Mets_), Intent(In) :: Mets    ! Collection of met module instance states.
  ! Function result:
  Logical :: MetPresent ! Indicates whether the met module instance is present.

  MetPresent = .false.

  If (.false.) Then

  Else If (                               &
    Mets%nPrototypeMets /= 0 .and.             &
    iMetMod == Mets%PrototypeMets(1)%C%iMetMod &
  ) Then

    MetPresent = .true.

  Else If (                               &
    Mets%nSingleSiteMets /= 0 .and.             &
    iMetMod == Mets%SingleSiteMets(1)%C%iMetMod &
  ) Then

    MetPresent = .true.

  Else If (                               &
    Mets%nNWPMets /= 0 .and.             &
    iMetMod == Mets%NWPMets(1)%C%iMetMod &
  ) Then

    MetPresent = .true.

  Else If (                               &
    Mets%nRadarMets /= 0 .and.             &
    iMetMod == Mets%RadarMets(1)%C%iMetMod &
  ) Then

    MetPresent = .true.

  Else If (                               &
    Mets%nAncillaryMets /= 0 .and.             &
    iMetMod == Mets%AncillaryMets(1)%C%iMetMod &
  ) Then

    MetPresent = .true.


  End If

End Function MetPresent

!-------------------------------------------------------------------------------------------------------------

Subroutine MetValid(Mets, iMetMod, iMet, Time, Valid, TValid, ValidAttribs)
! Returns information on the validity of a met module instance.

  Implicit None
  ! Argument list:
  Type(Mets_), Intent(In)            :: Mets
  Integer,     Intent(In)            :: iMetMod
  Integer,     Intent(In)            :: iMet
  Type(Time_), Intent(In)            :: Time
  Logical,     Intent(Out)           :: Valid
  Type(Time_), Intent(Out)           :: TValid
  Logical,     Intent(Out), Optional :: ValidAttribs(MaxMetAttribs)
  ! Mets         :: Collection of met module instance states.
  ! iMetMod      :: Met module index.
  ! iMet         :: Met module instance index.
  ! Time         :: Current time in non-fixed met cases and current time or the time of the met in fixed met
  !                 cases.
  ! Valid        :: Indicates the met module instance is known to be valid at the current time.
  ! TValid       :: Earliest time the validity of the met module instance might change.
  ! ValidAttribs ::
  ! Locals:
  Type(CommonMet_), Pointer :: CommonMet ! Abreviation for the part of met state
                                         ! common to all met modules.

  CommonMet => Mets%C(iMetMod, iMet)%P

# ifdef ExtraChecks
    If (                            &
      iMetMod > Mets%nMetMods .or.  &
      iMet    > Mets%nMets(iMetMod) &
    ) Then
      Call Message('UNEXPECTED FATAL ERROR in MetValid', 4)
    End If
    If (IsTimeInterval(Time)) Then
      Call Message('UNEXPECTED FATAL ERROR in MetValid', 4)
    End If
# endif

  If (Present(ValidAttribs)) Then
    ValidAttribs(:) = CommonMet%ValidAttribs(:) .and. Time < CommonMet%TValid
  End If
  Valid  = CommonMet%Valid .and. Time < CommonMet%TValid
  TValid = TMax(CommonMet%TValid, Time)

End Subroutine MetValid

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateMets(                      &
             Coords, Grids,                 &
             SameResultsWithUpdateOnDemand, &
             iCase, iMetCase,               &
             MetTime,                       &
             OverallTValid,                 &
             Mets,                          &
             Units                          &
           )
! Updates the state of the met module instances.


  Use OpenMPModule, only : LookaheadFileReadCompleteWait &
                          ,LookaheadFileReadRequest      &
                          ,OpenMPOpts_

  Implicit None

  ! Argument list:
  Type(Coords_), Intent(In)    :: Coords
  Type(Grids_),  Intent(In)    :: Grids
  Logical,       Intent(In)    :: SameResultsWithUpdateOnDemand
  Integer,       Intent(In)    :: iCase
  Integer,       Intent(In)    :: iMetCase
  Type(Time_),   Intent(In)    :: MetTime
  Type(Time_),   Intent(InOut) :: OverallTValid
  Type(Mets_),   Intent(InOut) :: Mets
  Type(Units_),  Intent(InOut) :: Units
  ! Coords                        :: Collection of coord systems.
  ! Grids                         :: Collection of grids.
  ! SameResultsWithUpdateOnDemand :: Indicates results should be made the same whether the met and flow module
  !                                  instances use update-on-demand or update-at-once (so it should really be
  !                                  called SameResultsWithOrWithoutUpdateOnDemand).
  ! iCase                         :: Number of case.
  ! iMetCase                      :: Number of the met realisation in the met ensemble.
  ! MetTime                       :: Time for which the met module instances are to be updated (this must be
  !                                  the current time unless we have a fixed met case, in which case it must
  !                                  be the time of the fixed met).
  ! OverallTValid                 :: Earliest time that the validity of any of the met module instances might
  !                                  change, assuming all the instances which have been prepared for
  !                                  update-on-demand are updated now. The value is that determined at the end
  !                                  of this routine (the actual time may be later).
  ! Mets                          :: Collection of met module instance states.
  ! Units                         :: Collection of information on input/output unit numbers.
  ! Locals:
  Type(CommonMet_), Pointer :: CommonMet
  Integer                   :: iMetMod
  Integer                   :: iMet
  Type(Time_)               :: TValid
  Logical                   :: UpdateNow
  Logical                   :: ParallelMetRead
  ! CommonMet :: Abbreviation for the common part of the met module instance.
  ! iMetMod   :: Met module index.
  ! iMet      :: Met module instance index.
  ! TValid    :: Earliest time that the validity (overall or for any single attribute) of the met module
  !              instance might change, assuming the met module instance is updated now. The value is that
  !              determined at the end of PrepareForUpdateMet (the actual time may be later).
  ! UpdateNow :: Indicates the met module instance must be updated now (even if update-on-demand is
  !              specified). If set, the value of TValid is not reliable.
  ! ParallelMetRead :: specifies whether (at least one) NWPMet module is updated by the parallel IO
  !                    thread

  ! Go through all Met modules and check whether one of them has OpenMPOpts%ParallelMetRead == .true.

  ParallelMetRead = .false.

  Do iMet = 1, Mets%nNWPMets
    ParallelMetRead = ParallelMetRead .or. Mets%NWPMets(iMet)%OpenMPOpts%ParallelMetRead
    ! Check if the ParallelMetRead flag is consistent between all NWPMet instances
    ! $$ At some point we might allow this flag to differ to only prefect selected NWPMet
    if ((iMet > 1) .and. (Mets%NWPMets(iMet)%OpenMPOpts%ParallelMetRead .neqv. ParallelMetRead)) Then
      Call Message('ERROR: The flag ParallelMetRead has to be identical for all NWPMet module instances',4)
    End If
  End Do

  ! Wait until any lookahead reads have completed before starting our Met Reads
  If (ParallelMetRead) Then
    Call LookaheadFileReadCompleteWait()
  End If

  ! Loop over met module instances.
  Do iMetMod = 1, Mets%nMetMods
  Do iMet    = 1, Mets%nMets(iMetMod)

    CommonMet => Mets%C(iMetMod, iMet)%P

    ! Check MetTime is the right type of time and is finite.
    If (iMetMod == 1 .and. iMet == 1) Then
      If (IsTimeInterval(MetTime)) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateMets', 4)
      End If
      If (IsInfFuture(MetTime) .or. IsInfPast(MetTime)) Then
        Call Message('UNEXPECTED FATAL ERROR in UpdateMets', 4)
      End If
    End If

    ! Is the validity of the met module instance improvable?
    If (MetTime >= CommonMet%TValid) Then

      ! Prepare to update met module instance.
      Call PrepareForUpdateMet( &
             Coords, Grids,     &
             iCase, iMetCase,   &
             iMetMod, iMet,     &
             MetTime,           &
             TValid, UpdateNow, &
             Mets,              &
             Units              &
           )

      If (UpdateNow) Then

        ! Try to update met module instance.
        Call UpdateMet(         &
               Coords, Grids,   &
               iCase, iMetCase, &
               iMetMod, iMet,   &
               MetTime,         &
               Mets,            &
               Units            &
             )
        ! Update OverallTValid.
        OverallTValid = TMin(OverallTValid, CommonMet%TValid)

      Else If (.not. CommonMet%UpdateOnDemand) Then

        ! Try to update met module instance.
        Call UpdateMet(         &
               Coords, Grids,   &
               iCase, iMetCase, &
               iMetMod, iMet,   &
               MetTime,         &
               Mets,            &
               Units            &
             )
        ! Update OverallTValid. For same results whether update-on-demand or update-at-once is used, update
        ! OverallTValid using the local variable TValid for update-at-once, as would be used for
        ! update-on-demand.
        If (.not. SameResultsWithUpdateOnDemand) Then
          OverallTValid = TMin(OverallTValid, CommonMet%TValid)
        Else
          OverallTValid = TMin(OverallTValid, TValid)
        End If

      Else

        ! Update OverallTValid.
        OverallTValid = TMin(OverallTValid, TValid)

      End If

    Else

      ! Update OverallTValid.
      OverallTValid = TMin(OverallTValid, CommonMet%TValid)

    End If

  End Do
  End Do

  ! Tell the IO Thread to prefetch any data
  If (ParallelMetRead) Then
    Call LookaheadFileReadRequest()
  End If

End Subroutine UpdateMets

!-------------------------------------------------------------------------------------------------------------

Subroutine PrepareForUpdateMet( &
             Coords, Grids,     &
             iCase, iMetCase,   &
             iMetMod, iMet,     &
             MetTime,           &
             TValid, UpdateNow, &
             Mets,              &
             Units              &
           )
! Prepares for updating the state of a met module instance.

  Implicit None
  ! Argument list:
  Type(Coords_), Intent(In)    :: Coords    ! Collection of coord systems.
  Type(Grids_),  Intent(In)    :: Grids     ! Collection of grids.
  Integer,       Intent(In)    :: iCase     ! Number of case.
  Integer,       Intent(In)    :: iMetCase  ! Number of the met realisation in the met ensemble.
  Integer,       Intent(In)    :: iMetMod   ! Met module index.
  Integer,       Intent(In)    :: iMet      ! Met module instance index.
  Type(Time_),   Intent(In)    :: MetTime   ! Time for which the met module instance is to be updated (this
                                            ! must be the current time unless we have a fixed met case, in
                                            ! which case it must be the time of the fixed met).
  Type(Time_),   Intent(Out)   :: TValid    ! Earliest time that the validity (overall or for any single
                                            ! attribute) of the met module instance might change, assuming the
                                            ! met module instance is updated now. The value is that determined
                                            ! at the end of this routine (the actual time may be later).
  Logical,       Intent(Out)   :: UpdateNow ! Indicates the met module instance must be updated now (even if
                                            ! update-on-demand is specified). If set, TValid need not be set
                                            ! to any particular time.
  Type(Mets_),   Intent(InOut) :: Mets      ! Collection of met module instance states.
  Type(Units_),  Intent(InOut) :: Units     ! Collection of information on input/output unit numbers.
  ! Locals:
  Type(CommonMet_), Pointer :: CommonMet
  ! CommonMet :: Abbreviation for the common part of the met module instance.

  CommonMet => Mets%C(iMetMod, iMet)%P

  If (.false.) Then

  Else If (                               &
    Mets%nPrototypeMets /= 0 .and.             &
    iMetMod == Mets%PrototypeMets(1)%C%iMetMod &
  ) Then

    Call PrepareForUpdatePrototypeMet( &
           Coords, Grids,         &
           iCase, iMetCase,       &
           MetTime,               &
           TValid, UpdateNow,     &
           Mets%PrototypeMets(iMet),   &
           Units                  &
         )

  Else If (                               &
    Mets%nSingleSiteMets /= 0 .and.             &
    iMetMod == Mets%SingleSiteMets(1)%C%iMetMod &
  ) Then

    Call PrepareForUpdateSingleSiteMet( &
           Coords, Grids,         &
           iCase, iMetCase,       &
           MetTime,               &
           TValid, UpdateNow,     &
           Mets%SingleSiteMets(iMet),   &
           Units                  &
         )

  Else If (                               &
    Mets%nNWPMets /= 0 .and.             &
    iMetMod == Mets%NWPMets(1)%C%iMetMod &
  ) Then

    Call PrepareForUpdateNWPMet( &
           Coords, Grids,         &
           iCase, iMetCase,       &
           MetTime,               &
           TValid, UpdateNow,     &
           Mets%NWPMets(iMet),   &
           Units                  &
         )

  Else If (                               &
    Mets%nRadarMets /= 0 .and.             &
    iMetMod == Mets%RadarMets(1)%C%iMetMod &
  ) Then

    Call PrepareForUpdateRadarMet( &
           Coords, Grids,         &
           iCase, iMetCase,       &
           MetTime,               &
           TValid, UpdateNow,     &
           Mets%RadarMets(iMet),   &
           Units                  &
         )

  Else If (                               &
    Mets%nAncillaryMets /= 0 .and.             &
    iMetMod == Mets%AncillaryMets(1)%C%iMetMod &
  ) Then

    Call PrepareForUpdateAncillaryMet( &
           Coords, Grids,         &
           iCase, iMetCase,       &
           MetTime,               &
           TValid, UpdateNow,     &
           Mets%AncillaryMets(iMet),   &
           Units                  &
         )


  Else
    Call Message('UNEXPECTED FATAL ERROR in PrepareForUpdateMet: met module not found', 4)
  End If

  ! Check TValid.
  If (.not. UpdateNow) Then

    ! Check TValid is the right type of time.
    If (IsTimeInterval(TValid)) Then
      Call Message('UNEXPECTED FATAL ERROR in PrepareForUpdateMet', 4)
    End If

    ! Check that TValid is in the future, and, for fixed met, is infinite.
    If (TValid <= MetTime) Then
      Call Message('UNEXPECTED FATAL ERROR in PrepareForUpdateMet', 4)
    End If
    If (CommonMet%FixedMet) Then
      If (.not. IsInfFuture(TValid)) Then
        Call Message('UNEXPECTED FATAL ERROR in PrepareForUpdateMet', 4)
      End If
    End If

  End If

  ! Update DueForUpdate.
  CommonMet%DueForUpdate = .true.

End Subroutine PrepareForUpdateMet

!-------------------------------------------------------------------------------------------------------------

Subroutine UpdateMet(         &
             Coords, Grids,   &
             iCase, iMetCase, &
             iMetMod, iMet,   &
             MetTime,         &
             Mets,            &
             Units            &
           )
! Updates the state of a met module instance.

  Implicit None
  ! Argument list:
  Type(Coords_), Intent(In)    :: Coords   ! Collection of coord systems.
  Type(Grids_),  Intent(In)    :: Grids    ! Collection of grids.
  Integer,       Intent(In)    :: iCase    ! Number of case.
  Integer,       Intent(In)    :: iMetCase ! Number of the met realisation in the met ensemble.
  Integer,       Intent(In)    :: iMetMod  ! Met module index.
  Integer,       Intent(In)    :: iMet     ! Met module instance index.
  Type(Time_),   Intent(In)    :: MetTime  ! Time for which the met module instance is to be updated (this
                                           ! must be the current time unless we have a fixed met case, in
                                           ! which case it must be the time of the fixed met).
  Type(Mets_),   Intent(InOut) :: Mets     ! Collection of met module instance states.
  Type(Units_),  Intent(InOut) :: Units    ! Collection of information on input/output unit numbers.
  ! Locals:
  Type(CommonMet_), Pointer :: CommonMet ! Abbreviation for the common part of the met module instance.

  CommonMet => Mets%C(iMetMod, iMet)%P

  If (.false.) Then

  Else If (                               &
    Mets%nPrototypeMets /= 0 .and.             &
    iMetMod == Mets%PrototypeMets(1)%C%iMetMod &
  ) Then

    Call UpdatePrototypeMet(         &
           Coords, Grids,       &
           iCase, iMetCase,     &
           MetTime,             &
           Mets%PrototypeMets(iMet), &
           Units                &
         )

  Else If (                               &
    Mets%nSingleSiteMets /= 0 .and.             &
    iMetMod == Mets%SingleSiteMets(1)%C%iMetMod &
  ) Then

    Call UpdateSingleSiteMet(         &
           Coords, Grids,       &
           iCase, iMetCase,     &
           MetTime,             &
           Mets%SingleSiteMets(iMet), &
           Units                &
         )

  Else If (                               &
    Mets%nNWPMets /= 0 .and.             &
    iMetMod == Mets%NWPMets(1)%C%iMetMod &
  ) Then

    Call UpdateNWPMet(         &
           Coords, Grids,       &
           iCase, iMetCase,     &
           MetTime,             &
           Mets%NWPMets(iMet), &
           Units                &
         )

  Else If (                               &
    Mets%nRadarMets /= 0 .and.             &
    iMetMod == Mets%RadarMets(1)%C%iMetMod &
  ) Then

    Call UpdateRadarMet(         &
           Coords, Grids,       &
           iCase, iMetCase,     &
           MetTime,             &
           Mets%RadarMets(iMet), &
           Units                &
         )

  Else If (                               &
    Mets%nAncillaryMets /= 0 .and.             &
    iMetMod == Mets%AncillaryMets(1)%C%iMetMod &
  ) Then

    Call UpdateAncillaryMet(         &
           Coords, Grids,       &
           iCase, iMetCase,     &
           MetTime,             &
           Mets%AncillaryMets(iMet), &
           Units                &
         )


  Else
    Call Message('UNEXPECTED FATAL ERROR in UpdateMet: met module not found', 4)
  End If

  ! Check TValid is the right type of time.
  If (IsTimeInterval(CommonMet%TValid)) Then
    Call Message('UNEXPECTED FATAL ERROR in UpdateMet', 4)
  End If

  ! Check that TValid is in the future, and, for fixed met, is infinite.
  If (CommonMet%TValid <= MetTime) Then
    Call Message('UNEXPECTED FATAL ERROR in UpdateMet', 4)
  End If
  If (CommonMet%FixedMet) Then
    If (.not. IsInfFuture(CommonMet%TValid)) Then
      Call Message('UNEXPECTED FATAL ERROR in UpdateMet', 4)
    End If
  End If

  ! Messages for whether met update was sucessful.
  If (CommonMet%Valid) Then
    If (CommonMet%FixedMet) Then
      Call Message(                        &
             'Met module "'             // &
             Trim(CommonMet%MetModName) // &
             '.'                        // &
             Trim(CommonMet%MetName)    // &
             '" now valid'                 &
           )
    Else
      Call Message(                                                   &
             'Met module "'                                        // &
             Trim(CommonMet%MetModName)                            // &
             '.'                                                   // &
             Trim(CommonMet%MetName)                               // &
             '" now valid until '                                  // &
             Trim(Time2Char(CommonMet%TValid, .false., 0, .true.))    &
           )
    End If
  Else
    Call Message(                             &
           'Unable to update met module "' // &
           Trim(CommonMet%MetModName)      // &
           '.'                             // &
           Trim(CommonMet%MetName)         // &
           '"'                                &
         )
  EndIf

  ! Update DueForUpdate.
  CommonMet%DueForUpdate = .false.

End Subroutine UpdateMet

!-------------------------------------------------------------------------------------------------------------

Subroutine ResetMets(Mets)
! Resets the state of the met module instances for a new realisation.

  Implicit None
  ! Argument list:
  Type(Mets_), Intent(InOut) :: Mets ! Collection of met module instance states.
  ! Locals:
  Integer :: iMetMod ! Met module index.
  Integer :: iMet    ! Met module instance index.

  Do iMetMod = 1, Mets%nMetMods
  Do iMet    = 1, Mets%nMets(iMetMod)

    Mets%C(iMetMod, iMet)%P%Valid        = .false.
    Mets%C(iMetMod, iMet)%P%TValid       = InfPastTime()
    Mets%C(iMetMod, iMet)%P%DueForUpdate = .false.

  End Do
  End Do

  Do iMet = 1, Mets%nSingleSiteMets
    Call ResetSingleSiteMet(Mets%SingleSiteMets(iMet))
  End Do
  Do iMet = 1, Mets%nNWPMets
    Call ResetNWPMet(Mets%NWPMets(iMet))
  End Do
  Do iMet = 1, Mets%nAncillaryMets
    Call ResetAncillaryMet(Mets%AncillaryMets(iMet))
  End Do

End Subroutine ResetMets

!-------------------------------------------------------------------------------------------------------------

End Module MetsModule
