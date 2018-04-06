! Module:  Particle Size Distribution Module

Module SizeDistModule

! This module provides code for handling particle size distributions.

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: SizeDist_          ! A particle size distribution.
Public :: SizeDists_         ! A collection of particle size distributions.
Public :: InitSizeDists      ! Initialises a collection of particle size distributions.
Public :: InitSizeDist       ! Initialises a particle size distribution.
Public :: AddSizeDist        ! Adds a particle size distribution to a collection of particle size 
                             ! distributions.
Public :: FindSizeDistIndex  ! Finds the index of a particle size distribution in a collection of particle 
                             ! size distributions.

!-------------------------------------------------------------------------------------------------------------

Type :: SizeDist_ ! A particle size distribution.
  Character(MaxCharLength) :: Name
  Integer                  :: nSizeRanges
  Real(Std)                :: DiameterRangeBoundary(MaxDiameterRanges + 1)
  Real(Std)                :: CumulativeFrac(MaxDiameterRanges + 1)
  Real(Std)                :: Density(MaxDiameterRanges)
  Real(Std)                :: ParticleShape(MaxDiameterRanges)
  Real(Std)                :: RepresentativeDiameter(MaxDiameterRanges)
  Logical                  :: CumulativeFracPresent
  Logical                  :: DensityPresent
  Logical                  :: ParticleShapePresent
  ! Name                   :: Name of particle size distribution.
  ! nSizeRanges            :: Number of particle size ranges.
  ! Diameter               :: Diameters defining the boundaries of the particle size ranges.
  ! CumulativeFrac         :: Cumulative fraction of particles up to each of the diameters.
  ! Density                :: Density of particles in each size range.
  ! ParticleShape          :: Shape of particles in each size range.
  ! RepresentativeDiameter :: Representative diameter for each size range.
  ! CumulativeFracPresent  :} Indicates whether CumulativeFrac, Density and ParticleShape, 
  ! DensityPresent         :} which are optional elements, 
  ! ParticleShapePresent   :} are present.

  ! The distribution is the distribution in terms of the material unit used for the
  ! species. Usually this will be a mass unit but could also be number of particles.
  ! The points at which the cumulative distribution is specified are joined together
  ! by means of straight lines on a log-log plot.

End Type SizeDist_

!-------------------------------------------------------------------------------------------------------------

Type :: SizeDists_ ! A collection of particle size distributions.
  Integer         :: nSizeDists              ! Number of particle size distributions.
  Type(SizeDist_) :: SizeDists(MaxSizeDists) ! Collection of particle size distributions.
End Type SizeDists_

!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of particle size distributions.
  Module Procedure SizeDistEq
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitSizeDists() Result(SizeDists)
! Initialises a collection of particle size distributions.

  Implicit None
  ! Function result:
  Type(SizeDists_) :: SizeDists ! Initialised collection of particle size distributions.

  SizeDists%nSizeDists = 0

End Function InitSizeDists

!-------------------------------------------------------------------------------------------------------------

Function InitSizeDist(                                                          &
           Name,                                                                &
           DiameterRangeBoundary,                                               &
           CumulativeFrac, Density, ParticleShape, RepresentativeDiameter,      &
           CumulativeFracPresent, DensityPresent,                               &
           ParticleShapePresent, RepresentativeDiameterPresent                  &
         )                                                                      &
Result(SizeDist)
! Initialises a particle size distribution.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name
  Real(Std),    Intent(In) :: DiameterRangeBoundary(:)
  Real(Std),    Intent(In) :: CumulativeFrac(:)
  Real(Std),    Intent(In) :: Density(:)
  Real(Std),    Intent(In) :: ParticleShape(:)
  Real(Std),    Intent(In) :: RepresentativeDiameter(:)
  Logical,      Intent(In) :: CumulativeFracPresent
  Logical,      Intent(In) :: DensityPresent
  Logical,      Intent(In) :: ParticleShapePresent
  Logical,      Intent(In) :: RepresentativeDiameterPresent
  ! Name                          :: Name of particle size distribution.
  ! DiameterRangeBoundary         :: Diameters defining the boundaries of the particle size ranges.
  ! CumulativeFrac                :: Cumulative fraction of particles up to each of the diameters.
  ! Density                       :: Density of particles in each size range.
  ! ParticleShape                 :: Shape of particles in each size range.
  ! RepresentativeDiameter        :: Representative diameter for each size range.
  ! CumulativeFracPresent         :} Indicates whether CumulativeFrac, Density, 
  ! DensityPresent                :} Particle Shape and RepresentativeDiameter  
  ! ParticleShapePresent          :} are present (i.e. defined).
  ! RepresentativeDiameterPresent :}
  ! Function result:
  Type(SizeDist_) :: SizeDist ! Initialised particle size distribution.
  ! Locals:
  Integer :: iSizeRange ! Index of particle size range.

  If (Name == ' ') Then
    Call Message(                                            &
           'FATAL ERROR in InitSizeDist: '                // &
           'name of particle size distribution is blank',    &
           3                                                 &
         )
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message(                                                &
           'FATAL ERROR in InitSizeDist: '                    // &
           'name of particle size distribution is given as "' // &
           Trim(Name)                                         // &
           '" and is too long',                                  &
           3                                                     &
         )
  End If

  If (Size(DiameterRangeBoundary) > MaxDiameterRanges + 1) Then
    Call Message(                                                        &
           'FATAL ERROR in InitSizeDist: '                            // &
           'too many diameter ranges in particle size distribution "' // &
           Trim(Name)                                                 // &
           '"',                                                          &
           3                                                             &
         )
  End If
  If (Size(DiameterRangeBoundary) < 2) Then
    Call Message(                                                     &
           'FATAL ERROR in InitSizeDist: '                         // &
           'particle size distribution "'                          // &
           Trim(Name)                                              // &
           ' should have at least two diameter range boundaries"',    &
           3                                                          &
         )
  End If

  ! $$ check DiameterRangeBoundarys and cumulative fractions > 0 and increasing / non-decreasing

  If (CumulativeFracPresent) Then
    If (Size(CumulativeFrac) /= Size(DiameterRangeBoundary)) Then
      Call Message(                              &
             'UNEXPECTED ERROR in InitSizeDist', &
             4                                   &
           )
    End If
    If (CumulativeFrac(1) /= 0.0) Then
      Call Message(                                                                   &
             'FATAL ERROR in InitSizeDist: fraction of particles in lowest range ' // &
             'of particle sizes must be zero',                                        &
             3                                                                        &
           )
    End If
    If (CumulativeFrac(Size(DiameterRangeBoundary)) /= 1.0) Then
      Call Message(                                                             &
             'FATAL ERROR in InitSizeDist: cumulative fraction of particle ' // &
             'sizes must add up to 1',                                          &
             3                                                                  &
           )
    End If
  End If
  
  If (DensityPresent) Then
    If (Size(Density) /= Size(DiameterRangeBoundary) - 1) Then
      Call Message(                              &
             'UNEXPECTED ERROR in InitSizeDist', &
             4                                   &
           )
    End If
    Do iSizeRange = 1, Size(Density)
      If (Density(iSizeRange) <= 0.0) Then
        Call Message(                                                    &
               'FATAL ERROR in InitSizeDist: Particle density is <= 0.', &
               3                                                         &
             )
      End If
    End Do
  End If

  If (ParticleShapePresent) Then
    If (Size(ParticleShape) /= Size(DiameterRangeBoundary) - 1) Then
      Call Message(                              &
             'UNEXPECTED ERROR in InitSizeDist', &
             4                                   &
           )
    End If
    Do iSizeRange = 1, Size(ParticleShape)
      If (ParticleShape(iSizeRange) <= 0.0 .or. ParticleShape(iSizeRange) > 1.0) Then
        Call Message(                                                                     &
               'FATAL ERROR in InitSizeDist: The Particle Shape value is inappropriate.', &
               3                                                                          &
             )
      End If
    End Do
  End If

  If (RepresentativeDiameterPresent) Then
    If (Size(RepresentativeDiameter) /= Size(DiameterRangeBoundary) - 1) Then
      Call Message(                              &
             'UNEXPECTED ERROR in InitSizeDist', &
             4                                   &
           )
    End If
    Do iSizeRange = 1, Size(RepresentativeDiameter)
      If (RepresentativeDiameter(iSizeRange) <= 0.0) Then
        Call Message(                                                           &
               'FATAL ERROR in InitSizeDist: Representative diameter is <= 0.', &
               3                                                                &
             )
      End If
    End Do
  End If

  SizeDist%Name                                              = Name
  SizeDist%nSizeRanges                                       = Size(DiameterRangeBoundary) - 1
  SizeDist%DiameterRangeBoundary(1:SizeDist%nSizeRanges + 1) = DiameterRangeBoundary(1:SizeDist%nSizeRanges + 1)
  SizeDist%CumulativeFracPresent                             = CumulativeFracPresent
  SizeDist%DensityPresent                                    = DensityPresent
  SizeDist%ParticleShapePresent                              = ParticleShapePresent
  If (CumulativeFracPresent) Then
    SizeDist%CumulativeFrac(1:SizeDist%nSizeRanges + 1) = CumulativeFrac(1:SizeDist%nSizeRanges + 1)
  End If
  If (DensityPresent) Then
    SizeDist%Density(1:SizeDist%nSizeRanges) = Density(1:SizeDist%nSizeRanges)
  End If
  If (ParticleShapePresent) Then
    SizeDist%ParticleShape(1:SizeDist%nSizeRanges) = ParticleShape(1:SizeDist%nSizeRanges)
  End If
  If (RepresentativeDiameterPresent) Then
    SizeDist%RepresentativeDiameter(1:SizeDist%nSizeRanges) = RepresentativeDiameter(1:SizeDist%nSizeRanges)
  Else
    Do iSizeRange = 1, SizeDist%nSizeRanges
      SizeDist%RepresentativeDiameter(iSizeRange) =                                                       &
        Sqrt(SizeDist%DiameterRangeBoundary(iSizeRange) * SizeDist%DiameterRangeBoundary(iSizeRange + 1))
    End Do
  End If

End Function InitSizeDist

!-------------------------------------------------------------------------------------------------------------

Subroutine AddSizeDist(SizeDist, SizeDists)
! Adds a particle size distribution to a collection of particle size distributions.

  Implicit None
  ! Argument list:
  Type(SizeDist_),  Intent(In)    :: SizeDist  ! Particle size distribution.
  Type(SizeDists_), Intent(InOut) :: SizeDists ! Collection of particle size distributions.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, SizeDists%nSizeDists
    If (SizeDist%Name .CIEq. SizeDists%SizeDists(i)%Name) Then
      If (SizeDist == SizeDists%SizeDists(i)) Then
        Return
      Else
        Call Message(                                                             &
               'FATAL ERROR in adding the particle size distribution "'        // &
               Trim(SizeDist%Name)                                             // &
               '": a different particle size distribution with the same name ' // &
               'already exists',                                                  &
               3                                                                  &
             )
      End If
    End If
  End Do

  If (SizeDists%nSizeDists >= MaxSizeDists) Then
    Call Message(                                                      &
           'FATAL ERROR in adding the particle size distribution "' // &
           Trim(SizeDist%Name)                                      // &
           '": there are too many particle size distributions',        &
           3                                                           &
         )
  End If

  SizeDists%nSizeDists                      = SizeDists%nSizeDists + 1
  SizeDists%SizeDists(SizeDists%nSizeDists) = SizeDist

End Subroutine AddSizeDist

!-------------------------------------------------------------------------------------------------------------

Function FindSizeDistIndex(Name, SizeDists)
! Finds the index of a particle size distribution in a collection of particle size
! distributions.

  Implicit None
  ! Argument list:
  Character(*),     Intent(In) :: Name      ! Name of particle size distribution.
  Type(SizeDists_), Intent(In) :: SizeDists ! Collection of particle size distributions.
  ! Function result:
  Integer :: FindSizeDistIndex ! Index of the particle size distribution.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, SizeDists%nSizeDists
    If (Name == SizeDists%SizeDists(i)%Name) Then
      FindSizeDistIndex = i
      Return
    End If
  End Do

  Call Message(                                         &
         'FATAL ERROR: particle size distribution "' // &
         Trim(Name)                                  // &
         '" not found',                                 &
         3                                              &
       )

End Function FindSizeDistIndex

!-------------------------------------------------------------------------------------------------------------

Function SizeDistEq(SizeDist1, SizeDist2)
! Tests for equality of particle size distributions.

  Implicit None
  ! Argument list:
  Type(SizeDist_), Intent(In) :: SizeDist1 !} The two particle size distributions.
  Type(SizeDist_), Intent(In) :: SizeDist2 !}
  ! Function result:
  Logical :: SizeDistEq ! Indicates if particle size distributions are equal.
  ! Locals:
  Integer :: i ! Loop index.

  SizeDistEq = (SizeDist1%Name                  .CIEq. SizeDist2%Name                 ) .and. &
                SizeDist1%nSizeRanges             ==   SizeDist2%nSizeRanges            .and. &
               (SizeDist1%CumulativeFracPresent .eqv.  SizeDist2%CumulativeFracPresent) .and. &
               (SizeDist1%DensityPresent        .eqv.  SizeDist2%DensityPresent       ) .and. &
               (SizeDist1%ParticleShapePresent  .eqv.  SizeDist2%ParticleShapePresent )
               
  If (.not. SizeDistEq) Return
  
  Do i = 1, SizeDist1%nSizeRanges + 1
    SizeDistEq = SizeDistEq .and. SizeDist1%DiameterRangeBoundary(i) == SizeDist2%DiameterRangeBoundary(i)
  End Do
  If (SizeDist1%CumulativeFracPresent) Then
    Do i = 1, SizeDist1%nSizeRanges + 1
      SizeDistEq = SizeDistEq .and. SizeDist1%CumulativeFrac(i) == SizeDist2%CumulativeFrac(i)
    End Do
  End If
  If (SizeDist1%DensityPresent) Then
    Do i = 1, SizeDist1%nSizeRanges
      SizeDistEq = SizeDistEq .and. SizeDist1%Density(i) == SizeDist2%Density(i)
    End Do
  End If
  If (SizeDist1%ParticleShapePresent) Then
    Do i = 1, SizeDist1%nSizeRanges
      SizeDistEq = SizeDistEq .and. SizeDist1%ParticleShape(i) == SizeDist2%ParticleShape(i)
    End Do
  End If
  Do i = 1, SizeDist1%nSizeRanges
    SizeDistEq = SizeDistEq .and. SizeDist1%RepresentativeDiameter(i) == SizeDist2%RepresentativeDiameter(i)
  End Do

End Function SizeDistEq

!-------------------------------------------------------------------------------------------------------------

End Module SizeDistModule
