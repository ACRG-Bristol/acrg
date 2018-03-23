! Module:  Grid And Domain Module

Module GridAndDomainModule

! This module provides code for handling grids and domains.

! Module overview
! ---------------

! $$
! Each grid also has an averaging region of size > 0 associated with each grid point. This can be specified
! explicitly or, if not, is e.g. taken to be dZ or determined from the z(i). This is assumed in Output.f90
!
! $$ Need to add grids defined in terms of BL to support bl av of met.

! Module use
! ----------

! $$

! Module call tree
! ----------------

! $$

!-------------------------------------------------------------------------------------------------------------

Use GlobalParametersModule
Use ErrorAndMessageModule
Use MathsModule
Use PhysicsModule
Use StringModule
Use TimeModule
Use CoordinateSystemModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: Locations_
Public  :: Locationses_
Public  :: HGrid_
Public  :: ZGrid_
Public  :: TGrid_
Public  :: DGrid_
Public  :: Grids_
Public  :: GHCoeffs_
Public  :: HCoeffs_
Public  :: ZCoeffs_
Public  :: TCoeffs_
Public  :: Domain_
Public  :: Domains_
Public  :: Operator(==)
Public  :: InitLocationses
Public  :: InitLocations
Public  :: AddLocations
Public  :: FindLocationsIndex
Public  :: FindLocationIndex
Public  :: InitGrids
Public  :: InitHGrid
Public  :: InitZGrid
Public  :: InitTGrid
Public  :: InitDGrid
Public  :: AddHGrid
Public  :: AddZGrid
Public  :: AddTGrid
Public  :: AddDGrid
Public  :: FindHGridIndex
Public  :: FindZGridIndex
Public  :: FindTGridIndex
Public  :: FindDGridIndex
Public  :: SetUpGrids_CoordsEtc
Public  :: SetUpGrids_iCoords
Public  :: TInTGrid
Public  :: FirstTAfterT
Public  :: LastTBeforeT
Public  :: nDInDGrid
Public  :: FloatingDGrid
Public  :: DInDGrid
Public  :: ShiftIndexInDGrid
Public  :: InitGHCoeffs
Public  :: GetGHCoeffs
Public  :: GInterpXY
Public  :: GInterpXYZ
Public  :: GetHCoeffs
Public  :: GetZCoeffs
Public  :: GetTCoeffs
Public  :: InterpXYZT
Public  :: InterpXYT
Public  :: InterpCXYT
Public  :: InterpXY
Public  :: CalcArea
Public  :: InitDomains
Public  :: InitDomain
Public  :: AddDomain
Public  :: FindDomainIndex
Public  :: SetUpDomains_CoordsEtc
Public  :: SetUpDomains_iCoords
Public  :: StartTimeOfDomain
Public  :: EndTimeOfDomain
Public  :: DomainStart
Public  :: DomainEnd
Public  :: XInDomain
Public  :: HInDomain
Public  :: ZInDomain
Public  :: SInDomain
Public  :: TInDomain
Public  :: CalcZenithAngle

!-------------------------------------------------------------------------------------------------------------

Type :: Locations_ ! A set of locations.
  Character(MaxCharLength)          :: Name
  Integer                           :: nLocations ! $$ Needed?
  Character(MaxCharLength), Pointer :: Names(:)
  Character(MaxCharLength), Pointer :: HCoordNames(:)
  Real(Std),                Pointer :: X(:)
  Real(Std),                Pointer :: Y(:)
End Type Locations_

!-------------------------------------------------------------------------------------------------------------

Type :: Locationses_ ! A collection of sets of locations.
  Integer          :: nLocationses
  Type(Locations_) :: Locationses(MaxLocationses)
End Type Locationses_

!-------------------------------------------------------------------------------------------------------------

Type :: HGrid_ ! A horizontal grid.
  Character(MaxCharLength)          :: Name
  Character(MaxCharLength)          :: HCoordName
  Integer                           :: iHCoord
  Logical                           :: Unstructured
  Logical                           :: Variable
  Logical                           :: Named
  Logical                           :: Wrap
  Integer                           :: nX
  Integer                           :: nY
  Real(Std)                         :: dX
  Real(Std)                         :: dY
  Real(Std)                         :: X0
  Real(Std)                         :: Y0
  Real(Std)                         :: XCentre
  Real(Std)                         :: YCentre
  Real(Std),                Pointer :: X(:)
  Real(Std),                Pointer :: Y(:)
  Real(Std)                         :: XCycle
  Real(Std)                         :: YCycle
  Character(MaxCharLength), Pointer :: Names(:)
  Character(MaxCharLength)          :: LocationsName
  Character(MaxCharLength)          :: CentreName
  ! Name          :: Name of grid.
  ! HCoordName    :: Name of coord system used to define the grid.
  ! iHCoord       :: Index (in Coords) of coord system used to define the grid. 0
  !                  indicates value unknown.
  ! Unstructured  :: Indicates whether the grid is structured, with grid coords of the
  !                  form (X(i), Y(j)), i = 1, nX, j = 1, nY, or unstructured, with
  !                  grid coords of the form (X(i), Y(i)), i = 1, nX.
  ! Variable      :: For structured grids, indicates that the grid spacing is variable.
  ! Named         :: For unstructured grids, indicates that the points are named. For
  !                  structured grids is set to false.
  ! Wrap          :: For structured grids, indicates that the grid wraps round in the
  !                  X-direction. Only relevant if interpolation is to be done on the
  !                  grid.
  ! nX            :} For structured grids, nX and nY are the number of grid points in
  ! nY            :} the coord directions. For unstructured grids, nX = total number of
  !                  grid points and nY = 1.
  ! dX            :] For structured non-variable grids, the grid coords are
  ! dY            :] (X0 + (i - 1) * dX, Y0 + (j - 1) * dY), i = 1, nX, j = 1, nY. Also
  ! X0            :] for structured non-variable grids and unstructured grids, dX and
  ! Y0            :] dY indicate an area to be associated with the grid points.
  ! XCentre       :} grid centre for structured non-variable grids
  ! YCentre       :}
  ! X             :] Coords of the grid.
  ! Y             :]
  ! XCycle        :} Coordinate change in cycling the underlying coordinate for cyclic
  ! YCycle        :} coords. 0 otherwise.
  ! Names         :: For unstructured named grids, gives the names of locations.
  ! LocationsName :} Named centre of structured non-variable grid. Blank indicates
  ! CentreName    :} centre not named or defined.
End Type HGrid_

!-------------------------------------------------------------------------------------------------------------

Type :: ZGrid_ ! A vertical grid.
  Character(MaxCharLength)         :: Name
  Character(MaxCharLength)         :: ZCoordName
  Integer                          :: iZCoord
  Logical                          :: Variable
  Integer                          :: nZ
  Real(Std)                        :: dZ
  Real(Std)                        :: Z0
  Real(Std),               Pointer :: Z(:)
  Real(Std),               Pointer :: AvZ(:)
  Logical                          :: UseiZ
  Integer,                 Pointer :: iZ(:)
  Integer                          :: IncreasingValue
  Integer                          :: IncreasingHeight
  ! Name             :: Name of grid.
  ! ZCoordName       :: Name of coord system used to define the grid.
  ! iZCoord          :: Index (in Coords) of coord system used to define the grid. 0
  !                     indicates value unknown.
  ! Variable         :: Indicates that the grid spacing is variable.
  ! nZ               :: Number of grid points.
  ! dZ               :} If Variable is false, then the grid coords are Z0 + (i - 1) * dZ,
  ! Z0               :} i = 1, nZ. If Variable is true, then dZ and Z0 are undefined.
  ! Z                :: Coords of the grid.
  ! AvZ              :: Heights over which vertical averages are to be calculated.
  ! UseiZ            :: Indicates array iZ is used.
  ! iZ               :: Indices used to label the points of the grid as an alternative to 1, ..., nZ.
  ! IncreasingValue  :: 1 indicates that value increases with index (i.e. Z(n + 1) >
  !                     Z(n)), -1 indicates that value decreases with index, and 0
  !                     indicates that value is non-monotonic (including the case nZ =
  !                     1).
  ! IncreasingHeight :: 1 indicates that height increases with index (i.e. Z(n + 1) is
  !                     higher in the atmosphere than Z(n)), -1 indicates that height
  !                     decreases with index, and 0 indicates that height is
  !                     non-monotonic (including the case nZ = 1) or that its
  !                     monotonicity hasn't yet been determined.
End Type ZGrid_

!-------------------------------------------------------------------------------------------------------------

Type :: TGrid_ ! A time grid.
  Character(MaxCharLength)         :: Name
  Logical                          :: Variable
  Integer                          :: nT
  Type(Time_)                      :: dT
  Type(Time_)                      :: T0
  Type(ShortTime_)                 :: sdT
  Type(ShortTime_)                 :: sT0
  Type(ShortTime_),        Pointer :: T_Off(:)
  ! Name     :: Name of grid.
  ! Variable :: Indicates that the grid spacing is variable.
  ! nT       :: Number of grid points.
  ! dT       :} If Variable is false, then the grid coords are T0 + (i - 1) * dT,
  ! T0       :} i = 1, nT. If Variable is true, then dT and T0 are undefined.
  ! T        :: Coords of the grid.
End Type TGrid_

!-------------------------------------------------------------------------------------------------------------

Type :: DGrid_ ! A data grid.
  Private
  Character(MaxCharLength)         :: Name
  Logical                          :: Variable
  Logical                          :: Floating
  Logical                          :: Geometric
  Logical                          :: Zero
  Integer                          :: nD
  Real(Std)                        :: dD
  Real(Std)                        :: D0
  Real(Std),               Pointer :: D(:)
  ! Name      :: Name of grid.
  ! Variable  :} Defines points of grid as indicated below.
  ! Geometric :}
  ! Floating  :}
  ! nD        :}
  ! dD        :}
  ! D0        :}
  ! D         :}

  ! Variable  Floating  Geometric  Zero   Grid values
  ! --------  --------  ---------  ----   -----------
  ! False     False     False      False     D0 + dD *  (i - 1),     i = 1, nD
  ! False     False     True       False     D0 * dD ** (i - 1),     i = 1, nD
  ! False     False     True       True   0, D0 * dD ** (i - 2),     i = 2, nD
  ! False     True      False      False     D0 + dD *  (i + s - 1), i = 1, nD where s can vary
  ! False     True      True       False     D0 * dD ** (i + s - 1), i = 1, nD where s can vary
  ! False     True      True       True   0, D0 * dD ** (i + s - 2), i = 2, nD where s can vary
  ! True      False     False      False  D(i),                      i = 1, nD
End Type DGrid_

!-------------------------------------------------------------------------------------------------------------

Type :: Grids_ ! A collection of grids.
  Integer      :: nHGrids           ! Number of horizontal grids.
  Type(HGrid_) :: HGrids(MaxHGrids) ! Collection of horizontal grids.
  Integer      :: nZGrids           ! Number of vertical grids.
  Type(ZGrid_) :: ZGrids(MaxZGrids) ! Collection of vertical grids.
  Integer      :: nTGrids           ! Number of time grids.
  Type(TGrid_) :: TGrids(MaxTGrids) ! Collection of time grids.
  Integer      :: nDGrids           ! Number of data grids.
  Type(DGrid_) :: DGrids(MaxDGrids) ! Collection of data grids.
End Type Grids_

!-------------------------------------------------------------------------------------------------------------

Type :: GHCoeffs_ ! A set of interpolation coefficients for interpolating from one
                  ! horizontal grid to another.
  Integer            :: n        !} Interpolation coefficients, defined so that a
  Integer,   Pointer :: iX(:, :) !} field F can be interpolated from FromHGrid to
  Integer,   Pointer :: iY(:, :) !} ToHGrid using F_ToHGrid(i, j) = sum_{m = 1}^n
  Real(Std), Pointer :: X(:, :)  !} F_FromHGrid(iX(i, m), iY(j, m)) * X(i, m) *
  Real(Std), Pointer :: Y(:, :)  !} Y(j, m).
  Logical            :: Valid    ! Indicates the information is valid.
End Type GHCoeffs_

!-------------------------------------------------------------------------------------------------------------

Type :: HCoeffs_ ! A set of interpolation coefficients for a horizontal grid.
  Integer   :: iX1      ! Index of neighbouring x-grid point with smaller index.
  Integer   :: iX2      ! Index of neighbouring x-grid point with larger index.
  Integer   :: iY1      ! Index of neighbouring y-grid point with smaller index.
  Integer   :: iY2      ! Index of neighbouring y-grid point with larger index.
  Real(Std) :: X1       ! Interpolation coefficient multiplying iX2 value.
  Real(Std) :: X2       ! Interpolation coefficient multiplying iX1 value.
  Real(Std) :: Y1       ! Interpolation coefficient multiplying iY2 value.
  Real(Std) :: Y2       ! Interpolation coefficient multiplying iY1 value.
  Integer   :: iXNear   ! Nearest X index.
  Integer   :: iYNear   ! Nearest Y index.
  Integer   :: XOutside ! -1, 0 or 1 depending on whether the point is before the first x-grid point, within
                        ! the x-grid, or after the last x-grid point.
  Integer   :: YOutside ! -1, 0 or 1 depending on whether the point is before the first y-grid point, within
                        ! the y-grid, or after the last y-grid point.
  Logical   :: Valid
  ! Indices range from 1 to the nX and nY values of the corresponding horizontal grid
  ! and the results always have iX1 /= iX2 and iY1 /= iY2 (in case the coefficients
  ! are to be modified to calculate derivatives). If the point is outside the grid,
  ! values corresponding to the nearest boundary point are given.
End Type HCoeffs_

!-------------------------------------------------------------------------------------------------------------

Type :: ZCoeffs_ ! A set of interpolation coefficients for a vertical grid.
  Integer   :: iZ1      ! Index of neighbouring z-grid point with larger index. $$ swap ?
  Integer   :: iZ2      ! Index of neighbouring z-grid point with smaller index.$$
  Real(Std) :: Z1       ! Interpolation coefficient multiplying iZ2 value.
  Real(Std) :: Z2       ! Interpolation coefficient multiplying iZ1 value.
  Integer   :: iZNear   ! Nearest Z index.
  Integer   :: ZOutside ! -1, 0 or 1 depending on whether the point is before the first z-grid point, within
                        ! the z-grid, or after the last z-grid point.
  Logical   :: Valid
  ! Indices range from 1 to the nX and nY values of the corresponding vertical grid
  ! and the results always have iZ1 /= iZ2 (in case the coefficients are to be
  ! modified to calculate derivatives). If the point is outside the grid, values
  ! corresponding to the nearest boundary point are given.

  ! $$ The ordering must be iZ2 < iZ1 (see NWPFlow.EtaToZ). $$ check always true even
  ! if grid decreses with height (and swap?).
End Type ZCoeffs_

!-------------------------------------------------------------------------------------------------------------

Type :: TCoeffs_ ! A set of time interpolation coefficients.
  Integer   :: iT1    ! Index of earlier time.
  Integer   :: iT2    ! Index of later time.
  Real(Std) :: T1     ! Interpolation coefficient multiplying iT2 value.
  Real(Std) :: T2     ! Interpolation coefficient multiplying iT1 value.
  Integer   :: iTNear ! Nearest T index.
  Logical   :: Valid
End Type TCoeffs_

!-------------------------------------------------------------------------------------------------------------

Type :: Domain_ ! A space-time domain.
  Character(MaxCharLength) :: Name
  Logical                  :: XUnbounded
  Logical                  :: YUnbounded
  Character(MaxCharLength) :: HCoordName
  Integer                  :: iHCoord
  Real(Std)                :: XMin
  Real(Std)                :: XMax
  Real(Std)                :: YMin
  Real(Std)                :: YMax
  Logical                  :: ZUnbounded
  Character(MaxCharLength) :: ZCoordName
  Integer                  :: iZCoord
  Real(Std)                :: ZMax
  Type(Time_)              :: StartTime
  Type(Time_)              :: EndTime   ! StopTime ?
  Type(Time_)              :: MaxTravelTime
  Type(ShortTime_)         :: sStartTime
  Type(ShortTime_)         :: sEndTime
  Type(ShortTime_)         :: sMaxTravelTime
  Character(MaxCharLength) :: LocationsName
  Character(MaxCharLength) :: CentreName
  ! XUnbounded    :: Indicates the domain is unbounded horizontally in the X direction.
  ! YUnbounded    :: Indicates the domain is unbounded horizontally in the Y direction.
  ! HCoordName    :: Name of horizontal coord system used to define the domain
  !                  (defined only if XUnbounded or YUnbounded is false).
  ! iHCoord       :: Index in Coords of horizontal coord system used to define the
  !                  domain. 0 indicates domain unbounded horizontally.
  ! XMin          :} Max and min values of the horizontal coords (defined only if
  ! XMax          :} XUnbounded or YUnbounded respectively are false).
  ! YMin          :}
  ! YMax          :}
  ! ZUnbounded    :: Indicates the domain is unbounded vertically.
  ! ZCoordName    :: Name of vertical coord system used to define the domain (defined
  !                  only if ZUnbounded is false).
  ! iZCoord       :: Index in Coords of vertical coord system used to define the
  !                  domain. 0 indicates domain unbounded vertically.
  ! ZMax          :: Top of domain (defined only if ZUnbounded is false).
  ! StartTime     :: Start time of domain.
  ! EndTime       :: End time of domain.
  ! MaxTravelTime :: Maximum travel time for which the domain is to be used.
  ! LocationsName :} Named centre of structured non-variable grid. Blank indicates
  ! CentreName    :} centre not named or defined.

  ! add short times.

  ! Domains cannot be null in time and so cannot have EndTime < StartTime or StartTime = infinite future or
  ! EndTime = infinite past

End Type Domain_

!-------------------------------------------------------------------------------------------------------------

Type :: Domains_ ! A collection of domains.
  Integer       :: nDomains
  Type(Domain_) :: Domains(MaxDomains)
End Type Domains_

!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of coord systems, grids and eta definitions.
  Module Procedure LocationsEq
  Module Procedure HGridEq
  Module Procedure ZGridEq
  Module Procedure TGridEq
  Module Procedure DGridEq
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function InitLocationses() Result (Locationses)
! Initialises a collection of sets of locations.

  Implicit None
  ! Function result:
  Type(Locationses_) :: Locationses ! Initialised collection of sets of locations.

  Locationses%nLocationses = 0

End Function InitLocationses

!-------------------------------------------------------------------------------------------------------------

Function InitLocations(Name, Names, HCoordNames, X, Y) Result (Locations)
! Initialises a set of locations.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name           !
  Character(*), Intent(In) :: Names(:)       !
  Character(*), Intent(In) :: HCoordNames(:) !
  Real(Std),    Intent(In) :: X(:)           !
  Real(Std),    Intent(In) :: Y(:)           !
  ! Function result:
  Type(Locations_) :: Locations ! Initialised set of locations.
  ! Locals:
  Integer :: i, j

  If (Name == ' ') Then
    Call Message(                                                          &
           'FATAL ERROR in InitLocations: set of locations name is blank', &
           3                                                               &
         )
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message(                                   &
           'FATAL ERROR in InitLocations: '      // &
           'set of locations name is given as "' // &
           Trim(Name)                            // &
           '" and is too long',                     &
           3                                        &
         )
  End If

  ! Check dimensions consistent, non-zero $$ - UNEXPECTED FATAL ERROR
  ! Check Names all blank or all non-blank $$ - FATAL ERROR (also name lengths)
  ! Check HCoordNames all non-blank $$ - FATAL ERROR (also name lengths)
  ! Check Names all different $$

  ! This is a rough go at checking for identical names.
  Do j = 2, Size(HCoordNames)
    Do i = 1, j - 1
      If (Names(j) .CIEq. Names(i)) Then
!        If (Location(j) == Location(i)) Then ! Restore this bit when we have an eq operator for locations $$
!          Exit
!        Else
          Call Message(                                                  &
                 'FATAL ERROR in initialising the set of locations "' // &
                 Trim(Name)                                           // &
                 '": two or more locations have the same name',          &
                 3                                                       &
               )
!        End If
      End If
    End Do
  End Do

  Allocate(Locations%Names      (Size(HCoordNames)))
  Allocate(Locations%HCoordNames(Size(HCoordNames)))
  Allocate(Locations%X          (Size(HCoordNames)))
  Allocate(Locations%Y          (Size(HCoordNames)))

  Locations%Name                                = Name
  Locations%nLocations                          = Size(Names)
  Locations%Names      (1:Locations%nLocations) = Names      (1:Locations%nLocations)
  Locations%HCoordNames(1:Locations%nLocations) = HCoordNames(1:Locations%nLocations)
  Locations%X          (1:Locations%nLocations) = X          (1:Locations%nLocations)
  Locations%Y          (1:Locations%nLocations) = Y          (1:Locations%nLocations)

End Function InitLocations

!-------------------------------------------------------------------------------------------------------------

Subroutine AddLocations(Locations, Locationses)
! Adds a set of locations to a collection of sets of locations.

  Implicit None
  ! Argument list:
  Type(Locations_),   Intent(In)    :: Locations   ! Set of locations to be added.
  Type(Locationses_), Intent(InOut) :: Locationses ! Collection of sets of locations.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Locationses%nLocationses
    If (Locations%Name .CIEq. Locationses%Locationses(i)%Name) Then
      If (Locations == Locationses%Locationses(i)) Then
        Return
      Else
        Call Message(                                              &
               'FATAL ERROR in adding the set of locations "'   // &
               Trim(Locations%Name)                             // &
               '": a different set of locations with the same ' // &
               'name already exists',                              &
               3                                                   &
             )
      End If
    End If
  End Do

  If (Locationses%nLocationses >= MaxLocationses) Then
    Call Message(                                            &
           'FATAL ERROR in adding the set of locations "' // &
           Trim(Locations%Name)                           // &
           '": there are too many sets of locations',        &
           3                                                 &
         )
  End If

  Locationses%nLocationses                          = Locationses%nLocationses + 1
  Locationses%Locationses(Locationses%nLocationses) = Locations

End Subroutine AddLocations

!-------------------------------------------------------------------------------------------------------------

Function FindLocationsIndex(Name, Locationses)
! Finds the index of a set of locations.

  Implicit None
  ! Argument list:
  Character(*),       Intent(In) :: Name        ! Name of set of locations.
  Type(Locationses_), Intent(In) :: Locationses ! Collection of sets of locations.
  ! Function result:
  Integer :: FindLocationsIndex ! Index of the set of locations.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Locationses%nLocationses
    If (Name .CIEq. Locationses%Locationses(i)%Name) Then
      FindLocationsIndex = i
      Return
    End If
  End Do

  Call Message(                               &
         'FATAL ERROR: set of locations "' // &
         Trim(Name)                        // &
         '" not found',                       &
         3                                    &
       )

End Function FindLocationsIndex

!-------------------------------------------------------------------------------------------------------------

Function LocationsEq(Locations1, Locations2)
! Tests for equality of sets of locations.

  Implicit None
  ! Argument list:
  Type(Locations_), Intent(In) :: Locations1 !} The two sets of locations.
  Type(Locations_), Intent(In) :: Locations2 !}
  ! Function result:
  Logical :: LocationsEq ! Indicates if sets of locations are equal.
  ! Locals:
  Integer :: i ! Loop index.

  LocationsEq = (Locations1%Name       .CIEq. Locations2%Name      ) .and. &
                 Locations1%nLocations   ==   Locations2%nLocations
  Do i = 1, Min(Locations1%nLocations, Locations2%nLocations)
    LocationsEq =  LocationsEq                                     .and. &
                  (Locations1%Names(i) .CIEq. Locations2%Names(i)) .and. &
                   Locations1%X(i)       ==   Locations2%X(i)      .and. &
                   Locations1%Y(i)       ==   Locations2%Y(i)
  End Do

End Function LocationsEq

!-------------------------------------------------------------------------------------------------------------

Function FindLocationIndex(Name, Locations)
! Finds the index of a location.

  Implicit None
  ! Argument list:
  Character(*),     Intent(In) :: Name      ! Name of location.
  Type(Locations_), Intent(In) :: Locations ! Set of locations.
  ! Function result:
  Integer :: FindLocationIndex ! Index of location.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Locations%nLocations
    If (Name .CIEq. Locations%Names(i)) Then
      FindLocationIndex = i
      Return
    End If
  End Do

  Call Message(                       &
         'FATAL ERROR: location "' // &
         Trim(Name)                // &
         '" not found',               & ! $$ 'in the set of locations " "'
         3                            &
       )

End Function FindLocationIndex

!-------------------------------------------------------------------------------------------------------------

Function InitGrids() Result(Grids)
! Initialises a collection of grids.

  Implicit None
  ! Function result:
  Type(Grids_) :: Grids !

  Grids%nHGrids = 0
  Grids%nZGrids = 0
  Grids%nTGrids = 0
  Grids%nDGrids = 0

End Function InitGrids

!-------------------------------------------------------------------------------------------------------------

Function InitHGrid(                                                         &
           Name, HCoordName, Wrap, nX, nY, dX, dY, X0, Y0, X, Y, Locations, &
           LocationsName, CentreName                                        &
         )                                                                  &
Result (HGrid)
! Initialises a horizontal grid.

  Implicit None
  ! Argument list:
  Character(*),     Intent(In)           :: Name
  Character(*),     Intent(In), Optional :: HCoordName
  Logical,          Intent(In), Optional :: Wrap
  Integer,          Intent(In), Optional :: nX
  Integer,          Intent(In), Optional :: nY
  Real(Std),        Intent(In), Optional :: dX
  Real(Std),        Intent(In), Optional :: dY
  Real(Std),        Intent(In), Optional :: X0
  Real(Std),        Intent(In), Optional :: Y0
  Real(Std),        Intent(In), Optional :: X(:)
  Real(Std),        Intent(In), Optional :: Y(:)
  Type(Locations_), Intent(In), Optional :: Locations
  Character(*),     Intent(In), Optional :: LocationsName
  Character(*),     Intent(In), Optional :: CentreName
  ! Function result:
  Type(HGrid_) :: HGrid !
  ! Locals:
  Integer :: i              ! Loop index.
  Logical :: CentreLocation ! Indicates the grid centre is specified by giving a location.

  If (Len_Trim(Name) > MaxCharLength .or. Name == ' ') Then
    Call Message('FATAL ERROR in InitHGrid', 3)
  End If

  ! Unstructured grids.
  If (                                &
    .not.Present(HCoordName)    .and. &
    .not.Present(Wrap)          .and. &
    .not.Present(nX)            .and. &
    .not.Present(nY)            .and. &
         Present(dX)            .and. &
         Present(dY)            .and. &
    .not.Present(X0)            .and. &
    .not.Present(Y0)            .and. &
    .not.Present(X)             .and. &
    .not.Present(Y)             .and. &
         Present(Locations)     .and. &
    .not.Present(LocationsName) .and. &
    .not.Present(CentreName)          &
  ) Then

    HGrid%Unstructured = .true.

  ! Structured variable grids.
  Else If (                           &
         Present(HCoordName)    .and. &
         Present(Wrap)          .and. &
    .not.Present(nX)            .and. &
    .not.Present(nY)            .and. &
    .not.Present(dX)            .and. &
    .not.Present(dY)            .and. &
    .not.Present(X0)            .and. &
    .not.Present(Y0)            .and. &
         Present(X)             .and. &
         Present(Y)             .and. &
    .not.Present(Locations)     .and. &
    .not.Present(LocationsName) .and. &
    .not.Present(CentreName)          &
  ) Then

    HGrid%Unstructured = .false.
    HGrid%Variable     = .true.

  ! Structured non-variable grids, specified by X0, Y0.
  Else If (                           &
         Present(HCoordName)    .and. &
         Present(Wrap)          .and. &
         Present(nX)            .and. &
         Present(nY)            .and. &
         Present(dX)            .and. &
         Present(dY)            .and. &
         Present(X0)            .and. &
         Present(Y0)            .and. &
    .not.Present(X)             .and. &
    .not.Present(Y)             .and. &
    .not.Present(Locations)     .and. &
    .not.Present(LocationsName) .and. &
    .not.Present(CentreName)          &
  ) Then

    HGrid%Unstructured = .false.
    HGrid%Variable     = .false.
    CentreLocation     = .false.

  ! Structured non-variable grids, specified by LocationsName, CentreName.
  Else If (                           &
         Present(HCoordName)    .and. &
         Present(Wrap)          .and. &
         Present(nX)            .and. &
         Present(nY)            .and. &
         Present(dX)            .and. &
         Present(dY)            .and. &
    .not.Present(X0)            .and. &
    .not.Present(Y0)            .and. &
    .not.Present(X)             .and. &
    .not.Present(Y)             .and. &
    .not.Present(Locations)     .and. &
         Present(LocationsName) .and. &
         Present(CentreName)          &
  ) Then

    HGrid%Unstructured = .false.
    HGrid%Variable     = .false.
    CentreLocation     = .true.

  Else

    Call Message('FATAL ERROR in InitHGrid', 4)

  End If

  ! Unstructured grids.
  If (HGrid%Unstructured) Then

    HGrid%Name = Name
    Do i = 2, Locations%nLocations
      If (.not.(Locations%HCoordNames(i) .CIEq. Locations%HCoordNames(1))) Then
        Call Message(                                                           &
               'FATAL ERROR in initialising the horizontal grid "'           // &
               Trim(Name)                                                    // &
               '": the locations do not all use the same coordinate system',    &
               3                                                                &
             )
      End If
    End Do
    HGrid%HCoordName = Locations%HCoordNames(1)
    HGrid%iHCoord    = 0
    HGrid%Named      = Locations%Names(1) /= ' '
    HGrid%nX         = Locations%nLocations
    HGrid%nY         = 1
    If (HGrid%nX <= 0 .or. HGrid%nY <= 0) Then
      Call Message('FATAL ERROR in InitHGrid', 4)
    End If
    HGrid%dX                        = dX
    HGrid%dY                        = dY
    Allocate(HGrid%X(Locations%nLocations))
    Allocate(HGrid%Y(Locations%nLocations))
    HGrid%X(1:Locations%nLocations) = Locations%X(1:Locations%nLocations)
    HGrid%Y(1:Locations%nLocations) = Locations%Y(1:Locations%nLocations)
    If (HGrid%Named) Then
      Allocate(HGrid%Names(Locations%nLocations))
      HGrid%Names(1:Locations%nLocations) = Locations%Names(1:Locations%nLocations)
    End If
    HGrid%LocationsName = ' '
    HGrid%CentreName    = ' '

    HGrid%XCentre = 0.0
    HGrid%YCentre = 0.0
    ! $$ should adjust locations mod xcycle etc to be close to 0.0 (and store old
    ! locations?) This might make output calculation easier (or might not). Otherwise
    ! no good reason to initialise XCentre, YCentre at all.

  ! Structured variable grids.
  Else If (HGrid%Variable) Then

    HGrid%Name = Name
    If (Len_Trim(HCoordName) > MaxCharLength .or. HCoordName == ' ') Then
      Call Message('FATAL ERROR in InitHGrid', 3)
    End If
    HGrid%HCoordName = HCoordName
    HGrid%iHCoord    = 0
    HGrid%Named      = .false.
    HGrid%Wrap       = Wrap
    HGrid%nX         = Size(X)
    HGrid%nY         = Size(Y)
    If (HGrid%nX <= 0 .or. HGrid%nY <= 0) Then
      Call Message('FATAL ERROR in InitHGrid', 3)
    End If
    Allocate(HGrid%X(HGrid%nX))
    Allocate(HGrid%Y(HGrid%nY))
    HGrid%X(1:HGrid%nX) = X(1:HGrid%nX)
    HGrid%Y(1:HGrid%nY) = Y(1:HGrid%nY)
    HGrid%LocationsName = ' '
    HGrid%CentreName    = ' '

  ! Structured non-variable grids, specified by X0, Y0.
  Else If (.not.CentreLocation) Then

    HGrid%Name = Name
    If (Len_Trim(HCoordName) > MaxCharLength .or. HCoordName == ' ') Then
      Call Message('FATAL ERROR in InitHGrid', 3)
    End If
    HGrid%HCoordName = HCoordName
    HGrid%iHCoord    = 0
    HGrid%Named      = .false.
    HGrid%Wrap       = Wrap
    HGrid%nX         = nX
    HGrid%nY         = nY
    If (HGrid%nX <= 0 .or. HGrid%nY <= 0) Then
      Call Message('FATAL ERROR in InitHGrid', 4)
    End If
    HGrid%dX = dX
    HGrid%dY = dY
    HGrid%X0 = X0
    HGrid%Y0 = Y0
    HGrid%XCentre = X0 + Real(nX - 1, Std) * 0.5 * dX
    HGrid%YCentre = Y0 + Real(nY - 1, Std) * 0.5 * dY
    Allocate(HGrid%X(HGrid%nX))
    Allocate(HGrid%Y(HGrid%nY))
    Do i = 1, nX
      HGrid%X(i) = HGrid%X0 + (i - 1)*HGrid%dX
    End Do
    Do i = 1, nY
      HGrid%Y(i) = HGrid%Y0 + (i - 1)*HGrid%dY
    End Do
    HGrid%LocationsName = ' '
    HGrid%CentreName    = ' '

  ! Structured non-variable grids, specified by LocationsName, CentreName.
  Else

    HGrid%Name = Name
    If (Len_Trim(HCoordName) > MaxCharLength .or. HCoordName == ' ') Then
      Call Message('FATAL ERROR in InitHGrid', 4)
    End If
    HGrid%HCoordName = HCoordName
    HGrid%iHCoord    = 0
    HGrid%Named      = .false.
    HGrid%Wrap       = Wrap
    HGrid%nX         = nX
    HGrid%nY         = nY
    If (HGrid%nX <= 0 .or. HGrid%nY <= 0) Then
      Call Message('FATAL ERROR in InitHGrid', 4)
    End If
    HGrid%dX = dX
    HGrid%dY = dY
    Allocate(HGrid%X(HGrid%nX))
    Allocate(HGrid%Y(HGrid%nY))
    If (Len_Trim(LocationsName) > MaxCharLength .or. LocationsName == ' ') Then
      Call Message('FATAL ERROR in InitHGrid', 4)
    End If
    If (Len_Trim(CentreName) > MaxCharLength .or. CentreName == ' ') Then
      Call Message('FATAL ERROR in InitHGrid', 4)
    End If
    HGrid%LocationsName = LocationsName
    HGrid%CentreName    = CentreName

  End If

! $$ check regular cyclic grids don't overlap. gap must be >= dx - epsilon for particle
! output at least
! $$ need YWrap too (e.g. polar coords)

End Function InitHGrid

!-------------------------------------------------------------------------------------------------------------

Function InitZGrid(Name, ZCoordName, nZ, dZ, Z0, Z, AvZ, iZ, ZOnBoundaries)
! Initialises a vertical grid.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: Name
  Character(*), Intent(In)           :: ZCoordName
  Integer,      Intent(In)           :: nZ
  Real(Std),    Intent(In), Optional :: dZ
  Real(Std),    Intent(In), Optional :: Z0
  Real(Std),    Intent(In), Optional :: Z(:)
  Real(Std),    Intent(In), Optional :: AvZ(:)
  Integer,      Intent(In), Optional :: iZ(:)
  Logical,      Intent(In), Optional :: ZOnBoundaries
  ! Function result:
  Type(ZGrid_) :: InitZGrid !
  ! Locals:
  Type(ZGrid_) :: ZGrid !
  Integer      :: i

  If (Len_Trim(Name) > MaxCharLength .or. Len_Trim(Name) == 0) Then
    Call Message('FATAL ERROR in InitZGrid', 4)
  End If
  If (Len_Trim(ZCoordName) > MaxCharLength .or. Len_Trim(ZCoordName) == 0) Then
    Call Message('FATAL ERROR in InitZGrid', 4)
  End If
  If (nZ <= 0) Then
    Call Message('FATAL ERROR in InitZGrid', 4)
  End If
  If (                          &
         Present(dZ)      .and. &
         Present(Z0)      .and. &
    .not.Present(Z)       .and. &
    .not.Present(AvZ)     .and. &
    .not.Present(ZOnBoundaries) &
  ) Then
    ZGrid%Variable = .false.
  Else If (                     &
    .not.Present(dZ)      .and. &
    .not.Present(Z0)      .and. &
         Present(Z)       .and. &
         Present(AvZ)     .and. &
         Present(ZOnBoundaries) &
  ) Then
    ZGrid%Variable = .true.
  Else
    Call Message('Error in InitZGrid', 4)
  End If

  ZGrid%Name       = Name
  ZGrid%ZCoordName = ZCoordName
  ZGrid%iZCoord    = 0
  If (ZGrid%Variable .and. ZOnBoundaries) Then
    ZGrid%nZ = nZ - 1
    Allocate(ZGrid%Z(nZ - 1))
    Allocate(ZGrid%AvZ(nZ - 1))
  Else
    ZGrid%nZ = nZ
    Allocate(ZGrid%Z(nZ))
    Allocate(ZGrid%AvZ(nZ))
  End If

  ! $$ Note that dZ and Z0 are not set for a variable grid.
  ! $$ check size of input arrays
  If (ZGrid%Variable) Then
    If (ZOnBoundaries) Then
      Do i = 1, ZGrid%nZ
        ZGrid%Z(i)   = (Z(i) + Z(i+1)) / 2.0
        ZGrid%AvZ(i) = Z(i+1) - Z(i)
      End Do
    Else
      ZGrid%Z(1:nZ)   = Z(1:nZ)
      ZGrid%AvZ(1:nZ) = AvZ(1:nZ)
    End If
  Else
    ZGrid%dZ = dZ
    ZGrid%Z0 = Z0
    Do i = 1, ZGrid%nZ
      ZGrid%Z(i)   = ZGrid%Z0 + (i - 1)*ZGrid%dZ
      ZGrid%AvZ(i) = ZGrid%dZ
    End Do
  End If

  ! iZ. $$ check iZ monotone.
  If (Present(iZ)) Then
    If (Size(iZ) /= nZ) Call Message('UNEXPECTED FATAL ERROR in InitZGrid', 4)
    ZGrid%UseiZ = .true.
    Allocate(ZGrid%iZ(nZ))
    ZGrid%iZ(:) = iZ(:)
  Else
    ZGrid%UseiZ = .false.
  End If

  ! IncreasingValue should be determined here (see SetUpGrids_iCoords). $$ (0 for nZ=1)
  ZGrid%IncreasingValue = 0

  ! The monotonicity of height with index cannot be determined here and so
  ! IncreasingHeight is set to zero.
  ZGrid%IncreasingHeight = 0

  InitZGrid = ZGrid

End Function InitZGrid

!-------------------------------------------------------------------------------------------------------------

Function InitTGrid(Name, nT, dT, T0, T) Result(TGrid)
! Initialises a time grid.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: Name
  Integer,      Intent(In)           :: nT
  Type(Time_),  Intent(In), Optional :: dT
  Type(Time_),  Intent(In), Optional :: T0
  Type(Time_),  Intent(In), Optional :: T(:)
  ! Function result:
  Type(TGrid_) :: TGrid !
  ! Locals:
  Integer :: i

  If (Len_Trim(Name) > MaxCharLength .or. Len_Trim(Name) == 0) Then
    Call Message('FATAL ERROR in InitTGrid', 4)
  End If
  If (nT <= 0) Then
    Call Message('FATAL ERROR in InitTGrid', 4)
  End If
  If (                     &
         Present(dT) .and. &
         Present(T0) .and. &
    .not.Present(T)        &
  ) Then
    TGrid%Variable = .false.
  Else If (                &
    .not.Present(dT) .and. &
    .not.Present(T0) .and. &
         Present(T)        &
  ) Then
    TGrid%Variable = .true.
  Else
    Call Message('FATAL ERROR in InitTGrid', 4)
  End If

 !$$ check time types re interval or clock times. dT = interval time

  TGrid%Name = Name
  TGrid%nT   = nT

  !Allocate(TGrid%T(nT))

  If (TGrid%Variable) Then
  !  Do i = 1, nT
  !    TGrid%T(i) = Time2ShortTime(T(i))
  !  End Do
  Else
    TGrid%dT = dT
    TGrid%T0 = T0
    TGrid%sdT = Time2ShortTime(dT)
    TGrid%sT0 = Time2ShortTime(T0)
  !  TGrid%T(1) = Time2ShortTime(TGrid%T0)
  !  Do i = 2, TGrid%nT
  !    TGrid%T(i) = TGrid%T(i - 1) + Time2ShortTime(TGrid%dT)
  !  End Do
  End If

End Function InitTGrid

!-------------------------------------------------------------------------------------------------------------

Function InitDGrid(Name, nD, dD, D0, D, Floating, Geometric, Zero) Result(DGrid)
! Initialises a data grid.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: Name
  Integer,      Intent(In)           :: nD
  Real(Std),    Intent(In), Optional :: dD
  Real(Std),    Intent(In), Optional :: D0
  Real(Std),    Intent(In), Optional :: D(:)
  Logical,      Intent(In), Optional :: Floating
  Logical,      Intent(In), Optional :: Geometric
  Logical,      Intent(In), Optional :: Zero
  ! Function result:
  Type(DGrid_) :: DGrid !
  ! Locals:
  Integer :: i

  If (Len_Trim(Name) > MaxCharLength .or. Len_Trim(Name) == 0) Then
    Call Message('FATAL ERROR in InitDGrid', 4)
  End If
  If (nD <= 0) Then
    Call Message('FATAL ERROR in InitDGrid', 4)
  End If
  If (                            &
         Present(Floating)  .and. &
         Present(Geometric) .and. &
         Present(Zero)      .and. &
         Present(dD)        .and. &
         Present(D0)        .and. &
    .not.Present(D)               &
  ) Then
    DGrid%Variable = .false.
  Else If (                       &
    .not.Present(Floating)  .and. &
    .not.Present(Geometric) .and. &
    .not.Present(Zero)      .and. &
    .not.Present(dD)        .and. &
    .not.Present(D0)        .and. &
         Present(D)               &
  ) Then
    DGrid%Variable = .true.
  Else
    Call Message('FATAL ERROR in InitDGrid', 4)
  End If

  DGrid%Name = Name
  DGrid%nD   = nD ! $$ check nD >= 1 - all grids must have one point at least

  If (DGrid%Variable) Then
    DGrid%Floating  = .false.
    DGrid%Geometric = .false.
    DGrid%Zero      = .false.
    Allocate(DGrid%D(nD))
    Do i = 1, nD
      DGrid%D(i) = D(i)
    End Do
  Else
    If (.not.Geometric .and. Zero) Then
      Call Message('FATAL ERROR in InitDGrid', 4)
    End If
    DGrid%Floating  = Floating
    DGrid%Geometric = Geometric
    DGrid%Zero      = Zero
    DGrid%dD        = dD
    DGrid%D0        = D0
  End If

End Function InitDGrid

!-------------------------------------------------------------------------------------------------------------

Subroutine AddHGrid(HGrid, Grids)
! Adds a horizontal grid to a collection of grids.

  Implicit None
  ! Argument list:
  Type(HGrid_),  Intent(In)    :: HGrid
  Type(Grids_),  Intent(InOut) :: Grids
  ! Locals:
  Integer :: i !

  Do i = 1, Grids%nHGrids
    If (HGrid%Name .CIEq. Grids%HGrids(i)%Name) Then
      If (HGrid == Grids%HGrids(i)) Then
        Return
      Else
        Call Message('FATAL ERROR: two different horizontal grids with the same name', 3)
      End If
    End If
  End Do
  If (Grids%nHGrids >= MaxHGrids) Then
    Call Message('FATAL ERROR: too many horizontal grids', 3)
  End If
  Grids%nHGrids               = Grids%nHGrids + 1
  Grids%HGrids(Grids%nHGrids) = HGrid

End Subroutine AddHGrid

!-------------------------------------------------------------------------------------------------------------

Subroutine AddZGrid(ZGrid, Grids)
! Adds a vertical grid to a collection of grids.

  Implicit None
  ! Argument list:
  Type(ZGrid_),  Intent(In)    :: ZGrid
  Type(Grids_),  Intent(InOut) :: Grids
  ! Locals:
  Integer :: i !

  Do i = 1, Grids%nZGrids
    If (ZGrid%Name .CIEq. Grids%ZGrids(i)%Name) Then
      If (ZGrid == Grids%ZGrids(i)) Then
        Return
      Else
        Call Message('FATAL ERROR: two different vertical grids with the same name', 3)
      End If
    End If
  End Do
  If (Grids%nZGrids >= MaxZGrids) Then
    Call Message('FATAL ERROR: too many vertical grids', 3)
  End If
  Grids%nZGrids               = Grids%nZGrids + 1
  Grids%ZGrids(Grids%nZGrids) = ZGrid

End Subroutine AddZGrid

!-------------------------------------------------------------------------------------------------------------

Subroutine AddTGrid(TGrid, Grids)
! Adds a time grid to a collection of grids.

  Implicit None
  ! Argument list:
  Type(TGrid_),  Intent(In)    :: TGrid
  Type(Grids_),  Intent(InOut) :: Grids
  ! Locals:
  Integer :: i !

  Do i = 1, Grids%nTGrids
    If (TGrid%Name .CIEq. Grids%TGrids(i)%Name) Then
      If (TGrid == Grids%TGrids(i)) Then
        Return
      Else
        Call Message('FATAL ERROR: two different time grids with the same name', 3)
      End If
    End If
  End Do
  If (Grids%nTGrids >= MaxTGrids) Then
    Call Message('FATAL ERROR: too many time grids', 3)
  End If
  Grids%nTGrids               = Grids%nTGrids + 1
  Grids%TGrids(Grids%nTGrids) = TGrid

End Subroutine AddTGrid

!-------------------------------------------------------------------------------------------------------------

Subroutine AddDGrid(DGrid, Grids)
! Adds a data grid to a collection of grids.

  Implicit None
  ! Argument list:
  Type(DGrid_),  Intent(In)    :: DGrid
  Type(Grids_),  Intent(InOut) :: Grids
  ! Locals:
  Integer :: i !

  Do i = 1, Grids%nDGrids
    If (DGrid%Name .CIEq. Grids%DGrids(i)%Name) Then
      If (DGrid == Grids%DGrids(i)) Then
        Return
      Else
        Call Message('FATAL ERROR: two different data grids with the same name', 3)
      End If
    End If
  End Do
  If (Grids%nDGrids >= MaxDGrids) Then
    Call Message('FATAL ERROR: too many data grids', 3)
  End If
  Grids%nDGrids               = Grids%nDGrids + 1
  Grids%DGrids(Grids%nDGrids) = DGrid

End Subroutine AddDGrid

!-------------------------------------------------------------------------------------------------------------

Function HGridEq(HGrid1, HGrid2)
! Tests for equality of horizontal grids.

  Implicit None
  ! Argument list:
  Type(HGrid_), Intent(In) :: HGrid1 !
  Type(HGrid_), Intent(In) :: HGrid2 !
  ! Function result:
  Logical :: HGridEq !

  HGridEq = (HGrid1%Name       .CIEq. HGrid2%Name      ) .and. &
            (HGrid1%HCoordName .CIEq. HGrid2%HCoordName) .and. &
            (HGrid1%Variable   .Eqv.  HGrid2%Variable  ) .and. &
            (HGrid1%Wrap       .Eqv.  HGrid2%Wrap      ) .and. &
             HGrid1%nX           ==   HGrid2%nX          .and. &
             HGrid1%nY           ==   HGrid2%nY

  If (HGrid1%Variable) Then ! $$ potential for array exception here if grids are different.
                            ! check similar places in other grid types.
    HGridEq = HGridEq .and. All(HGrid1%X(1:HGrid1%nX) == HGrid2%X(1:HGrid1%nX))
    HGridEq = HGridEq .and. All(HGrid1%Y(1:HGrid1%nY) == HGrid2%Y(1:HGrid1%nY))
  Else
    HGridEq = HGridEq                .and. &
              HGrid1%dX == HGrid2%dX .and. &
              HGrid1%dY == HGrid2%dY .and. &
              HGrid1%X0 == HGrid2%X0 .and. &
              HGrid1%Y0 == HGrid2%Y0
  End If

End Function HGridEq

!-------------------------------------------------------------------------------------------------------------

Function ZGridEq(ZGrid1, ZGrid2)
! Tests for equality of vertical grids.

  Implicit None
  ! Argument list:
  Type(ZGrid_), Intent(In) :: ZGrid1 !
  Type(ZGrid_), Intent(In) :: ZGrid2 !
  ! Function result:
  Logical :: ZGridEq !

  ZGridEq = (ZGrid1%Name       .CIEq. ZGrid2%Name      ) .and. &
            (ZGrid1%ZCoordName .CIEq. ZGrid2%ZCoordName) .and. &
            (ZGrid1%Variable   .Eqv.  ZGrid2%Variable  ) .and. &
            (ZGrid1%nZ           ==   ZGrid2%nZ        )

  If (ZGrid1%Variable) Then
    ZGridEq = ZGridEq .and. All(ZGrid1%Z(1:ZGrid1%nZ) == ZGrid2%Z(1:ZGrid1%nZ))
  Else
    ZGridEq = ZGridEq                .and. &
              ZGrid1%dZ == ZGrid2%dZ .and. &
              ZGrid1%Z0 == ZGrid2%Z0
  End If

  If (ZGrid1%UseiZ) Then
    ZGridEq = ZGridEq .and. All(ZGrid1%iZ(1:ZGrid1%nZ) == ZGrid2%iZ(1:ZGrid1%nZ))
  End If

End Function ZGridEq

!-------------------------------------------------------------------------------------------------------------

Function TGridEq(TGrid1, TGrid2)
! Tests for equality of time grids.

  Implicit None
  ! Argument list:
  Type(TGrid_), Intent(In) :: TGrid1 !
  Type(TGrid_), Intent(In) :: TGrid2 !
  ! Function result:
  Logical :: TGridEq !

  TGridEq = (TGrid1%Name     .CIEq. TGrid2%Name    ) .and. &
            (TGrid1%Variable .Eqv.  TGrid2%Variable) .and. &
            (TGrid1%nT         ==   TGrid2%nT      )

  If (TGrid1%Variable) Then
  !  TGridEq = TGridEq .and. All(TGrid1%T(1:TGrid1%nT) == TGrid2%T(1:TGrid1%nT))
  Else
    TGridEq = TGridEq                .and. &
              TGrid1%dT == TGrid2%dT .and. &
              TGrid1%T0 == TGrid2%T0
  End If

End Function TGridEq

!-------------------------------------------------------------------------------------------------------------

Function DGridEq(DGrid1, DGrid2)
! Tests for equality of data grids.

  Implicit None
  ! Argument list:
  Type(DGrid_), Intent(In) :: DGrid1 !
  Type(DGrid_), Intent(In) :: DGrid2 !
  ! Function result:
  Logical :: DGridEq !

  DGridEq = (DGrid1%Name      .CIEq. DGrid2%Name     ) .and. &
            (DGrid1%Variable  .Eqv.  DGrid2%Variable ) .and. &
            (DGrid1%Floating  .Eqv.  DGrid2%Floating ) .and. &
            (DGrid1%Geometric .Eqv.  DGrid2%Geometric) .and. &
            (DGrid1%Zero      .Eqv.  DGrid2%Zero     ) .and. &
            (DGrid1%nD          ==   DGrid2%nD       )

  If (DGrid1%Variable) Then
    DGridEq = DGridEq .and. All(DGrid1%D(1:DGrid1%nD) == DGrid2%D(1:DGrid1%nD))
  Else
    DGridEq = DGridEq                .and. &
              DGrid1%dD == DGrid2%dD .and. &
              DGrid1%D0 == DGrid2%D0
  End If

End Function DGridEq

!-------------------------------------------------------------------------------------------------------------

Function FindHGridIndex(Name, Grids)
! .

  Implicit None
  ! Argument list:
  Character(*), Intent(In)  :: Name
  Type(Grids_), Intent(In)  :: Grids
  ! Function result:
  Integer :: FindHGridIndex
  ! Locals:
  Integer :: i

  Do i = 1, Grids%nHGrids
    If (Name .CIEq. Grids%HGrids(i)%Name) Then
      FindHGridIndex = i
      Return
    End If
  End Do

  Call Message(                              &
         'FATAL ERROR: horizontal grid "' // &
         Trim(Name)                       // &
         '" not found',                      &
         3                                   &
       )

End Function FindHGridIndex

!-------------------------------------------------------------------------------------------------------------

Function FindZGridIndex(Name, Grids)
! .

  Implicit None
  ! Argument list:
  Character(*), Intent(In)  :: Name
  Type(Grids_), Intent(In)  :: Grids
  ! Function result:
  Integer :: FindZGridIndex
  ! Locals:
  Integer :: i

  Do i = 1, Grids%nZGrids
    If (Name .CIEq. Grids%ZGrids(i)%Name) Then
      FindZGridIndex = i
      Return
    End If
  End Do

  Call Message(                            &
         'FATAL ERROR: vertical grid "' // &
         Trim(Name)                     // &
         '" not found',                    &
         3                                 &
       )

End Function FindZGridIndex

!-------------------------------------------------------------------------------------------------------------

Function FindTGridIndex(Name, Grids)
! .

  Implicit None
  ! Argument list:
  Character(*), Intent(In)  :: Name
  Type(Grids_), Intent(In)  :: Grids
  ! Function result:
  Integer :: FindTGridIndex
  ! Locals:
  Integer :: i

  Do i = 1, Grids%nTGrids
    If (Name .CIEq. Grids%TGrids(i)%Name) Then
      FindTGridIndex = i
      Return
    End If
  End Do

  Call Message(                        &
         'FATAL ERROR: time grid "' // &
         Trim(Name)                 // &
         '" not found',                &
         3                             &
       )

End Function FindTGridIndex

!-------------------------------------------------------------------------------------------------------------

Function FindDGridIndex(Name, Grids)
! .

  Implicit None
  ! Argument list:
  Character(*), Intent(In)  :: Name
  Type(Grids_), Intent(In)  :: Grids
  ! Function result:
  Integer :: FindDGridIndex
  ! Locals:
  Integer :: i

  Do i = 1, Grids%nDGrids
    If (Name .CIEq. Grids%DGrids(i)%Name) Then
      FindDGridIndex = i
      Return
    End If
  End Do

  Call Message(                        &
         'FATAL ERROR: data grid "' // &
         Trim(Name)                 // &
         '" not found',                &
         3                             &
       )

End Function FindDGridIndex

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpGrids_CoordsEtc(Locationses, Coords, Grids)
! Sets up %IncreasingValue in Grids and grid centre and X/YCycle.

  Implicit None
  ! Argument list:
  Type(Locationses_), Intent(In)            :: Locationses
  Type(Coords_),      Intent(In)            :: Coords
  Type(Grids_),       Intent(InOut), Target :: Grids
  ! Locals:
  Integer               :: i
  Integer               :: j
  Integer               :: iLocations
  Integer               :: iCentre
  Integer               :: iHCoord1
  Integer               :: iHCoord2
  Integer               :: iZCoord
  Real(Std)             :: Point(2)
  Type(HGrid_), Pointer :: HGrid

  Do i = 1, Grids%nHGrids

    HGrid => Grids%HGrids(i)

    If (HGrid%CentreName /= ' ') Then

      iLocations = FindLocationsIndex(HGrid%LocationsName, Locationses)
      iCentre    = FindLocationIndex(                    &
                     HGrid%CentreName,                   &
                     Locationses%Locationses(iLocations) &
                   )
      iHCoord1   = FindHCoordIndex(                                            &
                     Locationses%Locationses(iLocations)%HCoordNames(iCentre), &
                     Coords                                                    &
                   )
      iHCoord2   = FindHCoordIndex(HGrid%HCoordName, Coords)
      Point      = ConvertH(                                           &
                     Coords%HCoords(iHCoord1),                         &
                     Coords%HCoords(iHCoord2),                         &
                     (/                                                &
                       Locationses%Locationses(iLocations)%X(iCentre), &
                       Locationses%Locationses(iLocations)%Y(iCentre)  &
                     /)                                                &
                   )

      HGrid%X0 = Point(1) - Real(HGrid%nX - 1, Std)*HGrid%dX/2.0
      HGrid%Y0 = Point(2) - Real(HGrid%nY - 1, Std)*HGrid%dY/2.0
      HGrid%XCentre = Point(1)
      HGrid%YCentre = Point(2)
      Do j = 1, HGrid%nX
        HGrid%X(j) = HGrid%X0 + (j - 1)*HGrid%dX
      End Do
      Do j = 1, HGrid%nY
        HGrid%Y(j) = HGrid%Y0 + (j - 1)*HGrid%dY
      End Do

    End If

    ! XCycle and YCycle.
    iHCoord2   = FindHCoordIndex(HGrid%HCoordName, Coords)
    HGrid%XCycle = Coords%HCoords(iHCoord2)%XCycle
    HGrid%YCycle = Coords%HCoords(iHCoord2)%YCycle

  End Do

  Do i = 1, Grids%nZGrids
    iZCoord = FindZCoordIndex(Grids%ZGrids(i)%ZCoordName, Coords)
    ! $$ IncreasingValue should have been set in InitZCoord and IncreasingHeight should
    ! be determined here from IncreasingValue and %IncreasingUpwards. Instead we
    ! assume that height and index always increase together. $$
    Grids%ZGrids(i)%IncreasingHeight = 1
    If (Coords%ZCoords(iZCoord)%IncreasingUpwards) Then
      Grids%ZGrids(i)%IncreasingValue = Grids%ZGrids(i)%IncreasingHeight
    Else
      Grids%ZGrids(i)%IncreasingValue = - Grids%ZGrids(i)%IncreasingHeight
    End If
  End Do

End Subroutine SetUpGrids_CoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpGrids_iCoords(Coords, Grids)
! Sets up indices in Grids for referring to coord systems.

  Implicit None
  ! Argument list:
  Type(Coords_), Intent(In)    :: Coords
  Type(Grids_),  Intent(InOut) :: Grids
  ! Locals:
  Integer :: i

  Do i = 1, Grids%nHGrids
    Grids%HGrids(i)%iHCoord = FindHCoordIndex(Grids%HGrids(i)%HCoordName, Coords)
  End Do

  Do i = 1, Grids%nZGrids
    Grids%ZGrids(i)%iZCoord = FindZCoordIndex(Grids%ZGrids(i)%ZCoordName, Coords)
  End Do

End Subroutine SetUpGrids_iCoords

!-------------------------------------------------------------------------------------------------------------

Function TInTGrid(TGrid, i)
! Calculates the i-th time T(i) = T0 + dT*(i-1), i = 1,...,n in the time grid TGrid.

  Implicit None
  ! Argument list:
  Type(TGrid_), Intent(In) :: TGrid
  Integer,      Intent(In) :: i
  ! Function result:
  Type(ShortTime_) :: TInTGrid ! The i-th time in TGrid.

 !$$ check Time Types here

  If (i < 1 .or. i > TGrid%nT) Then
    Call Message('FATAL ERROR in TInTGrid: time index outside TGrid bounds', 4)
  Else
    TInTGrid = TGrid%sT0 + TGrid%sdT * (i-1)
  End If

End Function TInTGrid

!-------------------------------------------------------------------------------------------------------------

Subroutine FirstTAfterT(TGrid, T, Strict, FirstT, iT)
! Calculates the first time in a TGrid at or after a given time (inf for an infinite grid).

  Implicit None
  ! Argument list:
  Type(TGrid_),     Intent(In)            :: TGrid
  Type(ShortTime_), Intent(In)            :: T      ! Must not be infinite future.
  Logical,          Intent(In)            :: Strict ! means first time > T. Otherwise >=.
  Type(ShortTime_), Intent(Out)           :: FirstT ! The first time in TGrid after T. inf if no such time
  Integer,          Intent(Out), Optional :: iT     ! Index in grid. huge if no such time

  If (IsInfPast(T)) Then
    FirstT = TGrid%sT0 ! $$ is dt always > 0 (in run direction) and nT >= 1
    If (Present(iT)) iT = 1
  Else If (IsInfFuture(T)) Then
    FirstT = InfFutureShortTime(IsTimeInterval(TGrid%sT0))
    If (Present(iT)) iT = Huge(iT)/3 ! $$ make Huge/3 a function or param? Use elsewhere?
  Else
    FirstT = Round(T, TGrid%sT0, TGrid%sdT, Up = .true.)
    If (Strict) Then
      If (FirstT == T) FirstT = FirstT + TGrid%sdT
    End If
    If (FirstT < TGrid%sT0) FirstT = TGrid%sT0
    If (Present(iT)) Then
      iT = (FirstT - TGrid%sT0) / TGrid%sdT + 1
    End If
    If (FirstT > TGrid%sT0 + TGrid%sdT * (TGrid%nT - 1)) Then
      FirstT = InfFutureShortTime(IsTimeInterval(TGrid%sT0))
      If (Present(iT)) iT = Huge(iT)/3
    End If
  End If

End Subroutine FirstTAfterT

!-------------------------------------------------------------------------------------------------------------

Subroutine LastTBeforeT(TGrid, T, Strict, LastT, iT)
! Calculates the last time in a TGrid at or before a given time (sup for an infinite grid).

  Implicit None
  ! Argument list:
  Type(TGrid_),     Intent(In)            :: TGrid
  Type(ShortTime_), Intent(In)            :: T      ! Must not be infinite past.
  Logical,          Intent(In)            :: Strict ! means last time < T. Otherwise <=.
  Type(ShortTime_), Intent(Out)           :: LastT  ! The last time in TGrid before T. -inf if no such time
  Integer,          Intent(Out), Optional :: iT     ! Index in grid. -Huge if no such time.

  If (IsInfFuture(T)) Then
    LastT = TGrid%sT0 + TGrid%sdT * (TGrid%nT - 1) ! $$ is dt always > 0 (in run direction)
    If (Present(iT)) iT = TGrid%nT
  Else If (IsInfPast(T)) Then
    LastT = InfPastShortTime(IsTimeInterval(TGrid%sT0))
    If (Present(iT)) iT = -Huge(iT)/3
  Else
    LastT = Round(T, TGrid%sT0, TGrid%sdT, Up = .false.)
    If (Strict) Then
      If (LastT == T) LastT = LastT - TGrid%sdT
    End If
    If (LastT > TGrid%sT0 + TGrid%sdT * (TGrid%nT - 1)) LastT = TGrid%sT0 + TGrid%sdT * (TGrid%nT - 1)
    If (Present(iT)) Then
      iT = (LastT - TGrid%sT0) / TGrid%sdT + 1
    End If
    If (LastT < TGrid%sT0) Then
      LastT = InfPastShortTime(IsTimeInterval(TGrid%sT0))
      If (Present(iT)) iT = -Huge(iT)/3
    End If
  End If

End Subroutine LastTBeforeT

!-------------------------------------------------------------------------------------------------------------

Function nDInDGrid(DGrid)
!

  Implicit None
  ! Argument list:
  Type(DGrid_), Intent(In) :: DGrid
  ! Function result:
  Integer :: nDInDGrid !

  nDInDGrid = DGrid%nD

End Function nDInDGrid

!-------------------------------------------------------------------------------------------------------------

Function FloatingDGrid(DGrid)
!

  Implicit None
  ! Argument list:
  Type(DGrid_), Intent(In) :: DGrid
  ! Function result:
  Logical :: FloatingDGrid !

  FloatingDGrid = DGrid%Floating

End Function FloatingDGrid

!-------------------------------------------------------------------------------------------------------------

Function ZeroDGrid(DGrid) ! $$ Needed?
!

  Implicit None
  ! Argument list:
  Type(DGrid_), Intent(In) :: DGrid
  ! Function result:
  Logical :: ZeroDGrid !

  ZeroDGrid = DGrid%Zero

End Function ZeroDGrid

!-------------------------------------------------------------------------------------------------------------

Function DInDGrid(DGrid, i, s)
! Calculates the i-th time T(i) = T0 + dT*(i-1), i = 1,...,n in the time grid TGrid.

  Implicit None
  ! Argument list:
  Type(DGrid_), Intent(In) :: DGrid
  Integer,      Intent(In) :: i
  Integer,      Intent(In) :: s
  ! Function result:
  Real(Std) :: DInDGrid ! The i-th time in TGrid.

  If (.not.DGrid%Floating .and. s /= 0) Then
    Call Message('FATAL ERROR in DInDGrid: index outside DGrid bounds', 3)
  End If
  If (i < 1 .or. i > DGrid%nD) Then
    Call Message('FATAL ERROR in DInDGrid: index outside DGrid bounds', 3)
  End If

  If (DGrid%Variable) Then
    DInDGrid = DGrid%D(i)
  Else
    If (DGrid%Geometric) Then
      If (DGrid%Zero) Then
        If (i == 1) Then
          DInDGrid = 0.0
        Else
          DInDGrid = DGrid%D0 * DGrid%dD ** (i + s - 2)
        End If
      Else
        DInDGrid = DGrid%D0 * DGrid%dD ** (i + s - 1)
      End If
    Else
      DInDGrid = DGrid%D0 + DGrid%dD * (i + s - 1)
    End If
  End If

End Function DInDGrid

!-------------------------------------------------------------------------------------------------------------

Function ShiftIndexInDGrid(DGrid, D) Result(s)
! Finds the shift value s needed to represent D.

  Implicit None
  ! Argument list:
  Type(DGrid_), Intent(In)  :: DGrid
  Real(Std),    Intent(In)  :: D     !
  ! Function result:
  Integer :: s

  If (.not.DGrid%Floating) Then
    Call Message('UNEXPECTED FATAL ERROR in ShiftIndexInDGrid', 3)
  End If

  If (DGrid%Geometric) Then
    If (D <= 0.0) Then
      s = -Huge(s)
    Else
      s = Ceiling(Log(D / DGrid%D0) / Log(DGrid%dD)) - DGrid%nD + 1
      If (DGrid%Zero) s = s + 1
    End If
  Else
    s = Ceiling(   (D - DGrid%D0) /     DGrid%dD ) - DGrid%nD + 1
  End If

  ! Rounding check
  If (DInDGrid(DGrid, DGrid%nD, s) < D) s = s + 1

End Function ShiftIndexInDGrid

!-------------------------------------------------------------------------------------------------------------

Subroutine FirstDAfterD(DGrid, s, D, FirstD, iD)
! Calculates the first value in a DGrid at or after a given value.

! $$ this is just a rough sketch for future development - not yet used.

  Implicit None
  ! Argument list:
  Type(DGrid_), Intent(In)            :: DGrid
  Integer,      Intent(In)            :: s
  Real(Std),    Intent(In)            :: D      !
  Real(Std),    Intent(Out)           :: FirstD ! The first time in TGrid after T. inf if no such time
  Integer,      Intent(Out), Optional :: iD     ! Index in grid. huge if no such value
  ! Locals:
  Integer :: i

  If (DGrid%Variable) Then
    i = 1 ! $$
  Else
    If (DGrid%Geometric) Then
      i = Ceiling(Log(D / DGrid%D0) / Log(DGrid%dD)) - s + 1
    Else
      i = Ceiling(   (D - DGrid%D0) /     DGrid%dD ) - s + 1
    End If
  End If

  FirstD = DInDGrid(DGrid, i, s)

  ! Rounding check
  If (FirstD < D) Then
    i = i + 1
    FirstD = DInDGrid(DGrid, i, s)
  End If

  If (FirstD < DGrid%D0) FirstD = DGrid%D0
  If (Present(iD)) Then
    iD = (FirstD - DGrid%D0) / DGrid%dD + 1
  End If
  If (FirstD > DGrid%D0 + DGrid%dD * (DGrid%nD - 1)) Then
    If (Present(iD)) iD = Huge(iD)
  End If

End Subroutine FirstDAfterD

!-------------------------------------------------------------------------------------------------------------

Subroutine InitGHCoeffs(GHCoeffs)

! Note this must be called before using interpolation coefficients for first time.

  Implicit None
  ! Argument List:
  Type(GHCoeffs_), Intent(Out) :: GHCoeffs ! Interpolation coefficients.

  ! Set pointers to null.
  GHCoeffs%iX => Null()
  GHCoeffs%iY => Null()
  GHCoeffs%X  => Null()
  GHCoeffs%Y  => Null()

End Subroutine InitGHCoeffs

!-------------------------------------------------------------------------------------------------------------

Subroutine GetGHCoeffs(FromHGrid, ToHGrid, Derivative, GHCoeffs)
! Calculates the interpolation coefficients needed to evaluate a quantity or its
! derivative when the quantity is defined on one horizontal grid and the result is
! required at points on another horizontal grid.

! Will work for points of ToHGrid outside FromHGrid, giving the values corresponding
! to the nearest boundary point. Assumes (i) that the grids have the same spacing,
! (ii) that, at every edge, the grids are aligned or off-set by 1/2 the grid spacing,
! (iii) that both grids are wrapped or both are not, and (iv) that the grids have at
! least two points in each direction.

  Implicit None
  ! Argument List:
  Type(HGrid_),    Intent(In)    :: FromHGrid  ! Grid on which the quantity in question is defined.
  Type(HGrid_),    Intent(In)    :: ToHGrid    ! Grid on which the result is required.
  Character(1),    Intent(In)    :: Derivative ! blank, X or Y - calculates derivative
                                               ! in direction of increasing grid point index
                                               ! with respect to the coordinate
  Type(GHCoeffs_), Intent(InOut) :: GHCoeffs   ! Interpolation coefficients.
  ! Locals:
  Integer   :: ii(ToHGrid%nX, 2)
  Integer   :: jj(ToHGrid%nY, 2)
  Real(Std) :: X(ToHGrid%nX, 2)
  Real(Std) :: Y(ToHGrid%nY, 2)
  Integer   :: i
  Integer   :: j
  Integer   :: nX
  Integer   :: nY
  Integer   :: n

  Do i = 1, ToHGrid%nX

    If (FromHGrid%X0 - ToHGrid%X0 == 0.0) Then ! $$ Precision problems ?
      If (Derivative == 'X') Then
        ii(i, 1) = i - 1
        ii(i, 2) = i + 1
        nX       = 2
      Else
        ii(i, 1) = i
        nX       = 1
      End If
    Else If ((FromHGrid%X0 - ToHGrid%X0) * Sign(1.0, ToHGrid%dX) > 0.0) Then
      ii(i, 1) = i - 1
      ii(i, 2) = i
      nX       = 2
    Else If ((FromHGrid%X0 - ToHGrid%X0) * Sign(1.0, ToHGrid%dX) < 0.0) Then
      ii(i, 1) = i
      ii(i, 2) = i + 1
      nX       = 2
    End If

    Do n = 1, nX
      If (ii(i, n) < 1) Then
        If (FromHGrid%Wrap) Then
          ii(i, n) = FromHGrid%nX
        Else
          ii(i, n) = 1
        End If
      Else If (ii(i, n) > FromHGrid%nX) Then
        If (FromHGrid%Wrap) Then
          ii(i, n) = 1
        Else
          ii(i, n) = FromHGrid%nX
        End If
      End If
    End Do

    If (Derivative == 'X' .and. ii(i, 1) == ii(i, 2)) Then
      If (ii(i, 1) == 1) Then
        ii(i, 2) = 2
      Else If (ii(i, 1) == FromHGrid%nX) Then
        ii(i, 1) = FromHGrid%nX - 1
      End If
    End If

    If (Derivative == 'X') Then
      X(i, 1) = -1.0 / (FromHGrid%dX * Real(ii(i, 2) - ii(i, 1), Std))
      X(i, 2) =  1.0 / (FromHGrid%dX * Real(ii(i, 2) - ii(i, 1), Std))
    Else
      X(i, 1:nX) = 1.0 / Real(nX, Std)
    End If

  End Do

  Do j = 1, ToHGrid%nY

    If (FromHGrid%Y0 - ToHGrid%Y0 == 0.0) Then
      If (Derivative == 'Y') Then
        jj(j, 1) = j - 1
        jj(j, 2) = j + 1
        nY       = 2
      Else
        jj(j, 1) = j
        nY       = 1
      End If
    Else If ((FromHGrid%Y0 - ToHGrid%Y0) * Sign(1.0, ToHGrid%dY) > 0.0) Then
      jj(j, 1) = j - 1
      jj(j, 2) = j
      nY       = 2
    Else If ((FromHGrid%Y0 - ToHGrid%Y0) * Sign(1.0, ToHGrid%dY) < 0.0) Then
      jj(j, 1) = j
      jj(j, 2) = j + 1
      nY       = 2
    End If

    Do n = 1, nY
      If (jj(j, n) < 1) Then
        jj(j, n) = 1
      Else If (jj(j, n) > FromHGrid%nY) Then
        jj(j, n) = FromHGrid%nY
      End If
    End Do

    If (Derivative == 'Y' .and. jj(j, 1) == jj(j, 2)) Then
      If (jj(j, 1) == 1) Then
        jj(j, 2) = 2
      Else If (jj(j, 1) == FromHGrid%nX) Then
        jj(j, 1) = FromHGrid%nX - 1
      End If
    End If

    If (Derivative == 'Y') Then
      Y(j, 1) = -1.0 / (FromHGrid%dY * Real(jj(j, 2) - jj(j, 1), Std))
      Y(j, 2) =  1.0 / (FromHGrid%dY * Real(jj(j, 2) - jj(j, 1), Std))
    Else
      Y(j, 1:nY) = 1.0 / Real(nY, Std)
    End If

  End Do

  If (Associated(GHCoeffs%iX)) Then
    If (Size(GHCoeffs%iX, 1) /= ToHGrid%nX) Deallocate(GHCoeffs%iX)
  End If
  If (.not.Associated(GHCoeffs%iX)) Then
    Allocate(GHCoeffs%iX(ToHGrid%nX, 4))
  End If

  If (Associated(GHCoeffs%iY)) Then
    If (Size(GHCoeffs%iY, 1) /= ToHGrid%nY) Deallocate(GHCoeffs%iY)
  End If
  If (.not.Associated(GHCoeffs%iY)) Then
    Allocate(GHCoeffs%iY(ToHGrid%nY, 4))
  End If

  If (Associated(GHCoeffs%X)) Then
    If (Size(GHCoeffs%X, 1) /= ToHGrid%nX) Deallocate(GHCoeffs%X)
  End If
  If (.not.Associated(GHCoeffs%X)) Then
    Allocate(GHCoeffs%X(ToHGrid%nX, 4))
  End If

  If (Associated(GHCoeffs%Y)) Then
    If (Size(GHCoeffs%Y, 1) /= ToHGrid%nY) Deallocate(GHCoeffs%Y)
  End If
  If (.not.Associated(GHCoeffs%Y)) Then
    Allocate(GHCoeffs%Y(ToHGrid%nY, 4))
  End If

  If (nX == 2 .and. nY == 2) Then
    GHCoeffs%iX(:, 1) = ii(:, 1)
    GHCoeffs%iX(:, 2) = ii(:, 1)
    GHCoeffs%iX(:, 3) = ii(:, 2)
    GHCoeffs%iX(:, 4) = ii(:, 2)
    GHCoeffs%iY(:, 1) = jj(:, 1)
    GHCoeffs%iY(:, 2) = jj(:, 2)
    GHCoeffs%iY(:, 3) = jj(:, 1)
    GHCoeffs%iY(:, 4) = jj(:, 2)
    GHCoeffs%X (:, 1) = X (:, 1)
    GHCoeffs%X (:, 2) = X (:, 1)
    GHCoeffs%X (:, 3) = X (:, 2)
    GHCoeffs%X (:, 4) = X (:, 2)
    GHCoeffs%Y (:, 1) = Y (:, 1)
    GHCoeffs%Y (:, 2) = Y (:, 2)
    GHCoeffs%Y (:, 3) = Y (:, 1)
    GHCoeffs%Y (:, 4) = Y (:, 2)
  Else If (nX == 2) Then
    GHCoeffs%iX(:, 1) = ii(:, 1)
    GHCoeffs%iX(:, 2) = ii(:, 2)
    GHCoeffs%iY(:, 1) = jj(:, 1)
    GHCoeffs%iY(:, 2) = jj(:, 1)
    GHCoeffs%X (:, 1) = X (:, 1)
    GHCoeffs%X (:, 2) = X (:, 2)
    GHCoeffs%Y (:, 1) = Y (:, 1)
    GHCoeffs%Y (:, 2) = Y (:, 1)
  Else If (nY == 2) Then
    GHCoeffs%iX(:, 1) = ii(:, 1)
    GHCoeffs%iX(:, 2) = ii(:, 1)
    GHCoeffs%iY(:, 1) = jj(:, 1)
    GHCoeffs%iY(:, 2) = jj(:, 2)
    GHCoeffs%X (:, 1) = X (:, 1)
    GHCoeffs%X (:, 2) = X (:, 1)
    GHCoeffs%Y (:, 1) = Y (:, 1)
    GHCoeffs%Y (:, 2) = Y (:, 2)
  Else
    GHCoeffs%iX(:, 1) = ii(:, 1)
    GHCoeffs%iY(:, 1) = jj(:, 1)
    GHCoeffs%X (:, 1) = X (:, 1)
    GHCoeffs%Y (:, 1) = Y (:, 1)
  End If
  GHCoeffs%n     = nX * nY
  GHCoeffs%Valid = .true.

End Subroutine GetGHCoeffs

!-------------------------------------------------------------------------------------------------------------

Function GInterpXY(Field, i, j, GHCoeffs)

  Implicit None
  ! Argument List:
  Real(Std),       Intent(In) :: Field(:, :)
  Integer,         Intent(In) :: i
  Integer,         Intent(In) :: j
  Type(GHCoeffs_), Intent(In) :: GHCoeffs
  ! Function result:
  Real(Std) :: GInterpXY
  ! Locals:
  Integer :: m

  GInterpXY = 0.0
  Do m = 1, GHCoeffs%n
    GInterpXY = GInterpXY +                                   &
                Field(GHCoeffs%iX(i, m), GHCoeffs%iY(j, m)) * &
                GHCoeffs%X(i, m) * GHCoeffs%Y(j, m)
  End Do

End Function GInterpXY

!-------------------------------------------------------------------------------------------------------------

Function GInterpXYZ(Field, i, j, GHCoeffs, ZCoeffs)

  Implicit None
  ! Argument List:
  Real(Std),       Intent(In) :: Field(:, :, :)
  Integer,         Intent(In) :: i
  Integer,         Intent(In) :: j
  Type(GHCoeffs_), Intent(In) :: GHCoeffs
  Type(ZCoeffs_),  Intent(In) :: ZCoeffs
  ! Function result:
  Real(Std) :: GInterpXYZ
  ! Locals:
  Integer :: m

  GInterpXYZ = 0.0
  Do m = 1, GHCoeffs%n
    GInterpXYZ = GInterpXYZ +                                               &
                 Field(GHCoeffs%iX(i, m), GHCoeffs%iY(j, m), ZCoeffs%iZ2) * &
                 GHCoeffs%X(i, m) * GHCoeffs%Y(j, m) * ZCoeffs%Z1
  End Do
  Do m = 1, GHCoeffs%n
    GInterpXYZ = GInterpXYZ +                                               &
                 Field(GHCoeffs%iX(i, m), GHCoeffs%iY(j, m), ZCoeffs%iZ1) * &
                 GHCoeffs%X(i, m) * GHCoeffs%Y(j, m) * ZCoeffs%Z2
  End Do

End Function GInterpXYZ

!-------------------------------------------------------------------------------------------------------------

Subroutine GetHCoeffs(X, Y, HGrid, HCoeffs)
! Calculates the interpolation coefficients needed to evaluate a quantity defined on
! the grid HGrid at the point (X, Y).

  Implicit None
  ! Arguments List:
  Real(Std),      Intent(In)  :: X       !} X and Y coords of location of interest.
  Real(Std),      Intent(In)  :: Y       !}
  Type(HGrid_),   Intent(In)  :: HGrid   ! Horizontal grid of interest.
  Type(HCoeffs_), Intent(Out) :: HCoeffs ! Interpolation coefficients.
  ! Locals:
  Real(Std) :: XL          !} X and Y expressed in grid index coords.
  Real(Std) :: YL          !}
  Real(Std) :: XGridCentre
  Real(Std) :: XCycle

  ! Convert X and Y into grid index coords. Note grid indices are (1, 1) not (0, 0) at
  ! the grid point with coords (HGrid%X0, HGrid%Y0). Note also there is a need to
  ! adjust cyclic coords.

  XL = (X - HGrid%X0)/HGrid%dX + 1.0
  YL = (Y - HGrid%Y0)/HGrid%dY + 1.0

  If (HGrid%XCycle /= 0.0) Then
    XGridCentre = Real(HGrid%nX - 1, Std) * 0.5 + 1.0
    XCycle      = Abs(HGrid%XCycle / HGrid%dX)
    Do While (XL < XGridCentre - XCycle * 0.5)
      XL = XL + XCycle
    End Do
    Do While (XL > XGridCentre + XCycle * 0.5)
      XL = XL - XCycle
    End Do
  End If

  ! $$ Should check grids don't overlap in cyclic coord directions when grid set up $$
  ! $$ Also that wrapped grids do join up reasonably accurately.

  ! Calculate X coefficients.

  HCoeffs%iX1 = Floor(XL)
  HCoeffs%iX2 = HCoeffs%iX1 + 1

  HCoeffs%X1 = XL - Float(HCoeffs%iX1)
  HCoeffs%X2 = 1.0 - HCoeffs%X1

  HCoeffs%XOutside = 0

  ! Correct X coefficients for points outside grid, taking account of whether the grid
  ! is wrapped.

  If (HGrid%Wrap) Then
    If (HCoeffs%iX1 < 1) Then
      HCoeffs%iX1 = HGrid%nX
      HCoeffs%iX2 = 1
      HCoeffs%X1  = (XL + XCycle - HGrid%nX) / (XCycle - Real(HGrid%nX - 1, Std))
      HCoeffs%X2  = 1.0 - HCoeffs%X1
    Else If (HCoeffs%iX2 > HGrid%nX) Then
      HCoeffs%iX1 = HGrid%nX
      HCoeffs%iX2 = 1
      HCoeffs%X1  = (XL - HGrid%nX) / (XCycle - Real(HGrid%nX - 1, Std))
      HCoeffs%X2  = 1.0 - HCoeffs%X1
    End If
  Else
    If (HCoeffs%iX1 < 1) Then
      HCoeffs%iX1      = 1
      HCoeffs%iX2      = 2
      HCoeffs%X1       = 0.0
      HCoeffs%X2       = 1.0
      HCoeffs%XOutside = -1
    Else If (HCoeffs%iX2 > HGrid%nX) Then
      HCoeffs%iX1      = HGrid%nX - 1
      HCoeffs%iX2      = HGrid%nX
      HCoeffs%X1       = 1.0
      HCoeffs%X2       = 0.0
      HCoeffs%XOutside = 1
    End If
  End If

  ! Calculate Y coefficients.

  HCoeffs%iY1 = Floor(YL)
  HCoeffs%iY2 = HCoeffs%iY1 + 1

  HCoeffs%Y1 = YL - Float(HCoeffs%iY1)
  HCoeffs%Y2 = 1.0 - HCoeffs%Y1

  HCoeffs%YOutside = 0

  ! Correct Y coefficients for points outside grid.

  If (HCoeffs%iY1 < 1) Then
    HCoeffs%iY1      = 1
    HCoeffs%iY2      = 2
    HCoeffs%Y1       = 0.0
    HCoeffs%Y2       = 1.0
    HCoeffs%YOutside = -1
  End If
  If (HCoeffs%iY2 > HGrid%nY) Then
    HCoeffs%iY1      = HGrid%nY - 1
    HCoeffs%iY2      = HGrid%nY
    HCoeffs%Y1       = 1.0
    HCoeffs%Y2       = 0.0
    HCoeffs%YOutside = 1
  End If

  ! Indices of nearest grid point.

  If (HCoeffs%X1 > 0.5) Then
    HCoeffs%iXNear = HCoeffs%iX2
  Else
    HCoeffs%iXNear = HCoeffs%iX1
  End If

  If (HCoeffs%Y1 > 0.5) Then
    HCoeffs%iYNear = HCoeffs%iY2
  Else
    HCoeffs%iYNear = HCoeffs%iY1
  End If

  HCoeffs%Valid = .true.

End Subroutine GetHCoeffs

!-------------------------------------------------------------------------------------------------------------

Subroutine GetZCoeffs(Z, ZGrid, ZCoeffs, Z0, H, dZdZ, ZS)
! Calculates the interpolation coefficients needed to evaluate a quantity defined on
! the grid ZGrid at the point Z.

  Implicit None
  ! Arguments List:
  Real(Std),      Intent(In)           :: Z       ! Z coord of location of interest.
  Type(ZGrid_),   Intent(In)           :: ZGrid   ! Vertical grid of interest.
  Type(ZCoeffs_), Intent(Out)          :: ZCoeffs ! Interpolation coefficients.
  Real(Std),      Intent(In), Optional :: Z0      ! Roughness length (m).
  Real(Std),      Intent(In), Optional :: H       ! Boundary layer depth (m).
  Real(Std),      Intent(In), Optional :: dZdZ    ! Rate of change of height agl (m) with the Z coord.
  Real(Std),      Intent(In), Optional :: ZS      ! Surface value of the Z coord.
  ! Locals:
  Integer   :: kMin    !} Indices of levels bracketing location.
  Integer   :: kMax    !}
  Integer   :: k       ! Index of a level between kMin and kMax.
  Real(Std) :: ZetaMin !} Natural logarithm of height agl (m) + roughness length at levels bracketing location
  Real(Std) :: ZetaMax !} and at the location.
  Real(Std) :: Zeta    !}

  ! Check the grid is monotonic.
# ifdef ExtraChecks
    If (ZGrid%IncreasingValue == 0) Then
      Call Message(                                            &
             'FATAL ERROR in CoordinateSystem.GetZCoeffs: ' // &
             'grid '                                        // &
             Trim(ZGrid%Name)                               // &
             ' is non-monotonic',                              &
             3                                                 &
           )
    End If ! $$ check optional variables all present or none.
# endif

  ! Find levels surrounding point.
  kMax = ZGrid%nZ + 1
  kMin = 0
  Do
    k = kMax - kMin
    If (k <= 1) Exit
    k = kMin + k/2
    If ((Z - ZGrid%Z(k)) * ZGrid%IncreasingValue > 0.0) Then
      kMin = k
    Else
      kMax = k
    End If
  End Do

  ! Calculate coefficients (ensuring that iZ1 /= iZ2 in case the coefficients are to
  ! be modified to calculate derivatives).
  ZCoeffs%iZ1 = kMax
  ZCoeffs%iZ2 = kMin
  If (ZCoeffs%iZ2 == ZGrid%nZ) Then
    ZCoeffs%Z1  = 0.0
    ZCoeffs%Z2  = 1.0
    ZCoeffs%iZ1 = ZGrid%nZ
    ZCoeffs%iZ2 = ZGrid%nZ - 1
  Else If (ZCoeffs%iZ1 == 1) Then
    ZCoeffs%Z1  = 1.0
    ZCoeffs%Z2  = 0.0
    ZCoeffs%iZ1 = 2
    ZCoeffs%iZ2 = 1
  ! Log coefficients for horizontal wind speed and direction
  Else If (Present(Z0)) Then
    If (H >= (ZGrid%Z(kMax) - ZS) * dZdZ) Then ! $$ could kMin be higher if ZGrid%IncreasingValue < 0?
      ZetaMax    = Log(dZdZ * (ZGrid%Z(ZCoeffs%iZ1) - ZS) + Z0)
      ZetaMin    = Log(dZdZ * (ZGrid%Z(ZCoeffs%iZ2) - ZS) + Z0)
      Zeta       = Log(dZdZ * (Z                    - ZS) + Z0)
      ZCoeffs%Z1 = (Zeta - ZetaMax) / (ZetaMin - ZetaMax)
      ZCoeffs%Z2 = 1.0 - ZCoeffs%Z1
    Else
      ZCoeffs%Z1 = (Z                    - ZGrid%Z(ZCoeffs%iZ1)) / &
                   (ZGrid%Z(ZCoeffs%iZ2) - ZGrid%Z(ZCoeffs%iZ1))
      ZCoeffs%Z2 = 1.0 - ZCoeffs%Z1
    End If
  Else
    ZCoeffs%Z1 = (Z                    - ZGrid%Z(ZCoeffs%iZ1)) / &
                 (ZGrid%Z(ZCoeffs%iZ2) - ZGrid%Z(ZCoeffs%iZ1))
    ZCoeffs%Z2 = 1.0 - ZCoeffs%Z1
  End If

  If ((Z - ZGrid%Z(ZGrid%nZ)) * ZGrid%IncreasingValue > 0.0) Then
    ZCoeffs%ZOutside = 1
  Else If ((Z - ZGrid%Z(1)) * ZGrid%IncreasingValue < 0.0) Then
    ZCoeffs%ZOutside = -1
  Else
    ZCoeffs%ZOutside = 0
  End If

  ! Index of nearest grid point.
  If (ZCoeffs%Z1 > 0.5) Then
    ZCoeffs%iZNear = ZCoeffs%iZ2
  Else
    ZCoeffs%iZNear = ZCoeffs%iZ1
  End If

  ! Validity.
  ZCoeffs%Valid = .true.

End Subroutine GetZCoeffs

!-------------------------------------------------------------------------------------------------------------

Subroutine GetTCoeffs(Time, OldTime, DtSec, iOld, iNew, TCoeffs)
! Calculates the interpolation coefficients needed to evaluate a quantity defined at
! two times ('old' and 'new') at the time Time.

  Implicit None
  ! Arguments List:
  Type(ShortTime_), Intent(In)  :: Time    ! Time of interest.
  Type(ShortTime_), Intent(In)  :: OldTime ! Old time.
  Integer,          Intent(In)  :: DtSec   ! Time between old and new times in seconds.
  Integer,          Intent(In)  :: iOld    !} Indices to be used to refer to the old and new times.
  Integer,          Intent(In)  :: iNew    !}
  Type(TCoeffs_),   Intent(Out) :: TCoeffs ! Interpolation coefficients.

  ! Calculate coefficients, ensuring that T1 and T2 lie in [0, 1] even if Time lies
  ! outside [OldTime, OldTime + DtSec].

  TCoeffs%iT1 = iOld
  TCoeffs%iT2 = iNew

  TCoeffs%T1 = ShortTime2RealTime(Time - OldTime)/Float(DtSec)
  TCoeffs%T1 = Min(Max(TCoeffs%T1, 0.0), 1.0)
  TCoeffs%T2 = 1.0 - TCoeffs%T1

  ! Index of nearest time.

  If (TCoeffs%T1 > 0.5) Then
    TCoeffs%iTNear = iNew
  Else
    TCoeffs%iTNear = iOld
  End If

  TCoeffs%Valid = .true.

End Subroutine GetTCoeffs

!-------------------------------------------------------------------------------------------------------------

Subroutine InterpXYZT(H, Z, T, Field, Value)
! Interpolates a gridded three-dimensional field with two time levels to a particular
! location and time using multi-linear interpolation.

  Implicit None
  ! Argument list:
  Type(HCoeffs_), Intent(In)  :: H                 ! Horizontal interpolation coefficients.
  Type(ZCoeffs_), Intent(In)  :: Z                 ! Vertical interpolation coefficients.
  Type(TCoeffs_), Intent(In)  :: T                 ! Time interpolation coefficients.
  Real(Std),      Intent(In)  :: Field(:, :, :, :) ! Field.
  Real(Std),      Intent(Out) :: Value             ! Interpolated Value.

    Value = T%T2*(Z%Z1*(H%Y1*(H%X2*Field(H%iX1, H%iY2, Z%iZ2, T%iT1)    + &
                              H%X1*Field(H%iX2, H%iY2, Z%iZ2, T%iT1))   + &
                        H%Y2*(H%X2*Field(H%iX1, H%iY1, Z%iZ2, T%iT1)    + &
                              H%X1*Field(H%iX2, H%iY1, Z%iZ2, T%iT1)))  + &
                  Z%Z2*(H%Y1*(H%X2*Field(H%iX1, H%iY2, Z%iZ1, T%iT1)    + &
                              H%X1*Field(H%iX2, H%iY2, Z%iZ1, T%iT1))   + &
                        H%Y2*(H%X2*Field(H%iX1, H%iY1, Z%iZ1, T%iT1)    + &
                              H%X1*Field(H%iX2, H%iY1, Z%iZ1, T%iT1)))) + &
            T%T1*(Z%Z1*(H%Y1*(H%X2*Field(H%iX1, H%iY2, Z%iZ2, T%iT2)    + &
                              H%X1*Field(H%iX2, H%iY2, Z%iZ2, T%iT2))   + &
                        H%Y2*(H%X2*Field(H%iX1, H%iY1, Z%iZ2, T%iT2)    + &
                              H%X1*Field(H%iX2, H%iY1, Z%iZ2, T%iT2)))  + &
                  Z%Z2*(H%Y1*(H%X2*Field(H%iX1, H%iY2, Z%iZ1, T%iT2)    + &
                              H%X1*Field(H%iX2, H%iY2, Z%iZ1, T%iT2))   + &
                        H%Y2*(H%X2*Field(H%iX1, H%iY1, Z%iZ1, T%iT2)    + &
                              H%X1*Field(H%iX2, H%iY1, Z%iZ1, T%iT2))))

End Subroutine InterpXYZT

!-------------------------------------------------------------------------------------------------------------

Subroutine InterpXYT(H, T, Field, Value)
! Interpolates a gridded two-dimensional field with two time levels to a particular
! location and time using multi-linear interpolation.

  Implicit None
  ! Argument List:
  Type(HCoeffs_), Intent(In)  :: H              ! Horizontal interpolation coefficients.
  Type(TCoeffs_), Intent(In)  :: T              ! Time interpolation coefficients.
  Real(Std),      Intent(In)  :: Field(:, :, :) ! Field.
  Real(Std),      Intent(Out) :: Value          ! Interpolated Value.

  Value = T%T2*(H%Y1*(H%X2*Field(H%iX1, H%iY2, T%iT1)   + &
                      H%X1*Field(H%iX2, H%iY2, T%iT1))  + &
                H%Y2*(H%X2*Field(H%iX1, H%iY1, T%iT1)   + &
                      H%X1*Field(H%iX2, H%iY1, T%iT1))) + &
          T%T1*(H%Y1*(H%X2*Field(H%iX1, H%iY2, T%iT2)   + &
                      H%X1*Field(H%iX2, H%iY2, T%iT2))  + &
                H%Y2*(H%X2*Field(H%iX1, H%iY1, T%iT2)   + &
                      H%X1*Field(H%iX2, H%iY1, T%iT2)))

End Subroutine InterpXYT

!-------------------------------------------------------------------------------------------------------------

Subroutine InterpCXYT(H, T, Field, Value)
! Interpolates a gridded two-dimensional field with two time levels to a particular
! location and time using multi-linear interpolation. The field may not be defined 
! at all values of the array (e.g. cloud base or cloud top) and hence the multi-linear 
! interpolation may not use all 8 array points.

  Implicit None
  ! Argument List:
  Type(HCoeffs_), Intent(In)  :: H              ! Horizontal interpolation coefficients.
  Type(TCoeffs_), Intent(In)  :: T              ! Time interpolation coefficients.
  Real(Std),      Intent(In)  :: Field(:, :, :) ! Field.
  Real(Std),      Intent(Out) :: Value          ! Interpolated Value.

  ! Locals:
  Real(Std)                   :: C(2,2,2)       ! Array denoting use field value
  Real(Std)                   :: InterpDenom    ! Denominator of interpolation formula
  
  ! Initialise C
  C(:, :, :) = 1.0

  ! Set C to be zero if field not defined  
  If (Field(H%iX1, H%iY2, T%iT1) < 0.0) C(1,2,1) = 0.0
  If (Field(H%iX2, H%iY2, T%iT1) < 0.0) C(2,2,1) = 0.0
  If (Field(H%iX1, H%iY1, T%iT1) < 0.0) C(1,1,1) = 0.0
  If (Field(H%iX2, H%iY1, T%iT1) < 0.0) C(2,1,1) = 0.0
  If (Field(H%iX1, H%iY2, T%iT2) < 0.0) C(1,2,2) = 0.0
  If (Field(H%iX2, H%iY2, T%iT2) < 0.0) C(2,2,2) = 0.0
  If (Field(H%iX1, H%iY1, T%iT2) < 0.0) C(1,1,2) = 0.0
  If (Field(H%iX2, H%iY1, T%iT2) < 0.0) C(2,1,2) = 0.0
  
  ! Calculate denominator of interpolation formula
  InterpDenom = T%T2*(H%Y1*(H%X2*C(1,2,1) + H%X1*C(2,2,1))    + &
                      H%Y2*(H%X2*C(1,1,1) + H%X1*C(2,1,1)))   + &
                T%T1*(H%Y1*(H%X2*C(1,2,2) + H%X1*C(2,2,2))    + &
                      H%Y2*(H%X2*C(1,1,2) + H%X1*C(2,1,2)))
                      

  If (InterpDenom == 0.0) Then
  
    ! If interpolation denominator is zero then field missing 
    Value = -1.0
    
  Else  
    
    Value = (T%T2*(H%Y1*(H%X2*C(1,2,1)*Field(H%iX1, H%iY2, T%iT1)   + &
                         H%X1*C(2,2,1)*Field(H%iX2, H%iY2, T%iT1))  + &
                   H%Y2*(H%X2*C(1,1,1)*Field(H%iX1, H%iY1, T%iT1)   + &
                         H%X1*C(2,1,1)*Field(H%iX2, H%iY1, T%iT1))) + &
             T%T1*(H%Y1*(H%X2*C(1,2,2)*Field(H%iX1, H%iY2, T%iT2)   + &
                         H%X1*C(2,2,2)*Field(H%iX2, H%iY2, T%iT2))  + &
                   H%Y2*(H%X2*C(1,1,2)*Field(H%iX1, H%iY1, T%iT2)   + &
                         H%X1*C(2,1,2)*Field(H%iX2, H%iY1, T%iT2))))/ &
             InterpDenom

  End If
                               

End Subroutine InterpCXYT
!-------------------------------------------------------------------------------------------------------------

Function InterpXY(H, Field) Result(Value)
! Interpolates a gridded two-dimensional field to a particular
! location using multi-linear interpolation.

  Implicit None
  ! Argument List:
  Type(HCoeffs_), Intent(In)  :: H           ! Horizontal interpolation coefficients.
  Real(Std),      Intent(In)  :: Field(:, :) ! Field.
  ! Function result:
  Real(Std) Value ! Interpolated Value.

  Value = H%Y1*(H%X2*Field(H%iX1, H%iY2)   + &
                H%X1*Field(H%iX2, H%iY2))  + &
          H%Y2*(H%X2*Field(H%iX1, H%iY1)   + &
                H%X1*Field(H%iX2, H%iY1))

End Function InterpXY

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcArea(HCoord, HGrid, iX, iY, Area)
! Calculates area of a grid box for latitude-longitude and cartesian horizontal coord systems.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In)  :: HCoord   ! Coord system in which the horizontal grid is defined.
  Type(HGrid_),  Intent(In)  :: HGrid    ! Horizontal grid for which a grid box area is to be calculated.
  Integer,       Intent(In)  :: iX       !} Indices of grid box in the horizontal grid.
  Integer,       Intent(In)  :: iY       !}
  Real(Std),     Intent(Out) :: Area     ! Area of grid box.

  If (HCoord%CoordType == H_LatLong) Then

    Area = (EarthRadius ** 2) * HCoord%Unit(1) * HGrid%dX *                               &
             Abs(Sin((HGrid%Y(iY) + 0.5 * HGrid%dY) * HCoord%Unit(2) + HCoord%Origin(2))  &
               - Sin((HGrid%Y(iY) - 0.5 * HGrid%dY) * HCoord%Unit(2) + HCoord%Origin(2)))

  Else If (HCoord%CoordType == H_PSCartesian) Then

    Area = HGrid%dX * HGrid%dY * HCoord%Unit(1) * HCoord%Unit(2)

  Else If (HCoord%CoordType == H_TMCartesian) Then

    Area = HGrid%dX * HGrid%dY * HCoord%Unit(1) * HCoord%Unit(2)

  Else

    Call Message('Error in CalcArea', 4)

  End If

End Subroutine CalcArea

!-------------------------------------------------------------------------------------------------------------

Function InitDomains()
! Initialises a collection of domains.

  Implicit None
  ! Function result:
  Type(Domains_) :: InitDomains !
  ! Locals:
  Type(Domains_) :: Domains !

  Domains%nDomains = 0

  InitDomains = Domains

End Function InitDomains

!-------------------------------------------------------------------------------------------------------------

Function InitDomain(                                                   &
           Name,                                                       &
           HCoordName, XMin, XMax, YMin, YMax, XUnbounded, YUnbounded, &
           ZCoordName, ZMax, ZUnbounded,                               &
           StartTime, EndTime,                                         &
           MaxTravelTime,                                              &
           LocationsName, CentreName                                   &
         )
! Initialises a domain.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name
  Character(*), Intent(In) :: HCoordName
  Real(Std),    Intent(In) :: XMin
  Real(Std),    Intent(In) :: XMax
  Real(Std),    Intent(In) :: YMin
  Real(Std),    Intent(In) :: YMax
  Logical,      Intent(In) :: XUnbounded
  Logical,      Intent(In) :: YUnbounded
  Character(*), Intent(In) :: ZCoordName
  Real(Std),    Intent(In) :: ZMax
  Logical,      Intent(In) :: ZUnbounded
  Type(Time_),  Intent(In) :: StartTime
  Type(Time_),  Intent(In) :: EndTime
  Type(Time_),  Intent(In) :: MaxTravelTime
  Character(*), Intent(In) :: LocationsName ! Blank => not used. If used XMin/Max must,
  Character(*), Intent(In) :: CentreName    ! at this point, be +/- XRange / 2
  ! Function result:
  Type(Domain_) :: InitDomain !
  ! Locals:
  Type(Domain_) :: Domain !

  If (Len_Trim(Name) > MaxCharLength .or. Len_Trim(Name) == 0) Then
    Call Message('FATAL ERROR in InitDomain', 4)
  End If
  If (Len_Trim(HCoordName) > MaxCharLength) Then
    Call Message('FATAL ERROR in InitDomain', 4)
  End If
  If (Len_Trim(HCoordName) == 0 .and. (.not.XUnbounded .or. .not.YUnbounded)) Then
    Call Message('FATAL ERROR in InitDomain: H-Coord needs to be given if ' // &
                 'X is bounded or if Y is bounded', 4)
  End If
  If (Len_Trim(ZCoordName) > MaxCharLength) Then
    Call Message('FATAL ERROR in InitDomain', 4)
  End If
  If (Len_Trim(ZCoordName) == 0 .and. .not.ZUnbounded) Then
    Call Message('FATAL ERROR in InitDomain', 4)
  End If

 ! Check Time Types.

  If (IsTimeInterval(StartTime) .or. IsTimeInterval(EndTime)) Then
    Call Message(                                                             &
           'FATAL ERROR in InitDomain: the start and/or end times of the ' // &
           'domain are time intervals',                                       &
           3                                                                  &
         )
  End If

  If (StartTime > EndTime) Then ! $$ check StartTime /= infinite future and
  ! EndTime /= infinite past too - domains cannot be null in time.
    Call Message(                                                        &
           'FATAL ERROR in InitDomain: the start time of the domain ' // &
           'is after the end time',                                      &
           3                                                             &
         )
  End If

  If (.not.IsTimeInterval(MaxTravelTime)) Then
    Call Message(                                                      &
           'FATAL ERROR in InitDomain: the max travel time of the ' // &
           'domain is not a time interval',                            &
           3                                                           &
         )
  End If

  Domain%Name           = Name
  Domain%HCoordName     = HCoordName
  Domain%XMin           = XMin
  Domain%XMax           = XMax
  Domain%YMin           = YMin
  Domain%YMax           = YMax
  Domain%XUnbounded     = XUnbounded
  Domain%YUnbounded     = YUnbounded
  Domain%ZCoordName     = ZCoordName
  Domain%ZMax           = ZMax
  Domain%ZUnbounded     = ZUnbounded
  Domain%StartTime      = StartTime
  Domain%EndTime        = EndTime
  Domain%MaxTravelTime  = MaxTravelTime
  Domain%sStartTime     = Time2ShortTime(StartTime)
  Domain%sEndTime       = Time2ShortTime(EndTime)
  Domain%sMaxTravelTime = Time2ShortTime(MaxTravelTime)
  Domain%LocationsName  = LocationsName
  Domain%CentreName     = CentreName

  InitDomain = Domain

End Function InitDomain

!-------------------------------------------------------------------------------------------------------------

Subroutine AddDomain(Domain, Domains)
! Adds a horizontal Domain to a collection of Domains.

  Implicit None
  ! Argument list:
  Type(Domain_),  Intent(In)    :: Domain
  Type(Domains_), Intent(InOut) :: Domains
  ! Locals:
  Integer :: i !

  Do i = 1, Domains%nDomains
    If (Domain%Name .CIEq. Domains%Domains(i)%Name) Then
      If (Domain%HCoordName == Domains%Domains(i)%HCoordName) Then ! $$ need eq
                                                                   ! function for
                                                                   ! domains
        Return
      Else
        Call Message('FATAL ERROR in AddDomain: two different domains with the same name', 3)
      End If
    End If
  End Do
  If (Domains%nDomains >= MaxDomains) Then
    Call Message('FATAL ERROR in AddDomain: too many horizontal domains', 3)
  End If
  Domains%nDomains                  = Domains%nDomains + 1
  Domains%Domains(Domains%nDomains) = Domain

End Subroutine AddDomain

!-------------------------------------------------------------------------------------------------------------

Function FindDomainIndex(Name, Domains, Error)
! .

  Implicit None
  ! Argument list:
  Character(*),   Intent(In)            :: Name
  Type(Domains_), Intent(In)            :: Domains
  Logical,        Intent(Out), Optional :: Error
  ! Function result:
  Integer :: FindDomainIndex
  ! Locals:
  Integer :: i

  Do i = 1, Domains%nDomains
    If (Name .CIEq. Domains%Domains(i)%Name) Then
      FindDomainIndex = i
      If (Present(Error)) Error = .false.
      Return
    End If
  End Do

  If (Present(Error)) Then
    FindDomainIndex = 0
    Error = .true.
  Else
    Call Message(                     &
           'FATAL ERROR: domain "' // &
           Trim(Name)              // &
           '" not found',             &
           3                          &
         )
  End If

End Function FindDomainIndex

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpDomains_CoordsEtc(Locationses, Coords, Domains)
! Sets up Domain using info in Locationses and coord systems.

  Implicit None
  ! Argument list:
  Type(Locationses_), Intent(In)            :: Locationses
  Type(Coords_),      Intent(In),    Target :: Coords
  Type(Domains_),     Intent(InOut), Target :: Domains
  ! Locals:
  Type(Domain_), Pointer :: Domain
  Integer                :: i
  Integer                :: iLocations
  Integer                :: iCentre
  Integer                :: iHCoord1
  Integer                :: iZCoord
  Real(Std)              :: Point(2)
  Type(HCoord_), Pointer :: HCoord

  Do i = 1, Domains%nDomains

    Domain => Domains%Domains(i)

    If (.not. Domain%XUnbounded .or. .not. Domain%YUnbounded) Then

      HCoord => Coords%HCoords(FindHCoordIndex(Domain%HCoordName, Coords))

      If (Domain%CentreName /= ' ') Then

        iLocations = FindLocationsIndex(Domain%LocationsName, Locationses)
        iCentre    = FindLocationIndex(                    &
                       Domain%CentreName,                  &
                       Locationses%Locationses(iLocations) &
                     )
        iHCoord1   = FindHCoordIndex(                                            &
                       Locationses%Locationses(iLocations)%HCoordNames(iCentre), &
                       Coords                                                    &
                     )
        Point      = ConvertH(                                           &
                       Coords%HCoords(iHCoord1),                         &
                       HCoord,                                           &
                       (/                                                &
                         Locationses%Locationses(iLocations)%X(iCentre), &
                         Locationses%Locationses(iLocations)%Y(iCentre)  &
                       /)                                                &
                     )

        If (.not. Domain%XUnbounded) Then
          Domain%XMin = Domain%XMin + Point(1)
          Domain%XMax = Domain%XMax + Point(1)
        End If
        If (.not. Domain%YUnbounded) Then
          Domain%YMin = Domain%YMin + Point(2)
          Domain%YMax = Domain%YMax + Point(2)
        End If

      End If

      ! Check X Min <= X Max etc unless coord cyclic. ($$ Z too)
      If (.not. Domain%XUnbounded .and. HCoord%XCycle == 0.0) Then
        If (Domain%XMin > Domain%XMax) Then
          Call Message(                                             &
                 'FATAL ERROR in reading domain "'               // &
                 Trim(Domain%Name)                               // &
                 '" from block "Domains":'                       // &
                 'X Min must be <= X Max for non-cyclic coords',    & ! mention XRange etc $$
                 3                                                  &
               )
        End If
      End If
      If (.not. Domain%YUnbounded .and. HCoord%YCycle == 0.0) Then
        If (Domain%YMin > Domain%YMax) Then
          Call Message(                                             &
                 'FATAL ERROR in reading domain "'               // &
                 Trim(Domain%Name)                               // &
                 '" from block "Domains":'                       // &
                 'Y Min must be <= Y Max for non-cyclic coords',    &
                 3                                                  &
               )
        End If
      End If

      ! Check for domain boundaries outside +/- XCycle. ! $$ y similarly
      If (.not. Domain%XUnbounded .and. HCoord%XCycle /= 0.0) Then
        If (                                 & ! $$ is this OK? Could generate warning with
          Domain%XMax >   HCoord%XCycle .or. & ! legitimate centre/range or location/range
          Domain%XMax < - HCoord%XCycle .or. & ! input?.
          Domain%XMin >   HCoord%XCycle .or. &
          Domain%XMin < - HCoord%XCycle      &
        ) Then
          Call Message(                                                                 &
                 'WARNING in processing domain "'                                    // &
                 Trim(Domain%Name)                                                   // &
                 '": X Min and/or X Max are larger than the coordinate periodicity ' // & ! mention
                                                                                          ! XRange etc $$
                 Trim(Std2Char(HCoord%XCycle)),                                         &
                 1                                                                      &
               )
        End If
        If (Abs(Domain%XMax - Domain%XMin) >= HCoord%XCycle) Then
          Call Message(                                                                &
                 'FATAL ERROR in processing domain "'                               // &
                 Trim(Domain%Name)                                                  // &
                 '": |X Max - X Min| must be less than the coordinate periodicity ' // & ! mention
                                                                                         ! XRange etc $$
                 Trim(Std2Char(HCoord%XCycle))                                      // &
                 ' (use "X Unbounded" option for unbounded domains)',                  &
                 3                                                                     &
               )
        End If
      End If ! $$ If pass this, could check that domain size is roughly right
             ! (e.g. check for Xmin = -0.1 XMax =359.89
             ! but v small domain after reducing 359.99 inaccurately to -0.09

      ! Ensure XMin and XMax lie in +/- HCoord%XCycle / 2. ! $$ y similarly
      If (.not. Domain%XUnbounded .and. HCoord%XCycle /= 0.0) Then
        Do While (Domain%XMin < - HCoord%XCycle * 0.5)
          Domain%XMin = Domain%XMin + HCoord%XCycle
        End Do
        Do While (Domain%XMin > HCoord%XCycle * 0.5)
          Domain%XMin = Domain%XMin - HCoord%XCycle
        End Do
        Do While (Domain%XMax < - HCoord%XCycle * 0.5)
          Domain%XMax = Domain%XMax + HCoord%XCycle
        End Do
        Do While (Domain%XMax > HCoord%XCycle * 0.5)
          Domain%XMax = Domain%XMax - HCoord%XCycle
        End Do
      End If

    End If

  End Do

End Subroutine SetUpDomains_CoordsEtc

!-------------------------------------------------------------------------------------------------------------

Subroutine SetUpDomains_iCoords(Coords, Domains)
! Sets up indices in Domain for referring to coord systems.

  Implicit None
  ! Argument list:
  Type(Coords_),  Intent(In)            :: Coords
  Type(Domains_), Intent(InOut), Target :: Domains
  ! Locals:
  Type(Domain_), Pointer :: Domain
  Integer                :: i

  Do i = 1, Domains%nDomains

    Domain => Domains%Domains(i)

    If (Domain%XUnbounded .and. Domain%YUnbounded) Then
      Domain%iHCoord = 0
    Else
      Domain%iHCoord = FindHCoordIndex(Domain%HCoordName, Coords)
    End If
    If (Domain%ZUnbounded) Then
      Domain%iZCoord = 0
    Else
      Domain%iZCoord = FindZCoordIndex(Domain%ZCoordName, Coords)
    End If

  End Do

End Subroutine SetUpDomains_iCoords

!-------------------------------------------------------------------------------------------------------------

Function StartTimeOfDomain(Domain)
! Returns the lower time limit of a domain.

  Implicit None
  ! Argument list:
  Type(Domain_), Intent(In) :: Domain !
  ! Function result:
  Type(Time_) :: StartTimeOfDomain !

  StartTimeOfDomain = Domain%StartTime

End Function StartTimeOfDomain

!-------------------------------------------------------------------------------------------------------------

Function EndTimeOfDomain(Domain)
! Returns the lower time limit of a domain.

  Implicit None
  ! Argument list:
  Type(Domain_), Intent(In) :: Domain !
  ! Function result:
  Type(Time_) :: EndTimeOfDomain !

  EndTimeOfDomain = Domain%EndTime

End Function EndTimeOfDomain

!-------------------------------------------------------------------------------------------------------------

Function DomainStart(Domain)
! Returns the lower time limit of a domain.

  Implicit None
  ! Argument list:
  Type(Domain_), Intent(In) :: Domain !
  ! Function result:
  Type(ShortTime_) :: DomainStart !

  DomainStart = Domain%sStartTime

End Function DomainStart

!-------------------------------------------------------------------------------------------------------------

Function DomainEnd(Domain)
! Returns the lower time limit of a domain.

  Implicit None
  ! Argument list:
  Type(Domain_), Intent(In) :: Domain !
  ! Function result:
  Type(ShortTime_) :: DomainEnd !

  DomainEnd = Domain%sEndTime

End Function DomainEnd

!-------------------------------------------------------------------------------------------------------------

Function XInDomain(Position, Domain, Coords)
! Finds out whether a point lies in a domain, when the coords are known in the coord
! systems in which the domain is defined.

  Implicit None
  ! Argument list:
  Type(Position_), Intent(In)         :: Position ! Coords of the point in various coord systems in Coords,
                                                  ! with flags to indicate whether the values are valid. The
                                                  ! coords corresponding to the coord system in which Domain
                                                  ! is defined must be valid.
  Type(Domain_),   Intent(In)         :: Domain   ! Domain.
  Type(Coords_),   Intent(In), Target :: Coords
  ! Function result:
  Logical :: XInDomain ! Indicates if the point lies in the domain.

  XInDomain = ZInDomain(Position, Domain, Coords)

  ! Test for point in horizontal extent of domain.
  If (XInDomain) Then
    XInDomain = HInDomain(Position, Domain, Coords)
  End If

End Function XInDomain

!-------------------------------------------------------------------------------------------------------------

Function HInDomain(Position, Domain, Coords)
! Finds out whether a point lies in a domain, when the coords are known in the coord
! systems in which the domain is defined.

  Implicit None
  ! Argument list:
  Type(Position_), Intent(In)         :: Position ! Coords of the point in various coord systems in Coords,
                                                  ! with flags to indicate whether the values are valid. The
                                                  ! coords corresponding to the coord system in which Domain
                                                  ! is defined must be valid.
  Type(Domain_),   Intent(In)         :: Domain   ! Domain.
  Type(Coords_),   Intent(In), Target :: Coords
  ! Function result:
  Logical :: HInDomain ! Indicates if the point lies in the domain.
  ! Locals:
  Real(Std)              :: Point(2) !
  Real(Std)              :: XMin
  Real(Std)              :: XMax
  Type(HCoord_), Pointer :: HCoord

  HInDomain = .true.

  ! Test for point in horizontal extent of domain.
  If (.not.Domain%XUnbounded .or. .not.Domain%YUnbounded) Then

    If (.not.Position%XYValid(Domain%iHCoord)) Then
      Call Message('FATAL ERROR in XInDomain', 4)
    End If

    Point(:) = Position%XY(:, Domain%iHCoord)

    HCoord => Coords%HCoords(Domain%iHCoord)

    If (.not.Domain%YUnbounded) Then
      If (Point(2) < Domain%YMin .or. Point(2) > Domain%YMax) Then
        HInDomain = .false.
        Return
      End If
    End If

    If (.not.Domain%XUnbounded) Then

      If (HCoord%XCycle /= 0.0) Then ! $$ y similarly

        Do While (Point(1) < - HCoord%XCycle * 0.5)
          Point(1) = Point(1) + HCoord%XCycle
        End Do
        Do While (Point(1) > HCoord%XCycle * 0.5)
          Point(1) = Point(1) - HCoord%XCycle
        End Do

        If (                                                                   &
          .not.(                                                               &
            (Domain%XMin <= Point(1)    .and. Point(1)    <= Domain%XMax) .or. &
            (Point(1)    <= Domain%XMax .and. Domain%XMax <= Domain%XMin) .or. &
            (Domain%XMax <= Domain%XMin .and. Domain%XMin <= Point(1)   )      &
          )                                                                    &
        ) Then
          HInDomain = .false.
          Return
        End If

      Else

        If (Point(1) < Domain%XMin .or. Point(1) > Domain%XMax) Then
          HInDomain = .false.
          Return
        End If

      End If

    End If

  End If

End Function HInDomain

!-------------------------------------------------------------------------------------------------------------

Function ZInDomain(Position, Domain, Coords)
! Finds out whether a point lies in a domain, when the coords are known in the coord
! systems in which the domain is defined.

  Implicit None
  ! Argument list:
  Type(Position_), Intent(In)         :: Position ! Coords of the point in various coord systems in Coords,
                                                  ! with flags to indicate whether the values are valid. The
                                                  ! coords corresponding to the coord system in which Domain
                                                  ! is defined must be valid.
  Type(Domain_),   Intent(In)         :: Domain   ! Domain.
  Type(Coords_),   Intent(In), Target :: Coords
  ! Function result:
  Logical :: ZInDomain ! Indicates if the point lies in the domain.

  ZInDomain = .true.

  If (.not.Domain%ZUnbounded) Then

    If (.not.Position%ZValid(Domain%iZCoord)) Then
      Call Message('FATAL ERROR in XInDomain', 4)
    End If

    If (                                                 &
      (                                                  &
        Position%Z(Domain%iZCoord) > Domain%ZMax   .and. &
        ZLike(Coords%ZCoords(Domain%iZCoord))            &
      )                                                  &
      .or.                                               &
      (                                                  &
        Position%Z(Domain%iZCoord) < Domain%ZMax   .and. &
        .not.ZLike(Coords%ZCoords(Domain%iZCoord))       &
      )                                                  &
    ) Then
      ZInDomain = .false.
      Return
    End If

  End If

End Function ZInDomain

!-------------------------------------------------------------------------------------------------------------

Function SInDomain(TravelTime, Domain)
! Finds out whether a travel time lies in a domain.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: TravelTime ! Travel Time.
  Type(Domain_),    Intent(In) :: Domain     ! Domain.
  ! Function result:
  Logical :: SInDomain ! Indicates if the point lies in the domain.

  SInDomain = .true.

  If (TravelTime >= Domain%sMaxTravelTime) Then
    SInDomain = .false.
  End If

End Function SInDomain

!-------------------------------------------------------------------------------------------------------------

Function TInDomain(Time, Domain)
! Finds out whether a time lies in a domain.

  Implicit None
  ! Argument list:
  Type(ShortTime_), Intent(In) :: Time   ! Time.
  Type(Domain_),    Intent(In) :: Domain ! Domain.
  ! Function result:
  Logical :: TInDomain ! Indicates if the point lies in the domain.

  TInDomain = .true.

  If (Time >= Domain%sEndTime) Then
    TInDomain = .false.
  Else If (Time < Domain%sStartTime) Then
    TInDomain = .false.
  End If

End Function TInDomain

!-------------------------------------------------------------------------------------------------------------

Function DistanceToHEdge(Point, HCoord, Domain)
! Computes distance to edge of domain. $$ Not used currently.

  Implicit None
  ! Argument list:
  Real(Std),     Intent(In) :: Point(2) !
  Type(HCoord_), Intent(In) :: HCoord   !
  Type(Domain_), Intent(In) :: Domain   !
  ! Function result:
  Real(Std) :: DistanceToHEdge !
  ! Locals:
  Integer   :: i  !
  Real(Std) :: d  !
  Real(Std) :: x  !
  Real(Std) :: y  !
  Real(Std) :: x1 !
  Real(Std) :: y1 !
  Real(Std) :: x2 !
  Real(Std) :: y2 !

  x = Point(1)
  y = Point(2)
  DistanceToHEdge = 1.0E20_Std  ! $$ use Huge
  If (HCoord%CoordType == H_LatLong) Then
  ! $$
  ElseIf (HCoord%CoordType == H_PSCartesian) Then
    Do i = 1,4
      If (i == 1) Then !$$ This is not a sensible approach for square domains - but
                       !   currently not used anyway
        x1 = Domain%XMin
        y1 = Domain%YMin
        x2 = Domain%XMax
        y2 = Domain%YMin
      Else If (i == 2) Then
        x1 = Domain%XMax
        y1 = Domain%YMin
        x2 = Domain%XMax
        y2 = Domain%YMax
      Else If (i == 3) Then
        x1 = Domain%XMin
        y1 = Domain%YMax
        x2 = Domain%XMax
        y2 = Domain%YMax
      Else If (i == 4) Then
        x1 = Domain%XMin
        y1 = Domain%YMin
        x2 = Domain%XMin
        y2 = Domain%YMax
      End If
      d = ((x - x1)*(y1 - y2) - (y - y1)*(x1 - x2))/Sqrt((y1 - y2)**2 + (x1 - x2)**2)
      DistanceToHEdge = Min(DistanceToHEdge, d)
    EndDo
  End If

End Function DistanceToHEdge

!-------------------------------------------------------------------------------------------------------------

Function CalcZenithAngle(Coords, Time, Position)
! Computes an approximate solar zenith angle (in radians) at the given position/time.
! (adapted from the subroutine ZENITH in NAME V8.08).
!$$ Precision/accuracy of this calculation??
! $$ Not clear where this routine should live.
! $$ review approach for 360 day year calendar

  Implicit None
  ! Argument list:
  Type(Coords_),   Intent(In)    :: Coords    ! Collection of coord systems.
  Type(Time_),     Intent(In)    :: Time      ! Current time.
  Type(Position_), Intent(InOut) :: Position  ! Position of interest.
  ! Function result:
  Real(Std) :: CalcZenithAngle ! Solar zenith angle (in radians).
  ! Locals:
  Type(WildTime_) :: WildTime ! Wild card time used to construct TimeRef.
  Type(Time_) :: TimeRef  ! Reference time (00:00 UTC, 21/6/YY where YY = Time%Year).
  Real(Std)   :: TimeDiff ! Time difference between the current and reference times.
  Integer     :: iLatLong ! Index (in Coords) of the standard lat-long coord system.
  Integer     :: iZCoord  ! Index (in Coords) of a valid vertical coord system.
  Integer     :: i        ! Loop index.
  Real(Std)   :: X(3)     ! Lat-long-height coords of position (lat-long in degrees).
  Real(Std)   :: HrAngle  ! Local hour angle.
  Real(Std)   :: Decl     ! Solar declination.
  Real(Std)   :: CosZen   ! Cosine of the zenith angle.

  ! Calculate time difference (in seconds) between the current time
  ! and the reference time (at midnight on the summer solstice).
  WildTime = Char2WildTime('21/6/* 00:00') ! $$ doesn't work in relative time frame
  TimeRef = SubstituteWildCards(Time, WildTime)

  TimeDiff = ShortTime2RealTime(Time2ShortTime(Time - TimeRef))

  ! Find index of the standard lat-long coord system.
  iLatLong = FindHCoordIndex('Lat-Long', Coords)

  ! Calculate lat-long of the position (if not already available).
  If (.not.Position%XYValid(iLatLong)) Then
    Call ConvertToH(Coords, iLatLong, Position)
  End If

  ! Determine a vertical coord system in which the position is defined. Any valid
  ! vertical coord system is sufficient here since height is not used any further.
  Do i = 1, Coords%nZCoords
    If (Position%ZValid(i)) Then
      iZCoord = i
      Exit
    Else If (i == Coords%nZCoords) Then
      Call Message('FATAL ERROR in CalcZenithAngle: undefined vertical coordinate', 4)
    End If
  End Do

  ! Retrieve lat-long of the position.
  X = Position2X(Coords, Position, iLatLong, iZCoord)

  ! Compute solar zenith angle based on lat-long location and current time.
  HrAngle = (1.0 + TimeDiff/4.32E4 + X(1)/180.0) * Pi
  Decl    = 0.4142 * Cos(2.0*Pi*TimeDiff/3.1536E+7)
  CosZen  = Cos(HrAngle)*Cos(X(2)*Pi/180.0)*Cos(Decl) + Sin(X(2)*Pi/180.0)*Sin(Decl)
  CalcZenithAngle = ACos(CosZen)

End Function CalcZenithAngle

!-------------------------------------------------------------------------------------------------------------

End Module GridAndDomainModule
