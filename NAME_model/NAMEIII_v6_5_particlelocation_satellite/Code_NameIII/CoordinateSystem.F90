! Module:  Coordinate System Module

Module CoordinateSystemModule

! This module provides code for handling coord systems.

!-------------------------------------------------------------------------------------------------------------

Use GlobalParametersModule
Use ErrorAndMessageModule
Use MathsModule
Use PhysicsModule
Use StringModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public  :: H_LatLong
Public  :: H_PSCartesian
Public  :: H_PSPolar
Public  :: H_TMCartesian
Public  :: H_TMPolar
Public  :: Z_AboveGround
Public  :: Z_AboveSea
Public  :: Z_P
Public  :: Z_PAsZ
Public  :: Z_PAsEta
Public  :: Z_ZAsEta
Public  :: EtaDefn_
Public  :: EtaDefns_
Public  :: HCoord_
Public  :: ZCoord_
Public  :: Coords_
Public  :: Position_
Public  :: Operator(==)
Public  :: Operator(.equiv.)
Public  :: ZLike
Public  :: InitEtaDefns
Public  :: InitEtaDefn
Public  :: AddEtaDefn
Public  :: FindEtaDefnIndex
Public  :: InitHCoord
Public  :: HCoord_LatLong
Public  :: HCoord_LatLongRadians
Public  :: HCoord_EMEP50kmGrid
Public  :: HCoord_EMEP150kmGrid
Public  :: HCoord_UKNationalGridM
Public  :: HCoord_UKNationalGrid100M
Public  :: ConvertH
Public  :: MetricCoeffs
Public  :: CalcDistanceBetweenTwoPoints
Public  :: CalcHAngle
Public  :: CalcdZdXZBased
Public  :: CalcdZdXPBased
Public  :: CalcdZdZZBased
Public  :: CalcdZdZPBased
Public  :: InitZCoord
Public  :: ZCoord_m_agl
Public  :: ZCoord_m_asl
Public  :: ZCoord_Pa
Public  :: ZCoord_FL
Public  :: ZBasedToZBased
Public  :: PBasedToPBased
Public  :: InitCoords
Public  :: AddHCoord
Public  :: FindHCoordIndex
Public  :: AddZCoord
Public  :: FindZCoordIndex
Public  :: X2Position
Public  :: Position2X
Public  :: Position2XUnknownCoord
Public  :: ConvertToH
Public  :: AboveSeaToPAsZICAO
Public  :: PAsZToAboveSeaICAO
Public  :: RationaliseHPosition

!-------------------------------------------------------------------------------------------------------------

! Codes for types of horizontal coord system.
Integer, Parameter :: H_LatLong     = 1 ! Lat-long coord system.
Integer, Parameter :: H_PSCartesian = 2 ! Cartesian coord system in PS projection.
Integer, Parameter :: H_PSPolar     = 3 ! Polar coord system in PS projection.
Integer, Parameter :: H_TMCartesian = 4 ! Cartesian coord system in TM projection.
Integer, Parameter :: H_TMPolar     = 5 ! Polar coord system in TM projection.

! Codes for types of vertical coord system.
Integer, Parameter :: Z_AboveGround = 1 ! Height above ground.
Integer, Parameter :: Z_AboveSea    = 2 ! Height above sea.
Integer, Parameter :: Z_P           = 3 ! Pressure.
Integer, Parameter :: Z_PAsZ        = 4 ! Pressure, converted to height above sea
                                        ! level using the ICAO standard atmosphere.
Integer, Parameter :: Z_PAsEta      = 5 ! Pressure-based eta.
Integer, Parameter :: Z_ZAsEta      = 6 ! Height-based eta.

!-------------------------------------------------------------------------------------------------------------

Type :: EtaDefn_ ! A set of arrays, used for defining vertical pressure-based eta
                 ! coord systems (at least only pressure-based at present).
  Character(MaxCharLength)         :: Name
  Character(1)                     :: EtaDefnType ! Exactly two of Eta, A and B must be given. 'E' = Eta
                                                  ! missing, 'A' = A missing, and 'B' = B missing.
  Integer                          :: nEtaLevels ! $$ Needed?
  Real(Std),               Pointer :: Eta(:)
  Real(Std),               Pointer :: A(:)
  Real(Std),               Pointer :: B(:)
End Type EtaDefn_

!-------------------------------------------------------------------------------------------------------------

Type :: EtaDefns_ ! A collection of eta definitions.
  Integer        :: nEtaDefns
  Type(EtaDefn_) :: EtaDefns(MaxEtaDefns)
End Type EtaDefns_

!-------------------------------------------------------------------------------------------------------------

Type :: HCoord_ ! A horizontal coord system.
  Character(MaxCharLength) :: Name
  Integer                  :: CoordType
  Real(Std)                :: Pole(2)
  Real(Std)                :: Angle
  Real(Std)                :: ScaleFactor
  Real(Std)                :: Origin(2)
  Real(Std)                :: ThetaOrigin
  Real(Std)                :: Unit(2)
  Real(Std)                :: XCycle
  Real(Std)                :: YCycle
  Logical                  :: YEdge
  Real(Std)                :: YEdgeMax
  Real(Std)                :: YEdgeMin
  ! Name        :: Name of coord system.
  ! CoordType   :: 1 = coord system based on a latitude-longitude system with
  !                    arbitrary orientation, i.e. arbitrary location for the
  !                    latitude-longitude system's north pole and for the 'third'
  !                    Euler angle representing a rotation about the system's north
  !                    pole ('coord system based on a latitude-longitude system'
  !                    rather than 'latitude-longitude system' because of the
  !                    possibility of a non-zero origin offset and non-standard choice
  !                    of units),
  !                2 = Cartesian coord system using a polar stereographic (PS)
  !                    projection,
  !                3 = polar coord system using a polar stereographic (PS) projection,
  !                4 = Cartesian coord system using a transverse Mercator (TM)
  !                    projection,
  !                5 = polar coord system using a transverse Mercator (TM) projection.
  ! Pole        :: For coord systems of type 1 (lat-long):
  !                    position of the coord system's north pole in the standard
  !                    latitude-longitude coord system (in radians), but with latitude
  !                    replaced by angle from the true north pole.
  !                For coord systems of type 2 and 3 (PS):
  !                    position of the tangent point of the PS projection in the
  !                    standard latitude-longitude coord system (in radians), but with
  !                    latitude replaced by angle from the true north pole.
  !                For coord systems of type 4 and 5 (TM):
  !                    position of the true origin of the TM projection in the
  !                    standard latitude-longitude coord system (in radians).
  !                The reason for using angle from true north pole in 1-3 rather than
  !                latitude is so that a standard latitude-longitude coord system is
  !                represented exactly.
  ! Angle       :: For coord systems of type 1 (lat-long):
  !                    rotation of the coord system relative to the system with the
  !                    same north pole location and with the zero longitude line
  !                    passing through the true south pole (in radians).
  !                For coord systems of type 2 (PS Cartesian):
  !                    angle between the x axis and the Pi/2 longitude line of a type
  !                    1 coord system with the same Pole but with a zero value of
  !                    Angle (in radians).
  !                For coord systems of type 3 (PS polar):
  !                    angle between the x axis of the underlying Cartesian coord
  !                    system and the Pi/2 longitude line of a type 1 coord system
  !                    with the same Pole but with a zero value of Angle (in radians).
  !                For coord systems of type 4 (TM Cartesian):
  !                    Angle is not used but is always set to 0. The coord system is
  !                    always east-north for type 4.
  !                For coord systems of type 5 (TM polar):
  !                    Angle is not used but is always set to 0. The underlying
  !                    Cartesian coord system is always east-north for type 5.
  !                In each case a positive value means that, standing at the north
  !                pole or tangent point and looking down, the system is rotated
  !                anticlockwise relative to its orientation for a zero value.
  ! ScaleFactor :: For coord systems of type 1, 2 and 3 (lat-long and PS):
  !                    ScaleFactor is not used but is always set to 1.
  !                For coord systems of type 4 and 5 (TM):
  !                    Scale adjustment used by the TM projection. Must be > 0. A
  !                    scale factor of 1.0 preserves distance along the central
  !                    meridian of the projection (with map points away from the
  !                    meridian being stretched). A scale factor of less than 1 is
  !                    normally applied to shrink the map projection to improve the
  !                    overall representation of distance. ScaleFactor plays a
  !                    different role from Unit, because it is applied to the basic
  !                    projection before Origin (see below) is considered, and hence
  !                    affects the way Origin should be interpreted. Unit affects only
  !                    the location relative to Origin.
  ! Origin      :: Offset of the origin of the coordinate system from:
  !                    long = lat = 0 for type 1 coord systems (lat-long);
  !                    the tangent point for type 2 and 3 coord systems (PS);
  !                    the true origin of the TM projection for type 4 and 5 coord
  !                    systems (TM).
  !                For polar coordinate systems (type 3 and 5 coord systems) the
  !                offset is a Cartesian offset in the underlying Cartesian coord
  !                system. Values are in radians for type 1 coord systems and in
  !                metres for the other types.
  ! ThetaOrigin :: For coord systems of type 1, 2 and 4 (lat-long and Cartesian):
  !                    ThetaOrigin is not used, but is always set to 0.
  !                For coord systems of type 3 and 5 (polar):
  !                    Angle between the line theta = 0 and the x axis of the
  !                    underlying Cartesian coord system. A positive value means that,
  !                    standing at r = 0 and looking down, the line theta = 0 is
  !                    rotated anticlockwise from the x axis.
  ! Unit        :: Units used for the two coordinates relative to:
  !                    radians for type 1 coord systems (lat-long);
  !                    metres for type 2 and 4 coord systems (Cartesian);
  !                    metres and radians for type 3 and 5 coord systems (polar).
  !                Values greater than 1 indicate that the units are bigger and the
  !                values smaller.
  ! XCycle       :] Coordinate change in cycling the coordinate for cyclic coords. 0
  ! YCycle       :] otherwise.
  ! YEdge        :: Indicates the coordinate system has an edge in the y-direction.
  ! YEdgeMax     ::} Location of coordinate system edge.
  ! YEdgeMin     ::}

  ! The above refers to the direction of the x axis. This means the direction in which
  ! x would increase if we had Unit(1) > 0.

  ! For polar coordinates, the above refers to the 'underlying Cartesian coord
  ! system'. This is the coord system using the same projection and the same values of
  ! Pole, Angle and ScaleFactor, but with no origin offset and values in metres (i.e.
  ! Origin = 0 and Unit = 1).

  ! For lat-long systems, coord 1 and coord 2 are often referred to as lambda and phi
  ! with lambda related to longitude and phi related to latitude. Lambda and phi
  ! increase to the east and north (when Unit > 0).

  ! For Cartesian coord systems, coord 1 and coord 2 are often referred to as x and y.
  ! The positive y axis is rotated 90 degrees anticlockwise from the positive x axis
  ! (when Unit > 0).

  ! For polar coord systems, coord 1 and coord 2 are often referred to as r and theta
  ! with r related to distance from the origin and theta related to the angle about
  ! the origin. Theta increases anticlockwise (when Unit(2) > 0).

  ! In two element arrays such as Pole, Origin, Unit, or any given pairs of coord
  ! values, the first element refers to coord 1.

End Type HCoord_

!-------------------------------------------------------------------------------------------------------------

Type :: ZCoord_ ! A vertical coord system.
  Character(MaxCharLength)         :: Name
  Integer                          :: CoordType
  Real(Std)                        :: Unit
  Integer                          :: nEtaLevels
  Real(Std),               Pointer :: Eta(:)
  Real(Std),               Pointer :: A(:)
  Real(Std),               Pointer :: B(:)
  Real(Std)                        :: PStarMin
  Real(Std)                        :: ModelTopHeight
  Real(Std)                        :: InterfaceHeight
  Real(Std)                        :: InterfaceEta
  Real(Std)                        :: MaxOrogHeight
  Logical                          :: IncreasingUpwards
  ! Name       :: Name of coord system.
  ! CoordType  :: 1 = height above ground (metres),
  !               2 = height above sea (metres),
  !               3 = pressure (Pa),
  !               4 = pressure, converted to height above sea level using the ICAO
  !                   standard atmosphere (metres),
  !               5 = pressure based eta coord (dimensionless),
  !               6 = height based eta coord (dimensionless).

  ! The height based eta coord is a simplified special case of a composite
  ! linear/quadratic transformation, as used by the New Dynamics formulation of
  ! the Unified Model (UM 5.2 onwards). This eta coord
  ! system has two parameters: height of the model top, and height of the
  ! interface between the lower (quadratic) region and upper (linear) region.
  ! The general case, which is not necessarily linear over the sea, would also require
  ! the eta value of the interface.

  ! Unit       :: Units used for the coord relative to above units. Values greater
  !               than 1 indicate that the units are bigger and the values smaller.
  !               Unit is not used for type 5 and 6 coord systems (eta coords) but
  !               its value is always set to 1.0.
  ! nEtaLevels :: Number of eta levels (defined for type 5 coord systems only).
  ! Eta        :} At a given level, eta = A/p0 + B and p = A + B*pstar,
  ! A          :} where p is pressure, p0 = 10E5 Pa and pstar = surface
  ! B          :} pressure (in Pa). In between levels, p and eta are
  ! PStarMin   :} linearly related. This is also the case above the top
  !               level, p and eta reaching zero together. The first
  !               level must be on the ground, i.e. eta = 1, A = 0, B =
  !               1, and eta and B must decrease with increasing level
  !               number. Usually B = 0 near the top levels so that eta
  !               is a function of p only. If A = 0, eta = p/pstar. p
  !               should decrease with level number, and this will be
  !               the case provided pstar exceeds PStarMin. These
  !               variables are defined for type 5 coord systems only.
  ! (Note p = eta * p0 + B * (PStar - P0) - so B decreasing is necessary for effect of
  ! pstar to decrease with height. If B increased one would have a PStarMax too).
  ! ModelTopHeight  :: Height above sea level of the model top (in metres); this
  !                    corresponds to eta = 1.
  ! InterfaceHeight :: Height above sea level of the linear/quadratic interface
  !                    (in metres); that is, the lowest constant-height eta level.
  ! InterfaceEta    :: Eta value on the linear/quadratic interface
  !                    (= InterfaceHeight / ModelTopHeight).
  ! MaxOrogHeight   :: Maximum orography height (in metres) allowable for a monotonic
  !                    height-based eta transformation (= 1/2 * InterfaceHeight).
  !                    The above four variables are defined for type 6 coord systems
  !                    only.
! $$ store p0 with coord (and call PRef).
  ! IncreasingUpwards :: Indicates the coord values increase upwards.

! $$ Note and ensure the following is true (its assumed in nwpflow.PToEta etc)
  ! Note that, above and below
  ! ZGrid, z-like eta is defined as having d eta / dz constant and equal to the grid
  ! edge value, and p-like eta is defined as proportional to pressure (i.e. defined as
  ! having d log eta / d log p constant and equal to 1).

  ! Note that coord systems are defined below ground. Here
  ! z-like eta is defined as having d eta / dz constant and equal to the ground
  ! value, and p-like eta is defined as proportional to pressure (i.e. defined as
  ! having d log eta / d log p constant and equal to 1).
  ! PAsZ coords are defined assuming the ICAO lapse rate at the ground remains constant
  ! below ground (above top of ICAO region?).

  ! No conversion between p-based and z-based systems is supported here. Any code
  ! using these routines must make its own assumptions.

End Type ZCoord_

!-------------------------------------------------------------------------------------------------------------

Type :: Coords_ ! A collection of coord systems.
  Integer       :: nHCoords            ! Number of horizontal coord systems.
  Type(HCoord_) :: HCoords(MaxHCoords) ! Collection of horizontal coord systems.
  Integer       :: nZCoords            ! Number of horizontal coord systems.
  Type(ZCoord_) :: ZCoords(MaxZCoords) ! Collection of vertical coord systems.
End Type Coords_

!-------------------------------------------------------------------------------------------------------------

Type :: Position_ ! A position in a variety of coord systems.
  Real(Std) :: XY(2, MaxHCoords)   !
  Real(Std) :: Z(MaxZCoords)       !
  Real(Std) :: Topog               !
  Real(Std) :: PS                  !
  Real(Std) :: Rho                 !
  Logical   :: XYValid(MaxHCoords) !
  Logical   :: ZValid(MaxZCoords)  !
  Logical   :: TopogValid          !
  Logical   :: PSValid             !
  Logical   :: RhoValid            !
End Type Position_

!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of coord systems, grids and eta definitions.
  Module Procedure EtaDefnEq
  Module Procedure HCoordEq
  Module Procedure ZCoordEq
End Interface

!-------------------------------------------------------------------------------------------------------------

Interface Operator(.equiv.) ! Equivalence of coord systems (i.e. equality except
                            ! possibly of their names).
  Module Procedure HCoordEquiv
  Module Procedure ZCoordEquiv
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function ZLike(ZCoord)
! Indicates that the size of ZCoord increases without bound as height increases. If
! false, then ZCoord approaches zero as height increases.

! Note that although coord systems can be complex functions of height, they are all
! defined so that at large heights or below ground they are, at least approximately,
! affine functions of height (ZLike = .true.) or proportional to pressure (ZLike =
! .false.) with constants which may vary with position.

! $$ Need Unit > 0

  Implicit None
  ! Argument list:
  Type(ZCoord_) :: ZCoord !
  ! Function result:
  Logical :: ZLike !

  Select Case (ZCoord%CoordType)
    Case (Z_AboveGround, Z_AboveSea, Z_PAsZ, Z_ZAsEta)
      ZLike = .true.
    Case (Z_P, Z_PAsEta)
      ZLike = .false.
  End Select

End Function ZLike

!-------------------------------------------------------------------------------------------------------------

Function InitEtaDefns() Result (EtaDefns)
! Initialises a collection of eta definitions.

  Implicit None
  ! Function result:
  Type(EtaDefns_) :: EtaDefns ! Initialised collection of eta definitions.

  EtaDefns%nEtaDefns = 0

End Function InitEtaDefns

!-------------------------------------------------------------------------------------------------------------

Function InitEtaDefn(Name, nEtaLevels, Eta, A, B) Result (EtaDefn)
! Initialises an eta definition.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: Name       ! Name of eta definition.
  Integer,      Intent(In)           :: nEtaLevels ! Number of eta levels.
  Real(Std),    Intent(In), Optional :: Eta(:)     ! Eta values.
  Real(Std),    Intent(In), Optional :: A(:)       ! A values.
  Real(Std),    Intent(In), Optional :: B(:)       ! B values.
  ! Function result:
  Type(EtaDefn_) :: EtaDefn ! Initialised eta defn.

  If (Name == ' ') Then
    Call Message('FATAL ERROR in InitEtaDefn: name of eta definition is blank', 3)
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message(                                                                &
           'FATAL ERROR in InitEtaDefn: name of eta definition is given as "' // &
           Trim(Name)                                                         // &
           '" and is too long',                                                  &
           3                                                                     &
         )
  End If

  ! $$ Check dimensions consistent etc

  If (.not.Present(Eta) .and. Present(A) .and. Present(B)) Then
    Allocate(EtaDefn%A(nEtaLevels))
    Allocate(EtaDefn%B(nEtaLevels))
    EtaDefn%Name        = Name
    EtaDefn%EtaDefnType = 'E'
    EtaDefn%nEtaLevels  = nEtaLevels
    EtaDefn%A(:)        = A(:)
    EtaDefn%B(:)        = B(:)
  Else If (Present(Eta) .and. .not.Present(A) .and. Present(B)) Then
    Allocate(EtaDefn%Eta(nEtaLevels))
    Allocate(EtaDefn%B  (nEtaLevels))
    EtaDefn%Name        = Name
    EtaDefn%EtaDefnType = 'A'
    EtaDefn%nEtaLevels  = nEtaLevels
    EtaDefn%Eta(:)      = Eta(:)
    EtaDefn%B(:)        = B(:)
  Else If (Present(Eta) .and. Present(A) .and. .not.Present(B)) Then
    Allocate(EtaDefn%Eta(nEtaLevels))
    Allocate(EtaDefn%A  (nEtaLevels))
    EtaDefn%Name        = Name
    EtaDefn%EtaDefnType = 'B'
    EtaDefn%nEtaLevels  = nEtaLevels
    EtaDefn%Eta(:)      = Eta(:)
    EtaDefn%A(:)        = A(:)
  Else
    Call Message('UNEXPECTED FATAL ERROR in InitEtaDefn', 4)
  End If

End Function InitEtaDefn

!-------------------------------------------------------------------------------------------------------------

Subroutine AddEtaDefn(EtaDefn, EtaDefns)
! Adds an eta definition to a collection of eta definitions.

  Implicit None
  ! Argument list:
  Type(EtaDefn_),  Intent(In)    :: EtaDefn  ! Eta definition to be added.
  Type(EtaDefns_), Intent(InOut) :: EtaDefns ! Collection of eta definitions.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, EtaDefns%nEtaDefns
    If (EtaDefn%Name .CIEq. EtaDefns%EtaDefns(i)%Name) Then
      If (EtaDefn == EtaDefns%EtaDefns(i)) Then
        Return
      Else
        Call Message(                                                               &
               'FATAL ERROR in adding the eta definition "'                      // &
               Trim(EtaDefn%Name)                                                // &
               '": a different eta definition with the same name already exists',   &
               3                                                                    &
             )
      End If
    End If
  End Do

  If (EtaDefns%nEtaDefns >= MaxEtaDefns) Then
    Call Message(                                          &
           'FATAL ERROR in adding the eta definition "' // &
           Trim(EtaDefn%Name)                           // &
           '": there are too many eta definitions',        &
           3                                               &
         )
  End If

  EtaDefns%nEtaDefns                    = EtaDefns%nEtaDefns + 1
  EtaDefns%EtaDefns(EtaDefns%nEtaDefns) = EtaDefn

End Subroutine AddEtaDefn

!-------------------------------------------------------------------------------------------------------------

Function FindEtaDefnIndex(Name, EtaDefns)
! Finds the index of an eta definition.

  Implicit None
  ! Argument list:
  Character(*),    Intent(In) :: Name     ! Name of eta definition.
  Type(EtaDefns_), Intent(In) :: EtaDefns ! Collection of eta definitions.
  ! Function result:
  Integer :: FindEtaDefnIndex ! Index of the eta definition.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, EtaDefns%nEtaDefns
    If (Name .CIEq. EtaDefns%EtaDefns(i)%Name) Then
      FindEtaDefnIndex = i
      Return
    End If
  End Do

  Call Message(                             &
         'FATAL ERROR: eta definition "' // &
         Trim(Name)                      // &
         '" not found',                     &
         3                                  &
       )

End Function FindEtaDefnIndex

!-------------------------------------------------------------------------------------------------------------

Function EtaDefnEq(EtaDefn1, EtaDefn2)
! Tests for equality of eta definitions.

  Implicit None
  ! Argument list:
  Type(EtaDefn_), Intent(In) :: EtaDefn1 !} The two eta definitions.
  Type(EtaDefn_), Intent(In) :: EtaDefn2 !}
  ! Function result:
  Logical :: EtaDefnEq ! Indicates if eta definitions are equal.
  ! Locals:
  Integer :: i ! Loop index.

  EtaDefnEq = (EtaDefn1%Name        .CIEq. EtaDefn2%Name       ) .and. &
               EtaDefn1%nEtaLevels    ==   EtaDefn2%nEtaLevels   .and. &
               EtaDefn1%EtaDefnType   ==   EtaDefn2%EtaDefnType
  Do i = 1, Min(EtaDefn1%nEtaLevels, EtaDefn2%nEtaLevels)
    If (EtaDefn1%EtaDefnType == 'E') Then
      EtaDefnEq = EtaDefnEq                          .and. &
                  EtaDefn1%A(i)   == EtaDefn2%A(i)   .and. &
                  EtaDefn1%B(i)   == EtaDefn2%B(i)
    Else If (EtaDefn1%EtaDefnType == 'A') Then
      EtaDefnEq = EtaDefnEq                          .and. &
                  EtaDefn1%Eta(i) == EtaDefn2%Eta(i) .and. &
                  EtaDefn1%B(i)   == EtaDefn2%B(i)
    Else If (EtaDefn1%EtaDefnType == 'B') Then
      EtaDefnEq = EtaDefnEq                          .and. &
                  EtaDefn1%Eta(i) == EtaDefn2%Eta(i) .and. &
                  EtaDefn1%A(i)   == EtaDefn2%A(i)
    End If
  End Do

End Function EtaDefnEq

!-------------------------------------------------------------------------------------------------------------

Function InitHCoord(Name, CoordType, Pole, Angle, Origin, Unit, ThetaOrigin, &
                    ScaleFactor)
! Initialises a horizontal coord system.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: Name
  Integer,      Intent(In) :: CoordType
  Real(Std),    Intent(In) :: Pole(2)
  Real(Std),    Intent(In) :: Angle       ! $$ could be optional ?
  Real(Std),    Intent(In) :: Origin(2)
  Real(Std),    Intent(In) :: Unit(2)
  Real(Std),    Intent(In) :: ThetaOrigin ! $$ could be optional ?
  Real(Std),    Intent(In) :: ScaleFactor ! $$ could be optional ?
  ! Name        :: Name of coord system.
  ! CoordType   :: H_LatLong, H_PSCartesian, H_PSPolar, H_TMCartesian or H_TMPolar
  ! Pole(2)     :: In degrees
  ! Angle       :: In degrees
  ! Origin(2)   :: In degrees or m
  ! Unit(2)     :: In degrees or m
  ! ThetaOrigin :: In degrees
  ! ScaleFactor ::
  ! Function result:
  Type(HCoord_) :: InitHCoord ! Initialised horizontal coord system.
  ! Locals:
  Type(HCoord_) :: HCoord ! Local copy of function result.

  If (CoordType /= H_LatLong     .and. &
      CoordType /= H_PSCartesian .and. &
      CoordType /= H_PSPolar     .and. &
      CoordType /= H_TMCartesian .and. &
      CoordType /= H_TMPolar) Then
    Call Message('Error in InitHCoord: invalid CoordType', 4)
  Else If (Unit(1) == 0.0 .or.  &
           Unit(2) == 0.0) Then
    Call Message('Error in InitHCoord: Unit should be non-zero', 4)
  Else If (ScaleFactor <= 0.0) Then
    Call Message('Error in InitHCoord: ScaleFactor should be positive', 4)
  Else If (Len_Trim(Name) > MaxCharLength) Then
    Call Message('Error in InitHCoord: Name is too long', 4)
  End If

  HCoord%Name        = Name
  HCoord%CoordType   = CoordType
  HCoord%Pole        = Pole
  HCoord%Angle       = Angle
  HCoord%Origin      = Origin
  HCoord%Unit        = Unit
  HCoord%ThetaOrigin = ThetaOrigin
  HCoord%ScaleFactor = ScaleFactor

  ! Process angles as radians (input is in degrees).

  If (CoordType == H_LatLong     .or. &
      CoordType == H_PSCartesian .or. &
      CoordType == H_PSPolar) Then
    ! Convert pole latitude to angle from north pole.
    HCoord%Pole(2) = 90.0 - HCoord%Pole(2)
    ! Convert both pole and angle to radians.
    HCoord%Pole  = HCoord%Pole  * Pi / 180.0
    HCoord%Angle = HCoord%Angle * Pi / 180.0
  Else If (CoordType == H_TMCartesian .or. &
           CoordType == H_TMPolar) Then
    ! Convert pole (true origin) to radians.
    HCoord%Pole  = HCoord%Pole  * Pi / 180.0
  End If

  If (CoordType == H_LatLong) Then
    ! Convert origin offset and units to radians.
    HCoord%Origin = HCoord%Origin * Pi / 180.0
    HCoord%Unit   = HCoord%Unit   * Pi / 180.0
  Else If (CoordType == H_PSPolar .or. &
           CoordType == H_TMPolar) Then
    ! Convert theta unit and theta offset to radians.
    HCoord%Unit(2)     = HCoord%Unit(2) * Pi / 180.0
    HCoord%ThetaOrigin = HCoord%ThetaOrigin * Pi / 180.0
  End If

  ! Set: Angle = 0, ThetaOrigin = 0 and ScaleFactor = 1
  !  for coord systems where these components are redundant.
  If (CoordType == H_LatLong) Then
    HCoord%ThetaOrigin = 0.0 ; HCoord%ScaleFactor = 1.0
  Else If (CoordType == H_PSCartesian) Then
    HCoord%ThetaOrigin = 0.0 ; HCoord%ScaleFactor = 1.0
  Else If (CoordType == H_PSPolar) Then
    HCoord%ScaleFactor = 1.0
  Else If (CoordType == H_TMCartesian) Then
    HCoord%Angle = 0.0 ; HCoord%ThetaOrigin = 0.0
  Else If (CoordType == H_TMPolar) Then
    HCoord%Angle = 0.0
  End If

  If (HCoord%CoordType == H_LatLong) Then
    HCoord%XCycle = 2.0 * Pi / HCoord%Unit(1)
    HCoord%YCycle = 0.0
  Else If (HCoord%CoordType == H_PSCartesian) Then
    HCoord%XCycle = 0.0
    HCoord%YCycle = 0.0
  Else If (HCoord%CoordType == H_PSPolar) Then
    HCoord%XCycle = 0.0
    HCoord%YCycle = 2.0 * Pi / HCoord%Unit(2)
  Else If (HCoord%CoordType == H_TMCartesian) Then
    HCoord%XCycle = 0.0
    HCoord%YCycle = 0.0
  Else If (HCoord%CoordType == H_TMPolar) Then
    HCoord%XCycle = 0.0
    HCoord%YCycle = 2.0 * Pi / HCoord%Unit(2)
  End If

  If (HCoord%CoordType == H_LatLong) Then
    HCoord%YEdge = .true.
    HCoord%YEdgeMax = ( ( Pi / 2.0) - HCoord%Origin(2) ) / HCoord%Unit(2)
    HCoord%YEdgeMin = ( (-Pi / 2.0) - HCoord%Origin(2) ) / HCoord%Unit(2)
  Else
    HCoord%YEdge = .false.
    HCoord%YEdgeMax = 0.0
    HCoord%YEdgeMin = 0.0
  End If

  InitHCoord = HCoord

End Function InitHCoord

!-------------------------------------------------------------------------------------------------------------

Function HCoordEq(HCoord1, HCoord2)
! Tests for equality of horizontal coord systems (including their names).
  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoord1 ! First horizontal coord system.
  Type(HCoord_), Intent(In) :: HCoord2 ! Second horizontal coord system.
  ! Function result:
  Logical :: HCoordEq ! Indicates if coord systems are equal.

  HCoordEq = (HCoord1%Name      .CIEq. HCoord2%Name)        .and. &
              HCoord1%CoordType   ==   HCoord2%CoordType    .and. &
              HCoord1%Pole(1)     ==   HCoord2%Pole(1)      .and. &
              HCoord1%Pole(2)     ==   HCoord2%Pole(2)      .and. &
              HCoord1%Angle       ==   HCoord2%Angle        .and. &
              HCoord1%Origin(1)   ==   HCoord2%Origin(1)    .and. &
              HCoord1%Origin(2)   ==   HCoord2%Origin(2)    .and. &
              HCoord1%Unit(1)     ==   HCoord2%Unit(1)      .and. &
              HCoord1%Unit(2)     ==   HCoord2%Unit(2)      .and. &
              HCoord1%ThetaOrigin ==   HCoord2%ThetaOrigin  .and. &
              HCoord1%ScaleFactor ==   HCoord2%ScaleFactor

End Function HCoordEq

!-------------------------------------------------------------------------------------------------------------

Function HCoordEquiv(HCoord1, HCoord2)
! Tests for equivalence of horizontal coord systems (i.e. equality except possibly of
! their names).

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoord1 ! First horizontal coord system.
  Type(HCoord_), Intent(In) :: HCoord2 ! Second horizontal coord system.
  ! Function result:
  Logical :: HCoordEquiv ! Indicates if coord systems are equal.

  HCoordEquiv = HCoord1%CoordType   == HCoord2%CoordType    .and. &
                HCoord1%Pole(1)     == HCoord2%Pole(1)      .and. &
                HCoord1%Pole(2)     == HCoord2%Pole(2)      .and. &
                HCoord1%Angle       == HCoord2%Angle        .and. &
                HCoord1%Origin(1)   == HCoord2%Origin(1)    .and. &
                HCoord1%Origin(2)   == HCoord2%Origin(2)    .and. &
                HCoord1%Unit(1)     == HCoord2%Unit(1)      .and. &
                HCoord1%Unit(2)     == HCoord2%Unit(2)      .and. &
                HCoord1%ThetaOrigin == HCoord2%ThetaOrigin  .and. &
                HCoord1%ScaleFactor == HCoord2%ScaleFactor

End Function HCoordEquiv

!-------------------------------------------------------------------------------------------------------------

Function HCoord_LatLong()
! Returns the standard latitude-longitude coord system (with units in degrees).

  Implicit None
  ! Function result:
  Type(HCoord_) :: HCoord_LatLong ! Standard latitude-longitude coord
                                  ! system (with units in degrees).

  HCoord_LatLong = InitHCoord(                                    &
                     Name        = 'Lat-Long',                    &
                     CoordType   = H_LatLong,                     &
                     Pole        = (/ 0.0_Std, 90.0_Std /),       &
                     Angle       = 0.0_Std,                       &
                     Origin      = (/ 0.0_Std, 0.0_Std /),        &
                     Unit        = (/ 1.0_Std, 1.0_Std /),        &
                     ThetaOrigin = 0.0_Std,                       &
                     ScaleFactor = 1.0_Std                        &
                   )

End Function HCoord_LatLong

!-------------------------------------------------------------------------------------------------------------

Function HCoord_LatLongRadians()
! Returns the standard latitude-longitude coord system (with units in radians).

  Implicit None
  ! Function result:
  Type(HCoord_) :: HCoord_LatLongRadians ! Standard latitude-longitude coord
                                         ! system (with units in radians).

  HCoord_LatLongRadians  = InitHCoord(                                    &
                             Name        = 'Lat-Long (Radians)',          &
                             CoordType   = H_LatLong,                     &
                             Pole        = (/ 0.0_Std, 90.0_Std /),       &
                             Angle       = 0.0_Std,                       &
                             Origin      = (/ 0.0_Std, 0.0_Std /),        &
                             Unit        = (/ 180.0/Pi, 180.0/Pi /),      &
                             ThetaOrigin = 0.0_Std,                       &
                             ScaleFactor = 1.0_Std                        &
                           )
          ! $$ ensure unit exactly 1.0

End Function HCoord_LatLongRadians

!-------------------------------------------------------------------------------------------------------------

Function HCoord_EMEP50kmGrid()
! Returns the polar stereographic Cartesian coord system for the EMEP 50km grid.

  Implicit None
  ! Function result:
  Type(HCoord_) :: HCoord_EMEP50kmGrid ! EMEP coordinate system with 50km resolution.
  ! Locals:
  Real(Std) :: GridUnit ! unit length of EMEP grid x_50km = 53 589.84 metres.

  GridUnit = 53589.84_Std
  HCoord_EMEP50kmGrid = InitHCoord(                                           &
                          Name        = 'EMEP 50km Grid',                     &
                          CoordType   = H_PSCartesian,                        &
                          Pole        = (/ 0.0_Std, 90.0_Std /),              &
                          Angle       = -32.0_Std,                            &
                          Origin      = (/ -8.0*GridUnit, -110.0*GridUnit /), &
                          Unit        = (/ GridUnit, GridUnit /),             &
                          ThetaOrigin = 0.0_Std,                              &
                          ScaleFactor = 1.0_Std                               &
                        )

End Function HCoord_EMEP50kmGrid

!-------------------------------------------------------------------------------------------------------------

Function HCoord_EMEP150kmGrid()
! Returns the polar stereographic Cartesian coord system for the EMEP 150km grid.

  Implicit None
  ! Function result:
  Type(HCoord_) :: HCoord_EMEP150kmGrid ! EMEP coordinate system with 150km resolution.
  ! Locals:
  Real(Std) :: GridUnit ! unit length of EMEP grid x_150km = 160 769.52 metres.

  GridUnit = 160769.52_Std
  HCoord_EMEP150kmGrid = InitHCoord(                                          &
                           Name        = 'EMEP 150km Grid',                   &
                           CoordType   = H_PSCartesian,                       &
                           Pole        = (/ 0.0_Std, 90.0_Std /),             &
                           Angle       = -32.0_Std,                           &
                           Origin      = (/ -3.0*GridUnit, -37.0*GridUnit /), &
                           Unit        = (/ GridUnit, GridUnit /),            &
                           ThetaOrigin = 0.0_Std,                             &
                           ScaleFactor = 1.0_Std                              &
                         )

End Function HCoord_EMEP150kmGrid

!-------------------------------------------------------------------------------------------------------------

Function HCoord_UKNationalGridM()
! Returns the UK National Grid coord system (with units in metres).

  Implicit None
  ! Function result:
  Type(HCoord_) :: HCoord_UKNationalGridM ! UK National Grid system (metre units).

  HCoord_UKNationalGridM = InitHCoord(                                     &
                             Name        = 'UK National Grid (m)',         &
                             CoordType   = H_TMCartesian,                  &
                             Pole        = (/ -2.0_Std, 49.0_Std /),       &
                             Angle       = 0.0_Std,                        &
                             Origin      = (/ -4.0E5_Std, 1.0E5_Std /),    &
                             Unit        = (/ 1.0_Std, 1.0_Std /),         &
                             ThetaOrigin = 0.0_Std,                        &
                             ScaleFactor = 0.9996012717_Std                &
                           )

End Function HCoord_UKNationalGridM

!-------------------------------------------------------------------------------------------------------------

Function HCoord_UKNationalGrid100M()
! Returns the UK National Grid coord system (with 100 metre units).

  Implicit None
  ! Function result:
  Type(HCoord_) :: HCoord_UKNationalGrid100M ! UK National Grid system (100 metre units).

  HCoord_UKNationalGrid100M = InitHCoord(                                         &
                                Name        = 'UK National Grid (100m)',          &
                                CoordType   = H_TMCartesian,                      &
                                Pole        = (/ -2.0_Std, 49.0_Std /),           &
                                Angle       = 0.0_Std,                            &
                                Origin      = (/ -4.0E5_Std, 1.0E5_Std /),        &
                                Unit        = (/ 100.0_Std, 100.0_Std /),         &
                                ThetaOrigin = 0.0_Std,                            &
                                ScaleFactor = 0.9996012717_Std                    &
                              )

End Function HCoord_UKNationalGrid100M

!-------------------------------------------------------------------------------------------------------------

Function LatLongToSLatLong(HCoordIn, HCoordOut, PointIn)
! Converts between two lat-long coord systems, the second one being the standard
! latitude-longitude coord system (apart from a possible origin offset and choice of
! units).

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   !
  Type(HCoord_), Intent(In) :: HCoordOut  !
  Real(Pos),     Intent(In) :: PointIn(2) !
  ! Function result:
  Real(Pos) :: LatLongToSLatLong(2) !
  ! Locals:
  Real(Std), Parameter :: Small = 1.0E-15_Std ! For cos(lat) less than this, treat as
                                              ! if at pole and set longitude = 0.
  Real(Std) :: LatIn       ! Input latitude in radians.
  Real(Std) :: LongIn      ! Input longitude in radians.
  Real(Std) :: LatPole     ! Latitude of input pole in radians.
  Real(Std) :: LongPole    ! Longitude of input pole in radians.
  Real(Pos) :: PointOut(2) ! Output coords.
  Real(Std) :: XTemp       !} Coords in cartesian system centred on the centre of the
  Real(Std) :: YTemp       !} Earth.
  Real(Std) :: ZTemp       !}

# ifdef ExtraChecks
    If (HCoordIn%CoordType  /= H_LatLong .or. &
        HCoordOut%CoordType /= H_LatLong .or. &
        HCoordOut%Pole(1)   /= 0.0       .or. &
        HCoordOut%Pole(2)   /= 0.0       .or. &
        HCoordOut%Angle     /= 0.0) Then
      Call Message('Error in LatLongToSLatLong', 4)
    End If
# endif

  ! Scaling.
  LongIn   = PointIn(1) * HCoordIn%Unit(1) + HCoordIn%Origin(1) + HCoordIn%Angle
  LatIn    = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%Origin(2)
  LongPole = HCoordIn%Pole(1)
  LatPole  = Pi/2.0 - HCoordIn%Pole(2)

  ! Coords in cartesian system centred on the centre of the Earth.
  XTemp = Cos(LongIn)*Cos(LatIn)*Sin(LatPole) + Sin(LatIn)*Cos(LatPole)
  YTemp = Sin(LongIn)*Cos(LatIn)
  ZTemp = Sin(LatIn)*Sin(LatPole) - Cos(LatIn)*Cos(LatPole)*Cos(LongIn)

  ! Latitude calculation.
  PointOut(2) = ATan2(ZTemp, Sqrt(XTemp**2 + YTemp**2))

  ! Longitude calculation.
  If ((Pi/2.0 - Abs(PointOut(2))) < Small) Then
    PointOut(1) = 0.0
  Else
    PointOut(1) = ATan2(YTemp, XTemp) + LongPole
    If (PointOut(1) >   Pi) PointOut(1) = PointOut(1) - 2.0*Pi
    If (PointOut(1) < - Pi) PointOut(1) = PointOut(1) + 2.0*Pi
  End If

  ! Scaling.
  PointOut(1) = (PointOut(1) - HCoordOut%Origin(1)) / HCoordOut%Unit(1)
  PointOut(2) = (PointOut(2) - HCoordOut%Origin(2)) / HCoordOut%Unit(2)

  LatLongToSLatLong = PointOut

End Function LatLongToSLatLong

!-------------------------------------------------------------------------------------------------------------

Function SLatLongToLatLong(HCoordIn, HCoordOut, PointIn)
! Converts between two lat-long coord systems, the first one being the standard
! latitude-longitude coord system (apart from a possible scale factor).

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   !
  Type(HCoord_), Intent(In) :: HCoordOut  !
  Real(Pos),     Intent(In) :: PointIn(2) !
  ! Function result:
  Real(Pos) :: SLatLongToLatLong(2) !
  ! Locals:
  Real(Std), Parameter :: Small = 1.0E-15_Std ! For cos(lat) less than this, treat as
                                              ! if at pole and set longitude = 0.
  Real(Std) :: LatIn       ! Input latitude in radians.
  Real(Std) :: LongIn      ! Input longitude in radians.
  Real(Std) :: LatPole     ! Latitude of output pole in radians.
  Real(Std) :: LongPole    ! Longitude of output pole in radians.
  Real(Pos) :: PointOut(2) ! Output coords.
  Real(Std) :: XTemp       !} Coords in cartesian system centred on the centre of the
  Real(Std) :: YTemp       !} Earth.
  Real(Std) :: ZTemp       !}

# ifdef ExtraChecks
    If (HCoordOut%CoordType /= H_LatLong .or. &
        HCoordIn%CoordType  /= H_LatLong .or. &
        HCoordIn%Pole(1)    /= 0.0       .or. &
        HCoordIn%Pole(2)    /= 0.0       .or. &
        HCoordIn%Angle      /= 0.0) Then
      Call Message('Error in SLatLongToLatLong', 4)
    End If
# endif

  ! Scaling.
  LongIn   = PointIn(1) * HCoordIn%Unit(1) + HCoordIn%Origin(1)
  LatIn    = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%Origin(2)
  LongPole = HCoordOut%Pole(1)
  LatPole  = Pi/2.0 - HCoordOut%Pole(2)

  ! Coords in cartesian system centred on the centre of the Earth.
  XTemp = Cos(LongIn - LongPole)*Cos(LatIn)*Sin(LatPole) - Sin(LatIn)*Cos(LatPole)
  YTemp = Sin(LongIn - LongPole)*Cos(LatIn)
  ZTemp = Sin(LatIn)*Sin(LatPole) + Cos(LatIn)*Cos(LatPole)*Cos(LongIn - LongPole)

  ! Latitude calculation.
  PointOut(2) = ATan2(ZTemp, Sqrt(XTemp**2 + YTemp**2))

  ! Longitude calculation.
  If ((Pi/2.0 - Abs(PointOut(2))) < Small) Then
    PointOut(1) = 0.0
  Else
    PointOut(1) = ATan2(YTemp, XTemp)
  End If

  ! Angle.
  PointOut(1) = PointOut(1) - HCoordOut%Angle

  If (PointOut(1) >   Pi) PointOut(1) = PointOut(1) - 2.0*Pi
  If (PointOut(1) < - Pi) PointOut(1) = PointOut(1) + 2.0*Pi

  ! Scaling.
  PointOut(1) = (PointOut(1) - HCoordOut%Origin(1)) / HCoordOut%Unit(1)
  PointOut(2) = (PointOut(2) - HCoordOut%Origin(2)) / HCoordOut%Unit(2)

  SLatLongToLatLong = PointOut

End Function SLatLongToLatLong

!-------------------------------------------------------------------------------------------------------------

Function LatLongToLatLong(HCoordIn, HCoordOut, PointIn)
! Converts between two lat-long coord systems.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: LatLongToLatLong(2) ! Output coords.
  ! Locals:
  Real(Std), Parameter :: Small = 1.0E-15_Std ! For cos(lat) less than Small, treat
                                              ! as if at pole and set longitude = 0.
  Real(Std) :: LatIn       ! Input latitude relative to input pole (in radians).
  Real(Std) :: LongIn      ! Input longitude relative to input pole (in radians).
  Real(Std) :: LatOut      ! Output latitude relative to output pole (in radians).
  Real(Std) :: LongOut     ! Output longitude relative to output pole (in radians).
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: X1, Y1, Z1  !} Coords in various intermediate 3-d cartesian systems
  Real(Std) :: X2, Y2, Z2  !} centred on the centre of the Earth.
  Real(Std) :: X3, Y3, Z3  !}
  Real(Std) :: X4, Y4, Z4  !}
  Real(Std) :: Rot1        !] Rotation angles in the coordinate transformation.
  Real(Std) :: Rot2        !]
  Real(Std) :: Rot3        !]

# ifdef ExtraChecks
    If (HCoordIn%CoordType  /= H_LatLong .or. &
        HCoordOut%CoordType /= H_LatLong) Then
      Call Message('Error in LatLongToLatLong', 4)
    End If
# endif

  ! Processing the input coords to lat/long in radians.
  LongIn   = PointIn(1) * HCoordIn%Unit(1) + HCoordIn%Origin(1) + HCoordIn%Angle
  LatIn    = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%Origin(2)

  ! Defining the rotation angles of the pole transformation.
  Rot1 = HCoordIn%Pole(2)
  Rot2 = HCoordOut%Pole(1) - HCoordIn%Pole(1)
  Rot3 = HCoordOut%Pole(2)

  ! Coords in 3-d cartesian system aligned with the input pole.
  X1 = Cos(LatIn) * Cos(LongIn)
  Y1 = Cos(LatIn) * Sin(LongIn)
  Z1 = Sin(LatIn)

  ! Coords in cartesian system rotated to true north and longitude of input pole.
  X2 = Cos(Rot1) * X1 + Sin(Rot1) * Z1
  Y2 = Y1
  Z2 = Cos(Rot1) * Z1 - Sin(Rot1) * X1

  ! Coords in cartesian system rotated to true north and longitude of output pole.
  X3 = Cos(Rot2) * X2 + Sin(Rot2) * Y2
  Y3 = Cos(Rot2) * Y2 - Sin(Rot2) * X2
  Z3 = Z2

  ! Coords in 3-d cartesian system aligned with the output pole.
  X4 = Cos(Rot3) * X3 - Sin(Rot3) * Z3
  Y4 = Y3
  Z4 = Sin(Rot3) * X3 + Cos(Rot3) * Z3

  ! Latitude calculation.
  LatOut = ATan2(Z4, Sqrt(X4**2 + Y4**2))

  ! Longitude calculation.
  If ((Pi/2.0 - Abs(LatOut)) < Small) Then
    LongOut = 0.0
  Else
    LongOut = ATan2(Y4, X4)
  End If

  ! Apply angle offset of output coord system.
  LongOut = LongOut - HCoordOut%Angle

  If (LongOut >   Pi) LongOut = LongOut - 2.0*Pi
  If (LongOut < - Pi) LongOut = LongOut + 2.0*Pi

  ! Scaling the lat/long position in radians to the output coord system.
  PointOut(1) = (LongOut - HCoordOut%Origin(1)) / HCoordOut%Unit(1)
  PointOut(2) = (LatOut  - HCoordOut%Origin(2)) / HCoordOut%Unit(2)

  LatLongToLatLong = PointOut

End Function LatLongToLatLong

!-------------------------------------------------------------------------------------------------------------

Function QuickLatLongToLatLong(HCoordIn, HCoordOut, PointIn)
! Converts between two lat-long coord systems where the coord systems have the same
! value for Pole.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: QuickLatLongToLatLong(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: Point(2)    ! True coords in each lat-long system (in radians).

# ifdef ExtraChecks
    If (HCoordIn%CoordType  /= H_LatLong         .or. &
        HCoordOut%CoordType /= H_LatLong         .or. &
        HCoordIn%Pole(1)    /= HCoordOut%Pole(1) .or. &
        HCoordIn%Pole(2)    /= HCoordOut%Pole(2)) Then
      Call Message('Error in QuickLatLongToLatLong', 4)
    End If
# endif

  ! Scaling of input coords.
  Point(1) = PointIn(1) * HCoordIn%Unit(1) + HCoordIn%Origin(1)
  Point(2) = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%Origin(2)

  ! Rotation angle adjustment.
  Point(1) = Point(1) + HCoordIn%Angle - HCoordOut%Angle
  If (Point(1) >   Pi) Point(1) = Point(1) - 2.0*Pi
  If (Point(1) < - Pi) Point(1) = Point(1) + 2.0*Pi

  ! Scaling of output coords.
  PointOut(1) = (Point(1) - HCoordOut%Origin(1)) / HCoordOut%Unit(1)
  PointOut(2) = (Point(2) - HCoordOut%Origin(2)) / HCoordOut%Unit(2)

  QuickLatLongToLatLong = PointOut

End Function QuickLatLongToLatLong

!-------------------------------------------------------------------------------------------------------------

Function LatLongToPSCartesian(HCoordIn, HCoordOut, PointIn)
! Converts from a lat-long coord system to a polar stereographic Cartesian coord
! system where the coord systems have the same value for Pole.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: LatLongToPSCartesian(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: U           !
  Real(Std) :: V           !
  Real(Std) :: UPrime      !
  Real(Std) :: Point(2)    !

# ifdef ExtraChecks
    If (HCoordIn%CoordType  /= H_LatLong         .or. &
        HCoordOut%CoordType /= H_PSCartesian     .or. &
        HCoordIn%Pole(1)    /= HCoordOut%Pole(1) .or. &
        HCoordIn%Pole(2)    /= HCoordOut%Pole(2)) Then
      Call Message('Error in LatLongToPSCartesian', 4)
    End If
# endif

  ! Scaling of input coords.
  Point(1) = PointIn(1) * HCoordIn%Unit(1) + HCoordIn%Origin(1)
  Point(2) = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%Origin(2)

  ! Rotation of reference angle.
  Point(1) = Point(1) + HCoordIn%Angle - HCoordOut%Angle

  ! Calculation of coordinate position on tangent plane.
  U           = Cos(Point(2))
  V           = Sin(Point(2))
  UPrime      = 2.0 * EarthRadius * U / (1.0 + V)
  PointOut(1) =   UPrime * Sin(Point(1))
  PointOut(2) = - UPrime * Cos(Point(1))

  ! Scaling of output coords.
  PointOut(1) = (PointOut(1) - HCoordOut%Origin(1)) / HCoordOut%Unit(1)
  PointOut(2) = (PointOut(2) - HCoordOut%Origin(2)) / HCoordOut%Unit(2)

  LatLongToPSCartesian = PointOut

End Function LatLongToPSCartesian

!-------------------------------------------------------------------------------------------------------------

Function PSCartesianToLatLong(HCoordIn, HCoordOut, PointIn)
! Converts from a polar stereographic Cartesian coord system to a lat-long coord
! system where the coord systems have the same value for Pole.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: PSCartesianToLatLong(2) ! Output coords.
  ! Locals:
  Real(Std), Parameter :: Small = 1.0E-15_Std ! If angle to pole is less than this,
                                              ! treat as at pole and set longitude = 0
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: U           !
  Real(Std) :: V           !
  Real(Std) :: Point(2)    !

# ifdef ExtraChecks
    If (HCoordOut%CoordType /= H_LatLong         .or. &
        HCoordIn%CoordType  /= H_PSCartesian     .or. &
        HCoordIn%Pole(1)    /= HCoordOut%Pole(1) .or. &
        HCoordIn%Pole(2)    /= HCoordOut%Pole(2)) Then
      Call Message('Error in PSCartesianToLatLong', 4)
    End If
# endif

  ! Scaling of input coords.
  Point(1) = PointIn(1) * HCoordIn%Unit(1) + HCoordIn%Origin(1)
  Point(2) = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%Origin(2)

  V           = Point(1)**2 + Point(2)**2
  V           = V / (2.0 * EarthRadius)**2
  V           = (1.0 - V) / (1.0 + V)
  U           = Sqrt(1.0 - V**2)

  ! Latitude calculation.
  PointOut(2) = ATan2(V, U)

  ! Longitude calculation.
  If ((Pi/2.0 - Abs(PointOut(2))) < Small) Then
    PointOut(1) = 0.0
  Else
    PointOut(1) = ATan2(Point(1), - Point(2))
    ! Rotation of reference angle.
    PointOut(1) = PointOut(1) + HCoordIn%Angle - HCoordOut%Angle
    If (PointOut(1) >  Pi) PointOut(1) = PointOut(1) - 2.0*Pi
    If (PointOut(1) < -Pi) PointOut(1) = PointOut(1) + 2.0*Pi
  End If

  ! Scaling of output coords.
  PointOut(1) = (PointOut(1) - HCoordOut%Origin(1)) / HCoordOut%Unit(1)
  PointOut(2) = (PointOut(2) - HCoordOut%Origin(2)) / HCoordOut%Unit(2)

  PSCartesianToLatLong = PointOut

End Function PSCartesianToLatLong

!-------------------------------------------------------------------------------------------------------------

Function QuickPSCartesianToPSCartesian(HCoordIn, HCoordOut, PointIn)
! Converts between two polar stereographic Cartesian coord systems where the
! coord systems have the same value for Pole.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: QuickPSCartesianToPSCartesian(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: Point(2)    ! Intermediate coords.
  Real(Std) :: Angle       ! Rotation angle between the two cartesian systems.

# ifdef ExtraChecks
    If (HCoordIn%CoordType  /= H_PSCartesian     .or. &
        HCoordOut%CoordType /= H_PSCartesian     .or. &
        HCoordIn%Pole(1)    /= HCoordOut%Pole(1) .or. &
        HCoordIn%Pole(2)    /= HCoordOut%Pole(2)) Then
      Call Message('Error in QuickPSCartesianToPSCartesian', 4)
    End If
# endif

  ! Scaling of input coords.
  Point(1) = PointIn(1) * HCoordIn%Unit(1) + HCoordIn%Origin(1)
  Point(2) = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%Origin(2)

  ! Rotation angle adjustment.
  If (HCoordOut%Angle /= HCoordIn%Angle) Then
    Angle       = HCoordOut%Angle - HCoordIn%Angle
    PointOut(1) = Point(1) * Cos(Angle) + Point(2) * Sin(Angle)
    PointOut(2) = Point(2) * Cos(Angle) - Point(1) * Sin(Angle)
  Else
    PointOut = Point
  End If

  ! Scaling of output coords.
  PointOut(1) = (PointOut(1) - HCoordOut%Origin(1)) / HCoordOut%Unit(1)
  PointOut(2) = (PointOut(2) - HCoordOut%Origin(2)) / HCoordOut%Unit(2)

  QuickPSCartesianToPSCartesian = PointOut

End Function QuickPSCartesianToPSCartesian

!-------------------------------------------------------------------------------------------------------------

Function LatLongToPSPolar(HCoordIn, HCoordOut, PointIn)
! Converts from a lat-long coord system to a polar system in a polar stereographic
! projection where the coord systems have the same value for Pole.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: LatLongToPSPolar(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: Point(2)    ! Lat-long of point.
  Real(Std) :: U           ! Cos(lat).
  Real(Std) :: V           ! Sin(lat).
  Real(Std) :: R           ! Radial distance from tangent point.
  Real(Std) :: Theta       ! Direction from tangent point.
  Real(Std) :: X           ! } Cartesian coords relative to polar origin.
  Real(Std) :: Y           ! }
  Real(Std) :: RP          ! Radial distance from polar origin.
  Real(Std) :: ThetaP      ! Direction from polar origin.

# ifdef ExtraChecks
    If (HCoordIn%CoordType  /= H_LatLong         .or. &
        HCoordOut%CoordType /= H_PSPolar         .or. &
        HCoordIn%Pole(1)    /= HCoordOut%Pole(1) .or. &
        HCoordIn%Pole(2)    /= HCoordOut%Pole(2)) Then
      Call Message('Error in LatLongToPSPolar', 4)
    End If
# endif

  ! Scaling of input coords.
  Point(1) = PointIn(1) * HCoordIn%Unit(1) + HCoordIn%Origin(1)
  Point(2) = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%Origin(2)

  ! Compute polar stereographic projection on the tangent plane.
  U     = Cos(Point(2))
  V     = Sin(Point(2))
  R     = 2.0 * EarthRadius * U / (1.0 + V)
  Theta = Point(1) - Pi/2.0

  ! Change of reference angle between coord systems.
  Theta = Theta + HCoordIn%Angle - HCoordOut%Angle

  ! Compute coords in associated Cartesian system (relative to polar origin).
  X = R * Cos(Theta) - HCoordOut%Origin(1)
  Y = R * Sin(Theta) - HCoordOut%Origin(2)

  ! Compute polar coords relative to polar origin.
  RP     = Sqrt(X**2 + Y**2)
  ThetaP = ATan2ZeroTest(Y, X) - HCoordOut%ThetaOrigin
  If (ThetaP >   Pi) ThetaP = ThetaP - 2.0*Pi
  If (ThetaP < - Pi) ThetaP = ThetaP + 2.0*Pi

  ! Scaling of output coords.
  PointOut(1) =  RP / HCoordOut%Unit(1)
  PointOut(2) =  ThetaP / HCoordOut%Unit(2)

  LatLongToPSPolar = PointOut

End Function LatLongToPSPolar

!-------------------------------------------------------------------------------------------------------------

Function PSPolarToLatLong(HCoordIn, HCoordOut, PointIn)
! Converts from a polar coord system in a polar stereographic projection to
! a lat-long system where the coord systems have the same value for Pole.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: PSPolarToLatLong(2) ! Output coords.
  ! Locals:
  Real(Std), Parameter :: Small = 1.0E-15_Std ! If angle to pole is less than this,
                                              ! treat as at pole and set longitude = 0
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: RP          ! Radial distance from polar origin.
  Real(Std) :: ThetaP      ! Direction from polar origin.
  Real(Std) :: X           ! } Cartesian coords relative to tangent point.
  Real(Std) :: Y           ! }
  Real(Std) :: R2          ! Radial distance (squared) from tangent point.
  Real(Std) :: U           ! Cos(lat).
  Real(Std) :: V           ! Sin(lat).

# ifdef ExtraChecks
    If (HCoordOut%CoordType /= H_LatLong         .or. &
        HCoordIn%CoordType  /= H_PSPolar         .or. &
        HCoordIn%Pole(1)    /= HCoordOut%Pole(1) .or. &
        HCoordIn%Pole(2)    /= HCoordOut%Pole(2)) Then
      Call Message('Error in PSPolarToLatLong', 4)
    End If
# endif

  ! Scaling of input coords.
  RP     = PointIn(1) * HCoordIn%Unit(1)
  ThetaP = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%ThetaOrigin

  ! Compute coords in associated Cartesian system (relative to tangent point).
  X = RP * Cos(ThetaP) + HCoordIn%Origin(1)
  Y = RP * Sin(ThetaP) + HCoordIn%Origin(2)

  ! Radial distance (squared) from tangent point.
  R2 = X**2 + Y**2

  ! Compute lat-long of point.
  V = R2 / (2.0 * EarthRadius)**2
  V = (1.0 - V) / (1.0 + V)       ! Sine of latitude
  U = Sqrt(1.0 - V**2)            ! Cosine of latitude

  ! Latitude calculation.
  PointOut(2) = ATan2(V, U)

  ! Longitude calculation.
  If ((Pi/2.0 - Abs(PointOut(2))) < Small) Then
    PointOut(1) = 0.0
  Else
    PointOut(1) = ATan2(X, -Y)
  ! Rotation of reference angle between coord systems.
    PointOut(1) = PointOut(1) + HCoordIn%Angle - HCoordOut%Angle
    If (PointOut(1) >  Pi) PointOut(1) = PointOut(1) - 2.0*Pi
    If (PointOut(1) < -Pi) PointOut(1) = PointOut(1) + 2.0*Pi
  End If

  ! Scaling of output coords.
  PointOut(1) = (PointOut(1) - HCoordOut%Origin(1)) / HCoordOut%Unit(1)
  PointOut(2) = (PointOut(2) - HCoordOut%Origin(2)) / HCoordOut%Unit(2)

  PSPolarToLatLong = PointOut

End Function PSPolarToLatLong

!-------------------------------------------------------------------------------------------------------------

Function PSCartesianToPSPolar(HCoordIn, HCoordOut, PointIn)
! Converts from a cartesian coord system to a polar system in a polar stereographic
! projection where the coord systems have the same value for Pole.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: PSCartesianToPSPolar(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: Point(2)    ! } Intermediate coords.
  Real(Std) :: Point2(2)   ! }
  Real(Std) :: Angle       ! Rotation angle.

# ifdef ExtraChecks
    If (HCoordOut%CoordType /= H_PSPolar         .or. &
        HCoordIn%CoordType  /= H_PSCartesian     .or. &
        HCoordIn%Pole(1)    /= HCoordOut%Pole(1) .or. &
        HCoordIn%Pole(2)    /= HCoordOut%Pole(2)) Then
      Call Message('Error in PSCartesianToPSPolar', 4)
    End If
# endif

  ! Scaling of input coords.
  Point(1) = PointIn(1) * HCoordIn%Unit(1) + HCoordIn%Origin(1)
  Point(2) = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%Origin(2)

  ! Rotation of reference angle between coord systems.
  If (HCoordOut%Angle /= HCoordIn%Angle) Then
    Angle     = HCoordOut%Angle - HCoordIn%Angle
    Point2(1) = Point(1) * Cos(Angle) + Point(2) * Sin(Angle)
    Point2(2) = Point(2) * Cos(Angle) - Point(1) * Sin(Angle)
    Point     = Point2
  End If

  ! Conversion to polar coords relative to polar origin.
  Point(1)    = Point(1) - HCoordOut%Origin(1)
  Point(2)    = Point(2) - HCoordOut%Origin(2)
  PointOut(1) = Sqrt(Point(1)**2 + Point(2)**2)
  PointOut(2) = ATan2ZeroTest(Point(2), Point(1)) - HCoordOut%ThetaOrigin

  ! Scaling of output coords.
  PointOut(1) = PointOut(1) / HCoordOut%Unit(1)
  PointOut(2) = PointOut(2) / HCoordOut%Unit(2)

  PSCartesianToPSPolar = PointOut

End Function PSCartesianToPSPolar

!-------------------------------------------------------------------------------------------------------------

Function PSPolarToPSCartesian(HCoordIn, HCoordOut, PointIn)
! Converts from a polar coord system to a cartesian system in a polar stereographic
! projection where the coord systems have the same value for Pole.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: PSPolarToPSCartesian(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: R           ! } Scaled input coords.
  Real(Std) :: Theta       ! }
  Real(Std) :: Point(2)    ! ] Intermediate coords.
  Real(Std) :: Point2(2)   ! ]
  Real(Std) :: Angle       ! Rotation angle.

# ifdef ExtraChecks
    If (HCoordOut%CoordType /= H_PSCartesian     .or. &
        HCoordIn%CoordType  /= H_PSPolar         .or. &
        HCoordIn%Pole(1)    /= HCoordOut%Pole(1) .or. &
        HCoordIn%Pole(2)    /= HCoordOut%Pole(2)) Then
      Call Message('Error in PSPolarToPSCartesian', 4)
    End If
# endif

  ! Scaling of input coords.
  R     = PointIn(1) * HCoordIn%Unit(1)
  Theta = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%ThetaOrigin

  ! Conversion to cartesian coords relative to tangent point.
  Point(1) = R * Cos(Theta) + HCoordIn%Origin(1)
  Point(2) = R * Sin(Theta) + HCoordIn%Origin(2)

  ! Rotation of reference angle between coord systems.
  If (HCoordOut%Angle /= HCoordIn%Angle) Then
    Angle     = HCoordOut%Angle - HCoordIn%Angle
    Point2(1) = Point(1) * Cos(Angle) + Point(2) * Sin(Angle)
    Point2(2) = Point(2) * Cos(Angle) - Point(1) * Sin(Angle)
    Point     = Point2
  End If

  ! Scaling of output coords.
  PointOut(1) = (Point(1) - HCoordOut%Origin(1)) / HCoordOut%Unit(1)
  PointOut(2) = (Point(2) - HCoordOut%Origin(2)) / HCoordOut%Unit(2)

  PSPolarToPSCartesian = PointOut

End Function PSPolarToPSCartesian

!-------------------------------------------------------------------------------------------------------------

Function QuickPSPolarToPSPolar(HCoordIn, HCoordOut, PointIn)
! Converts between two polar coord systems in a polar stereographic projection
! where the coord systems have the same value for Pole, Angle and Origin.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn    ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut   ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2)  ! Input coords.
  ! Function result:
  Real(Pos) :: QuickPSPolarToPSPolar(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: Point(2)    ! Scaled coords in standard polar units.

# ifdef ExtraChecks
    If (HCoordIn%CoordType  /= H_PSPolar           .or. &
        HCoordOut%CoordType /= H_PSPolar           .or. &
        HCoordIn%Pole(1)    /= HCoordOut%Pole(1)   .or. &
        HCoordIn%Pole(2)    /= HCoordOut%Pole(2)   .or. &
        HCoordIn%Angle      /= HCoordOut%Angle     .or. &
        HCoordIn%Origin(1)  /= HCoordOut%Origin(1) .or. &
        HCoordIn%Origin(2)  /= HCoordOut%Origin(2)) Then
      Call Message('Error in QuickPSPolarToPSPolar', 4)
    End If
# endif

  ! Scaling of input coords.
  Point(1) = PointIn(1) * HCoordIn%Unit(1)
  Point(2) = PointIn(2) * HCoordIn%Unit(2)

  ! Rotation of theta offset between coord systems.
  Point(2) = Point(2) + HCoordIn%ThetaOrigin - HCoordOut%ThetaOrigin
  If (Point(2) >   Pi) Point(2) = Point(2) - 2.0*Pi
  If (Point(2) < - Pi) Point(2) = Point(2) + 2.0*Pi

  ! Scaling of output coords.
  PointOut(1) = Point(1) / HCoordOut%Unit(1)
  PointOut(2) = Point(2) / HCoordOut%Unit(2)

  QuickPSPolarToPSPolar = PointOut

End Function QuickPSPolarToPSPolar

!-------------------------------------------------------------------------------------------------------------

Function SLatLongToTMCartesian(HCoordIn, HCoordOut, PointIn)
! Converts from the standard lat-long coord system (in radians) to a Cartesian coord
! system in a transverse Mercator projection.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: SLatLongToTMCartesian(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: Lat1        ! Latitude in offset system.
  Real(Std) :: Long1       ! Longitude in offset system.
  Real(Std) :: S           ! Intermediate expression equal to sin(pseudo-latitude).
  Real(Std) :: C           ! Intermediate expression equal to cos(pseudo-latitude).
  Real(Std) :: PLat        ! Pseudo-latitude (relative to rotated N pole).
  Real(Std) :: PLong       ! Pseudo-longitude (relative to rotated N pole).
  Real(Std) :: X           ! Normalised distance in the East direction.
  Real(Std) :: Y           ! Normalised distance in the North direction.
  Real(Std) :: DistanceCon ! Conversion constant for grid distances.

# ifdef ExtraChecks
    If (HCoordOut%CoordType /= H_TMCartesian .or. &
        HCoordIn%CoordType  /= H_LatLong     .or. &
        HCoordIn%Pole(1)    /= 0.0           .or. &
        HCoordIn%Pole(2)    /= 0.0           .or. &
        HCoordIn%Angle      /= 0.0           .or. &
        HCoordIn%Origin(1)  /= 0.0           .or. &
        HCoordIn%Origin(2)  /= 0.0           .or. &
        HCoordIn%Unit(1)    /= 1.0           .or. &
        HCoordIn%Unit(2)    /= 1.0) Then
      Call Message('Error in SLatLongToTMCartesian', 4)
    End If
# endif

  ! Latitude and longitude in the offset lat-long system.
  Long1 = PointIn(1) - HCoordOut%Pole(1)
  Lat1  = PointIn(2)

  ! Calculate pseudo-latitude and pseudo-longitude (relative to rotated pole).
  S     = Cos(Lat1) * Sin(Long1)
  C     = (Cos(Lat1) * Cos(Long1))**2 + Sin(Lat1)**2
  C     = Sqrt(C)
  PLat  = ATan2(S, C)
  PLong = ATan2ZeroTest( - Sin(Lat1), Cos(Lat1) * Cos(Long1))

  ! Calculate normalised distances from true origin.
  X = Log(Tan(PLat/2.0 + Pi/4.0))
  Y = - (PLong + HCoordOut%Pole(2))
  If (Y >   Pi) Y = Y - 2.0*Pi
  If (Y < - Pi) Y = Y + 2.0*Pi

  ! Conversion constant for calculating actual grid distances.
  DistanceCon = HCoordOut%ScaleFactor * EarthRadius

  ! The Easting and Northing grid coordinate.
  PointOut(1) = (X * DistanceCon - HCoordOut%Origin(1))/HCoordOut%Unit(1)
  PointOut(2) = (Y * DistanceCon - HCoordOut%Origin(2))/HCoordOut%Unit(2)

  SLatLongToTMCartesian = PointOut

End Function SLatLongToTMCartesian

!-------------------------------------------------------------------------------------------------------------

Function TMCartesianToSLatLong(HCoordIn, HCoordOut, PointIn)
! Converts from a Cartesian coord system in a transverse Mercator projection to the
! standard lat-long coord system (in radians).

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: TMCartesianToSLatLong(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: DistanceCon ! Conversion constant for grid distances.
  Real(Std) :: X           ! Normalised distance in the East direction.
  Real(Std) :: Y           ! Normalised distance in the North direction.
  Real(Std) :: PLat        ! Pseudo-latitude (relative to rotated N pole).
  Real(Std) :: PLong       ! Pseudo-longitude (relative to rotated N pole).
  Real(Std) :: S           ! Intermediate expression equal to sin(latitude).
  Real(Std) :: C           ! Intermediate expression equal to cos(latitude).

# ifdef ExtraChecks
    If (HCoordIn%CoordType  /= H_TMCartesian .or. &
        HCoordOut%CoordType /= H_LatLong     .or. &
        HCoordOut%Pole(1)   /= 0.0           .or. &
        HCoordOut%Pole(2)   /= 0.0           .or. &
        HCoordOut%Angle     /= 0.0           .or. &
        HCoordOut%Origin(1) /= 0.0           .or. &
        HCoordOut%Origin(2) /= 0.0           .or. &
        HCoordOut%Unit(1)   /= 1.0           .or. &
        HCoordOut%Unit(2)   /= 1.0) Then
      Call Message('Error in TMCartesianToSLatLong', 4)
    End If
# endif

  ! Conversion constant for calculating normalised grid distances.
  DistanceCon = HCoordIn%ScaleFactor * EarthRadius

  ! Calculate normalised distances from true origin.
  X = (PointIn(1)*HCoordIn%Unit(1) + HCoordIn%Origin(1)) / DistanceCon
  Y = (PointIn(2)*HCoordIn%Unit(2) + HCoordIn%Origin(2)) / DistanceCon

  ! Pseudo-latitude and pseudo-longitude (relative to rotated pole).
  PLat  = 2.0 * ATan(Exp(X)) - Pi/2.0
  PLong = - (Y + HCoordIn%Pole(2))

  ! Calculate the true latitude and longitude.
  S           = - Cos(PLat) * Sin(PLong)
  C           = (Cos(PLat) * Cos(PLong))**2 + Sin(PLat)**2
  C           = Sqrt(C)
  PointOut(1) = ATan2ZeroTest(Sin(PLat), Cos(PLat) * Cos(PLong)) + HCoordIn%Pole(1)
  PointOut(2) = ATan2(S, C)

  If (PointOut(1) >   Pi) PointOut(1) = PointOut(1) - 2.0*Pi
  If (PointOut(1) < - Pi) PointOut(1) = PointOut(1) + 2.0*Pi

  TMCartesianToSLatLong = PointOut

End Function TMCartesianToSLatLong

!-------------------------------------------------------------------------------------------------------------

Function QuickTMCartesianToTMCartesian(HCoordIn, HCoordOut, PointIn)
! Converts between two cartesian coordinate systems (Easting/Northing grids) in a
! transverse Mercator projection, where the projection has the same value of Pole
! (true origin) and ScaleFactor.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: QuickTMCartesianToTMCartesian(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: Point(2)    ! Coords relative to true origin (in metres).

# ifdef ExtraChecks
    If (HCoordIn%CoordType   /= H_TMCartesian        .or. &
        HCoordOut%CoordType  /= H_TMCartesian        .or. &
        HCoordIn%Pole(1)     /= HCoordOut%Pole(1)    .or. &
        HCoordIn%Pole(2)     /= HCoordOut%Pole(2)    .or. &
        HCoordIn%ScaleFactor /= HCoordOut%ScaleFactor) Then
      Call Message('Error in QuickTMCartesianToTMCartesian', 4)
    End If
# endif

  ! If the two cartesian systems have the same origin, it is only necessary
  ! to scale the coordinates; otherwise the origin offset is also required.
  If (HCoordIn%Origin(1) == HCoordOut%Origin(1) .and. &
      HCoordIn%Origin(2) == HCoordOut%Origin(2)) Then
    PointOut(1) = PointIn(1) * HCoordIn%Unit(1) / HCoordOut%Unit(1)
    PointOut(2) = PointIn(2) * HCoordIn%Unit(2) / HCoordOut%Unit(2)
  Else
    ! Scaling of input coords.
    Point(1) = PointIn(1) * HCoordIn%Unit(1) + HCoordIn%Origin(1)
    Point(2) = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%Origin(2)
    ! Scaling of output coords.
    PointOut(1) = (Point(1) - HCoordOut%Origin(1)) / HCoordOut%Unit(1)
    PointOut(2) = (Point(2) - HCoordOut%Origin(2)) / HCoordOut%Unit(2)
  End If

  QuickTMCartesianToTMCartesian = PointOut

End Function QuickTMCartesianToTMCartesian

!-------------------------------------------------------------------------------------------------------------

Function TMCartesianToTMPolar(HCoordIn, HCoordOut, PointIn)
! Converts from a cartesian coord system to a polar system in a transverse Mercator
! projection where the coord systems have the same value for Pole and ScaleFactor.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: TMCartesianToTMPolar(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: Point(2)    ! Intermediate coords.

# ifdef ExtraChecks
    If (HCoordOut%CoordType  /= H_TMPolar                   .or. &
        HCoordIn%CoordType   /= H_TMCartesian               .or. &
        HCoordIn%Pole(1)     /= HCoordOut%Pole(1)           .or. &
        HCoordIn%Pole(2)     /= HCoordOut%Pole(2)           .or. &
        HCoordIn%ScaleFactor /= HCoordOut%ScaleFactor) Then
      Call Message('Error in TMCartesianToTMPolar', 4)
    End If
# endif

  ! Scaling of input coords.
  Point(1) = PointIn(1) * HCoordIn%Unit(1) + HCoordIn%Origin(1)
  Point(2) = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%Origin(2)

  ! Conversion to polar coords relative to polar origin and theta offset.
  Point(1)    = Point(1) - HCoordOut%Origin(1)
  Point(2)    = Point(2) - HCoordOut%Origin(2)
  PointOut(1) = Sqrt(Point(1)**2 + Point(2)**2)
  PointOut(2) = ATan2ZeroTest(Point(2), Point(1)) - HCoordOut%ThetaOrigin
  If (PointOut(2) >   Pi) PointOut(2) = PointOut(2) - 2.0*Pi
  If (PointOut(2) < - Pi) PointOut(2) = PointOut(2) + 2.0*Pi

  ! Scaling of output coords.
  PointOut(1) = PointOut(1) / HCoordOut%Unit(1)
  PointOut(2) = PointOut(2) / HCoordOut%Unit(2)

  TMCartesianToTMPolar = PointOut

End Function TMCartesianToTMPolar

!-------------------------------------------------------------------------------------------------------------

Function TMPolarToTMCartesian(HCoordIn, HCoordOut, PointIn)
! Converts from a polar coord system to a cartesian system in a transverse Mercator
! projection where the coord systems have the same value for Pole and ScaleFactor.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: TMPolarToTMCartesian(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: R           ! } Scaled input coords.
  Real(Std) :: Theta       ! }
  Real(Std) :: Point(2)    ! Coords relative to true origin (in metres).

# ifdef ExtraChecks
    If (HCoordOut%CoordType  /= H_TMCartesian               .or. &
        HCoordIn%CoordType   /= H_TMPolar                   .or. &
        HCoordIn%Pole(1)     /= HCoordOut%Pole(1)           .or. &
        HCoordIn%Pole(2)     /= HCoordOut%Pole(2)           .or. &
        HCoordIn%ScaleFactor /= HCoordOut%ScaleFactor) Then
      Call Message('Error in TMPolarToTMCartesian', 4)
    End If
# endif

  ! Scaling of input coords.
  R     = PointIn(1) * HCoordIn%Unit(1)
  Theta = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%ThetaOrigin

  ! Conversion to cartesian coords relative to true origin.
  Point(1) = R * Cos(Theta) + HCoordIn%Origin(1)
  Point(2) = R * Sin(Theta) + HCoordIn%Origin(2)

  ! Scaling of output coords.
  PointOut(1) = (Point(1) - HCoordOut%Origin(1)) / HCoordOut%Unit(1)
  PointOut(2) = (Point(2) - HCoordOut%Origin(2)) / HCoordOut%Unit(2)

  TMPolarToTMCartesian = PointOut

End Function TMPolarToTMCartesian

!-------------------------------------------------------------------------------------------------------------

Function QuickTMPolarToTMPolar(HCoordIn, HCoordOut, PointIn)
! Converts between two polar coord systems in a transverse Mercator projection
! where the coord systems have the same value for Pole, ScaleFactor and Origin.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn    ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut   ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2)  ! Input coords.
  ! Function result:
  Real(Pos) :: QuickTMPolarToTMPolar(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: Point(2)    ! Scaled coords in standard polar units.

# ifdef ExtraChecks
    If (HCoordIn%CoordType   /= H_TMPolar              .or. &
        HCoordOut%CoordType  /= H_TMPolar              .or. &
        HCoordIn%Pole(1)     /= HCoordOut%Pole(1)      .or. &
        HCoordIn%Pole(2)     /= HCoordOut%Pole(2)      .or. &
        HCoordIn%ScaleFactor /= HCoordOut%ScaleFactor  .or. &
        HCoordIn%Origin(1)   /= HCoordOut%Origin(1)    .or. &
        HCoordIn%Origin(2)   /= HCoordOut%Origin(2)) Then
      Call Message('Error in QuickTMPolarToTMPolar', 4)
    End If
# endif

  ! Scaling of input coords.
  Point(1) = PointIn(1) * HCoordIn%Unit(1)
  Point(2) = PointIn(2) * HCoordIn%Unit(2)

  ! Rotation of theta offset between coord systems.
  Point(2) = Point(2) + HCoordIn%ThetaOrigin - HCoordOut%ThetaOrigin
  If (Point(2) >   Pi) Point(2) = Point(2) - 2.0*Pi
  If (Point(2) < - Pi) Point(2) = Point(2) + 2.0*Pi

  ! Scaling of output coords.
  PointOut(1) = Point(1) / HCoordOut%Unit(1)
  PointOut(2) = Point(2) / HCoordOut%Unit(2)

  QuickTMPolarToTMPolar = PointOut

End Function QuickTMPolarToTMPolar

!-------------------------------------------------------------------------------------------------------------

Function QuickTMCartesianToTMPolar(HCoordIn, HCoordOut, PointIn)
! Converts from a cartesian coord system to a polar system in a transverse Mercator
! projection, where the systems have the same value for Pole, ScaleFactor and Origin.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn    ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut   ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2)  ! Input coords.
  ! Function result:
  Real(Pos) :: QuickTMCartesianToTMPolar(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: Point(2)    ! Scaled input coords in metres.
  Real(Std) :: R           ! Radial distance from origin in metres.
  Real(Std) :: Theta       ! Direction from origin in radians.

# ifdef ExtraChecks
    If (HCoordIn%CoordType   /= H_TMCartesian          .or. &
        HCoordOut%CoordType  /= H_TMPolar              .or. &
        HCoordIn%Pole(1)     /= HCoordOut%Pole(1)      .or. &
        HCoordIn%Pole(2)     /= HCoordOut%Pole(2)      .or. &
        HCoordIn%ScaleFactor /= HCoordOut%ScaleFactor  .or. &
        HCoordIn%Origin(1)   /= HCoordOut%Origin(1)    .or. &
        HCoordIn%Origin(2)   /= HCoordOut%Origin(2)) Then
      Call Message('Error in QuickTMCartesianToTMPolar', 4)
    End If
# endif

  ! Scaling of input coords.
  Point(1) = PointIn(1) * HCoordIn%Unit(1)
  Point(2) = PointIn(2) * HCoordIn%Unit(2)

  ! Conversion to polar coords centred at origin.
  R     = Sqrt(Point(1)**2 + Point(2)**2)
  Theta = ATan2ZeroTest(Point(2), Point(1))

  ! Rotation to theta offset of output coord system.
  Theta = Theta - HCoordOut%ThetaOrigin
  If (Theta >   Pi) Theta = Theta - 2.0*Pi
  If (Theta < - Pi) Theta = Theta + 2.0*Pi

  ! Scaling of output coords.
  PointOut(1) = R     / HCoordOut%Unit(1)
  PointOut(2) = Theta / HCoordOut%Unit(2)

  QuickTMCartesianToTMPolar = PointOut

End Function QuickTMCartesianToTMPolar

!-------------------------------------------------------------------------------------------------------------

Function QuickTMPolarToTMCartesian(HCoordIn, HCoordOut, PointIn)
! Converts from a polar coord system to a cartesian system in a transverse Mercator
! projection, where the systems have the same value for Pole, ScaleFactor and Origin.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn    ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut   ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2)  ! Input coords.
  ! Function result:
  Real(Pos) :: QuickTMPolarToTMCartesian(2) ! Output coords.
  ! Locals:
  Real(Pos) :: PointOut(2) ! Local copy of output coords.
  Real(Std) :: R           ! Radial distance from origin in metres.
  Real(Std) :: Theta       ! Direction from origin in radians (relative to East).
  Real(Std) :: Point(2)    ! Cartesian coords centred at origin in metres.

# ifdef ExtraChecks
    If (HCoordIn%CoordType   /= H_TMPolar              .or. &
        HCoordOut%CoordType  /= H_TMCartesian          .or. &
        HCoordIn%Pole(1)     /= HCoordOut%Pole(1)      .or. &
        HCoordIn%Pole(2)     /= HCoordOut%Pole(2)      .or. &
        HCoordIn%ScaleFactor /= HCoordOut%ScaleFactor  .or. &
        HCoordIn%Origin(1)   /= HCoordOut%Origin(1)    .or. &
        HCoordIn%Origin(2)   /= HCoordOut%Origin(2)) Then
      Call Message('Error in QuickTMPolarToTMCartesian', 4)
    End If
# endif

  ! Scaling of input coords.
  R     = PointIn(1) * HCoordIn%Unit(1)
  Theta = PointIn(2) * HCoordIn%Unit(2) + HCoordIn%ThetaOrigin

  ! Conversion to cartesian coords centred at origin.
  Point(1) = R * Cos(Theta)
  Point(2) = R * Sin(Theta)

  ! Scaling of output coords.
  PointOut(1) = Point(1) / HCoordOut%Unit(1)
  PointOut(2) = Point(2) / HCoordOut%Unit(2)

  QuickTMPolarToTMCartesian = PointOut

End Function QuickTMPolarToTMCartesian

!-------------------------------------------------------------------------------------------------------------

Function ConvertH(HCoordIn, HCoordOut, PointIn)
! Converts horizontal coords between horizontal coord systems.
!
! Note: the design philosophy adopted with intermediate coord systems is:
!        + to have no origin offset (i.e. Origin = (0, 0)) and
!          to use the standard units (i.e. Units = (1, 1))
!          for all intermediate systems, *except*
!        + when converting between TMCartesian <-> TMPolar using the
!          'Quick' conversion routines (where origin offsets are required,
!          but the above constraint is still used for scaling the units).

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In) :: HCoordIn   ! Coord system of input coords.
  Type(HCoord_), Intent(In) :: HCoordOut  ! Coord system of output coords.
  Real(Pos),     Intent(In) :: PointIn(2) ! Input coords.
  ! Function result:
  Real(Pos) :: ConvertH(2) ! Output coords.
  ! Locals:
  Real(Pos)     :: PointOut(2) ! Local copy of output coords.
  Type(HCoord_) :: HCoord1     !} Intermediate coord systems.
  Type(HCoord_) :: HCoord2     !}
  Type(HCoord_) :: HCoord3     !}
  Type(HCoord_) :: HCoord4     !}
  Real(Pos)     :: Point1(2)   !] Coords in intermediate coord systems.
  Real(Pos)     :: Point2(2)   !]
  Real(Pos)     :: Point3(2)   !]
  Real(Pos)     :: Point4(2)   !]

  Select Case (HCoordIn%CoordType)

    Case (H_LatLong)

      Select Case (HCoordOut%CoordType)

        Case (H_LatLong)      ! IN: H_LatLong      OUT: H_LatLong

          If (HCoordIn%Pole(1) == HCoordOut%Pole(1) .and. &
              HCoordIn%Pole(2) == HCoordOut%Pole(2)) Then
            PointOut = QuickLatLongToLatLong(HCoordIn, HCoordOut, PointIn)
          Else
            PointOut = LatLongToLatLong(HCoordIn, HCoordOut, PointIn)
          End If

        Case (H_PSCartesian)  ! IN: H_LatLong      OUT: H_PSCartesian

          If (HCoordIn%Pole(1) == HCoordOut%Pole(1) .and. &
              HCoordIn%Pole(2) == HCoordOut%Pole(2)) Then
            PointOut = LatLongToPSCartesian(HCoordIn, HCoordOut, PointIn)
          Else
            HCoord1 = HCoordOut ; HCoord1%CoordType = H_LatLong
            HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
            Point1   = LatLongToLatLong(HCoordIn, HCoord1, PointIn)
            PointOut = LatLongToPSCartesian(HCoord1, HCoordOut, Point1)
          End If

        Case (H_PSPolar)      ! IN: H_LatLong      OUT: H_PSPolar

          If (HCoordIn%Pole(1) == HCoordOut%Pole(1) .and. &
              HCoordIn%Pole(2) == HCoordOut%Pole(2)) Then
            HCoord1 = HCoordOut ; HCoord1%CoordType = H_PSCartesian
            HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
            Point1   = LatLongToPSCartesian(HCoordIn, HCoord1, PointIn)
            PointOut = PSCartesianToPSPolar(HCoord1, HCoordOut, Point1)
          Else
            HCoord1 = HCoordOut ; HCoord1%CoordType = H_LatLong
            HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
            HCoord2 = HCoordOut ; HCoord2%CoordType = H_PSCartesian
            HCoord2%Origin = (/ 0.0, 0.0 /) ; HCoord2%Unit = (/ 1.0, 1.0 /)
            Point1   = LatLongToLatLong(HCoordIn, HCoord1, PointIn)
            Point2   = LatLongToPSCartesian(HCoord1, HCoord2, Point1)
            PointOut = PSCartesianToPSPolar(HCoord2, HCoordOut, Point2)
          End If

        Case (H_TMCartesian)  ! IN: H_LatLong      OUT: H_TMCartesian

          HCoord1 = HCoord_LatLongRadians()
          If (HCoordIn .equiv. HCoord1) Then
            PointOut = SLatLongToTMCartesian(HCoordIn, HCoordOut, PointIn)
          Else If (HCoordIn%Pole(1) == 0.0  .and. &
                   HCoordIn%Pole(2) == 0.0) Then
            Point1   = QuickLatLongToLatLong(HCoordIn, HCoord1, PointIn)
            PointOut = SLatLongToTMCartesian(HCoord1, HCoordOut, Point1)
          Else
            Point1   = LatLongToSLatLong(HCoordIn, HCoord1, PointIn)
            PointOut = SLatLongToTMCartesian(HCoord1, HCoordOut, Point1)
          End If

        Case (H_TMPolar)      ! IN: H_LatLong      OUT: H_TMPolar

          HCoord1 = HCoord_LatLongRadians()
          HCoord2 = HCoordOut ; HCoord2%CoordType = H_TMCartesian
          HCoord2%Unit = (/ 1.0, 1.0 /)
          If (HCoordIn .equiv. HCoord1) Then
            Point2   = SLatLongToTMCartesian(HCoordIn, HCoord2, PointIn)
            PointOut = QuickTMCartesianToTMPolar(HCoord2, HCoordOut, Point2)
          Else If (HCoordIn%Pole(1) == 0.0  .and. &
                   HCoordIn%Pole(2) == 0.0) Then
            Point1   = QuickLatLongToLatLong(HCoordIn, HCoord1, PointIn)
            Point2   = SLatLongToTMCartesian(HCoord1, HCoord2, Point1)
            PointOut = QuickTMCartesianToTMPolar(HCoord2, HCoordOut, Point2)
          Else
            Point1   = LatLongToSLatLong(HCoordIn, HCoord1, PointIn)
            Point2   = SLatLongToTMCartesian(HCoord1, HCoord2, Point1)
            PointOut = QuickTMCartesianToTMPolar(HCoord2, HCoordOut, Point2)
          End If

        Case Default

          Call Message('Error in ConvertH', 4)

      End Select

    Case (H_PSCartesian)

      Select Case (HCoordOut%CoordType)

        Case (H_LatLong)      ! IN: H_PSCartesian  OUT: H_LatLong

          If (HCoordIn%Pole(1) == HCoordOut%Pole(1) .and. &
              HCoordIn%Pole(2) == HCoordOut%Pole(2)) Then
            PointOut = PSCartesianToLatLong(HCoordIn, HCoordOut, PointIn)
          Else
            HCoord1 = HCoordIn ; HCoord1%CoordType = H_LatLong
            HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
            Point1   = PSCartesianToLatLong(HCoordIn, HCoord1, PointIn)
            PointOut = LatLongToLatLong(HCoord1, HCoordOut, Point1)
          End If

        Case (H_PSCartesian)  ! IN: H_PSCartesian  OUT: H_PSCartesian

          If (HCoordIn%Pole(1) == HCoordOut%Pole(1) .and. &
              HCoordIn%Pole(2) == HCoordOut%Pole(2)) Then
            PointOut = QuickPSCartesianToPSCartesian(HCoordIn, HCoordOut, PointIn)
          Else
            HCoord1 = HCoordIn  ; HCoord1%CoordType = H_LatLong
            HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
            HCoord2 = HCoordOut ; HCoord2%CoordType = H_LatLong
            HCoord2%Origin = (/ 0.0, 0.0 /) ; HCoord2%Unit = (/ 1.0, 1.0 /)
            Point1   = PSCartesianToLatLong(HCoordIn, HCoord1, PointIn)
            Point2   = LatLongToLatLong(HCoord1, HCoord2, Point1)
            PointOut = LatLongToPSCartesian(HCoord2, HCoordOut, Point2)
          End If

        Case (H_PSPolar)      ! IN: H_PSCartesian  OUT: H_PSPolar

          If (HCoordIn%Pole(1) == HCoordOut%Pole(1) .and. &
              HCoordIn%Pole(2) == HCoordOut%Pole(2)) Then
            PointOut = PSCartesianToPSPolar(HCoordIn, HCoordOut, PointIn)
          Else
            HCoord1 = HCoordIn  ; HCoord1%CoordType = H_LatLong
            HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
            HCoord2 = HCoordOut ; HCoord2%CoordType = H_LatLong
            HCoord2%Origin = (/ 0.0, 0.0 /) ; HCoord2%Unit = (/ 1.0, 1.0 /)
            HCoord3 = HCoordOut ; HCoord3%CoordType = H_PSCartesian
            HCoord3%Origin = (/ 0.0, 0.0 /) ; HCoord3%Unit = (/ 1.0, 1.0 /)
            Point1   = PSCartesianToLatLong(HCoordIn, HCoord1, PointIn)
            Point2   = LatLongToLatLong(HCoord1, HCoord2, Point1)
            Point3   = LatLongToPSCartesian(HCoord2, HCoord3, Point2)
            PointOut = PSCartesianToPSPolar(HCoord3, HCoordOut, Point3)
          End If

        Case (H_TMCartesian)  ! IN: H_PSCartesian  OUT: H_TMCartesian

          HCoord2 = HCoord_LatLongRadians()
          If (HCoordIn%Pole(1) == 0.0  .and. &
              HCoordIn%Pole(2) == 0.0) Then
            Point2   = PSCartesianToLatLong(HCoordIn, HCoord2, PointIn)
            PointOut = SLatLongToTMCartesian(HCoord2, HCoordOut, Point2)
          Else
            HCoord1 = HCoordIn ; HCoord1%CoordType = H_LatLong
            HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
            Point1   = PSCartesianToLatLong(HCoordIn, HCoord1, PointIn)
            Point2   = LatLongToSLatLong(HCoord1, HCoord2, Point1)
            PointOut = SLatLongToTMCartesian(HCoord2, HCoordOut, Point2)
          End If

        Case (H_TMPolar)      ! IN: H_PSCartesian  OUT: H_TMPolar

          HCoord2 = HCoord_LatLongRadians()
          HCoord3 = HCoordOut ; HCoord3%CoordType = H_TMCartesian
          HCoord3%Unit = (/ 1.0, 1.0 /)
          If (HCoordIn%Pole(1) == 0.0  .and. &
              HCoordIn%Pole(2) == 0.0) Then
            Point2   = PSCartesianToLatLong(HCoordIn, HCoord2, PointIn)
            Point3   = SLatLongToTMCartesian(HCoord2, HCoord3, Point2)
            PointOut = QuickTMCartesianToTMPolar(HCoord3, HCoordOut, Point3)
          Else
            HCoord1 = HCoordIn ; HCoord1%CoordType = H_LatLong
            HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
            Point1   = PSCartesianToLatLong(HCoordIn, HCoord1, PointIn)
            Point2   = LatLongToSLatLong(HCoord1, HCoord2, Point1)
            Point3   = SLatLongToTMCartesian(HCoord2, HCoord3, Point2)
            PointOut = QuickTMCartesianToTMPolar(HCoord3, HCoordOut, Point3)
          End If

        Case Default

          Call Message('Error in ConvertH', 4)

      End Select

    Case (H_PSPolar)

      Select Case (HCoordOut%CoordType)

        Case (H_LatLong)      ! IN: H_PSPolar      OUT: H_LatLong

          HCoord1 = HCoordIn ; HCoord1%CoordType = H_PSCartesian
          HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
          Point1 = PSPolarToPSCartesian(HCoordIn, HCoord1, PointIn)
          If (HCoordIn%Pole(1) == HCoordOut%Pole(1) .and. &
              HCoordIn%Pole(2) == HCoordOut%Pole(2)) Then
            PointOut = PSCartesianToLatLong(HCoord1, HCoordOut, Point1)
          Else
            HCoord2 = HCoordIn ; HCoord2%CoordType = H_LatLong
            HCoord2%Origin = (/ 0.0, 0.0 /) ; HCoord2%Unit = (/ 1.0, 1.0 /)
            Point2   = PSCartesianToLatLong(HCoord1, HCoord2, Point1)
            PointOut = LatLongToLatLong(HCoord2, HCoordOut, Point2)
          End If

        Case (H_PSCartesian)  ! IN: H_PSPolar      OUT: H_PSCartesian

          If (HCoordIn%Pole(1) == HCoordOut%Pole(1)  .and. &
              HCoordIn%Pole(2) == HCoordOut%Pole(2)) Then
            PointOut = PSPolarToPSCartesian(HCoordIn, HCoordOut, PointIn)
          Else
            HCoord1 = HCoordIn  ; HCoord1%CoordType = H_PSCartesian
            HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
            HCoord2 = HCoordIn  ; HCoord2%CoordType = H_LatLong
            HCoord2%Origin = (/ 0.0, 0.0 /) ; HCoord2%Unit = (/ 1.0, 1.0 /)
            HCoord3 = HCoordOut ; HCoord3%CoordType = H_LatLong
            HCoord3%Origin = (/ 0.0, 0.0 /) ; HCoord3%Unit = (/ 1.0, 1.0 /)
            Point1   = PSPolarToPSCartesian(HCoordIn, HCoord1, PointIn)
            Point2   = PSCartesianToLatLong(HCoord1, HCoord2, Point1)
            Point3   = LatLongToLatLong(HCoord2, HCoord3, Point2)
            PointOut = LatLongToPSCartesian(HCoord3, HCoordOut, Point3)
          End If

        Case (H_PSPolar)      ! IN: H_PSPolar      OUT: H_PSPolar

          If (HCoordIn%Pole(1) == HCoordOut%Pole(1)  .and. &
              HCoordIn%Pole(2) == HCoordOut%Pole(2)) Then
            If (HCoordIn%Angle     == HCoordOut%Angle      .and. &
                HCoordIn%Origin(1) == HCoordOut%Origin(1)  .and. &
                HCoordIn%Origin(2) == HCoordOut%Origin(2)) Then
              PointOut = QuickPSPolarToPSPolar(HCoordIn, HCoordOut, PointIn)
            Else
              HCoord1 = HCoordIn ; HCoord1%CoordType = H_PSCartesian
              HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
              Point1   = PSPolarToPSCartesian(HCoordIn, HCoord1, PointIn)
              PointOut = PSCartesianToPSPolar(HCoord1, HCoordOut, Point1)
            End If
          Else
            HCoord1 = HCoordIn  ; HCoord1%CoordType = H_PSCartesian
            HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
            HCoord2 = HCoordIn  ; HCoord2%CoordType = H_LatLong
            HCoord2%Origin = (/ 0.0, 0.0 /) ; HCoord2%Unit = (/ 1.0, 1.0 /)
            HCoord3 = HCoordOut ; HCoord3%CoordType = H_LatLong
            HCoord3%Origin = (/ 0.0, 0.0 /) ; HCoord3%Unit = (/ 1.0, 1.0 /)
            HCoord4 = HCoordOut ; HCoord4%CoordType = H_PSCartesian
            HCoord4%Origin = (/ 0.0, 0.0 /) ; HCoord4%Unit = (/ 1.0, 1.0 /)
            Point1   = PSPolarToPSCartesian(HCoordIn, HCoord1, PointIn)
            Point2   = PSCartesianToLatLong(HCoord1, HCoord2, Point1)
            Point3   = LatLongToLatLong(HCoord2, HCoord3, Point2)
            Point4   = LatLongToPSCartesian(HCoord3, HCoord4, Point3)
            PointOut = PSCartesianToPSPolar(HCoord4, HCoordOut, Point4)
          End If

        Case (H_TMCartesian)  ! IN: H_PSPolar      OUT: H_TMCartesian

          HCoord1 = HCoordIn ; HCoord1%CoordType = H_PSCartesian
          HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
          Point1 = PSPolarToPSCartesian(HCoordIn, HCoord1, PointIn)
          HCoord3 = HCoord_LatLongRadians()
          If (HCoordIn%Pole(1) == 0.0  .and. &
              HCoordIn%Pole(2) == 0.0) Then
            Point3   = PSCartesianToLatLong(HCoord1, HCoord3, Point1)
            PointOut = SLatLongToTMCartesian(HCoord3, HCoordOut, Point3)
          Else
            HCoord2 = HCoordIn ; HCoord2%CoordType = H_LatLong
            HCoord2%Origin = (/ 0.0, 0.0 /) ; HCoord2%Unit = (/ 1.0, 1.0 /)
            Point2   = PSCartesianToLatLong(HCoord1, HCoord2, Point1)
            Point3   = LatLongToSLatLong(HCoord2, HCoord3, Point2)
            PointOut = SLatLongToTMCartesian(HCoord3, HCoordOut, Point3)
          End If

        Case (H_TMPolar)      ! IN: H_PSPolar      OUT: H_TMPolar

          HCoord1 = HCoordIn ; HCoord1%CoordType = H_PSCartesian
          HCoord1%Origin = (/ 0.0, 0.0 /) ; HCoord1%Unit = (/ 1.0, 1.0 /)
          Point1 = PSPolarToPSCartesian(HCoordIn, HCoord1, PointIn)
          HCoord3 = HCoord_LatLongRadians()
          HCoord4 = HCoordOut ; HCoord4%CoordType = H_TMCartesian
          HCoord4%Unit = (/ 1.0, 1.0 /)
          If (HCoordIn%Pole(1) == 0.0  .and. &
              HCoordIn%Pole(2) == 0.0) Then
            Point3   = PSCartesianToLatLong(HCoord1, HCoord3, Point1)
            Point4   = SLatLongToTMCartesian(HCoord3, HCoord4, Point3)
            PointOut = QuickTMCartesianToTMPolar(HCoord4, HCoordOut, Point4)
          Else
            HCoord2 = HCoordIn ; HCoord2%CoordType = H_LatLong
            HCoord2%Origin = (/ 0.0, 0.0 /) ; HCoord2%Unit = (/ 1.0, 1.0 /)
            Point2   = PSCartesianToLatLong(HCoord1, HCoord2, Point1)
            Point3   = LatLongToSLatLong(HCoord2, HCoord3, Point2)
            Point4   = SLatLongToTMCartesian(HCoord3, HCoord4, Point3)
            PointOut = QuickTMCartesianToTMPolar(HCoord4, HCoordOut, Point4)
          End If

        Case Default

          Call Message('Error in ConvertH', 4)

      End Select

    Case (H_TMCartesian)

      Select Case (HCoordOut%CoordType)

        Case (H_LatLong)      ! IN: H_TMCartesian  OUT: H_LatLong

          HCoord1 = HCoord_LatLongRadians()
          If (HCoordOut .equiv. HCoord1) Then
            PointOut = TMCartesianToSLatLong(HCoordIn, HCoordOut, PointIn)
          Else If (HCoordOut%Pole(1) == 0.0  .and. &
                   HCoordOut%Pole(2) == 0.0) Then
            Point1   = TMCartesianToSLatLong(HCoordIn, HCoord1, PointIn)
            PointOut = QuickLatLongToLatLong(HCoord1, HCoordOut, Point1)
          Else
            Point1   = TMCartesianToSLatLong(HCoordIn, HCoord1, PointIn)
            PointOut = SLatLongToLatLong(HCoord1, HCoordOut, Point1)
          End If

        Case (H_PSCartesian)  ! IN: H_TMCartesian  OUT: H_PSCartesian

          HCoord1 = HCoord_LatLongRadians()
          Point1  = TMCartesianToSLatLong(HCoordIn, HCoord1, PointIn)
          If (HCoordOut%Pole(1) == 0.0  .and. &
              HCoordOut%Pole(2) == 0.0) Then
            PointOut = LatLongToPSCartesian(HCoord1, HCoordOut, Point1)
          Else
            HCoord2 = HCoordOut ; HCoord2%CoordType = H_LatLong
            HCoord2%Origin = (/ 0.0, 0.0 /) ; HCoord2%Unit = (/ 1.0, 1.0 /)
            Point2   = SLatLongToLatLong(HCoord1, HCoord2, Point1)
            PointOut = LatLongToPSCartesian(HCoord2, HCoordOut, Point2)
          End If

        Case (H_PSPolar)      ! IN: H_TMCartesian  OUT: H_PSPolar

          HCoord1 = HCoord_LatLongRadians()
          Point1  = TMCartesianToSLatLong(HCoordIn, HCoord1, PointIn)
          If (HCoordOut%Pole(1) == 0.0  .and. &
              HCoordOut%Pole(2) == 0.0) Then
            HCoord2 = HCoordOut ; HCoord2%CoordType = H_PSCartesian
            HCoord2%Origin = (/ 0.0, 0.0 /) ; HCoord2%Unit = (/ 1.0, 1.0 /)
            Point2   = LatLongToPSCartesian(HCoord1, HCoord2, Point1)
            PointOut = PSCartesianToPSPolar(HCoord2, HCoordOut, Point2)
          Else
            HCoord2 = HCoordOut ; HCoord2%CoordType = H_LatLong
            HCoord2%Origin = (/ 0.0, 0.0 /) ; HCoord2%Unit = (/ 1.0, 1.0 /)
            HCoord3 = HCoordOut ; HCoord3%CoordType = H_PSCartesian
            HCoord3%Origin = (/ 0.0, 0.0 /) ; HCoord3%Unit = (/ 1.0, 1.0 /)
            Point2   = SLatLongToLatLong(HCoord1, HCoord2, Point1)
            Point3   = LatLongToPSCartesian(HCoord2, HCoord3, Point2)
            PointOut = PSCartesianToPSPolar(HCoord3, HCoordOut, Point3)
          End If

        Case (H_TMCartesian)  ! IN: H_TMCartesian  OUT: H_TMCartesian

          If (HCoordIn%Pole(1)     == HCoordOut%Pole(1)      .and. &
              HCoordIn%Pole(2)     == HCoordOut%Pole(2)      .and. &
              HCoordIn%ScaleFactor == HCoordOut%ScaleFactor) Then
            PointOut = QuickTMCartesianToTMCartesian(HCoordIn, HCoordOut, PointIn)
          Else
            HCoord1  = HCoord_LatLongRadians()
            Point1   = TMCartesianToSLatLong(HCoordIn, HCoord1, PointIn)
            PointOut = SLatLongToTMCartesian(HCoord1, HCoordOut, Point1)
          End If

        Case (H_TMPolar)      ! IN: H_TMCartesian  OUT: H_TMPolar

          If (HCoordIn%Pole(1)     == HCoordOut%Pole(1)      .and. &
              HCoordIn%Pole(2)     == HCoordOut%Pole(2)      .and. &
              HCoordIn%ScaleFactor == HCoordOut%ScaleFactor) Then
            If (HCoordIn%Origin(1) == HCoordOut%Origin(1)  .and. &
                HCoordIn%Origin(2) == HCoordOut%Origin(2)) Then
              PointOut = QuickTMCartesianToTMPolar(HCoordIn, HCoordOut, PointIn)
            Else
              PointOut = TMCartesianToTMPolar(HCoordIn, HCoordOut, PointIn)
            End If
          Else
            HCoord1 = HCoord_LatLongRadians()
            HCoord2 = HCoordOut ; HCoord2%CoordType = H_TMCartesian
            HCoord2%Unit = (/ 1.0, 1.0 /)
            Point1   = TMCartesianToSLatLong(HCoordIn, HCoord1, PointIn)
            Point2   = SLatLongToTMCartesian(HCoord1, HCoord2, Point1)
            PointOut = QuickTMCartesianToTMPolar(HCoord2, HCoordOut, Point2)
          End If

        Case Default

          Call Message('Error in ConvertH', 4)

      End Select

    Case (H_TMPolar)

      Select Case (HCoordOut%CoordType)

        Case (H_LatLong)      ! IN: H_TMPolar      OUT: H_LatLong

          HCoord1 = HCoordIn ; HCoord1%CoordType = H_TMCartesian
          HCoord1%Unit = (/ 1.0, 1.0 /)
          Point1 = QuickTMPolarToTMCartesian(HCoordIn, HCoord1, PointIn)
          HCoord2 = HCoord_LatLongRadians()
          If (HCoordOut .equiv. HCoord2) Then
            PointOut = TMCartesianToSLatLong(HCoord1, HCoordOut, Point1)
          Else If (HCoordOut%Pole(1) == 0.0  .and. &
                   HCoordOut%Pole(2) == 0.0) Then
            Point2   = TMCartesianToSLatLong(HCoord1, HCoord2, Point1)
            PointOut = QuickLatLongToLatLong(HCoord2, HCoordOut, Point2)
          Else
            Point2   = TMCartesianToSLatLong(HCoord1, HCoord2, Point1)
            PointOut = SLatLongToLatLong(HCoord2, HCoordOut, Point2)
          End If

        Case (H_PSCartesian)  ! IN: H_TMPolar      OUT: H_PSCartesian

          HCoord1 = HCoordIn ; HCoord1%CoordType = H_TMCartesian
          HCoord1%Unit = (/ 1.0, 1.0 /)
          Point1 = QuickTMPolarToTMCartesian(HCoordIn, HCoord1, PointIn)
          HCoord2 = HCoord_LatLongRadians()
          Point2  = TMCartesianToSLatLong(HCoord1, HCoord2, Point1)
          If (HCoordOut%Pole(1) == 0.0  .and. &
              HCoordOut%Pole(2) == 0.0) Then
            PointOut = LatLongToPSCartesian(HCoord2, HCoordOut, Point2)
          Else
            HCoord3 = HCoordOut ; HCoord3%CoordType = H_LatLong
            HCoord3%Origin = (/ 0.0, 0.0 /) ; HCoord3%Unit = (/ 1.0, 1.0 /)
            Point3   = SLatLongToLatLong(HCoord2, HCoord3, Point2)
            PointOut = LatLongToPSCartesian(HCoord3, HCoordOut, Point3)
          End If

        Case (H_PSPolar)      ! IN: H_TMPolar      OUT: H_PSPolar

          HCoord1 = HCoordIn ; HCoord1%CoordType = H_TMCartesian
          HCoord1%Unit = (/ 1.0, 1.0 /)
          Point1 = QuickTMPolarToTMCartesian(HCoordIn, HCoord1, PointIn)
          HCoord2 = HCoord_LatLongRadians()
          Point2  = TMCartesianToSLatLong(HCoord1, HCoord2, Point1)
          If (HCoordOut%Pole(1) == 0.0  .and. &
              HCoordOut%Pole(2) == 0.0) Then
            HCoord3 = HCoordOut ; HCoord3%CoordType = H_PSCartesian
            HCoord3%Origin = (/ 0.0, 0.0 /) ; HCoord3%Unit = (/ 1.0, 1.0 /)
            Point3   = LatLongToPSCartesian(HCoord2, HCoord3, Point2)
            PointOut = PSCartesianToPSPolar(HCoord3, HCoordOut, Point3)
          Else
            HCoord3 = HCoordOut ; HCoord3%CoordType = H_LatLong
            HCoord3%Origin = (/ 0.0, 0.0 /) ; HCoord3%Unit = (/ 1.0, 1.0 /)
            HCoord4 = HCoordOut ; HCoord4%CoordType = H_PSCartesian
            HCoord4%Origin = (/ 0.0, 0.0 /) ; HCoord4%Unit = (/ 1.0, 1.0 /)
            Point3   = SLatLongToLatLong(HCoord2, HCoord3, Point2)
            Point4   = LatLongToPSCartesian(HCoord3, HCoord4, Point3)
            PointOut = PSCartesianToPSPolar(HCoord4, HCoordOut, Point4)
          End If

        Case (H_TMCartesian)  ! IN: H_TMPolar      OUT: H_TMCartesian

          If (HCoordIn%Pole(1)     == HCoordOut%Pole(1)      .and. &
              HCoordIn%Pole(2)     == HCoordOut%Pole(2)      .and. &
              HCoordIn%ScaleFactor == HCoordOut%ScaleFactor) Then
            If (HCoordIn%Origin(1) == HCoordOut%Origin(1)  .and. &
                HCoordIn%Origin(2) == HCoordOut%Origin(2)) Then
              PointOut = QuickTMPolarToTMCartesian(HCoordIn, HCoordOut, PointIn)
            Else
              PointOut = TMPolarToTMCartesian(HCoordIn, HCoordOut, PointIn)
            End If
          Else
            HCoord1 = HCoordIn ; HCoord1%CoordType = H_TMCartesian
            HCoord1%Unit = (/ 1.0, 1.0 /)
            HCoord2  = HCoord_LatLongRadians()
            Point1   = QuickTMPolarToTMCartesian(HCoordIn, HCoord1, PointIn)
            Point2   = TMCartesianToSLatLong(HCoord1, HCoord2, Point1)
            PointOut = SLatLongToTMCartesian(HCoord2, HCoordOut, Point2)
          End If

        Case (H_TMPolar)      ! IN: H_TMPolar      OUT: H_TMPolar

          If (HCoordIn%Pole(1)     == HCoordOut%Pole(1)      .and. &
              HCoordIn%Pole(2)     == HCoordOut%Pole(2)      .and. &
              HCoordIn%ScaleFactor == HCoordOut%ScaleFactor) Then
            If (HCoordIn%Origin(1) == HCoordOut%Origin(1)  .and. &
                HCoordIn%Origin(2) == HCoordOut%Origin(2)) Then
              PointOut = QuickTMPolarToTMPolar(HCoordIn, HCoordOut, PointIn)
            Else
              HCoord1 = HCoordOut ; HCoord1%CoordType = H_TMCartesian
              HCoord1%Unit = (/ 1.0, 1.0 /)
              Point1   = TMPolarToTMCartesian(HCoordIn, HCoord1, PointIn)
              PointOut = QuickTMCartesianToTMPolar(HCoord1, HCoordOut, Point1)
            End If
          Else
            HCoord1 = HCoordIn  ; HCoord1%CoordType = H_TMCartesian
            HCoord1%Unit = (/ 1.0, 1.0 /)
            HCoord2 = HCoord_LatLongRadians()
            HCoord3 = HCoordOut ; HCoord3%CoordType = H_TMCartesian
            HCoord3%Unit = (/ 1.0, 1.0 /)
            Point1   = QuickTMPolarToTMCartesian(HCoordIn, HCoord1, PointIn)
            Point2   = TMCartesianToSLatLong(HCoord1, HCoord2, Point1)
            Point3   = SLatLongToTMCartesian(HCoord2, HCoord3, Point2)
            PointOut = QuickTMCartesianToTMPolar(HCoord3, HCoordOut, Point3)
          End If

        Case Default

          Call Message('Error in ConvertH', 4)

      End Select

    Case Default

      Call Message('Error in ConvertH', 4)

  End Select

  ConvertH = PointOut

End Function ConvertH

!-------------------------------------------------------------------------------------------------------------

Subroutine MetricCoeffs(HCoord, Point, HMax, H1, H2)
! Calculates metric coefficients for lat-long and cartesian horizontal coord systems.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In)  :: HCoord   ! Coord system for which metric coefficients
                                         ! are required.
  Real(Pos),     Intent(In)  :: Point(2) ! Coords of location for which metric
                                         ! coefficients are required.
  Real(Std),     Intent(Out) :: HMax     ! Max value of H1 and H2 for this coord
                                         ! system and for any point.
  Real(Std),     Intent(Out) :: H1       ! Distance per unit change in x coord.
  Real(Std),     Intent(Out) :: H2       ! Distance per unit change in y coord.

  If (HCoord%CoordType == H_LatLong) Then

    H1 = EarthRadius * HCoord%Unit(1) *                  &
         Cos(Point(2)*HCoord%Unit(2) + HCoord%Origin(2))
    H2 = EarthRadius * HCoord%Unit(2)
    HMax = Max(EarthRadius * HCoord%Unit(1), EarthRadius * HCoord%Unit(2))

  Else If (HCoord%CoordType == H_PSCartesian) Then

    H1 = HCoord%Unit(1) ! $$
    H2 = HCoord%Unit(2) ! $$
    HMax = Max(H1, H2)  ! $$

  Else If (HCoord%CoordType == H_TMCartesian) Then

    H1 = HCoord%Unit(1) ! $$
    H2 = HCoord%Unit(2) ! $$
    HMax = Max(H1, H2)  ! $$

  Else

    Call Message('Error in MetricCoeffs', 4)

  End If

End Subroutine MetricCoeffs

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcDistanceBetweenTwoPoints( &
             HCoord,                     &
             Point1, Point2,             &
             Distance                    &
           )

! Calculates distance between two points for lat-long and cartesian horizontal coord systems.

! This subroutine assumes that the vertical coord system is 'm agl' and the horizontal coord system is one
! of the above types - for efficiency reasons no additional checks are carried out in this subroutine.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In)  :: HCoord    ! Horizontal coord system.
  Real(Std),     Intent(In)  :: Point1(3) ! Coords of first location.
  Real(Std),     Intent(In)  :: Point2(3) ! Coords of second location.
  Real(Std)                  :: Distance  ! Distance between the two locations (in metres).
  ! Locals:
  Real(Std)     :: H1       ! Distance in metres per unit change in x coord.
  Real(Std)     :: H2       ! Distance in metres per unit change in y coord.
  Real(Std)     :: HMax     ! Dummy value not used further in this subroutine.

  ! Calculate metric terms at second location
  Call MetricCoeffs(HCoord, Point2, HMax, H1, H2)

  ! Calculate distance between the two points
  Distance = (H1 * (Point2(1) - Point1(1)))**2 + (H2 * (Point2(2) - Point1(2)))**2  &
             + (Point2(3) - Point1(3))**2
  Distance = Sqrt(Distance)

End Subroutine CalcDistanceBetweenTwoPoints

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcHAngle(HCoord1, HCoord2, PointIn1, Angle)
! Calculates angle between two horizontal coord systems at a given position (for
! lat-long and cartesian type horizontal coord systems).

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In)  :: HCoord1     ! First horizontal coord system.
  Type(HCoord_), Intent(In)  :: HCoord2     ! Second horizontal coord system.
  Real(Std),     Intent(In)  :: PointIn1(2) ! Coords of position in HCoord1 system.
  Real(Std),     Intent(Out) :: Angle       ! Rotation angle of HCoord2's x-axis (for
                                            ! cartesian) or easterly (for lat-long)
                                            ! from HCoord1's x-axis or easterly
                                            ! (in radians, and positive for
                                            ! anticlockwise rotations).
  ! Locals:
  Real(Std), Parameter :: Small = 1.0E-15_Std ! If latitude difference from pole less
                                              ! than this, treat position as at pole.
  Type(HCoord_) :: HCoord    !  Standard lat-long coord system.
  Type(HCoord_) :: HCoord0   !  Rotated lat-long coord system (for TM projection).
  Type(HCoord_) :: HCoordA   !} Intermediate horizontal coord systems.
  Type(HCoord_) :: HCoordB   !}
  Real(Std)     :: Point(2)  !  Coords of point in standard lat-long coord system.
  Real(Std)     :: Point0(2) !  Coords of point in rotated lat-long coord system.
  Real(Std)     :: PointA(2) !] Coords of point in intermediate coord systems.
  Real(Std)     :: PointB(2) !]
  Real(Std)     :: Lon          !} Values used in correction angle calculations.
  Real(Std)     :: Lat          !}
  Real(Std)     :: Lon0         !}
  Real(Std)     :: Lat0         !}
  Real(Std)     :: LonPrime     !}
  Real(Std)     :: LatPrime     !}
  Real(Std)     :: LonStar      !}
  Real(Std)     :: PolePoint(2) !}
  Real(Std)     :: EdgeA        !}
  Real(Std)     :: EdgeB        !}
  Real(Std)     :: EdgeC        !}
  Real(Std)     :: AngleA       !}
  Real(Std)     :: AngleB       !}
  Real(Std)     :: S            !}
  Real(Std)     :: C            !}
  Real(Std)     :: CorrAngle    !}

# ifdef ExtraChecks
    If (HCoord1%CoordType == H_PSPolar .or. &
        HCoord1%CoordType == H_TMPolar) Then
      Call Message(                                                        &
             'Error in CalcHAngle: polar coords in HCoord1 not supported', &
             4                                                             &
           )
    End If
    If (HCoord2%CoordType == H_PSPolar .or. &
        HCoord2%CoordType == H_TMPolar) Then
      Call Message(                                                        &
             'Error in CalcHAngle: polar coords in HCoord2 not supported', &
             4                                                             &
           )
    End If
# endif

  If (HCoord1%CoordType == HCoord2%CoordType  .and. &
      HCoord1%Pole(1)   == HCoord2%Pole(1)    .and. &
      HCoord1%Pole(2)   == HCoord2%Pole(2)) Then

    Select Case (HCoord1%CoordType)

      Case (H_LatLong)

        Angle = 0.0 ! Easterly orientation is same in both lat-long systems.

      Case (H_PSCartesian)

        Angle = HCoord2%Angle - HCoord1%Angle ! Rotation angle between systems.

      Case (H_TMCartesian)

        Angle = 0.0 ! Orientation of grid Easting is same in both systems.

      Case Default

        Call Message('Error in CalcHAngle', 4)

    End Select

  Else

    ! Process PS/TM input systems: convert from HCoord1 to a lat-long system HCoordA.
    Select Case (HCoord1%CoordType)

      Case (H_LatLong) ! No conversion required.

        HCoordA = HCoord1
        PointA  = PointIn1
        Angle   = 0.0

      Case (H_PSCartesian) ! Convert to a lat-long system with same pole.

        HCoordA = HCoord1 ; HCoordA%CoordType = H_LatLong
        HCoordA%Origin = (/ 0.0, 0.0 /) ; HCoordA%Unit = (/ 1.0, 1.0 /)
        PointA = ConvertH(HCoord1, HCoordA, PointIn1)
        Angle = PointA(1) ! Longitude of point in intermediate coord system.

      Case (H_TMCartesian) ! Convert to the standard lat-long system.

        HCoordA = HCoord_LatLongRadians()
        PointA  = ConvertH(HCoord1, HCoordA, PointIn1)
        ! Set up rotated lat-long coord system corresponding to TM projection.
        ! Note: (0, 0) of HCoord0 is positioned at true origin of TM projection.
        Lon0 = HCoord1%Pole(1) - Pi/2.0
        Lat0 = Pi/2.0 ! Actually angle from true north pole.
        HCoord0 = HCoordA
        HCoord0%Pole  = (/ Lon0, Lat0 /)
        HCoord0%Angle = HCoord1%Pole(2) + Pi/2.0
        Point0 = ConvertH(HCoordA, HCoord0, PointA)
        ! Use spherical trigonometric formulae to calculate angle correction;
        ! see documentation for details.
        ! Construct spherical triangle.
        EdgeA  = Pi/2.0 - Point0(2)
        EdgeB  = Pi/2.0 - PointA(2)
        EdgeC  = Pi/2.0
        AngleA = PointA(1) - Lon0
        AngleB = (Pi/2.0 - HCoord1%Pole(2)) - Point0(1)
        ! Calculate angle correction.
        S = Sin(EdgeC) * Sin(AngleA) / Sin(EdgeA)
        C = (Sin(AngleA) * Sin(AngleB) * Cos(EdgeC)) - (Cos(AngleA) * Cos(AngleB))
        Angle = Pi/2.0 - ATan2(S, C) ! Angle from grid north to true north.

      Case Default

        Call Message('Error in CalcHAngle', 4)

    End Select

    HCoord = HCoord_LatLongRadians()

    ! Convert to the standard lat-long coord system.
    ! If the north pole of the intermediate coord system HCoordA is not at the
    ! standard north pole then an angle correction is required. Note that
    ! when the two lat-long coord systems are parallel (i.e. same pole) then
    ! no angle correction is necessary.
    If (HCoordA%Pole(2) > Small) Then
      ! Use spherical trigonometric formulae to calculate angle correction;
      ! see documentation for details.
      ! Point P.
      LonPrime = PointA(1) * HCoordA%Unit(1) + HCoordA%Origin(1) + HCoordA%Angle
      LatPrime = PointA(2) * HCoordA%Unit(2) + HCoordA%Origin(2)
      Point    = ConvertH(HCoordA, HCoord, PointA)
      Lon      = Point(1)
      Lat      = Point(2)
      ! Points N and N'.
      PolePoint = (/ 0.0, Pi/2.0 /)
      PolePoint = ConvertH(HCoord, HCoordA, PolePoint)
      LonStar   = PolePoint(1) * HCoordA%Unit(1) + HCoordA%Origin(1) + HCoordA%Angle
      Lon0      = HCoordA%Pole(1)
      Lat0      = HCoordA%Pole(2) ! Actually angle from true north pole.
      ! Construct spherical triangle.
      EdgeA  = Pi/2.0 - LatPrime
      EdgeB  = Pi/2.0 - Lat
      EdgeC  = Lat0
      AngleA = Lon - Lon0
      AngleB = LonStar - LonPrime
      ! Calculate angle correction.
      S = Sin(EdgeC) * Sin(AngleA) / Sin(EdgeA)
      C = (Sin(AngleA) * Sin(AngleB) * Cos(EdgeC)) - (Cos(AngleA) * Cos(AngleB))
      Angle = Angle - ATan2(S, C) ! Note: negative rotation for N' -> N.
    Else
      Point = ConvertH(HCoordA, HCoord, PointA)
    End If

    ! Identify the lat-long coord system HCoordB appropriate to HCoord2.
    Select Case (HCoord2%CoordType)

      Case (H_LatLong) ! HCoord2 directly.

        HCoordB = HCoord2

      Case (H_PSCartesian) ! A lat-long system with same pole as HCoord2.

        HCoordB = HCoord2 ; HCoordB%CoordType = H_LatLong
        HCoordB%Origin = (/ 0.0, 0.0 /) ; HCoordB%Unit = (/ 1.0, 1.0 /)

      Case (H_TMCartesian) ! The standard lat-long system.

        HCoordB = HCoord_LatLongRadians()

      Case Default

        Call Message('Error in CalcHAngle', 4)

    End Select

    ! Convert from the standard lat-long coord system to HCoordB.
    ! If the north pole of the intermediate coord system HCoordB is not at the
    ! standard north pole then an angle correction has to be applied. Note that
    ! when the two lat-long coord systems are parallel (i.e. same pole)
    ! then no angle correction is necessary.
    If (HCoordB%Pole(2) > Small) Then
      ! Use spherical trigonometric formulae to calculate angle correction;
      ! see documentation for details.
      ! Point P.
      Lon      = Point(1)
      Lat      = Point(2)
      PointB   = ConvertH(HCoord, HCoordB, Point)
      LonPrime = PointB(1) * HCoordB%Unit(1) + HCoordB%Origin(1) + HCoordB%Angle
      LatPrime = PointB(2) * HCoordB%Unit(2) + HCoordB%Origin(2)
      ! Points N and N'.
      PolePoint = (/ 0.0, Pi/2.0 /)
      PolePoint = ConvertH(HCoord, HCoordB, PolePoint)
      LonStar   = PolePoint(1) * HCoordB%Unit(1) + HCoordB%Origin(1) + HCoordB%Angle
      Lon0      = HCoordB%Pole(1)
      Lat0      = HCoordB%Pole(2) ! Actually angle from true north pole.
      ! Construct spherical triangle.
      EdgeA  = Pi/2.0 - LatPrime
      EdgeB  = Pi/2.0 - Lat
      EdgeC  = Lat0
      AngleA = Lon - Lon0
      AngleB = LonStar - LonPrime
      ! Calculate angle correction.
      S = Sin(EdgeC) * Sin(AngleA) / Sin(EdgeA)
      C = (Sin(AngleA) * Sin(AngleB) * Cos(EdgeC)) - (Cos(AngleA) * Cos(AngleB))
      Angle = Angle + ATan2(S, C) ! Note: positive rotation for N -> N'.
    Else
      PointB = ConvertH(HCoord, HCoordB, Point)
    End If

    ! Process PS/TM output systems: convert from lat-long system HCoordB to HCoord2.
    If (HCoord2%CoordType == H_PSCartesian) Then
      Angle = Angle - PointB(1) ! Subtract longitude of point in coord system HCoordB.
    Else If (HCoord2%CoordType == H_TMCartesian) Then
        ! Set up rotated lat-long coord system corresponding to TM projection.
        ! Note: (0, 0) of HCoord0 is positioned at true origin of TM projection.
        Lon0 = HCoord2%Pole(1) - Pi/2.0
        Lat0 = Pi/2.0 ! Actually angle from true north pole.
        HCoord0 = HCoordB
        HCoord0%Pole  = (/ Lon0, Lat0 /)
        HCoord0%Angle = HCoord2%Pole(2) + Pi/2.0
        Point0 = ConvertH(HCoordB, HCoord0, PointB)
        ! Use spherical trigonometric formulae to calculate angle correction;
        ! see documentation for details.
        ! Construct spherical triangle.
        EdgeA  = Pi/2.0 - Point0(2)
        EdgeB  = Pi/2.0 - PointB(2)
        EdgeC  = Pi/2.0
        AngleA = PointB(1) - Lon0
        AngleB = (Pi/2.0 - HCoord2%Pole(2)) - Point0(1)
        ! Calculate angle correction.
        S = Sin(EdgeC) * Sin(AngleA) / Sin(EdgeA)
        C = (Sin(AngleA) * Sin(AngleB) * Cos(EdgeC)) - (Cos(AngleA) * Cos(AngleB))
        CorrAngle = ATan2(S, C) - Pi/2.0 ! Angle from true north to grid north.
      Angle = Angle + CorrAngle
    End If

    ! Check that output angle is in the range [-pi, pi].
    If (Angle >  Pi) Angle = Angle - 2.0*Pi
    If (Angle < -Pi) Angle = Angle + 2.0*Pi

  End If

End Subroutine CalcHAngle

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcdZdXZBased(ZCoord1, ZCoord2, ZIn1, Topog, dTopogdX, dZdX)
! Calculates dZ_1/d(X,Y,Z_2) for two height-based coord systems Z_1 and Z_2
! with horizontal units of metres.
! $$ Should this be horizontal units of dTopogdX and dZdX need to agree?
! $$ Rename CalcdZdXZBased
! Note this routine can also calculate dZ_1/dT (at constant X, Y and Z_2) if dTopogdX(1) is the time
! derivative d Topog / dT. The result is then returned in dZdX(1).

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In)  :: ZCoord1     ! Primary vertical coord (dependent variable).
  Type(ZCoord_), Intent(In)  :: ZCoord2     ! Secondary vertical coord (independent variable).
  Real(Std),     Intent(In)  :: ZIn1        ! Height in coord system ZCoord1.
  Real(Std),     Intent(In)  :: Topog       ! Topographic height in metres above sea level.
  Real(Std),     Intent(In)  :: dTopogdX(2) ! Horizontal gradient of topographic height above sea level.
  Real(Std),     Intent(Out) :: dZdX(3)     ! dZ_1/d(X,Y,Z_2)
  ! Locals:
  Real(Std) :: Z    ! Height in metres above sea level.
  Real(Std) :: H    ! Model top height in a height-based eta coord system.
  Real(Std) :: Eta  ! Height in an eta coord system.
  Real(Std) :: EtaI ! Eta value on the linear-quadratic interface of an eta system.
  Real(Std) :: R    ! Orography damping component for an eta system (= 1 - eta/etaI).
  Real(Std) :: SZ   !] Expressions used in eta calculations.
  Real(Std) :: SXY  !]

  Select Case (ZCoord1%CoordType)

    Case (Z_AboveGround)

      Select Case (ZCoord2%CoordType)

        Case (Z_AboveGround)

          dZdX(1) = 0.0
          dZdX(2) = 0.0
          dZdX(3) = ZCoord2%Unit / ZCoord1%Unit

        Case (Z_AboveSea)

          ! Height gradients in metres per unit dX.
          dZdX(1) = -dTopogdX(1)
          dZdX(2) = -dTopogdX(2)
          dZdX(3) = ZCoord2%Unit
          ! Gradients scaled by the Z_1 coord unit.
          dZdX = dZdX / ZCoord1%Unit

        Case (Z_ZAsEta)

          ! Calculate height in metres above sea level.
          Z = ZIn1 * ZCoord1%Unit + Topog
          H = ZCoord2%ModelTopHeight

          If (Z < ZCoord2%InterfaceHeight) Then
            ! Lower subdomain calculation.
            Eta  = AboveGroundToZAsEta(ZCoord1, ZCoord2, Topog, ZIn1)
            EtaI = ZCoord2%InterfaceEta
            R    = 1.0 - Eta/EtaI
            ! Height gradients in metres per unit dX.
            dZdX(1) = (R * R - 1) * dTopogdX(1)
            dZdX(2) = (R * R - 1) * dTopogdX(2)
            dZdX(3) = H - 2.0 * (R / EtaI) * Topog
            ! Gradients scaled by the Z_1 coord unit.
            dZdX = dZdX / ZCoord1%Unit
          Else
            ! Upper subdomain calculation.
            ! Height gradients in metres per unit dX.
            dZdX(1) = -dTopogdX(1)
            dZdX(2) = -dTopogdX(2)
            dZdX(3) = H
            ! Gradients scaled by the Z_1 coord unit.
            dZdX = dZdX / ZCoord1%Unit
          End If

        Case Default

          Call Message('Error in CalcdZdXZBased', 4)

      End Select

    Case (Z_AboveSea)

      Select Case (ZCoord2%CoordType)

        Case (Z_AboveGround)

          ! Height gradients in metres per unit dX.
          dZdX(1) = dTopogdX(1)
          dZdX(2) = dTopogdX(2)
          dZdX(3) = ZCoord2%Unit
          ! Gradients scaled by the Z_1 coord unit.
          dZdX = dZdX / ZCoord1%Unit

        Case (Z_AboveSea)

          dZdX(1) = 0.0
          dZdX(2) = 0.0
          dZdX(3) = ZCoord2%Unit / ZCoord1%Unit

        Case (Z_ZAsEta)

          ! Calculate height in metres above sea level.
          Z = ZIn1 * ZCoord1%Unit
          H = ZCoord2%ModelTopHeight

          If (Z < ZCoord2%InterfaceHeight) Then
            ! Lower subdomain calculation.
            Eta  = AboveSeaToZAsEta(ZCoord1, ZCoord2, Topog, ZIn1)
            EtaI = ZCoord2%InterfaceEta
            R    = 1.0 - Eta/EtaI
            ! Height gradients in metres per unit dX.
            dZdX(1) = R * R * dTopogdX(1)
            dZdX(2) = R * R * dTopogdX(2)
            dZdX(3) = H - 2.0 * (R / EtaI) * Topog
            ! Gradients scaled by the Z_1 coord unit.
            dZdX = dZdX / ZCoord1%Unit
          Else
            ! Upper subdomain calculation.
            dZdX(1) = 0.0
            dZdX(2) = 0.0
            dZdX(3) = H / ZCoord1%Unit
          End If

        Case Default

          Call Message('Error in CalcdZdXZBased', 4)

      End Select

    Case (Z_ZAsEta)

      Select Case (ZCoord2%CoordType)

        Case (Z_AboveGround)

          ! Set up information on eta coord system.
          Eta  = ZIn1
          EtaI = ZCoord1%InterfaceEta
          H    = ZCoord1%ModelTopHeight

          If (Eta < EtaI) Then
            ! Lower subdomain calculation.
            R   = 1.0 - Eta/EtaI
            SZ  = H - 2.0 * (R / EtaI) * Topog
            SXY = (1.0 - R * R) / SZ
            dZdX(1) = SXY * dTopogdX(1)
            dZdX(2) = SXY * dTopogdX(2)
            dZdX(3) = ZCoord2%Unit / SZ
          Else
            ! Upper subdomain calculation.
            dZdX(1) = dTopogdX(1) / H
            dZdX(2) = dTopogdX(2) / H
            dZdX(3) = ZCoord2%Unit / H
          End If

        Case (Z_AboveSea)

          ! Set up information on eta coord system.
          Eta  = ZIn1
          EtaI = ZCoord1%InterfaceEta
          H    = ZCoord1%ModelTopHeight

          If (Eta < EtaI) Then
            ! Lower subdomain calculation.
            R   = 1.0 - Eta/EtaI
            SZ  = H - 2.0 * (R / EtaI) * Topog
            SXY = -(R * R / SZ)
            dZdX(1) = SXY * dTopogdX(1)
            dZdX(2) = SXY * dTopogdX(2)
            dZdX(3) = ZCoord2%Unit / SZ
          Else
            ! Upper subdomain calculation.
            dZdX(1) = 0.0
            dZdX(2) = 0.0
            dZdX(3) = ZCoord2%Unit / H
          End If

        Case (Z_ZAsEta)

          Call Message('Error in CalcdZdXZBased: eta -> eta option not supported', 4)

        Case Default

          Call Message('Error in CalcdZdXZBased', 4)

      End Select

    Case Default

      Call Message('Error in CalcdZdXZBased', 4)

  End Select

End Subroutine CalcdZdXZBased

!-------------------------------------------------------------------------------------------------------------

Subroutine CalcdZdXPBased(ZCoord1, ZCoord2, ZIn1, PStar, dPStardX, dZdX)
! Calculates dZ_1/d(X,Y,Z_2) for two pressure-based coord systems Z_1 and Z_2.
! $$ Currently only handles coord systems Z_P and Z_PAsEta.
! Note this routine can also calculate dZ_1/dT (at constant X, Y and Z_2) if dPStardX(1) is the time
! derivative d PStar / dT. The result is then returned in dZdX(1).

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In)  :: ZCoord1     ! Primary vertical coord (dependent variable).
  Type(ZCoord_), Intent(In)  :: ZCoord2     ! Secondary vertical coord (independent variable).
  Real(Std),     Intent(In)  :: ZIn1        ! Height in coord system ZCoord1.
  Real(Std),     Intent(In)  :: PStar       ! Surface pressure in Pa.
  Real(Std),     Intent(In)  :: dPStardX(2) ! Horizontal gradient of surface pressure.
  Real(Std),     Intent(Out) :: dZdX(3)     ! dZ_1/d(X,Y,Z_2).
  ! Locals:
  Real(Std) :: P        ! Pressure (in Pa).
  Real(Std) :: Eta      ! Height in an eta coord system.
  Integer   :: EtaLevel ! Eta level index defining the base of a model layer.
  Integer   :: i        ! Loop index.
  Real(Std) :: PBelow   !] Pressures and eta values at base/top of a model layer.
  Real(Std) :: PAbove   !]
  Real(Std) :: EtaBelow !]
  Real(Std) :: EtaAbove !]
  Real(Std) :: B        ! Interpolated value of B coefficients.

  Select Case (ZCoord1%CoordType)

    Case (Z_P)

      Select Case (ZCoord2%CoordType)

        Case (Z_P)

          dZdX(1) = 0.0
          dZdX(2) = 0.0
          dZdX(3) = ZCoord2%Unit / ZCoord1%Unit

        Case (Z_PAsEta)

          ! Calculate pressure in Pa.
          P = ZIn1 * ZCoord1%Unit

          ! Determine eta level of the layer, i.e. z in (EtaLevel, EtaLevel + 1].
          EtaLevel = ZCoord2%nEtaLevels
          Do i = 1, ZCoord2%nEtaLevels          ! $$ Use bisection
            If (P >= ZCoord2%A(i) + ZCoord2%B(i)*PStar) Then
              EtaLevel = i - 1
              Exit
            End If
          End Do
          ! If at surface, use the first and second eta levels for calculating vertical gradient.
          If (EtaLevel == 0) EtaLevel = 1

          ! For horizontal gradient: interpolate B to pressure P.
          ! For vertical gradient: calculate P and Eta on neighbouring eta levels.
          If (EtaLevel == ZCoord2%nEtaLevels) Then
            ! Above highest eta level
            PBelow   = ZCoord2%A(EtaLevel) + ZCoord2%B(EtaLevel)*PStar
            PAbove   = 0.0
            EtaBelow = ZCoord2%Eta(EtaLevel)
            EtaAbove = 0.0
            B        = ZCoord2%B(EtaLevel) * P / PBelow
          Else
            PBelow   = ZCoord2%A(EtaLevel    ) + ZCoord2%B(EtaLevel    )*PStar
            PAbove   = ZCoord2%A(EtaLevel + 1) + ZCoord2%B(EtaLevel + 1)*PStar
            EtaBelow = ZCoord2%Eta(EtaLevel    )
            EtaAbove = ZCoord2%Eta(EtaLevel + 1)
            B        = ZCoord2%B(EtaLevel) + (ZCoord2%B(EtaLevel + 1) - ZCoord2%B(EtaLevel)) * &
                         (PBelow - P)/(PBelow - PAbove)
          End If

          ! Compute gradients in Pa per unit dX.
          dZdX(1) = B * dPStardX(1)
          dZdX(2) = B * dPStardX(2)
          dZdX(3) = (PBelow - PAbove)/(EtaBelow - EtaAbove)
          ! Gradients scaled by the Z_1 coord unit.
          dZdX(:) = dZdX(:) / ZCoord1%Unit

        Case Default

          Call Message('FATAL ERROR in CalcdZdXPBased: unsupported option', 4)

      End Select

    Case (Z_PAsEta)

      Select Case (ZCoord2%CoordType)

        Case (Z_P)

          ! Note: first calculate dp/d(x, y, eta) on eta (model) levels. Then invert matrix
          !       to get the required d eta/d(x, y, p) on pressure surfaces.

          ! Height in eta coord system.
          Eta = ZIn1

          ! Determine eta level of the layer, i.e. z in (EtaLevel, EtaLevel + 1].
          EtaLevel = ZCoord1%nEtaLevels
          Do i = 1, ZCoord1%nEtaLevels          ! $$ Use bisection
            If (Eta >= ZCoord1%Eta(i)) Then
              EtaLevel = i - 1
              Exit
            End If
          End Do
          ! If at surface, use the first and second eta levels for calculating vertical gradient.
          If (EtaLevel == 0) EtaLevel = 1

          ! For horizontal gradient: interpolate B to pressure P.
          ! For vertical gradient: calculate P and Eta on neighbouring eta levels.
          If (EtaLevel == ZCoord1%nEtaLevels) Then
            ! Above highest eta level
            PBelow   = ZCoord1%A(EtaLevel) + ZCoord1%B(EtaLevel)*PStar
            PAbove   = 0.0
            EtaBelow = ZCoord1%Eta(EtaLevel)
            EtaAbove = 0.0
            B        = ZCoord1%B(EtaLevel) * Eta / EtaBelow
          Else
            PBelow   = ZCoord1%A(EtaLevel    ) + ZCoord1%B(EtaLevel    )*PStar
            PAbove   = ZCoord1%A(EtaLevel + 1) + ZCoord1%B(EtaLevel + 1)*PStar
            EtaBelow = ZCoord1%Eta(EtaLevel    )
            EtaAbove = ZCoord1%Eta(EtaLevel + 1)
            B        = ZCoord1%B(EtaLevel) + (ZCoord1%B(EtaLevel + 1) - ZCoord1%B(EtaLevel)) * &
                         (EtaBelow - Eta)/(EtaBelow - EtaAbove)
          End If

          ! Compute gradients in Pa per unit dX.
          dZdX(1) = B * dPStardX(1)
          dZdX(2) = B * dPStardX(2)
          dZdX(3) = (PBelow - PAbove)/(EtaBelow - EtaAbove)
          ! Gradients scaled by the Z_2 coord unit.
          dZdX(:) = dZdX(:) / ZCoord2%Unit

          ! Invert matrix to obtain d eta/d(x, y, p) on pressure surfaces.
          dZdX(1) = - dZdX(1) / dZdX(3)
          dZdX(2) = - dZdX(2) / dZdX(3)
          dZdX(3) = 1.0 / dZdX(3)

        Case (Z_PAsEta)

          Call Message('FATAL ERROR in CalcdZdXPBased: eta -> eta option not supported', 4)

        Case Default

          Call Message('FATAL ERROR in CalcdZdXPBased: unsupported option', 4)

      End Select

    Case Default

      Call Message('FATAL ERROR in CalcdZdXPBased: unsupported option', 4)

  End Select

End Subroutine CalcdZdXPBased

!-------------------------------------------------------------------------------------------------------------

Function CalcdZdZZBased(ZCoord1, ZCoord2, ZIn1, Topog) Result(dZdZ)
! Calculates dZCoord1/dZCoord2 for two height-based coord systems dZCoord1 and
! dZCoord2.

! If one is eta, better for this to be ZCoord1.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In)           :: ZCoord1 !} The vertical coord systems.
  Type(ZCoord_), Intent(In)           :: ZCoord2 !}
  Real(Std),     Intent(In), Optional :: ZIn1    ! Height in coord system ZCoord1.
  Real(Std),     Intent(In), Optional :: Topog   ! Topographic height in metres above
                                                 ! sea level.
  ! ZIn1 and Topog are only needed if one of the coord systems is eta.
  ! Result.
  Real(Std) :: dZdZ ! dZ_1/dZ_2.
  ! Locals:
  Real(Std) :: Z    ! Height in metres above sea level.
  Real(Std) :: H    ! Model top height in a height-based eta coord system.
  Real(Std) :: Eta  ! Height in an eta coord system.
  Real(Std) :: EtaI ! Eta value on the linear-quadratic interface of an eta system.
  Real(Std) :: R    ! Orography damping component for an eta system (= 1 - eta/etaI).
  Real(Std) :: SZ   !] Expressions used in eta calculations.
  Real(Std) :: SXY  !]

  Select Case (ZCoord1%CoordType)

    Case (Z_AboveGround)

      Select Case (ZCoord2%CoordType)

        Case (Z_AboveGround)

          dZdZ = ZCoord2%Unit / ZCoord1%Unit

        Case (Z_AboveSea)

          dZdZ = ZCoord2%Unit / ZCoord1%Unit

        Case (Z_ZAsEta)

          If (.not.Present(ZIn1) .or. .not.Present(Topog)) Then
            Call Message('UNEXPECTED FATAL ERROR in CalcdZdZZBased', 4)
          End If

          ! Calculate height in metres above sea level.
          Z = ZIn1 * ZCoord1%Unit + Topog
          H = ZCoord2%ModelTopHeight

          If (Z < ZCoord2%InterfaceHeight) Then
            ! Lower subdomain calculation.
            Eta  = AboveGroundToZAsEta(ZCoord1, ZCoord2, Topog, ZIn1)
            EtaI = ZCoord2%InterfaceEta
            R    = 1.0 - Eta/EtaI
            ! Height gradients in metres per unit dX.
            dZdZ = H - 2.0 * (R / EtaI) * Topog
            ! Gradients scaled by the Z_1 coord unit.
            dZdZ = dZdZ / ZCoord1%Unit
          Else
            ! Upper subdomain calculation.
            ! Height gradients in metres per unit dX.
            dZdZ = H
            ! Gradients scaled by the Z_1 coord unit.
            dZdZ = dZdZ / ZCoord1%Unit
          End If

        Case Default

          Call Message('Error in CalcdZdZZBased', 4)

      End Select

    Case (Z_AboveSea)

      Select Case (ZCoord2%CoordType)

        Case (Z_AboveGround)

          dZdZ = ZCoord2%Unit / ZCoord1%Unit

        Case (Z_AboveSea)

          dZdZ = ZCoord2%Unit / ZCoord1%Unit

        Case (Z_ZAsEta)

          If (.not.Present(ZIn1) .or. .not.Present(Topog)) Then
            Call Message('UNEXPECTED FATAL ERROR in CalcdZdZZBased', 4)
          End If

          ! Calculate height in metres above sea level.
          Z = ZIn1 * ZCoord1%Unit
          H = ZCoord2%ModelTopHeight

          If (Z < ZCoord2%InterfaceHeight) Then
            ! Lower subdomain calculation.
            Eta  = AboveSeaToZAsEta(ZCoord1, ZCoord2, Topog, ZIn1)
            EtaI = ZCoord2%InterfaceEta
            R    = 1.0 - Eta/EtaI
            ! Height gradients in metres per unit dX.
            dZdZ = H - 2.0 * (R / EtaI) * Topog
            ! Gradients scaled by the Z_1 coord unit.
            dZdZ = dZdZ / ZCoord1%Unit
          Else
            ! Upper subdomain calculation.
            dZdZ = H / ZCoord1%Unit
          End If

        Case Default

          Call Message('Error in CalcdZdZZBased', 4)

      End Select

    Case (Z_ZAsEta)

      If (.not.Present(ZIn1) .or. .not.Present(Topog)) Then
        Call Message('UNEXPECTED FATAL ERROR in CalcdZdZZBased', 4)
      End If

      Select Case (ZCoord2%CoordType)

        Case (Z_AboveGround)

          ! Set up information on eta coord system.
          Eta  = ZIn1
          EtaI = ZCoord1%InterfaceEta
          H    = ZCoord1%ModelTopHeight

          If (Eta < EtaI) Then
            ! Lower subdomain calculation.
            R   = 1.0 - Eta/EtaI
            SZ  = H - 2.0 * (R / EtaI) * Topog
            SXY = (1.0 - R * R) / SZ
            dZdZ = ZCoord2%Unit / SZ
          Else
            ! Upper subdomain calculation.
            dZdZ = ZCoord2%Unit / H
          End If

        Case (Z_AboveSea)

          ! Set up information on eta coord system.
          Eta  = ZIn1
          EtaI = ZCoord1%InterfaceEta
          H    = ZCoord1%ModelTopHeight

          If (Eta < EtaI) Then
            ! Lower subdomain calculation.
            R   = 1.0 - Eta/EtaI
            SZ  = H - 2.0 * (R / EtaI) * Topog
            SXY = -(R * R / SZ)
            dZdZ = ZCoord2%Unit / SZ
          Else
            ! Upper subdomain calculation.
            dZdZ = ZCoord2%Unit / H
          End If

        Case (Z_ZAsEta)

          Call Message('Error in CalcdZdZZBased: eta -> eta option not supported', 4)

        Case Default

          Call Message('Error in CalcdZdZZBased', 4)

      End Select

    Case Default

      Call Message('Error in CalcdZdZZBased', 4)

  End Select

End Function CalcdZdZZBased

!-------------------------------------------------------------------------------------------------------------

Function CalcdZdZPBased(ZCoord1, ZCoord2, ZIn1, PS) Result(dZdZ)
! Calculates dZCoord1/dZCoord2 for two pressure-based coord systems ZCoord1 and ZCoord2.

! If one is eta, better for this to be ZCoord1.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In)           :: ZCoord1 !} The vertical coord systems.
  Type(ZCoord_), Intent(In)           :: ZCoord2 !}
  Real(Std),     Intent(In), Optional :: ZIn1    ! Height in coord system ZCoord1.
  Real(Std),     Intent(In), Optional :: PS      ! Surface pressure.
  ! ZIn1 and PS are needed if one of the coord systems is eta.
  ! ZIn1 is needed if one of the coord systems is PAsZ.
  ! Function result.
  Real(Std) :: dZdZ ! dZ_1/dZ_2.
  ! Locals:
  Real(Std) :: P   ! Pressure (Pa).
  Real(Std) :: Z   ! Height above sea level (m).
  Real(Std) :: T   ! Temperature (K).
  Real(Std) :: Rho ! Density.

  Select Case (ZCoord1%CoordType)

    Case (Z_P)

      Select Case (ZCoord2%CoordType)

        Case (Z_P)

          dZdZ = ZCoord2%Unit / ZCoord1%Unit

        Case (Z_PAsZ)

          If (.not.Present(ZIn1)) Then
            Call Message('UNEXPECTED FATAL ERROR in CalcdZdZPBased', 4)
          End If
          Call ICAOAtP(ZIn1 * ZCoord1%Unit, Z, T, Rho)
          dZdZ = Rho * Gravity * ZCoord2%Unit / ZCoord1%Unit

        Case (Z_PAsEta)

          Call Message('Error in CalcdZdZPBased: eta option not supported yet', 4)

        Case Default

          Call Message('Error in CalcdZdZPBased', 4)

      End Select

    Case (Z_PAsZ)

      Select Case (ZCoord2%CoordType)

        Case (Z_P)

          If (.not.Present(ZIn1)) Then
            Call Message('UNEXPECTED FATAL ERROR in CalcdZdZPBased', 4)
          End If
          Call ICAOAtZ(ZIn1 * ZCoord1%Unit, P, T, Rho)
          dZdZ = ZCoord2%Unit / (ZCoord1%Unit * Rho * Gravity)

        Case (Z_PAsZ)

          dZdZ = ZCoord2%Unit / ZCoord1%Unit

        Case (Z_PAsEta)

          Call Message('Error in CalcdZdZPBased: eta option not supported yet', 4)

        Case Default

          Call Message('Error in CalcdZdZPBased', 4)

      End Select

    Case (Z_PAsEta)

      Select Case (ZCoord2%CoordType)

        Case (Z_P)

          Call Message('Error in CalcdZdZPBased: eta option not supported yet', 4)

        Case (Z_PAsZ)

          Call Message('Error in CalcdZdZPBased: eta option not supported yet', 4)

        Case (Z_PAsEta)

          Call Message('Error in CalcdZdZPBased: eta -> eta option not supported', 4)

        Case Default

          Call Message('Error in CalcdZdZPBased', 4)

      End Select

    Case Default

      Call Message('Error in CalcdZdZPBased', 4)

  End Select

End Function CalcdZdZPBased

!-------------------------------------------------------------------------------------------------------------

Function InitZCoord(Name, CoordType, Unit, EtaDefn, ModelTopHeight, InterfaceHeight)
! Initialises a vertical coord system.

  Implicit None
  ! Argument list:
  Character(*),   Intent(In)           :: Name            ! Name of coord system.
  Integer,        Intent(In)           :: CoordType       !
  Real(Std),      Intent(In)           :: Unit            !
  Type(EtaDefn_), Intent(In), Optional :: EtaDefn         ! Eta definition.
  Real(Std),      Intent(In), Optional :: ModelTopHeight  !} Parameters in height-based
  Real(Std),      Intent(In), Optional :: InterfaceHeight !} eta coord system.
  ! Function result:
  Type(ZCoord_) :: InitZCoord !
  ! Locals parameters:
  Real(Std), Parameter :: P0 = 100000.0_Std ! Reference pressure for pressure-based eta definition (Pa).
  ! Locals:
  Real(Std)     :: P0Min  ! Minimum surface pressure for which p decreases with level number (Pa).
  Type(ZCoord_) :: ZCoord ! Local copy of function result.
  Integer       :: i      !

  If (CoordType < 1 .or. CoordType > 6) Then
    Call Message('UNEXPECTED FATAL ERROR in InitZCoord: invalid coord type', 4)
  End If
  If (Unit == 0.0) Then
    Call Message('FATAL ERROR in InitZCoord: unit should be non-zero', 3)
  End If
  If (Len_Trim(Name) > MaxCharLength) Then
    Call Message('FATAL ERROR in InitZCoord: name is too long', 3)
  End If
  If (CoordType == Z_PAsEta .and. .not.Present(EtaDefn)) Then
    Call Message('UNEXPECTED FATAL ERROR in InitZCoord: missing pressure-eta information', 4)
  End If
  If (CoordType == Z_ZAsEta .and. (.not.Present(ModelTopHeight) .or. .not.Present(InterfaceHeight))) Then
    Call Message('UNEXPECTED FATAL ERROR in InitZCoord: missing height-eta information', 4)
  End If

  ZCoord%Name      = Name
  ZCoord%CoordType = CoordType
  ZCoord%Unit      = Unit

  ! Pressure-based eta system.
  If (CoordType == Z_PAsEta) Then

    Allocate(ZCoord%Eta(EtaDefn%nEtaLevels))
    Allocate(ZCoord%A  (EtaDefn%nEtaLevels))
    Allocate(ZCoord%B  (EtaDefn%nEtaLevels))

    ZCoord%nEtaLevels = EtaDefn%nEtaLevels

    If (ZCoord%nEtaLevels <= 2) Then
      Call Message('FATAL ERROR in InitZCoord: invalid pressure-eta information', 3)
    End If

    If (EtaDefn%EtaDefnType == 'E') Then
      ZCoord%Eta(:) = EtaDefn%A(:) / P0 + EtaDefn%B(:)
      ZCoord%A  (:) = EtaDefn%A(:)
      ZCoord%B  (:) = EtaDefn%B(:)
    Else If (EtaDefn%EtaDefnType == 'A') Then
      ZCoord%Eta(:) = EtaDefn%Eta(:)
      ZCoord%A  (:) = (EtaDefn%Eta(:) - EtaDefn%B(:)) * P0
      ZCoord%B  (:) = EtaDefn%B(:)
    Else If (EtaDefn%EtaDefnType == 'B') Then
      ZCoord%Eta(:) = EtaDefn%Eta(:)
      ZCoord%A  (:) = EtaDefn%A(:)
      ZCoord%B  (:) = EtaDefn%Eta(:) - EtaDefn%A(:) / P0
    Else
      Call Message('UNEXPECTED FATAL ERROR in InitZCoord', 4)
    End If

    If (ZCoord%Eta(1) /= 1.0 .or. ZCoord%A(1) /= 0.0 .or. ZCoord%B(1) /= 1.0) Then
      Call Message('FATAL ERROR in InitZCoord: invalid pressure-eta information', 3)
    End If

    ! Check that the Unit argument has value unity. ($$ note handled differently to
    !                                                non-optional, but unused, args in
    !                                                InitHCoord)
    If (Unit /= 1.0_Std) Then
      ZCoord%Unit = 1.0_Std
      Call Message(                                                       &
             'Warning: Unit not used for eta system; user input ignored', &
             1                                                            &
           )
    End If
    ! Check that eta and B decrease with increasing level number.
    Do i = 1, ZCoord%nEtaLevels - 1
      If (ZCoord%Eta(i) <= ZCoord%Eta(i + 1) .or. ZCoord%B(i) < ZCoord%B(i + 1)) Then
        Call Message(                                                                       &
               'FATAL ERROR in InitZCoord: invalid pressure-eta information'             // &
               'Eta must decrease and B must not increase with increasing level number',    &
               3                                                                            &
             )
      End If
    End Do
    ! Calculate minimum surface pressure for which pressure decreases with level
    ! number.
    ZCoord%PStarMin = 0.0
    Do i = 1, ZCoord%nEtaLevels - 1 ! $$ B(i) = B(i+1) case
      ZCoord%PStarMin = Max((ZCoord%A(i + 1) - ZCoord%A(i))/(ZCoord%B(i) - ZCoord%B(i + 1)), ZCoord%PStarMin)
    End Do
  End If

  ! Height-based eta system.
  If (CoordType == Z_ZAsEta) Then
    ! Check that the Unit argument has value unity.
    If (Unit /= 1.0_Std) Then
      ZCoord%Unit = 1.0_Std
      Call Message(                                                       &
             'Warning: Unit not used for eta system; user input ignored', &
             1                                                            &
           )
    End If
    ! Check that ModelTopHeight > InterfaceHeight > 0, and set
    ! the values of InterfaceEta and MaxOrogHeight.
    If (ModelTopHeight <= InterfaceHeight .or. InterfaceHeight <= 0.0) Then
      Call Message('Error in InitZCoord: invalid height-eta information', 4)
    End If
    ZCoord%ModelTopHeight  = ModelTopHeight
    ZCoord%InterfaceHeight = InterfaceHeight
    ZCoord%InterfaceEta    = InterfaceHeight / ModelTopHeight
    ZCoord%MaxOrogHeight   = 0.5 * InterfaceHeight
  End If

  Select Case (CoordType) ! $$ check - e.g. do we allow units to be < 0.0
    Case (Z_AboveGround, Z_AboveSea, Z_ZAsEta, Z_PAsZ)
      ZCoord%IncreasingUpwards = .true.
    Case (Z_P, Z_PAsEta)
      ZCoord%IncreasingUpwards = .false.
  End Select

  InitZCoord = ZCoord

End Function InitZCoord

!-------------------------------------------------------------------------------------------------------------

Function ZCoordEq(ZCoord1, ZCoord2)
! Tests for equality of vertical coord systems (including their names).

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoord1 ! First vertical coord system.
  Type(ZCoord_), Intent(In) :: ZCoord2 ! Second vertical coord system.
  ! Function result:
  Logical :: ZCoordEq ! Indicates if coord systems are equal.
  ! Locals:
  Integer :: i ! Loop variable.

  ZCoordEq = (ZCoord1%Name      .CIEq. ZCoord2%Name)     .and. &
              ZCoord1%CoordType   ==   ZCoord2%CoordType .and. &
              ZCoord1%Unit        ==   ZCoord2%Unit
  If (ZCoord1%CoordType == Z_PAsEta) Then
    ZCoordEq = ZCoordEq .and. ZCoord1%nEtaLevels == ZCoord2%nEtaLevels
    Do i = 1, ZCoord1%nEtaLevels
      ZCoordEq = ZCoordEq .and. ZCoord1%A(i) == ZCoord2%A(i)
      ZCoordEq = ZCoordEq .and. ZCoord1%B(i) == ZCoord2%B(i)
    End Do
  End If
  If (ZCoord1%CoordType == Z_ZAsEta) Then
    ZCoordEq = ZCoordEq                                           .and. &
               ZCoord1%ModelTopHeight  == ZCoord2%ModelTopHeight  .and. &
               ZCoord1%InterfaceHeight == ZCoord2%InterfaceHeight
  End If

End Function ZCoordEq

!-------------------------------------------------------------------------------------------------------------

Function ZCoordEquiv(ZCoord1, ZCoord2)
! Tests for equivalence of vertical coord systems (i.e. equality except possibly of
! their names).

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoord1 ! First vertical coord system.
  Type(ZCoord_), Intent(In) :: ZCoord2 ! Second vertical coord system.
  ! Function result:
  Logical :: ZCoordEquiv ! Indicates if coord systems are equivalent.
  ! Locals:
  Integer :: i ! Loop variable.

  ZCoordEquiv = ZCoord1%CoordType == ZCoord2%CoordType .and. &
                ZCoord1%Unit      == ZCoord2%Unit
  If (ZCoord1%CoordType == Z_PAsEta) Then
    ZCoordEquiv = ZCoordEquiv .and. ZCoord1%nEtaLevels == ZCoord2%nEtaLevels
    Do i = 1, ZCoord1%nEtaLevels
      ZCoordEquiv = ZCoordEquiv .and. ZCoord1%A(i) == ZCoord2%A(i)
      ZCoordEquiv = ZCoordEquiv .and. ZCoord1%B(i) == ZCoord2%B(i)
    End Do
  End If
  If (ZCoord1%CoordType == Z_ZAsEta) Then
    ZCoordEquiv = ZCoordEquiv                                        .and. &
                  ZCoord1%ModelTopHeight  == ZCoord2%ModelTopHeight  .and. &
                  ZCoord1%InterfaceHeight == ZCoord2%InterfaceHeight
  End If

End Function ZCoordEquiv

!-------------------------------------------------------------------------------------------------------------

Function ZCoord_m_agl()
! Returns the vertical coord system giving the height above ground in metres.

  Implicit None
  ! Function result:
  Type(ZCoord_) :: ZCoord_m_agl ! Required coord system.

  ZCoord_m_agl = InitZCoord('m agl', Z_AboveGround, 1.0_Std)

End Function ZCoord_m_agl

!-------------------------------------------------------------------------------------------------------------

Function ZCoord_m_asl()
! Returns the vertical coord system giving the height above sea in metres.

  Implicit None
  ! Function result:
  Type(ZCoord_) :: ZCoord_m_asl ! Required coord system.

  ZCoord_m_asl = InitZCoord('m asl', Z_AboveSea, 1.0_Std)

End Function ZCoord_m_asl

!-------------------------------------------------------------------------------------------------------------

Function ZCoord_Pa()
! Returns the vertical coord system giving pressure in Pascals.

  Implicit None
  ! Function result:
  Type(ZCoord_) :: ZCoord_Pa ! Required coord system.

  ZCoord_Pa = InitZCoord('Pa', Z_P, 1.0_Std)

End Function ZCoord_Pa

!-------------------------------------------------------------------------------------------------------------

Function ZCoord_FL()
! Returns the flight level vertical coord system (flight level is pressure,
! converted to height above sea level using the ICAO standard atmosphere and
! expressed in hundreds of feet).

  Implicit None
  ! Function result:
  Type(ZCoord_) :: ZCoord_FL ! Required coord system.

  ZCoord_FL = InitZCoord('FL', Z_PAsZ, 100.0*Foot)

End Function ZCoord_FL

!-------------------------------------------------------------------------------------------------------------

Function ZSameToSame(ZCoordIn, ZCoordOut, ZIn)
! Converts coords between coord systems of same type (but not eta coord systems).

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn  !
  Type(ZCoord_), Intent(In) :: ZCoordOut !
  Real(Std),     Intent(In) :: ZIn       !
  ! Function result:
  Real(Std) :: ZSameToSame !

# ifdef ExtraChecks
    If (ZCoordIn%CoordType /= ZCoordOut%CoordType .or. &
        ZCoordIn%CoordType == Z_PAsEta            .or. &
        ZCoordIn%CoordType == Z_ZAsEta) Then
      Call Message('Error in ZSameToSame', 4)
    End If
# endif

  ZSameToSame = (ZIn*ZCoordIn%Unit)/ZCoordOut%Unit

End Function ZSameToSame

!-------------------------------------------------------------------------------------------------------------

Function AboveGroundToAboveSea(ZCoordIn, ZCoordOut, GroundHeight, ZIn)
! Converts coords between coord systems.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn     !
  Type(ZCoord_), Intent(In) :: ZCoordOut    !
  Real(Std),     Intent(In) :: GroundHeight ! Height of ground.
  Real(Std),     Intent(In) :: ZIn          !
  ! Function result:
  Real(Std) :: AboveGroundToAboveSea !

# ifdef ExtraChecks
    If (ZCoordIn%CoordType  /= Z_AboveGround .or. &
        ZCoordOut%CoordType /= Z_AboveSea) Then
      Call Message('Error in AboveGroundToAboveSea', 4)
    End If
# endif

  AboveGroundToAboveSea = (ZIn*ZCoordIn%Unit + GroundHeight)/ZCoordOut%Unit

End Function AboveGroundToAboveSea

!-------------------------------------------------------------------------------------------------------------

Function AboveSeaToAboveGround(ZCoordIn, ZCoordOut, GroundHeight, ZIn)
! Converts coords between coord systems.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn     !
  Type(ZCoord_), Intent(In) :: ZCoordOut    !
  Real(Std),     Intent(In) :: GroundHeight ! Height of ground.
  Real(Std),     Intent(In) :: ZIn          !
  ! Function result:
  Real(Std) :: AboveSeaToAboveGround !

# ifdef ExtraChecks
    If (ZCoordIn%CoordType  /= Z_AboveSea .or.     &
        ZCoordOut%CoordType /= Z_AboveGround) Then
      Call Message('Error in AboveSeaToAboveGround', 4)
    End If
# endif

  AboveSeaToAboveGround = (ZIn*ZCoordIn%Unit - GroundHeight)/ZCoordOut%Unit

End Function AboveSeaToAboveGround

!-------------------------------------------------------------------------------------------------------------

Function AboveGroundToZAsEta(ZCoordIn, ZCoordOut, GroundHeight, ZIn)
! Converts coords between coord systems.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn     ! Input vertical coord system.
  Type(ZCoord_), Intent(In) :: ZCoordOut    ! Output vertical coord system.
  Real(Std),     Intent(In) :: GroundHeight ! Height of ground in metres.
  Real(Std),     Intent(In) :: ZIn          ! Input vertical coord.
  ! Function result:
  Real(Std) :: AboveGroundToZAsEta ! Height as eta.
  ! Locals:
  Real(Std) :: eta  ! Local copy of function result.
  Real(Std) :: z1   ! Height above ground level in metres.
  Real(Std) :: z2   ! Height above sea level in metres.
  Real(Std) :: z3   ! Expression used in z -> eta conversion   (= zI - 2h).

# ifdef ExtraChecks
    If (ZCoordIn%CoordType  /= Z_AboveGround  .or. &
        ZCoordOut%CoordType /= Z_ZAsEta       .or. &
        GroundHeight         > ZCoordOut%MaxOrogHeight) Then
      Call Message('Error in AboveGroundToZAsEta', 4)
    End If
# endif

  z1 = ZIn * ZCoordIn%Unit
  z2 = z1 + GroundHeight

  If (z1 < 0.0) Then

    eta = z1 / (ZCoordOut%ModelTopHeight - 2.0 * GroundHeight / ZCoordOut%InterfaceEta)

  Else If (z2 < ZCoordOut%InterfaceHeight) Then

    z3  = ZCoordOut%InterfaceHeight - 2.0 * GroundHeight
    eta = (2.0*ZCoordOut%InterfaceEta*z1) / (Sqrt(z3*z3 + 4.0*GroundHeight*z1) + z3)

  Else

    eta = z2 / ZCoordOut%ModelTopHeight

  End If

  AboveGroundToZAsEta = eta

End Function AboveGroundToZAsEta

!-------------------------------------------------------------------------------------------------------------

Function ZAsEtaToAboveGround(ZCoordIn, ZCoordOut, GroundHeight, ZIn)
! Converts coords between coord systems.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn     ! Input vertical coord system.
  Type(ZCoord_), Intent(In) :: ZCoordOut    ! Output vertical coord system.
  Real(Std),     Intent(In) :: GroundHeight ! Height of ground in metres.
  Real(Std),     Intent(In) :: ZIn          ! Input vertical coord.
  ! Function result:
  Real(Std) :: ZAsEtaToAboveGround ! Height above ground.
  ! Locals:
  Real(Std) :: z1   ! Height above ground level in metres.
  Real(Std) :: z2   ! Height above sea level in metres.
  Real(Std) :: r    ! Orography damping component (= 1 - eta/etaI).

# ifdef ExtraChecks
    If (ZCoordIn%CoordType  /= Z_ZAsEta       .or. &
        ZCoordOut%CoordType /= Z_AboveGround  .or. &
        GroundHeight         > ZCoordIn%MaxOrogHeight) Then
      Call Message('Error in ZAsEtaToAboveGround', 4)
    End If
# endif

  If (ZIn < 0.0) Then

    z1 = ZIn * (ZCoordIn%ModelTopHeight - 2.0 * GroundHeight / ZCoordIn%InterfaceEta)
    z2 = z1 + GroundHeight

  Else If (ZIn < ZCoordIn%InterfaceEta) Then

    r  = 1.0 - ZIn / ZCoordIn%InterfaceEta
    z2 = (ZIn * ZCoordIn%ModelTopHeight) + (r * r * GroundHeight)

  Else

    z2 = ZIn * ZCoordIn%ModelTopHeight

  End If

  z1 = z2 - GroundHeight

  ZAsEtaToAboveGround = z1 / ZCoordOut%Unit

End Function ZAsEtaToAboveGround

!-------------------------------------------------------------------------------------------------------------

Function AboveSeaToZAsEta(ZCoordIn, ZCoordOut, GroundHeight, ZIn)
! Converts coords between coord systems.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn     ! Input vertical coord system.
  Type(ZCoord_), Intent(In) :: ZCoordOut    ! Output vertical coord system.
  Real(Std),     Intent(In) :: GroundHeight ! Height of ground in metres.
  Real(Std),     Intent(In) :: ZIn          ! Input vertical coord.
  ! Function result:
  Real(Std) :: AboveSeaToZAsEta ! Height as eta.
  ! Locals:
  Real(Std) :: eta  ! Local copy of function result.
  Real(Std) :: z1   ! Height above ground level in metres.
  Real(Std) :: z2   ! Height above sea level in metres.
  Real(Std) :: z3   ! Expression used in z -> eta conversion   (= zI - 2h).

# ifdef ExtraChecks
    If (ZCoordIn%CoordType  /= Z_AboveSea     .or. &
        ZCoordOut%CoordType /= Z_ZAsEta       .or. &
        GroundHeight         > ZCoordOut%MaxOrogHeight) Then
      Call Message('Error in AboveSeaToZAsEta', 4)
    End If
# endif

  z2 = ZIn * ZCoordIn%Unit
  z1 = z2 - GroundHeight

  If (z1 < 0.0) Then

    eta = z1 / (ZCoordOut%ModelTopHeight - 2.0 * GroundHeight / ZCoordOut%InterfaceEta)

  Else If (z2 < ZCoordOut%InterfaceHeight) Then

    z3  = ZCoordOut%InterfaceHeight - 2.0 * GroundHeight
    eta = (2.0*ZCoordOut%InterfaceEta*z1) / (Sqrt(z3*z3 + 4.0*GroundHeight*z1) + z3)

  Else

    eta = z2 / ZCoordOut%ModelTopHeight

  End If

  AboveSeaToZAsEta = eta

End Function AboveSeaToZAsEta

!-------------------------------------------------------------------------------------------------------------

Function ZAsEtaToAboveSea(ZCoordIn, ZCoordOut, GroundHeight, ZIn)
! Converts coords between coord systems.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn     ! Input vertical coord system.
  Type(ZCoord_), Intent(In) :: ZCoordOut    ! Output vertical coord system.
  Real(Std),     Intent(In) :: GroundHeight ! Height of ground in metres.
  Real(Std),     Intent(In) :: ZIn          ! Input vertical coord.
  ! Function result:
  Real(Std) :: ZAsEtaToAboveSea ! Height above ground.
  ! Locals:
  Real(Std) :: z1   ! Height above ground level in metres.
  Real(Std) :: z2   ! Height above sea level in metres.
  Real(Std) :: r    ! Orography damping component (= 1 - eta/etaI).

# ifdef ExtraChecks
    If (ZCoordIn%CoordType  /= Z_ZAsEta       .or. &
        ZCoordOut%CoordType /= Z_AboveSea     .or. &
        GroundHeight         > ZCoordIn%MaxOrogHeight) Then
      Call Message('Error in ZAsEtaToAboveSea', 4)
    End If
# endif

  If (ZIn < 0.0) Then

    z1 = ZIn * (ZCoordIn%ModelTopHeight - 2.0 * GroundHeight / ZCoordIn%InterfaceEta)
    z2 = z1 + GroundHeight

  Else If (ZIn < ZCoordIn%InterfaceEta) Then

    r  = 1.0 - ZIn / ZCoordIn%InterfaceEta
    z2 = ZIn * ZCoordIn%ModelTopHeight + r * r * GroundHeight

  Else

    z2 = ZIn * ZCoordIn%ModelTopHeight

  End If

  z1 = z2 - GroundHeight ! Not used here - but retained for consistency with
                         ! other routines.

  ZAsEtaToAboveSea = z2 / ZCoordOut%Unit

End Function ZAsEtaToAboveSea

!-------------------------------------------------------------------------------------------------------------

Function AboveSeaToPAsZICAO(ZCoordIn, ZCoordOut, ZIn)
! Converts vertical coord from coord systems type 2 (height above sea) to type 4
! (flight level) assuming the ICAO standard atmosphere.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn  ! Coord system of input vertical coord.
  Type(ZCoord_), Intent(In) :: ZCoordOut ! Coord system of desired vertical coord.
  Real(Std),     Intent(In) :: ZIn       ! Vertical coord in coord system ZCoordIn.
  ! Function result:
  Real(Std) :: AboveSeaToPAsZICAO ! Output coord.

# ifdef ExtraChecks
    If (ZCoordIn%CoordType  /= Z_AboveSea .or. &
        ZCoordOut%CoordType /= Z_PAsZ) Then
      Call Message('Error in AboveSeaToPAsZICAO', 4)
    End If
# endif

  AboveSeaToPAsZICAO = (ZIn*ZCoordIn%Unit)/ZCoordOut%Unit

End Function AboveSeaToPAsZICAO

!-------------------------------------------------------------------------------------------------------------

Function PAsZToAboveSeaICAO(ZCoordIn, ZCoordOut, ZIn)
! Converts vertical coord from coord systems type 4 (flight level) to type 2 (height
! above sea) assuming the ICAO standard atmosphere.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn  ! Coord system of input vertical coord.
  Type(ZCoord_), Intent(In) :: ZCoordOut ! Coord system of desired vertical coord.
  Real(Std),     Intent(In) :: ZIn       ! Vertical coord in coord system ZCoordIn.
  ! Function result:
  Real(Std) :: PAsZToAboveSeaICAO ! Output coord.

# ifdef ExtraChecks
    If (ZCoordIn%CoordType  /= Z_PAsZ .or. &
        ZCoordOut%CoordType /= Z_AboveSea) Then
      Call Message('Error in PAsZToAboveSeaICAO', 4)
    End If
# endif

  PAsZToAboveSeaICAO = (ZIn*ZCoordIn%Unit)/ZCoordOut%Unit

End Function PAsZToAboveSeaICAO

!-------------------------------------------------------------------------------------------------------------

Function PToPAsZ(ZCoordIn, ZCoordOut, ZIn)
! Converts vertical coord from coord systems type 3 (pressure) to type 4 (flight
! level).

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn  ! Coord system of input vertical coord.
  Type(ZCoord_), Intent(In) :: ZCoordOut ! Coord system of desired vertical coord.
  Real(Std),     Intent(In) :: ZIn       ! Vertical coord in coord system ZCoordIn.
  ! Function result:
  Real(Std) :: PToPAsZ ! Output coord.
  ! Locals:
  Real(Std) :: P ! Pressure (Pa).
  Real(Std) :: Z ! Height above sea level as read by altimeter callibrated to the ICAO
                 ! standard atmosphere (m).
  Real(Std) :: T ! ICAO Temperature (K).
  Real(Std) :: Rho ! ICAO Density.

# ifdef ExtraChecks
    If (ZCoordIn%CoordType /= Z_P .or. ZCoordOut%CoordType /= Z_PAsZ) Then
      Call Message('Error in PToPAsZ', 4)
    End If
# endif

  P = ZIn * ZCoordIn%Unit

  Call ICAOAtP(P, Z, T, Rho)

  PToPAsZ = Z / ZCoordOut%Unit

End Function PToPAsZ

!-------------------------------------------------------------------------------------------------------------

Function PAsZToP(ZCoordIn, ZCoordOut, ZIn)
! Converts vertical coord from coord systems type 4 (flight level) to type 3
! (pressure).

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn  ! Coord system of input vertical coord.
  Type(ZCoord_), Intent(In) :: ZCoordOut ! Coord system of desired vertical coord.
  Real(Std),     Intent(In) :: ZIn       ! Vertical coord in coord system ZCoordIn.
  ! Function result:
  Real(Std) :: PAsZToP ! Output coord.
  ! Locals:
  Real(Std) :: P       ! Pressure at location of interest (Pa).
  Real(Std) :: Z       ! Height above sea level as read by altimeter callibrated to
                       ! the ICAO standard atmosphere.
  Real(Std) :: T ! ICAO Temperature (K).
  Real(Std) :: Rho ! ICAO Density.

# ifdef ExtraChecks
    If (ZCoordIn%CoordType /= Z_PAsZ .or. ZCoordOut%CoordType /= Z_P) Then
      Call Message('Error in PAsZToP', 4)
    End If
# endif

  Z = ZIn * ZCoordIn%Unit

  Call ICAOAtZ(Z, P, T, Rho)

  PAsZToP = P / ZCoordOut%Unit

End Function PAsZToP

!-------------------------------------------------------------------------------------------------------------

Function PToPAsEta(ZCoordIn, ZCoordOut, PStar, ZIn)
! Converts coords between coord systems.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn  !
  Type(ZCoord_), Intent(In) :: ZCoordOut !
  Real(Std),     Intent(In) :: PStar     ! Surface pressure (Pa).
  Real(Std),     Intent(In) :: ZIn       !
  ! Function result:
  Real(Std) :: PToPAsEta !
  ! Locals:
  Integer   :: i        !
  Integer   :: EtaLevel ! Number of eta level below ZIn.
  Real(Std) :: P        !
  Real(Std) :: PAbove   !
  Real(Std) :: PBelow   !

# ifdef ExtraChecks
    If (ZCoordIn%CoordType  /= Z_P .or. &
        ZCoordOut%CoordType /= Z_PAsEta) Then
      Call Message('Error in PToPAsEta', 4)
    End If
    If (PStar < ZCoordOut%PStarMin) Then
      Call Message(                                 &
             'Warning in PToPAsEta: PStar = '    // &
             Trim(Std2Char(PStar))               // &
             ' < PStarMin ='                     // &
             Trim(Std2Char(ZCoordOut%PStarMin)),    &
             1                                      &
           )
    End If
# endif

  P = ZIn*ZCoordIn%Unit

# ifdef ExtraChecks
    If (P < 0.0) Then
      Call Message('Error in PToPAsEta', 4)
    End If
# endif

  If (P < 0.0)   P = 0.0

  EtaLevel = ZCoordOut%nEtaLevels
  Do i = 1, ZCoordOut%nEtaLevels          ! Use bisection $$
    If (P >= ZCoordOut%A(i) + ZCoordOut%B(i)*PStar) Then
      EtaLevel = i - 1
      Exit
    End If
  End Do

  ! Interpolate.
  If (EtaLevel == ZCoordOut%nEtaLevels) Then
    PBelow = ZCoordOut%A(EtaLevel) + ZCoordOut%B(EtaLevel)*PStar
    PToPAsEta = ZCoordOut%Eta(EtaLevel) * P / PBelow
  Else If (EtaLevel == 0) Then
    PAbove = ZCoordOut%A(EtaLevel + 1) + ZCoordOut%B(EtaLevel + 1)*PStar
    PToPAsEta = ZCoordOut%Eta(EtaLevel + 1) * P / PAbove
  Else
    PBelow = ZCoordOut%A(EtaLevel    ) + ZCoordOut%B(EtaLevel    )*PStar
    PAbove = ZCoordOut%A(EtaLevel + 1) + ZCoordOut%B(EtaLevel + 1)*PStar
    PToPAsEta = ZCoordOut%Eta(EtaLevel) +                                &
                (ZCoordOut%Eta(EtaLevel + 1) - ZCoordOut%Eta(EtaLevel))* &
                (PBelow - P)/                                            &
                (PBelow - PAbove)
  End If

End Function PToPAsEta

!-------------------------------------------------------------------------------------------------------------

Function PAsEtaToP(ZCoordIn, ZCoordOut, PStar, ZIn)
! Converts coords between coord systems.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn  !
  Type(ZCoord_), Intent(In) :: ZCoordOut !
  Real(Std),     Intent(In) :: PStar     ! Surface pressure (Pa).
  Real(Std),     Intent(In) :: ZIn       !
  ! Function result:
  Real(Std) :: PAsEtaToP !
  ! Locals:
  Integer   :: i        !
  Integer   :: EtaLevel ! Number of eta level below ZIn.
  Real(Std) :: Eta      !
  Real(Std) :: PAbove   !
  Real(Std) :: PBelow   !

# ifdef ExtraChecks
    If (ZCoordIn%CoordType  /= Z_PAsEta .or. &
        ZCoordOut%CoordType /= Z_P) Then
      Call Message('Error in PAsEtaToP', 4)
    End If
    If (PStar < ZCoordIn%PStarMin) Then
      Call Message(                                &
             'Warning in PAsEtaToP: PStar = '   // &
             Trim(Std2Char(PStar))              // &
             ' < PStarMin ='                    // &
             Trim(Std2Char(ZCoordIn%PStarMin)),    &
             1                                     &
           )
    End If
# endif

  Eta = ZIn

  EtaLevel = ZCoordIn%nEtaLevels
  Do i = 1, ZCoordIn%nEtaLevels ! Use bisection $$
    If (Eta >= ZCoordIn%Eta(i)) Then
      EtaLevel = i - 1
      Exit
    End If
  End Do

  ! Interpolate.
  If (EtaLevel == ZCoordIn%nEtaLevels) Then
    PBelow = ZCoordIn%A(EtaLevel) + ZCoordIn%B(EtaLevel)*PStar
    PAsEtaToP = PBelow * Eta / ZCoordIn%Eta(EtaLevel)
  Else If (EtaLevel == 0) Then
    PAbove = ZCoordIn%A(EtaLevel + 1) + ZCoordIn%B(EtaLevel + 1)*PStar
    PAsEtaToP = PAbove * Eta / ZCoordIn%Eta(EtaLevel + 1)
  Else
    PBelow = ZCoordIn%A(EtaLevel    ) + ZCoordIn%B(EtaLevel    )*PStar
    PAbove = ZCoordIn%A(EtaLevel + 1) + ZCoordIn%B(EtaLevel + 1)*PStar
    PAsEtaToP = PBelow +                                              &
                (PAbove - PBelow)*                                    &
                (ZCoordIn%Eta(EtaLevel) - Eta)/                       &
                (ZCoordIn%Eta(EtaLevel) - ZCoordIn%Eta(EtaLevel + 1))
  End If

  PAsEtaToP = PAsEtaToP/ZCoordOut%Unit

End Function PAsEtaToP

!-------------------------------------------------------------------------------------------------------------

Function ZBasedToZBased(ZCoordIn, ZCoordOut, Topog, ZIn)
! Converts vertical coords between two height based coord systems.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn  ! Input vertical coordinate system.
  Type(ZCoord_), Intent(In) :: ZCoordOut ! Output vertical coordinate system.
  Real(Std),     Intent(In) :: Topog     ! Height of topography above sea level. Value
                                         ! irrelevant if the coord systems are of the
                                         ! same type and are not eta coord systems.
  Real(Std),     Intent(In) :: ZIn       ! Vertical coord in input coord system.
  ! Function result:
  Real(Std) :: ZBasedToZBased ! Output coord.
  ! Locals:
  Type(ZCoord_) :: ZCoord !

  Select Case (ZCoordIn%CoordType)

    Case (Z_AboveGround)

      Select Case (ZCoordOut%CoordType)

        Case (Z_AboveGround)

          ZBasedToZBased = ZSameToSame(ZCoordIn, ZCoordOut, ZIn)

        Case (Z_AboveSea)

          ZBasedToZBased = AboveGroundToAboveSea(ZCoordIn, ZCoordOut, Topog, ZIn)

        Case (Z_ZAsEta)

          ZBasedToZBased = AboveGroundToZAsEta(ZCoordIn, ZCoordOut, Topog, ZIn)

        Case Default

          Call Message('Error in ZBasedToZBased', 4)

      End Select

    Case (Z_AboveSea)

      Select Case (ZCoordOut%CoordType)

        Case (Z_AboveGround)

          ZBasedToZBased = AboveSeaToAboveGround(ZCoordIn, ZCoordOut, Topog, ZIn)

        Case (Z_AboveSea)

          ZBasedToZBased = ZSameToSame(ZCoordIn, ZCoordOut, ZIn)

        Case (Z_ZAsEta)

          ZBasedToZBased = AboveSeaToZAsEta(ZCoordIn, ZCoordOut, Topog, ZIn)

        Case Default

          Call Message('Error in ZBasedToZBased', 4)

      End Select

    Case (Z_ZAsEta)

      Select Case (ZCoordOut%CoordType)

        Case (Z_AboveGround)

          ZBasedToZBased = ZAsEtaToAboveGround(ZCoordIn, ZCoordOut, Topog, ZIn)

        Case (Z_AboveSea)

          ZBasedToZBased = ZAsEtaToAboveSea(ZCoordIn, ZCoordOut, Topog, ZIn)

        Case (Z_ZAsEta)

          ZCoord = ZCoord_m_agl()
          ZBasedToZBased = ZAsEtaToAboveGround(ZCoordIn, ZCoord,    Topog, ZIn)
          ZBasedToZBased = AboveGroundToZAsEta(ZCoord,   ZCoordOut, Topog, ZBasedToZBased)

        Case Default

          Call Message('Error in ZBasedToZBased', 4)

      End Select

    Case Default

      Call Message('Error in ZBasedToZBased', 4)

  End Select

End Function ZBasedToZBased

!-------------------------------------------------------------------------------------------------------------

Function PBasedToPBased(ZCoordIn, ZCoordOut, PS, ZIn)
! Converts vertical coords between two pressure based coord systems.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In) :: ZCoordIn  !
  Type(ZCoord_), Intent(In) :: ZCoordOut !
  Real(Std),     Intent(In) :: PS        ! Surface pressure. Value
                                         ! irrelevant for same type of coord system.
                                         ! or pressure and p as z.
  Real(Std),     Intent(In) :: ZIn       !
  ! Function result:
  Real(Std) :: PBasedToPBased !
  ! Local:
  Type(ZCoord_) :: ZCoord

  Select Case (ZCoordIn%CoordType)

    Case (Z_P)

      Select Case (ZCoordOut%CoordType)

        Case (Z_P)

          PBasedToPBased = ZSameToSame(ZCoordIn, ZCoordOut, ZIn)

        Case (Z_PAsZ)

          PBasedToPBased = PToPAsZ(ZCoordIn, ZCoordOut, ZIn)

        Case (Z_PAsEta)

          PBasedToPBased = PToPAsEta(ZCoordIn, ZCoordOut, PS, ZIn)

        Case Default

          Call Message('Error in PBasedToPBased', 4)

      End Select

    Case (Z_PAsZ)

      Select Case (ZCoordOut%CoordType)

        Case (Z_P)

          PBasedToPBased = PAsZToP(ZCoordIn, ZCoordOut, ZIn)

        Case (Z_PAsZ)

          PBasedToPBased = ZSameToSame(ZCoordIn, ZCoordOut, ZIn)

        Case (Z_PAsEta)

          ZCoord         = ZCoord_Pa()
          PBasedToPBased = PAsZToP(ZCoordIn, ZCoord, ZIn)
          PBasedToPBased = PToPAsEta(ZCoord, ZCoordOut, PS, PBasedToPBased)

        Case Default

          Call Message('Error in PBasedToPBased', 4)

      End Select

    Case (Z_PAsEta)

      Select Case (ZCoordOut%CoordType)

        Case (Z_P)

          PBasedToPBased = PAsEtaToP(ZCoordIn, ZCoordOut, PS, ZIn)

        Case (Z_PAsZ)

          ZCoord         = ZCoord_Pa()
          PBasedToPBased = PAsEtaToP(ZCoordIn, ZCoord, PS, ZIn)
          PBasedToPBased = PToPAsZ(ZCoord, ZCoordOut, PBasedToPBased)

        Case (Z_PAsEta)

          ZCoord         = ZCoord_Pa()
          PBasedToPBased = PAsEtaToP(ZCoordIn, ZCoord, PS, ZIn)
          PBasedToPBased = PToPAsEta(ZCoord, ZCoordOut, PS, PBasedToPBased)

        Case Default

          Call Message('Error in PBasedToPBased', 4)

      End Select

    Case Default

      Call Message('Error in PBasedToPBased', 4)

  End Select

End Function PBasedToPBased

!-------------------------------------------------------------------------------------------------------------

Function InitCoords()
! Initialises a collection of coord systems.

  Implicit None
  ! Function result:
  Type(Coords_) :: InitCoords !
  ! Locals:
  Type(Coords_) :: Coords !

  Coords%nHCoords = 0
  Coords%nZCoords = 0

  InitCoords = Coords

End Function InitCoords

!-------------------------------------------------------------------------------------------------------------

Subroutine AddHCoord(HCoord, Coords)
! Adds a horizontal coord system to a collection of coord systems.

  Implicit None
  ! Argument list:
  Type(HCoord_), Intent(In)    :: HCoord !
  Type(Coords_), Intent(InOut) :: Coords !
  ! Locals:
  Integer :: i !

  Do i = 1, Coords%nHCoords
    If (HCoord%Name .CIEq. Coords%HCoords(i)%Name) Then
      If (HCoord == Coords%HCoords(i)) Then
        Return
      Else
        Call Message(                                               &
               'Error in AddHCoord: '                            // &
               'two different coord systems with the same name',    &
               3                                                    &
             )
      End If
    End If
  End Do
  If (Coords%nHCoords >= MaxHCoords) Then
    Call Message('Error in AddHCoord', 3)
  End If
  Coords%nHCoords                 = Coords%nHCoords + 1
  Coords%HCoords(Coords%nHCoords) = HCoord

End Subroutine AddHCoord

!-------------------------------------------------------------------------------------------------------------

Function FindHCoordIndex(Name, Coords)
! .

  Implicit None
  ! Argument list:
  Character(*),  Intent(In) :: Name
  Type(Coords_), Intent(In) :: Coords
  ! Function result:
  Integer :: FindHCoordIndex !
  ! Locals:
  Integer :: i

  Do i = 1, Coords%nHCoords
    If (Name .CIEq. Coords%HCoords(i)%Name) Then
      FindHCoordIndex = i
      Return
    End If
  End Do

  Call Message(                                           &
         'FATAL ERROR: horizontal coordinate system "' // &
         Trim(Name)                                    // &
         '" not found',                                   &
         3                                                &
       )

End Function FindHCoordIndex

!-------------------------------------------------------------------------------------------------------------

Subroutine AddZCoord(ZCoord, Coords)
! Adds a vertical coord system to a collection of coord systems.

  Implicit None
  ! Argument list:
  Type(ZCoord_), Intent(In)    :: ZCoord !
  Type(Coords_), Intent(InOut) :: Coords !
  ! Locals:
  Integer :: i !

  Do i = 1, Coords%nZCoords
    If (ZCoord%Name .CIEq. Coords%ZCoords(i)%Name) Then
      If (ZCoord == Coords%ZCoords(i)) Then
        Return
      Else
        Call Message(                                               &
               'Error in AddZCoord: '                            // &
               'two different coord systems with the same name',    &
               3                                                    &
             )
      End If
    End If
  End Do
  If (Coords%nZCoords >= MaxZCoords) Then
    Call Message('Error in AddZCoord', 3)
  End If
  Coords%nZCoords                 = Coords%nZCoords + 1
  Coords%ZCoords(Coords%nZCoords) = ZCoord

End Subroutine AddZCoord

!-------------------------------------------------------------------------------------------------------------

Function FindZCoordIndex(Name, Coords)
! .

  Implicit None
  ! Argument list:
  Character(*),  Intent(In) :: Name
  Type(Coords_), Intent(In) :: Coords
  ! Function result:
  Integer :: FindZCoordIndex !
  ! Locals:
  Integer :: i

  Do i = 1, Coords%nZCoords
    If (Name .CIEq. Coords%ZCoords(i)%Name) Then
      FindZCoordIndex = i
      Return
    End If
  End Do

  Call Message(                                         &
         'FATAL ERROR: vertical coordinate system "' // &
         Trim(Name)                                  // &
         '" not found',                                 &
         3                                              &
       )

End Function FindZCoordIndex

!-------------------------------------------------------------------------------------------------------------

Function X2Position(Coords, X, iHCoord, iZCoord) Result(Position)
! .

  Implicit None
  ! Argument list:
  Type(Coords_), Intent(In)  :: Coords
  Real(Std),     Intent(In)  :: X(3)
  Integer,       Intent(In)  :: iHCoord
  Integer,       Intent(In)  :: iZCoord
  ! Function result:
  Type(Position_) :: Position

# ifdef ExtraChecks
    If (iHCoord > Coords%nHCoords .or. &
        iHCoord < 1               .or. &
        iZCoord > Coords%nZCoords .or. &
        iZCoord < 1) Then
      Call Message('Error in X2Position', 4)
    End If
# endif

  Position%XYValid(1:Coords%nHCoords) = .false.
  Position%ZValid (1:Coords%nZCoords) = .false.
  Position%TopogValid                 = .false.
  Position%PSValid                    = .false.
  Position%RhoValid                   = .false.
  Position%XY     (:,iHCoord)         = X(1:2)
  Position%Z      (  iZCoord)         = X(3)
  Position%XYValid(  iHCoord)         = .true.
  Position%ZValid (  iZCoord)         = .true.

End Function X2Position

!-------------------------------------------------------------------------------------------------------------

Function Position2X(Coords, Position, iHCoord, iZCoord) Result(X)
! .

  Implicit None
  ! Argument list:
  Type(Coords_),   Intent(In)  :: Coords
  Type(Position_), Intent(In)  :: Position
  Integer,         Intent(In)  :: iHCoord
  Integer,         Intent(In)  :: iZCoord
  ! Function result:
  Real(Std) :: X(3)

# ifdef ExtraChecks
    If (iHCoord > Coords%nHCoords .or. &
        iHCoord < 1               .or. &
        iZCoord > Coords%nZCoords .or. &
        iZCoord < 1) Then
      Call Message('Error in Position2X', 4)
    End If
    If (.not.Position%XYValid(iHCoord) .or. &
        .not.Position%ZValid (iZCoord)) Then
      Call Message('Error in Position2X', 4)
    End If
# endif

  X(1:2) = Position%XY(:,iHCoord)
  X(3)   = Position%Z (  iZCoord)

End Function Position2X

!-------------------------------------------------------------------------------------------------------------

Subroutine Position2XUnknownCoord(Coords, Position, iHCoord, iZCoord, X)
! .

  Implicit None
  ! Argument list:
  Type(Coords_),   Intent(In)  :: Coords
  Type(Position_), Intent(In)  :: Position
  Integer,         Intent(Out) :: iHCoord
  Integer,         Intent(Out) :: iZCoord
  Real(Std),       Intent(Out) :: X(3)
  ! Locals:
  Integer :: i

  iHCoord = 0
  iZCoord = 0

  Do i = 1, Coords%nHCoords
    If (Position%XYValid(i)) Then
      iHCoord = i
      X(1:2)  = Position%XY(:, iHCoord)
      Exit
    End If
  End Do

  Do i = 1, Coords%nZCoords
    If (Position%ZValid(i)) Then
      iZCoord = i
      X(3)    = Position%Z(iZCoord)
      Exit
    End If
  End Do

  If (iHCoord == 0 .or. iZCoord == 0) Then
    Call Message('Error in Position2XUnknownCoord', 4)
  End If

End Subroutine Position2XUnknownCoord

!-------------------------------------------------------------------------------------------------------------

Subroutine ConvertToH(Coords, iHCoord, Position)
! Converts horizontal coords to a particular horizontal coord system using information
! on the value of the coords in a number of horizontal coord systems.

  Implicit None
  ! Argument list:
  Type(Coords_),   Intent(In)    :: Coords   !
  Integer,         Intent(In)    :: iHCoord  !
  Type(Position_), Intent(InOut) :: Position !
  ! Locals:
  Integer :: iHCoordUse !
  Integer :: i          !

# ifdef ExtraChecks
    If (                        &
      iHCoord < 1 .or.          &
      Coords%nHCoords < iHCoord &
    ) Then
      Call Message('Error in ConvertToH', 4)
    End If
# endif

  If (.not.Position%XYValid(iHCoord)) Then

    Do i = 1, Coords%nHCoords ! $$ replace by efficient selection
      If (Position%XYValid(i)) Then
        iHCoordUse = i
        Exit
      End If
    End Do

    Position%XY(:,iHCoord) = ConvertH(Coords%HCoords(iHCoordUse), &
                             Coords%HCoords(iHCoord),             &
                             Position%XY(:,iHCoordUse))
    Position%XYValid(iHCoord) = .true.

  End If

End Subroutine ConvertToH

!-------------------------------------------------------------------------------------------------------------

Subroutine RationaliseHPosition(XY, HCoord)
! If particle position is beyond the edge of the coordinate system, this
! subroutine will place it back into the domain.
! Also, if X is beyond edge of wrapped domain, places X back in domain.

  Implicit None
  ! Argument list:
  Real(Std),     Intent(InOut) :: XY(2)  ! X = XY(1), Y = XY(2).
  Type(HCoord_), Intent(In)    :: HCoord

  ! Ensure Y is within range of the coord system, correcting X if necessary.
  If (HCoord%YEdge) Then
    If (HCoord%CoordType == H_LatLong) Then
      If (XY(2) > HCoord%YEdgeMax) Then  ! $$ need do-while if a long way outside range.
        XY(2) = 2.0 * HCoord%YEdgeMax - XY(2)
        XY(1) = XY(1) + HCoord%XCycle / 2.0
      Else If (XY(2) < HCoord%YEdgeMin) Then
        XY(2) = 2.0 * HCoord%YEdgeMin - XY(2)
        XY(1) = XY(1) + HCoord%XCycle / 2.0
      End If
    End If
  End If

  ! Ensure X is between +/- XCycle.
  If (HCoord%XCycle /= 0.0) Then
    Do While (XY(1) > Abs(HCoord%XCycle) / 2.0)
      XY(1) = XY(1) - Abs(HCoord%XCycle)
    End Do
    Do While (XY(1) < - Abs(HCoord%XCycle) / 2.0)
      XY(1) = XY(1) + Abs(HCoord%XCycle)
    End Do
  End If

! $$ needs tidying up for other coord systems than lat-long

End Subroutine RationaliseHPosition

!-------------------------------------------------------------------------------------------------------------

End Module CoordinateSystemModule
