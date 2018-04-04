! Module:  Physical Units Module

Module PhysicalUnitsModule

! This module defines a framework for handling physical units.
!
! The Type MaterialUnit_ is used to describe the unit a species/source is measured in.
! The conversion factor describes how the numerical value of a quantity changes when it is converted between
! different units. Consider, for example, mass. The reference unit is 'g' (grams) with 
! a conversion factor of 1.0. The conversion factor of the unit 'mg' (milligrams) is 1e-3, i.e. 
! 1 mg contains 1.0e-3 g and to convert a quantity measured in mg to the unit g, it has to be multiplied
! by 1.0e-3. To convert from mg to kg we would need to multiply the quantity by 1.0e-3/(1.0e3) = 1.0e-6.
! Only units with the same UnitType can be converted into each other.
!
! The supported units are stored in the global variable MaterialUnits, the index of a unit with a given 
! name can be extracted with the function FindMaterialUnitIndex.
! MaterialUnits is initialied in the subroutine MaterialUnitsInitialise.

Use StringModule
Use GlobalParametersModule
Use ErrorAndMessageModule

Private

Public  :: MaterialUnitsInitialise       ! Initialise material units
Public  :: AddMaterialUnit               ! Add material unit to the collection
Public  :: GetMaterialUnit               ! Return the material unit with a given
Public  :: FindMaterialUnitIndex
Public  :: MaterialUnit_
Public  :: MaterialUnits_
Public  :: MassUnitType
Public  :: VolumetricUnitType
Public  :: ActivityUnitType
Public  :: DobsonUnitType
Public  :: UnknownUnitType               ! This type is assigned to all unknown units
Public  :: MaterialUnitEq
Public  :: MaterialUnitNotEq
Public  :: Operator(.eq.)                ! Comparison operators for Material Units
Public  :: Operator(.ne.)

!-------------------------------------------------------------------------------------------------------------

Type :: MaterialUnit_
  Character(MaxCharLength) :: Name              ! name of unit, e.g. 'kg'
  Character(MaxCharLength) :: Description       ! } longer description, e.g. 'kilogram', which can be
                                                ! } used in NAME input file.
  Integer                  :: UnitType          ! } Type of unit 1: mass unit, 2: activity unit. Unknown
                                                ! } units which are added at runtime should have the negative
                                                ! } type UnknownUnitType
  Real                     :: Conversionfactor  ! } Described the number of time the reference unit 
                                                ! } fits into this unit (see comment on top).
End Type MaterialUnit_

! Collection of material units

Integer, Parameter :: MaxMaterialUnits = 100    ! Maximal number of material units

Type :: MaterialUnits_
  Integer :: nMaterialUnits                               ! Current number of material units
  Type(MaterialUnit_) :: MaterialUnits(MaxMaterialUnits)  ! List with material units
End Type

! Different predefined unit types
Integer, Parameter :: MassUnitType       = 1
Integer, Parameter :: VolumetricUnitType = 2
Integer, Parameter :: ActivityUnitType   = 3
Integer, Parameter :: DobsonUnitType     = 4
Integer, Parameter :: UnknownUnitType    = -9

!-------------------------------------------------------------------------------------------------------------

Interface Operator(.eq.) ! Equality of Material Units
  Module Procedure MaterialUnitEq
End Interface

Interface Operator(.ne.) ! Equality of Material Units
  Module Procedure MaterialUnitNotEq
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function MaterialUnitEq(MaterialUnit1, MaterialUnit2)

  Implicit None

  Logical                         :: MaterialUnitEq
  Type(MaterialUnit_), Intent(In) :: MaterialUnit1
  Type(MaterialUnit_), Intent(In) :: MaterialUnit2

  MaterialUnitEq = ( Trim(MaterialUnit1%Name)        .eq. Trim(MaterialUnit2%Name)                ) .and. &
                   ( Trim(MaterialUnit1%Description)  ==  Trim(MaterialUnit2%Description)         ) .and. &
                   ( MaterialUnit1%UnitType           ==  MaterialUnit2%UnitType                  ) .and. &
                   ( Abs(MaterialUnit1%ConversionFactor - MaterialUnit2%ConversionFactor) < 1e-12 )
  
End Function MaterialUnitEq

!-------------------------------------------------------------------------------------------------------------

Function MaterialUnitNotEq(MaterialUnit1, MaterialUnit2)

  Implicit None

  Logical                         :: MaterialUnitNotEq
  Type(MaterialUnit_), Intent(In) :: MaterialUnit1
  Type(MaterialUnit_), Intent(In) :: MaterialUnit2

  MaterialUnitNotEq = .not. MaterialUnitEq(MaterialUnit1,MaterialUnit2)
  
End Function MaterialUnitNotEq

!-------------------------------------------------------------------------------------------------------------

SubRoutine MaterialUnitsInitialise(MaterialUnits)

! initialise structure MaterialUnits and add standard units.
   
  Implicit None

  Type(MaterialUnits_), Intent(InOut) :: MaterialUnits
  Real :: p_std = 1.013e5   ! standard pressure [Pa]
  Real :: T_std = 273.0     ! standard temperature [K]
  Real :: Mr_Ozone = 47.998 ! molar mass of ozone [g mol^{-1}]
  Real :: R_ideal = 8.3145  ! ideal gas constant [J K^{-1} mol^{-1}]

! Mass units. Reference = 1 g (gram)
  Call AddMaterialUnit('t'   ,'tonne'    , MassUnitType, 1.0e6       , MaterialUnits )
  Call AddMaterialUnit('kg'  ,'kilogram' , MassUnitType, 1.0e3       , MaterialUnits )
  Call AddMaterialUnit('g'   ,'gram'     , MassUnitType, 1.0         , MaterialUnits )
  Call AddMaterialUnit('mg'  ,'milligram', MassUnitType, 1.0e-3      , MaterialUnits )
  Call AddMaterialUnit('mcg' ,'microgram', MassUnitType, 1.0e-6      , MaterialUnits )
  Call AddMaterialUnit('lb'  ,'pound'    , MassUnitType, 453.59237   , MaterialUnits )
  Call AddMaterialUnit('oz'  ,'ounce'    , MassUnitType, 28.349523125, MaterialUnits )

! Volumetric units. Reference = 1 ppm (parts per million)
  Call AddMaterialUnit('ppm', 'parts per million', VolumetricUnitType, 1.0  , MaterialUnits )
  Call AddMaterialUnit('ppb', 'parts per billion', VolumetricUnitType, 1.0e3, MaterialUnits )
  
! Activity units. Reference = 1 Bq (Becquerel)
  Call AddMaterialUnit('Bq'   , 'Becquerel'     , ActivityUnitType, 1.0   , MaterialUnits )
  Call AddMaterialUnit('mBq'  , 'milliBecquerel', ActivityUnitType, 1.0e-3, MaterialUnits )
  Call AddMaterialUnit('mcBq' , 'microBecquerel', ActivityUnitType, 1.0e-6, MaterialUnits )
  Call AddMaterialUnit('Ci'   , 'Curie'         , ActivityUnitType, 3.7e10, MaterialUnits )

! Dobson units. The conversion factor is chosen for the case that the species is ozone and 
!               is given in g (note that in the subroutine ProcessFields() [Output.F90] the conversion
!               factor is given as Species%ConversionFactor/FieldReq%ConversionFactor).
  Call AddMaterialUnit('DU'   , 'Dobson units' , DobsonUnitType, 1.0e-5*p_std*Mr_Ozone/(R_ideal*T_std), &
                       MaterialUnits                                                                    &
       )
  
End SubRoutine MaterialUnitsInitialise

!-------------------------------------------------------------------------------------------------------------

Subroutine AddMaterialUnit(Name,Description,UnitType,ConversionFactor,MaterialUnits)

! Add a unit to the structure MaterialUnits.
   
  Implicit None

  Character(*),         Intent(In)    :: Name             !} name of unit, e.g. 'kg', which can be used in 
                                                          !} input file
  Character(*),         Intent(In)    :: Description      ! more detailled description, e.g. 'kilogram'
  Integer,              Intent(In)    :: UnitType         ! Type of unit 1: mass unit, 2: activity unit 
  Real,                 Intent(In)    :: Conversionfactor ! Conversion factor relative to a standard unit
  Type(MaterialUnits_), Intent(InOut) :: MaterialUnits    ! set of Material Units
  Integer                  :: i,n

  If (MaterialUnits%nMaterialUnits == MaxMaterialUnits) Then
    Call Message('ERROR: Can not have more than '//Int2Char(MaxMaterialUnits)//' material units',4)
  End If

  Do i = 1, MaterialUnits%nMaterialUnits
    If (Trim(MaterialUnits%MaterialUnits(i)%Name) .EQ. Trim(Name)) Then
      Call Message('ERROR: Material unit with name '//Trim(Name)//' already exists.',4)
    End If
  End Do
  
  MaterialUnits%nMaterialUnits = MaterialUnits%nMaterialUnits + 1
  n = MaterialUnits%nMaterialUnits
  MaterialUnits%MaterialUnits(n)%Name             = Trim(Name)
  MaterialUnits%MaterialUnits(n)%Description      = Trim(Description)
  MaterialUnits%MaterialUnits(n)%UnitType         = UnitType
  MaterialUnits%MaterialUnits(n)%ConversionFactor = ConversionFactor
  
End Subroutine AddMaterialUnit

!-------------------------------------------------------------------------------------------------------------

Function FindMaterialUnitIndex(Name, MaterialUnits)

! Return the index in the MaterialUnits structure of unit with a specific name.

  Implicit None

  Integer                          :: FindMaterialUnitIndex
  Character(*),         Intent(In) :: Name
  Type(MaterialUnits_), Intent(In) :: MaterialUnits
  Integer                          :: i
  Integer                          :: Index

  Index = -1
  Do i = 1, MaterialUnits%nMaterialUnits
    If (Trim(MaterialUnits%MaterialUnits(i)%Name) .EQ. Trim(Name)) Index = i
  End Do

  FindMaterialUnitIndex = Index
   
End Function FindMaterialUnitIndex

!-------------------------------------------------------------------------------------------------------------

End Module PhysicalUnitsModule

