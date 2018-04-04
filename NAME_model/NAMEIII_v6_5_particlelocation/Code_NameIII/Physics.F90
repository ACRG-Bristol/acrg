! Module: Physics Module

Module PhysicsModule

! This module contains physical constants and functions.

!-------------------------------------------------------------------------------------------------------------

Use GlobalParametersModule
Use ErrorAndMessageModule
Use StringModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private

! General physical constants:
Public :: Gravity       ! Acceleration due to gravity.
Public :: Avogadro      ! Avogadro's number (molecules per mole).
Public :: GasConstant   ! Gas constant for dry air.
Public :: Cp            ! Specific heat capacity for dry air.
Public :: EarthRadius   ! Earth radius.
Public :: VK            ! Von Karman's constant.
Public :: MoleMassAir   ! Molecular mass of dry air (g/mole).
Public :: MoleMassWater ! Molecular mass of water (g/mole).
Public :: TKAtTCEq0     ! Thermodynamic temperature at 0 degrees Celsius (K).
Public :: c_pv          ! Specific heat capacity of water vapour.
Public :: c_pw          ! Specific heat capacity of liquid water.
Public :: c_ps          ! Specific heat capacity of volcanic ash (Mastin 2007).  
Public :: RhoWater      ! Density of liquid water.
Public :: RVapour       ! Gas constant for water vapour.

! Shape Scheme
Public :: ShapeScheme_Ganser
Public :: ShapeScheme_WH
Public :: ShapeScheme_White

! Constants for the ICAO standard atmosphere:
Public :: TAt0km            ! Temperature at 0km above mean sea level (K).
Public :: TAt11km           ! Temperature at 11km above mean sea level (K).
Public :: TAt20km           ! Temperature at 20km above mean sea level (K).
Public :: LapseRate0To11km  ! Lapse rate from 0 to 11km above mean sea level.
Public :: LapseRate20kmPlus ! Lapse rate at more than 20km above mean sea level.
Public :: PAt0km            ! Pressure at mean sea level (Pa).

! Reference values used in definitions:
Public :: PRef ! Reference pressure for potential temperature (Pa).

! Non-SI units (in SI units):
Public :: Foot ! One foot.

! Properties of the ICAO standard atmosphere:
Public :: ICAOAtP ! Calculates properties of the ICAO standard atmosphere at a given pressure. This routine
                  ! works below ground and assumes the tropospheric lapse rate continues below ground.
Public :: ICAOAtZ ! Calculates properties of the ICAO standard atmosphere at a given height. This routine
                  ! works below ground and assumes the tropospheric lapse rate continues below ground.

! Humidity functions:
Public :: CalcQ           ! Calculates specific humidity from relative humidity, temperature and pressure.
Public :: CalcRH          ! Calculates relative humidity from specific humidity, temperature and pressure.
Public :: CalcSatVapP     ! Calculates saturation vapour pressure from temperature using Wexler's formula
                          ! (1976, Journal of Research of the National Bureau of Standards - A Physics and
                          ! Chemistry, 80A, 775-785).
Public :: CalcMixingRatio ! Calculates humidity mixing ratio from pressure and vapour pressure.
Public :: CalcLatentHeat  ! Calculates latent heat of vaporization of water.

! Pasquill stability:
Public :: CalcPasquill ! Calculates Pasquill stability using the continuous equations described in TDN 206 but
                       ! with Nielsen et al's modified cloud amount replaced by cloud amount.

! Particle terminal velocity.
Public :: TerminalVelocity ! Calculates particle terminal velocity.

!-------------------------------------------------------------------------------------------------------------

! Physics parameters (in SI units):
Real(Std), Parameter :: Gravity       =    9.80665_Std ! Acceleration due to gravity.
Real(Std), Parameter :: Avogadro      = 6.02214E23_Std ! Avogadro's number (molecules per mole).
Real(Std), Parameter :: GasConstant   =     287.05_Std ! Gas constant for dry air.
Real(Std), Parameter :: Cp            =     1004.6_Std ! Specific heat capacity for dry air.
Real(Std), Parameter :: EarthRadius   =  6371229.0_Std ! Earth radius.
Real(Std), Parameter :: VK            =        0.4_Std ! Von Karman's constant.
Real(Std), Parameter :: MoleMassAir   =     28.966_Std ! Molecular mass of dry air (g/mole).
Real(Std), Parameter :: MoleMassWater =     18.016_Std ! Molecular mass of water (g/mole).
Real(Std), Parameter :: TKAtTCEq0     =     273.15_Std ! Thermodynamic temperature at 0 degrees Celsius (K).
Real(Std), Parameter :: c_pv          =     1859.0_Std ! Specific heat capacity of water vapour
Real(Std), Parameter :: c_pw          =     4183.0_Std ! Specific heat capacity of liquid water
Real(Std), Parameter :: c_ps          =     1100.0_Std ! Specific heat capacity of volcanic ash (Mastin 2007)  
Real(Std), Parameter :: RhoWater      =     1000.0_Std ! Density of liquid water
Real(Std), Parameter :: RVapour       =      462.0_Std ! Gas constant for water vapour

! Parameters for the ICAO standard atmosphere:
Real(Std), Parameter :: TAt0km            =   288.15_Std ! Temperature at 0km above mean sea level (K).
Real(Std), Parameter :: TAt11km           =   216.65_Std ! Temperature at 11km above mean sea level (K).
Real(Std), Parameter :: TAt20km           =   216.65_Std ! Temperature at 20km above mean sea level (K).
Real(Std), Parameter :: LapseRate0To11km  =   0.0065_Std ! Lapse rate from 0 to 11km above mean sea level.
Real(Std), Parameter :: LapseRate20kmPlus =   -0.001_Std ! Lapse rate at more than 20km above mean sea level.
Real(Std), Parameter :: PAt0km            = 101325.0_Std ! Pressure at mean sea level (Pa).

! Definition of particle Shape Scheme:
Integer,   Parameter :: ShapeScheme_Ganser = 1 ! }
Integer,   Parameter :: ShapeScheme_WH     = 2 ! } Associate shape scheme with an integer
Integer,   Parameter :: ShapeScheme_White  = 3 ! }

! Reference values used in definitions (in SI units):
Real(Std), Parameter :: PRef = 100000.0 ! Reference pressure for potential temperature.

! Non-SI units (in SI units):
Real(Std), Parameter :: Foot = 0.3048_Std ! One foot.

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine ICAOAtP(P, Z, T, Rho)
! Calculates properties of the ICAO standard atmosphere at a given pressure. This routine works below ground
! and assumes the tropospheric lapse rate continues below ground.

  Implicit None
  ! Argument list:
  Real(Std), Intent(In)  :: P   ! Pressure (Pa).
  Real(Std), Intent(Out) :: Z   ! Height above sea level (m).
  Real(Std), Intent(Out) :: T   ! Temperature (K).
  Real(Std), Intent(Out) :: Rho ! Density.
  ! Locals:
  Real(Std) :: PAt11km ! Pressure at height 11km (Pa).
  Real(Std) :: PAt20km ! Pressure at height 20km (Pa).

  PAt11km = PAt0km *                                                                                    &
            (1.0 - LapseRate0To11km * 11000.0 / TAt0km) ** (Gravity / (GasConstant * LapseRate0To11km))

  PAt20km = PAt11km * Exp( - Gravity * 9000.0 / (GasConstant * TAt11km))

  If (P > PAt11km) Then

    Z = (1.0 - (P / PAt0km) ** (GasConstant * LapseRate0To11km / Gravity)) * TAt0km / LapseRate0To11km
    T = TAt0km - LapseRate0To11km * Z

  Else If (P > PAt20km) Then

    Z = 11000.0 + Log(PAt11km / P) * (GasConstant * TAt11km) / Gravity
    T = TAt11km

  Else

    Z = 20000.0 +                                                                                          &
        (1.0 - (P / PAt20km) ** (GasConstant * LapseRate20kmPlus / Gravity)) * TAt20km / LapseRate20kmPlus
    T = TAt20km - LapseRate20kmPlus * (Z - 20000.0)

  End If

  Rho = P / (GasConstant * T)

End Subroutine ICAOAtP

!-------------------------------------------------------------------------------------------------------------

Subroutine ICAOAtZ(Z, P, T, Rho)
! Calculates properties of the ICAO standard atmosphere at a given height. This routine works below ground and
! assumes the tropospheric lapse rate continues below ground.

  Implicit None
  ! Argument list:
  Real(Std), Intent(In)  :: Z   ! Height above sea level (m).
  Real(Std), Intent(Out) :: P   ! Pressure (Pa).
  Real(Std), Intent(Out) :: T   ! Temperature (K).
  Real(Std), Intent(Out) :: Rho ! Density.
  ! Locals:
  Real(Std) :: PAt11km ! Pressure at height 11km (Pa).
  Real(Std) :: PAt20km ! Pressure at height 20km (Pa).

  PAt11km = PAt0km *                                                                                    &
            (1.0 - LapseRate0To11km * 11000.0 / TAt0km) ** (Gravity / (GasConstant * LapseRate0To11km))

  PAt20km = PAt11km * Exp( - Gravity * 9000.0 / (GasConstant * TAt11km))

  If (Z < 11000.0) Then

    P = PAt0km *                                                                              &
        (1.0 - LapseRate0To11km * Z / TAt0km) ** (Gravity / (GasConstant * LapseRate0To11km))
    T = TAt0km - LapseRate0To11km * Z

  Else If (Z < 20000.0) Then

    P = PAt11km * Exp( - Gravity * (Z - 11000.0) / (GasConstant * TAt11km))
    T = TAt11km

  Else

    P = PAt20km *                                                                                            &
        (1.0 - LapseRate20kmPlus * (Z - 20000.0) / TAt20km) ** (Gravity / (GasConstant * LapseRate20kmPlus))
    T = TAt20km - LapseRate20kmPlus * (Z - 20000.0)

  End If

  Rho = P / (GasConstant * T)

End Subroutine ICAOAtZ

!-------------------------------------------------------------------------------------------------------------

Function CalcQ(RH, T, P) Result(Q)
! Calculates specific humidity from relative humidity, temperature and pressure.

  Implicit None
  ! Argument list:
  Real(Std), Intent(In) :: RH ! Relative humidity (%).
  Real(Std), Intent(In) :: T  ! Temperature (K).
  Real(Std), Intent(In) :: P  ! Pressure (Pa).
  ! Function result:
  Real(Std) :: Q ! Specific humidity.
  ! Locals:
  Real(Std) :: SatVP ! Saturated vapour pressure (Pa).
  Real(Std) :: SatMR ! Saturated mixing ratio.
  Real(Std) :: MR    ! Mixing ratio.

  SatVP = CalcSatVapP(T)
  SatMR = CalcMixingRatio(P, SatVP)
  MR    = SatMR * RH / 100.0
  Q     = MR / (1.0 + MR)

End Function CalcQ

!-------------------------------------------------------------------------------------------------------------

Function CalcRH(Q, T, P) Result(RH)
! Calculates relative humidity from specific humidity, temperature and pressure.

  Implicit None
  ! Argument list:
  Real(Std), Intent(In) :: Q ! Specific humidity.
  Real(Std), Intent(In) :: T ! Temperature (K).
  Real(Std), Intent(In) :: P ! Pressure (Pa).
  ! Function result:
  Real(Std) :: RH ! Relative humidity (%).
  ! Locals:
  Real(Std) :: SatVP ! Saturated vapour pressure (Pa).
  Real(Std) :: SatMR ! Saturated mixing ratio.
  Real(Std) :: MR    ! Mixing ratio.

  SatVP  = CalcSatVapP(T)
  SatMR  = CalcMixingRatio(P, SatVP)
  MR     = Q / (1.0 - Q)
  RH     = Min(MR * 100.0 / SatMR, 100.0)

End Function CalcRH

!-------------------------------------------------------------------------------------------------------------

Function CalcSatVapP(T) Result(SatVapP)
! Calculates saturation vapour pressure from temperature using Wexler's formula (1976, Journal of Research of
! the National Bureau of Standards - A Physics and Chemistry, 80A, 775-785).

  Implicit None
  ! Argument list:
  Real(Std), Intent(In) :: T ! Temperature (K).
  ! Function result:
  Real(Std) :: SatVapP ! Saturation vapour pressure (Pa).
  ! Local parameters:
  Real(Std), Parameter :: G1 = -2.9912729E3   !} Constants.
  Real(Std), Parameter :: G2 = -6.0170128E3   !}
  Real(Std), Parameter :: G3 = 18.87643854    !}
  Real(Std), Parameter :: G4 = -2.8354721E-2  !}
  Real(Std), Parameter :: G5 =  1.7838301E-5  !}
  Real(Std), Parameter :: G6 = -8.4150417E-10 !}
  Real(Std), Parameter :: G7 =  4.4412543E-13 !}
  Real(Std), Parameter :: G8 =  2.858487      !}
  ! Locals:
  Real(Std) :: Temp ! Temporary variable.

# ifdef ExtraChecks
    ! Check for T <= 0. !$$ perhaps should be unexpected error, if check T when read in.
    If (T <= 0.0) Then
      Call Message('FATAL ERROR in CalcSatVapP: T = ' // Trim(Std2Char(T)) // ' <= 0.', 3)
    End If
# endif

  Temp = G1 * T**(- 2) + &
         G2 * T**(- 1) + &
         G3 * T**0     + &
         G4 * T**1     + &
         G5 * T**2     + &
         G6 * T**3     + &
         G7 * T**4     + &
         G8 * ALog(T)
  SatVapP = Exp(Temp)

End Function CalcSatVapP

!-------------------------------------------------------------------------------------------------------------

Function CalcMixingRatio(P, VapP) Result(MixingRatio)
! Calculates mixing ratio from pressure and vapour pressure.

  Implicit None
  ! Argument list:
  Real(Std), Intent(In) :: P    ! Pressure.
  Real(Std), Intent(In) :: VapP ! Vapour pressure.
  ! Function result:
  Real(Std) :: MixingRatio ! Mixing ratio.

# ifdef ExtraChecks
    If (VapP >= P) Then
      Call Message(                                                  &
             'UNEXPECTED FATAL ERROR in CalcMixingRatio. VapP = ' // &
             Trim(Std2Char(VapP))                                 // &
             ' >= P = '                                           // &
             Trim(Std2Char(P)),                                      &
             4                                                       &
           )
    End If
# endif

  MixingRatio = (MoleMassWater/MoleMassAir)*VapP/(P - VapP)

End Function CalcMixingRatio

!-------------------------------------------------------------------------------------------------------------

Function CalcLatentHeat(T0) Result(LatentHeat)
! Calculates latent heat of vaporization of water.

  Implicit None
  ! Argument list:
  Real(Std), Intent(In) :: T0 ! Temperature (K).
  ! Function result:
  Real(Std) :: LatentHeat ! Latent heat of vaporization of water.

  LatentHeat = 2500800.0_Std - 2300.0 * (T0 - TKAtTCEq0)

End Function CalcLatentHeat

!-------------------------------------------------------------------------------------------------------------

Function CalcPasquill(WindSpeed, HeatFlux, Cloud) Result(P)
! Calculates Pasquill stability using the continuous equations described in TDN 206 but with Nielsen et al's
! modified cloud amount replaced by cloud amount.

  Implicit None
  ! Argument list:
  Real(Std), Intent(In) :: WindSpeed ! Wind speed.
  Real(Std), Intent(In) :: HeatFlux  ! Surface sensible heat flux.
  Real(Std), Intent(In) :: Cloud     ! Cloud amount (fraction).
  ! Function result:
  Real(Std) :: P ! Pasquill stability (continuously varying rather than a category, as in Smith's extensions
                 ! of Pasquill's approach).
  ! Locals:
  Real(Std) :: UHat       ! Min(WindSpeed, 8).
  Real(Std) :: CloudOktas ! Cloud amount (oktas).

  UHat = Min(WindSpeed, 8.0)

  If (HeatFlux > 0.0) Then

    P = 7.0 -                                                                        &
        (2.26 + 0.019 * (UHat - 5.6)**2) *                                           &
        (0.1 * HeatFlux + 2.0 + 0.4 * UHat**1.5) ** (0.28 - 0.004 * (UHat - 2.0)**2)
    If (P < 0.0) P = 0.0

  Else

    CloudOktas = Cloud * 8.0
    P = 3.6 + Exp(- 3.0 * UHat / 8.0) * (120.0 - 13.3 * CloudOktas) / (27.0 - 2.0 * CloudOktas)
    If (P > 7.0) P = 7.0

  End If

End Function CalcPasquill

!-------------------------------------------------------------------------------------------------------------

Function TerminalVelocity(Diameter, Density, ParticleShape, ShapeSchemeCode, T, P, Rho) Result(WSed)
! Calculates particle terminal velocity.

  Implicit None
  ! Argument list:
  Real(Std), Intent(In) :: Diameter        ! Particle diameter (um).
  Real(Std), Intent(In) :: Density         ! Particle density.
  Real(Std), Intent(In) :: ParticleShape   ! Particle shape.
  Real(Std), Intent(In) :: T               ! Temperature (K).
  Real(Std), Intent(In) :: P               ! Pressure.
  Real(Std), Intent(In) :: Rho             ! Air density.
  Integer,   Intent(In) :: ShapeSchemeCode ! Shape Scheme.
  ! Function result:
  Real(Std) :: WSed ! Particle terminal velocity.
  ! Locals:
  Real(Std) :: DynViscosity   ! Dynamic viscosity.
  Real(Std) :: ReynoldsNumber ! Reynolds number.
  Real(Std) :: LambdaA        ! Mean free path.
  Real(Std) :: CCF            ! Cunningham correction factor.
  Real(Std) :: CD             ! Drag coefficient.
  Real(Std) :: K1             ! For Ganser drag coefficient.
  Real(Std) :: K2             ! For Ganser drag coefficient.
  Real(Std) :: Drag100        ! For Wilson + Huang drag coefficient. 
  Integer   :: i              ! Loop index.  

  If (Diameter > 0.0) Then

    If (T >= 273.15) Then
      DynViscosity = (1.718 + 0.0049 * (T - 273.15)) * 1.0E-5
    Else
      DynViscosity = (1.718 + 0.0049 * (T - 273.15) -                 &
                      1.2E-5 * (T - 273.15) * (T - 273.15)) * 1.0E-5
    End If
     
    ! Stokes regime (Re < 1).
    WSed = Diameter * 1.0E-6 * Diameter * 1.0E-6 * Gravity * (Density - Rho) &
             / (18.0 * DynViscosity)

    ! Iteration using Stokes value as initial value.
    Do i = 1, 10
    
      ReynoldsNumber = WSed * Diameter * 1.0E-6 * Rho / DynViscosity
      
      ! Choose scheme to calculate drag coefficient
      ! If non-spherical particles, then either the Ganser (1993)
      ! or the Wilson and Huang (1979) scheme must be selected. 
      ! If Spherical particles the White (1974) scheme is used (TDN 244)
      
      ! Ganser.
      If  (ShapeSchemeCode == ShapeScheme_Ganser) Then
        
        K1 = 3.0 / (1.0 + 2.0 * ParticleShape ** (-0.5))
        
        K2 = 10.0 ** (1.84148 * ((-1.0 * log10(ParticleShape)) ** 0.5743))
        
        CD = ( (24.0 / (ReynoldsNumber * K1)) *                                &
               (1.0 + 0.1118 * ((ReynoldsNumber * K1 * K2) ** 0.6567)) )  +    &
             ( (0.4305 * K2) / (1.0 + (3305.0 / (ReynoldsNumber * K1 * K2))) )
      
      ! Wilson and Huang.
      Else If  (ShapeSchemeCode == ShapeScheme_WH) Then
      
        Drag100 = ((24.0 / 100.0) * ParticleShape ** (-0.828) ) + &
                   (2.0 * Sqrt(1.07 - ParticleShape))

        If (ReynoldsNumber > 100 .AND. ReynoldsNumber < 1000) Then
                          
           CD = ( (1.0 - Drag100) / (1000.0 - 100.0) ) * (ReynoldsNumber - 100.0) + Drag100

        Else If (ReynoldsNumber <= 100) Then

           CD = ( (24.0 / ReynoldsNumber) * ParticleShape ** (-0.828) ) + &
                (2.0 * Sqrt(1.07 - ParticleShape))  
                              
        Else 
                   
           CD = 1

        End If
      
      ! White (Original Scheme).
      Else If (ShapeSchemeCode == ShapeScheme_White) Then
      
         CD  = 0.25 + 24.0 / ReynoldsNumber + 6.0 / (1.0 + Sqrt(ReynoldsNumber))
            
      End If 

      WSed = Sqrt(4.0 * Diameter * 1.0E-6 * Gravity * (Density - Rho) / (3.0 * CD * Rho)) 
   
    End Do

    ! Apply slip flow 'Cunningham' Correction Factor - important for small particles
    LambdaA = 6.6E-8 * (101325.0 / P) * (T / 293.15)
    CCF     = 1.0 + 2.0 * LambdaA / (Diameter * 1.0E-6) *                         &
               (1.257 + 0.4 * Exp(- 1.1 * Diameter * 1.0E-6 / (2.0 * LambdaA)))
    WSed    = CCF * WSed

  Else

    WSed = 0.0

  End If

End Function TerminalVelocity

!-------------------------------------------------------------------------------------------------------------

End Module PhysicsModule
