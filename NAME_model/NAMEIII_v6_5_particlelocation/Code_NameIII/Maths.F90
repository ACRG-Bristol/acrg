! Module:  Maths Module

Module MathsModule

! This module contains mathematical constants and functions.

!-------------------------------------------------------------------------------------------------------------

Use GlobalParametersModule
Use ErrorAndMessageModule
Use StringModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: Pi                    ! Pi.
Public :: Pi_Pos                ! Pi.
Public :: RandomState_          ! A state of the random number generator.
Public :: Erf                   ! Calculates the error function using an approximation given by Abramowitz and
                                ! Stegun, Handbook of mathematical functions, Dover Publications (1965).
Public :: ATan2ZeroTest         ! Two-argument form of the ArcTan function with zero returned if both
                                ! arguments
                                ! are zero.
Public :: GaussA                ! Generates a Gaussian random number with mean zero and variance 1 by
                                ! adding 12
                                ! independent uniform random numbers.
Public :: GaussB                ! Generates a Gaussian random number with mean zero and variance 1 by
                                ! transforming
                                ! two uniform random numbers into two Gaussian random numbers
                                ! following the method
                                ! in Numerical Recipes, 2nd ed., page 279.
Public :: SetRandomState        ! Sets the state of the random number generator.
Public :: WriteRandomState      ! Writes the state of the random number generator.
Public :: ReadRandomState       ! Reads the state of the random number generator.
Public :: GetRandomNumber       ! Wrapper subroutine for random number generator.
Public :: InitialiseRandomSeeds ! Initialises the array of random seeds.

!-------------------------------------------------------------------------------------------------------------

Real(Std), Parameter :: Pi     = 3.141592653589793238462643_Std ! Pi.
Real(Pos), Parameter :: Pi_Pos = 3.141592653589793238462643_Pos ! Pi.

!-------------------------------------------------------------------------------------------------------------

Type :: RandomState_ ! A state of the random number generator.
  ! If changing this type, remember WriteRandomState and ReadRandomState routines.
  Integer,  Pointer :: IntrinsicSeed(:) => Null()
  Logical           :: Set
  Real(Std)         :: Gauss2
  Character(1)      :: Mode
  Integer,  Pointer :: Seeds(:, :) => Null()
  ! IntrinsicSeed :: The seed of the intrinsic random number generator Random_Number.
  ! Set           :} Used in the routine GaussB which generates two Gaussian random variables from two uniform
  ! Gauss2        :} random variables. Set indicates if a second Gaussian variable has been generated but not
  !                  used and Gauss2 is its value.
  ! Mode          :: Random or fixed seeds, complete array or last 2 used
  ! Seeds         :: Seeds for random number generator.

! $$ Preinit routine to set to null?
End Type

!-------------------------------------------------------------------------------------------------------------

Type(RandomState_), Public, Save :: GlobalRandomState

!-------------------------------------------------------------------------------------------------------------

Integer(I64), Parameter :: Multipliers(2) = (/40014, 40692/)
                                      ! Multipliers for random number generator.
Integer(I64), Parameter :: Moduli(2)      = (/2147483563, 2147483399/)
                                      ! Moduli for random number generator.
Integer(I64), Parameter :: PowerOf2       = 34
                                      ! Spacing between random seeds of different particles/puffs
                                      ! Allows up to 2**27 particles/puffs

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Function Erf(X)
! Calculates the error function using an approximation given by Abramowitz and Stegun, Handbook of
! mathematical functions, Dover Publications (1965).

  Implicit None
  ! Argument list:
  Real(Std), Intent(In) :: X ! Argument of error function.
  ! Function result:
  Real(Std) :: Erf ! Value of error function at X.
  ! Local parameters:
  Real(Std), Parameter :: P  =  0.3275911_Std   !} Constants used in numerical approximation to the error
  Real(Std), Parameter :: A1 =  0.254829592_Std !} function.
  Real(Std), Parameter :: A2 = -0.284496736_Std !}
  Real(Std), Parameter :: A3 =  1.421413741_Std !}
  Real(Std), Parameter :: A4 = -1.453152027_Std !}
  Real(Std), Parameter :: A5 =  1.061405429_Std !}
  ! Locals:
  Real(Std) :: Temp ! Temporary variable.

  Temp = 1.0 / (1.0 + P * Abs(X))
  Erf = 1.0 - Exp( - X*X) * (A1*Temp + A2*Temp**2 + A3*Temp**3 + A4*Temp**4 + A5*Temp**5)
  If (X < 0.0) Erf = -Erf

End Function Erf

!-------------------------------------------------------------------------------------------------------------

Function ATan2ZeroTest(Y, X)
! Two-argument form of the ArcTan function with zero returned if both arguments are zero.

  Implicit None
  ! Argument list:
  Real(Pos), Intent(In) :: Y !} Components of the vector whose angle is required.
  Real(Pos), Intent(In) :: X !}
  ! Function result:
  Real(Pos) :: ATan2ZeroTest ! Angle of vector (x,y) to the x-axis, with anticlockwise values positive. Lies
                             ! in [-Pi, Pi].

  If (X == 0.0 .and. Y == 0.0) Then
    ATan2ZeroTest = 0.0
  Else
    ATan2ZeroTest = ATan2(Y, X)
  End If

End Function ATan2ZeroTest

!-------------------------------------------------------------------------------------------------------------

Function GaussA(N)
! Generates a Gaussian random number with mean zero and variance 1 by adding independent uniform random
! numbers.

  Implicit None
  ! Argument list:
  Integer :: N ! Index of random seeds.
  ! Function result:
  Real(Std) :: GaussA ! Gaussian Random number.
  ! Locals:
  Integer   :: i    ! Loop index.
  Real(Std) :: Temp ! Temporary variable.

  GaussA = - 6.0
  Do i = 1, 12
    Call GetRandomNumber(Temp, N)
    GaussA = GaussA + Temp
  End Do

End Function GaussA

!-------------------------------------------------------------------------------------------------------------

Function GaussB(N)
! Generates a Gaussian random number with mean zero and variance 1 by transforming two uniform random numbers
! into two Gaussian random numbers folowing the method in Numerical Recipes, 2nd ed., page 279.

  Implicit None
  ! Argument list:
  Integer :: N ! Index of random seeds.
  ! Function result:
  Real(Std) :: GaussB ! Gaussian Random number.
  ! Locals:
  Real(Std) :: Temp(2) ! Uniform random variables on [-1, 1].
  Real(Std) :: RSq     ! | (Temp(1), Temp(2)) |^2.
  Real(Std) :: Factor  ! Transformation factor.

  If (.not. GlobalRandomState%Set) Then

    Do

      Call GetRandomNumber(Temp(1), N)
      Call GetRandomNumber(Temp(2), N)

      Temp(1) = 2.0 * Temp(1) - 1.0
      Temp(2) = 2.0 * Temp(2) - 1.0
      RSq     = Temp(1)**2 + Temp(2)**2

      If (RSq < 1.0 .and. RSq /= 0.0) Exit

    End Do

    Factor                   = Sqrt(- 2.0 * Log(RSq) / RSq)
    GaussB                   = Temp(1)*Factor
    GlobalRandomState%Gauss2 = Temp(2)*Factor
    If (.not. Associated(GlobalRandomState%Seeds)) GlobalRandomState%Set = .true.

  Else

    GaussB                = GlobalRandomState%Gauss2
    GlobalRandomState%Set = .false.

  End If

End Function GaussB

!-------------------------------------------------------------------------------------------------------------

Subroutine SetRandomState(Mode, RandomState, Seed)
! Sets the state of the random number generator.

  Implicit None
  ! Argument list:
  Character(1),       Intent(In)            :: Mode
  Type(RandomState_), Intent(Out)           :: RandomState
  Integer,            Intent(In),  Optional :: Seed(:)
  ! Mode        :: Mode of initialisation of the random number generator. 'F', 'R' and 'G' indicate that the
  !                random number generator's state is to be initialised with fixed, random and given values
  !                respectively. For option 'G', Seed must be given. In the latter case
  !                the remaining elements of the state of the random number generator are set to default
  !                values.
  ! RandomState :: The state of the random number generator.
  ! Seed        :: The given value for the seed of the intrinsic random number generator Random_Number.
  ! Locals:
  Integer :: SeedSize ! Size of the seed of the intrinsic random number generator Random_Number.
  Integer :: Error    ! Error code.
  Integer :: i

  RandomState%Mode = Mode

  ! Checks on consistency of input.
  If (Mode == 'F') Then
    If (Present(Seed)) Then
      Call Message('UNEXPECTED ERROR in SetRandomState', 4)
    End If
  Else If (Mode == 'R') Then
    If (Present(Seed)) Then
      Call Message('UNEXPECTED ERROR in SetRandomState', 4)
    End If
  Else If (Mode == 'G') Then
    If (.not.Present(Seed)) Then
      Call Message('UNEXPECTED ERROR in SetRandomState', 4)
    End If
  Else If (Mode == 'P') Then
    Return
  Else If (Mode == 'I') Then
    Return
  Else If (Mode == 'A') Then
    Return
  Else
    Call Message('UNEXPECTED ERROR in SetRandomState', 4)
  End If

  If (Mode == 'F' .or. Mode == 'R' .or. Mode == 'G') Then
    Call Random_Seed(Size = SeedSize)
    Allocate(RandomState%IntrinsicSeed(SeedSize), Stat = Error)
    If (Error /= 0) Call Message('UNEXPECTED FATAL ERROR: Unable to allocate array for random number seed', 4)
  End If

  ! Set random seed to fixed value.
  If (Mode == 'F') Then

    Call Random_Seed(Size = SeedSize)

    ! Compaq Win compiler. These values are the compiler's defaults for starting the random number generator
    ! if no call to Random_Seed is made. Starting from (1, 1) seems to generate the sequence obtained by
    ! subtracting the default sequence from 1.
#   ifdef CompaqPCCompiler
      If (SeedSize /= 2) Call Message('UNEXPECTED ERROR in SetRandomState', 4)
      RandomState%IntrinsicSeed(:) = (/ 2147483562, 2147483398 /)
      Call Random_Seed(Put = RandomState%IntrinsicSeed)
#   endif

    ! Intel Linux compiler. These values are the compiler's defaults for starting the random number generator
    ! if no call to Random_Seed is made.
#   ifdef IntelLinCompiler
      If (SeedSize /= 2) Call Message('UNEXPECTED ERROR in SetRandomState', 4)
      RandomState%IntrinsicSeed(:) = (/ 2147483562, 2147483398 /)
      Call Random_Seed(Put = RandomState%IntrinsicSeed)
#   endif

    RandomState%Set    = .false.
    RandomState%Gauss2 = 0.0

  ! Set random seed to random value.
  Else If (Mode == 'R') Then

    Call Random_Seed()

    Call Random_Seed(Get = RandomState%IntrinsicSeed)
    Call Message(' ')
    Call Message('The random number seed or seed array for this run is as follows:')
    Do i = 1, Size(RandomState%IntrinsicSeed)
      Call Message(Int2Char(RandomState%IntrinsicSeed(i)))
    End Do

    RandomState%Set    = .false.
    RandomState%Gauss2 = 0.0

  ! Set random seed to given value.
  Else If (Mode == 'G') Then

    Call Random_Seed(Size = SeedSize)

    If (Size(Seed) /= SeedSize) Then
      Call Message(                                                               &
             'ERROR in setting the state of the random number generator: the ' // &
             'seed array size is not correct for the random rumber generator',    &
             3                                                                    &
           )
    End If

    Call Random_Seed(Put = Seed)

    RandomState%Set    = .false.
    RandomState%Gauss2 = 0.0

  End If

! $$
! Something like the following may be needed for random seed on different compilers
!  Integer :: DateTime(8)
!  Allocate(Seed(SeedSize))
!
!  ! Random. Taken from current time.
!  If (Mode = 'R') Then
!
!    Call Date_and_Time(Values = DateTime)
!
!    Seed = DateTime(8) + DateTime(7)*1000 + DateTime(6)*60*1000 &
!           + DateTime(5)*60*60*1000 + (DateTime(3) + DateTime(2) + DateTime(1))*1000
!    Call Random_Seed(Put = Seed)  ! Sets random number seed
!
!  End If

End Subroutine SetRandomState

!-------------------------------------------------------------------------------------------------------------

Subroutine WriteRandomState(Unit, RandomState)
! Writes the state of the random number generator.

  Implicit None
  ! Argument List:
  Integer,            Intent(In)    :: Unit        ! Input/output unit number.
  Type(RandomState_), Intent(InOut) :: RandomState ! The state of the random number generator. This is InOut
                                                   ! because it is updated using the part of the random state
                                                   ! hidden inside the intrinsic random number generator.

  Write (Unit) RandomState%Mode

  Write (Unit) Associated(RandomState%IntrinsicSeed)
  If (Associated(RandomState%IntrinsicSeed)) Then
    Call Random_Seed(Get = RandomState%IntrinsicSeed)
    Write (Unit) Shape(RandomState%IntrinsicSeed)
    Write (Unit) RandomState%IntrinsicSeed(:)
  End If

  Write (Unit) RandomState%Set

  Write (Unit) RandomState%Gauss2

  Write (Unit) Associated(RandomState%Seeds)
  If (Associated(RandomState%Seeds)) Then
    Write (Unit) Shape(RandomState%Seeds)
    Write (Unit) RandomState%Seeds(:,:)
  End If

End Subroutine WriteRandomState

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadRandomState(Unit, RandomState)
! Reads the state of the random number generator.

  Implicit None
  ! Argument List:
  Integer,            Intent(In)    :: Unit        ! Input/output unit number.
  Type(RandomState_), Intent(InOut) :: RandomState ! The state of the random number generator. This is InOut
                                                   ! because it carries the allocation information in.
  ! Locals:
  Logical :: IsAssociated  ! Indicates the allocatable array should be allocated.
  Integer :: ArrayShape(2) ! Shape of array to be allocated.
  Integer :: Error         ! Error code.

  Read (Unit, End = 1, Err = 2) RandomState%Mode

  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(1)
    If (.not.Associated(RandomState%IntrinsicSeed)) Then ! $$ needed?
      Allocate(RandomState%IntrinsicSeed(ArrayShape(1)), Stat = Error)
      If (Error /= 0) Call Message('FATAL ERROR: Unable to allocate array for random number seed', 3)
    Else If (Size(RandomState%IntrinsicSeed) /= ArrayShape(1)) Then
      Call Message('UNEXPECTED ERROR in ReadRandomState', 4)
    End If
    Read (Unit, End = 1, Err = 2) RandomState%IntrinsicSeed(:)
    Call Random_Seed(Put = RandomState%IntrinsicSeed)
  End If

  Read (Unit, End = 1, Err = 2) RandomState%Set

  Read (Unit, End = 1, Err = 2) RandomState%Gauss2

  Read (Unit, End = 1, Err = 2) IsAssociated
  If (IsAssociated) Then
    Read (Unit, End = 1, Err = 2) ArrayShape(:)
    If (.not.Associated(RandomState%Seeds)) Then ! $$ needed?
      Allocate(RandomState%Seeds(ArrayShape(1), 0:ArrayShape(2) - 1), Stat = Error)
      If (Error /= 0) Call Message('FATAL ERROR: Unable to allocate array for random number seeds', 3)
    End If
    Read (Unit, End = 1, Err = 2) RandomState%Seeds(:, :)
  End If

  Return

1 Continue
  Call Message('FATAL ERROR: The restart file is shorter than expected', 3)

2 Continue
  Call Message('FATAL ERROR: An error occurred when the reading restart file', 3)

End Subroutine ReadRandomState

!-------------------------------------------------------------------------------------------------------------

Subroutine GetRandomNumber(Ran, N)
! Wrapper subroutine to call random number generator.
! Random Number Generator.  Returns a uniform random number between 0 and 1 from 2 linear
! congruential random number generators.  Moduli and Multipliers global variables in maths module.

! $$

  Implicit None
  ! Argument List:
  Real(Std), Intent(Out)          :: Ran        ! Random number to be returned.
  Integer,   Intent(In)           :: N          ! Index of particle or puff
  ! Locals:
  Integer      :: Seed(2)
  Integer(I64) :: Seed64Bit(2)

  If (N < 0) Call Message('Particle or Puff index fault', 4) ! $$

  If (Associated(GlobalRandomState%Seeds)) Then

    Seed(1) = GlobalRandomState%Seeds(1, N)
    Seed(2) = GlobalRandomState%Seeds(2, N)

    Seed64Bit = Seed
    Seed(1) = Mod(Multipliers(1) * Seed64Bit(1), Moduli(1))
    Seed(2) = Mod(Multipliers(2) * Seed64Bit(2), Moduli(2))
    Ran = Seed(1) - Seed(2)
    Do While (Ran < 0.5)
      Ran = Ran + (Real(Moduli(1)) - 1.0)
    End Do
    Do While (Ran > Moduli(1) - 0.5)
      Ran = Ran - (Real(Moduli(1)) - 1.0)
    End Do
    Ran = Ran / Real(Moduli(1))

    GlobalRandomState%Seeds(1, N) = Seed(1)
    GlobalRandomState%Seeds(2, N) = Seed(2)

  Else

    Call Random_Number(Ran)

  End If

End Subroutine GetRandomNumber

!-------------------------------------------------------------------------------------------------------------

Subroutine InitialiseRandomSeeds(Max, Mode)
! Initialises the array of random seeds.  Has four mode options:
! P - Program defined array, based on Intel compiler standard initial random seeds
! I - Initial seeds taken from standard input
! S - Saved array from restart file
! A - Array based on system time defined seeds, obtained using Random_Seed()

  Implicit None
  !Argument List:
  Integer,   Intent(In)           :: Max     ! Maximum number of particles/puffs in run.
  Character, Intent(In)           :: Mode    ! Fixed or random seed
  ! Locals:
  Integer      :: GetSeed(2)
  Integer      :: i
  Integer(I64) :: Seed64Bit(2) ! 64-bit integers to prevent overflow on 32-bit machines
  Integer(I64) :: Multi(2)     ! Multipliers to separate inital seeds for each particle/puff
  Integer :: Error             ! Error code.

  If (Mode == 'P') Then

    GetSeed = Moduli - 1

  Else If (Mode == 'I') Then

    Write(*,*) 'Input 1st random seed'
    Read(*,*) GetSeed(1)
    If (GetSeed(1) == 0) GetSeed(1) = Moduli(1) - 1
    Write(*,*) 'Input 2nd random seed'
    Read(*,*) GetSeed(2)
    If (GetSeed(2) == 0) GetSeed(2) = Moduli(2) - 1

  Else If (Mode == 'A') Then

    Call Random_Seed()
    Call Random_Seed(GET = GetSeed)
    Write(LogFileUnit,*) 'First pair of random seeds is ',GetSeed

  Else

    Return

  End If

  If (.not.Associated(GlobalRandomState%Seeds)) Then
    Allocate (GlobalRandomState%Seeds(1:2,0:Max), Stat = Error)
    If (Error /= 0) Then
      Call Message('FATAL ERROR: Unable to allocate array for random number seeds', 3)
    End If
    GlobalRandomState%Set = .false.
    Call Message('Random seed array allocated')
  Else
    Return
  End If

  GlobalRandomState%Seeds(1,0) = GetSeed(1)
  GlobalRandomState%Seeds(2,0) = GetSeed(2)

  Multi = Multipliers
  Do i=1,PowerOf2
    Multi(1) = Mod(Multi(1)*Multi(1),Moduli(1))
    Multi(2) = Mod(Multi(2)*Multi(2),Moduli(2))
  End Do

  Do i=1,Max
    Seed64Bit(1) = GlobalRandomState%Seeds(1,i-1)
    Seed64Bit(2) = GlobalRandomState%Seeds(2,i-1)
    GlobalRandomState%Seeds(1,i) = Mod(Seed64Bit(1)*Multi(1),Moduli(1))
    GlobalRandomState%Seeds(2,i) = Mod(Seed64Bit(2)*Multi(2),Moduli(2))
  End Do

End Subroutine InitialiseRandomSeeds

!-------------------------------------------------------------------------------------------------------------

End Module MathsModule
