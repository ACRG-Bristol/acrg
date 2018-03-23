! Module:  OpenMP Module

Module OpenMPModule

! This module provides subroutines for synchronising the IO thread and the worker threads
! during the parallel reading of MetData.
! 
! The synchronisation is realised with a set of locks (array ParallelReadLocks) and subroutines for 
! setting and unsetting them: 
!
!   LookaheadFileReadRequestWait()  <--- [waits for ] --- LookaheadFileReadRequest()
!   LookaheadFileReadCompleteWait() <--- [waits for ] --- LookaheadFileReadComplete()
! 
! The module also contains routines for initialising the parallel IO configuration, the setup is controlled
! by the options in the structure Type(OpenMPOpts_). These are read from the input file at runtime.
!
!------------------------------------------------------------------------------------------------------------

Use TimeModule, Only : ShortTime_

!$ Use omp_lib

Use StringModule, Only : Int2Char
Use ErrorAndMessageModule

!------------------------------------------------------------------------------------------------------------

Implicit None

!------------------------------------------------------------------------------------------------------------

Private

! public subroutines & functions
Public ::    PreInitOpenMPOpts             &
            ,InitOpenMPOpts                &
            ,CheckOpenMPOpts               &
            ,FinaliseOpenMPConfiguration   &
            ,InitIOSynchronisation         &
            ,LookaheadFileReadRequest      &
            ,LookaheadFileReadRequestWait  &
            ,LookaheadFileReadComplete     &
            ,LookaheadFileReadCompleteWait &
            ,GetNumberOfParallelIOThreads  &
            ,OpenMPMaxThreads
            
! public variables
Public ::    IOLookaheadComplete           &    ! } Prefetch termination test. Set to .false. if there is no 
                                                ! } more data to prefetch
            ,IOMaxLookahead                &    ! } Maximal lookahead time. This is set in MainName.F90
                                                ! } $$ Could be refined
            ,WorkerThreadID                &    ! } id of worker and IO thread
            ,IOThreadID                         ! }

! public data types
Public :: OpenMPOpts_


! Data type for OpenMP options
Type :: OpenMPOpts_
  Logical        :: Initialised                 ! Have the OpenMP options already been initialised?
  Logical        :: UseOpenMP                   ! Do we use OpenMP at all? 
  Integer        :: nParticleThreads            ! number of threads in the particle loop 
  Integer        :: nParticleUpdateThreads      ! number of threads in the particle update loop
  Integer        :: nChemistryThreads           ! number of threads in the chemistry loop
  Integer        :: nOutputGroupThreads         ! number of threads in the output group loop 
  Integer        :: nOutputProcessThreads       ! number of threads in the output processing loop 
  Logical        :: ParallelMetRead             ! read metdata with separate thread?
  Logical        :: ParallelMetProcess          ! process metdata with separate IO thread?
End Type OpenMPOpts_

Type(ShortTime_)    :: IOMaxLookahead           ! Maximal lookahead time. This is calculated in MainName.F90
                                                ! $$ Could be refined

Integer             :: NumberOfParallelIOThreads=1

! Locks for thread synchronisation

!$ Integer(kind=OMP_lock_kind) :: ParallelReadLocks(3)

! Index of current lock, private to each thread
Integer             :: ParallelReadLockIndex
!$OMP THREADPRIVATE(ParallelReadLockIndex)

Integer, Parameter  :: WorkerThreadID = 0
Integer, Parameter  :: IOThreadID     = 1
Logical             :: IOLookaheadComplete=.false.

Integer             :: NRequests = 0

! The following variables are used for reading environment variables:
Character(len=64)   :: VALUE
Integer             :: LENGTH
Integer             :: STATUS

! Maximal allowed number of OpenMP threads
Integer, Parameter :: OpenMPMaxThreads = 32

!------------------------------------------------------------------------------------------------------------

Contains

!------------------------------------------------------------------------------------------------------------

Function PreInitOpenMPOpts()
! Pre-Initialise OpenMP configuration.

  Implicit None

  Type(OpenMPOpts_) :: PreInitOpenMPOpts

  PreInitOpenMPOpts%Initialised = .true.

  ! Set to default values
  PreInitOpenMPOpts%UseOpenMP              = .false.
  PreInitOpenMPOpts%nParticleThreads       = 1
  PreInitOpenMPOpts%nParticleUpdateThreads = 1
  PreInitOpenMPOpts%nChemistryThreads      = 1
  PreInitOpenMPOpts%nOutputGroupThreads    = 1
  PreInitOpenMPOpts%nOutputProcessThreads  = 1
  PreInitOpenMPOpts%ParallelMetRead        = .false.
  PreInitOpenMPOpts%ParallelMetprocess     = .false.
  
End Function PreInitOpenMPOpts

!------------------------------------------------------------------------------------------------------------


Subroutine InitOpenMPOpts(OpenMPOpts)
! Initialise OpenMP configuration.
!
! Print out parallel model and number of threads. Options are passed via the structure OpenMPOpts, which is 
! of type(OpenMPOpts_).

  Implicit None

  Type(OpenMPOpts_), Intent(InOut) :: OpenMPOpts

  ! Number of threads
  Integer :: nThreads

  ! Loop variable
  Integer :: i

  OpenMPOpts%Initialised = .true.
  
  nThreads = -1

  !$OMP PARALLEL
  !$ nThreads = omp_get_num_threads()
  !$OMP END PARALLEL
  
  If (nThreads .EQ. -1) Then
    If (OpenMPOpts%UseOpenMP) Then
      Call Message('WARNING: Running in serial mode, OpenMP has been disabled',1)
      OpenMPOpts%UseOpenMP=.false.
      OpenMPOpts%nParticleThreads=1
      OpenMPOpts%nParticleUpdateThreads=1
      OpenMPOpts%nChemistryThreads=1
      OpenMPOpts%nOutputGroupThreads=1
      OpenMPOpts%nOutputProcessThreads=1
      OpenMPOpts%ParallelMetRead=.false.
      OpenMPOpts%ParallelMetProcess=.false.
    End If
  Else
    !$ Call omp_set_nested(.true.)
    If (OpenMPOpts%UseOpenMP) Then
      If (OpenMPOpts%nParticleThreads < 1) Then
        Call Message('WARNING: Number of particle threads has to be positive, adjusted to 1',1)
        OpenMPOpts%nParticleThreads=1
      End If
      If (OpenMPOpts%nParticleThreads > OpenMPMaxThreads) Then
        Call Message('WARNING: Number of particle threads has to be smaller than ' // &
                     Trim(Int2Char(OpenMPMaxThreads)) //                              &
                     ', adjusted to ' // Trim(Int2Char(OpenMPMaxThreads)) // '.',     &
                     1                                                                &
              )
        OpenMPOpts%nParticleThreads=OpenMPMaxThreads
      End If
      If (OpenMPOpts%nParticleUpdateThreads < 1) Then
        Call Message('WARNING: Number of particle update threads has to be positive, adjusted to 1',1)
        OpenMPOpts%nParticleUpdateThreads=1
      End If
      If (OpenMPOpts%nParticleUpdateThreads > OpenMPMaxThreads) Then
        Call Message('WARNING: Number of particle update threads has to be smaller than ' // &
                     Trim(Int2Char(OpenMPMaxThreads)) //                                     &
                     ', adjusted to ' // Trim(Int2Char(OpenMPMaxThreads)) // '.',            &
                     1                                                                       &
              )
        OpenMPOpts%nParticleUpdateThreads=OpenMPMaxThreads
      End If
      If (OpenMPOpts%nChemistryThreads < 1) Then
        Call Message('WARNING: Number of chemistry threads has to be positive, adjusted to 1',1)
        OpenMPOpts%nChemistryThreads=1
      End If
      If (OpenMPOpts%nChemistryThreads > OpenMPMaxThreads) Then
        Call Message('WARNING: Number of chemistry threads has to be smaller than ' // &
                     Trim(Int2Char(OpenMPMaxThreads)) //                               &
                     ', adjusted to ' // Trim(Int2Char(OpenMPMaxThreads)) // '.',      &
                     1                                                                 &
              )
        OpenMPOpts%nChemistryThreads=OpenMPMaxThreads
      End If
      If (OpenMPOpts%nOutputGroupThreads < 1) Then
        Call Message('WARNING: Number of output group threads has to be positive, adjusted to 1',1)
        OpenMPOpts%nOutputGroupThreads=1
      End If
      If (OpenMPOpts%nOutputGroupThreads > OpenMPMaxThreads) Then
        Call Message('WARNING: Number of output group threads has to be smaller than ' // &
                     Trim(Int2Char(OpenMPMaxThreads)) //                                  &
                     ', adjusted to ' // Trim(Int2Char(OpenMPMaxThreads)) // '.',         &
                     1                                                                    &
              )
        OpenMPOpts%nOutputGroupThreads=OpenMPMaxThreads
      End If
      If (OpenMPOpts%nOutputProcessThreads < 1) Then
        Call Message('WARNING: Number of output process threads has to be positive, adjusted to 1',1)
        OpenMPOpts%nOutputProcessThreads=1
      End If
      If (OpenMPOpts%nOutputProcessThreads > OpenMPMaxThreads) Then
        Call Message('WARNING: Number of output process threads has to be smaller than ' // &
                     Trim(Int2Char(OpenMPMaxThreads)) //                                  &
                     ', adjusted to ' // Trim(Int2Char(OpenMPMaxThreads)) // '.',         &
                     1                                                                    &
              )
        OpenMPOpts%nOutputProcessThreads=OpenMPMaxThreads
      End If
      If ((OpenMPopts%ParallelMetProcess) .and. (.not. OpenMPOpts%ParallelMetRead)) Then
        Call Message('WARNING: Parallel MetRead has to be .true. if Parallel MetProcess is .true. ' // &
                     'Parallel MetRead has been adjusted to .true.',1)
        OpenMPOpts%ParallelMetRead=.true.
      End If
    Else
      OpenMPOpts%nParticleThreads=1
      OpenMPOpts%nParticleUpdateThreads=1
      OpenMPOpts%nChemistryThreads=1
      OpenMPOpts%nOutputGroupThreads=1
      OpenMPOpts%nOutputProcessThreads=1
      OpenMPOpts%ParallelMetRead=.false.
      OpenMPOpts%ParallelMetProcess=.false.
    End If
    
! Initialise (create) synchronisation locks
    If (OpenMPOpts%ParallelMetRead) Then
      Do i = 1, 3
        !$ Call omp_init_lock(ParallelReadLocks(i))
      End Do
    End If

  End If
   
End Subroutine InitOpenMPOpts

!------------------------------------------------------------------------------------------------------------


subroutine To_Upper(str)
! Convert a string to upper case

  Implicit None
  
  Character(*), Intent(InOut) :: str
  Integer                     :: i
 
  Do i = 1, len(str)
    Select Case(str(i:i))
      Case("a":"z")
        str(i:i) = AChar(IAChar(str(i:i))-32)
    End Select
  End Do 

End Subroutine To_Upper

!------------------------------------------------------------------------------------------------------------

Subroutine InitIOSynchronisation(OpenMPOpts)
! Initialise synchronisation between IO- and worker- thread by setting the appropriate locks

  Type(OpenMPOpts_), Intent(In) :: OpenMPOpts

  Integer :: myid

  If (OpenMPOpts%ParallelMetRead) Then
  myid = -1
  !$ myid=omp_get_thread_num()
   If (myid .EQ. WorkerThreadID) Then
     !$ Call omp_set_lock(ParallelReadLocks(2))
     ParallelReadLockIndex = 1
   End If
   If (myid .EQ. IOThreadID) Then
     !$ Call omp_set_lock(ParallelReadLocks(1))
     !$ Call omp_set_lock(ParallelReadLocks(3))
     ParallelReadLockIndex = 1
   End If
  End If

End Subroutine InitIOSynchronisation

!------------------------------------------------------------------------------------------------------------

Subroutine CheckOpenMPOpts(OpenMPOpts)
! Check for right number of threads in OpenMP configuration and print out (fatal) error messages if necessary

  Implicit None

  Type(OpenMPOpts_), Intent(InOut) :: OpenMPOpts

  ! Number of threads
  Integer :: nThreads

  ! Loop variable
  Integer :: i

  If (.not. OpenMPOpts%Initialised) Then
    Call Message('FATAL ERROR: The OpenMP options have not been initialised', 3)
  End If
 
  nThreads = -1

  !$OMP PARALLEL
  !$ nThreads = omp_get_num_threads()
  !$OMP END PARALLEL

  ! Print out OpenMP options

  If (nThreads .EQ. -1) Then
    Call Message('Code compiled in serial mode.',0)
    ! $$ Print out warning if we are trying to use OpenMP
  Else
    Call Message('Code compiled in parallel mode.',0)
    Call Message(                                           &
           'Number of Particle Threads          : '      // &
           Int2Char(OpenMPOpts%nParticleThreads,3),         &
           0                                                &
         )
    Call Message(                                           &
           'Number of Particle Update Threads   : '      // &
           Int2Char(OpenMPOpts%nParticleUpdateThreads,3),   &
           0                                                &
         )
    Call Message(                                           &
           'Number of Chemistry Threads         : '      // &
           Int2Char(OpenMPOpts%nChemistryThreads,3),        &
           0                                                &
         )
    Call Message(                                           &
           'Number of Output Group Threads      : '      // &
           Int2Char(OpenMPOpts%nOutputGroupThreads,3),      &
           0                                                &
         )
    Call Message(                                           &
           'Number of Output Process Threads    : '      // &
           Int2Char(OpenMPOpts%nOutputProcessThreads,3),    &
           0                                                &
         )
    If (OpenMPOpts%ParallelMetRead) Then
      NumberOfParallelIOThreads = 2
      If (OpenMPOpts%ParallelMetprocess) Then
        Call Message('Reading and processing data with parallel IO thread',0)
      Else   
        Call Message('Reading data with parallel IO thread',0)
      End If
    Else
      NumberOfParallelIOThreads = 1
      Call Message('Parallel MetRead disabled',0)
    End If

  End If
  
End Subroutine CheckOpenMPOpts

!------------------------------------------------------------------------------------------------------------

Function GetNumberOfParallelIOThreads()
! Return number of parallel IO threads

   Implicit None

   Integer :: GetNumberOfParallelIOThreads

   GetNumberOfParallelIOThreads = NumberOfParallelIOThreads

End Function GetNumberOfParallelIOThreads

!------------------------------------------------------------------------------------------------------------

Subroutine FinaliseOpenMPConfiguration(OpenMPOpts)
! Finalise OpenMP module by destroying the synchronisation locks

  Implicit None
  
  Type(OpenMPOpts_), Intent(In) :: OpenMPOpts

  Integer :: i

  If (OpenMPOpts%ParallelMetRead) Then
    Do i = 1, 3
      !$ Call omp_destroy_lock(ParallelReadLocks(i))
    End Do
  End If

End Subroutine FinaliseOpenMPConfiguration

!------------------------------------------------------------------------------------------------------------

Subroutine LookaheadFileReadRequest()
  
  Implicit None 
  
  !$ Call omp_unset_lock(ParallelReadLocks(ParallelReadLockIndex))
  ParallelReadLockIndex = MOD(ParallelReadLockIndex,3)+1
  
End Subroutine LookaheadFileReadRequest

!------------------------------------------------------------------------------------------------------------

Subroutine LookaheadFileReadRequestWait()
  
  Implicit None
  
  !$ Call omp_set_lock(ParallelReadLocks(ParallelReadLockIndex))
  ParallelReadLockIndex = MOD(ParallelReadLockIndex,3)+1
  
End Subroutine LookaheadFileReadRequestWait

!------------------------------------------------------------------------------------------------------------

Subroutine LookaheadFileReadComplete()

  Implicit None
  
  !$ Call omp_unset_lock(ParallelReadLocks(ParallelReadLockIndex))
  ParallelReadLockIndex = MOD(ParallelReadLockIndex,3)+1

End Subroutine LookaheadFileReadComplete

!------------------------------------------------------------------------------------------------------------

Subroutine LookaheadFileReadCompleteWait()

  Implicit None
  
  !$ Call omp_set_lock(ParallelReadLocks(ParallelReadLockIndex))
  ParallelReadLockIndex = MOD(ParallelReadLockIndex,3)+1

End Subroutine LookaheadFileReadCompleteWait

End Module OpenMPModule
