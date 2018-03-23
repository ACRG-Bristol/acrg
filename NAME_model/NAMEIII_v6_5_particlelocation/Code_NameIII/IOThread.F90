! Module:  IOThread Module

Module IOThreadModule

! This module controls the IO thread for parallel MetData reads

!-------------------------------------------------------------------------------------------------------------

Use ServiceModule
Use CommonMetModule, Only : CommonMet_
Use MetsModule,      Only : Mets_
Use OpenMPModule,    Only : IOLookaheadComplete           &
                           ,LookaheadFileReadComplete     &
                           ,LookaheadFileReadRequestWait  &
                           ,OpenMPOpts_
Use NWPMetModule,    Only : ReadNWPMet,                   &
                            ProcessNWPMet,                &
                            NWPMet_

Use TimerModule
Use ErrorAndMessageModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

Private

Public :: StartIOThread
Public :: IOThreadModuleTimerSummary
Public :: IOThreadModuleTimerInitialise

Type(Timer_), Save :: IONWPMetReadTimer,   &        ! Timers
                      IONWPMetProcessTimer
Logical      :: TimersInitialised = .false.

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

! Initialise timers

Subroutine IOThreadModuleTimerInitialise(TimerOpts)

  Implicit None

  Type(TimerOpts_) :: TimerOpts

  ! Create timers

  If (.not. TimersInitialised) Then
    Call TimerCreate(IONWPMetReadTimer,   "IOThread_NWPMetRead"   ,TimerOpts)
    Call TimerCreate(IONWPMetProcessTimer,"IOThread_NWPMetProcess",TimerOpts)
    TimersInitialised = .true.
  End If

End Subroutine IOThreadModuleTimerInitialise

!-------------------------------------------------------------------------------------------------------------

! Output timer summary information

Subroutine IOThreadModuleTimerSummary()

  Implicit None

  If (TimersInitialised) Then
    Call TimerWriteSummary(IONWPMetReadTimer)
    Call TimerWriteSummary(IONWPMetProcessTimer)
  End If

End Subroutine IOThreadModuleTimerSummary

!-------------------------------------------------------------------------------------------------------------

! Main IO thread routine.
! As long as IOLookaheadComplete is .false., the IO thread loops over the different MetModules and
! checks whether one of them has data marked for prefetching. It the reads the data into the prefetch
! buffer and processes it, depending on the flag ProcessDataOnIOThread (which is set in OpenMPModule).

Subroutine StartIOThread(Coords,Grids,Mets,Units,OpenMPOpts)

  Implicit None

  Type(Coords_), Intent(In)         :: Coords      ! Collection of coord systems.
  Type(Grids_),  Intent(In)         :: Grids       ! Collection of grids.
  Type(Mets_),   Intent(In), Target :: Mets        ! Collection of met module instance states.
  Type(Units_),  Intent(InOut)      :: Units       ! Collection of information on input/output unit numbers.
  Type(OpenMPOpts_), Intent(In)     :: OpenMPOpts  ! OpenMP options

  Integer                           :: iMetMod     ! MetModule type index
  Integer                           :: iMet        ! MetModule index
  Integer                           :: NRequests=0
  Type(CommonMet_), Pointer         :: CommonMet   ! } Abbreviation for the common part of the Met
                                                   ! } module instance.
  Integer                           :: iCase
  Integer                           :: iMetCase
  Type(Time_)                       :: MetTime     ! } This is the Metdata time in the MetModule,
                                                   ! } NWPMet%PrefetchTime will be set to MetTime after
                                                   ! } the data has been read
  Integer                           :: PrefetchBufferIndex  ! Index of the prefetch data buffer
  Integer                           :: NewBufferIndex       ! Index of the new data buffer
  Logical                           :: Error       ! Error flag
  Type(NWPMet_), Pointer            :: NWPMet      ! Pointer to Met module instance

  Logical                           :: AllowSkip = .false.    ! Do not allow skipping of MetData
  Logical                           :: Skip


  ! Loop until the worker thread set IOLookaheadComplete to .true.
  Do While (.not.IOLookaheadComplete)

    ! Synchronisation: Tell the worker thread that any required data has been read
    Call LookaheadFileReadComplete()

    ! Synchronisation: Wait for next prefetch request
    Call LookaheadFileReadRequestWait()

    If (.not. IOLookaheadComplete) Then

      ! Loop over Met modules and check whether one of them has been marked for prefetch
      NRequests=0
      Do iMetMod = 1, Mets%nMetMods
        Do iMet    = 1, Mets%nMets(iMetMod)
          CommonMet => Mets%C(iMetMod, iMet)%P
          If (CommonMet%Prefetch) Then
            ! Check whether UpdateOnDemand is used for this met module. If so, then abort with a
            ! fatal error
            If (CommonMet%UpdateOnDemand) Then
              Call Message('ERROR: UpdateOnDemand can not be used together with parallel MetRead',4)
            End If

            NRequests=NRequests+1
            iCase=CommonMet%iCase
            iMetCase=CommonMet%iMetCase
            MetTime=CommonMet%MetTime
            PrefetchBufferIndex=CommonMet%PrefetchBufferIndex
            NewBufferIndex=CommonMet%NewBufferIndex

            ! Check this is an NWPMet
            If ((Mets%nNWPMets /= 0) .and. (iMetMod == Mets%NWPMets(1)%C%iMetMod)) Then

              NWPMet=>Mets%NWPMets(iMet)
!             Read the NWPMet data

              Call TimerOn(IONWPMetReadTimer)
              Call ReadNWPMet(MetTime, iCase, iMetCase, Coords, Grids,   &
                              AllowSkip, Error, Skip, NWPMet, Units,     &
                              PrefetchBufferIndex)
              Call TimerOff(IONWPMetReadTimer)
              Call TimerWriteLast(IONWPMetReadTimer)

              ! Process data, if necessary
              If (OpenMPOpts%ParallelMetProcess) Then
                Call TimerOn(IONWPMetProcessTimer)
                Call ProcessNWPMet(Coords, Grids, NWPMet, NewBufferIndex, PrefetchBufferIndex)
                Call TimerOff(IONWPMetProcessTimer)
                Call TimerWriteLast(IONWPMetProcessTimer)
              End If

              ! Update Prefetch time to signal that buffer has been read
              NWPMet%PrefetchTime=Time2ShortTime(MetTime)
            Else
              Call Message(                                                            &
                           'Error: IO thread received a prefetch that is not NWPMet.', &
                           2                                                           &
                   )
              ! abort here
            End If
            CommonMet%Error=Error
!           Reset Flag
            CommonMet%Prefetch=.false.
          End If
        End Do
      End Do
      Call Message(                                                          &
                   Trim(Int2Char(NRequests,Length=3))                    //  &
                   ' Met read requests found and completed by IO thread',     &
                   0                                                         &
           )

    End If

  End Do

End Subroutine StartIOThread

!-------------------------------------------------------------------------------------------------------------

End Module IOThreadModule
