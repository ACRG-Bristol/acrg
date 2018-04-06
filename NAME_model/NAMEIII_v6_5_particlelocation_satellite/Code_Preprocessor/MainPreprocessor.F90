! Module:  Main Preprocessor Module

!-------------------------------------------------------------------------------------------------------------

Program Preprocessor
! Main routine for preprocessing source code for NAME III using 'dollar' directives.
  
  Use PreprocessorModule
  Use SystemModule
   
  Implicit None
  ! Locals: 
  Integer       :: nArguments
  Character(80) :: Arguments(2)

  ! Get command line arguments.
  Call GetCommandLineArguments(nArguments, Arguments)
  If (nArguments /= 2) Then
    Write(6,*) 'Error: incorrect number of command line arguments' 
    Stop
  End If

  Call PreProcessFile(Arguments(1), Arguments(2))
    
End Program Preprocessor

!-------------------------------------------------------------------------------------------------------------
