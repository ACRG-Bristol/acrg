! Module:  Preprocessor Module

Module PreprocessorModule

! Code for preprocessing a Fortran source file according to various 'dollar' 
! directives.
!
! Three directives are catered for.
! 1) Begin block directive
!    Format: '$ BeginBlock blockname'
! 2) End block directive
!    Format: '$ EndBlock blockname'
! 3) Use block directive
!    Format: '$ UseBlock blockname string1 string2'
! Each directive must start in column 1 and must be the only material on the line.
! Extra spaces are permitted between each component and after the last component, 
! while the space(s) after the $ can be omitted. Blockname, string1 and string2 cannot 
! be longer than 32 characters, be blank, or contain spaces.
!
! A begin block and end block directive, which must give the same block name, are used 
! to identify a block of code. This code is not included directly in the output file.
! A use block statment is used to include a previously identified block of code in the
! output file with all occurrances of string1 replaced by string2. Currently, only the 
! most recently defined block can be used. 
  
!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Type :: DollarData_ ! Information about a line of code.
  Logical        Test       !} These have different meanings when the type is used to
  Character(132) CharString !} return data from IsDollar, IsBeginBlock, IsEndBlock and
                            !} IsUseBlock. 
                            !} 1) IsDollar: Test indicates if the line is a dollar 
                            !}    line (a line with a dollar in column 1) and 
                            !}    CharString contains the left-justified text 
                            !}    following the $.
                            !} 2) IsBeginBlock: Test indicates if the line is a begin
                            !}    block directive and (if it is) CharString gives the 
                            !}    block name. 
                            !} 3) IsEndBlock: Test indicates if the line is an end
                            !}    block directive and (if it is) CharString gives the
                            !}    block name. 
                            !} 4) IsUseBlock: Test indicates if the line is a use  
                            !}    block directive and (if it is) CharString gives the 
                            !}    block name, string 1 and string 2.  
End Type DollarData_

!-------------------------------------------------------------------------------------------------------------
    
Contains

!-------------------------------------------------------------------------------------------------------------

Function IsDollar(TextLine)
! Tests whether TextLine is a 'dollar line' (a line of text that has a $ character in  
! the first position). Also gives an error if the input text line is greater than 132 
! characters in length.
      
  Implicit None      
  ! Argument list:      
  Character(*), Intent(In) :: TextLine ! Line of text to tested to see if its a dollar
                                       ! line.      
  ! Function result:      
  Type(DollarData_) IsDollar ! Indicates if the line is a dollar line and gives the
                             ! left justified text following the '$'.      
  ! Locals:       
  Integer        nChars      
  Character(132) CharString 
  Logical        Test
      
  nChars = Len(TextLine)
      
  If (nChars == 0) Then
    Test       = .false.
	CharString = ''
  Else If (nChars > 132) Then
    Write (6, *) 'Error: Line in input file exceeds the maximum permitted ' // &
                 'length of 132 characters'
    Stop
  Else If (TextLine(1:1) == '$') Then
    Test       = .true.
    CharString = AdjustL(TextLine(2:))
  Else
    Test       = .false.
	CharString = ''
  End If
      
  IsDollar = DollarData_(Test, CharString)
      
End Function IsDollar

!-------------------------------------------------------------------------------------------------------------
    
Function IsBeginBlock(TextLine)
! Tests whether TextLine is a begin block directive. Also gives an error if the block 
! name is greater than 32 characters in length, is blank, or contains spaces.
      
  Implicit None
  ! Argument list:      
  Character(*), Intent(In) :: TextLine ! Line of text to tested to see if its a begin
                                       ! block directive.     
  ! Function result:     
  Type(DollarData_) IsBeginBlock ! Indicates whether the line is a begin block 
                                 ! directive and (if it is) gives the block name.   
  ! Locals:      
  Character(132)    CharString
  Character(132)    BlockName      
  Logical           Test     
  Type(DollarData_) IsDollarLine
      
  Test      = .false.
  BlockName = ''
      
  IsDollarLine = IsDollar(TextLine)
      
  If (IsDollarLine%Test) Then
    CharString = IsDollarLine%CharString
    If (CharString(1:11) == 'BeginBlock ') Then
      Test      = .true.
      BlockName = AdjustL(CharString(12:))
      If (Len_Trim(BlockName) > 32) Then
        Write (6, *) 'Error: Blockname "' //                                   & 
                     Trim(BlockName)      //                                   &
                     '" exceeds the maximum permitted length of 32 characters' 
        Stop
      Else If (Len_Trim(BlockName) == 0) Then
        Write (6, *) 'Error: Blockname is absent in begin block directive' 
        Stop
      Else If (Index(BlockName, ' ') /= Len_Trim(BlockName) + 1) Then
        Write (6, *) 'Error: Blockname "' // Trim(BlockName) // '" contains spaces' 
        Stop
      End If
    End If
  End If
      
  IsBeginBlock = DollarData_(Test, BlockName)
      
End Function IsBeginBlock

!-------------------------------------------------------------------------------------------------------------
    
Function IsEndBlock(TextLine)
! Tests whether TextLine is an end block directive. Also gives an error if the block 
! name is greater than 32 characters in length, is blank, or contains spaces.
      
  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: TextLine ! Line of text to tested to see if its an end
                                       ! block directive. 
  ! Function result:
  Type(DollarData_) IsEndBlock ! Indicates whether the line is an end block directive 
                               ! and (if it is) gives the block name.  
  ! Locals:
  Character(132)    CharString
  Character(132)    BlockName
  Logical           Test    
  Type(DollarData_) IsDollarLine
      
  Test      = .false.
  BlockName = ''
      
  IsDollarLine = IsDollar(TextLine)
      
  If (IsDollarLine%Test) Then
    CharString = IsDollarLine%CharString
	If (CharString(1:9) == 'EndBlock ') Then
	  Test      = .true.
	  BlockName = AdjustL(CharString(10:))
	  If (Len_Trim(BlockName) > 32) Then
	    Write(6, *) 'Error: Blockname "' //                                   &
                    Trim(BlockName)      //                                   &
                    '" exceeds the maximum permitted length of 32 characters'
	    Stop
      Else If (Len_Trim(BlockName) == 0) Then
        Write (6, *) 'Error: Blockname is absent in end block directive' 
        Stop
      Else If (Index(BlockName, ' ') /= Len_Trim(BlockName) + 1) Then
        Write (6, *) 'Error: Blockname "' // Trim(BlockName) // '" contains spaces' 
        Stop
	  End If
	End If
  End If
      
  IsEndBlock = DollarData_(Test, BlockName)
      
End Function IsEndBlock

!-------------------------------------------------------------------------------------------------------------
    
Function IsUseBlock(TextLine)
! Tests whether TextLine is a use block directive. 
      
  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: TextLine ! Line of text to tested to see if its a use
                                       ! block directive. 
  ! Function result:
  Type(DollarData_) IsUseBlock ! Indicates whether the line is a use block directive 
                               ! and (if it is) gives the block name, string 1 and 
                               ! string 2.  
  ! Locals: 
  Character(132)    CharString
  Character(132)    BlockNameAndStrings
  Logical           Test
  Type(DollarData_) IsDollarLine
      
  Test                = .false.
  BlockNameAndStrings = ''
      
  IsDollarLine = IsDollar(TextLine)
      
  If (IsDollarLine%Test) Then
    CharString = IsDollarLine%CharString
	If (CharString(1:9) == 'UseBlock ') Then
	  Test                = .true.
	  BlockNameAndStrings = AdjustL(CharString(10:))
	End If
  End If
      
  IsUseBlock = DollarData_(Test, BlockNameAndStrings)
      
End Function IsUseBlock

!-------------------------------------------------------------------------------------------------------------
        
Subroutine ParseString(CharString, BlockName, String1, String2)
! Extracts the arguments blockname, string1 and string2 from a use block directive. 
! Also gives an error if the block name extracted differs from that of the last block 
! read, if the strings cannot be found, or if additional text is found.
      
  Implicit None
  ! Argument list:
  Character(*),  Intent(In)  :: CharString ! String obtained from a use block 
                                           ! directive by IsUseBlock. It should 
                                           ! containing the block name, string 1 and 
                                           ! string 2.
  Character(32), Intent(In)  :: BlockName  ! Name of last block read.
  Character(32), Intent(Out) :: String1    ! String 1 extracted from use block 
                                           ! directive.
  Character(32), Intent(Out) :: String2    ! String 2 extracted from use block 
                                           ! directive.
  ! Local:
  Integer        Pos
  Character(32)  BlockName2
  Character(132) CharString2
      
  CharString2 = CharString
  Pos         = Index(CharString2, ' ')
  BlockName2  = CharString2(1:Pos - 1)
  If (BlockName /= BlockName2) Then
    Write (6, *) 'Error: invalid block name in use block directive'
    Stop
  End If

  CharString2 = AdjustL(CharString2(Pos:))    
  Pos         = Index(CharString2, ' ')
  String1     = CharString2(1:Pos - 1)
  If (Pos == 1) Then
    Write (6, *) 'Error: string 1 not found in use block statement'
    Stop
  End If

  CharString2 = AdjustL(CharString2(Pos:))
  Pos         = Index(CharString2, ' ')
  String2     = CharString2(1:Pos - 1)
  If (Pos == 1) Then
    Write (6, *) 'Error: string 2 not found in use block statement'
    Stop
  End If
      
  If (Len_Trim(CharString2(Pos:)) > 0) Then
    Write (6,*) 'Error: additional text found in UseBlock statement'
    Stop
  End If
      
End Subroutine ParseString

!-------------------------------------------------------------------------------------------------------------
    
Subroutine SubstituteString(String1, String2, TextLine)
! Substitutes every occurrence of one string by another string in a given line of 
! text. Also gives an error if the line becomes longer than 132 characters or if
! recursive substitution occurs. 

  Implicit None
  ! Argument list:
  Character(*),   Intent(In)    :: String1  ! String to be      } Strings are 
                                            ! replaced.         } interpreted
  Character(*),   Intent(In)    :: String2  ! String to replace } strictly, including    
                                            ! String2.          } trailing spaces.
  Character(132), Intent(InOut) :: TextLine ! Line of text to be processed.
  ! Locals:
  Integer nChars          ! Trimmed length of TextLine (which varies during 
                          ! execution).
  Integer Len1            ! Length of String1.
  Integer Len2            ! Length of String2. 
  Integer Position        ! Position in TextLine where next substitution will occur.
  Integer LastSubPosition ! Position in TextLine just after where last substitution 
                          ! was made.
      
  nChars          = Len_Trim(TextLine)
  Len1            = Len(String1)
  Len2            = Len(String2)
  LastSubPosition = 1
      
  Do
    Position = Index(TextLine, String1)
	If (Position == 0) Exit
    If (Position < LastSubPosition) Then
      Write (6, *) 'Error: recursive substitution occurring'
      Stop
    End If
    If (nChars + Len2 - Len1 > 132) Then
      Write (6, *) 'Error: text line becomes longer than 132 characters'
      Stop
    End If
    TextLine = TextLine(1:Position - 1) // String2 // TextLine(Position + Len1:nChars)
    nChars = nChars + Len2 - Len1
    LastSubPosition = Position + Len1
  End Do
      
End Subroutine SubstituteString

!-------------------------------------------------------------------------------------------------------------
    
Subroutine PreprocessFile(FileIn, FileOut)
! Preprocesses a Fortran source file according to the various dollar directives.
     
  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: FileIn  ! Name of file containing code for 
                                      ! pre-processing.
  Character(*), Intent(In) :: FileOut ! Name of file to contain pre-processed code.
  ! Declaration of parameters
  Integer, Parameter :: ReadUnit        = 10  ! Unit number for FileIn.
  Integer, Parameter :: WriteUnit       = 11  ! Unit number for FileOut.
  Integer, Parameter :: MaxLinesInBlock = 128
  ! Locals:
  Integer           OpenStatus
  Integer           ReadStatus
  Integer           WriteStatus
  Integer           CloseStatus
  Integer           TabCheck
  Integer           LineCounter
  Integer           LineNum
  Character(1024)   TextLine
  Character(132)    CharString
  Character(32)     CurrentBlockName
  Character(32)     BlockName
  Character(32)     String1
  Character(32)     String2
  Character(132)    BlockLines(MaxLinesInBlock)
  Logical           ReadingBlock                ! Indicates a block of data is being
                                                ! read.
  Type(DollarData_) IsDollarLine
  Type(DollarData_) IsBeginLine
  Type(DollarData_) IsEndLine
  Type(DollarData_) IsUseLine
      
  ! Open the input text file for processing.      
  Open (Unit     = ReadUnit,	&
        File     = FileIn,		&
        Action   = 'Read',	    &
        Form     = 'Formatted', &
        IOStat   = OpenStatus,  &
        Position = 'Rewind',    &
	    Status   = 'Old')
  If (OpenStatus /= 0) Then
    Write (6, *) 'Error: unable to open input file ' // FileIn
	Stop
  End If
      
  ! Open the output file for writing out the processed source code.      
  Open (Unit     = WriteUnit,	&
        File     = FileOut,		&
        Action   = 'Write',	    &
        Form     = 'Formatted', &
	    IOStat   = OpenStatus,  &
	    Position = 'Rewind',    &
	    Status   = 'Replace')
  If (OpenStatus /= 0) Then
    Write (6, *) 'Error: unable to open output file ' // FileOut
	Stop 
  End If
      
  ! Initialise ReadingBlock to false.
  ReadingBlock = .false.
      
  ! Read one line at a time from the input file and process in the appropriately.      

  Do
        
    ! Read a line into the character variable TextLine.
    Read (Unit = ReadUnit, Fmt = '(A)', IOStat = ReadStatus) TextLine
	
    ! Terminate the loop when the end of file is encountered (or if there is an error 
    ! in reading the input file).
	If (ReadStatus /= 0) Exit
	
	! Check that the input line does not contain any TABs.
	TabCheck = Index(TextLine, '	')
	If (TabCheck /= 0) Then
	  Write (6, *) 'Warning: input text line contains a TAB character. Text line is:'
      Write (6, *) Trim(TextLine)
	End If
	
	! Test if the input line is a dollar statement line.	
	IsDollarLine = IsDollar(Trim(TextLine))
	
	! Process the input line in the appropriate manner.
    	
	If (ReadingBlock) Then
	  
	  ! Processing text in block-reading mode.	  
	  If (IsDollarLine%Test) Then	    
	    ! A dollar line has been found in block-reading mode - only an EndBlock 
        ! statement is valid here. Test for BeginBlock, EndBlock and UseBlock 
        ! statements.
	    IsBeginLine = IsBeginBlock(Trim(TextLine))
	    IsEndLine   = IsEndBlock  (Trim(TextLine))
	    IsUseLine   = IsUseBlock  (Trim(TextLine))
	    If (IsEndLine%Test) Then
	      ! An EndBlock statement has been found and so block-reading mode is switched 
          ! off (ReadingBlock is set to false). If the blockname doesn't match the 
          ! current block then an error is generated. Also if the block length is zero
          ! a warning is given.
	      BlockName = IsEndLine%CharString(1:32)
	      ReadingBlock = .false.
	      If (BlockName /= CurrentBlockName) Then
	        Write (6, *) 'Error: begin/end blocknames do not match'
	        Stop
	      End If	      
	      If (LineCounter == 0) Write (6, *) 'Warning: empty block found'
	    Else If (IsBeginLine%Test) Then	      
	      ! A BeginBlock directive has been found and implies a syntax error.	      
	      Write (6, *) 'Error: BeginBlock directive found while reading a block'
	      Stop	      
	    Else If (IsUseLine%Test) Then
	      ! A UseBlock directive has been found and implies a syntax error.	      
	      Write (6, *) 'Error: UseBlock directive found while reading a block'
	      Stop	      
	    Else	      
	      ! An invalid dollar line has been found and implies a syntax error.	      
	      Write (6, *) 'Error: invalid directive'
	      Stop	      
	    End If	    
	  Else	    
	    ! A standard line has been found in block-reading mode. Increment the line 
        ! counter by 1, check that the line counter does not exceed the maximum number 
        ! of lines allowed in a block, and add the text line to the block.
	    LineCounter = LineCounter + 1
	    If (LineCounter > MaxLinesInBlock) Then
	      Write (6, *) 'Error: number of lines in block exceeds the maximum ' // &
	                   'permitted value'
	      Write (6, *) 'Reset the program parameter MaxLinesInBlock'
	      Stop
	    End If	    
	    BlockLines(LineCounter) = TextLine	    
	  End If
	  
	Else
	  
	  ! Processing text in non-block-reading mode.	  
	  If (IsDollarLine%Test) Then	    
	    ! A dollar line has been found in non-block-reading mode - only a BeginBlock 
        ! or UseBlock statement is valid here. Test for BeginBlock, EndBlock and 
        ! UseBlock statements.
	    IsBeginLine = IsBeginBlock(Trim(TextLine))
	    IsEndLine   = IsEndBlock  (Trim(TextLine))
	    IsUseLine   = IsUseBlock  (Trim(TextLine))	    
	    If (IsBeginLine%Test) Then	      
	      ! A BeginBlock statement has been found and so block-reading mode is 
          ! switched on. (ReadingBlock is set to true, LineCounter is reset to zero, 
          ! and CurrentBlockName is set to the block name).
	      ReadingBlock     = .true.
	      LineCounter      = 0
	      CurrentBlockName = IsBeginLine%CharString(1:32)
	    Else If (IsEndLine%Test) Then
	      ! An EndBlock statement has been found and implies a syntax error
	      Write (6, *) 'Error: EndBlock directive found while not reading block'
	      Stop
	    Else If (IsUseLine%Test) Then
	      ! A UseBlock statement has been found and so the contents of the most recent
	      ! block are writen to the output file with the prescribed string 
          ! substitutions. In more detail, we isolate the arguments of the UseBlock 
          ! directive and check that the block name in the directive matches the 
          ! current block name, and then write out each line in turn with the required 
          ! substitutions.		
	      Call ParseString(IsUseLine%CharString, CurrentBlockName, String1, String2)
	      Do LineNum = 1, LineCounter		
            CharString = BlockLines(LineNum)		  
		    Call SubstituteString(Trim(String1), Trim(String2), CharString)		  
		    Write (Unit = WriteUnit, Fmt = '(A)', IOStat = WriteStatus) &
                    Trim(CharString)		    
		    If (WriteStatus /= 0) Then
	          Write (6, *) 'Error: problem writing line to output file'
	          Stop
	        End If		  
		  End Do		  
	    Else
	      ! An invalid dollar line has been found and implies a syntax error.	      
	      Write (6, *) 'Error: invalid directive'
	      Stop
	    End If
	  Else	    
	    ! A standard line has been found in non-block-reading mode. Write the text 
        ! line to the output file.	    
	    Write (Unit = WriteUnit, Fmt = '(A)', IOStat = WriteStatus) Trim(TextLine)
	    If (WriteStatus /= 0) Then
	      Write (6, *) 'Error: problem writing line to output file'
	      Stop
	    End If	    
	  End If
	  
	End If
	
  End Do
      
  ! Close the input file      
  Close (Unit = ReadUnit, IOStat = CloseStatus)      
  If (CloseStatus /= 0) Then
    Write (6, *) 'Error: unable to close input text file ' // FileIn
	Stop
  End If
      
  ! Close the output file      
  Close (Unit = WriteUnit, IOStat = CloseStatus)      
  If (CloseStatus /= 0) Then
    Write (6, *) 'Error: unable to close output text file ' // FileOut
	Stop
  End If

  ! Write success message.
  Write (6, *) 'File '                                       // &
               Trim(FileIn)                                  // &
               ' successfully preprocessed to produce file ' // &
               Trim(FileOut)                                 // &
               '.'
      
End Subroutine PreprocessFile

!-------------------------------------------------------------------------------------------------------------
    
End Module PreProcessorModule
