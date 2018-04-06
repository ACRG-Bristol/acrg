! Module:  Headed File Module

Module HeadedFileModule

! This module provides code for reading headed files.

! Module overview
! ---------------

! Headed files have a format in which the only significant content takes the form of a number of blocks. Each
! block has the form:
!
! Block keyword followed, for certain block keywords, by a block name
! Comma separated list of column keywords - this can be repeated several times
! Comma separated list of data items - this can be repeated several times
! Blank line or end of file.
!
! Each comma separated list in a block must have the same number of entries, so that the block in effect
! consists of a set of columns. Blank or zero length column keywords or data items are possible - but any such
! entry which is last on a line must be followed by a trailing comma. Otherwise, trailing commas at the end of
! lines are optional. Blank or zero length block keywords or block names are not allowed. Spaces can be
! included before, between or after any of the elements (block keyword, block name, column keywords, data
! items, comma). As a result spaces at the start or end of any element are not significant. Block keywords and
! column keywords are not case sensitive, but block names and data items may be, depending on how they are
! interpreted by the calling code.
!
! Outside the blocks any material can appear in the file and will have no effect, provided no line in this
! material starts with a block keyword or a block keyword preceeded by spaces.
!
! The block keyword says what sort of information each block contains and the column keywords say what
! information each block column contains. Block names serve to allow the possibility of distinguishing
! different blocks with the same block keyword, and are usually used when each block is presenting information
! on a particular instance of a certain type of entity.
!
! Column keywords starting with '!' are ignored together with the associated column. This provides a way of
! putting comments in a column.
!
! The module can be set up to read an ordered set of files. A single block cannot be split between the files,
! but otherwise the effect is the same as if the files were concatenated in order with blank lines added
! between the files.
!
! This module makes no assumptions about what files to read, what block keywords and column keywords to search
! for, whether block names are needed, or the number of lines of column keywords required. These choices need
! to be made by calling InitAndAddHFBlockForm and InitHFState before the file(s) are read by calling ReadHFs.
! The number of lines of data items does not need to be specified when the module is set up and is determined
! only by what's in the file.
!
! The call to InitAndAddHFBlockForm also enables information on the column specifications, i.e. how to
! identify particular columns of data, to be set. Often this is just a question of looking for a column with a
! particular column keyword, but, if there are multiple lines of column keywords, things can be more
! complicated with the need to look for certain combinations of column keywords in the various lines.
!
! Not all the column specifications to search for have to be actually matched by column keyword combinations
! in a block. The effect of this is the same as if the column specification was matched but with the
! corresponding data items being blank. If default data values are provided by the calling routine, the
! defaults are used when a data item is blank. It is possible to define the column specifications so that some
! are functionally identical. In this case the first column to match the column specification will be deemed
! to match the first of the functionally identical column specifications, the second to match the second etc,
! until either the columns in the block or the column specifications are exhausted (with the latter giving an
! error - see below).
!
! The block keyword names to search for do not have to be present in the file or can occur more than once. If
! not present, this means no data is provided for this block keyword. The effect is the same if the block
! occurs one or more times with no lines of data items in the block(s).
!
! Unrecognised block keywords or column keyword combinations will be ignored together with the rest of the
! block or column and there is an option to give a warning message for unrecognised column keyword
! combinations.
!
! $$ do something about unrecognised block keywords ?
!
! More than one column matching a column specification (or more columns matching than there are functionally
! identical column specifications) will result in an error message. A column matching two or more
! (non-identical) column specifications will also result in an error message. If the flag DisjointColSpecs is
! set, the column specifications are tested at initialisation to ensure the latter possibility cannot occur.
!
! The order of the blocks, the way the data items for a particular block keyword (and, if relevant, a
! particular block name) are distributed between the blocks with this block keyword (and, if relevant, with
! this block name), and the order of data items within the blocks may be significant depending on how the
! module is used. To understand what options are possible it is important to understand the way the data items
! are returned to the calling program (see below). The calling program can make use of whatever information
! this provides.

! Module use
! ----------

! A variable of type HFBlockForms_ is used to store the information defining the expected forms of the blocks.
! This is initialised with InitHFBlockForms and then information defining the expected form of the blocks is
! recorded by calling InitAndAddHFBlockForm once for each block keyword. A variable of type HFState_ stores
! the names of the files to read and other information needed for internal purposes within the module. It is
! initialised and the file names are recorded using InitHFState. Finally ReadHFs is used to read the files
! with results being returned in a variable of type Tokens_. Information on all data lines in all data blocks
! which match the list of block keywords that are being searched for are returned in the order they occur in
! the files together with the block keyword and, if relevant, block name. Successive calls to ReadHFs return
! information on a single data line except where the description of the expected block form for a particular
! block keyword has TwoD = true; in this case the information in any block with this keyword is returned in
! one go (in a 2-D array) instead of line by line (in a 1-D array), but with the line order within the block
! still clear.
!
! BlockNameLengthTest and TokenLengthTest can be used to test if returned block names or tokens have an
! acceptable (trimmed) length.
!
! Token2Int, Token2Std, Token2Log and Token2Time can be used to convert tokens to integers, reals of kind Std,
! logicals and times.
!
! Types Array_ and Arrays_ (input arrays and collections of input arrays) can be used to store information
! read in from a block. InitArrays is used to initialise a collection of input arrays. Then InitArray is used
! to initialise an input array with the required information and AddArray is used to store the input array in
! the collection of input arrays. FindArrayIndex can be used to find the index of an input array in the
! collection of input arrays.
!
! ParseListChar, ParseListInt, ParseListStd and ParseListLog can be used to parse a token containing a
! delimited list and return the list elements as character strings or converted to integers, reals of kind
! Std and logicals.

! Module call tree
! ----------------

! InitHFBlockForms---(ErrorAndMessage.Message
!
! InitAndAddHFBlockForm---(ErrorAndMessage.Message
!                         (String.CaseInsensitiveEq(.CIEq.)
!
! InitHFState---(ErrorAndMessage.Message
!               (String.ConvertFileName
!
! ReadHFs---(OpenCloseHFs---(Unit.OpenFile
!           (               (Unit.CloseUnit
!           (ReadHeaders---(ReadLine---(ErrorAndMessage.Message
!           (              (ErrorAndMessage.Message
!           (              (String.CaseInsensitiveEq(.CIEq.)
!           (              (String.Int2Char
!           (              (String.Int2Ordinal
!           (ReadData---(ReadLine---(ErrorAndMessage.Message
!                       (ErrorAndMessage.Message
!
! BlockNameLengthTest---(ErrorAndMessage.Message
!                       (String.Int2Char
!
! TokenLengthTest---(ErrorAndMessage.Message
!                   (String.Int2Char
!                   (String.Int2Ordinal
!
! Token2Int---(ErrorAndMessage.Message
!             (String.Int2Char
!             (String.Char2Int
!             (String.Int2Ordinal
!
! Token2Std---(ErrorAndMessage.Message
!             (String.Int2Char
!             (String.Char2Std
!             (String.Int2Ordinal
!
! Token2Log---(ErrorAndMessage.Message
!             (String.Int2Char
!             (String.Char2Log
!             (String.Int2Ordinal
!
! Token2Time---(ErrorAndMessage.Message
!              (String.Int2Char
!              (String.Int2Ordinal
!              (Time.Char2Time
!
! InitArrays
!
! InitArray---(BlockNameLengthTest---(ErrorAndMessage.Message
!             (                      (String.Int2Char
!             (TokenLengthTest---(ErrorAndMessage.Message
!             (                  (String.Int2Char
!             (                  (String.Int2Ordinal
!             (ErrorAndMessage.Message
!             (String.Int2Char
!             (String.CaseInsensitiveEq(.CIEq.)
!
! AddArray---(ArrayEq(==)
!            (ErrorAndMessage.Message
!            (String.CaseInsensitiveEq(.CIEq.)
!
! FindArrayIndex---(ErrorAndMessage.Message
!                  (String.CaseInsensitiveEq(.CIEq.)
!
! ParseListChar---(TokenLengthTest---(ErrorAndMessage.Message
!                 (                  (String.Int2Char
!                 (                  (String.Int2Ordinal
!                 (ErrorAndMessage.Message
!
! ParseListInt---(Token2Int---(ErrorAndMessage.Message
!                (            (String.Int2Char
!                (            (String.Char2Int
!                (            (String.Int2Ordinal
!                (ErrorAndMessage.Message
!
! ParseListStd---(Token2Std---(ErrorAndMessage.Message
!                (            (String.Int2Char
!                (            (String.Char2Std
!                (            (String.Int2Ordinal
!                (ErrorAndMessage.Message
!
! ParseListLog---(Token2Log---(ErrorAndMessage.Message
!                (            (String.Int2Char
!                (            (String.Char2Log
!                (            (String.Int2Ordinal
!                (ErrorAndMessage.Message
!
! Calls to routines in other modules are included (indicated by ModuleName.RoutineName), but subsequent calls
! from the called routine are not.

!-------------------------------------------------------------------------------------------------------------

Use GlobalParametersModule
Use ErrorAndMessageModule
Use StringModule
Use UnitModule
Use TimeModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Private
Public :: HFBlockForms_         ! A collection of block form descriptions.
Public :: HFState_              ! Information on a set of headed files and on what is currently being read
                                ! from the files.
Public :: Tokens_               ! A set of tokens as read from a set of headed files. This might contain data
                                ! from a single line or from an entire block.
Public :: Array_                ! An input array.
Public :: Arrays_               ! A collection of input arrays.
Public :: InitHFBlockForms      ! Initialises a collection of block form descriptions.
Public :: InitAndAddHFBlockForm ! Initialises a block form description and adds it to a collection of block
                                ! form descriptions.
Public :: InitHFState           ! Initialises information on a set of headed files and on what is currently
                                ! being read from the files.
Public :: ReadHFs               ! Reads the next set of tokens (if any) from a set of headed files.
Public :: BlockNameLengthTest   ! Checks a block name has an acceptable (trimmed) length.
Public :: TokenLengthTest       ! Checks a token has an acceptable (trimmed) length.
Public :: Token2Int             ! Converts a token to an integer.
Public :: Token2Std             ! Converts a token to a real of kind Std.
Public :: Token2Log             ! Converts a token to a logical variable.
Public :: Token2Time            ! Converts a token to a time.
Public :: InitArrays            ! Initialises a collection of arrays.
Public :: InitArray             ! Initialises an input array.
Public :: AddArray              ! Adds an input array to a collection of input arrays.
Public :: FindArrayIndex        ! Finds the index of an input array.
Public :: ParseListChar         ! Parses a list of character values separated by some delimiter into a
                                ! character array.
Public :: ParseListInt          ! Parses a list of integer values separated by some delimiter into an integer
                                ! array.
Public :: ParseListStd          ! Parses a list of real values separated by some delimiter into a real array.
Public :: ParseListLog          ! Parses a list of logical values separated by some delimiter into a logical
                                ! array.

!-------------------------------------------------------------------------------------------------------------

! $$ make more use of allocation, assumed shape etc to avoid needing to define so many parameters in global
! params.

Type :: HFBlockForm_ ! A block form description.
  Private
  Character(MaxTokenLength) :: BlockKey
  Logical                   :: NamedBlock
  Logical                   :: MatchMultipleLines
  Integer                   :: nHeaderLines
  Integer                   :: nColumnKeys
  Integer                   :: nSetsOfColumnKeys
  Integer                   :: HeaderLinesToMatch(MaxSetsOfColumnKeys)
  Character(MaxTokenLength) :: ColumnKeys(MaxColumnKeys, MaxSetsOfColumnKeys)
  Character(MaxTokenLength) :: Defaults(MaxColumnKeys)
  Logical                   :: TwoD
  Logical                   :: UnrecognisedMessages
  Logical                   :: SameColSpec(MaxColumnKeys, MaxColumnKeys)
  ! BlockKey             :: Block keyword.
  ! NamedBlock           :: Indicates whether the block keyword corresponds to a named or unnamed block.
  ! MatchMultipleLines   :: Indicates that columns are identified by matching a number of header lines.
  ! nHeaderLines         :: Number of column header lines for this block keyword.
  ! nColumnKeys          :: Number of column keywords for this block keyword.
  ! nSetsOfColumnKeys    :: Number of sets of column keywords for this block keyword.
  ! HeaderLinesToMatch   :: If MatchMultipleLines is true, this contains the indices of the header lines to
  !                         match. The active part of the array has size nSetsOfColumnKeys.
  ! ColumnKeys           :: Column keywords for this block keyword. The active part of the array has size
  !                         (nColumnKeys, nSetsOfColumnKeys).
  ! Defaults             :: Default values for the data items. The active part of the array has size
  !                         nColumnKeys.
  ! TwoD                 :: Indicates the block is read as a single '2-D' block of data rather than a line at
  !                         a time.
  ! UnrecognisedMessages :: Indicates warning messages for unrecognised columns will be given.
  ! SameColSpec          :: Indicates which pairs of column specifications are functionally identical.

  ! This information gives rise to the column specifications to be used in reading the headed files. The
  ! criteria for a column in a headed file to match the ith column specification are as follows:
  !     Case 1: MatchMultipleLines = true:
  !     There is a match if, for each j, the HeaderLinesToMatch(j)-th header line matches ColumnKeys(i, j),
  !     with anything regarded as matching ColumnKeys(i, j) if ColumnKeys(i, j) = '*'.
  !     Case 2: MatchMultipleLines = false:
  !     There is a match if any of the header lines matches ColumnKeys(i, j) for any j.
  ! If a column specification is repeated in HFBlockForm_, then only the ith column to match by the above
  ! criteria is deemed to match the ith repetition of the column specification.
  !
  ! Case 1 requires matches to specific header lines while case 2 provides a way of giving alternative sets
  ! of column keywords. Its not possible to do both.

End Type HFBlockForm_

!-------------------------------------------------------------------------------------------------------------

Type :: HFBlockForms_ ! A collection of block form descriptions.
  Private
  Integer            :: nHFBlockForms               ! Number of block form descriptions.
  Type(HFBlockForm_) :: HFBlockForms(MaxBlockForms) ! A block form description.
End Type HFBlockForms_

!-------------------------------------------------------------------------------------------------------------

Type :: HFState_ ! Information on a set of headed files and on what is currently being read from the files.
  Private
  Integer                      :: nFiles
  Character(MaxFileNameLength) :: FileNames(MaxFilesPerSet)
  Integer                      :: FileUnit
  Integer                      :: iFile
  Character(MaxTokenLength)    :: BlockKey
  Integer                      :: BlockKeyCode
  Character(MaxTokenLength)    :: BlockName
  Integer                      :: nColumns
  Integer                      :: ColumnCodes(MaxColumns)
  Logical                      :: ReadingFile
  Logical                      :: ReadingBlock
  ! nFiles       :: Number of files in the set.
  ! FileNames    :: Names of the files in the set.
  ! FileUnit     :: Unit number used for reading the files.
  ! iFile        :: Index in FileNames of the last file opened. Zero indicates that no files have yet been
  !                 opened.
  ! BlockKey     :: Block keyword of the block being read.
  ! BlockKeyCode :: Index in a collection of block form descriptions of the block form description
  !                 corresponding to the block being read.
  ! BlockName    :: Name of named block being read or blank if an unnamed block is being read.
  ! nColumns     :: Number of columns in the block being read.
  ! ColumnCodes  :: Indices in a block form description of the column specifications corresponding to the
  !                 columns in the block being read. Zero indicates that the column does not match any of the
  !                 column specifications.
  ! ReadingFile  :: Indicates a file is being read.
  ! ReadingBlock :: Indicates a block is being read.

  ! Note 'being read' and 'is being read' here means 'being read or just read' and 'is being read or has just
  ! been read'.

End Type HFState_

!-------------------------------------------------------------------------------------------------------------

Type :: Tokens_ ! A set of tokens as read from a set of headed files. This might contain data from a single
                ! line or from an entire block.
  Character(MaxTokenLength) :: BlockKey
  Character(MaxTokenLength) :: BlockName
  Logical                   :: TwoD
  Integer                   :: nLines
  Character(MaxTokenLength) :: Tokens(MaxColumnKeys)
  Character(MaxTokenLength) :: Tokens2d(MaxColumnKeysTwoD, MaxLinesTwoD)
  ! BlockKey  :: Block keyword of the block just read.
  ! BlockName :: Name of named block just read or blank if an unnamed block has just been read.
  ! TwoD      :: Indicates the block was read as a single '2-D' block of data rather than a line at a time.
  ! nLines    :: Number of lines read (1 if TwoD is false, unless there is no more data to be read).
  ! Tokens    :: Tokens read, in order of column specifications. Used when TwoD is false.
  ! Tokens2d  :: 2-d array of tokens read, in order of column specifications / lines. Used when TwoD is true.

  ! Note that, in the description of Tokens and Tokens2d, the order of the column specifications is the order
  ! defined in the corresponding block form description.

End Type Tokens_

!-------------------------------------------------------------------------------------------------------------

Type :: Array_ ! An input array.
  Character(MaxCharLength) :: Name                  ! Name of input array.
  Integer                  :: n                     ! Number of elements in input array.
  Character(MaxCharLength) :: Array(MaxArrayLength) ! The input array.
End Type Array_

!-------------------------------------------------------------------------------------------------------------

Type :: Arrays_ ! A collection of input arrays.
  Integer      :: nArrays           ! Number of input arrays.
  Type(Array_) :: Arrays(MaxArrays) ! Collection of input arrays.
End Type Arrays_

!-------------------------------------------------------------------------------------------------------------

Interface Operator(==) ! Equality of input arrays.
  Module Procedure ArrayEq
End Interface

!-------------------------------------------------------------------------------------------------------------

Contains

!-------------------------------------------------------------------------------------------------------------

Subroutine InitHFBlockForms(HFBlockForms)
! Initialises a collection of block form descriptions.

  Implicit None
  ! Argument List:
  Type(HFBlockForms_), Intent(Out) :: HFBlockForms ! Initialised collection of block form descriptions.

  If (MaxColumnKeysTwoD >  MaxColumnKeys ) Call Message('UNEXPECTED FATAL ERROR in InitHFBlockForms', 4)
  If (MaxLineLength     <= MaxTokenLength) Call Message('UNEXPECTED FATAL ERROR in InitHFBlockForms', 4)

  HFBlockForms%nHFBlockForms = 0

End Subroutine InitHFBlockForms

!-------------------------------------------------------------------------------------------------------------

Subroutine InitAndAddHFBlockForm(                          &
             BlockKey, NamedBlock, MatchMultipleLines,     &
             nHeaderLines, HeaderLinesToMatch,             &
             ColumnKeys, Defaults,                         &
             TwoD, UnrecognisedMessages, DisjointColSpecs, &
             HFBlockForms                                  &
           )
! Initialises a block form description and adds it to a collection of block form descriptions.

  Implicit None
  ! Argument List:
  Character(*),        Intent(In)              :: BlockKey
  Logical,             Intent(In)              :: NamedBlock
  Logical,             Intent(In)              :: MatchMultipleLines
  Integer,             Intent(In)              :: nHeaderLines
  Integer,             Intent(In),    Optional :: HeaderLinesToMatch(:)
  Character(*),        Intent(In)              :: ColumnKeys(:, :)
  Character(*),        Intent(In),    Optional :: Defaults(:)
  Logical,             Intent(In)              :: TwoD
  Logical,             Intent(In)              :: UnrecognisedMessages
  Logical,             Intent(In)              :: DisjointColSpecs
  Type(HFBlockForms_), Intent(InOut), Target   :: HFBlockForms
  ! BlockKey             :: Block keyword.
  ! NamedBlock           :: Indicates whether the block keyword corresponds to a named or unnamed block.
  ! nHeaderLines         :: Number of column header lines for this block keyword.
  ! ColumnKeys           :: Column keywords for this block keyword.
  ! Defaults             :: Default values for the data items.
  ! TwoD                 :: Indicates the block is read as a single '2-D' block of data rather than a line at
  !                         a time.
  ! UnrecognisedMessages :: Indicates warning messages for unrecognised columns will be given.
  ! DisjointColSpecs     :: Indicates the column specifications are tested to ensure any two are either
  !                         functionally identical or are disjoint (i.e. cannot be satisfied by the same
  !                         column).
  ! HFBlockForms         :: A collection of block form descriptions.
  ! Locals:
  Type(HFBlockForm_), Pointer :: HFBlockForm
  Logical                     :: EqualColumnKeys(MaxSetsOfColumnKeys)
  Logical                     :: EqualColumnKeys2(MaxSetsOfColumnKeys, MaxSetsOfColumnKeys)
  Integer                     :: i
  Integer                     :: j
  Integer                     :: k
  Integer                     :: m
  ! HFBlockForm      :: Abbreviation for block form description.
  ! EqualColumnKeys  :} 1-D and 2-D arrays indicating if various column keys are equal.
  ! EqualColumnKeys2 :}
  ! i                :] Loop indices.
  ! j                :]
  ! k                :]
  ! m                :]

  ! Find index of block form description and set HFBlockForm.
  If (HFBlockForms%nHFBlockForms >= MaxBlockForms) Then
    Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)
  End If

  HFBlockForms%nHFBlockForms = HFBlockForms%nHFBlockForms + 1

  HFBlockForm => HFBlockForms%HFBlockForms(HFBlockForms%nHFBlockForms)

  ! BlockKey.
  If (Len(BlockKey) == 0) Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)

  If (                                        &
    Len_Trim(BlockKey) >  MaxTokenLength .or. &
    BlockKey(1:1)      == ' '                 &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)
  End If

  HFBlockForm%BlockKey = BlockKey

  ! NamedBlock and MatchMultipleLines.
  HFBlockForm%NamedBlock         = NamedBlock
  HFBlockForm%MatchMultipleLines = MatchMultipleLines

  ! nHeaderLines.
  If (nHeaderLines < 1) Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)

  HFBlockForm%nHeaderLines = nHeaderLines

  ! nColumnKeys and nSetsOfColumnKeys.
  If (                                                         &
     Size(ColumnKeys, 1) == 0                             .or. &
     Size(ColumnKeys, 1) >  MaxColumnKeys                 .or. &
    (Size(ColumnKeys, 1) >  MaxColumnKeysTwoD .and. TwoD) .or. &
     Size(ColumnKeys, 2) == 0                             .or. &
     Size(ColumnKeys, 2) >  MaxSetsOfColumnKeys                &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)
  End If

  HFBlockForm%nColumnKeys       = Size(ColumnKeys, 1)
  HFBlockForm%nSetsOfColumnKeys = Size(ColumnKeys, 2)

  ! HeaderLinesToMatch.
  If (MatchMultipleLines) Then

    If (.not.Present(HeaderLinesToMatch)) Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)

    If (Size(HeaderLinesToMatch) /= Size(ColumnKeys, 2)) Then
      Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)
    End If

    If (                                                               &
      HeaderLinesToMatch(1)                        < 1            .or. &
      HeaderLinesToMatch(Size(HeaderLinesToMatch)) > nHeaderLines      &
    ) Then
      Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)
    End If

    Do i = 1, Size(HeaderLinesToMatch) - 1
      If (HeaderLinesToMatch(i) >= HeaderLinesToMatch(i + 1)) Then
        Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)
      End If
    End Do

    HFBlockForm%HeaderLinesToMatch(1:Size(HeaderLinesToMatch)) = HeaderLinesToMatch

  Else

    If (Present(HeaderLinesToMatch)) Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)

  End If

  ! ColumnKeys.
  If (Len(ColumnKeys(1, 1)) == 0) Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)

  Do i = 1, Size(ColumnKeys, 1)
    Do j = 1, Size(ColumnKeys, 2)
      If (                                                       &
        Len_Trim(ColumnKeys(i, j))   >    MaxTokenLength    .or. &
        AdjustL(ColumnKeys(i, j))    /=   ColumnKeys(i, j)  .or. &
        (ColumnKeys(i, j)(1:1)     .CIEq. '!')                   &
      ) Then
        Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)
      End If
    End Do
  End Do

  HFBlockForm%ColumnKeys(1:Size(ColumnKeys, 1), 1:Size(ColumnKeys, 2)) = ColumnKeys(:, :)

  ! Defaults.
  If (Present(Defaults)) Then

    If (Size(ColumnKeys, 1) /= Size(Defaults)) Then
      Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)
    End If

    Do i = 1, Size(Defaults)
      If (                                           &
        Len_Trim(Defaults(i)) >  MaxTokenLength .or. &
        AdjustL(Defaults(i))  /= Defaults(i)         &
      ) Then
        Call Message('UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm', 4)
      End If
    End Do

    HFBlockForm%Defaults(1:Size(Defaults)) = Defaults(:)

  Else

    HFBlockForm%Defaults(1:Size(ColumnKeys, 1)) = ' '

  End If

  ! TwoD and UnrecognisedMessages.
  HFBlockForm%TwoD                 = TwoD
  HFBlockForm%UnrecognisedMessages = UnrecognisedMessages

  ! Identify functionally identical column specifications and, if DisjointColSpecs set, ensure any two column
  ! specifications are either functionally identical or are disjoint (i.e. cannot be satisfied by the same
  ! column).
  Do i = 1,     HFBlockForm%nColumnKeys
  Do j = i + 1, HFBlockForm%nColumnKeys

    If (MatchMultipleLines) Then

      Do k = 1, HFBlockForm%nSetsOfColumnKeys
        EqualColumnKeys(k) = ColumnKeys(i, k) .CIEq. ColumnKeys(j, k)
      End Do
      HFBlockForm%SameColSpec(i, j) = All(EqualColumnKeys(1:HFBlockForm%nSetsOfColumnKeys))

      If (DisjointColSpecs) Then
        If (.not. HFBlockForm%SameColSpec(i, j)) Then
          If (                                                                   &
            .not. Any(                                                           &
                    .not. EqualColumnKeys(1:HFBlockForm%nSetsOfColumnKeys) .and. &
                    ColumnKeys(i, 1:HFBlockForm%nSetsOfColumnKeys) /= '*'  .and. &
                    ColumnKeys(j, 1:HFBlockForm%nSetsOfColumnKeys) /= '*'        &
                  )                                                              &
          ) Then
            Call Message(                                                    &
                   'UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm: '    // &
                   'Column specifications are not identical or disjoint',    &
                   4                                                         &
                 )
          End If
        End If
      End If

    Else

      Do k = 1, HFBlockForm%nSetsOfColumnKeys
        Do m = 1, HFBlockForm%nSetsOfColumnKeys
          EqualColumnKeys2(k, m) = ColumnKeys(i, k) .CIEq. ColumnKeys(j, m)
        End Do
      End Do
      HFBlockForm%SameColSpec(i, j) =                                                           &
        All(                                                                                    &
          Any(                                                                                  &
            EqualColumnKeys2(1:HFBlockForm%nSetsOfColumnKeys, 1:HFBlockForm%nSetsOfColumnKeys), &
            Dim = 2                                                                             &
          )                                                                                     &
        )                                                                                       &
        .and.                                                                                   &
        All(                                                                                    &
          Any(                                                                                  &
            EqualColumnKeys2(1:HFBlockForm%nSetsOfColumnKeys, 1:HFBlockForm%nSetsOfColumnKeys), &
            Dim = 1                                                                             &
          )                                                                                     &
        )

      If (DisjointColSpecs) Then
        If (.not. HFBlockForm%SameColSpec(i, j)) Then
          If (Any(EqualColumnKeys2(1:HFBlockForm%nSetsOfColumnKeys, 1:HFBlockForm%nSetsOfColumnKeys))) Then
            Call Message(                                                    &
                   'UNEXPECTED FATAL ERROR in InitAndAddHFBlockForm: '    // &
                   'Column specifications are not identical or disjoint',    &
                   4                                                         &
                 )
          End If
        End If
      End If

    End If

    HFBlockForm%SameColSpec(j, i) = HFBlockForm%SameColSpec(i, j)

  End Do
  End Do

End Subroutine InitAndAddHFBlockForm

!-------------------------------------------------------------------------------------------------------------

Subroutine InitHFState(FileNames, HFState)
! Initialises information on a set of headed files and on what is currently being read from the files.

  Implicit None
  ! Argument List:
  Character(*),   Intent(In)  :: FileNames(:) ! File names.
  Type(HFState_), Intent(Out) :: HFState      ! Variable to be initialised.
  ! Locals:
  Integer :: i ! Loop index.

  If (                                     &
    Size(FileNames) == 0              .or. &
    Size(FileNames) >  MaxFilesPerSet      &
  ) Then
    Call Message('UNEXPECTED FATAL ERROR in InitHFState', 4)
  End If
  Do i = 1, Size(FileNames)
    If (                                               &
      Len_Trim(FileNames(i)) == 0                 .or. &
      Len_Trim(FileNames(i)) >  MaxFileNameLength      &
    ) Then
      Call Message('UNEXPECTED FATAL ERROR in InitHFState', 4)
    End If
  End Do

  HFState%nFiles       = Size(FileNames)
  HFState%iFile        = 0
  HFState%ReadingFile  = .false.
  HFState%ReadingBlock = .false.
  Do i = 1, Size(FileNames)
    HFState%FileNames(i) = ConvertFileName(FileNames(i))
  End Do

End Subroutine InitHFState

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadHFs(HFBlockForms, Tokens, NoMoreData, Units, HFState)
! Reads the next set of tokens (if any) from a set of headed files.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(In)    :: HFBlockForms
  Type(Tokens_),       Intent(Out)   :: Tokens
  Logical,             Intent(Out)   :: NoMoreData
  Type(Units_),        Intent(InOut) :: Units
  Type(HFState_),      Intent(InOut) :: HFState
  ! HFBlockForms :: A collection of block form descriptions.
  ! Tokens       :: A set of tokens as read from the set of headed files.
  ! NoMoreData   :: Indicates there is no more data to be read from the set of headed files.
  ! Units        :: Collection of information on input/output unit numbers.
  ! HFState      :: Information on the set of headed files and on what is currently being read from the files.

  NoMoreData = .false.

  Do

    ! If not reading a file, try to open a file. Exit if failure with NoMoreData set to true.
    If (.not.HFState%ReadingFile) Then

      Call OpenCloseHFs(Units, HFState)
      If (.not.HFState%ReadingFile) Then
        NoMoreData = .true.
        Exit
      End If

    ! If reading a file but not reading a block of data, search for block headers and read headers.
    Else If (.not.HFState%ReadingBlock) Then

      Call ReadHeaders(HFBlockForms, HFState)

    ! If reading a file and reading a block of data, try to read the next line of data items from the block or
    ! the entire block of data items (depending on HFBlockForms%HFBlockForms(HFState%BlockKeyCode)%TwoD). Exit
    ! if successful.
    Else

      Call ReadData(HFBlockForms%HFBlockForms(HFState%BlockKeyCode), Tokens, HFState)
      If (Tokens%nLines /= 0) Exit

    End If

  End Do

End Subroutine ReadHFs

!-------------------------------------------------------------------------------------------------------------

Subroutine OpenCloseHFs(Units, HFState)
! Open next headed file (if any) and close last headed file (if any).

  Implicit None
  ! Argument list:
  Type(Units_),   Intent(InOut) :: Units   ! Collection of information on input/output unit numbers.
  Type(HFState_), Intent(InOut) :: HFState ! Information on the set of headed files and on what is currently
                                           ! being read from the files.
  ! Locals:
  Integer :: IOStat ! Error code for open statement.

  If (HFState%iFile > 0) Call CloseUnit(HFState%FileUnit, Units)

  HFState%iFile = HFState%iFile + 1

  If (HFState%iFile <= HFState%nFiles) Then

    HFState%FileUnit = OpenFile(                           &
                         HFState%FileNames(HFState%iFile), &
                         Units,                            &
                         Status          = 'Old',          &
                         Action          = 'Read',         &
                         FileDescription = 'input file'    &
                       )
    HFState%ReadingFile = .true.

  End If

End Subroutine OpenCloseHFs

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadHeaders(HFBlockForms, HFState)
! Reads the next block header.

  Implicit None
  ! Argument list:
  Type(HFBlockForms_), Intent(In),   Target :: HFBlockForms
  Type(HFState_),      Intent(InOut)        :: HFState
  ! HFBlockForms :: A collection of block form descriptions.
  ! HFState      :: Information on the set of headed files and on what is currently being read from the files.
  ! Locals:
  Type(HFBlockForm_),      Pointer      :: HFBlockForm
  Integer                               :: IOStat
  Character(MaxLineLength)              :: Line
  Integer                               :: k1
  Integer                               :: k2
  Integer                               :: s
  Character(MaxTokenLength)             :: TempString
  Logical,                  Allocatable :: VarMatchesList(:, :, :, :)
  Logical                               :: ColMatchesList(MaxColumns, MaxColumnKeys)
  Logical,                  Allocatable :: VarMatchesComment(:, :)
  Logical                               :: ColMatchesComment(MaxColumns)
  Logical                               :: ColumnsFound(MaxColumnKeys)
  Integer                               :: nAllowed
  Integer                               :: i
  Integer                               :: j
  Integer                               :: k
  Integer                               :: m
  ! HFBlockForm       :: Abbreviation for block form description.
  ! IOStat            :: Error code for read statement.
  ! Line              :: One line as read from a headed file.
  ! k1                :} Positions of commas spanning a column keyword (k1 = 0 for first column keyword and k2
  ! k2                :} = Len_Trim(Line) + 1 for last column keyword if no comma after last column keyword).
  ! s                 :: Variable used in scanning for commas.
  ! TempString        :: Left adjusted version of a column keyword.
  ! VarMatchesList    :: Indicates which column keywords in the column header lines match which column
  !                      keywords in HFBlockForm.
  ! ColMatchesList    :: Indicates which columns match which column specifications in HFBlockForm.
  ! VarMatchesComment :: Indicates which column keywords in the column header lines are comments.
  ! ColMatchesComment :: Indicates which columns are comments.
  ! ColumnsFound      :: Indicates which columns specifications in HFBlockForm have been matched in the block
  !                      being read.
  ! nAllowed          :: Counts the number of times a column specification in HFBlockForm may be matched, i.e.
  !                      the number of times the column specification is defined in HFBlockForm.
  ! i                 :} Loop indices.
  ! j                 :}
  ! k                 :}
  ! m                 :}

  !---------------------------!
  ! Search for block keyword. !
  !---------------------------!

  OuterLoop: Do

    Call ReadLine(HFState, Line, IOStat)

    If (IOStat < 0) Then
      HFState%ReadingFile = .false.
      Return
    End If

    Line = AdjustL(Line)

    Do i = 1, HFBlockForms%nHFBlockForms

      HFBlockForm => HFBlockForms%HFBlockForms(i)

      If (Line(1:Len_Trim(HFBlockForm%BlockKey)) .CIEq. HFBlockForm%BlockKey) Then

        HFState%BlockKey     = HFBlockForm%BlockKey
        HFState%BlockKeyCode = i
        HFState%ReadingBlock = .true.

        Line = AdjustL(Line(Len_Trim(HFBlockForm%BlockKey) + 1: ))

        If (HFBlockForm%NamedBlock) Then
          If (Len_Trim(Line) > MaxTokenLength) Then
            Call Message(                                                 &
                   'FATAL ERROR in reading file '                      // &
                   Trim(HFState%FileNames(HFState%iFile))              // &
                   ': the block name "'                                // &
                   Trim(Line)                                          // &
                   '" which has been given for a block with keyword "' // &
                   Trim(HFBlockForm%BlockKey)                          // &
                   '" is too long',                                       &
                   3                                                      &
                 )
          End If
          If (Len_Trim(Line) == 0) Then
            Call Message(                                         &
                   'FATAL ERROR in reading file '              // &
                   Trim(HFState%FileNames(HFState%iFile))      // &
                   ': no block name '                          // &
                   'has been given for a block with keyword "' // &
                   Trim(HFBlockForm%BlockKey)                  // &
                   '"',                                           &
                   3                                              &
                 )
          End If
        Else
          If (Len_Trim(Line) > 0) Then
            Call Message(                                                 &
                   'FATAL ERROR in reading file '                      // &
                   Trim(HFState%FileNames(HFState%iFile))              // &
                   ': the block name "'                                // &
                   Trim(Line)                                          // &
                   '" which has been given for a block with keyword "' // &
                   Trim(HFBlockForm%BlockKey)                          // &
                   '" should not be present',                             &
                   3                                                      &
                 )
          End If
        End If

        HFState%BlockName = Line

        Exit OuterLoop

      End If

    End Do

  End Do OuterLoop

  HFBlockForm => HFBlockForms%HFBlockForms(HFState%BlockKeyCode)

  !--------------------------------------------------------------------------------!
  ! Read column keywords and identify matches with column keywords in HFBlockForm. !
  !--------------------------------------------------------------------------------!

  Allocate(                         &
    VarMatchesList(                 &
      MaxColumns,                   &
      HFBlockForm%nHeaderLines,     &
      HFBlockForm%nColumnKeys,      &
      HFBlockForm%nSetsOfColumnKeys &
    )                               &
  )
  Allocate(                    &
    VarMatchesComment(         &
      MaxColumns,              &
      HFBlockForm%nHeaderLines &
    )                          &
  )

  Do j = 1, HFBlockForm%nHeaderLines

    Call ReadLine(HFState, Line, IOStat)

    If (IOStat < 0 .or. Line == ' ') Then
      Call Message(                                      &
             'FATAL ERROR in reading file '           // &
             Trim(HFState%FileNames(HFState%iFile))   // &
             ': A block with keyword "'               // &
             Trim(HFState%BlockKey)                   // &
             '" ends before the expected '            // &
             Trim(Int2Char(HFBlockForm%nHeaderLines)) // &
             ' line(s) of column headers',               &
             3                                           &
           )
    End If

    k1 = 0

    Do i = 1, MaxColumns + 1

      If (k1 >= Len_Trim(Line)) Then
        If (j == 1) Then
          HFState%nColumns = i - 1
        Else
          If (HFState%nColumns == i) Then
            Call Message(                                                                             &
                   'FATAL ERROR in reading file '                                                  // &
                   Trim(HFState%FileNames(HFState%iFile))                                          // &
                   ': In a block with keyword "'                                                   // &
                   Trim(HFState%BlockKey)                                                          // &
                   '", the last column keyword in the '                                            // &
                   Trim(Int2Char(j)) // Int2Ordinal(j)                                             // &
                   ' line of column headers is either missing or is blank with no trailing comma',    &
                   3                                                                                  &
                 )
          Else If (HFState%nColumns /= i - 1) Then
            Call Message(                                                          &
                   'FATAL ERROR in reading file '                               // &
                   Trim(HFState%FileNames(HFState%iFile))                       // &
                   ': In a block with keyword "'                                // &
                   Trim(HFState%BlockKey)                                       // &
                   '", the '                                                    // &
                   Trim(Int2Char(j)) // Int2Ordinal(j)                          // &
                   ' line of column headers has a different number of columns',    &
                   3                                                               &
                 )
          End If
        End If
        Exit
      End If

      If (i > MaxColumns) Then
        Call Message(                                              &
               'FATAL ERROR in reading file '                   // &
               Trim(HFState%FileNames(HFState%iFile))           // &
               ': In a block with keyword "'                    // &
               Trim(HFState%BlockKey)                           // &
               '", the '                                        // &
               Trim(Int2Char(j)) // Int2Ordinal(j)              // &
               '" line of column headers has too many columns',    &
               3                                                   &
             )
      End If

      If (k1 + 1 <= MaxLineLength) Then
        s = Scan(Line(k1 + 1: ), ',')
      Else
        s = 0
      End If
      If (s == 0) Then
        k2 = Len_Trim(Line) + 1
      Else
        k2 = k1 + s
      End If

      If (Len_Trim(AdjustL(Line(k1 + 1:k2 - 1))) > MaxTokenLength) Then
        Call Message(                                              &
               'FATAL ERROR in reading file '                   // &
               Trim(HFState%FileNames(HFState%iFile))           // &
               ': In a block with keyword "'                    // &
               Trim(HFState%BlockKey)                           // &
               '", the column keyword "'                        // &
               Trim(AdjustL(Line(k1 + 1:k2 - 1)))               // &
               '" is too long',                                    &
               3                                                   &
             )
      End If

      TempString = AdjustL(Line(k1 + 1:k2 - 1))

      Do k = 1, HFBlockForm%nColumnKeys
      Do m = 1, HFBlockForm%nSetsOfColumnKeys

        VarMatchesList(i, j, k, m) =                                                      &
          (TempString .CIEq. HFBlockForm%ColumnKeys(k, m))                           .or. &
          (HFBlockForm%MatchMultipleLines .and. HFBlockForm%ColumnKeys(k, m) == '*')

      End Do
      End Do

      VarMatchesComment(i, j) = TempString(1:1) == '!'

      ! Message for unrecognised column keyword in the case of 1 header line.
      If (HFBlockForm%nHeaderLines == 1 .and. HFBlockForm%UnrecognisedMessages) Then
        If (                                                                                                &
          .not. Any(VarMatchesList(i, j, 1:HFBlockForm%nColumnKeys, 1:HFBlockForm%nSetsOfColumnKeys)) .and. &
          .not. VarMatchesComment(i, j)                                                                     &
        ) Then
          Call Message(                                    &
                 'WARNING in reading file '             // &
                 Trim(HFState%FileNames(HFState%iFile)) // &
                 ': In a block with keyword "'          // &
                 Trim(HFState%BlockKey)                 // &
                 '", the column keyword "'              // &
                 Trim(AdjustL(Line(k1 + 1:k2 - 1)))     // &
                 '" has not been recognised',              &
                 1                                         &
               )
        End If
      End If

      k1 = k2

    End Do

  End Do

  !------------------------------------------------------------------------!
  ! Identify matches of columns with column specifications in HFBlockForm. !
  !------------------------------------------------------------------------!

  Do i = 1, HFState%nColumns

    If (HFBlockForm%MatchMultipleLines) Then

      Do k = 1, HFBlockForm%nColumnKeys
        ColMatchesList(i, k) = .true.
        Do m = 1, HFBlockForm%nSetsOfColumnKeys
          ColMatchesList(i, k) = ColMatchesList(i, k)                                       .and. &
                                 VarMatchesList(i, HFBlockForm%HeaderLinesToMatch(m), k, m)
        End Do
      End Do

    Else

      Do k = 1, HFBlockForm%nColumnKeys
        ColMatchesList(i, k) = Any(                                                                  &
                                 VarMatchesList(                                                     &
                                   i, 1:HFBlockForm%nHeaderLines, k, 1:HFBlockForm%nSetsOfColumnKeys &
                                 )                                                                   &
                               )
      End Do

    End If

    If (All(VarMatchesComment(i, 1:HFBlockForm%nHeaderLines))) Then
      ColMatchesComment(i) = .true.
    Else If (.not. Any(VarMatchesComment(i, 1:HFBlockForm%nHeaderLines))) Then
      ColMatchesComment(i) = .false.
    Else
      Call Message(                                                                         &
             'FATAL ERROR in reading file '                                              // &
             Trim(HFState%FileNames(HFState%iFile))                                      // &
             ': In a block with keyword "'                                               // &
             Trim(HFState%BlockKey)                                                      // &
             '", the '                                                                   // &
             Trim(Int2Char(i)) // Int2Ordinal(i)                                         // &
             ' column has some header lines starting with "!" indicating a comment and ' // &
             'some not starting with "!". Consistency is required here',                    &
             3                                                                              &
           )
    End If

  End Do

  ! Message for unrecognised column keyword combination in the case of more than 1 header line.
  If (HFBlockForm%nHeaderLines > 1 .and. HFBlockForm%UnrecognisedMessages) Then

    Do i = 1, HFState%nColumns
      If (                                                            &
        .not. Any(ColMatchesList(i, 1:HFBlockForm%nColumnKeys)) .and. &
        .not. ColMatchesComment(i)                                    &
      ) Then
        Call Message(                                    &
               'WARNING in reading file '             // &
               Trim(HFState%FileNames(HFState%iFile)) // &
               ': In a block with keyword "'          // &
               Trim(HFState%BlockKey)                 // &
               '", the '                              // &
               Trim(Int2Char(i)) // Int2Ordinal(i)    // &
               ' column has not been recognised',        &
               1                                         &
             )
      End If
    End Do

  End If

  !-----------------------------------------------------------------------------------!
  ! Error message for a column matching more than one different column specification. !
  !-----------------------------------------------------------------------------------!

  Do i = 1, HFState%nColumns

    Do k = 1, HFBlockForm%nColumnKeys
      If (ColMatchesList(i, k)) Then
        j = k
        Exit
      End If
    End Do

    Do m = j + 1, HFBlockForm%nColumnKeys
      If (ColMatchesList(i, m) .and. .not. HFBlockForm%SameColSpec(j, m)) Then
        Call Message(                                                             &
               'FATAL ERROR in reading file '                                  // &
               Trim(HFState%FileNames(HFState%iFile))                          // &
               ': In a block with keyword "'                                   // &
               Trim(HFState%BlockKey)                                          // &
               '", the '                                                       // &
               Trim(Int2Char(i)) // Int2Ordinal(i)                             // &
               ' column matches more than one different column specification',    &
               3                                                                  &
             )
      End If
    End Do

  End Do

  !----------------------------------------------------------------------------------------------------------!
  ! Assign codes to each column in the block and give error message for a column specification matching more !
  ! than the allowed number of columns.                                                                      !
  !----------------------------------------------------------------------------------------------------------!

  ColumnsFound(:)        = .false.
  HFState%ColumnCodes(:) = 0

  Do i = 1, HFState%nColumns

    nAllowed = 0

    Do k = 1, HFBlockForm%nColumnKeys

      If (ColMatchesList(i, k)) Then
        nAllowed = nAllowed + 1
        If (.not.ColumnsFound(k)) Then
          ColumnsFound(k)        = .true.
          HFState%ColumnCodes(i) = k
          Exit
        End If
      End If

    End Do

    If (HFState%ColumnCodes(i) == 0 .and. Any(ColMatchesList(i, 1:HFBlockForm%nColumnKeys))) Then
      Call Message(                                                       &
             'FATAL ERROR in reading file '                            // &
             Trim(HFState%FileNames(HFState%iFile))                    // &
             ': in a block with keyword "'                             // &
             Trim(HFState%BlockKey)                                    // &
             '", a column specification is matched by more columns '   // &
             'than the number allowed ('                               // &
             Trim(Int2Char(nAllowed))                                  // &
             ') - the '                                                // &
             Trim(Int2Char(i)) // Int2Ordinal(i)                       // &
             ' column is the '                                         // &
             Trim(Int2Char(nAllowed + 1)) // Int2Ordinal(nAllowed + 1) // &
             ' such column',                                              &
             3                                                            &
           )
    End If

  End Do

  Deallocate(VarMatchesList)
  Deallocate(VarMatchesComment)

End Subroutine ReadHeaders

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadData(HFBlockForm, Tokens, HFState)
! Reads the data within a block. If HFBlockForm%TwoD is true an entire block of data is read; otherwise a
! single line of data is read.

  Implicit None
  ! Argument list:
  Type(HFBlockForm_), Intent(In)    :: HFBlockForm ! A block form description.
  Type(Tokens_),      Intent(Out)   :: Tokens      ! A set of tokens as read from the set of headed files.
  Type(HFState_),     Intent(InOut) :: HFState     ! Information on the set of headed files and on what is
                                                   ! currently being read from the files.
  ! Locals:
  Integer                  :: IOStat ! Error code for read statement.
  Character(MaxLineLength) :: Line   ! One line as read from a headed file.
  Integer                  :: iLine  ! Loop index for the lines of the block.
  Integer                  :: k1     !} Positions of commas spanning token (k1 = 0 for first token and k2 =
  Integer                  :: k2     !} Len_Trim(Line) + 1 for last token if no comma after last token).
  Integer                  :: s      ! Variable used in scanning for commas.
  Integer                  :: i      ! Loop index.

  ! Set Tokens%BlockKey, Tokens%BlockName and Tokens%TwoD.
  Tokens%BlockKey  = HFState%BlockKey
  Tokens%BlockName = HFState%BlockName
  Tokens%TwoD      = HFBlockForm%TwoD

  ! Initialise Tokens%Tokens or Tokens%Tokens2d as appropriate.
  If (HFBlockForm%TwoD) Then
    Do i = 1, MaxLinesTwoD
      Tokens%Tokens2d(:, i) = HFBlockForm%Defaults(1:MaxColumnKeysTwoD)
    End Do
  Else
    Tokens%Tokens(:) = HFBlockForm%Defaults(:)
  End If

  ! Loop over lines of data in the named block (one-time loop if HFBlockForm%TwoD is false).
  Do iLine = 1, MaxLinesTwoD + 1

    ! Read line of data.
    Call ReadLine(HFState, Line, IOStat)

    ! End of file reached.
    If (IOStat < 0) Then

      HFState%ReadingFile  = .false.
      HFState%ReadingBlock = .false.
      Tokens%nLines        = iLine - 1
      Exit

    ! End of block reached.
    Else If (Line == ' ') Then

      HFState%ReadingBlock = .false.
      Tokens%nLines        = iLine - 1
      Exit

    ! Too many lines of data in block.
    Else If (HFBlockForm%TwoD .and. iLine > MaxLinesTwoD) Then

      Call Message(                                                      &
             'FATAL ERROR in reading file '                           // &
             Trim(HFState%FileNames(HFState%iFile))                   // &
             ': Too many lines of data items in block with keyword "' // &
             Trim(HFState%BlockKey)                                   // &
             '"',                                                        &
             3                                                           &
           )

    ! Line of data read. Assign token using column codes.
    Else

      k1 = 0
      Do i = 1, HFState%nColumns

        If (k1 + 1 <= MaxLineLength) Then
          s = Scan(Line(k1 + 1: ), ',')
        Else
          s = 0
        End If

        If (s == 0) Then
          k2 = Len_Trim(Line) + 1
          If (i /= HFState%nColumns) Then
            Call Message(                                                        &
                   'FATAL ERROR in reading file '                             // &
                   Trim(HFState%FileNames(HFState%iFile))                     // &
                   ': too few data items on a line in a block with keyword "' // &
                   Trim(HFState%BlockKey)                                     // &
                   '"',                                                          &
                   3                                                             &
                 )
          Else If (k1 == Len_Trim(Line)) Then
            Call Message(                                                        &
                   'FATAL ERROR in reading file '                             // &
                   Trim(HFState%FileNames(HFState%iFile))                     // &
                   ': the last data item on a line in a block with keyword "' // &
                   Trim(HFState%BlockKey)                                     // &
                   '" is either missing or is blank with no trailing comma',     &
                   3                                                             &
                 )
          End If
        Else
          k2 = k1 + s
          If (i == HFState%nColumns .and. Len_Trim(Line) > k2) Then
            Call Message(                                                         &
                   'FATAL ERROR in reading file '                              // &
                   Trim(HFState%FileNames(HFState%iFile))                      // &
                   ': too many data items on a line in a block with keyword "' // &
                   Trim(HFState%BlockKey)                                      // &
                   '"',                                                           &
                   3                                                              &
                 )
          End If
        End If

        If (HFState%ColumnCodes(i) /= 0) Then
          If (k1 + 1 <= k2 - 1 .and. Line(k1 + 1:k2 - 1) /= ' ') Then
            If (Len_Trim(AdjustL(Line(k1 + 1:k2 - 1))) > MaxTokenLength) Then
              Call Message(                                     &
                     'FATAL ERROR in reading file '          // &
                     Trim(HFState%FileNames(HFState%iFile))  // &
                     ': the data item "'                     // &
                     Trim(AdjustL(Line(k1 + 1:k2 - 1)))      // &
                     '" on a line in a block with keyword "' // &
                     Trim(HFState%BlockKey)                  // &
                     '" is too long',                           &
                     3                                          &
                   )
            End If
            If (HFBlockForm%TwoD) Then
              Tokens%Tokens2d(HFState%ColumnCodes(i), iLine) = AdjustL(Line(k1 + 1:k2 - 1))
            Else
              Tokens%Tokens(HFState%ColumnCodes(i)) = AdjustL(Line(k1 + 1:k2 - 1))
            End If
          End If
        End If

        k1 = k2

      End Do

    End If

    If (.not. HFBlockForm%TwoD) Then
      Tokens%nLines = iLine
      Exit
    End If

  End Do

End Subroutine ReadData

!-------------------------------------------------------------------------------------------------------------

Subroutine ReadLine(HFState, Line, IOStat)
! Reads a line.

  Implicit None
  ! Argument list:
  Type(HFState_), Intent(In)  :: HFState
  Character(*),   Intent(Out) :: Line
  Integer,        Intent(Out) :: IOStat
  ! HFState :: Information on the set of headed files and on what is currently being read from the files.
  ! Line    :: One line as read from a headed file.
  ! IOStat  :: Error code for read statement.
  ! Locals:
  Character(1) :: C ! Used to trigger end of record if record length <= Len(Line) rather than if record length < Len(Line).

  ! We use a non-advancing read statement in order to tell whether we have successfully captured the whole line or whether it has 
  ! been truncated.
  Read (HFState%FileUnit, '(2A)', IOStat = IOStat, Advance='No') Line, C

  ! End of file or line read without truncation (do nothing).
  If (Is_IOStat_End(IOStat) .or. Is_IOStat_Eor(IOStat)) Then

    If (.not. Is_IOStat_End(IOStat)) IOStat = 0 ! $$ this not necessary if later use "Is_IOStat_End(IOStat)" not "IOStat < 0".

  ! Line read with truncation (raise fatal error).
  Else If (IOStat == 0) Then

    Call Message(                                            &
           'FATAL Error: Line truncated in reading file ' // &
           Trim(HFState%FileNames(HFState%iFile)),           &
           3                                                 &
         )

  ! Line read with error (raise fatal error).
  Else

    Call Message(                                               &
           'FATAL Error: An error occurred in reading file ' // &
           Trim(HFState%FileNames(HFState%iFile)),              & 
           3                                                    &
         )

  End If

End Subroutine ReadLine

!-------------------------------------------------------------------------------------------------------------

Subroutine BlockNameLengthTest(BlockName, Length, BlockKey)
! Checks a block name has an acceptable (trimmed) length.

  Implicit None
  ! Argument list:
  Character(*), Intent(In) :: BlockName
  Integer,      Intent(In) :: Length
  Character(*), Intent(In) :: BlockKey
  ! BlockName :: Block name.
  ! Length    :: Maximum acceptable length.
  ! BlockKey  :: Input block keyword. A blank value produces an unexpected fatal error in the case of
  !              unacceptable length.

  If (BlockName == ' ' .or. Len_Trim(BlockName) > Length) Then

    If (BlockKey == ' ') Then
      Call Message('UNEXPECTED FATAL ERROR in BlockNameLengthTest.', 4)
    Else If (BlockName == ' ') Then
      Call Message(                                               &
             'FATAL ERROR in reading a block with keyword "'   // &
             Trim(BlockKey)                                    // &
             '": no value has been given for the block name.',    &
             3                                                    &
           )
    Else
      Call Message(                                             &
             'FATAL ERROR in reading a block with keyword "' // &
             Trim(BlockKey)                                  // &
             '": the value "'                                // &
             Trim(BlockName)                                 // &
             '" given for the block name has length > '      // &
             Trim(Int2Char(Length))                          // &
             ' and is too long.',                               &
             3                                                  &
           )
    End If

  End If

End Subroutine BlockNameLengthTest

!-------------------------------------------------------------------------------------------------------------

Subroutine TokenLengthTest(C, Length, Zero, BlockKey, Item, ColumnKey, List, i)
! Checks a token has an acceptable (trimmed) length.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: C
  Integer,      Intent(In)           :: Length
  Logical,      Intent(In)           :: Zero
  Character(*), Intent(In)           :: BlockKey
  Character(*), Intent(In)           :: Item
  Character(*), Intent(In)           :: ColumnKey
  Logical,      Intent(In), Optional :: List
  Integer,      Intent(In), Optional :: i
  ! C         :: Character string.
  ! Length    :: Maximum acceptable length.
  ! Zero      :: Indicates zero length is not acceptable.
  ! BlockKey  :: Block keyword.
  ! Item      :: Name of item read (block name or name associated with an input line). If blank, any message
  !              will not refer to the item read.
  ! ColumnKey :: Column keyword.
  ! List      :: Indicates that C is a list element extracted from a token.
  ! i         :: Indicates that this is the ith data item for the column keyword in the block. This is best
  !              used when reading a block with the TwoD option (otherwise the value of i won't be known) and
  !              with Item omitted or equal to the block name (if i is used and Item refers to an input line,
  !              then i and Item are both referring to the same thing which may make any message confusing).
  ! Locals:
  Character(MaxCharLength) :: ListAsChar ! Value of List converted to a form suitable for including in the
                                         ! error message.
  Character(MaxCharLength) :: iAsChar    ! Value of i converted to a form suitable for including in the error
                                         ! message.

  If ((Zero .and. C == ' ') .or. Len_Trim(C) > Length) Then

    If (Present(List)) Then
      If (List) Then
        ListAsChar = ' a list element within'
      Else
        ListAsChar = ' '
      End If
    Else
      ListAsChar = ' '
    End If

    If (Present(i)) Then
      iAsChar = ' the ' // Trim(Int2Char(i)) // Trim(Int2Ordinal(i)) // ' instance in the block of'
    Else
      iAsChar = ' an instance of'
    End If

    If (BlockKey == ' ' .or. ColumnKey == ' ') Then
      Call Message('UNEXPECTED FATAL ERROR in TokenLengthTest.', 4)
    Else If (Zero .and. C == ' ') Then
      If (Item == ' ') Then
        Call Message(                                             &
               'FATAL ERROR in reading a block with keyword "' // &
               Trim(BlockKey)                                  // &
               '": no value has been given for'                // &
               Trim(ListAsChar)                                // &
               Trim(iAsChar)                                   // &
               ' "'                                            // &
               Trim(ColumnKey)                                 // &
               '".',                                              &
               3                                                  &
             )
      Else
        Call Message(                              &
               'FATAL ERROR in reading item "'  // &
               Trim(Item)                       // &
               '" from a block with keyword "'  // &
               Trim(BlockKey)                   // &
               '": no value has been given for' // &
               Trim(ListAsChar)                 // &
               Trim(iAsChar)                    // &
               ' "'                             // &
               Trim(ColumnKey)                  // &
               '".',                               &
               3                                   &
             )
      End If
    Else
      If (Item == ' ') Then
        Call Message(                                             &
               'FATAL ERROR in reading a block with keyword "' // &
               Trim(BlockKey)                                  // &
               '": the value "'                                // &
               Trim(C)                                         // &
               '" given for'                                   // &
               Trim(ListAsChar)                                // &
               Trim(iAsChar)                                   // &
               ' "'                                            // &
               Trim(ColumnKey)                                 // &
               '" has length > '                               // &
               Trim(Int2Char(Length))                          // &
               ' and is too long.',                               &
               3                                                  &
             )
      Else
        Call Message(                             &
               'FATAL ERROR in reading item "' // &
               Trim(Item)                      // &
               '" from a block with keyword "' // &
               Trim(BlockKey)                  // &
               '": the value "'                // &
               Trim(C)                         // &
               '" given for'                   // &
               Trim(ListAsChar)                // &
               Trim(iAsChar)                   // &
               ' "'                            // &
               Trim(ColumnKey)                 // &
               '" has length > '               // &
               Trim(Int2Char(Length))          // &
               ' and is too long.',               &
               3                                  &
             )
      End If
    End If

  End If

End Subroutine TokenLengthTest

!-------------------------------------------------------------------------------------------------------------

Function Token2Int(C, BlockKey, Item, ColumnKey, List, i)
! Converts a token to an integer.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: C
  Character(*), Intent(In)           :: BlockKey
  Character(*), Intent(In)           :: Item
  Character(*), Intent(In)           :: ColumnKey
  Logical,      Intent(In), Optional :: List
  Integer,      Intent(In), Optional :: i
  ! C         :: Character string to be converted. Acceptable values for the character string are those
  !              accepted by a list directed read of an integer.
  ! BlockKey  :: Block keyword.
  ! Item      :: Name of item read (block name or name associated with an input line). If blank, any message
  !              will not refer to the item read.
  ! ColumnKey :: Column keyword.
  ! List      :: Indicates that C is a list element extracted from a token.
  ! i         :: Indicates that this is the ith data item for the column keyword in the block. This is best
  !              used when reading a block with the TwoD option (otherwise the value of i won't be known) and
  !              with Item omitted or equal to the block name (if i is used and Item refers to an input line,
  !              then i and Item are both referring to the same thing which may make any message confusing).
  ! Function result:
  Integer :: Token2Int ! The integer result of the conversion.
  ! Locals:
  Character(MaxCharLength) :: ListAsChar ! Value of List converted to a form suitable for including in the
                                         ! error message.
  Character(MaxCharLength) :: iAsChar    ! Value of i converted to a form suitable for including in the error
                                         ! message.
  Integer                  :: ErrorCode  ! Error code.

  Token2Int = Char2Int(C, ErrorCode)

  If (ErrorCode /= 0) Then

    If (Present(List)) Then
      If (List) Then
        ListAsChar = ' a list element within'
      Else
        ListAsChar = ' '
      End If
    Else
      ListAsChar = ' '
    End If

    If (Present(i)) Then
      iAsChar = ' the ' // Trim(Int2Char(i)) // Trim(Int2Ordinal(i)) // ' instance in the block of'
    Else
      iAsChar = ' an instance of'
    End If

    If (BlockKey == ' ' .or. ColumnKey == ' ') Then
      Call Message('UNEXPECTED FATAL ERROR in Token2Int.', 4)
    Else If (C == ' ') Then
      If (Item == ' ') Then
        Call Message(                                             &
               'FATAL ERROR in reading a block with keyword "' // &
               Trim(BlockKey)                                  // &
               '": no value has been given for'                // &
               Trim(ListAsChar)                                // &
               Trim(iAsChar)                                   // &
               ' "'                                            // &
               Trim(ColumnKey)                                 // &
               '".',                                              &
               3                                                  &
             )
      Else
        Call Message(                              &
               'FATAL ERROR in reading item "'  // &
               Trim(Item)                       // &
               '" from a block with keyword "'  // &
               Trim(BlockKey)                   // &
               '": no value has been given for' // &
               Trim(ListAsChar)                 // &
               Trim(iAsChar)                    // &
               ' "'                             // &
               Trim(ColumnKey)                  // &
               '".',                               &
               3                                   &
             )
      End If
    Else
      If (Item == ' ') Then
        Call Message(                                             &
               'FATAL ERROR in reading a block with keyword "' // &
               Trim(BlockKey)                                  // &
               '": the value "'                                // &
               Trim(C)                                         // &
               '" given for'                                   // &
               Trim(ListAsChar)                                // &
               Trim(iAsChar)                                   // &
               ' "'                                            // &
               Trim(ColumnKey)                                 // &
               '" is not an integer.',                            &
               3                                                  &
             )
      Else
        Call Message(                             &
               'FATAL ERROR in reading item "' // &
               Trim(Item)                      // &
               '" from a block with keyword "' // &
               Trim(BlockKey)                  // &
               '": the value "'                // &
               Trim(C)                         // &
               '" given for'                   // &
               Trim(ListAsChar)                // &
               Trim(iAsChar)                   // &
               ' "'                            // &
               Trim(ColumnKey)                 // &
               '" is not an integer.',            &
               3                                  &
             )
      End If
    End If

  End If

End Function Token2Int

!-------------------------------------------------------------------------------------------------------------

Function Token2Std(C, BlockKey, Item, ColumnKey, List, i)
! Converts a token to a real of kind Std.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: C
  Character(*), Intent(In)           :: BlockKey
  Character(*), Intent(In)           :: Item
  Character(*), Intent(In)           :: ColumnKey
  Logical,      Intent(In), Optional :: List
  Integer,      Intent(In), Optional :: i
  ! C         :: Character string to be converted. Acceptable values for the character string are those
  !              accepted by a list directed read of a real of kind Std.
  ! BlockKey  :: Block keyword.
  ! Item      :: Name of item read (block name or name associated with an input line). If blank, any message
  !              will not refer to the item read.
  ! ColumnKey :: Column keyword.
  ! List      :: Indicates that C is a list element extracted from a token.
  ! i         :: Indicates that this is the ith data item for the column keyword in the block. This is best
  !              used when reading a block with the TwoD option (otherwise the value of i won't be known) and
  !              with Item omitted or equal to the block name (if i is used and Item refers to an input line,
  !              then i and Item are both referring to the same thing which may make any message confusing).
  ! Function result:
  Real(Std) :: Token2Std ! The real result of the conversion.
  ! Locals:
  Character(MaxCharLength) :: ListAsChar ! Value of List converted to a form suitable for including in the
                                         ! error message.
  Character(MaxCharLength) :: iAsChar    ! Value of i converted to a form suitable for including in the error
                                         ! message.
  Integer                  :: ErrorCode  ! Error code.

  Token2Std = Char2Std(C, ErrorCode)

  If (ErrorCode /= 0) Then

    If (Present(List)) Then
      If (List) Then
        ListAsChar = ' a list element within'
      Else
        ListAsChar = ' '
      End If
    Else
      ListAsChar = ' '
    End If

    If (Present(i)) Then
      iAsChar = ' the ' // Trim(Int2Char(i)) // Trim(Int2Ordinal(i)) // ' instance in the block of'
    Else
      iAsChar = ' an instance of'
    End If

    If (BlockKey == ' ' .or. ColumnKey == ' ') Then
      Call Message('UNEXPECTED FATAL ERROR in Token2Std.', 4)
    Else If (C == ' ') Then
      If (Item == ' ') Then
        Call Message(                                             &
               'FATAL ERROR in reading a block with keyword "' // &
               Trim(BlockKey)                                  // &
               '": no value has been given for'                // &
               Trim(ListAsChar)                                // &
               Trim(iAsChar)                                   // &
               ' "'                                            // &
               Trim(ColumnKey)                                 // &
               '".',                                              &
               3                                                  &
             )
      Else
        Call Message(                              &
               'FATAL ERROR in reading item "'  // &
               Trim(Item)                       // &
               '" from a block with keyword "'  // &
               Trim(BlockKey)                   // &
               '": no value has been given for' // &
               Trim(ListAsChar)                 // &
               Trim(iAsChar)                    // &
               ' "'                             // &
               Trim(ColumnKey)                  // &
               '".',                               &
               3                                   &
             )
      End If
    Else
      If (Item == ' ') Then
        Call Message(                                             &
               'FATAL ERROR in reading a block with keyword "' // &
               Trim(BlockKey)                                  // &
               '": the value "'                                // &
               Trim(C)                                         // &
               '" given for'                                   // &
               Trim(ListAsChar)                                // &
               Trim(iAsChar)                                   // &
               ' "'                                            // &
               Trim(ColumnKey)                                 // &
               '" is not a real.',                                &
               3                                                  &
             )
      Else
        Call Message(                             &
               'FATAL ERROR in reading item "' // &
               Trim(Item)                      // &
               '" from a block with keyword "' // &
               Trim(BlockKey)                  // &
               '": the value "'                // &
               Trim(C)                         // &
               '" given for'                   // &
               Trim(ListAsChar)                // &
               Trim(iAsChar)                   // &
               ' "'                            // &
               Trim(ColumnKey)                 // &
               '" is not a real.',                &
               3                                  &
             )
      End If
    End If

  End If

End Function Token2Std

!-------------------------------------------------------------------------------------------------------------

Function Token2Log(C, BlockKey, Item, ColumnKey, List, i)
! Converts a token to a logical variable.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: C
  Character(*), Intent(In)           :: BlockKey
  Character(*), Intent(In)           :: Item
  Character(*), Intent(In)           :: ColumnKey
  Logical,      Intent(In), Optional :: List
  Integer,      Intent(In), Optional :: i
  ! C         :: Character string to be converted. Acceptable values for the character string are those
  !              accepted by a list directed read of a logical (i.e. those begining T, .T, t, .t, F, .F, f or
  !              .f, possibly preceeded by blanks) and those begining Y, y, N or n (again possibly preceeded
  !              by blanks).
  ! BlockKey  :: Block keyword.
  ! Item      :: Name of item read (block name or name associated with an input line). If blank, any message
  !              will not refer to the item read.
  ! ColumnKey :: Column keyword.
  ! List      :: Indicates that C is a list element extracted from a token.
  ! i         :: Indicates that this is the ith data item for the column keyword in the block. This is best
  !              used when reading a block with the TwoD option (otherwise the value of i won't be known) and
  !              with Item omitted or equal to the block name (if i is used and Item refers to an input line,
  !              then i and Item are both referring to the same thing which may make any message confusing).
  ! Function result:
  Logical :: Token2Log ! The logical result of the conversion.
  ! Locals:
  Character(MaxCharLength) :: ListAsChar ! Value of List converted to a form suitable for including in the
                                         ! error message.
  Character(MaxCharLength) :: iAsChar    ! Value of i converted to a form suitable for including in the error
                                         ! message.
  Integer                  :: ErrorCode  ! Error code.

  Token2Log = Char2Log(C, ErrorCode)

  If (ErrorCode /= 0) Then

    If (Present(List)) Then
      If (List) Then
        ListAsChar = ' a list element within'
      Else
        ListAsChar = ' '
      End If
    Else
      ListAsChar = ' '
    End If

    If (Present(i)) Then
      iAsChar = ' the ' // Trim(Int2Char(i)) // Trim(Int2Ordinal(i)) // ' instance in the block of'
    Else
      iAsChar = ' an instance of'
    End If

    If (BlockKey == ' ' .or. ColumnKey == ' ') Then
      Call Message('UNEXPECTED FATAL ERROR in Token2Log.', 4)
    Else If (C == ' ') Then
      If (Item == ' ') Then
        Call Message(                                             &
               'FATAL ERROR in reading a block with keyword "' // &
               Trim(BlockKey)                                  // &
               '": no value has been given for'                // &
               Trim(ListAsChar)                                // &
               Trim(iAsChar)                                   // &
               ' "'                                            // &
               Trim(ColumnKey)                                 // &
               '".',                                              &
               3                                                  &
             )
      Else
        Call Message(                              &
               'FATAL ERROR in reading item "'  // &
               Trim(Item)                       // &
               '" from a block with keyword "'  // &
               Trim(BlockKey)                   // &
               '": no value has been given for' // &
               Trim(ListAsChar)                 // &
               Trim(iAsChar)                    // &
               ' "'                             // &
               Trim(ColumnKey)                  // &
               '".',                               &
               3                                   &
             )
      End If
    Else
      If (Item == ' ') Then
        Call Message(                                             &
               'FATAL ERROR in reading a block with keyword "' // &
               Trim(BlockKey)                                  // &
               '": the value "'                                // &
               Trim(C)                                         // &
               '" given for'                                   // &
               Trim(ListAsChar)                                // &
               Trim(iAsChar)                                   // &
               ' "'                                            // &
               Trim(ColumnKey)                                 // &
               '" is not "yes" or "no".',                         &
               3                                                  &
             )
      Else
        Call Message(                             &
               'FATAL ERROR in reading item "' // &
               Trim(Item)                      // &
               '" from a block with keyword "' // &
               Trim(BlockKey)                  // &
               '": the value "'                // &
               Trim(C)                         // &
               '" given for'                   // &
               Trim(ListAsChar)                // &
               Trim(iAsChar)                   // &
               ' "'                            // &
               Trim(ColumnKey)                 // &
               '" is not "yes" or "no".',         &
               3                                  &
             )
      End If
    End If

  End If

End Function Token2Log

!-------------------------------------------------------------------------------------------------------------

Function Token2Time(C, BlockKey, Item, ColumnKey, List, i, Interval, Clock)
! Converts a token to a time.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: C
  Character(*), Intent(In)           :: BlockKey
  Character(*), Intent(In)           :: Item
  Character(*), Intent(In)           :: ColumnKey
  Logical,      Intent(In), Optional :: List
  Integer,      Intent(In), Optional :: i
  Logical,      Intent(In), Optional :: Interval
  Logical,      Intent(In), Optional :: Clock
  ! C              :: Character string to be converted. Acceptable values for the character string are those
  !                   accepted by the time module routine Char2Time.
  ! BlockKey       :: Block keyword.
  ! Item           :: Name of item read (block name or name associated with an input line). If blank, any
  !                   message will not refer to the item read.
  ! ColumnKey      :: Column keyword.
  ! List           :: Indicates that C is a list element extracted from a token.
  ! i              :: Indicates that this is the ith data item for the column keyword in the block. This is
  !                   best used when reading a block with the TwoD option (otherwise the value of i won't be
  !                   known) and with Item omitted or equal to the block name (if i is used and Item refers to
  !                   an input line, then i and Item are both referring to the same thing which may make any
  !                   message confusing).
  ! Interval       :} Options for Char2Time routine.
  ! Clock          :}
  ! Function result:
  Type(Time_) :: Token2Time ! The time result of the conversion.
  ! Locals:
  Character(MaxCharLength) :: ListAsChar ! Value of List converted to a form suitable for including in the
                                         ! error message.
  Character(MaxCharLength) :: iAsChar    ! Value of i converted to a form suitable for including in the error
                                         ! message.
  Integer                  :: ErrorCode  ! Error code.

  Token2Time = Char2Time(C, Interval = Interval, Clock = Clock, ErrorCode = ErrorCode)

  If (ErrorCode /= 0) Then

    If (Present(List)) Then
      If (List) Then
        ListAsChar = ' a list element within'
      Else
        ListAsChar = ' '
      End If
    Else
      ListAsChar = ' '
    End If

    If (Present(i)) Then
      iAsChar = ' the ' // Trim(Int2Char(i)) // Trim(Int2Ordinal(i)) // ' instance in the block of'
    Else
      iAsChar = ' an instance of'
    End If

    If (BlockKey == ' ' .or. ColumnKey == ' ') Then
      Call Message('UNEXPECTED FATAL ERROR in Token2Time.', 4)
    Else If (C == ' ') Then
      If (Item == ' ') Then
        Call Message(                                             &
               'FATAL ERROR in reading a block with keyword "' // &
               Trim(BlockKey)                                  // &
               '": no value has been given for'                // &
               Trim(ListAsChar)                                // &
               Trim(iAsChar)                                   // &
               ' "'                                            // &
               Trim(ColumnKey)                                 // &
               '".',                                              &
               3                                                  &
             )
      Else
        Call Message(                              &
               'FATAL ERROR in reading item "'  // &
               Trim(Item)                       // &
               '" from a block with keyword "'  // &
               Trim(BlockKey)                   // &
               '": no value has been given for' // &
               Trim(ListAsChar)                 // &
               Trim(iAsChar)                    // &
               ' "'                             // &
               Trim(ColumnKey)                  // &
               '".',                               &
               3                                   &
             )
      End If
    Else
      If (Item == ' ') Then
        Call Message(                                             &
               'FATAL ERROR in reading a block with keyword "' // &
               Trim(BlockKey)                                  // &
               '": the value "'                                // &
               Trim(C)                                         // &
               '" given for'                                   // &
               Trim(ListAsChar)                                // &
               Trim(iAsChar)                                   // &
               ' "'                                            // &
               Trim(ColumnKey)                                 // &
               '" is not an appropriately formatted time.',       &
               3                                                  &
             )
      Else
        Call Message(                                          &
               'FATAL ERROR in reading item "'              // &
               Trim(Item)                                   // &
               '" from a block with keyword "'              // &
               Trim(BlockKey)                               // &
               '": the value "'                             // &
               Trim(C)                                      // &
               '" given for'                                // &
               Trim(ListAsChar)                             // &
               Trim(iAsChar)                                // &
               ' "'                                         // &
               Trim(ColumnKey)                              // &
               '" is not an appropriately formatted time.',    &
               3                                               &
             )
      End If
    End If

  End If

End Function Token2Time

!-------------------------------------------------------------------------------------------------------------

Function InitArrays() Result (Arrays)
! Initialises a collection of input arrays.

  Implicit None
  ! Function result:
  Type(Arrays_) :: Arrays ! Initialised collection of input arrays.

  Arrays%nArrays = 0

End Function InitArrays

!-------------------------------------------------------------------------------------------------------------

Function InitArray(Name, ArrayElements, BlockKey) Result(Array)
! Initialises an input array.

! Note there are two types of error messages determined by BlockKey.
! BlockKey = Array : array assumed input in such a block
! BlockKey = blank : array assumed internally generated.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)  :: Name             ! Name of input array.
  Character(*), Intent(In)  :: ArrayElements(:) ! Elements of input array.
  Character(*), Intent(In)  :: BlockKey         ! Input block keyword. Blank if internally generated array.
  ! Function result:
  Type(Array_) :: Array ! Initialised input array.
  ! Locals:
  Integer :: i ! Loop index.

  Call BlockNameLengthTest(BlockName = Name, Length = MaxCharLength, BlockKey = BlockKey)

  If (Size(ArrayElements) > MaxArrayLength) Then
    If (BlockKey .CIEq. 'Array') Then
      Call Message('FATAL ERROR in reading an array: array "' // Trim(Name) // '" is too long.', 3)
    Else
      Call Message('UNEXPECTED FATAL ERROR in InitArray.', 4)
    End If
  End If

  Do i = 1, Size(ArrayElements)
    Call TokenLengthTest(                &
           C         = ArrayElements(i), &
           Length    = MaxCharLength,    &
           Zero      = .true.,           &
           BlockKey  = BlockKey,         &
           Item      = Name,             &
           ColumnKey = 'Array Values',   &
           i         = i                 &
         )
  End Do

  Array%Name                         = Name
  Array%n                            = Size(ArrayElements)
  Array%Array(1:Size(ArrayElements)) = ArrayElements(1:Size(ArrayElements))

End Function InitArray

!-------------------------------------------------------------------------------------------------------------

Subroutine AddArray(Array, Arrays)
! Adds an input array to a collection of input arrays.

  Implicit None
  ! Argument list:
  Type(Array_),  Intent(In)    :: Array  ! Input array.
  Type(Arrays_), Intent(InOut) :: Arrays ! Collection of input arrays.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Arrays%nArrays
    If (Array%Name .CIEq. Arrays%Arrays(i)%Name) Then
      If (Array == Arrays%Arrays(i)) Then
        Return
      Else
        Call Message(                                                            &
               'FATAL ERROR in adding the input array "'                      // &
               Trim(Array%Name)                                               // &
               '": a different input array with the same name already exists',   &
               3                                                                 &
             )
      End If
    End If
  End Do

  If (Arrays%nArrays >= MaxArrays) Then
    Call Message(                                       &
           'FATAL ERROR in adding the input array "' // &
           Trim(Array%Name)                          // &
           '": there are too many input arrays',        &
           3                                            &
         )
  End If

  Arrays%nArrays                = Arrays%nArrays + 1
  Arrays%Arrays(Arrays%nArrays) = Array

End Subroutine AddArray

!-------------------------------------------------------------------------------------------------------------

Function FindArrayIndex(Name, Arrays)
! Finds the index of an input array.

  Implicit None
  ! Argument list:
  Character(*),  Intent(In) :: Name   ! Name of input array.
  Type(Arrays_), Intent(In) :: Arrays ! Collection of input arrays.
  ! Function result:
  Integer :: FindArrayIndex ! Index of the input array.
  ! Locals:
  Integer :: i ! Loop index.

  Do i = 1, Arrays%nArrays
    If (Name .CIEq. Arrays%Arrays(i)%Name) Then
      FindArrayIndex = i
      Return
    End If
  End Do

  Call Message(                          &
         'FATAL ERROR: input array "' // &
         Trim(Name)                   // &
         '" not found',                  &
         3                               &
       )

End Function FindArrayIndex

!-------------------------------------------------------------------------------------------------------------

Function ArrayEq(Array1, Array2)
! Tests for equality of input arrays.

  Implicit None
  ! Argument list:
  Type(Array_), Intent(In) :: Array1 !} The two input arrays.
  Type(Array_), Intent(In) :: Array2 !}
  ! Function result:
  Logical :: ArrayEq ! Indicates if arrays are equal.
  ! Locals:
  Integer :: i ! Loop index.

  ArrayEq = (Array1%Name .CIEq. Array2%Name) .and. &
             Array1%n      ==   Array2%n

  If (ArrayEq) Then
    Do i = 1, Array1%n
      ArrayEq = ArrayEq .and. (Array1%Array(i) .CIEq. Array2%Array(i))
    End Do
  End If

End Function ArrayEq

!-------------------------------------------------------------------------------------------------------------

Subroutine ParseListChar(C, Delim, BlockKey, Item, ColumnKey, i, nValues, Values)
! Parses a list of character values separated by some delimiter into a character array.

! Note that, if there are n delimiters, then there are usually n + 1 elements returned, even though some have
! zero length. An exception is when the last of the n + 1 elements is blank when it is ignored and only n
! elements are returned. Equivalently, if the last element is intentionally blank, it must be followed by a
! delimiter.

! $$ Allow list of possible delimiters? Allow option of giving list as Array_.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: C
  Character(1), Intent(In)           :: Delim
  Character(*), Intent(In)           :: BlockKey
  Character(*), Intent(In)           :: Item
  Character(*), Intent(In)           :: ColumnKey
  Integer,      Intent(In), Optional :: i
  Integer,      Intent(Out)          :: nValues
  Character(*), Intent(Out)          :: Values(:)
  ! C         :: Character string containing the list.
  ! Delim     :: List delimiter. Cannot be a blank space.
  ! BlockKey  :: Block keyword.
  ! Item      :: Name of item read (block name or name associated with an input line). If blank, any message
  !              will not refer to the item read.
  ! ColumnKey :: Column keyword.
  ! i         :: Indicates that this is the ith data item for the column keyword in the block. This is best
  !              used when reading a block with the TwoD option (otherwise the value of i won't be known) and
  !              with Item omitted or equal to the block name (if i is used and Item refers to an input line,
  !              then i and Item are both referring to the same thing which may make any message confusing).
  ! nValues   :: Number of list elements.
  ! Values    :: List elements.
  ! Locals:
  Integer :: j1 !} Position of delimiters for list element (j1 = 0 and j2 = Len(C) + 1 for first and last
  Integer :: j2 !} elements).
  Integer :: k1 !] Position of start and end of list element (k1 > k2 for zero length elements).
  Integer :: k2 !]
  Integer :: k  ! Temporary variable.

# ifdef ExtraxChecks
    If (Delim == ' ') Then
      Call Message('UNEXPECTED FATAL ERROR in ParseListChar.', 4)
    End If
# endif

  nValues = 0
  j2      = 0

  Do

    j1 = j2

    k1 = j1 + 1

    If (j1 == Len(C)) Then
      j2 = j1 + 1
    Else
      j2 = Scan(C(k1:), Delim)
      If (j2 == 0) j2 = Len(C(k1:)) + 1
      j2 = j1 + j2
    End If

    k2 = j2 - 1

    ! Left justify element.
    If (k1 < k2) Then
      k = Verify(C(k1:k2), ' ')
      If (k /= 0) k1 = j1 + k
    End If

    ! Ignore last element if blank.
    If (j2 > Len(C)) Then
      If (k1 > k2) Exit
      If (C(k1:k2) == ' ') Exit
    End If

    nValues = nValues + 1
    If (nValues > Size(Values)) Then
      Call Message(                                                                &
             'FATAL ERROR reading input. Too many elements in the input list "' // &
             Trim(C)                                                            // &
             '". Number allowed is '                                            // &
             Int2Char(Size(Values)),                                               &
             3                                                                     &
           )
    End If

    If (k1 > k2) Then
      Values(nValues) = ' '
    Else
      Call TokenLengthTest(              &
             C         = C(k1:k2),       &
             Length    = Len(Values(1)), &
             Zero      = .false.,        &
             BlockKey  = BlockKey,       &
             Item      = Item,           &
             ColumnKey = ColumnKey,      &
             List      = .true.,         &
             i         = i               &
           )
      Values(nValues) = C(k1:k2)
    End If

    If (j2 > Len(C)) Exit

  End Do

End Subroutine ParseListChar

!-------------------------------------------------------------------------------------------------------------

Subroutine ParseListInt(C, Delim, BlockKey, Item, ColumnKey, i, nValues, Values)
! Parses a list of integer values separated by some delimiter into an integer array.

! Note that, if there are n delimiters, then there are usually n + 1 elements returned, even though some have
! zero length. An exception is when the last of the n + 1 elements is blank when it is ignored and only n
! elements are returned. Equivalently, if the last element is intentionally blank, it must be followed by a
! delimiter.

! $$ Allow list of possible delimiters? Allow option of giving list as Array_.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: C
  Character(1), Intent(In)           :: Delim
  Character(*), Intent(In)           :: BlockKey
  Character(*), Intent(In)           :: Item
  Character(*), Intent(In)           :: ColumnKey
  Integer,      Intent(In), Optional :: i
  Integer,      Intent(Out)          :: nValues
  Integer,      Intent(Out)          :: Values(:)
  ! C         :: Character string containing the list.
  ! Delim     :: List delimiter. Cannot be a blank space.
  ! BlockKey  :: Block keyword.
  ! Item      :: Name of item read (block name or name associated with an input line). If blank, any message
  !              will not refer to the item read.
  ! ColumnKey :: Column keyword.
  ! i         :: Indicates that this is the ith data item for the column keyword in the block. This is best
  !              used when reading a block with the TwoD option (otherwise the value of i won't be known) and
  !              with Item omitted or equal to the block name (if i is used and Item refers to an input line,
  !              then i and Item are both referring to the same thing which may make any message confusing).
  ! nValues   :: Number of list elements.
  ! Values    :: List elements.
  ! Locals:
  Integer :: j1 !} Position of delimiters for list element (j1 = 0 and j2 = Len(C) + 1 for first and last
  Integer :: j2 !} elements).
  Integer :: k1 !] Position of start and end of list element (k1 > k2 for zero length elements).
  Integer :: k2 !]
  Integer :: k  ! Temporary variable.

# ifdef ExtraxChecks
    If (Delim == ' ') Then
      Call Message('UNEXPECTED FATAL ERROR in ParseListChar.', 4)
    End If
# endif

  nValues = 0
  j2      = 0

  Do

    j1 = j2

    k1 = j1 + 1

    If (j1 == Len(C)) Then
      j2 = j1 + 1
    Else
      j2 = Scan(C(k1:), Delim)
      If (j2 == 0) j2 = Len(C(k1:)) + 1
      j2 = j1 + j2
    End If

    k2 = j2 - 1

    ! Left justify element.
    If (k1 < k2) Then
      k = Verify(C(k1:k2), ' ')
      If (k /= 0) k1 = j1 + k
    End If

    ! Ignore last element if blank.
    If (j2 > Len(C)) Then
      If (k1 > k2) Exit
      If (C(k1:k2) == ' ') Exit
    End If

    nValues = nValues + 1
    If (nValues > Size(Values)) Then
      Call Message(                                                                &
             'FATAL ERROR reading input. Too many elements in the input list "' // &
             Trim(C)                                                            // &
             '". Number allowed is '                                            // &
             Int2Char(Size(Values)),                                               &
             3                                                                     &
           )
    End If

    If (k1 > k2) Then
      Values(nValues) = Token2Int(               &
                          C         = ' ',       &
                          BlockKey  = BlockKey,  &
                          Item      = Item,      &
                          ColumnKey = ColumnKey, &
                          List      = .true.,    &
                          i         = i          &
                        )
    Else
      Values(nValues) = Token2Int(               &
                          C         = C(k1:k2),  &
                          BlockKey  = BlockKey,  &
                          Item      = Item,      &
                          ColumnKey = ColumnKey, &
                          List      = .true.,    &
                          i         = i          &
                        )
    End If

    If (j2 > Len(C)) Exit

  End Do

End Subroutine ParseListInt

!-------------------------------------------------------------------------------------------------------------

Subroutine ParseListStd(C, Delim, BlockKey, Item, ColumnKey, i, nValues, Values)
! Parses a list of real values separated by some delimiter into a real array.

! Note that, if there are n delimiters, then there are usually n + 1 elements returned, even though some have
! zero length. An exception is when the last of the n + 1 elements is blank when it is ignored and only n
! elements are returned. Equivalently, if the last element is intentionally blank, it must be followed by a
! delimiter.

! $$ Allow list of possible delimiters? Allow option of giving list as Array_.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: C
  Character(1), Intent(In)           :: Delim
  Character(*), Intent(In)           :: BlockKey
  Character(*), Intent(In)           :: Item
  Character(*), Intent(In)           :: ColumnKey
  Integer,      Intent(In), Optional :: i
  Integer,      Intent(Out)          :: nValues
  Real(Std),    Intent(Out)          :: Values(:)
  ! C         :: Character string containing the list.
  ! Delim     :: List delimiter. Cannot be a blank space.
  ! BlockKey  :: Block keyword.
  ! Item      :: Name of item read (block name or name associated with an input line). If blank, any message
  !              will not refer to the item read.
  ! ColumnKey :: Column keyword.
  ! i         :: Indicates that this is the ith data item for the column keyword in the block. This is best
  !              used when reading a block with the TwoD option (otherwise the value of i won't be known) and
  !              with Item omitted or equal to the block name (if i is used and Item refers to an input line,
  !              then i and Item are both referring to the same thing which may make any message confusing).
  ! nValues   :: Number of list elements.
  ! Values    :: List elements.
  ! Locals:
  Integer :: j1 !} Position of delimiters for list element (j1 = 0 and j2 = Len(C) + 1 for first and last
  Integer :: j2 !} elements).
  Integer :: k1 !] Position of start and end of list element (k1 > k2 for zero length elements).
  Integer :: k2 !]
  Integer :: k  ! Temporary variable.

# ifdef ExtraxChecks
    If (Delim == ' ') Then
      Call Message('UNEXPECTED FATAL ERROR in ParseListChar.', 4)
    End If
# endif

  nValues = 0
  j2      = 0

  Do

    j1 = j2

    k1 = j1 + 1

    If (j1 == Len(C)) Then
      j2 = j1 + 1
    Else
      j2 = Scan(C(k1:), Delim)
      If (j2 == 0) j2 = Len(C(k1:)) + 1
      j2 = j1 + j2
    End If

    k2 = j2 - 1

    ! Left justify element.
    If (k1 < k2) Then
      k = Verify(C(k1:k2), ' ')
      If (k /= 0) k1 = j1 + k
    End If

    ! Ignore last element if blank.
    If (j2 > Len(C)) Then
      If (k1 > k2) Exit
      If (C(k1:k2) == ' ') Exit
    End If

    nValues = nValues + 1
    If (nValues > Size(Values)) Then
      Call Message(                                                                &
             'FATAL ERROR reading input. Too many elements in the input list "' // &
             Trim(C)                                                            // &
             '". Number allowed is '                                            // &
             Int2Char(Size(Values)),                                               &
             3                                                                     &
           )
    End If

    ! $$ remove in due course.
    If (k1 > k2 .or. C(k1:Max(k1,k2)) == ' ') Then
      Call Message(                                                                                &
             'FATAL ERROR: This is probably because the input quantity '                        // &
             Trim(ColumnKey)                                                                    // &
             ' now (v5.4) needs to be a semi-colon separated list not a space separated list.',    &
             3                                                                                     &
           )
    End If

    If (k1 > k2) Then
      Values(nValues) = Token2Std(               &
                          C         = ' ',       &
                          BlockKey  = BlockKey,  &
                          Item      = Item,      &
                          ColumnKey = ColumnKey, &
                          List      = .true.,    &
                          i         = i          &
                        )
    Else
      Values(nValues) = Token2Std(               &
                          C         = C(k1:k2),  &
                          BlockKey  = BlockKey,  &
                          Item      = Item,      &
                          ColumnKey = ColumnKey, &
                          List      = .true.,    &
                          i         = i          &
                        )
    End If

    If (j2 > Len(C)) Exit

  End Do

End Subroutine ParseListStd

!-------------------------------------------------------------------------------------------------------------

Subroutine ParseListLog(C, Delim, BlockKey, Item, ColumnKey, i, nValues, Values)
! Parses a list of logical values separated by some delimiter into a logical array.

! Note that, if there are n delimiters, then there are usually n + 1 elements returned, even though some have
! zero length. An exception is when the last of the n + 1 elements is blank when it is ignored and only n
! elements are returned. Equivalently, if the last element is intentionally blank, it must be followed by a
! delimiter.

! $$ Allow list of possible delimiters? Allow option of giving list as Array_.

  Implicit None
  ! Argument list:
  Character(*), Intent(In)           :: C
  Character(1), Intent(In)           :: Delim
  Character(*), Intent(In)           :: BlockKey
  Character(*), Intent(In)           :: Item
  Character(*), Intent(In)           :: ColumnKey
  Integer,      Intent(In), Optional :: i
  Integer,      Intent(Out)          :: nValues
  Logical,      Intent(Out)          :: Values(:)
  ! C         :: Character string containing the list.
  ! Delim     :: List delimiter. Cannot be a blank space.
  ! BlockKey  :: Block keyword.
  ! Item      :: Name of item read (block name or name associated with an input line). If blank, any message
  !              will not refer to the item read.
  ! ColumnKey :: Column keyword.
  ! i         :: Indicates that this is the ith data item for the column keyword in the block. This is best
  !              used when reading a block with the TwoD option (otherwise the value of i won't be known) and
  !              with Item omitted or equal to the block name (if i is used and Item refers to an input line,
  !              then i and Item are both referring to the same thing which may make any message confusing).
  ! nValues   :: Number of list elements.
  ! Values    :: List elements.
  ! Locals:
  Integer :: j1 !} Position of delimiters for list element (j1 = 0 and j2 = Len(C) + 1 for first and last
  Integer :: j2 !} elements).
  Integer :: k1 !] Position of start and end of list element (k1 > k2 for zero length elements).
  Integer :: k2 !]
  Integer :: k  ! Temporary variable.

# ifdef ExtraxChecks
    If (Delim == ' ') Then
      Call Message('UNEXPECTED FATAL ERROR in ParseListChar.', 4)
    End If
# endif

  nValues = 0
  j2      = 0

  Do

    j1 = j2

    k1 = j1 + 1

    If (j1 == Len(C)) Then
      j2 = j1 + 1
    Else
      j2 = Scan(C(k1:), Delim)
      If (j2 == 0) j2 = Len(C(k1:)) + 1
      j2 = j1 + j2
    End If

    k2 = j2 - 1

    ! Left justify element.
    If (k1 < k2) Then
      k = Verify(C(k1:k2), ' ')
      If (k /= 0) k1 = j1 + k
    End If

    ! Ignore last element if blank.
    If (j2 > Len(C)) Then
      If (k1 > k2) Exit
      If (C(k1:k2) == ' ') Exit
    End If

    nValues = nValues + 1
    If (nValues > Size(Values)) Then
      Call Message(                                                                &
             'FATAL ERROR reading input. Too many elements in the input list "' // &
             Trim(C)                                                            // &
             '". Number allowed is '                                            // &
             Int2Char(Size(Values)),                                               &
             3                                                                     &
           )
    End If

    If (k1 > k2) Then
      Values(nValues) = Token2Log(               &
                          C         = ' ',       &
                          BlockKey  = BlockKey,  &
                          Item      = Item,      &
                          ColumnKey = ColumnKey, &
                          List      = .true.,    &
                          i         = i          &
                        )
    Else
      Values(nValues) = Token2Log(               &
                          C         = C(k1:k2),  &
                          BlockKey  = BlockKey,  &
                          Item      = Item,      &
                          ColumnKey = ColumnKey, &
                          List      = .true.,    &
                          i         = i          &
                        )
    End If

    If (j2 > Len(C)) Exit

  End Do

End Subroutine ParseListLog

!-------------------------------------------------------------------------------------------------------------

End Module HeadedFileModule
