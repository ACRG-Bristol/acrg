! Module:  Compiler Information Module

Module CompilerInfoModule

! This module provides Compiler information

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Public

#ifdef COMPILETIME
  Character(200) :: CompileTime=COMPILETIME
#else
  Character(200) :: CompileTime="unknown"
#endif

#ifdef COMPILERVERSION
  Character(200) :: CompilerVersion=COMPILERVERSION
#else
  Character(200) :: CompilerVersion="unknown"
#endif

!-------------------------------------------------------------------------------------------------------------

End Module CompilerInfoModule
