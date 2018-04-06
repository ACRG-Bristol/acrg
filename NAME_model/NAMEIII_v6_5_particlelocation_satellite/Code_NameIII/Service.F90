! Module: Service Module

Module ServiceModule

! This module provides a route through which basic support services from a number of modules are made
! available.

!-------------------------------------------------------------------------------------------------------------

Use GlobalParametersModule
Use ErrorAndMessageModule
Use MathsModule, Gauss => GaussB ! Select GaussA or GaussB Gassian random number generator.
Use PhysicsModule
Use StringModule
Use SystemModule
Use ScreenModule
Use UnitModule
Use HeadedFileModule
Use ErrorAndMessageIIModule
Use TimeModule
Use CoordinateSystemModule
Use GridAndDomainModule
Use SortModule
Use PhysicalUnitsModule

!-------------------------------------------------------------------------------------------------------------

Implicit None

!-------------------------------------------------------------------------------------------------------------

Public

!-------------------------------------------------------------------------------------------------------------

End Module ServiceModule
