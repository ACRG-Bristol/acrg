! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: how to set a bitmap in a grib message 
!
!
!  Author: Enrico Fucile 
!
!
program set_bitmap
  use grib_api
  implicit none
  integer                         :: infile,outfile
  integer                         :: igrib, iret
  integer                         :: numberOfValues
  real, dimension(:), allocatable :: values
  real                            :: missingValue
  logical                         :: grib1Example

  grib1Example=.true.

  if (grib1Example) then
    ! GRIB 1 example
    call grib_open_file(infile,'../../data/regular_latlon_surface.grib1','r')
  else
    ! GRIB 2 example
    call grib_open_file(infile,'../../data/regular_latlon_surface.grib2','r')
  end if
  
  call grib_open_file(outfile,'out.grib','w')
  
  !     a new grib message is loaded from file
  !     igrib is the grib id to be used in subsequent calls
  call grib_new_from_file(infile,igrib)
  
  ! The missingValue is not coded in the message. 
  ! It is a value we define as a placeholder for a missing value
  ! in a point of the grid.
  ! It should be choosen in a way that it cannot be confused 
  ! with a valid field value
  missingValue=9999
  call grib_set(igrib, 'missingValue',missingValue)
  write(*,*) 'missingValue=',missingValue

  ! get the size of the values array
  call grib_get_size(igrib,'values',numberOfValues)
  write(*,*) 'numberOfValues=',numberOfValues
  
  allocate(values(numberOfValues), stat=iret)

  ! get data values
  call grib_get(igrib,'values',values)
  
  ! enable bitmap 
  call grib_set(igrib,"bitmapPresent",1)

  ! some values are missing
  values(1:10) = missingValue

  ! set the values (the bitmap will be automatically built)
  call grib_set(igrib,'values', values)

  !  write modified message to a file
  call grib_write(igrib,outfile)
  
  ! FREE MEMORY
  call grib_release(igrib)
  
  call grib_close_file(infile)

  call grib_close_file(outfile)

  deallocate(values)

end program set_bitmap
