! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: how to set key values.
!
!
!  Author: Anne Fouilloux                            
!
!
program set
  use grib_api
  implicit none
  integer(kind = 4)    :: centre, date
  integer              :: infile,outfile
  integer              :: igrib


  centre = 80
  call current_date(date)
  call grib_open_file(infile, &
       '../../data/regular_latlon_surface_constant.grib1','r')

  call grib_open_file(outfile, &
       'out.grib1','w')

  !     a new grib message is loaded from file
  !     igrib is the grib id to be used in subsequent calls
  call grib_new_from_file(infile,igrib)

  call grib_set(igrib,'date',date)
  !     set centre as a integer */
  call grib_set(igrib,'centre',centre)

! check if it is correct in the actual GRIB message

  call check_settings(igrib)

  !     write modified message to a file
  call grib_write(igrib,outfile)

  call grib_release(igrib)

  call grib_close_file(infile)

  call grib_close_file(outfile)

contains

!======================================
subroutine current_date(date)
integer, intent(out) :: date

integer              :: val_date(8)
call date_and_time ( values = val_date)

date = val_date(1)* 10000 + val_date(2)*100 + val_date(3) 
end subroutine current_date
!======================================
subroutine check_settings(gribid)
  use grib_api
  implicit none
  integer, intent(in) :: gribid
  
  integer(kind = 4)    :: int_value
  character(len = 10)  :: string_value

  !     get centre as a integer
  call grib_get(gribid,'centre',int_value)
  write(*,*) "get centre as a integer - centre = ",int_value
  
  !     get centre as a string
  call grib_get(gribid,'centre',string_value)
  write(*,*) "get centre as a string  - centre = ",string_value
  
  !     get date as a string
  call grib_get(gribid,'dataDate',string_value)
  write(*,*) "get date as a string    - date = ",string_value
  
end subroutine check_settings
end program set
