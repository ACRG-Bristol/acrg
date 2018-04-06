! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: how to copy a namespace from a message to another.
!               
!
program copy_namespace
  use grib_api
  implicit none  
  integer  :: file1, file2, file3
  integer  :: igrib1,igrib2,igrib3

  call grib_open_file(file1,'../../data/reduced_latlon_surface.grib2','r')
  call grib_open_file(file2,'../../data/regular_latlon_surface.grib1','r')
  call grib_open_file(file3,'out.grib','w')

  call grib_new_from_file(file1,igrib1)
  call grib_new_from_file(file2,igrib2)

  call grib_clone(igrib2,igrib3)

  call grib_copy_namespace(igrib1,'geography',igrib3)

   call grib_write(igrib3,file3)

  call grib_close_file(file1)
  call grib_close_file(file2)
  call grib_close_file(file3)

end program copy_namespace
