! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: how to create a multi field message in memory and write
!               it in a file. The multi field messages can be created
!               only in grib edition 2.
!
!
!
program multi_write
  use grib_api
  implicit none
  integer                            :: infile,outfile
  integer                            :: in_gribid,iret
  integer                            :: multi_gribid
  integer                            :: step,startsection

  ! multi field messages can be created only in edition 2 
  call grib_open_file(infile,'../../data/sample.grib2','r')

  call grib_open_file(outfile,'multi_created.grib2','w')

  !     a grib message is loaded from file
  !     in_gribid is the grib id to be used in subsequent calls
  call grib_new_from_file(infile,in_gribid)

  startsection=4
  do step=0,240,12

    call grib_set(in_gribid,"step",step)
    call grib_multi_append(in_gribid,startsection,multi_gribid)

  enddo

!  write messages to a file
  call grib_multi_write(multi_gribid,outfile)

  call grib_close_file(infile)
  call grib_close_file(outfile)

end program multi_write
