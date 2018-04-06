! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: how to set missing a key value.
!
!
!  Author: Enrico Fucile
!
!
!
program set
  use grib_api
  implicit none
  integer              :: infile,outfile
  integer              :: igrib

  infile=5
  outfile=6

  call grib_open_file(infile, &
       '../../data/reduced_gaussian_pressure_level.grib2','r')

  call grib_open_file(outfile, &
       'out_surface_level.grib2','w')

  !     a new grib message is loaded from file
  !     igrib is the grib id to be used in subsequent calls
  call grib_new_from_file(infile,igrib)

  call grib_set(igrib,'typeOfFirstFixedSurface','sfc')
  call grib_set_missing(igrib,'scaleFactorOfFirstFixedSurface')
  call grib_set_missing(igrib,'scaledValueOfFirstFixedSurface')

  call grib_write(igrib,outfile)

  call grib_release(igrib)

  call grib_close_file(infile)

  call grib_close_file(outfile)

end program set
