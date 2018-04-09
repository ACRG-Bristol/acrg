! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: how to use grib_find_nearest and grib_get_element 
!
!
!  Author: Enrico Fucile 
!
!
!
program find
  use grib_api
  implicit none
  integer                                      :: npoints
  integer                                      :: infile
  integer                                      :: igrib, ios, i
  real(8), dimension(:), allocatable  :: lats, lons
  real(8), dimension(:), allocatable  :: nearest_lats, nearest_lons
  real(8), dimension(:), allocatable  :: distances, values, lsm_values
  integer(kind=kindOfInt), dimension(:), allocatable  :: indexes
  real(kind=8)                        :: value

! initialization
  open( unit=1, file="../../data/list_points",form="formatted",action="read")
  read(unit=1,fmt=*) npoints
  allocate(lats(npoints))
  allocate(lons(npoints))
  allocate(nearest_lats(npoints))
  allocate(nearest_lons(npoints))
  allocate(distances(npoints))
  allocate(lsm_values(npoints))
  allocate(values(npoints))
  allocate(indexes(npoints))
  do i=1,npoints
     read(unit=1,fmt=*, iostat=ios) lats(i), lons(i)
     if (ios /= 0) then
        npoints = i - 1
        exit
     end if
  end do
  close(unit=1)
  call grib_open_file(infile, &
       '../../data/reduced_gaussian_lsm.grib1','r')
  
  !     a new grib message is loaded from file
  !     igrib is the grib id to be used in subsequent calls
  call grib_new_from_file(infile,igrib)
  

  call grib_find_nearest(igrib, .true., lats, lons, nearest_lats, nearest_lons,lsm_values, distances, indexes)
  call grib_release(igrib)
  
  call grib_close_file(infile)

! will apply it to another GRIB
  call grib_open_file(infile, &
       '../../data/reduced_gaussian_pressure_level.grib1','r')
  call grib_new_from_file(infile,igrib)

  call grib_get_element(igrib,"values", indexes, values)
  call grib_release(igrib)
  call grib_close_file(infile)

  do i=1, npoints
     print*,lats(i), lons(i), nearest_lats(i), nearest_lons(i), distances(i), lsm_values(i), values(i)
  end do

  deallocate(lats)
  deallocate(lons)
  deallocate(nearest_lats)
  deallocate(nearest_lons)
  deallocate(distances)
  deallocate(lsm_values)
  deallocate(values)
  deallocate(indexes)

end program find
