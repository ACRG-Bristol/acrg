! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: how to get lat/lon/values.
!
!
!  Author: Enrico Fucile
!
!
program get_data
use grib_api
implicit none
  integer            :: ifile
  integer            :: iret,i
  real(kind=8),dimension(:),allocatable     :: lats,lons,values
  integer(4)        :: numberOfPoints
  real(8)  :: missingValue=9999
  integer           :: count=0
  character(len=256) :: filename

!     Message identifier.
  integer            :: igrib

  ifile=5

  call grib_open_file(ifile, &
       '../../data/reduced_latlon_surface.grib1','R')

! Loop on all the messages in a file.

  call grib_new_from_file(ifile,igrib,iret)

  do while (iret/=GRIB_END_OF_FILE)
    count=count+1
    print *, "===== Message #",count
    call grib_get(igrib,'numberOfPoints',numberOfPoints)
    call grib_set(igrib,'missingValue',missingValue)

    allocate(lats(numberOfPoints))
    allocate(lons(numberOfPoints))
    allocate(values(numberOfPoints))

    call grib_get_data(igrib,lats,lons,values)

    do i=1,numberOfPoints
      if (values(i) /= missingValue) then
        print *, lats(i),lons(i),values(i)
      end if
    enddo

    deallocate(lats)
    deallocate(lons)
    deallocate(values)

    call grib_release(igrib)
    call grib_new_from_file(ifile,igrib, iret)

  end do 


  call grib_close_file(ifile)

end program 
