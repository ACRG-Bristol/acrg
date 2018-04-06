! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: prints all the data contained in a grib file
!
!
!  Author: Anne Fouilloux
!
!
program print_data
use grib_api
implicit none
integer            :: ifile
integer            :: iret
integer            :: igrib
integer            :: i
real(kind=8), dimension(:), allocatable       :: values
integer(kind=4)    :: numberOfValues
real(kind=8)       :: average
real(kind=8)       :: max
real(kind=8)       :: min
character(len=256) :: error

call grib_open_file(ifile, &
           '../../data/constant_field.grib1','r')

!     a new grib message is loaded from file
!     igrib is the grib id to be used in subsequent calls
      call grib_new_from_file(ifile,igrib)


!     get the size of the values array
      call grib_get_size(igrib,'values',numberOfValues)

!     get data values
  print*, 'number of values ', numberOfValues
  allocate(values(numberOfValues), stat=iret)

  call grib_get(igrib,'values',values)

  do i=1,numberOfValues
    write(*,*)'  ',i,values(i)
  enddo


  write(*,*)numberOfValues,' values found '

  call grib_get(igrib,'max',max)
  write(*,*) 'max=',max
  call grib_get(igrib,'min',min)
  write(*,*) 'min=',min
  call grib_get(igrib,'average',average)
  write(*,*) 'average=',average

  call grib_release(igrib)

  call grib_close_file(ifile)

  deallocate(values)

end program print_data
