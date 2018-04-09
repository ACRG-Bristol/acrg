! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!
!  Author: Enrico Fucile
!
!
program read_from_file
use grib_api
  implicit none
  integer :: size,intsize
  parameter (intsize=5000000,size=intsize*4)
  integer            :: ifile
  integer            :: iret
  integer           :: count=0
  character(len=256) :: filename
  integer(kind=4),dimension(intsize) :: buffer
  !character,dimension(size) :: buffer
  integer :: len

!     Message identifier.
  integer            :: igrib

  ifile=5

  call grib_open_file(ifile,'../../data/collection.grib1','r')

! Loop on all the messages in a file.

  len=size
  call grib_read_from_file(ifile,buffer,len,iret)
 ! call grib_new_from_file(ifile,igrib, iret)

  do while (iret==GRIB_SUCCESS)
    count=count+1
    print *, '===== Message #',count,' len=',len
    len=size
    call grib_read_from_file(ifile,buffer,len,iret)
   ! call grib_new_from_file(ifile,igrib, iret)

  end do 

  if (iret/=GRIB_END_OF_FILE) then
    call grib_check(iret,'read_from_file','')
  endif 


  call grib_close_file(ifile)

end program 
