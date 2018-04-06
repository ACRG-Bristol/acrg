!
!  Description: how to get values using keys.
!
!  Author: Enrico Fucile <enrico.fucile@ecmwf.int>
!
!  Copyright: See COPYING file that comes with this distribution
!
!
program read_message
use grib_api
implicit none
  integer                                ::  ifile,ofile
  integer                                ::  iret,igrib
  integer , dimension(50000)   ::  buffer
  integer                                ::  len,step,level

  call grib_open_file(ifile,'../../data/index.grib','r')
  call grib_open_file(ofile,'out.grib','w')
 
! a grib message is read from file into buffer
  len=size(buffer)*4
  call  grib_read_from_file(ifile,buffer,len,iret) 

  do while (iret/=GRIB_END_OF_FILE)

!   a new grib message is created from buffer
    call grib_new_from_message(igrib,buffer)

!   get as a integer
    call grib_get(igrib,'step', step) 
    write(*,*) 'step=',step

    call grib_get(igrib,'level',level)
    write(*,*) 'level=',level

    call grib_release(igrib)

    call grib_write_bytes(ofile,buffer,len)

!   a grib message is read from file into buffer
    len=size(buffer)*4
    call  grib_read_from_file(ifile,buffer,len,iret) 

  enddo

  call grib_close_file(ifile)
  call grib_close_file(ofile)

end program read_message
