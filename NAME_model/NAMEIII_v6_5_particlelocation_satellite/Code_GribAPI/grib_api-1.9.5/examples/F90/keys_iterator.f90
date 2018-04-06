! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description:
!       How to use keys_iterator to get all the available
!       keys in a message.
!
!  Author: Enrico Fucile 
!
!
program keys_iterator
  use grib_api
  implicit none
  character(len=20)  :: name_space
  integer            :: kiter,ifile,igrib,iret
  character(len=256) :: key
  character(len=256) :: value
  character(len=512) :: all
  integer            :: grib_count
  
  call grib_open_file(ifile, &
       '../../data/regular_latlon_surface.grib1','r')
  
  ! Loop on all the messages in a file.
  
  call grib_new_from_file(ifile,igrib, iret)
  
  do while (iret /= GRIB_END_OF_FILE)

     grib_count=grib_count+1
     write(*,*) '-- GRIB N. ',grib_count,' --'
     
     ! valid name_spaces are ls and mars
     name_space='ls'
     
     call grib_keys_iterator_new(igrib,kiter,name_space)
     
     do
        call grib_keys_iterator_next(kiter, iret) 
        
        if (iret .ne. 1) exit
        
        call grib_keys_iterator_get_name(kiter,key)
        call grib_get(igrib,trim(key),value)
        all=trim(key)// ' = ' // trim(value)
        write(*,*) trim(all)
        
     end do
     
     call grib_keys_iterator_delete(kiter)
     call grib_release(igrib)
     call grib_new_from_file(ifile,igrib, iret)
  end do
  
  
  call grib_close_file(ifile)
  
end program keys_iterator

