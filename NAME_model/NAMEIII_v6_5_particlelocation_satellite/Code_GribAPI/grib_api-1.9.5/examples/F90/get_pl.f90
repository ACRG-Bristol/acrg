! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: how to get PL values.
!
!
!  Author: Anne Fouilloux
!
!
program get_pl
  use grib_api
  implicit none
  integer                         :: infile
  integer                         :: igrib
  integer                         :: PLPresent, nb_pl
  real, dimension(:), allocatable :: pl


  call grib_open_file(infile, &
       '../../data/reduced_gaussian_surface.grib1','r')
  
  !     a new grib message is loaded from file
  !     igrib is the grib id to be used in subsequent calls
  call grib_new_from_file(infile,igrib)
  
  !     set PVPresent as an integer 
  call grib_get(igrib,'PLPresent',PLPresent)
  print*, "PLPresent= ", PLPresent
  if (PLPresent == 1) then
     call grib_get_size(igrib,'pl',nb_pl)
     print*, "there are ", nb_pl, " PL values"
     allocate(pl(nb_pl))
     call grib_get(igrib,'pl',pl)
     print*, "pl = ", pl
     deallocate(pl)
  else
     print*, "There is no PL values in your GRIB message!"
  end if
  call grib_release(igrib)
     
  call grib_close_file(infile)

end program get_pl
