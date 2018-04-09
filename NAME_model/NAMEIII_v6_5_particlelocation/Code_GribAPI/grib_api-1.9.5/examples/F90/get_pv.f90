! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: how to get PV values.
!
!
!  Author: Anne Fouilloux
!
!
program get_pv
  use grib_api
  implicit none
  integer                         :: infile
  integer                         :: igrib
  integer                         :: PVPresent, nb_pv
  real, dimension(:), allocatable :: pv
  

  call grib_open_file(infile, &
       '../../data/reduced_gaussian_model_level.grib1','r')
  
  !     a new grib message is loaded from file
  !     igrib is the grib id to be used in subsequent calls
  call grib_new_from_file(infile,igrib)
  
  !     set PVPresent as an integer 
  call grib_get(igrib,'PVPresent',PVPresent)
  print*, "PVPresent = ", PVPresent
  if (PVPresent == 1) then
     call grib_get_size(igrib,'pv',nb_pv)
     print*, "There are ", nb_pv, " PV values"
     allocate(pv(nb_pv))
     call grib_get(igrib,'pv',pv)
     print*, "pv = ", pv
     deallocate(pv)
  else
     print*, "There is no PV values in your GRIB message!"
  end if
  call grib_release(igrib)
  
  call grib_close_file(infile)
  
end program get_pv
