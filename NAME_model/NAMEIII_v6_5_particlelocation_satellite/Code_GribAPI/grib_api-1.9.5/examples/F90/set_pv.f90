! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: how to set pv values.
!
!
!  Author: Anne Fouilloux
!
!
program set_pv
  use grib_api
  implicit none
  integer                         :: numberOfLevels
  integer                         :: numberOfCoefficients
  integer                         :: outfile, igrib
  integer                         :: i, ios
  real, dimension(:),allocatable  :: pv
  
  numberOfLevels=60
  numberOfCoefficients=2*(numberOfLevels+1)

  allocate(pv(numberOfCoefficients))

  ! read the model level coefficients from file
  open( unit=1, file="../../data/60_model_levels", &
                form="formatted",action="read")

  do i=1,numberOfCoefficients,2
     read(unit=1,fmt=*, iostat=ios) pv(i), pv(i+1)
     if (ios /= 0) then
        print *, "I/O error: ",ios
        exit
     end if
  end do
  
  ! print coefficients
  !do i=1,numberOfCoefficients,2
  !  print *,"  a=",pv(i)," b=",pv(i+1)
  !end do

  close(unit=1)

  call grib_open_file(outfile, 'out.grib1','w')
  
  !     a new grib message is loaded from file
  !     igrib is the grib id to be used in subsequent calls
  call grib_new_from_samples(igrib, "reduced_gg_sfc_grib1")

  !     set levtype to ml (model level)
  call grib_set(igrib,'levtype','ml')

  !     set level 
  call grib_set(igrib,'level',2)

  !     set PVPresent as an integer 
  call grib_set(igrib,'PVPresent',1)
  
  call grib_set(igrib,'pv',pv)
  
  !     write modified message to a file
  call grib_write(igrib,outfile)
  
  !  FREE MEMORY
  call grib_release(igrib)
  deallocate(pv)

  call grib_close_file(outfile)
  
end program set_pv
