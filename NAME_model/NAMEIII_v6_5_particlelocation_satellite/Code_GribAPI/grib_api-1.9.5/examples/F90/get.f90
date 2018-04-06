! Copyright 2005-2007 ECMWF
! 
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: how to get values using keys.
!
!  Author: Enrico Fucile 
!
!
program get
  use grib_api
  implicit none
  
  integer                            ::  ifile
  integer                            ::  iret
  integer                            ::  igrib
  real                               ::  latitudeOfFirstPointInDegrees
  real                               ::  longitudeOfFirstPointInDegrees
  real                               ::  latitudeOfLastPointInDegrees
  real                               ::  longitudeOfLastPointInDegrees
  integer                            ::  numberOfPointsAlongAParallel
  integer                            ::  numberOfPointsAlongAMeridian
  real, dimension(:), allocatable    ::  values
  integer                            ::  numberOfValues
  real                               ::  average,min_val, max_val
  integer                            ::  is_missing
  
  call grib_open_file(ifile, &
       '../../data/reduced_latlon_surface.grib1','r')
  
  ! Loop on all the messages in a file.
  
  !     a new grib message is loaded from file
  !     igrib is the grib id to be used in subsequent calls
  call  grib_new_from_file(ifile,igrib, iret) 
  
  LOOP: DO WHILE (iret /= GRIB_END_OF_FILE)

     !check if the value of the key is MISSING
     is_missing=0;
     call grib_is_missing(igrib,'numberOfPointsAlongAParallel', &
          is_missing);
     if ( is_missing /= 1 ) then
        !     get as a integer
        call grib_get(igrib,'numberOfPointsAlongAParallel', &
             numberOfPointsAlongAParallel) 
        write(*,*) 'numberOfPointsAlongAParallel=', &
             numberOfPointsAlongAParallel
     else
        write(*,*) 'numberOfPointsAlongAParallel is missing'
     endif     
     !     get as a integer
     call grib_get(igrib,'numberOfPointsAlongAMeridian', &
          numberOfPointsAlongAMeridian) 
     write(*,*) 'numberOfPointsAlongAMeridian=', &
          numberOfPointsAlongAMeridian
     
     !     get as a real
     call grib_get(igrib, 'latitudeOfFirstGridPointInDegrees', &
          latitudeOfFirstPointInDegrees) 
     write(*,*) 'latitudeOfFirstGridPointInDegrees=', &
          latitudeOfFirstPointInDegrees
     
     !     get as a real
     call grib_get(igrib, 'longitudeOfFirstGridPointInDegrees', &
          longitudeOfFirstPointInDegrees) 
     write(*,*) 'longitudeOfFirstGridPointInDegrees=', &
          longitudeOfFirstPointInDegrees
     
     !     get as a real
     call grib_get(igrib, 'latitudeOfLastGridPointInDegrees', &
          latitudeOfLastPointInDegrees) 
     write(*,*) 'latitudeOfLastGridPointInDegrees=', &
          latitudeOfLastPointInDegrees
     
     !     get as a real
     call grib_get(igrib, 'longitudeOfLastGridPointInDegrees', &
          longitudeOfLastPointInDegrees) 
     write(*,*) 'longitudeOfLastGridPointInDegrees=', &
          longitudeOfLastPointInDegrees
     
     
     !     get the size of the values array
     call grib_get_size(igrib,'values',numberOfValues)
     write(*,*) 'numberOfValues=',numberOfValues
     
     allocate(values(numberOfValues), stat=iret)
     !     get data values
     call grib_get(igrib,'values',values)
     call grib_get(igrib,'min',min_val) ! can also be obtained through minval(values)
     call grib_get(igrib,'max',max_val) ! can also be obtained through maxval(values)
     call grib_get(igrib,'average',average) ! can also be obtained through maxval(values)

     deallocate(values)
          
     write(*,*)'There are ',numberOfValues, &
          ' average is ',average, &
          ' min is ',  min_val, &
          ' max is ',  max_val
     
     call grib_release(igrib)
     
     call grib_new_from_file(ifile,igrib, iret)
     
  end do LOOP
  
  call grib_close_file(ifile)
  
end program get
