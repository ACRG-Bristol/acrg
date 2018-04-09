! Copyright 2005-2007 ECMWF
!
! Licensed under the GNU Lesser General Public License which
! incorporates the terms and conditions of version 3 of the GNU
! General Public License.
! See LICENSE and gpl-3.0.txt for details.
!
!
!  Description: count messages before processing
!
!  Author: Enrico Fucile
!
!
program get
  use grib_api
  implicit none

  integer                            ::  ifile
  integer                            ::  iret
  integer                            ::  n
  integer                            ::  i
  integer,dimension(:),allocatable   ::  igrib
  real                               ::  latitudeOfFirstPointInDegrees
  real                               ::  longitudeOfFirstPointInDegrees
  real                               ::  latitudeOfLastPointInDegrees
  real                               ::  longitudeOfLastPointInDegrees
  integer                            ::  numberOfPointsAlongAParallel
  integer                            ::  numberOfPointsAlongAMeridian
  real, dimension(:), allocatable    ::  values
  integer                            ::  numberOfValues
  real                               ::  average,min_val, max_val

  call grib_open_file(ifile, &
       '../../data/tigge_pf_ecmwf.grib2','r')

  ! count the messages in the file
  call grib_count_in_file(ifile,n)
  allocate(igrib(n))
  igrib=-1

  ! Load the messages from the file.
  DO i=1,n
     call grib_new_from_file(ifile,igrib(i), iret)
  END DO

  ! we can close the file
  call grib_close_file(ifile)

  ! Loop on all the messages in memory
  DO i=1,n
     write(*,*) 'processing message number ',i
     !     get as a integer
     call grib_get(igrib(i),'numberOfPointsAlongAParallel', &
          numberOfPointsAlongAParallel)
     write(*,*) 'numberOfPointsAlongAParallel=', &
          numberOfPointsAlongAParallel

     !     get as a integer
     call grib_get(igrib(i),'numberOfPointsAlongAMeridian', &
          numberOfPointsAlongAMeridian)
     write(*,*) 'numberOfPointsAlongAMeridian=', &
          numberOfPointsAlongAMeridian

     !     get as a real
     call grib_get(igrib(i), 'latitudeOfFirstGridPointInDegrees', &
          latitudeOfFirstPointInDegrees)
     write(*,*) 'latitudeOfFirstGridPointInDegrees=', &
          latitudeOfFirstPointInDegrees

     !     get as a real
     call grib_get(igrib(i), 'longitudeOfFirstGridPointInDegrees', &
          longitudeOfFirstPointInDegrees)
     write(*,*) 'longitudeOfFirstGridPointInDegrees=', &
          longitudeOfFirstPointInDegrees

     !     get as a real
     call grib_get(igrib(i), 'latitudeOfLastGridPointInDegrees', &
          latitudeOfLastPointInDegrees)
     write(*,*) 'latitudeOfLastGridPointInDegrees=', &
          latitudeOfLastPointInDegrees

     !     get as a real
     call grib_get(igrib(i), 'longitudeOfLastGridPointInDegrees', &
          longitudeOfLastPointInDegrees)
     write(*,*) 'longitudeOfLastGridPointInDegrees=', &
          longitudeOfLastPointInDegrees


     !     get the size of the values array
     call grib_get_size(igrib(i),'values',numberOfValues)
     write(*,*) 'numberOfValues=',numberOfValues

     allocate(values(numberOfValues), stat=iret)
     !     get data values
     call grib_get(igrib(i),'values',values)
     call grib_get(igrib(i),'min',min_val) ! can also be obtained through minval(values)
     call grib_get(igrib(i),'max',max_val) ! can also be obtained through maxval(values)
     call grib_get(igrib(i),'average',average) ! can also be obtained through maxval(values)

     write(*,*)'There are ',numberOfValues, &
          ' average is ',average, &
          ' min is ',  min_val, &
          ' max is ',  max_val
     write(*,*) '---------------------'
  END DO

  DO i=1,n
    call grib_release(igrib(i))
  END DO

  deallocate(values)
  deallocate(igrib)

end program get
