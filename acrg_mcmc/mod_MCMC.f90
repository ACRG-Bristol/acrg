Module mod_MCMC

 Implicit none

 real,parameter             :: pi = 3.14159265       		! a fixed constant

 Contains

 subroutine pdf_calc(x,pdf_param1,pdf_param2,pdf,p1)      	! Subroutine for calculating P1 for specified PDF

 Implicit none
 ! Args
 integer,intent(in)   :: pdf
 real,intent(in)      :: x, pdf_param1, pdf_param2
 real,intent(out)     :: p1
 integer              :: positive
 ! Locals

!  PDF .eq. 'LOGNORMAL'
  If(pdf .eq. 0) then                            ! Lognormal PDF = 0
     p1=(1./(x*pdf_param2*sqrt(2.*pi)))*exp( -1.*(alog(x) - alog(pdf_param1))**2 / 2./(pdf_param2**2))
     positive=1
  Endif

!  PDF .eq. 'EXPONENTIAL'
  If(pdf .eq. 1) then                     	! Exponential PDF = 1
     p1=pdf_param1*exp(-1.*pdf_param1*x)
     positive=1
  Endif

!  PDF .eq. 'GAUSSIAN'
  If(pdf .eq. 2) then                   	! Gaussian PDF = 2
     p1=(1./(pdf_param2*sqrt(2.*pi)))*exp(-1.*(x - pdf_param1)**2 / 2./(pdf_param2**2))
     positive=0
  Endif

!  PDF .eq. 'UNIFORM'
  If(pdf .eq. 3) then                   	! Uniform PDF = 3	
	if (x .lt. pdf_param1 .OR. x .gt. pdf_param2) then
	p1 = 0
	else
    	p1= 1./(pdf_param2 - pdf_param1)
	endif
     	positive=0
  Endif

  if((x .le. 0.) .and. (positive .eq. 1)) p1 = 0.
    
 ! return p1

end subroutine pdf_calc
  
FUNCTION random_normal() RESULT(fn_val)
IMPLICIT NONE
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

REAL :: fn_val

!     Local variables
REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
            r1 = 0.27597, r2 = 0.27846, half =0.5, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
RETURN

END FUNCTION random_normal

 subroutine init_random_seed()
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid, t(2), s, getpid
            integer(8) :: count, tms
           
            call random_seed(size = n)
            allocate(seed(n))

            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(count)
               if (count /= 0) then
                  t = transfer(count, t)
               else
                  call date_and_time(values=dt)
                  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
               pid = getpid() + 1099279 ! Add a prime
               s = ieor(s, pid)
               if (n >= 3) then
                  seed(1) = t(1) + 36269
                  seed(2) = t(2) + 72551
                  seed(3) = pid
                  if (n > 3) then
                     seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  end if
               else
                  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
               end if
            end if
            call random_seed(put=seed)
          end subroutine init_random_seed

function Ainv(A) 
  real, dimension(:,:), intent(in) :: A
  real, dimension(size(A,1),size(A,2)) :: Ainv

  real, dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external SGETRF
  external SGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call SGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call SGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function Ainv


function det(A)
  !Returns the log of the "square root of the determinant" of a postive definite matrix calculated from the Cholesky decomposition
  ! Uses LAPACK routines
  real, dimension(:,:), intent(in) :: A
  real, dimension(size(A,1),size(A,2)) :: Adum
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  real :: det
  integer :: n, info,j

  ! External procedures defined in LAPACK
  external SPOTRF

  ! Store A in Adum to prevent it from being overwritten by LAPACK
  Adum = A
  n = size(A,1)

  call SPOTRF('U',n,Adum,n,info)
  det = 0

  do j = 1,n
  det = det + alog(Adum(j,j))
  enddo

end function det

subroutine rkbesl ( x, alpha, nb, ize, k_arg, ncalc )

!*****************************************************************************80
!
!! RKBESL calculates K Bessel function with non-integer orders.
!
!  Discussion:
!
!    This routine calculates modified Bessel functions of the second
!    kind, K SUB(N+ALPHA) (X), for non-negative argument X, and
!    non-negative order N+ALPHA, with or without exponential scaling.
!
!    This program is based on a program written by J. B. Campbell
!    that computes values of the Bessel functions K of real
!    argument and real order.  Modifications include the addition
!    of non-scaled functions, parameterization of machine
!    dependencies, and the use of more accurate approximations
!    for SINH and SIN.
!
!    In case of an error, NCALC .NE. NB, and not all K's are
!    calculated to the desired accuracy.
!
!    NCALC < -1:  An argument is out of range. For example,
!    NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
!    XMAX.  In this case, the B-vector is not calculated,
!    and NCALC is set to MIN0(NB,0)-2  so that NCALC .NE. NB.
!
!    NCALC = -1:  Either  K(ALPHA,X) .GE. XINF  or
!    K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) .GE. XINF.  In this case,
!    the B-vector is not calculated.  Note that again
!    NCALC .NE. NB.
!
!    0 < NCALC < NB: Not all requested function values could
!    be calculated accurately.  BK(I) contains correct function
!    values for I <= NCALC, and contains the ratios
!    K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!    22 January 2014 A. Ganesan - added bk(1:nb) = 0 (ln 396) after  'ncalc = min ( nb, 0 ) - 2'
!    16 January 2014 A. Ganesan - only the final order of bk is output k_arg = bk(nb)
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    JB Campbell,
!    On Temme's Algorithm for the Modified Bessel Functions of the
!    Third Kind,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 581-586.
!
!    JB Campbell,
!    A FORTRAN IV Subroutine for the Modified Bessel Functions of
!    the Third Kind of Real Order and Real Argument,
!    Report NRC/ERB-925,
!    National Research Council, Canada.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the non-negative argument for which
!    K's or exponentially scaled K's (K*EXP(X))
!    are to be calculated.  If K's are to be calculated,
!    X must not be greater than XMAX.
!
!    Input, real ( kind = 8 ) ALPHA, the fractional part of order for which
!    K's or exponentially scaled K's (K*EXP(X)) are to be calculated.
!    0 <= ALPHA < 1.0.
!
!    Input, integer ( kind = 4 ) NB, the number of functions to be calculated, 
!    NB .GT. 0.  The first function calculated is of order ALPHA, and the
!    last is of order (NB - 1 + ALPHA).
!
!    Input, integer ( kind = 4 ) IZE, scaling option.
!    1, unscaled functions are to calculated,
!    2, exponentially scaled functions are to be calculated.
!
!    Output, real ( kind = 8 ) BK(NB), the results.  If the routine
!    terminates normally, with NCALC = NB, the vector BK contains the
!    functions K(ALPHA,X), ... , K(NB-1+ALPHA,X), or the corresponding
!    exponentially scaled functions.
!    If (0 < NCALC < NB), BK(I) contains correct function
!    values for I <= NCALC, and contains the ratios
!    K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
!
!    Output, integer ( kind = 4 ) NCALC, error indicator.  If NCALC = NB, then 
!    all the requested values were calculated to the desired accuracy.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) blpha
  real ( kind = 8 ) bk(1)
  real ( kind = 8 ) bk1
  real ( kind = 8 ) bk2
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dm
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) d3
  real ( kind = 8 ) enu
  real ( kind = 8 ) estf(7)
  real ( kind = 8 ) estm(6)
  real ( kind = 8 ) ex
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) ize
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mplus1
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) ncalc
  real ( kind = 8 ) p(8)
  real ( kind = 8 ) p0
  real ( kind = 8 ) q(7)
  real ( kind = 8 ) q0
  real ( kind = 8 ) r(5)
  real ( kind = 8 ) ratio
  real ( kind = 8 ) s(4)
  real ( kind = 8 ) sqxmin
  real ( kind = 8 ) t(6)
  real ( kind = 8 ) tinyx
  real ( kind = 8 ) twonu
  real ( kind = 8 ) twox
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) wminf
  real ( kind = 8 ) x
  real ( kind = 8 ) xinf
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) x2by4
  real ( kind = 8 ) k_arg
!
!  Mathematical constants
!    A = LOG(2.D0) - Euler's constant
!    D = SQRT(2.D0/PI)
!
  data tinyx / 1.0d-10/
  data a / 0.11593151565841244881d0/
  data d /0.797884560802865364d0/
!
!  Machine dependent parameters
!
  data sqxmin / 1.49d-154 /
  data xinf / 1.79d+308 /
  data xmin / 2.23d-308 /
  data xmax / 705.342d0 /
!
!  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA
!                                         + Euler's constant
!  Coefficients converted from hex to decimal and modified
!  by W. J. Cody, 2/26/82
!
!  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA)
!  T    - Approximation for SINH(Y)/Y
!
  data p/ 0.805629875690432845d00,    0.204045500205365151d02, &
          0.157705605106676174d03,    0.536671116469207504d03, &
          0.900382759291288778d03,    0.730923886650660393d03, &
          0.229299301509425145d03,    0.822467033424113231d00/
  data q/ 0.294601986247850434d02,    0.277577868510221208d03, &
          0.120670325591027438d04,    0.276291444159791519d04, &
          0.344374050506564618d04,    0.221063190113378647d04, &
          0.572267338359892221d03/
  data r/-0.48672575865218401848d+0,  0.13079485869097804016d+2, &
         -0.10196490580880537526d+3,  0.34765409106507813131d+3, &
          0.34958981245219347820d-3/
  data s/-0.25579105509976461286d+2,  0.21257260432226544008d+3, &
         -0.61069018684944109624d+3,  0.42269668805777760407d+3/
  data t/ 0.16125990452916363814d-9, 0.25051878502858255354d-7, &
          0.27557319615147964774d-5, 0.19841269840928373686d-3, &
          0.83333333333334751799d-2, 0.16666666666666666446d+0/
  data estm / 5.20583d1, 5.7607d0, 2.7782d0, 1.44303d1, 1.853004d2, &
            9.3715d0/
  data estf / 4.18341d1, 7.1075d0, 6.4306d0, 4.25110d1, 1.35633d0, &
            8.45096d1, 2.0d1/

  ex = x
  enu = alpha
  ncalc = min ( nb, 0 ) - 2
  !bk(1:nb) = 0

  if ( 0 < nb .and. &
    ( 0.0D+00 <= enu .and. enu < 1.0D+00 ) .and. &
    ( 1 <= ize .and. ize <= 2 ) .and. &
    ( ize /= 1 .or. ex <= xmax ) .and. &
    0.0D+00 < ex )  then

    k = 0
    if ( enu < sqxmin ) then
      enu = 0.0D+00
    end if

    if ( 0.5D+00 < enu ) then
      k = 1
      enu = enu - 1.0D+00
    end if

    twonu = enu + enu
    iend = nb + k - 1
    c = enu * enu
    d3 = -c

    if ( ex <= 1.0D+00 ) then
!
!  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA,
!                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA.
!
      d1 = 0.0D+00
      d2 = p(1)
      t1 = 1.0D+00
      t2 = q(1)

      do i = 2, 7, 2
        d1 = c * d1 + p(i)
        d2 = c * d2 + p(i+1)
        t1 = c * t1 + q(i)
        t2 = c * t2 + q(i+1)
      end do

      d1 = enu * d1
      t1 = enu * t1
      f1 = log ( ex )
      f0 = a + enu * ( p(8) &
         - enu * ( d1 + d2 ) / ( t1 + t2 ) ) - f1
      q0 = exp ( -enu * ( a - enu * &
         ( p(8) + enu * ( d1 - d2 ) / ( t1 - t2 ) ) - f1 ) )
      f1 = enu * f0
      p0 = exp ( f1 )
!
!  Calculation of F0.
!
      d1 = r(5)
      t1 = 1.0D+00
      do i = 1, 4
        d1 = c * d1 + r(i)
        t1 = c * t1 + s(i)
      end do

      if ( abs ( f1 ) <= 0.5D+00 ) then
        f1 = f1 * f1
        d2 = 0.0D+00
        do i = 1, 6
          d2 = f1 * d2 + t(i)
        end do
        d2 = f0 + f0 * f1 * d2
      else
        d2 = sinh ( f1 ) / enu
      end if

      f0 = d2 - enu * d1 / ( t1 * p0 )
!
!  X <= 1.0E-10.
!
!  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X).
!
      if ( ex <= tinyx ) then

        bk(1) = f0 + ex * f0

        if ( ize == 1 ) then
          bk(1) = bk(1) - ex * bk(1)
        end if

        ratio = p0 / f0
        c = ex * xinf
!
!  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
!  1/2 <= ALPHA.
!
        if ( k /= 0 ) then

          ncalc = -1

          if ( c / ratio <= bk(1) ) then
!	  k_arg = bk(nb)
            return
          end if

          bk(1) = ratio * bk(1) / ex
          twonu = twonu + 2.0D+00
          ratio = twonu

        end if

        ncalc = 1

        if ( nb == 1 ) then
	  k_arg = bk(nb)
          return
        end if
!
!  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1.
!
        ncalc = -1
        do i = 2, nb
          if ( c <= ratio ) then
!	  k_arg = bk(nb)
            return
          end if
          bk(i) = ratio / ex
          twonu = twonu + 2.0D+00
          ratio = twonu
        end do

        ncalc = 1
        j = ncalc + 1

        do i = j, nb
          if ( xinf / bk(i) <= bk(ncalc) ) then
	  k_arg = bk(nb)
            return
          end if
          bk(i) = bk(ncalc) * bk(i)
          ncalc = i
        end do
        k_arg = bk(nb)
        return
!
!  1.0E-10 < X <= 1.0.
!
      else

        c = 1.0D+00
        x2by4 = ex * ex / 4.0D+00
        p0 = 0.5D+00 * p0
        q0 = 0.5D+00 * q0
        d1 = - 1.0D+00
        d2 = 0.0D+00
        bk1 = 0.0D+00
        bk2 = 0.0D+00
        f1 = f0
        f2 = p0

  100       continue

        d1 = d1 + 2.0D+00
        d2 = d2 + 1.0D+00
        d3 = d1 + d3
        c = x2by4 * c / d2
        f0 = ( d2 * f0 + p0 + q0 ) / d3
        p0 = p0 / ( d2 - enu )
        q0 = q0 / ( d2 + enu )
        t1 = c * f0
        t2 = c * ( p0 - d2 * f0 )
        bk1 = bk1 + t1
        bk2 = bk2 + t2

        if ( epsilon ( t1 ) < abs ( t1 / ( f1 + bk1 ) ) .or. &
             epsilon ( t2 ) < abs ( t2 / ( f2 + bk2 ) ) )  then
          go to 100
        end if

        bk1 = f1 + bk1
        bk2 = 2.0D+00 * ( f2 + bk2 ) / ex

        if ( ize == 2 ) then
          d1 = exp ( ex )
          bk1 = bk1 * d1
          bk2 = bk2 * d1
        end if

        wminf = estf(1) * ex + estf(2)

      end if
!
!  1/EPS < X.
!
    else if ( 1.0D+00 < epsilon ( ex ) * ex ) then

      ncalc = nb
      bk1 = 1.0D+00 / ( d * sqrt ( ex ) )
      do i = 1, nb
        bk(i) = bk1
      end do
      k_arg = bk(nb)
      return

    else
!
!  1 < X.
!
      twox = ex + ex
      blpha = 0.0D+00
      ratio = 0.0D+00

      if ( ex <= 4.0D+00 ) then
!
!  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 <= X <= 4.0.
!
        d2 = aint ( estm(1) / ex + estm(2) )
        m = int ( d2 )
        d1 = d2 + d2
        d2 = d2 - 0.5D+00
        d2 = d2 * d2
        do i = 2, m
          d1 = d1 - 2.0D+00
          d2 = d2 - d1
          ratio = ( d3 + d2 ) / ( twox + d1 - ratio )
        end do
!
!  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
!  recurrence and K(ALPHA,X) from the Wronskian.
!
        d2 = aint ( estm(3) * ex + estm(4) )
        m = int ( d2 )
        c = abs ( enu )
        d3 = c + c
        d1 = d3 - 1.0D+00
        f1 = xmin
        f0 = ( 2.0D+00 * ( c + d2 ) / ex &
           + 0.5D+00 * ex / ( c + d2 + 1.0D+00 ) ) * xmin

        do i = 3, m
          d2 = d2 - 1.0D+00
          f2 = ( d3 + d2 + d2 ) * f0
          blpha = ( 1.0D+00 + d1 / d2 ) * ( f2 + blpha )
          f2 = f2 / ex + f1
          f1 = f0
          f0 = f2
        end do

        f1 = ( d3 + 2.0D+00 ) * f0 / ex + f1
        d1 = 0.0D+00
        t1 = 1.0D+00
        do i = 1, 7
          d1 = c * d1 + p(i)
          t1 = c * t1 + q(i)
        end do

        p0 = exp ( c * ( a + c * ( p(8) &
           - c * d1 / t1 ) - log ( ex ) ) ) / ex
        f2 = ( c + 0.5D+00 - ratio ) * f1 / ex
        bk1 = p0 + ( d3 * f0 - f2 + f0 + blpha ) &
          / ( f2 + f1 + f0 ) * p0

        if ( ize == 1 ) then
          bk1 = bk1 * exp ( - ex )
        end if

        wminf = estf(3) * ex + estf(4)

      else
!
!  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
!  recurrence, for 4 < X.
!
        dm = aint ( estm(5) / ex + estm(6) )
        m = int ( dm )
        d2 = dm - 0.5D+00
        d2 = d2 * d2
        d1 = dm + dm

        do i = 2, m
          dm = dm - 1.0D+00
          d1 = d1 - 2.0D+00
          d2 = d2 - d1
          ratio = ( d3 + d2 ) / ( twox + d1 - ratio )
          blpha = ( ratio + ratio * blpha ) / dm
        end do

        bk1 = 1.0D+00 / ( ( d + d * blpha ) * sqrt ( ex ) )

        if ( ize == 1 ) then
          bk1 = bk1 * exp ( - ex )
        end if

        wminf = estf(5) * ( ex - abs ( ex - estf(7) ) ) + estf(6)

      end if
!
!  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
!  K(ALPHA+1,X)/K(ALPHA,X).
!
      bk2 = bk1 + bk1 * ( enu + 0.5D+00 - ratio ) / ex

    end if
!
!  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1,
!  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1.
!
    ncalc = nb
    bk(1) = bk1

    if ( iend == 0 ) then
      k_arg = bk(nb)
      return
    end if

    j = 2 - k

    if ( 0 < j ) then
      bk(j) = bk2
    end if

    if ( iend == 1 ) then
      k_arg = bk(nb)
      return
    end if

    m = min ( int ( wminf - enu ), iend )

    do i = 2, m

      t1 = bk1
      bk1 = bk2
      twonu = twonu + 2.0D+00

      if ( ex < 1.0D+00 ) then

        if ( ( xinf / twonu ) * ex <= bk1 ) then
          exit
        end if

      else

        if ( xinf / twonu <= bk1 / ex ) then
          exit
        end if

      end if

      bk2 = twonu / ex * bk1 + t1
      itemp = i
      j = j + 1

      if ( 0 < j ) then
        bk(j) = bk2
      end if

    end do

    m = itemp

    if ( m == iend ) then
      k_arg = bk(nb)
      return
    end if

    ratio = bk2 / bk1
    mplus1 = m + 1
    ncalc = -1

    do i = mplus1, iend

      twonu = twonu + 2.0D+00
      ratio = twonu / ex + 1.0D+00 / ratio
      j = j + 1

      if ( 1 < j ) then
        bk(j) = ratio
      else
        if ( xinf / ratio <= bk2 ) then
	  k_arg = bk(nb)
          return
        end if
        bk2 = ratio * bk2
      end if

    end do

    ncalc = max ( mplus1 - k, 1 )

    if ( ncalc == 1 ) then
      bk(1) = bk2
    end if

    if ( nb == 1 ) then
      k_arg = bk(nb)
      return
    end if

    j = ncalc + 1

    do i = j, nb
      if ( xinf / bk(i) <= bk(ncalc) ) then
	k_arg = bk(nb)
        return
      end if
      bk(i) = bk(ncalc) * bk(i)
      ncalc = i
    end do

  end if
  k_arg = bk(nb)
  return
end subroutine rkbesl

subroutine gamma ( x, ga )
!*****************************************************************************80
!
!! GAMMA evaluates the Gamma function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into a user program that the copyright 
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!    X must not be 0, or any negative integer.
!
!    Output, real ( kind = 8 ) GA, the value of the Gamma function.
!
  implicit none

  real ( kind = 8 ), dimension ( 26 ) :: g = (/ &
    1.0D+00, &
    0.5772156649015329D+00, &
   -0.6558780715202538D+00, &
   -0.420026350340952D-01, &
    0.1665386113822915D+00, &
   -0.421977345555443D-01, &
   -0.96219715278770D-02, &
    0.72189432466630D-02, &
   -0.11651675918591D-02, &
   -0.2152416741149D-03, &
    0.1280502823882D-03, & 
   -0.201348547807D-04, &
   -0.12504934821D-05, &
    0.11330272320D-05, &
   -0.2056338417D-06, & 
    0.61160950D-08, &
    0.50020075D-08, &
   -0.11812746D-08, &
    0.1043427D-09, & 
    0.77823D-11, &
   -0.36968D-11, &
    0.51D-12, &
   -0.206D-13, &
   -0.54D-14, &
    0.14D-14, &
    0.1D-15 /)
  real ( kind = 8 ) ga
  real ( kind = 8 ) gr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) z

  if ( x == aint ( x ) ) then

    if ( 0.0D+00 < x ) then
      ga = 1.0D+00
      m1 = int ( x ) - 1
      do k = 2, m1
        ga = ga * k
      end do
    else
      ga = 1.0D+300
    end if

  else

    if ( 1.0D+00 < abs ( x ) ) then
      z = abs ( x )
      m = int ( z )
      r = 1.0D+00
      do k = 1, m
        r = r * ( z - real ( k, kind = 8 ) )
      end do
      z = z - real ( m, kind = 8 )
    else
      z = x
    end if

    gr = g(26)
    do k = 25, 1, -1
      gr = gr * z + g(k)
    end do

    ga = 1.0D+00 / ( gr * z )

    if ( 1.0D+00 < abs ( x ) ) then
      ga = ga * r
      if ( x < 0.0D+00 ) then
        ga = - pi / ( x* ga * sin ( pi * x ) )
      end if
    end if

  end if

  return

end subroutine gamma
 
subroutine kron(K, A, B)
! Calculate the Kronecker product of two matrices
    IMPLICIT NONE
    REAL, INTENT(IN)  :: A(:,:), B(:,:)
    REAL, INTENT(INOUT) :: K(:,:)
    INTEGER :: I, J, MA, NA, MB, NB
    MA = UBOUND(A, 1)
    NA = UBOUND(A, 2)
    MB = UBOUND(B, 1)
    NB = UBOUND(B, 2)
    IF (SIZE(K,1) /= MA*MB .OR. SIZE(K,2) /= NA*NB) THEN
        WRITE(*,*) 'K has invalid size'
        CALL ABORT
    END IF

    FORALL(I=1:MA, J=1:NA)
        K(MB*(I-1)+1:MB*I,NB*(J-1)+1:NB*J) = A(I,J)*B
    END FORALL

END subroutine kron

End module mod_MCMC

