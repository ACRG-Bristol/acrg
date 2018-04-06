/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/***************************************************************************
 *   Enrico Fucile  - 19.06.2007                                           *
 *                                                                         *
 ***************************************************************************/

unsigned long grib_decode_unsigned_long(const unsigned char* p, long *bitp, long nbits)
{
  long countOfLeftmostBits=0,leftmostBits=0;
  long startBit= *bitp ;
  long remainingBits = nbits;
  long *pp=(long*)p;
  unsigned long val=0;

  if (startBit >= max_nbits) {
    pp += startBit/max_nbits;
    startBit %= max_nbits;
  }

  countOfLeftmostBits = startBit + remainingBits;
  if (countOfLeftmostBits > max_nbits) {
    countOfLeftmostBits = max_nbits - startBit;
    remainingBits -= countOfLeftmostBits;
    leftmostBits=(VALUE(*pp,startBit,countOfLeftmostBits)) << remainingBits;
    startBit = 0;
    pp++;
  } else
    leftmostBits = 0;

  val=leftmostBits+(VALUE(*pp,startBit,remainingBits));

  *bitp += nbits;

  return val;
}

int grib_encode_unsigned_long(unsigned char* p, unsigned long val ,long *bitp, long nbits)
{
  long* destination = (long*)p;
  long countOfLeftmostBits=0,nextWord=0,startBit=0,remainingBits=0,rightmostBits=0;

  startBit = *bitp;
  remainingBits = nbits;

  if (startBit >= max_nbits) {
    nextWord = startBit / max_nbits;
    startBit %= max_nbits;
  } else
    nextWord = 0;

  countOfLeftmostBits = startBit + remainingBits;
  if (countOfLeftmostBits > max_nbits) {
    countOfLeftmostBits = max_nbits - startBit;
    startBit = max_nbits - remainingBits;
    remainingBits -= countOfLeftmostBits;
    destination[nextWord] =
      ((destination[nextWord] >> countOfLeftmostBits) << countOfLeftmostBits)
      + (VALUE(val,startBit,countOfLeftmostBits));
    startBit = 0;
    nextWord++;
  }

  rightmostBits = VALUE(val,max_nbits-remainingBits,remainingBits);
  destination[nextWord] =
    (destination[nextWord] & ~MASKVALUE(startBit,remainingBits))
    + (rightmostBits << max_nbits-(remainingBits+startBit));

  *bitp+=nbits;
  return GRIB_SUCCESS;
}

int grib_encode_unsigned_longb(unsigned char* p, unsigned long val ,long *bitp, long nbits)
{
  return grib_encode_unsigned_long(p,val ,bitp, nbits);
}

#if VECTOR

#include "grib_bits_fast_big_endian_vector.c"

#elif OMP

#include "grib_bits_fast_big_endian_omp.c"

#else

#include "grib_bits_fast_big_endian_simple.c"

#endif
