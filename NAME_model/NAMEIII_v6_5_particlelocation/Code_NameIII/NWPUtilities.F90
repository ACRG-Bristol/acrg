Module NWPUtilitiesModule

Use ServiceModule

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL   SUBROUTINE XPND_ALL -----------------
!LL
!LL   PURPOSE:   UNPACK FROM WGDOS FORMAT
!LL
!LL   (Note that for optimal performance the routines INSTIN and
!LL    EXTRIN are inline expanded (enable option 8 on fpp)
!LL
!LL
!LL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
!LL   STANDARD B, VERSION 2, DATED 18/01/90
!LL
!LL  Logical component number: S72
!LL
!LL   SYSTEM TASK: P7
!LL
!LLEND-------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owner's HTML page
!  This file belongs in section: STASH
Contains

   SUBROUTINE XPND_ALL(FIELD,ICOMP64,N,IX,IY,ISC,RMDI,ICODE,CMESSAGE)

!     Better vectorized version of XPND expanding the whole field at
!     once


      IMPLICIT NONE

!     Subroutine arguments

      INTEGER, INTENT(IN)  :: n, ix, iy, isc
      INTEGER(I64), INTENT(IN)  :: icomp64(n)
      REAL,    INTENT(IN)  :: rmdi
      REAL,    INTENT(OUT) :: field(ix,iy)
      INTEGER              :: icode
      CHARACTER(MaxCharLength)    :: cmessage

!     Local variables

      INTEGER(I64) :: i, j, nshft, num, iword, ioff, imask, ival, mant, iexp
      INTEGER(I64) :: i1, i2, nbits_bmap
      INTEGER(I64), DIMENSION(3*ix) :: itmp
      INTEGER(I64), DIMENSION(ix)   :: idx, imap
      INTEGER(I64), DIMENSION(iy)   :: istart, nop, ibase, nbits
      INTEGER(I64), DIMENSION(iy*(2*ix+2)+4) :: icomp

      REAL                 :: aprec
      REAL, DIMENSION (iy) :: base
      REAL, DIMENSION (ix) :: tmp

      LOGICAL, DIMENSION(iy) :: obtzer, obtmin, obtmis, obtmap

      INTEGER(I64), PARAMETER :: MASK16 = Z'FFFF'
      INTEGER(I64), PARAMETER :: MASK32 = Z'FFFFFFFF'
      INTEGER(I64), PARAMETER :: MASK_MANT_IBM = Z'00FFFFFF'
      INTEGER(I64), PARAMETER :: MASK_EXPT_IBM = Z'7F000000'
      INTEGER(I64), PARAMETER :: MASK_SIGN_IBM = Z'80000000'

      INTEGER(I64), PARAMETER :: ONE = 1
      INTEGER(I64), SAVE :: MASK_BITS(0:63), first

      DATA first /1/

      IF (first/=0) THEN
        DO i=0,63
          MASK_BITS(i) = ISHFT(ONE,63-i)
        END DO
        first = 0
      END IF

! Scale factor

      aprec = 2.**isc

! All lengths and alignments in WGDOS packing are for 32-bit words,
! so life gets much easier when we treat the packed data as 32-bit
! words.
! We split therefore the 64-bit compressed data into two 32 bit words

      num = ISHFT(icomp64(1),-32) ! Number of 32 bit words

      IF (num > SIZE(icomp)-2) THEN
        ICODE = 2
        CMESSAGE='COEX: Compressed data has too many elements'
        RETURN
      END IF

      DO i=1,(num+1)/2
        icomp(2*i-1) = IAND(ISHFT(icomp64(i),-32),MASK32)
        icomp(2*i)   = IAND(icomp64(i),MASK32)

      ENDDO
      ! The following word MUST be 0, it is used during decomposition!
      icomp(num+1) = 0
      icomp(num+2) = 0

! Get start word and length of every row

      istart(1) = 6
      nop(1) = IAND(icomp(5),MASK16)

      DO j=2,iy
        istart(j) = istart(j-1) + nop(j-1) + 2
        nop(j) = IAND(icomp(istart(j)-1),mask16)
        IF (istart(j)+nop(j)-1>num) THEN
          ICODE = 2
          CMESSAGE='COEX: Compressed data inconsistent'
          RETURN
        END IF
      END DO

! Get base (as a 32-bit IBM floating point number) and number of bits
! for every row and convert IBM floats to native floats
! The routine IBM2IEEE does a rather bad job, so we code it explicitly

      DO j=1,iy
        ibase(j) = icomp(istart(j)-2)
        nbits(j) = IAND(ISHFT(icomp(istart(j)-1),-16),mask16)

        mant = IAND(ibase(j),MASK_MANT_IBM)
        iexp = ISHFT(IAND(ibase(j),MASK_EXPT_IBM),-24)-64-6
        base(j) = 16.0**iexp*mant
        IF (IAND(ibase(j),MASK_SIGN_IBM) /= 0) base(j) = -base(j)

! Check if bitmaps are used

        obtzer(j) = IAND(nbits(j),128) /= 0
        obtmin(j) = IAND(nbits(j),64)  /= 0
        obtmis(j) = IAND(nbits(j),32)  /= 0
        obtmap(j) = obtzer(j) .OR. obtmin(j) .OR. obtmis(j)
        nbits(j)  = IAND(nbits(j),31)

! Decode data row by row


        ! Care about bitmaps

        imap(:) = 1 ! Data present indicator

        nbits_bmap = 0
        IF (obtmis(j)) nbits_bmap = nbits_bmap + ix
        IF (obtmin(j)) nbits_bmap = nbits_bmap + ix
        IF (obtzer(j)) nbits_bmap = nbits_bmap + ix

        IF (nbits_bmap > 0) THEN
          iword = istart(j)
          DO i1=1,nbits_bmap,64
            ival  = IOR(ISHFT(icomp(iword),32),icomp(iword+1))
            iword = iword+2
            DO i2=0,MIN(nbits_bmap-i1,63)
              itmp(i1+i2) = MERGE(1,0,IAND(ival,MASK_BITS(i2))/=0)
            END DO
          END DO
          istart(j) = istart(j) + (nbits_bmap+31)/32
        END IF

        nbits_bmap = 0

        ! Extract missing data bitmap

        IF (obtmis(j)) THEN
          WHERE(itmp(nbits_bmap+1:nbits_bmap+ix)/=0)
            field(:,j) = rmdi
            imap (:) = 0
          END WHERE
          nbits_bmap = nbits_bmap + ix
        END IF

        ! Extract minimum value bitmap

        IF(obtmin(j)) THEN
          WHERE(itmp(nbits_bmap+1:nbits_bmap+ix)/=0)
            field(:,j) = base(j)
            imap (:) = 0
          END WHERE
          nbits_bmap = nbits_bmap + ix
        END IF

        ! Extract zero value bitmap

        IF(obtzer(j)) THEN
          WHERE(itmp(nbits_bmap+1:nbits_bmap+ix)==0)
            field(:,j) = 0.
            imap (:) = 0
          END WHERE
          nbits_bmap = nbits_bmap + ix
        END IF

        IF(nbits(j)==0) THEN

          ! All points in row have same value

          IF(obtmap(j)) THEN
            WHERE(imap(:)/=0) field(:,j) = base(j)
          ELSE
            field(:,j) = base(j)
          END IF

        ELSE

          ! Get number [and index] of values to decode

          IF(obtmap(j)) THEN
            num = 0
            DO i=1,ix
              IF(imap(i) /= 0) THEN
                num = num+1
                idx(num) = i
              END IF
            END DO
          ELSE
            num = ix
          END IF

          imask = ISHFT(ONE,nbits(j))-1

          ! Decode data

          DO i=1,num

            ! Bit offset to value:
            ioff  = (i-1)*nbits(j)

            ! Number of word in icomp which contains first bit:
            iword = ISHFT(ioff,-5)+istart(j)

            ! We load this word and the following into ival,
            ! this way we don't have to care if a word boundary
            ! is crossed. This requires that ival is a 64 bit word!
            ival  = IOR(ISHFT(icomp(iword),32),icomp(iword+1))

            ! Number of bits we have to shift to the right:
            nshft = 64 - IAND(ioff,31) - nbits(j)

            ival = ISHFT(ival,-nshft)

            ! Mask ival and calculate decoded value:
            ival = IAND(ival,imask)
            tmp(i) = ival*aprec + base(j)
          END DO

          ! Write decoded values to field

          IF (obtmap(j)) THEN
            field(idx(1:num),j) = tmp(1:num)
          ELSE
            field(:,j) = tmp(:)
          END IF

        END IF

      END DO

      END SUBROUTINE XPND_ALL

End Module NWPUtilitiesModule
