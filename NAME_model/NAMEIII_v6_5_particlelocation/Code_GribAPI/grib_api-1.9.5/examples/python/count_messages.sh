#!/bin/sh

. ./include.sh

TEMP1=temp1
TEMP2=temp2

$PYTHON count_messages.py 2> $TEMP1 > $TEMP1
./count_messages ../../data/tigge_pf_ecmwf.grib2 2> $TEMP2 > $TEMP2

diff $TEMP1 $TEMP2
rm $TEMP1 $TEMP2 || true
