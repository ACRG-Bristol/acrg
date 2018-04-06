#!/bin/sh

. ./include.sh

REDIRECT=/dev/null
$PYTHON get.py 2> $REDIRECT > $REDIRECT
