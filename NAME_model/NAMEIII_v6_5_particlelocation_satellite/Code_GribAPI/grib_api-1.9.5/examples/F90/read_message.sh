	#!/bin/sh

	. ./include.sh

	in=../../data/index.grib
	out=out.grib

    rm -f $out | true

	${examples_dir}read_message > /dev/null

	${tools_dir}grib_compare $in $out

	rm -f $out

