	#!/bin/sh

	. ./include.sh

	${examples_dir}multi_write > /dev/null

    ${tools_dir}/grib_cmp ${data_dir}/multi_created.grib2 multi_created.grib2

#	rm -f multi_created.grib2

