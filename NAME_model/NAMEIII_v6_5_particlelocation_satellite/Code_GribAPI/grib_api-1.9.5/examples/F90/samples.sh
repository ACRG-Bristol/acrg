	#!/bin/sh

	. ./include.sh

    rm -f ${data_dir}/regular_latlon_surface.grib1.tmpl | true
    cp ${data_dir}/regular_latlon_surface.grib1 ${data_dir}/regular_latlon_surface.grib1.tmpl
    GRIB_SAMPLES_PATH=${data_dir}
	export GRIB_SAMPLES_PATH
	${examples_dir}samples > /dev/null
