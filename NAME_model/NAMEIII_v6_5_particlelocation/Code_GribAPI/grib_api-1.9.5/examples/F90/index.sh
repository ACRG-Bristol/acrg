	#!/bin/sh

	. ./include.sh

    # if an index file does not exist then
    # create index and save to disk
	${examples_dir}index > index_f90.out
	diff index_f90.out ${data_dir}/index_f90.ok

    # if an index file exists then load the index from it
	${examples_dir}index > index_f90.out
	diff index_f90.out ${data_dir}/index_f90.ok
