	#!/bin/sh

	. ./include.sh

	${examples_dir}multi > multi.out
	diff multi.out ${data_dir}/multi.ok
