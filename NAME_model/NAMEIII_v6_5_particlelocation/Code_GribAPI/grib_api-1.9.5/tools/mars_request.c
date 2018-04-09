/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

#include <stdio.h>
#include <string.h>

#include "grib_api.h"

void usage(const char *prog)
{
	fprintf(stderr,"%s: file\n",prog);
	exit(1);
}


int main(int argc, const char *argv[])
{
	FILE *in;
	int e;
	grib_handle *h;

	if(argc != 2) usage(argv[0]);

	in = fopen(argv[argc-1],"r");
	if(!in) {
		perror(argv[argc-1]);
		exit(1);
	}

	while( (h = grib_handle_new_from_file(NULL,in,&e)) != NULL )
	{
		grib_keys_iterator* i  = grib_keys_iterator_new(h,GRIB_KEYS_ITERATOR_ALL_KEYS,"mars");

		while(grib_keys_iterator_next(i))
		{
			char value[1024];
			size_t len = sizeof(value);
			const char* name;

			name=grib_keys_iterator_get_name(i);
			GRIB_CHECK(grib_get_string(h,name,value,&len),0);

			printf("%s = %s\n",name,value);

		}

		grib_keys_iterator_delete(i);
		grib_handle_delete(h);
	}

	GRIB_CHECK(e,argv[argc-1]);

	return 0;
}
