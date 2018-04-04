/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

#include <stdio.h>

#include "grib_api_internal.h"

int main(int argc, char *argv[])
{
	int i = 0;
	int fail = 0;

	grib_context* c = grib_context_get_default();
grib_action *a;

	for(i = 1; i < argc; i++){
		if(!(a = grib_parse_file(c,argv[i])))
		{
			fail++;
			exit(1);
		}
		else
		{
grib_xref_action_branch(stdout,a,argv[i]);
	
		}
	}
	return fail;
}
