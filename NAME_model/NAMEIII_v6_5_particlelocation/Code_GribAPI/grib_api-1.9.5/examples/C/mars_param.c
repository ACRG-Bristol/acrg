/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: set
 *
 * Description: how to set key values.
 *
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>

#include "grib_api.h"

int main(int argc, char** argv) {
	grib_string_list* list;
	grib_string_list* list2;

	list=grib_util_get_param_id("11.2");
	printf("mars.param=11.2 -> paramId= ");
	while (list) {
		printf("%s ",list->value);
	list2=grib_util_get_param_id("130.128");
	printf("mars.param=11.2 -> paramId= ");
	while (list2) {
		printf("%s ",list2->value);
		list2=list2->next;
	}
		list=list->next;
	}
	printf("\n");
	printf("\n");

	printf("paramId=130 -> mars.param= ");
	list=grib_util_get_mars_param("130");
	while (list) {
        printf("%s ",list->value);
        list=list->next;
    }
	printf("\n");

  return 0;
}
