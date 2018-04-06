/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */


#include "grib_api_internal.h"

/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = box
   IMPLEMENTS = destroy
   IMPLEMENTS = get_points
   IMPLEMENTS = get_values
   IMPLEMENTS = init
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "box.class" and rerun ./make_class.pl

*/


static void init_class              (grib_box_class*);

static int init               (grib_box* box,grib_handle* h,grib_arguments* args);
static grib_points* get_points(grib_box* box, double north, double west, double south,double east,int *err);
static int destroy            (grib_box* box);

typedef struct grib_box_gen{
  grib_box box;
/* Members defined in gen */
} grib_box_gen;


static grib_box_class _grib_box_class_gen = {
    0,                         /* super                     */
    "gen",                         /* name                      */
    sizeof(grib_box_gen),      /* size of instance          */
    0,                              /* inited */
    &init_class,                    /* init_class */
    &init,                          /* constructor               */
    &destroy,                       /* destructor                */
    &get_points,                    /* get points           */
};

grib_box_class* grib_box_class_gen = &_grib_box_class_gen;


static void init_class(grib_box_class* c)
{
}
/* END_CLASS_IMP */

static int init(grib_box* box,grib_handle* h,grib_arguments* args)
{
  return 0;
}

static int destroy(grib_box* box) {
  grib_points_delete(box->points);
  box->points=0;
  return GRIB_SUCCESS;
}

static grib_points* get_points(grib_box* box, double north, double west, double south,double east, int *err) {
	*err=GRIB_NOT_IMPLEMENTED;
	return NULL;
}



