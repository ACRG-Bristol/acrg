/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/**
* Author: Enrico Fucile
* date: 14/06/2010
*
*/

#include "grib_api_internal.h"


grib_points* grib_box_get_points(grib_box *box,double north, double west,double south,double east, int *err)
{
  grib_box_class *c = box->cclass;
  while(c)
  {
    grib_box_class *s = c->super ? *(c->super) : NULL;
    if(c->get_points) {
		return c->get_points(box,north,west,south,east,err);
	}
    c = s;
  }
  Assert(0);
  return 0;
}


/* For this one, ALL init are called */

static int init_box(grib_box_class* c,grib_box* i, grib_handle *h, grib_arguments* args)
{

  if(c) {
    int ret = GRIB_SUCCESS;
    grib_box_class *s = c->super ? *(c->super) : NULL;
    if(!c->inited)
    {
      if(c->init_class) c->init_class(c);
      c->inited = 1;
    }
    if(s) ret = init_box(s,i,h,args);

    if(ret != GRIB_SUCCESS) return ret;

    if(c->init) return c->init(i,h, args);
  }
  return GRIB_INTERNAL_ERROR;
}

int grib_box_init(grib_box* box, grib_handle *h, grib_arguments* args)
{
  return init_box(box->cclass,box,h,args);
}

/* For this one, ALL destroy are called */

int grib_box_delete(grib_box *box)
{
  grib_box_class *c = box->cclass;
  while(c)
  {
    grib_box_class *s = c->super ? *(c->super) : NULL;
    if(c->destroy) c->destroy(box);
    c = s;
  }
  return 0;
}

grib_points* grib_points_new(grib_context* c,size_t size) {
	grib_points* points=grib_context_malloc_clear(c,sizeof(grib_points));

	points->latitudes=grib_context_malloc_clear(c,sizeof(double)*size);
    points->longitudes=grib_context_malloc_clear(c,sizeof(double)*size);
	points->indexes=grib_context_malloc_clear(c,sizeof(double)*size);
	points->group_start=grib_context_malloc_clear(c,sizeof(double)*size);
	points->group_len=grib_context_malloc_clear(c,sizeof(double)*size);
	points->size=size;
	points->context=c;

	return points;
}

void grib_points_delete(grib_points *points) {

	grib_context* c;
	if (!points) return;
	
	c=points->context;
	grib_context_free(c,points->latitudes);
	grib_context_free(c,points->longitudes);
	grib_context_free(c,points->indexes);
	grib_context_free(c,points->group_start);
	grib_context_free(c,points->group_len);
	grib_context_free(c,points);

}
