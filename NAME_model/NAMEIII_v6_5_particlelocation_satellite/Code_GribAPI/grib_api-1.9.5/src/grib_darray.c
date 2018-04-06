/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/***************************************************************************
 *
 *   Enrico Fucile
 *
 ***************************************************************************/

#include "grib_api_internal.h"

grib_darray* grib_darray_new(grib_context* c,size_t size,size_t incsize) {
  grib_darray* v=NULL;
  if (!c) c=grib_context_get_default();
  v=(grib_darray*)grib_context_malloc(c,sizeof(grib_darray));
  if (!v) {
    grib_context_log(c,GRIB_LOG_ERROR,
          "grib_darray_new unable to allocate %d bytes\n",sizeof(grib_darray));
    return NULL;
  }
  v->size=size;
  v->n=0;
  v->incsize=incsize;
  v->v=(double*)grib_context_malloc(c,sizeof(double)*size);
  if (!v->v) {
    grib_context_log(c,GRIB_LOG_ERROR,
          "grib_darray_new unable to allocate %d bytes\n",sizeof(double)*size);
    return NULL;
  }
  return v;
}

grib_darray* grib_darray_resize(grib_context* c,grib_darray* v) {
  int newsize=v->incsize+v->size;

  if (!c) c=grib_context_get_default();

  v->v=grib_context_realloc(c,v->v,newsize*sizeof(double));
  v->size=newsize;
  if (!v->v) {
    grib_context_log(c,GRIB_LOG_ERROR,
          "grib_darray_resize unable to allocate %d bytes\n",sizeof(double)*newsize);
    return NULL;
  }
  return v;
}

grib_darray* grib_darray_push(grib_context* c,grib_darray* v,double val) {
  size_t start_size=100;
  size_t start_incsize=100;
  if (!v) v=grib_darray_new(c,start_size,start_incsize);

  if (v->n >= v->size) v=grib_darray_resize(c,v);
  v->v[v->n]=val;
  v->n++;
  return v;
}

void grib_darray_delete(grib_context* c,grib_darray* v) {
  if (!v) return;
  if (!c) grib_context_get_default();
  if (v->v) grib_context_free(c,v->v);
  grib_context_free(c,v);
}

