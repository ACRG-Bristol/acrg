/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

#include "grib_api_internal.h"
/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = accessor
   SUPER      = grib_accessor_class_bytes

   IMPLEMENTS = next_offset
   IMPLEMENTS = unpack_double;unpack_double_element
   IMPLEMENTS = unpack_string
   IMPLEMENTS = init;dump;update_size
   MEMBERS=const char* tableReference
   MEMBERS=const char* missing_value
   MEMBERS=const char* offsetbsec
   MEMBERS=const char* sLength
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int unpack_double(grib_accessor*, double* val,size_t *len);
static int unpack_string (grib_accessor*, char*, size_t *len);
static long next_offset(grib_accessor*);
static void dump(grib_accessor*, grib_dumper*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);
static void update_size(grib_accessor*,size_t);
static int unpack_double_element(grib_accessor*,size_t i, double* val);

typedef struct grib_accessor_bitmap {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in bytes */
/* Members defined in bitmap */
	const char* tableReference;
	const char* missing_value;
	const char* offsetbsec;
	const char* sLength;
} grib_accessor_bitmap;

extern grib_accessor_class* grib_accessor_class_bytes;

static grib_accessor_class _grib_accessor_class_bitmap = {
    &grib_accessor_class_bytes,                      /* super                     */
    "bitmap",                      /* name                      */
    sizeof(grib_accessor_bitmap),  /* size                      */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                       /* init                      */
    0,                  /* post_init                      */
    0,                    /* free mem                       */
    &dump,                       /* describes himself         */
    &next_offset,                /* get length of section     */
    0,                /* get number of values      */
    0,                 /* get number of bytes      */
    0,                /* get offset to bytes           */
    0,            /* get native type               */
    0,                /* get sub_section                */
    0,               /* grib_pack procedures long      */
    0,               /* grib_pack procedures long      */
    0,                  /* grib_pack procedures long      */
    0,                /* grib_unpack procedures long    */
    0,                /* grib_pack procedures double    */
    &unpack_double,              /* grib_unpack procedures double  */
    0,                /* grib_pack procedures string    */
    &unpack_string,              /* grib_unpack procedures string  */
    0,                 /* grib_pack procedures bytes     */
    0,               /* grib_unpack procedures bytes   */
    0,            /* pack_expression */
    0,              /* notify_change   */
    &update_size,                /* update_size   */
    0,            /* preferred_size   */
    0,                    /* resize   */
    0,      /* nearest_smaller_value */
    0,                       /* next accessor    */
    0,                    /* compare vs. another accessor   */
    &unpack_double_element,     /* unpack only ith value          */
    0,     /* unpack a subarray         */
    0,             		/* clear          */
};


grib_accessor_class* grib_accessor_class_bitmap = &_grib_accessor_class_bitmap;


static void init_class(grib_accessor_class* c)
{
	c->value_count	=	(*(c->super))->value_count;
	c->byte_count	=	(*(c->super))->byte_count;
	c->byte_offset	=	(*(c->super))->byte_offset;
	c->get_native_type	=	(*(c->super))->get_native_type;
	c->sub_section	=	(*(c->super))->sub_section;
	c->pack_missing	=	(*(c->super))->pack_missing;
	c->is_missing	=	(*(c->super))->is_missing;
	c->pack_long	=	(*(c->super))->pack_long;
	c->unpack_long	=	(*(c->super))->unpack_long;
	c->pack_double	=	(*(c->super))->pack_double;
	c->pack_string	=	(*(c->super))->pack_string;
	c->pack_bytes	=	(*(c->super))->pack_bytes;
	c->unpack_bytes	=	(*(c->super))->unpack_bytes;
	c->pack_expression	=	(*(c->super))->pack_expression;
	c->notify_change	=	(*(c->super))->notify_change;
	c->preferred_size	=	(*(c->super))->preferred_size;
	c->resize	=	(*(c->super))->resize;
	c->nearest_smaller_value	=	(*(c->super))->nearest_smaller_value;
	c->next	=	(*(c->super))->next;
	c->compare	=	(*(c->super))->compare;
	c->unpack_double_subarray	=	(*(c->super))->unpack_double_subarray;
	c->clear	=	(*(c->super))->clear;
}

/* END_CLASS_IMP */


static void compute_size(grib_accessor* a)
{
  long slen = 0;
  long off = 0;

  grib_accessor_bitmap* self = (grib_accessor_bitmap*)a;
  grib_get_long_internal(a->parent->h, self->offsetbsec,&off);
  grib_get_long_internal(a->parent->h, self->sLength, &slen);

  if(slen == 0)
  {
    grib_accessor* seclen;
    size_t size;
    /* Assume reparsing */
    Assert(a->parent->h->loader != 0);
    if (a->parent->h->loader != 0) {
      seclen = grib_find_accessor(a->parent->h, self->sLength);
      Assert(seclen);
      grib_get_block_length(seclen->parent,&size);
      slen = size;
    }
  }

#if 0
  printf("compute_size off=%ld slen=%ld a->offset=%ld\n",
    (long)off,(long)slen,(long)a->offset);
#endif

  a->length = off+(slen-a->offset);

  if(a->length < 0)
  {
    /* Assume reparsing */
    /*Assert(a->parent->h->loader != 0);*/
    a->length = 0;
  }


  Assert(a->length>=0);
}

static void init(grib_accessor* a, const long len , grib_arguments* arg )
{

  grib_accessor_bitmap* self = (grib_accessor_bitmap*)a;
  int n = 0;

  self->tableReference = grib_arguments_get_name(a->parent->h,arg,n++);
  self->missing_value  = grib_arguments_get_name(a->parent->h,arg,n++);
  self->offsetbsec     = grib_arguments_get_name(a->parent->h,arg,n++);
  self->sLength        = grib_arguments_get_name(a->parent->h,arg,n++);

  compute_size(a);

}

static long next_offset(grib_accessor* a){
  return grib_byte_offset(a) + grib_byte_count(a);
}

static void dump(grib_accessor* a, grib_dumper* dumper)
{
  long len = grib_value_count(a);
  char label[1024];
  sprintf(label,"Bitmap of %ld values",len);
  grib_dump_bytes(dumper, a,label);
}


static int unpack_double   (grib_accessor* a, double* val, size_t *len)
{
  long pos = a->offset*8;
  size_t tlen;
  long i;

  tlen = grib_value_count(a);

  if(*len < tlen)
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, "Wrong size for %s it contains %d values ", a->name , tlen );
    *len = 0;
    return GRIB_ARRAY_TOO_SMALL;
  }

  for(i=0;i<tlen;i++)
  {
    val[i] = (double)grib_decode_unsigned_long(a->parent->h->buffer->data, &pos,1);
  }
  *len = tlen;
  return GRIB_SUCCESS;

}

static int unpack_double_element   (grib_accessor* a, size_t idx, double* val)
{
  long pos = a->offset*8;

  pos+=idx;
  *val = (double)grib_decode_unsigned_long(a->parent->h->buffer->data, &pos,1);

  return GRIB_SUCCESS;

}


static void update_size(grib_accessor* a,size_t s)
{
  a->length = s;
}

static int unpack_string(grib_accessor* a, char* val, size_t *len)
{

  int i = 0;

  if(len[0] < (a->length))
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, "unpack_string: Wrong size (%d) for %s it contains %d values ", len[0], a->name , a->length );
    len[0] = 0;
    return GRIB_ARRAY_TOO_SMALL;
  }

  for ( i = 0; i < a->length; i++)
    val[i] = a->parent->h->buffer->data[a->offset+i];

  len[0] = a->length;

  return GRIB_SUCCESS;
}
