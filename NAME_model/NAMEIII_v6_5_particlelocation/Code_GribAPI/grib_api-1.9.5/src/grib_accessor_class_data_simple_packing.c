/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/********************************
 *   Enrico Fucile
 *******************************/

#include "grib_api_internal.h"

/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = accessor
   SUPER      = grib_accessor_class_values
   IMPLEMENTS = init
   IMPLEMENTS = unpack_double
   IMPLEMENTS = unpack_double_element
   IMPLEMENTS = unpack_double_subarray
   IMPLEMENTS = pack_double
   IMPLEMENTS = value_count
   MEMBERS=const char*  units_factor
   MEMBERS=const char*  units_bias  
   MEMBERS=const char*  changing_precision  
   MEMBERS=const char*  number_of_values
   MEMBERS=const char*  bits_per_value
   MEMBERS=const char*  reference_value
   MEMBERS=const char*  binary_scale_factor
   MEMBERS=const char*  decimal_scale_factor
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int pack_double(grib_accessor*, const double* val,size_t *len);
static int unpack_double(grib_accessor*, double* val,size_t *len);
static long value_count(grib_accessor*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);
static int unpack_double_element(grib_accessor*,size_t i, double* val);
static int unpack_double_subarray(grib_accessor*, double* val,size_t start,size_t len);

typedef struct grib_accessor_data_simple_packing {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in values */
	int  carg;
	const char* seclen;
	const char* offsetdata;
	const char* offsetsection;
	int dirty;
/* Members defined in data_simple_packing */
	const char*  units_factor;
	const char*  units_bias;
	const char*  changing_precision;
	const char*  number_of_values;
	const char*  bits_per_value;
	const char*  reference_value;
	const char*  binary_scale_factor;
	const char*  decimal_scale_factor;
} grib_accessor_data_simple_packing;

extern grib_accessor_class* grib_accessor_class_values;

static grib_accessor_class _grib_accessor_class_data_simple_packing = {
    &grib_accessor_class_values,                      /* super                     */
    "data_simple_packing",                      /* name                      */
    sizeof(grib_accessor_data_simple_packing),  /* size                      */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                       /* init                      */
    0,                  /* post_init                      */
    0,                    /* free mem                       */
    0,                       /* describes himself         */
    0,                /* get length of section     */
    &value_count,                /* get number of values      */
    0,                 /* get number of bytes      */
    0,                /* get offset to bytes           */
    0,            /* get native type               */
    0,                /* get sub_section                */
    0,               /* grib_pack procedures long      */
    0,               /* grib_pack procedures long      */
    0,                  /* grib_pack procedures long      */
    0,                /* grib_unpack procedures long    */
    &pack_double,                /* grib_pack procedures double    */
    &unpack_double,              /* grib_unpack procedures double  */
    0,                /* grib_pack procedures string    */
    0,              /* grib_unpack procedures string  */
    0,                 /* grib_pack procedures bytes     */
    0,               /* grib_unpack procedures bytes   */
    0,            /* pack_expression */
    0,              /* notify_change   */
    0,                /* update_size   */
    0,            /* preferred_size   */
    0,                    /* resize   */
    0,      /* nearest_smaller_value */
    0,                       /* next accessor    */
    0,                    /* compare vs. another accessor   */
    &unpack_double_element,     /* unpack only ith value          */
    &unpack_double_subarray,     /* unpack a subarray         */
    0,             		/* clear          */
};


grib_accessor_class* grib_accessor_class_data_simple_packing = &_grib_accessor_class_data_simple_packing;


static void init_class(grib_accessor_class* c)
{
	c->dump	=	(*(c->super))->dump;
	c->next_offset	=	(*(c->super))->next_offset;
	c->byte_count	=	(*(c->super))->byte_count;
	c->byte_offset	=	(*(c->super))->byte_offset;
	c->get_native_type	=	(*(c->super))->get_native_type;
	c->sub_section	=	(*(c->super))->sub_section;
	c->pack_missing	=	(*(c->super))->pack_missing;
	c->is_missing	=	(*(c->super))->is_missing;
	c->pack_long	=	(*(c->super))->pack_long;
	c->unpack_long	=	(*(c->super))->unpack_long;
	c->pack_string	=	(*(c->super))->pack_string;
	c->unpack_string	=	(*(c->super))->unpack_string;
	c->pack_bytes	=	(*(c->super))->pack_bytes;
	c->unpack_bytes	=	(*(c->super))->unpack_bytes;
	c->pack_expression	=	(*(c->super))->pack_expression;
	c->notify_change	=	(*(c->super))->notify_change;
	c->update_size	=	(*(c->super))->update_size;
	c->preferred_size	=	(*(c->super))->preferred_size;
	c->resize	=	(*(c->super))->resize;
	c->nearest_smaller_value	=	(*(c->super))->nearest_smaller_value;
	c->next	=	(*(c->super))->next;
	c->compare	=	(*(c->super))->compare;
	c->clear	=	(*(c->super))->clear;
}

/* END_CLASS_IMP */

static void init(grib_accessor* a,const long v, grib_arguments* args)
{
  grib_accessor_data_simple_packing *self =(grib_accessor_data_simple_packing*)a;
  self->units_factor  = grib_arguments_get_name(a->parent->h,args,self->carg++);
  self->units_bias  = grib_arguments_get_name(a->parent->h,args,self->carg++);
  self->changing_precision  = grib_arguments_get_name(a->parent->h,args,self->carg++);
  self->number_of_values  = grib_arguments_get_name(a->parent->h,args,self->carg++);
  self->bits_per_value  = grib_arguments_get_name(a->parent->h,args,self->carg++);
  self->reference_value = grib_arguments_get_name(a->parent->h,args,self->carg++);
  self->binary_scale_factor = grib_arguments_get_name(a->parent->h,args,self->carg++);
  self->decimal_scale_factor = grib_arguments_get_name(a->parent->h,args,self->carg++);
  a->flags |= GRIB_ACCESSOR_FLAG_DATA;
  self->dirty=1;

}

static long value_count(grib_accessor* a){
  grib_accessor_data_simple_packing* self =  (grib_accessor_data_simple_packing*)a;

  long number_of_values;

  if(grib_get_long_internal(a->parent->h,self->number_of_values,&number_of_values) != GRIB_SUCCESS)
       return 0;
  
  return number_of_values;
  
}

static int  unpack_double_element(grib_accessor* a, size_t idx, double* val)
{
  grib_accessor_data_simple_packing* self =  (grib_accessor_data_simple_packing*)a;

  size_t n_vals;
  int err = 0;

  double reference_value;
  long   binary_scale_factor;
  long   bits_per_value;
  long   decimal_scale_factor;
  unsigned char* buf = (unsigned char*)a->parent->h->buffer->data;
  double s = 0;
  double d = 0;
  long pos = 0;
  size_t o = 0;

  n_vals = grib_value_count(a);
  if(n_vals==0){
    return GRIB_NOT_FOUND;
  }

  if((err = grib_get_long_internal(a->parent->h,self->bits_per_value,&bits_per_value)) !=
       GRIB_SUCCESS)
    return err;

  self->dirty=0;

  if((err = grib_get_double_internal(a->parent->h,self->reference_value, &reference_value)) !=
       GRIB_SUCCESS)
    return err;

  if((err = grib_get_long_internal(a->parent->h,self->binary_scale_factor, &binary_scale_factor))
      != GRIB_SUCCESS)
    return err;

  if((err = grib_get_long_internal(a->parent->h,self->decimal_scale_factor, &decimal_scale_factor))
     != GRIB_SUCCESS)
    return err;

      /* Special case */

  if (bits_per_value == 0)  {
    *val=reference_value;
    return GRIB_SUCCESS;
  }

  s = grib_power(binary_scale_factor,2);
  d = grib_power(-decimal_scale_factor,10) ;

  grib_context_log(a->parent->h->context, GRIB_LOG_DEBUG,
      "grib_accessor_data_simple_packing : unpack_double : creating %s, %d values",
      a->name, n_vals);

  buf += grib_byte_offset(a);

  Assert(((bits_per_value*n_vals)/8) < (1<<29));
  /*ensure that the bit pointer is not overflown*/

  if(bits_per_value%8)
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_DEBUG,
        "unpack_double : calling outline function : bpv %d, rv : %g, sf : %d, dsf : %d ",
        bits_per_value,reference_value,binary_scale_factor, decimal_scale_factor);
    pos=idx*bits_per_value;
    *val= (double) (((
       grib_decode_unsigned_long(buf, &pos, bits_per_value)*s)+reference_value)*d);
    /* val[i] = grib_decode_unsigned_long(buf, &pos, bits_per_value);                                   */
    /* fprintf(stdout,"unpck uuu-o: %d vals %d bitspv buf %d by long \n", n_vals, bits_per_value, pos/8);*/
  }
  else
  {
    int bc;
    long lvalue = 0;
    int l = bits_per_value/8;

    pos=idx*l;
    buf+=pos;
    lvalue  = 0;
    lvalue  <<= 8;
    lvalue |= buf[o++] ;

    for ( bc=1; bc<l; bc++ ) {
        lvalue <<= 8;
        lvalue |= buf[o++] ;
    }
    *val = (double) (((lvalue*s)+reference_value)*d);

   }

  return err;
}

static int  _unpack_double(grib_accessor* a, double* val, size_t *len,unsigned char* buf,long pos, size_t n_vals)
{
  grib_accessor_data_simple_packing* self =  (grib_accessor_data_simple_packing*)a;

  size_t i = 0;
  int err = 0;

  double reference_value;
  long   binary_scale_factor;
  long   bits_per_value;
  long   decimal_scale_factor;
  double s = 0;
  double d = 0;
  double units_factor=1.0;
  double units_bias=0.0;


  if(*len < n_vals)
  {
    *len = (long) n_vals;
    return GRIB_ARRAY_TOO_SMALL;
  }

  if((err = grib_get_long_internal(a->parent->h,self->bits_per_value,&bits_per_value)) !=
       GRIB_SUCCESS)
    return err;

  /*/
   * check we don't decode bpv > max(ulong) as it is 
   * not currently supported by the algorithm 
   */
  if ( bits_per_value > (sizeof(long)*8) ) {
    return GRIB_INVALID_BPV;
  }

  if(self->units_factor &&
	 (grib_get_double_internal(a->parent->h,self->units_factor,&units_factor)== GRIB_SUCCESS)) {
		  grib_set_double_internal(a->parent->h,self->units_factor,1.0);
  }

  if(self->units_bias &&
	 (grib_get_double_internal(a->parent->h,self->units_bias,&units_bias)== GRIB_SUCCESS)) {
		  grib_set_double_internal(a->parent->h,self->units_bias,0.0);
	}

  if(n_vals==0){
    *len = 0;
    return GRIB_SUCCESS;
  }

  self->dirty=0;

  if((err = grib_get_double_internal(a->parent->h,self->reference_value, &reference_value)) !=
       GRIB_SUCCESS)
    return err;

  if((err = grib_get_long_internal(a->parent->h,self->binary_scale_factor, &binary_scale_factor))
      != GRIB_SUCCESS)
    return err;

  if((err = grib_get_long_internal(a->parent->h,self->decimal_scale_factor, &decimal_scale_factor))
     != GRIB_SUCCESS)
    return err;

      /* Special case */

  if(bits_per_value == 0)
  {
    for(i = 0; i < n_vals; i++)
      val[i] = reference_value;
    *len = n_vals;
    return GRIB_SUCCESS;
  }

  s = grib_power(binary_scale_factor,2);
  d = grib_power(-decimal_scale_factor,10) ;

  grib_context_log(a->parent->h->context, GRIB_LOG_DEBUG,
      "grib_accessor_data_simple_packing : unpack_double : creating %s, %d values",
      a->name, n_vals);

  buf += grib_byte_offset(a);

  Assert(((bits_per_value*n_vals)/8) < (1<<29));
  /*ensure that the bit pointer is not overflown*/

  grib_context_log(a->parent->h->context, GRIB_LOG_DEBUG,
        "unpack_double : calling outline function : bpv %d, rv : %g, sf : %d, dsf : %d ",
        bits_per_value,reference_value,binary_scale_factor, decimal_scale_factor);
  grib_decode_double_array(buf,&pos,bits_per_value,reference_value,s,d,n_vals,val);

  *len = (long) n_vals;

  if (units_factor != 1.0) {
	if (units_bias != 0.0)
		for (i=0;i<n_vals;i++) val[i]=val[i]*units_factor+units_bias;
	else 
		for (i=0;i<n_vals;i++) val[i]*=units_factor;
  } else if (units_bias != 0.0)
		for (i=0;i<n_vals;i++) val[i]+=units_bias;

  return err;
}

static int  unpack_double_subarray(grib_accessor* a, double* val, size_t start, size_t len)
{
  grib_accessor_data_simple_packing* self =  (grib_accessor_data_simple_packing*)a;
  unsigned char* buf = (unsigned char*)a->parent->h->buffer->data;
  size_t nvals = len;
  size_t *plen=&len;
  long bits_per_value=0;
  long pos;
  int err;

  if((err = grib_get_long_internal(a->parent->h,self->bits_per_value,&bits_per_value)) !=
       GRIB_SUCCESS)
    return err;

  buf+=(start*bits_per_value)/8;
  pos=start*bits_per_value%8;
  return _unpack_double(a,val,plen,buf,pos,nvals);
}

static int  unpack_double(grib_accessor* a, double* val, size_t *len) {
  unsigned char* buf = (unsigned char*)a->parent->h->buffer->data;
  size_t nvals = grib_value_count(a);
  long pos=0;
  return _unpack_double(a,val,len,buf,pos,nvals);
}


static int pack_double(grib_accessor* a, const double* val, size_t *len)
{
  grib_accessor_data_simple_packing* self =  (grib_accessor_data_simple_packing*)a;

  size_t i = 0;
  size_t n_vals = *len;
  int err = 0;
  int last;
  double reference_value = 0;
  long   binary_scale_factor = 0;
  long   bits_per_value = 0;
  long   decimal_scale_factor = 0;
  long   decimal_scale_factor_get = 0;
  double decimal = 1;
  double max = 0;
  double min = 0;
  double unscaled_max = 0;
  double unscaled_min = 0;
  long imin,imax;
  double f=0;
  double range=0;
  double minrange=0,maxrange=0;
  long changing_precision=0;
  grib_context* c=a->parent->h->context;

  decimal_scale_factor=0;

  if (*len ==0) return GRIB_NO_VALUES;

  if((err = grib_get_long_internal(a->parent->h,self->bits_per_value,&bits_per_value)) !=
      GRIB_SUCCESS)
    return err;

  if(*len == 0) return GRIB_SUCCESS;

  if((err = grib_get_long_internal(a->parent->h,self->decimal_scale_factor, &decimal_scale_factor_get))
     != GRIB_SUCCESS)
     return err;
  /*/
   * check we don't encode bpv > max(ulong)-1 as it is 
   * not currently supported by the algorithm 
   */
  if ( bits_per_value > (sizeof(long)*8-1) ) {
    return GRIB_INVALID_BPV;
  }

  self->dirty=1;

  max = val[0];
  min = max;
  for(i=1;i< n_vals;i++) {
    if (val[i] > max ) max = val[i];
    if (val[i] < min ) min = val[i];
  }

  /* constant field only reference_value is set and bits_per_value=0 */
  if(max==min) {
	if (grib_get_nearest_smaller_value(a->parent->h,self->reference_value,val[0],&reference_value)
		!=GRIB_SUCCESS) {
		grib_context_log(a->parent->h->context,GRIB_LOG_ERROR,
		 "unable to find nearest_smaller_value of %g for %s",min,self->reference_value);
		exit(GRIB_INTERNAL_ERROR);
	}
	if((err = grib_set_double_internal(a->parent->h,self->reference_value, reference_value)) !=
	  GRIB_SUCCESS)
	  return err;

	{
	  /* Make sure we can decode it again */
	  double ref = 1e-100;
	  grib_get_double_internal(a->parent->h,self->reference_value,&ref);
	  if (ref != reference_value) 
		printf("%.20e  !=  %.20e",ref,reference_value);
	  Assert(ref == reference_value);
	}

	if (c->gribex_mode_on==0) {
		bits_per_value=0;
		if((err = grib_set_long_internal(a->parent->h,self->bits_per_value, bits_per_value)) !=
		  GRIB_SUCCESS)
		  return err;

		 return GRIB_CONSTANT_FIELD;
	 } else {
		if((err = grib_set_long_internal(a->parent->h,self->binary_scale_factor, 0)) !=
		  GRIB_SUCCESS)
		  return err;

		if((err = grib_set_long_internal(a->parent->h,self->decimal_scale_factor, 0)) !=
		  GRIB_SUCCESS)
		  return err;

	 	return GRIB_SUCCESS;
	 }
  }

  if((err = grib_get_long_internal(a->parent->h,self->binary_scale_factor, &binary_scale_factor))
      != GRIB_SUCCESS)
    return err;

  if((err = grib_get_long_internal(a->parent->h,self->changing_precision, &changing_precision))
      != GRIB_SUCCESS)
    return err;

  /* the packing parameters are not properly defined
  this is a safe way of fixing the problem */
  if ( changing_precision==0 && bits_per_value==0  && decimal_scale_factor_get==0) {

    grib_context_log(a->parent->h->context,GRIB_LOG_WARNING,
                     "%s==0 and %s==0 (setting %s=24)",
                     self->bits_per_value,
                     self->decimal_scale_factor,
                     self->bits_per_value);
    
    bits_per_value=24;
    if((err = grib_set_long_internal(a->parent->h,self->bits_per_value,
        bits_per_value))!= GRIB_SUCCESS)
      return err;
  }

  if ( bits_per_value == 0 ) {
    /* decimal_scale_factor is given, binary_scale_factor=0
       and bits_per_value is computed */
    binary_scale_factor=0;
	decimal_scale_factor=decimal_scale_factor_get;
    decimal = grib_power(decimal_scale_factor,10) ;
    min*=decimal;
    max*=decimal;

    imin=(long)rint(min);
    imax=(long)rint(max);
    bits_per_value=(long)ceil(log((double)(imax-imin+1))/log(2.0));
    /*printf("bits_per_value=%ld\n",bits_per_value);*/
    if((err = grib_set_long_internal(a->parent->h,self->bits_per_value, bits_per_value)) !=
      GRIB_SUCCESS)
      return err;
    if (grib_get_nearest_smaller_value(a->parent->h,self->reference_value,min,&reference_value)
        !=GRIB_SUCCESS) {
        grib_context_log(a->parent->h->context,GRIB_LOG_ERROR,
         "unable to find nearest_smaller_value of %g for %s",min,self->reference_value);
        exit(GRIB_INTERNAL_ERROR);
    }
    /* divisor=1; */
  } else {
    last=127;
	if (c->gribex_mode_on) last=99;
    /* bits_per_value is given and decimal_scale_factor
       and binary_scale_factor are calcualated
    */
    if (max == min) {
      binary_scale_factor=0;
      /* divisor=1; */
      if (grib_get_nearest_smaller_value(a->parent->h,self->reference_value,min,&reference_value)
        !=GRIB_SUCCESS) {
        grib_context_log(a->parent->h->context,GRIB_LOG_ERROR,
         "unable to find nearest_smaller_value of %g for %s",min,self->reference_value);
        exit(GRIB_INTERNAL_ERROR);
      }
    } else {
     /* printf("max=%g reference_value=%g grib_power(-last,2)=%g decimal_scale_factor=%ld bits_per_value=%ld\n",
            max,reference_value,grib_power(-last,2),decimal_scale_factor,bits_per_value);*/
      /* last must be a parameter coming from the def file*/
      range=(max-min);
	  unscaled_min=min;
	  unscaled_max=max;
      f=(grib_power(bits_per_value,2)-1);
      minrange=grib_power(-last,2)*f;
      maxrange=grib_power(last,2)*f;

      while (range<minrange) {
        decimal_scale_factor+=1;
        decimal*=10;
        min=unscaled_min*decimal;
        max=unscaled_max*decimal;
        range=(max-min);
      }
      while (range>maxrange) {
        decimal_scale_factor-=1;
        decimal/=10;
        min=unscaled_min*decimal;
        max=unscaled_max*decimal;
        range=(max-min);
      }

      if (grib_get_nearest_smaller_value(a->parent->h,self->reference_value,
              min,&reference_value)!=GRIB_SUCCESS) {
        grib_context_log(a->parent->h->context,GRIB_LOG_ERROR,
           "unable to find nearest_smaller_value of %g for %s",min,self->reference_value);
        exit(GRIB_INTERNAL_ERROR);
      }

      binary_scale_factor = grib_get_binary_scale_fact(max,reference_value,bits_per_value,&err);
    }
  }


  if((err = grib_set_double_internal(a->parent->h,self->reference_value, reference_value)) !=
      GRIB_SUCCESS)
    return err;

  if((err = grib_set_long_internal(a->parent->h,self->changing_precision, 0)) !=
     GRIB_SUCCESS)
    return err;
  if((err = grib_set_long_internal(a->parent->h,self->binary_scale_factor, binary_scale_factor)) !=
     GRIB_SUCCESS)
    return err;
  if((err = grib_set_long_internal(a->parent->h,self->decimal_scale_factor, decimal_scale_factor))
     != GRIB_SUCCESS)
    return err;


  return GRIB_SUCCESS;

}

