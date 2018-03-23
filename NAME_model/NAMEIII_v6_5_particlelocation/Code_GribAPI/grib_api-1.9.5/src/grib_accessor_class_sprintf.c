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
   SUPER      = grib_accessor_class_ascii
   IMPLEMENTS = pack_string;unpack_string;value_count
   IMPLEMENTS = init
   MEMBERS= grib_arguments* args  
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int pack_string(grib_accessor*, const char*, size_t *len);
static int unpack_string (grib_accessor*, char*, size_t *len);
static long value_count(grib_accessor*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_sprintf {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in ascii */
/* Members defined in sprintf */
	grib_arguments* args;
} grib_accessor_sprintf;

extern grib_accessor_class* grib_accessor_class_ascii;

static grib_accessor_class _grib_accessor_class_sprintf = {
    &grib_accessor_class_ascii,                      /* super                     */
    "sprintf",                      /* name                      */
    sizeof(grib_accessor_sprintf),  /* size                      */
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
    0,                /* grib_pack procedures double    */
    0,              /* grib_unpack procedures double  */
    &pack_string,                /* grib_pack procedures string    */
    &unpack_string,              /* grib_unpack procedures string  */
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
    0,     /* unpack only ith value          */
    0,     /* unpack a subarray         */
    0,             		/* clear          */
};


grib_accessor_class* grib_accessor_class_sprintf = &_grib_accessor_class_sprintf;


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
	c->pack_double	=	(*(c->super))->pack_double;
	c->unpack_double	=	(*(c->super))->unpack_double;
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
	c->unpack_double_element	=	(*(c->super))->unpack_double_element;
	c->unpack_double_subarray	=	(*(c->super))->unpack_double_subarray;
	c->clear	=	(*(c->super))->clear;
}

/* END_CLASS_IMP */

static void init(grib_accessor* a,const long l, grib_arguments* c)
{
	grib_accessor_sprintf* self = (grib_accessor_sprintf*)a; 
	self->args = c;

}

static int pack_string(grib_accessor* a, const char* val, size_t *len){

	return GRIB_NOT_IMPLEMENTED;

}


static int    unpack_string(grib_accessor* a, char* val, size_t *len)
{   
	grib_accessor_sprintf* self = (grib_accessor_sprintf*)a; 

	char result[1024]  ;
	char sres[1024]  ;
	long ires = 0;
	double dres= 0;

	int  i = 0;
	size_t replen = 1024;

	int ret = GRIB_SUCCESS;

	int carg= 0;

	const char* uname = NULL;
	const char* tempname = NULL;


	uname = grib_arguments_get_string(a->parent->h,self->args,carg++);
	sprintf(result,"%s",""); 

	for(i=0;i<strlen(uname);i++)
	{

		if(uname[i]=='%'){
			i++;
			switch(uname[i]){

				case 'd':
					tempname = grib_arguments_get_name(a->parent->h,self->args,carg++);

					if((ret = grib_get_long_internal(a->parent->h,tempname,&ires)) != GRIB_SUCCESS)
						return ret;
					sprintf(result,"%s%ld",result, ires); 


					break;

				case 'g':
					tempname = grib_arguments_get_name(a->parent->h,self->args,carg++);
					if((ret = grib_get_double_internal(a->parent->h,tempname,&dres)) != GRIB_SUCCESS) 
						return ret;
					sprintf(result,"%s%g",result, dres); 


					break;

				case 's':
					tempname = grib_arguments_get_name(a->parent->h,self->args,carg++);
					if((ret = grib_get_string_internal(a->parent->h,tempname,sres, &replen )) != GRIB_SUCCESS)
						return ret;
					sprintf(result,"%s%s",result, sres); 
					replen = 1024;


					break;


			}
		} 
		else
			sprintf(result,"%s%c",result, uname[i]); 


	}

	replen = strlen(result)+1;

	if(*len < replen){
		*len = replen;
		return GRIB_ARRAY_TOO_SMALL;
	}
	*len = replen;

	sprintf(val,"%s",result);    

	return GRIB_SUCCESS;
}


static long value_count(grib_accessor* a)
{
	char result[1024]  ;
	size_t s = sizeof(result);

	unpack_string(a,result,&s);

    return s;
}
