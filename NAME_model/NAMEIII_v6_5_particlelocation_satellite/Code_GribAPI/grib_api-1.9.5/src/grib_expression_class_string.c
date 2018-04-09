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
   CLASS      = expression
   IMPLEMENTS = init_class
   IMPLEMENTS = native_type
   IMPLEMENTS = destroy
   IMPLEMENTS = evaluate_string
   IMPLEMENTS = print
   IMPLEMENTS = compile
   IMPLEMENTS = add_dependency
   MEMBERS    = char* value
   END_CLASS_DEF

 */
/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "expression.class" and rerun ./make_class.pl

*/

typedef const char* string; /* to keep make_class.pl happy */


static void init_class              (grib_expression_class*);

static void        destroy(grib_context*,grib_expression* e);

static void        print(grib_context*,grib_expression*,grib_handle*);
static void        compile(grib_expression*,grib_compiler*);
static void        add_dependency(grib_expression* e, grib_accessor* observer);

static int        native_type(grib_expression*,grib_handle*);

static string evaluate_string(grib_expression*,grib_handle*,char*,size_t*,int*);

typedef struct grib_expression_string{
  grib_expression base;
/* Members defined in string */
	char* value;
} grib_expression_string;


static grib_expression_class _grib_expression_class_string = {
    0,                    /* super                     */
    "string",                    /* name                      */
    sizeof(grib_expression_string),/* size of instance          */
    0,                           /* inited */
    &init_class,                 /* init_class */
    0,                     /* constructor               */
    &destroy,                  /* destructor                */
    &print,                 
    &compile,                 
    &add_dependency,       

	&native_type,
	0,

	0,
	0,
	&evaluate_string,
};

grib_expression_class* grib_expression_class_string = &_grib_expression_class_string;


static void init_class(grib_expression_class* c)
{
}
/* END_CLASS_IMP */

static const char* evaluate_string(grib_expression* g,grib_handle *h,char* buf,size_t* size,int* err)
{
	grib_expression_string* e = (grib_expression_string*)g;
	*err=0;
	return e->value;
}

static void print(grib_context* c,grib_expression* g,grib_handle* f)
{
	grib_expression_string* e = (grib_expression_string*)g;
	printf("string('%s')",e->value);
}

static void destroy(grib_context* c,grib_expression* g)
{
	grib_expression_string* e = (grib_expression_string*)g;
	grib_context_free_persistent(c,e->value);
}


static void  add_dependency(grib_expression* g, grib_accessor* observer){
	/* grib_expression_string* e = (grib_expression_string*)g; */
}


grib_expression* new_string_expression(grib_context* c,const char* value)
{
	grib_expression_string* e = grib_context_malloc_clear_persistent(c,sizeof(grib_expression_string));
	e->base.cclass                 = grib_expression_class_string;
	e->value               = grib_context_strdup_persistent(c,value);
	return (grib_expression*)e;
}

static void compile(grib_expression* g,grib_compiler* c)
{
	grib_expression_string* e = (grib_expression_string*)g;
    fprintf(c->out,"new_string_expression(ctx,\"%s\")",e->value);
}


static int native_type(grib_expression* g,grib_handle *h)
{
	return GRIB_TYPE_STRING;
}
