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
   IMPLEMENTS = destroy
   IMPLEMENTS = native_type
   IMPLEMENTS = evaluate_long
   IMPLEMENTS = evaluate_double
   IMPLEMENTS = print
   IMPLEMENTS = compile
   IMPLEMENTS = add_dependency
   MEMBERS    = double value
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

static int        evaluate_long(grib_expression*,grib_handle*,long*);
static int      evaluate_double(grib_expression*,grib_handle*,double*);

typedef struct grib_expression_double{
  grib_expression base;
/* Members defined in double */
	double value;
} grib_expression_double;


static grib_expression_class _grib_expression_class_double = {
    0,                    /* super                     */
    "double",                    /* name                      */
    sizeof(grib_expression_double),/* size of instance          */
    0,                           /* inited */
    &init_class,                 /* init_class */
    0,                     /* constructor               */
    &destroy,                  /* destructor                */
    &print,                 
    &compile,                 
    &add_dependency,       

	&native_type,
	0,

	&evaluate_long,
	&evaluate_double,
	0,
};

grib_expression_class* grib_expression_class_double = &_grib_expression_class_double;


static void init_class(grib_expression_class* c)
{
}
/* END_CLASS_IMP */

static int evaluate_long(grib_expression* g,grib_handle *h,long* lres)
{
	grib_expression_double* e = (grib_expression_double*)g;
	*lres=e->value;
	return GRIB_SUCCESS;
}

static int evaluate_double(grib_expression* g,grib_handle *h,double* dres)
{
	grib_expression_double* e = (grib_expression_double*)g;
	*dres=e->value;
	return GRIB_SUCCESS;
}

static void print(grib_context* c,grib_expression* g,grib_handle* f)
{
	grib_expression_double* e = (grib_expression_double*)g;
	printf("double(%g)",e->value);
}

static void destroy(grib_context* c,grib_expression* g)
{
	/* grib_expression_double* e = (grib_expression_double*)g; */
}


static void  add_dependency(grib_expression* g, grib_accessor* observer){
	/* grib_expression_double* e = (grib_expression_double*)g; */
}


grib_expression* new_double_expression(grib_context* c,double value)
{
	grib_expression_double* e = grib_context_malloc_clear_persistent(c,sizeof(grib_expression_double));
	e->base.cclass  = grib_expression_class_double;
	e->value = value;
	return (grib_expression*)e;
}

static void compile(grib_expression* g,grib_compiler* c)
{
	grib_expression_double* e = (grib_expression_double*)g;
    fprintf(c->out,"new_double_expression(ctx,%g)",e->value);
}


static int native_type(grib_expression* g,grib_handle *h)
{
	return GRIB_TYPE_DOUBLE;
}

