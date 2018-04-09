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
   IMPLEMENTS = print
   IMPLEMENTS = compile
   IMPLEMENTS = add_dependency
   MEMBERS    = char *name
   MEMBERS    = grib_arguments *args
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

typedef struct grib_expression_functor{
  grib_expression base;
/* Members defined in functor */
	char *name;
	grib_arguments *args;
} grib_expression_functor;


static grib_expression_class _grib_expression_class_functor = {
    0,                    /* super                     */
    "functor",                    /* name                      */
    sizeof(grib_expression_functor),/* size of instance          */
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
	0,
	0,
};

grib_expression_class* grib_expression_class_functor = &_grib_expression_class_functor;


static void init_class(grib_expression_class* c)
{
}
/* END_CLASS_IMP */

static int evaluate_long(grib_expression* g,grib_handle* h,long* lres)
{
	grib_expression_functor* e = (grib_expression_functor*)g;
  
	/*
	TODO: needs OO code here
	*/

	if(strcmp(e->name,"lookup") == 0) {
		return GRIB_SUCCESS;
	}
	
	if(strcmp(e->name,"new") == 0) {
		*lres=h->loader != NULL;
		return GRIB_SUCCESS;
	}
		
	if(strcmp(e->name,"missing") == 0)
	{
		const char *p = grib_arguments_get_name(h,e->args,0);

		if(p)
		{

			long val = 0;
			grib_get_long_internal(h,p,&val);
			*lres = (val == GRIB_MISSING_LONG);
			return GRIB_SUCCESS;
		}
		else
			*lres=GRIB_MISSING_LONG;
			return GRIB_SUCCESS;
	}

  if(strcmp(e->name,"defined") == 0)
  {
    const char *p = grib_arguments_get_name(h,e->args,0);

    if(p)  {
      grib_accessor* a=grib_find_accessor(h,p);
      *lres=a!=NULL ? 1 : 0;
      return GRIB_SUCCESS;
    }
    *lres=0;
    return GRIB_SUCCESS;
  }

  if(strcmp(e->name,"changed") == 0)
  {
    *lres=1;
    return GRIB_SUCCESS;
  }

  if(strcmp(e->name,"gribex_mode_on") == 0)
  {
    *lres= h->context->gribex_mode_on ? 1 : 0;
    return GRIB_SUCCESS;
  }
		
	return GRIB_NOT_IMPLEMENTED;
}

static void print(grib_context* c,grib_expression* g,grib_handle* f)
{
	grib_expression_functor* e = (grib_expression_functor*)g;
	printf("%s(",e->name);
	/*grib_expression_print(c,e->args,f);*/
	printf(")");
}

static void destroy(grib_context* c,grib_expression* g)
{
	grib_expression_functor* e = (grib_expression_functor*)g;
	grib_context_free_persistent(c,e->name);
	grib_arguments_free(c,e->args);
}


static void  add_dependency(grib_expression* g, grib_accessor* observer){
	grib_expression_functor* e = (grib_expression_functor*)g;
	if (strcmp(e->name,"defined")) 
		grib_dependency_observe_arguments(observer,e->args);
}

grib_expression* new_func_expression(grib_context* c,const char *name,grib_arguments* args)
{
	grib_expression_functor* e = grib_context_malloc_clear_persistent(c,sizeof(grib_expression_functor));
	e->base.cclass                 = grib_expression_class_functor;
	e->name                   = grib_context_strdup_persistent(c,name);
	e->args                  = args;
	return (grib_expression*)e;
}

static void compile(grib_expression* g,grib_compiler* c)
{
	grib_expression_functor* e = (grib_expression_functor*)g;
    fprintf(c->out,"new_func_expression(ctx,");
    fprintf(c->out,"\"%s\",",e->name);
    grib_compile_arguments(e->args,c);
    fprintf(c->out,")");
}

static int native_type(grib_expression* g,grib_handle *h)
{
	return GRIB_TYPE_LONG;
}
