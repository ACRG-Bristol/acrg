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
   IMPLEMENTS = print;compile
   IMPLEMENTS = add_dependency
   MEMBERS    = grib_expression *left
   MEMBERS = grib_expression *right
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

typedef struct grib_expression_string_compare{
  grib_expression base;
/* Members defined in string_compare */
	grib_expression *left;
	grib_expression *right;
} grib_expression_string_compare;


static grib_expression_class _grib_expression_class_string_compare = {
    0,                    /* super                     */
    "string_compare",                    /* name                      */
    sizeof(grib_expression_string_compare),/* size of instance          */
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

grib_expression_class* grib_expression_class_string_compare = &_grib_expression_class_string_compare;


static void init_class(grib_expression_class* c)
{
}
/* END_CLASS_IMP */

GRIB_INLINE static int grib_inline_strcmp(const char* a,const char* b) {
  if (*a != *b) return 1;
  while((*a!=0 && *b!=0) &&  *(a) == *(b) ) {a++;b++;}
  return (*a==0 && *b==0) ? 0 : 1;
}

static int evaluate_long(grib_expression *g,grib_handle* h,long* lres)
{
  int ret=0;
  char b1[1024]; size_t l1 = sizeof(b1);
  char b2[1024]; size_t l2 = sizeof(b2);
  const char *v1;
  const char *v2;

  grib_expression_string_compare* e = (grib_expression_string_compare*)g;

  v1 = grib_expression_evaluate_string(h,e->left,b1,&l1,&ret);
  if(!v1) {*lres=0;return ret;}

  v2 = grib_expression_evaluate_string(h,e->right,b2,&l2,&ret);
  if(!v2) {*lres=0;return ret;}

  *lres = (grib_inline_strcmp(v1,v2) == 0);
  return GRIB_SUCCESS;
}

static int evaluate_double(grib_expression *g,grib_handle* h,double* dres)
{
  long n;
  int ret = evaluate_long(g,h,&n);
  *dres = n;
  return ret;
}

static void print(grib_context* c,grib_expression* g,grib_handle* f)
{
  grib_expression_string_compare* e = (grib_expression_string_compare*)g;
  printf("string_compare(");
  grib_expression_print(c,e->left,f);
  printf(",");
  grib_expression_print(c,e->right,f);
  printf(")");
}

static void destroy(grib_context* c,grib_expression* g)
{
  grib_expression_string_compare* e = (grib_expression_string_compare*)g;
  grib_expression_free(c,e->left);
  grib_expression_free(c,e->right);
}


static void  add_dependency(grib_expression* g, grib_accessor* observer){
  grib_expression_string_compare* e = (grib_expression_string_compare*)g;
  grib_dependency_observe_expression(observer,e->left);
  grib_dependency_observe_expression(observer,e->right);
}

grib_expression* new_string_compare_expression(grib_context* c,
  grib_expression* left,grib_expression* right)
{
  grib_expression_string_compare* e = grib_context_malloc_clear_persistent(c,sizeof(grib_expression_string_compare));
  e->base.cclass                 = grib_expression_class_string_compare;
  e->left                = left;
  e->right               = right;
  return (grib_expression*)e;
}

static void compile(grib_expression* g,grib_compiler* c)
{
	grib_expression_string_compare* e = (grib_expression_string_compare*)g;
    fprintf(c->out,"new_string_compare_expression(ctx,");
    grib_expression_compile(e->left,c);
    fprintf(c->out,",");
    grib_expression_compile(e->right,c);
    fprintf(c->out,")");
}

static int native_type(grib_expression* g,grib_handle *h)
{
  return GRIB_TYPE_LONG ;
}
