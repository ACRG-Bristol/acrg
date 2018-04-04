/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/***************************************************************************
 * Enrico Fucile                                                                         *
 ***************************************************************************/
#include "grib_api_internal.h"
/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = action
   IMPLEMENTS = create_accessor
   IMPLEMENTS = dump;xref;compile
   IMPLEMENTS = destroy;notify_change
   MEMBERS    = grib_expression *expression
   MEMBERS    = grib_action     *block_true
   MEMBERS    = grib_action     *block_false
   MEMBERS    = int loop
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "action.class" and rerun ./make_class.pl

*/

static void init_class      (grib_action_class*);
static void dump            (grib_action* d, FILE*,int);
static void xref            (grib_action* d, FILE* f,const char* path);
static void compile         (grib_action* a, grib_compiler* compiler);
static void destroy         (grib_context*,grib_action*);
static int create_accessor(grib_section*,grib_action*,grib_loader*);
static int notify_change(grib_action* a, grib_accessor* observer,grib_accessor* observed);


typedef struct grib_action_when {
    grib_action          act;  
/* Members defined in when */
	grib_expression *expression;
	grib_action     *block_true;
	grib_action     *block_false;
	int loop;
} grib_action_when;


static grib_action_class _grib_action_class_when = {
    0,                              /* super                     */
    "action_class_when",                              /* name                      */
    sizeof(grib_action_when),            /* size                      */
    0,                                   /* inited */
    &init_class,                         /* init_class */
    0,                               /* init                      */
    &destroy,                            /* destroy */

    &dump,                               /* dump                      */
    &xref,                               /* xref                      */

    &create_accessor,             /* create_accessor*/

    &notify_change,                            /* notify_change */
    0,                            /* reparse */
    0,                            /* execute */
    &compile,                            /* compile */
};

grib_action_class* grib_action_class_when = &_grib_action_class_when;

static void init_class(grib_action_class* c)
{
}
/* END_CLASS_IMP */

grib_action* grib_action_create_when( grib_context* context,
    grib_expression* expression,
    grib_action* block_true,grib_action* block_false)
{
  char name[1024];

  grib_action_when* a ;
  grib_action_class* c   = grib_action_class_when;
  grib_action* act       = (grib_action*)grib_context_malloc_clear_persistent(context,c->size);
  act->op              = grib_context_strdup_persistent(context,"when");

  act->cclass       = c;
  a                 = (grib_action_when*)act;
  act->context      = context;

  a->expression  = expression;
  a->block_true       = block_true;
  a->block_false      = block_false;


  sprintf(name,"_when%p",(void*)expression);

  act->name      = grib_context_strdup_persistent(context,name);

  return act;
}

static void compile(grib_action* act, grib_compiler *compiler)
{
    grib_action_when* a  = (grib_action_when*)act;
    char t[80];
    char f[80];

    if(a->block_true)
        grib_compile_action_branch(a->block_true, compiler,t); 
    else
        strcpy(t,"NULL");

    if(a->block_false)
        grib_compile_action_branch(a->block_false, compiler,f);
    else
        strcpy(f,"NULL");

    fprintf(compiler->out,"%s = grib_action_create_when(ctx,",compiler->var);
    grib_compile_expression(a->expression, compiler);
    fprintf(compiler->out,",%s,%s);\n", t,f);
}


static int create_accessor(grib_section* p, grib_action* act,grib_loader *h)
{
  grib_action_when* self = (grib_action_when*)act;
  grib_accessor*  as = grib_accessor_factory(p, act,0,0);
  if(!as)return GRIB_INTERNAL_ERROR;

  grib_dependency_observe_expression(as,self->expression);

  grib_push_accessor(as,p->block);

  return GRIB_SUCCESS;
}


static void dump(grib_action* act, FILE* f, int lvl)
{
  grib_action_when* a = (grib_action_when*)act;
  int i = 0;

  for (i=0;i<lvl;i++)
    grib_context_print(act->context,f,"     ");

  printf("when(%s) { ",act->name);
  grib_expression_print(act->context,a->expression,0);
  printf("\n");

  grib_dump_action_branch(f,a->block_true,lvl+1);

  for (i=0;i<lvl;i++)
    grib_context_print(act->context,f,"     ");
  printf("}");
  
  if (a->block_false) {
    printf(" else { ");
    
    grib_dump_action_branch(f,a->block_true,lvl+1);
    
    for (i=0;i<lvl;i++)
      grib_context_print(act->context,f,"     ");
    printf("}");
  }
  printf("\n");
}

static int notify_change(grib_action* a, grib_accessor* observer,grib_accessor* observed)
{
  grib_action_when* self = (grib_action_when*) a;
  grib_action *b = NULL;
  int ret = GRIB_SUCCESS;
  long lres;

  if ((ret = grib_expression_evaluate_long(observed->parent->h, self->expression,&lres))
    != GRIB_SUCCESS) return ret;

  if(0 && self->loop)
  {
    printf("LOOP detected...\n");
    printf("WHEN triggered by %s %ld\n",observed->name,lres);
    grib_expression_print(observed->parent->h->context,self->expression,0);
    printf("\n");
    return ret;
  }

  self->loop = 1;

  if(lres) 
    b=self->block_true;
  else
    b=self->block_false;

  while(b) {
    ret = grib_action_execute(b,observed->parent->h);
    if(ret != GRIB_SUCCESS) {
      self->loop = 0;
      return ret;
    }
    b = b->next;
  }

  self->loop = 0;

  return GRIB_SUCCESS;
}

static void destroy(grib_context* context,grib_action* act)
{
  grib_action_when* self = (grib_action_when*) act;
  grib_action *t = self->block_true;


  while(t)
  {
    grib_action *nt = t->next;
    grib_free_action(context,t);
    t = nt;
  }

  t=self->block_false;
  while(t)
  {
    grib_action *nt = t->next;
    grib_free_action(context,t);
    t = nt;
  }

  grib_expression_free(context,self->expression);

  grib_context_free_persistent(context, act->name);
  grib_context_free_persistent(context, act->op);
}

static void xref(grib_action* d, FILE* f,const char* path)
{
}
