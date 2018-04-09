/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/***************************************************************************
 *   Jean Baptiste Filippi - 01.11.2005                                                           *
 *   Enrico Fucile
 *                                                                         *
 ***************************************************************************/
#include "grib_api_internal.h"
/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = action
   SUPER      = action_class_section
   IMPLEMENTS = create_accessor
   IMPLEMENTS = dump
   IMPLEMENTS = destroy
   IMPLEMENTS = xref   
   IMPLEMENTS = compile
   IMPLEMENTS = reparse;execute
   MEMBERS    = grib_expression *expression
   MEMBERS    = grib_action     *block_true
   MEMBERS    = grib_action     *block_false
   MEMBERS    = int transient
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
static grib_action* reparse(grib_action* a,grib_accessor* acc,int *doit);
static int execute(grib_action* a,grib_handle* h);


typedef struct grib_action_if {
    grib_action          act;  
/* Members defined in section */
/* Members defined in if */
	grib_expression *expression;
	grib_action     *block_true;
	grib_action     *block_false;
	int transient;
} grib_action_if;

extern grib_action_class* grib_action_class_section;

static grib_action_class _grib_action_class_if = {
    &grib_action_class_section,                              /* super                     */
    "action_class_if",                              /* name                      */
    sizeof(grib_action_if),            /* size                      */
    0,                                   /* inited */
    &init_class,                         /* init_class */
    0,                               /* init                      */
    &destroy,                            /* destroy */

    &dump,                               /* dump                      */
    &xref,                               /* xref                      */

    &create_accessor,             /* create_accessor*/

    0,                            /* notify_change */
    &reparse,                            /* reparse */
    &execute,                            /* execute */
    &compile,                            /* compile */
};

grib_action_class* grib_action_class_if = &_grib_action_class_if;

static void init_class(grib_action_class* c)
{
	c->notify_change	=	(*(c->super))->notify_change;
}
/* END_CLASS_IMP */

grib_action* grib_action_create_if( grib_context* context,
    grib_expression* expression,
    grib_action* block_true,grib_action* block_false,int transient)
{   char name[1024];
  grib_action_if* a ;
  grib_action_class* c   = grib_action_class_if;
  grib_action* act       = (grib_action*)grib_context_malloc_clear_persistent(context,c->size);
  act->op              = grib_context_strdup_persistent(context,"section");

  act->cclass       = c;
  a                 = (grib_action_if*)act;
  act->context      = context;

  a->expression  = expression;
  a->block_true  = block_true;
  a->block_false = block_false;
  a->transient   = transient;

  if (transient)
	  sprintf(name,"__if%p",(void*)a);
  else
	  sprintf(name,"_if%p",(void*)a);

  act->name      = grib_context_strdup_persistent(context,name);

  return act;
}

static void compile(grib_action* act, grib_compiler *compiler)
{
    grib_action_if* a  = (grib_action_if*)act;
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

    fprintf(compiler->out,"%s = grib_action_create_if(ctx,",compiler->var);
    grib_compile_expression(a->expression, compiler);
    fprintf(compiler->out,",%s,%s,%d);\n", t,f,a->transient);
}

static int create_accessor( grib_section* p, grib_action* act, grib_loader *h)
{
  grib_action_if* a = (grib_action_if*)act;
  grib_action* next = NULL;
  int ret = 0;
  long lres=0;

  grib_accessor* as = NULL;
  grib_section*  gs = NULL;

  as = grib_accessor_factory(p, act,0,NULL);
  if(!as)return GRIB_INTERNAL_ERROR;
  gs = as->sub_section;
  grib_push_accessor(as,p->block);

  if ((ret=grib_expression_evaluate_long(p->h,a->expression,&lres)) != GRIB_SUCCESS)
    return ret;

  if(lres)
    next = a->block_true;
  else
    next = a->block_false;

#if 0
if(p->h->context->debug > 1)
{
  printf("EVALUATE create_accessor_handle ");
grib_expression_print(p->h->context,a->expression,p->h);
printf(" [%d]\n", next == a->block_true);

  grib_dump_action_branch(stdout,next,5);
}
#endif

  gs->branch = next;
  grib_dependency_observe_expression(as,a->expression);

  while(next){

    ret = grib_create_accessor(gs, next, h);
    if(ret != GRIB_SUCCESS) return ret;
    next= next->next;
  }

  return GRIB_SUCCESS;
}

static int execute(grib_action* act, grib_handle *h)
{
  grib_action_if* a = (grib_action_if*)act;
  grib_action* next = NULL;
  int ret = 0;
  long lres=0;

  if ((ret=grib_expression_evaluate_long(h,a->expression,&lres)) != GRIB_SUCCESS) {
    if (ret == GRIB_NOT_FOUND) lres=0;
    else
     return ret;
  }

  if(lres)
    next = a->block_true;
  else
    next = a->block_false;

  while(next){

    ret = grib_action_execute(next, h);
    if(ret != GRIB_SUCCESS) return ret;
    next= next->next;
  }

  return GRIB_SUCCESS;
}

static void dump(grib_action* act, FILE* f, int lvl)
{
  grib_action_if* a = (grib_action_if*)act;
  int i = 0;

  for (i=0;i<lvl;i++)
    grib_context_print(act->context,f,"     ");

  printf("if(%s) { ",act->name);  grib_expression_print(act->context,a->expression,0);
    printf("\n");

  if(a->block_true){
    /*      grib_context_print(act->context,f,"IF \t TODO \n");  TODO */
    grib_dump_action_branch(f,a->block_true,lvl+1);
  }
  if(a->block_false){
    printf("}\n");
  for (i=0;i<lvl;i++)
    grib_context_print(act->context,f,"     ");
  printf("else(%s) { ",act->name);  grib_expression_print(act->context,a->expression,0);
    /*     grib_context_print(act->context,f,"ELSE \n" );*/
    grib_dump_action_branch(f,a->block_false,lvl+1);
  }
  for (i=0;i<lvl;i++)
    grib_context_print(act->context,f,"     ");
    printf("}\n");
}


static grib_action* reparse(grib_action* a,grib_accessor* acc,int* doit)
{
  int ret=0;
  long lres=0;
  grib_action_if* self = (grib_action_if*)a;

  /* printf("reparse %s %s\n",a->name,acc->name); */

  if((ret=grib_expression_evaluate_long(acc->parent->h,self->expression,&lres)) != GRIB_SUCCESS)
    grib_context_log(acc->parent->h->context,
      GRIB_LOG_ERROR,"if reparse  grib_expression_evaluate_long %s",
        grib_get_error_message(ret));

  if(lres)
    return self->block_true;
  else
    return self->block_false;

}

static void destroy(grib_context* context,grib_action* act)
{
  grib_action_if* a = (grib_action_if*) act;
  grib_action *t = a->block_true;
  grib_action *f = a->block_false;

  while(t)
  {
    grib_action *nt = t->next;
    grib_free_action(context,t);
    t = nt;
  }

  while(f)
  {
    grib_action *nf = f->next;
    grib_free_action(context,f);
    f = nf;
  }


  grib_expression_free(context,a->expression);

  grib_context_free_persistent(context, act->name);
  grib_context_free_persistent(context, act->op);
}

static void xref(grib_action* d, FILE* f,const char* path)
{
}


