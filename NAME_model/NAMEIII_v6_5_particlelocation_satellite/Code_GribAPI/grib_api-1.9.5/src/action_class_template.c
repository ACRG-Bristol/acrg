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
   IMPLEMENTS = reparse
   IMPLEMENTS = compile
   MEMBERS    = int nofail
   MEMBERS    = char*           arg
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
static void compile         (grib_action* a, grib_compiler* compiler);
static void destroy         (grib_context*,grib_action*);
static int create_accessor(grib_section*,grib_action*,grib_loader*);
static grib_action* reparse(grib_action* a,grib_accessor* acc,int *doit);


typedef struct grib_action_template {
    grib_action          act;  
/* Members defined in section */
/* Members defined in template */
	int nofail;
	char*           arg;
} grib_action_template;

extern grib_action_class* grib_action_class_section;

static grib_action_class _grib_action_class_template = {
    &grib_action_class_section,                              /* super                     */
    "action_class_template",                              /* name                      */
    sizeof(grib_action_template),            /* size                      */
    0,                                   /* inited */
    &init_class,                         /* init_class */
    0,                               /* init                      */
    &destroy,                            /* destroy */

    &dump,                               /* dump                      */
    0,                               /* xref                      */

    &create_accessor,             /* create_accessor*/

    0,                            /* notify_change */
    &reparse,                            /* reparse */
    0,                            /* execute */
    &compile,                            /* compile */
};

grib_action_class* grib_action_class_template = &_grib_action_class_template;

static void init_class(grib_action_class* c)
{
	c->xref	=	(*(c->super))->xref;
	c->notify_change	=	(*(c->super))->notify_change;
	c->execute	=	(*(c->super))->execute;
}
/* END_CLASS_IMP */

grib_action* grib_action_create_template( grib_context* context,int nofail,const char* name , const char* arg1)
{
  grib_action_template* a   ;
  grib_action_class* c   =  grib_action_class_template;
  grib_action* act       =  (grib_action*)grib_context_malloc_clear_persistent(context,c->size);
  act-> name             =  grib_context_strdup_persistent(context,name);
  act-> op               =  grib_context_strdup_persistent(context,"section");
  act-> cclass           =  c;
  act-> next             =  NULL;
  act->context           =  context;
  a                      =  (grib_action_template*)act;
  a->nofail=nofail;
  if (arg1) a->arg       =  grib_context_strdup_persistent(context,arg1);
  else  a->arg = NULL;

  return act;
}

static void compile(grib_action* act, grib_compiler *compiler)
{
    grib_action_template* a  = (grib_action_template*)act;
    fprintf(compiler->out,"%s = grib_action_create_template(ctx,", compiler->var);
    fprintf(compiler->out,"%d,",a->nofail);
    fprintf(compiler->out,"\"%s\",",act->name);
    if(a->arg) {
        fprintf(compiler->out,"\"%s\");",a->arg);
    }
    else
    {
        fprintf(compiler->out,"NULL);");
    }    
    fprintf(compiler->out,"\n");
}    

static void dump( grib_action* act, FILE* f, int lvl)
{
  grib_action_template* a = ( grib_action_template*)act;
  int i = 0;
  for (i=0;i<lvl;i++) grib_context_print(act->context,f,"     ");
  grib_context_print(act->context,f,"Template %s  %s\n",act->name , a->arg );
}

GRIB_INLINE grib_action* get_empty_template(grib_context* c,int *err) {
	char fname[]="empty_template.def";
	char* path=0;
	
	path=grib_context_full_path(c,fname);
	if (path) {
		*err=GRIB_SUCCESS;
		return grib_parse_file(c, path);
	} else {
		*err=GRIB_INTERNAL_ERROR;
		grib_context_log(c,GRIB_LOG_ERROR,"get_empty_template: unable to get template %s",fname);
		return NULL;
	}
	
}

static int  create_accessor(grib_section* p, grib_action* act, grib_loader *h ){
  int ret = GRIB_SUCCESS;
  grib_action_template* a = ( grib_action_template*)act;
  grib_action* la = NULL;
  grib_action* next = NULL;
  grib_accessor* as = NULL;
  grib_section*         gs = NULL;

  char fname[1024]={0,};
  char *fpath=0;

  as = grib_accessor_factory(p, act,0,NULL);

  if(!as) return GRIB_INTERNAL_ERROR;
  if(a->arg){
    ret = grib_recompose_name(p->h,as,a->arg,fname,1);

	if ((fpath=grib_context_full_path(p->h->context,fname))==NULL) {
      if (!a->nofail) {
        grib_context_log(p->h->context,GRIB_LOG_ERROR,
                         "Unable to find template %s from %s ",act->name,fname);
        return GRIB_FILE_NOT_FOUND;
      }
	  la = get_empty_template(p->h->context,&ret);
	  if (ret) return ret;
    } else 
      la = grib_parse_file(p->h->context, fpath);
  }
  as->flags |= GRIB_ACCESSOR_FLAG_HIDDEN;
  gs = as->sub_section;
  gs->branch = la; /* Will be used to prevent unecessary reparse */

  grib_push_accessor(as,p->block);

  if(la){
    next = la;

    while(next){
      ret = grib_create_accessor(gs, next,h);
      if(ret != GRIB_SUCCESS) {
      if(p->h->context->debug)
    {
      grib_context_log(p->h->context,GRIB_LOG_ERROR,
      "Error processing template %s: %s [%s] %04lx",
      fname,grib_get_error_message(ret),next->name,next->flags);
    }
      return ret;
    }
      next= next->next;
    }
  }
  return GRIB_SUCCESS;
}


static grib_action* reparse(grib_action* a,grib_accessor* acc,int *doit)
{
  grib_action_template* self = (grib_action_template*)a;
  char *fpath=0;

  if(self->arg){
    char fname[1024];
    grib_recompose_name(acc->parent->h,NULL,self->arg,fname,1);

    if ((fpath=grib_context_full_path(acc->parent->h->context,fname))==NULL) {
      if (!self->nofail) {
        grib_context_log(acc->parent->h->context,GRIB_LOG_ERROR,
                         "Unable to find template %s from %s ",a->name,fname);
        return NULL;
      } return a;
    }

    /* printf("REPARSE %s\n",fpath); */
    return grib_parse_file(acc->parent->h->context, fpath);
  }

  return NULL;

}


static void destroy(grib_context* context,grib_action* act)
{
  grib_action_template* a = (grib_action_template*)act;

  grib_context_free_persistent(context, a->arg);
  grib_context_free_persistent(context, act->name);
  grib_context_free_persistent(context, act->op);
}
