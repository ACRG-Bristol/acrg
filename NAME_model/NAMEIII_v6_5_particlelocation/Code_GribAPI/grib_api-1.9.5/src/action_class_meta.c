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
 *                                                                         *
 ***************************************************************************/
#include "grib_api_internal.h"

/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = action
   SUPER      = action_class_gen
   IMPLEMENTS = dump
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


typedef struct grib_action_meta {
    grib_action          act;  
/* Members defined in gen */
	long            len;
	grib_arguments* params;
/* Members defined in meta */
} grib_action_meta;

extern grib_action_class* grib_action_class_gen;

static grib_action_class _grib_action_class_meta = {
    &grib_action_class_gen,                              /* super                     */
    "action_class_meta",                              /* name                      */
    sizeof(grib_action_meta),            /* size                      */
    0,                                   /* inited */
    &init_class,                         /* init_class */
    0,                               /* init                      */
    0,                            /* destroy */

    &dump,                               /* dump                      */
    0,                               /* xref                      */

    0,             /* create_accessor*/

    0,                            /* notify_change */
    0,                            /* reparse */
    0,                            /* execute */
    0,                            /* compile */
};

grib_action_class* grib_action_class_meta = &_grib_action_class_meta;

static void init_class(grib_action_class* c)
{
	c->xref	=	(*(c->super))->xref;
	c->create_accessor	=	(*(c->super))->create_accessor;
	c->notify_change	=	(*(c->super))->notify_change;
	c->reparse	=	(*(c->super))->reparse;
	c->execute	=	(*(c->super))->execute;
	c->compile	=	(*(c->super))->compile;
}
/* END_CLASS_IMP */

grib_action* grib_action_create_meta( grib_context* context, const char* name, const char* op, 
		grib_arguments*   params,  grib_arguments*   default_value,unsigned long flags,const char* name_space)
{
	grib_action_meta*  a   =  (grib_action_meta*)grib_context_malloc_clear_persistent(context,sizeof(grib_action_meta));
	grib_action* act       =  (grib_action*)a;
	act->next              =  NULL;
	act->name              =  grib_context_strdup_persistent(context, name);
	act->op                =  grib_context_strdup_persistent(context, op);
	if(name_space)
	act->name_space                =  grib_context_strdup_persistent(context, name_space);
	act->cclass            =  grib_action_class_meta;
	act->context           = context;
	act->flags             = flags;
	a->params              =  params;
	act->default_value                = default_value;
	a->len                 = 0;
	

	/* grib_arguments_print(context,a->params,0); printf("\n"); */


	return act;
}


static void dump( grib_action* act, FILE* f, int lvl)
{

	int i =0;
	for (i=0;i<lvl;i++)
		grib_context_print(act->context,f,"     ");
	grib_context_print(act->context,f," meta %s \n", act->name );
}

