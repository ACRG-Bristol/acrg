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

static void init(grib_action_class *c)
{
    if(c && !c->inited)
    {
        init(c->super ? *(c->super) : NULL);
        c->init_class(c);
        c->inited = 1;
    }
}

void grib_dump(grib_action* a, FILE* f, int l)
{
    grib_action_class *c = a->cclass;
    init(c);

    while(c)
    {
        if(c->dump)
        {
            c->dump(a, f, l);
            return;
        }
        c = c->super ? *(c->super) : NULL;
    }
    Assert(0);
}

void grib_xref(grib_action* a, FILE* f,const char* path)
{
    grib_action_class *c = a->cclass;
    init(c);

    while(c)
    {
        if(c->xref)
        {
            c->xref(a, f,path);
            return;
        }
        c = c->super ? *(c->super) : NULL;
    }
    printf("xref not implemented for %s\n",a->cclass->name);
    Assert(0);
}


void grib_free_action(grib_context* context,grib_action* a)
{
    grib_action_class *c = a->cclass;
    init(c);
    while(c)
    {
        if(c->destroy)
            c->destroy(context, a);
        c = c->super ? *(c->super) : NULL;
    }
    grib_context_free_persistent(context, a);
}

int grib_create_accessor(grib_section* p, grib_action* a,  grib_loader* h)
{
    grib_action_class *c = a->cclass;
    init(c);
    while(c)
    {
        if(c->create_accessor)
            return c->create_accessor(p, a, h);
        c = c->super ? *(c->super) : NULL;
    }
    fprintf(stderr,"Cannot create accessor %s %s\n",a->name,a->cclass->name);
    Assert(0);
    return 0;
}

int grib_action_notify_change( grib_action* a, grib_accessor *observer, grib_accessor *observed)
{
    grib_action_class *c = a->cclass;
    init(c);
    while(c)
    {
        if(c->notify_change)
            return c->notify_change(a,observer,observed);
        c = c->super ? *(c->super) : NULL;
    }
    Assert(0);
    return 0;
}

grib_action* grib_action_reparse( grib_action* a, grib_accessor* acc,int* doit)
{
    grib_action_class *c = a->cclass;
    init(c);
    while(c)
    {
        if(c->reparse)
            return c->reparse(a,acc,doit);
        c = c->super ? *(c->super) : NULL;
    }
    Assert(0);
    return 0;
}

int grib_action_execute( grib_action* a, grib_handle* h)
{
    grib_action_class *c = a->cclass;
    init(c);
    while(c)
    {
        if(c->execute)
            return c->execute(a,h);
        c = c->super ? *(c->super) : NULL;
    }
    Assert(0);
    return 0;
}

void grib_dump_action_branch(FILE* out,grib_action* a, int decay)
{
    while(a)
    {
        grib_dump(a,out,decay);
        a=a->next;
    }
}

void grib_dump_action_tree( grib_context* ctx,FILE* out)
{
    grib_dump_action_branch( out, ctx->grib_reader->first->root ,0);
}

void grib_xref_action_branch(FILE* out,grib_action* a,const char* path)
{
    while(a)
    {
        grib_xref(a,out,path);
        a=a->next;
    }
}

void grib_compile(grib_action* a, grib_compiler* compiler)
{
    grib_action_class *c = a->cclass;
    init(c);
    if(c->compile) {
        c->compile(a,compiler);
    }
    else 
    {
        fprintf(stderr, "NO COMPILE METHOD '%s'\n", c->name);
        Assert(0);
    }
}
