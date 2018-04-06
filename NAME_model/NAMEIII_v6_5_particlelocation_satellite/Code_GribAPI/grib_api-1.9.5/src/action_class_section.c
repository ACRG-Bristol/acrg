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
   IMPLEMENTS = notify_change
   IMPLEMENTS = reparse
   IMPLEMENTS = xref
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "action.class" and rerun ./make_class.pl

*/

static void init_class      (grib_action_class*);
static void xref            (grib_action* d, FILE* f,const char* path);
static int notify_change(grib_action* a, grib_accessor* observer,grib_accessor* observed);
static grib_action* reparse(grib_action* a,grib_accessor* acc,int *doit);


typedef struct grib_action_section {
    grib_action          act;  
/* Members defined in section */
} grib_action_section;


static grib_action_class _grib_action_class_section = {
    0,                              /* super                     */
    "action_class_section",                              /* name                      */
    sizeof(grib_action_section),            /* size                      */
    0,                                   /* inited */
    &init_class,                         /* init_class */
    0,                               /* init                      */
    0,                            /* destroy */

    0,                               /* dump                      */
    &xref,                               /* xref                      */

    0,             /* create_accessor*/

    &notify_change,                            /* notify_change */
    &reparse,                            /* reparse */
    0,                            /* execute */
    0,                            /* compile */
};

grib_action_class* grib_action_class_section = &_grib_action_class_section;

static void init_class(grib_action_class* c)
{
}
/* END_CLASS_IMP */

static void check_sections(grib_section *s,grib_handle* h)
{
  grib_accessor *a = s?s->block->first:NULL;
  if(s) Assert(s->h == h);
  while(a)
  {
    Assert(a->parent->h == h);
    check_sections(a->sub_section,h);
    a = a->next;
  }
}


static int notify_change(grib_action* act, grib_accessor * notified,
                         grib_accessor* changed)
{

  grib_loader loader = { 0,};

  grib_section *old_section = NULL;
  grib_handle *h = notified->parent->h;
  size_t len = 0;
  size_t size = 0;
  int err=0;
  grib_handle* tmp_handle;
  int doit = 0;

  grib_action* la        = NULL;

  grib_context_log(h->context,
      GRIB_LOG_DEBUG,"------------- SECTION action %s (%s) is triggerred by [%s]",
      act->name, notified->name, changed->name);

  la = grib_action_reparse(act,notified,&doit);
  old_section = notified->sub_section;
  Assert(old_section);

  Assert(old_section->h == h);

  /* printf("old = %p\n",(void*)old_section->branch); */
  /* printf("new = %p\n",(void*)la); */

  grib_context_log(h->context,
      GRIB_LOG_DEBUG,"------------- DOIT %ld OLD %p NEW %p",
      doit,old_section->branch,la);


  if(!doit) {
    if(la != NULL || old_section->branch != NULL)
      if(la == old_section->branch)
      {
        grib_context_log(h->context,GRIB_LOG_DEBUG,"IGNORING TRIGGER action %s (%s) is triggerred %p", act->name, notified->name
            ,(void*)la);
        return GRIB_SUCCESS;
      }
  }

  loader.list_is_resized = (la == old_section->branch);

  if (!strcmp(changed->name,"GRIBEditionNumber")) loader.changing_edition=1;
  else loader.changing_edition=0;

  old_section->branch = la;

  tmp_handle = grib_new_handle(h->context);
  if(!tmp_handle)
    return GRIB_OUT_OF_MEMORY;

  tmp_handle->buffer = grib_create_growable_buffer(h->context);
  Assert(tmp_handle->buffer); /* FIXME */

  loader.data          = h;
  loader.lookup_long   = grib_lookup_long_from_handle;
  loader.init_accessor = grib_init_accessor_from_handle;


  Assert(h->kid == NULL);
  tmp_handle->loader = &loader;
  tmp_handle->main  = h;
  h->kid = tmp_handle;
  /* printf("tmp_handle- main %p %p\n",(void*)tmp_handle,(void*)h); */

  grib_context_log(h->context,GRIB_LOG_DEBUG,"------------- CREATE TMP BLOCK ", act->name, notified->name);
  tmp_handle->root  = grib_section_create(tmp_handle,NULL);

  tmp_handle->use_trie=1;

  err=grib_create_accessor(tmp_handle->root, act, &loader);

  grib_section_adjust_sizes(tmp_handle->root,1,0);

  grib_section_post_init(tmp_handle->root);

  /* grib_recompute_sections_lengths(tmp_handle->root); */
  grib_get_block_length(tmp_handle->root,&len);
  grib_context_log(h->context,GRIB_LOG_DEBUG,"-------------  TMP BLOCK IS sectlen=%d buffer=%d", len,  tmp_handle->buffer->ulength);

#if 0
  if(h->context->debug > 10)
    grib_dump_content(tmp_handle,stdout,NULL,0,NULL);
#endif

  Assert(tmp_handle->buffer->ulength == len);
  /* grib_empty_section(h->context,old_section); */

  grib_buffer_replace(notified, tmp_handle->buffer->data, tmp_handle->buffer->ulength,0,1);

  grib_swap_sections(old_section,
      tmp_handle->root->block->first->sub_section);

  Assert(tmp_handle->dependencies == NULL);
  /* printf("grib_handle_delete %p\n",(void*)tmp_handle); */


  grib_handle_delete(tmp_handle);

  h->use_trie = 1;
  h->trie_invalid=1;
  h->kid = NULL;

  grib_section_adjust_sizes(h->root,1,0);

  grib_section_post_init(h->root);

  grib_get_block_length(old_section,&size);

  grib_context_log(h->context,GRIB_LOG_DEBUG,"-------------   BLOCK SIZE %ld, buffer len=%ld", size,len);
  if(h->context->debug > 10)
    grib_dump_content(h,stdout,"debug",~0,NULL);

  Assert(size == len);

  grib_update_paddings(old_section);


  return err;
}


static grib_action* reparse(grib_action* a,grib_accessor* acc,int* doit)
{
  /* Should be inerited */
  printf("reparse should be inerited: %s\n",a->name);

  Assert(1==0);
  return 0;
}

static void xref(grib_action* d, FILE* f,const char* path)
{
}
