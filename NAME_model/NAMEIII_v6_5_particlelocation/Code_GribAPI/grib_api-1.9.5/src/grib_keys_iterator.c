/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: grib_keys_iterator
 *
 * Description:
 *
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */

#include "grib_api_internal.h"

GRIB_INLINE static int grib_inline_strcmp(const char* a,const char* b) {
  if (*a != *b) return 1;
  while((*a!=0 && *b!=0) &&  *(a) == *(b) ) {a++;b++;}
  return (*a==0 && *b==0) ? 0 : 1;
}

struct grib_keys_iterator{
  grib_handle   *handle;
  unsigned long filter_flags;     /** flags to filter out accessors */
  unsigned long accessor_flags;     /** flags to filter out accessors */
  grib_accessor *current;
  char    *name_space;
  int            at_start;
  int            match;
  grib_trie     *seen;
};


grib_keys_iterator*  grib_keys_iterator_new(grib_handle* h,unsigned long filter_flags, char* name_space) {

  grib_keys_iterator* ki=NULL;

  if (!h) return NULL;

  ki= grib_context_malloc_clear(h->context,sizeof(grib_keys_iterator));
  if (!ki) return NULL;

  ki->filter_flags = filter_flags;
  ki->handle       = h;
  ki->name_space   = NULL;

  if (name_space != NULL)
    ki->name_space   = grib_context_strdup(h->context,name_space);

  ki->at_start     = 1;
  ki->match        = 0;

  grib_keys_iterator_set_flags(ki,filter_flags);

  return ki;
}

int grib_keys_iterator_set_flags(grib_keys_iterator* ki,unsigned long flags) {
  int ret=0;
  grib_handle* h;

  if (!ki) return GRIB_INTERNAL_ERROR;

  h=ki->handle;

  if(flags & GRIB_KEYS_ITERATOR_SKIP_DUPLICATES && ki->seen!=NULL )
    ki->seen = grib_trie_new(h->context);

  if(flags & GRIB_KEYS_ITERATOR_SKIP_FUNCTION)
    ki->accessor_flags |= GRIB_ACCESSOR_FLAG_FUNCTION;

  if(flags & GRIB_KEYS_ITERATOR_SKIP_READ_ONLY)
    ki->accessor_flags |= GRIB_ACCESSOR_FLAG_READ_ONLY;

  if(flags & GRIB_KEYS_ITERATOR_SKIP_EDITION_SPECIFIC)
    ki->accessor_flags |= GRIB_ACCESSOR_FLAG_EDITION_SPECIFIC;

  return ret;
}

static void mark_seen(grib_keys_iterator* ki,const char* name)
{
  grib_trie_insert(ki->seen,name,(void*)name);
}

static int was_seen(grib_keys_iterator* ki,const char* name)
{
  return grib_trie_get(ki->seen,name) != NULL;
}


int grib_keys_iterator_rewind(grib_keys_iterator* ki)
{
  ki->at_start = 1;
  return GRIB_SUCCESS;
}

static int skip(grib_keys_iterator* kiter) {

  /* TODO: set the section to hidden, to speed up that */
  /* if(grib_get_sub_section(kiter->current)) */
  if(kiter->current->sub_section)
    return 1;

  if(kiter->current->flags & GRIB_ACCESSOR_FLAG_HIDDEN)
    return 1;

  if(kiter->current->flags &  kiter->accessor_flags)
    return 1;

  if((kiter->filter_flags & GRIB_KEYS_ITERATOR_SKIP_COMPUTED) && kiter->current->length == 0)
    return 1;

  if((kiter->filter_flags & GRIB_KEYS_ITERATOR_SKIP_CODED) && kiter->current->length != 0)
    return 1;

  if(kiter->name_space)
  {

    kiter->match = 0;

    while(kiter->match < MAX_ACCESSOR_NAMES)
    {
      if(kiter->current->all_name_spaces[kiter->match] != NULL)
        if(grib_inline_strcmp(kiter->current->all_name_spaces[kiter->match],kiter->name_space) == 0)
        {
          if(kiter->seen)
          {
            if(was_seen(kiter,kiter->current->all_names[kiter->match]))
              return 1;
            mark_seen(kiter,kiter->current->all_names[kiter->match]);
          }
          return 0;
        }

      kiter->match++;
    }

    return 1;

  }

  if(kiter->seen)
  {
    if(was_seen(kiter,kiter->current->name))
      return 1;
    mark_seen(kiter,kiter->current->name);
  }

  return 0;

}

int grib_keys_iterator_next(grib_keys_iterator* kiter){
  if(kiter->at_start) {
    kiter->current  = kiter->handle->root->block->first;
    kiter->at_start = 0;
  }
  else
  {
    kiter->current = grib_next_accessor(kiter->current);
  }

  while(kiter->current && skip(kiter))
    kiter->current = grib_next_accessor(kiter->current);

  return kiter->current != NULL;
}


const char* grib_keys_iterator_get_name(grib_keys_iterator* kiter)
{
  /* if(kiter->name_space) */
  return kiter->current->all_names[kiter->match];
}

grib_accessor* grib_keys_iterator_get_accessor(grib_keys_iterator* kiter)
{
  return kiter->current;
}

int grib_keys_iterator_delete( grib_keys_iterator* kiter){
  if (kiter) {
    if(kiter->seen)
      grib_trie_delete(kiter->seen);
    if (kiter->name_space)
     grib_context_free(kiter->handle->context,kiter->name_space);
    grib_context_free(kiter->handle->context,kiter);
  }
  return 0;
}

int grib_keys_iterator_get_long(grib_keys_iterator* kiter,long* v,size_t* len) {
  return grib_unpack_long( kiter->current,v,len);
}

int grib_keys_iterator_get_double(grib_keys_iterator* kiter,double* v,size_t* len) {
  return grib_unpack_double( kiter->current,v,len);
}

int grib_keys_iterator_get_string(grib_keys_iterator* kiter,char* v,size_t* len) {
  return grib_unpack_string( kiter->current,v,len);
}

int grib_keys_iterator_get_bytes(grib_keys_iterator* kiter,unsigned char* v,size_t* len) {
  return grib_unpack_bytes( kiter->current,v,len);
}

int grib_keys_iterator_get_native_type(grib_keys_iterator* kiter) {
  return grib_accessor_get_native_type(kiter->current);
}

