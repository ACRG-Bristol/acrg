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

void grib_get_buffer_ownership(const grib_context *c, grib_buffer *b)
{
  unsigned char* newdata;
  if(b->property == GRIB_MY_BUFFER)
    return;

  newdata = grib_context_malloc(c, b->length);
  memcpy(newdata, b->data, b->length);
  b->data = newdata;
  b->property = GRIB_MY_BUFFER;
}

grib_buffer* grib_create_growable_buffer(const grib_context* c)
{
  grib_buffer  *b =  (grib_buffer*)grib_context_malloc_clear(c,sizeof(grib_buffer));

  if(b == NULL)
  {
    grib_context_log(c,GRIB_LOG_ERROR,"grib_new_buffer: cannot allocate buffer");
    return NULL;
  }

  b->property = GRIB_MY_BUFFER;
  b->length   = 10240;
  b->ulength  = 0;
  b->data     = grib_context_malloc_clear(c,b->length);
  b->growable = 1;

  if(!b->data)
  {
    grib_context_log(c,GRIB_LOG_ERROR,"grib_new_buffer: cannot allocate buffer");
    grib_context_free(c,b);
    return NULL;
  }

  return b;

}


grib_buffer* grib_new_buffer(const grib_context* c,unsigned char* data,size_t buflen)
{
  grib_buffer  *b =  (grib_buffer*)grib_context_malloc_clear(c,sizeof(grib_buffer));

  if(b == NULL)
  {
    grib_context_log(c,GRIB_LOG_ERROR,"grib_new_buffer: cannot allocate buffer");
    return NULL;
  }


  b->property = GRIB_USER_BUFFER;
  b->length   = buflen;
  b->ulength  = buflen;
  b->data     = data;

  return b;
}

void grib_buffer_delete(const grib_context *c, grib_buffer *b)
{

  if(b->property == GRIB_MY_BUFFER)
    grib_context_free(c,b->data);
  b->length = 0;
  b->ulength = 0;
  grib_context_free(c,b);

}

static void grib_grow_buffer_to(const grib_context *c, grib_buffer *b, size_t ns)
{
  unsigned char* newdata;

  if(ns>b->length)
  {
    grib_get_buffer_ownership(c, b);
    newdata = grib_context_malloc_clear(c, ns);
    memcpy(newdata, b->data, b->length);
    grib_context_free(c,b->data);
    b->data = newdata;
    b->length = ns;
  }
}

void grib_grow_buffer(const grib_context *c, grib_buffer *b, size_t new_size)
{
  size_t len = ((new_size + 1023)/1024)*1024;
  grib_grow_buffer_to(c,b,len);
}

void grib_buffer_set_ulength(const grib_context *c, grib_buffer *b, size_t length)
{
  grib_grow_buffer_to(c,b,length);
  b->ulength = length;
}


static void update_offsets(grib_accessor* a,long len)
{
  while(a)
  {
    grib_section* s = a->sub_section;
    a->offset += len;
    grib_context_log(a->parent->h->context,GRIB_LOG_DEBUG, "::::: grib_buffer : accessor %s is moving by %d bytes to %ld",a->name ,len, a->offset);
    if(s) update_offsets(s->block->first,len);
    a = a->next;
  }
}


static void update_offsets_after(grib_accessor* a,long len)
{
  while(a)
  {
    update_offsets(a->next,len);
    a = a->parent->owner;
  }
}

void grib_recompute_sections_lengths(grib_section* s)
{
  if(s)
  {
    long   plen = 0;
    size_t  len = 1;

    grib_accessor* a = s->block->first;

    while(a)
    {
      /* grib_recompute_sections_lengths(grib_get_sub_section(a)); */
      grib_recompute_sections_lengths(a->sub_section);
      a = a->next;
    }

    if(s->aclength)
    {
      int ret;
      if(s->owner)
        plen = grib_get_next_position_offset(s->block->last) - s->owner->offset;
      else
        plen = grib_get_next_position_offset(s->block->last);

      if((ret = grib_pack_long(s->aclength, &plen, &len)) != GRIB_SUCCESS)
        ;

#if 0
      if(s->h->context->debug)
        printf("SECTION updating length %ld .. %s\n",plen,s->owner->name);
#endif

    }

  }

}

static void update_sections_lengths(grib_section* s)
{

  long   plen = 0;
  size_t  len = 1;

  if(!s) return;


  if(s->aclength)
  {
    int ret;
    if(s->owner)
      plen = grib_get_next_position_offset(s->block->last) - s->owner->offset;
    else
      plen = grib_get_next_position_offset(s->block->last);

    /* if(s->owner) */
      /* s->owner->length = plen; */

  /* if(s->aclength)  */
    if((ret = grib_pack_long(s->aclength, &plen, &len)) != GRIB_SUCCESS)
      ;

     if(s->h->context->debug)
     {
      printf("SECTION updating length %ld .. %s\n",plen,s->owner->name);
      printf("NEXT_POS = %ld, owner offset= %ld %s %s\n",
          grib_get_next_position_offset(s->block->last),
          s->owner ? s->owner->offset : 0L, s->owner->name,
          s->block->last->name);

     }

  }

  if(s->owner)
    update_sections_lengths(s->owner->parent);


}

void grib_buffer_replace( grib_accessor *a, const unsigned char* data,
                          size_t newsize,int update_lengths,int update_paddings)
{
  size_t offset   = a->offset;
  long   oldsize  = grib_get_next_position_offset(a)-offset;
  long   increase = (long)newsize - (long)oldsize;

  grib_buffer *buffer     = a->parent->h->buffer;
  size_t message_length   = buffer->ulength;

  grib_context_log(a->parent->h->context,GRIB_LOG_DEBUG,
     "grib_buffer_replace %s offset=%ld oldsize=%ld newsize=%ld message_length=%ld\n",
      a->name,(long)offset,oldsize,(long)newsize,(long)message_length);

    grib_buffer_set_ulength(a->parent->h->context,
        buffer,
        buffer->ulength+increase);

  /* move the end */
  if(increase)
    memmove(
      buffer->data + offset + newsize,
      buffer->data + offset + oldsize,
      message_length - offset - oldsize);

  /* copy new data */

  memcpy(buffer->data + offset, data, newsize);

  /* if(increase) */
  if(increase)
  {
    update_offsets_after(a,increase);
    if(update_lengths)
    {
      grib_update_size(a,newsize);
      grib_section_adjust_sizes(a->parent->h->root,1,0);
      if(update_paddings)
        grib_update_paddings(a->parent->h->root);
    }

  }


}




void grib_update_sections_lengths(grib_handle* h) {
  grib_section_adjust_sizes(h->root,2,0);
    grib_update_paddings(h->root);
}




