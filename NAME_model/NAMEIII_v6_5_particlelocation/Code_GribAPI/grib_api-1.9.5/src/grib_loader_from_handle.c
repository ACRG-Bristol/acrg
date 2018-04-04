/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/**************************************
 *  Enrico Fucile
 **************************************/


#include "grib_api_internal.h"

static int copy_values(grib_handle* h,grib_accessor* ga)
{
  int i,j,k;
  /* printf("copy_values stack is %ld\n",(long)h->values_stack);*/
  for(j = 0; j < h->values_stack; j++)
  {
    for(i = 0; i < h->values_count[j]; i++)
    {
      for(k = 0; (k < MAX_ACCESSOR_NAMES) && (ga->all_names[k] != NULL); k++)
      {
        /*printf("copy_values: %s %s\n",h->values[j][i].name,ga->all_names[k]);*/
        if(strcmp(h->values[j][i].name,ga->all_names[k]) == 0)
        {
          size_t len = 1;
          /*printf("SET VALUES %s\n",h->values[j][i].name);*/
          switch(h->values[j][i].type)
          {
            case GRIB_TYPE_LONG:
              return grib_pack_long(ga,&h->values[j][i].long_value,&len);
              break;

            case GRIB_TYPE_DOUBLE:
              return grib_pack_double(ga,&h->values[j][i].double_value,&len);
              break;

            case GRIB_TYPE_STRING:
              len = strlen(h->values[j][i].string_value);
              return grib_pack_string(ga,h->values[j][i].string_value,&len);
              break;
          }

        }
      }

    }
  }

  return GRIB_NOT_FOUND;
}


int grib_lookup_long_from_handle(grib_context* gc,grib_loader* loader,const char* name,long* value)
{
  grib_handle* h = (grib_handle*)loader->data;
  grib_accessor *b = grib_find_accessor(h,name);
  size_t len = 1;
  if(b) return grib_unpack_long(b,value,&len);

  /* TODO: fix me. For now, we don't fail on a lookup. */
#if 1
  *value = -1;
  return GRIB_SUCCESS;
#else
  return GRIB_NOT_FOUND;
#endif
}

int grib_init_accessor_from_handle(grib_loader* loader,grib_accessor* ga,grib_arguments* default_value)
{
  grib_handle* h = (grib_handle*)loader->data;
  int ret = GRIB_SUCCESS;
  size_t len = 0;
  char*   sval = NULL;
  unsigned char*   uval = NULL;
  long*   lval = NULL;
  double* dval = NULL;
  static int first  = 1;
  static const char* missing = 0;
  const char* name = NULL;
  int k = 0;
  grib_handle *g;
  int e;
  grib_context_log(h->context,GRIB_LOG_DEBUG, "XXXXX Copying  %s",   ga->name);

  if(default_value)
  {
    grib_context_log(h->context,GRIB_LOG_DEBUG, "Copying:  setting %s to default value",
        ga->name);
    grib_pack_expression(ga,grib_arguments_get_expression(h,default_value,0));
  }

  if( (ga->flags & GRIB_ACCESSOR_FLAG_NO_COPY) ||
       ((ga->flags & GRIB_ACCESSOR_FLAG_EDITION_SPECIFIC) &&
        loader->changing_edition ) ||
      (ga->flags & GRIB_ACCESSOR_FLAG_FUNCTION) ||
     ( (ga->flags & GRIB_ACCESSOR_FLAG_READ_ONLY) &&
      !(ga->flags & GRIB_ACCESSOR_FLAG_COPY_OK) ) )
  {
    grib_context_log(h->context,GRIB_LOG_DEBUG, "Copying  %s  ignored",   ga->name);
    return GRIB_SUCCESS;
  }

#if 0
  if(h->values)
    if(copy_values(h,ga) == GRIB_SUCCESS)
    {
      grib_context_log(h->context,GRIB_LOG_DEBUG, "Copying: setting %s to multi-set-value",   ga->name);
      return GRIB_SUCCESS;
    }
#endif

#if 0
  if(h->loader)
    h->loader->init_accessor(h->loader,ga,default_value);
#else
  /* COMEBACK here
     this is needed if we reparse during reparse....
   */

  g = h;
  while(g)
  {
    if(g->values) {
      if(copy_values(g,ga) == GRIB_SUCCESS) {
        grib_context_log(h->context,GRIB_LOG_DEBUG, "Copying: setting %s to multi-set-value",   ga->name);
        return GRIB_SUCCESS;
      }
    }
    g = g->main;
  }
#endif




  /* Check if the same name exists in the original message ... */
  k = 0;
  while( (k < MAX_ACCESSOR_NAMES)            &&
         ((name = ga->all_names[k]) != NULL) &&
       ((ret = grib_get_size(h,name,&len)) != GRIB_SUCCESS)) k++;

  if(ret != GRIB_SUCCESS)
  {
    name = ga->name;

    if(first)
    {
      missing = getenv("GRIB_PRINT_MISSING");
      first = 0;
    }

    grib_context_log(h->context,GRIB_LOG_DEBUG, "Copying [%s] failed: %s",
        name,grib_get_error_message(ret));

    if(missing)
    {
      fprintf(stdout,"REPARSE: no value for %s",name);
      if(default_value)
        fprintf(stdout," (default value)");
      fprintf(stdout,"\n");
    }

    return GRIB_SUCCESS;
  }

  if(len == 0)
  {
    grib_context_log(h->context,GRIB_LOG_DEBUG,
                     "Copying  %s  failed, length is 0",   name);
    return GRIB_SUCCESS;
  }

  if((ga->flags & GRIB_ACCESSOR_FLAG_CAN_BE_MISSING) && grib_is_missing(h,name,&e) && e == GRIB_SUCCESS && len == 1)
    grib_pack_missing(ga);


  switch(grib_accessor_get_native_type(ga))
  {
    case GRIB_TYPE_STRING:

      len = len > 1024 ? len : 1024;
      sval = grib_context_malloc(h->context,len*sizeof(char));

      if(grib_get_string_internal(h,name,sval,&len) == GRIB_SUCCESS)
      {
        grib_context_log(h->context,GRIB_LOG_DEBUG, "Copying string %s to %s",  sval, name);
        ret = grib_pack_string(ga,sval,&len);
      }
      grib_context_free(h->context,sval);

      break;

    case GRIB_TYPE_LONG:
      lval = grib_context_malloc(h->context,len*sizeof(long));

      if(grib_get_long_array_internal(h,name,lval,&len) == GRIB_SUCCESS)
      {
        grib_context_log(h->context,GRIB_LOG_DEBUG, "Copying %d long(s) %d to %s",  len, lval[0], name);
    if(ga->same)
    {
      ret = grib_set_long_array(ga->parent->h,ga->name,lval,len);

      /* Allow for lists to be resized */
      if((ret == GRIB_WRONG_ARRAY_SIZE || ret == GRIB_ARRAY_TOO_SMALL) && loader->list_is_resized)
        ret = GRIB_SUCCESS;
    }
    else
      ret = grib_pack_long(ga,lval,&len);
      }

      grib_context_free(h->context,lval);

      break;

    case GRIB_TYPE_DOUBLE:
      dval = grib_context_malloc(h->context,len*sizeof(double));
      if(grib_get_double_array_internal(h,name,dval,&len) == GRIB_SUCCESS)
      {
        grib_context_log(h->context,GRIB_LOG_DEBUG, "Copying %d double(s) %g to %s",  len, dval[0], name);
        if(ga->same)
        {
          ret = grib_set_double_array(ga->parent->h,ga->name,dval,len);

          /* Allow for lists to be resized */
          if((ret == GRIB_WRONG_ARRAY_SIZE || ret == GRIB_ARRAY_TOO_SMALL) && loader->list_is_resized)
            ret = GRIB_SUCCESS;
        }
        else ret = grib_pack_double(ga,dval,&len);
      }


      grib_context_free(h->context,dval);
      break;

    case GRIB_TYPE_BYTES:

      uval = grib_context_malloc(h->context,len*sizeof(char));

      if(grib_get_bytes_internal(h,name,uval,&len) == GRIB_SUCCESS)
      {
        grib_context_log(h->context,GRIB_LOG_DEBUG, "Copying %d byte(s) to %s",  len, name);
        ret = grib_pack_bytes(ga,uval,&len);
      }


      grib_context_free(h->context,uval);

      break;

    case GRIB_TYPE_LABEL:
      break;

    default:
      grib_context_log(h->context,GRIB_LOG_ERROR, "Copying %s, cannot establish type %d [%s]"
          , name,grib_accessor_get_native_type(ga),ga->creator->cclass->name);
      break;
  }

  return ret;

}
