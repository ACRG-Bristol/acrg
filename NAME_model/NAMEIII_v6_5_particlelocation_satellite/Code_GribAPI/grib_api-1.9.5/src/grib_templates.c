/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

#include "grib_api_internal.h"

#if HAVE_UNISTD_H
# include <unistd.h>
#endif

typedef struct grib_templates {
  const char*           name;
  const unsigned char* data;
  size_t               size;
} grib_templates;

#if 1
#include "grib_templates.h"

#define NUMBER(x) (sizeof(x)/sizeof(x[0]))

grib_handle* grib_internal_template(grib_context* c,const char* name)
{
  int i;
  for(i = 0; i < NUMBER(templates); i++)
    if(strcmp(name,templates[i].name) == 0)
      return grib_handle_new_from_message_copy(c,templates[i].data,templates[i].size);
  return NULL;
}
#else
grib_handle* grib_internal_template(grib_context* c,const char* name)
{
      return NULL;
}
#endif

static grib_handle* try_template(grib_context* c,const char* dir,const char* name)
{
  char path[1024];
  grib_handle *g = NULL;
  int err = 0;

  sprintf(path,"%s/%s.tmpl",dir,name);

  if(access(path,F_OK) == 0)
  {
    FILE* f = fopen(path,"r");
    if(!f)
    {
      grib_context_log(c,GRIB_LOG_PERROR,"cannot open %s",path);
      return NULL;
    }
    g = grib_handle_new_from_file(c,f,&err);
    fclose(f);
  }

  return g;
}

static char* try_template_path(grib_context* c,const char* dir,const char* name)
{
  char path[1024];

  sprintf(path,"%s/%s.tmpl",dir,name);

  if(access(path,R_OK) == 0)
  {
    return grib_context_strdup(c,path);
  }

  return NULL;
}

grib_handle* grib_external_template(grib_context* c,const char* name)
{
  const char *base = c->grib_templates_path;
  char buffer[1024];
  char *p = buffer;
  grib_handle *g = NULL;
  /* printf("GRIB_TEMPLATES_PATH=%s\n",base); */

  if(!base) return NULL;

  while(*base)
  {
    if(*base == ':')
    {
      *p = 0;
      g = try_template(c,buffer,name);
      if(g) return g;
      p = buffer;
    }
    *p++ = *base++;
  }

  *p = 0;
  return g = try_template(c,buffer,name);
}

char* grib_external_template_path(grib_context* c,const char* name)
{
  const char *base = c->grib_templates_path;
  char buffer[1024];
  char *p = buffer;
  char *g = NULL;

  if(!base) return NULL;

  while(*base)
  {
    if(*base == ':')
    {
      *p = 0;
      g = try_template_path(c,buffer,name);
      if(g) return g;
      p = buffer;
      base++;
    }
    *p++ = *base++;
  }

  *p = 0;
  return g = try_template_path(c,buffer,name);
}
