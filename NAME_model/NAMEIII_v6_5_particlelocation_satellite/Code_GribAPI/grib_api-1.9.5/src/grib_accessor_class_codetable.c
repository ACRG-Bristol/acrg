/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*****************************************
 *  Enrico Fucile
 ****************************************/

#include "grib_api_internal.h"
#include <ctype.h>

/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = accessor
   SUPER      = grib_accessor_class_unsigned
   IMPLEMENTS = init;dump;unpack_string;pack_expression;unpack_long
   IMPLEMENTS = value_count;pack_string; destroy; get_native_type;
   MEMBERS    =  const char* tablename
   MEMBERS    =  const char* masterDir
   MEMBERS    =  const char* localDir
   MEMBERS    =  grib_codetable* table
   END_CLASS_DEF

 */

/* START_CLASS_IMP */

/*

Don't edit anything between START_CLASS_IMP and END_CLASS_IMP
Instead edit values between START_CLASS_DEF and END_CLASS_DEF
or edit "accessor.class" and rerun ./make_class.pl

*/

static int  get_native_type(grib_accessor*);
static int pack_string(grib_accessor*, const char*, size_t *len);
static int pack_expression(grib_accessor*, grib_expression*);
static int unpack_long(grib_accessor*, long* val,size_t *len);
static int unpack_string (grib_accessor*, char*, size_t *len);
static long value_count(grib_accessor*);
static void destroy(grib_context*,grib_accessor*);
static void dump(grib_accessor*, grib_dumper*);
static void init(grib_accessor*,const long, grib_arguments* );
static void init_class(grib_accessor_class*);

typedef struct grib_accessor_codetable {
    grib_accessor          att;
/* Members defined in gen */
/* Members defined in long */
/* Members defined in unsigned */
	long nbytes;
	grib_arguments* arg;
/* Members defined in codetable */
	const char* tablename;
	const char* masterDir;
	const char* localDir;
	grib_codetable* table;
} grib_accessor_codetable;

extern grib_accessor_class* grib_accessor_class_unsigned;

static grib_accessor_class _grib_accessor_class_codetable = {
    &grib_accessor_class_unsigned,                      /* super                     */
    "codetable",                      /* name                      */
    sizeof(grib_accessor_codetable),  /* size                      */
    0,                           /* inited */
    &init_class,                 /* init_class */
    &init,                       /* init                      */
    0,                  /* post_init                      */
    &destroy,                    /* free mem                       */
    &dump,                       /* describes himself         */
    0,                /* get length of section     */
    &value_count,                /* get number of values      */
    0,                 /* get number of bytes      */
    0,                /* get offset to bytes           */
    &get_native_type,            /* get native type               */
    0,                /* get sub_section                */
    0,               /* grib_pack procedures long      */
    0,               /* grib_pack procedures long      */
    0,                  /* grib_pack procedures long      */
    &unpack_long,                /* grib_unpack procedures long    */
    0,                /* grib_pack procedures double    */
    0,              /* grib_unpack procedures double  */
    &pack_string,                /* grib_pack procedures string    */
    &unpack_string,              /* grib_unpack procedures string  */
    0,                 /* grib_pack procedures bytes     */
    0,               /* grib_unpack procedures bytes   */
    &pack_expression,            /* pack_expression */
    0,              /* notify_change   */
    0,                /* update_size   */
    0,            /* preferred_size   */
    0,                    /* resize   */
    0,      /* nearest_smaller_value */
    0,                       /* next accessor    */
    0,                    /* compare vs. another accessor   */
    0,     /* unpack only ith value          */
    0,     /* unpack a subarray         */
    0,             		/* clear          */
};


grib_accessor_class* grib_accessor_class_codetable = &_grib_accessor_class_codetable;


static void init_class(grib_accessor_class* c)
{
	c->next_offset	=	(*(c->super))->next_offset;
	c->byte_count	=	(*(c->super))->byte_count;
	c->byte_offset	=	(*(c->super))->byte_offset;
	c->sub_section	=	(*(c->super))->sub_section;
	c->pack_missing	=	(*(c->super))->pack_missing;
	c->is_missing	=	(*(c->super))->is_missing;
	c->pack_long	=	(*(c->super))->pack_long;
	c->pack_double	=	(*(c->super))->pack_double;
	c->unpack_double	=	(*(c->super))->unpack_double;
	c->pack_bytes	=	(*(c->super))->pack_bytes;
	c->unpack_bytes	=	(*(c->super))->unpack_bytes;
	c->notify_change	=	(*(c->super))->notify_change;
	c->update_size	=	(*(c->super))->update_size;
	c->preferred_size	=	(*(c->super))->preferred_size;
	c->resize	=	(*(c->super))->resize;
	c->nearest_smaller_value	=	(*(c->super))->nearest_smaller_value;
	c->next	=	(*(c->super))->next;
	c->compare	=	(*(c->super))->compare;
	c->unpack_double_element	=	(*(c->super))->unpack_double_element;
	c->unpack_double_subarray	=	(*(c->super))->unpack_double_subarray;
	c->clear	=	(*(c->super))->clear;
}

/* END_CLASS_IMP */

static int grib_load_codetable(grib_context* c,const char* filename,
           const char* recomposed_name,size_t size,grib_codetable* t); 
static void init(grib_accessor* a, const long len, grib_arguments* params) {
  int n=0;
  grib_accessor_codetable* self  = (grib_accessor_codetable*)a;
  grib_action* act=(grib_action*)(a->creator);

  self->tablename = grib_arguments_get_string(a->parent->h,params,n++);
  self->masterDir = grib_arguments_get_name(a->parent->h,params,n++);
  self->localDir = grib_arguments_get_name(a->parent->h,params,n++);

  /*if (a->flags & GRIB_ACCESSOR_FLAG_STRING_TYPE)
    printf("-------- %s type string (%ld)\n",a->name,a->flags);*/

  if (a->flags & GRIB_ACCESSOR_FLAG_TRANSIENT) {
    a->length = 0;
	if (!a->vvalue) 
		a->vvalue = grib_context_malloc_clear(a->parent->h->context,sizeof(grib_virtual_value));
    a->vvalue->type=grib_accessor_get_native_type(a);
    a->vvalue->length=len;
    if (act->default_value!=NULL) {
      const char* p = 0;
      size_t len = 1;
      long l;
      int ret=0;
      double d;
      char tmp[1024];
      grib_expression* expression=grib_arguments_get_expression(a->parent->h,act->default_value,0);
      int type = grib_expression_native_type(a->parent->h,expression);
      switch(type) {
        case GRIB_TYPE_DOUBLE:
          grib_expression_evaluate_double(a->parent->h,expression,&d);
          grib_pack_double(a,&d,&len);
          break;

        case GRIB_TYPE_LONG:
          grib_expression_evaluate_long(a->parent->h,expression,&l);
          grib_pack_long(a,&l,&len);
          break;

        default:
          len = sizeof(tmp);
          p = grib_expression_evaluate_string(a->parent->h,expression,tmp,&len,&ret);
          if (ret != GRIB_SUCCESS) {
			  grib_context_log(a->parent->h->context,GRIB_LOG_FATAL,
							   "unable to evaluate %s as string",a->name);
          }
          len = strlen(p)+1;
          pack_string(a,p,&len);
          break;
      }
    }
  } else
    a->length = len;

}

static grib_codetable* load_table(grib_accessor_codetable* self)
{
  size_t size = 0;
  grib_handle*    h = ((grib_accessor*)self)->parent->h;
  grib_context*   c = h->context;
  grib_codetable* t = NULL;
  grib_codetable* next=NULL ;
  grib_accessor* a=(grib_accessor*)self;
  char *filename=0;
  char name[1024]={0,};
  char recomposed[1024]={0,};
  char localRecomposed[1024]={0,};
  char *localFilename=0;
  char localName[1024]={0,};
  char masterDir[1024]={0,};
  char localDir[1024]={0,};
  size_t len=1024;

  if (self->masterDir != NULL)
    grib_get_string(h,self->masterDir,masterDir,&len);

  len=1024;
  if (self->localDir != NULL)
    grib_get_string(h,self->localDir,localDir,&len);

  if (*masterDir!=0) {
    sprintf(name,"%s/%s",masterDir,self->tablename);
    grib_recompose_name(h, NULL,name, recomposed,0);
    filename=grib_context_full_path(c,recomposed);
  } else {
    grib_recompose_name(h, NULL,self->tablename, recomposed,0);
    filename=grib_context_full_path(c,recomposed);
  }

  if (*localDir!=0) {
    sprintf(localName,"%s/%s",localDir,self->tablename);
    grib_recompose_name(h, NULL,localName, localRecomposed,0);
    localFilename=grib_context_full_path(c,localRecomposed);
  }
  
  next=c->codetable;
  while(next) {
    if((filename && next->filename[0] && strcmp(filename,next->filename[0]) == 0) &&
       ((localFilename==0 && next->filename[1]==NULL) ||
           ((localFilename!=0 && next->filename[1]!=NULL)
           && strcmp(localFilename,next->filename[1]) ==0)) )
      return next;
    next = next->next;
  }

  if (a->flags & GRIB_ACCESSOR_FLAG_TRANSIENT) {
	  Assert(a->vvalue!=NULL);
    size=a->vvalue->length*8;
  } else {
    size = grib_byte_count((grib_accessor*)self) * 8;
  }
  size = grib_power(size,2);

  t = (grib_codetable*)grib_context_malloc_clear_persistent(c,sizeof(grib_codetable) +
      (size-1)*sizeof(code_table_entry));

  if (filename!=0) grib_load_codetable(c,filename,recomposed,size,t);

  if (localFilename!=0) grib_load_codetable(c,localFilename,localRecomposed,size,t);

  if (t->filename[0]==NULL && t->filename[1]==NULL) {
    grib_context_free_persistent(c,t);
    return NULL;
  }
  
  return t;

}

static int grib_load_codetable(grib_context* c,const char* filename,
           const char* recomposed_name,size_t size,grib_codetable* t) {
  char line[1024];
  FILE *f = NULL;

  grib_context_log(c,GRIB_LOG_DEBUG,"Loading code table form %s",filename);

  f=fopen(filename, "r");
  if (!f) return GRIB_IO_PROBLEM;

  Assert(t!=NULL);

  if (t->filename[0] == NULL ){
    t->filename[0]  = grib_context_strdup_persistent(c,filename);
    t->recomposed_name[0]  = grib_context_strdup_persistent(c,recomposed_name);
    t->next      = c->codetable;
    t->size      = size;
    c->codetable = t;
  } else {
    t->filename[1]  = grib_context_strdup_persistent(c,filename);
    t->recomposed_name[1]  = grib_context_strdup_persistent(c,recomposed_name);
  }

  while(fgets(line,sizeof(line)-1,f))
  {
    char* p = line;
    int code    = 0;
    char abbreviation[1024] = {0,};
    char title[1024]={0,};
    char* q = abbreviation;
    char* r = title;
    char* units=0;
    char unknown[]="unknown";

    line[strlen(line)-1] = 0;

    while(*p != '\0' && isspace(*p)) p++;

    if(*p == '#')
      continue;

    while(*p != '\0' && isspace(*p)) p++;

    if( *p =='\0' ) continue;

    Assert(isdigit(*p));

    while(*p != '\0')
    {
      if(isspace(*p)) break;
      code *= 10;
      code += *p - '0';
      p++;
    }

    if(code <0 || code >= size)
    {
      grib_context_log(c,GRIB_LOG_WARNING,"code_table_entry: invalide code in %s: %d (table size=%d)",filename,code,size);
      continue;
    }

    while(*p != '\0' && isspace(*p)) p++;

    while(*p != '\0')
    {
      if(isspace(*p)) break;
      *q++ = *p++;
    }
    *q = 0;
    while(*p != '\0' && isspace(*p)) p++;

    while(*p != '\0')
    {
      if(*p == '(' ) break;
      *r++ = *p++;
    }
    *r = 0;

    while(*p != '\0' && isspace(*p)) p++;
    if (*p != '\0') {
      units=++p;
      while(*p != '\0' && *p != ')' ) p++;
      *p='\0';
    } else {
      units=unknown;
    }

    Assert(*abbreviation);
    Assert(*title);

    if(t->entries[code].abbreviation != NULL)
    {
      grib_context_log(c,GRIB_LOG_WARNING,"code_table_entry: duplicate code in %s: %d (table size=%d)",filename,code,size);
      continue;
    }

    Assert(t->entries[code].abbreviation == NULL);
    Assert(t->entries[code].title        == NULL);

    t->entries[code].abbreviation = grib_context_strdup_persistent(c,abbreviation);
    t->entries[code].title        = grib_context_strdup_persistent(c,title);
    t->entries[code].units        = grib_context_strdup_persistent(c,units);


  }

  fclose(f);

  return 0;

}


void grib_codetable_delete(grib_context* c) {
  grib_codetable* t = c->codetable;

  while(t)
  {
    grib_codetable* s = t->next;
    int i;

    for(i = 0; i < t->size; i++)
    {
      grib_context_free_persistent(c,t->entries[i].abbreviation);
      grib_context_free_persistent(c,t->entries[i].title);
    }
    grib_context_free_persistent(c,t->filename[0]);
    if(t->filename[1])
      grib_context_free_persistent(c,t->filename[1]);
    grib_context_free_persistent(c,t->recomposed_name[0]);
    if (t->recomposed_name[1])
      grib_context_free_persistent(c,t->recomposed_name[1]);
    grib_context_free_persistent(c,t);
    t = s;
  }

}



static void dump(grib_accessor* a, grib_dumper* dumper) {
  grib_accessor_codetable* self  = (grib_accessor_codetable*)a;
  char comment[2048];
  grib_codetable* table;

  size_t llen = 1;
  long value;

  if(!self->table) self->table = load_table(self);
  table=self->table;

  grib_unpack_long(a, &value,&llen);

  if(value == GRIB_MISSING_LONG)
  {
    if(a->length < 4)
    {
      value = (1L << a->length) - 1;
    }
  }

  if(table && value >= 0 && value < table->size)
  {
    if(table->entries[value].abbreviation)
    {
      int b = atol(table->entries[value].abbreviation);
      if(b == value)
        strcpy(comment,table->entries[value].title);
      else
        sprintf(comment,"%s", table->entries[value].title);

      if (table->entries[value].units!=NULL && strcmp(table->entries[value].units,"unknown")) {
        strcat(comment," (");
        strcat(comment,table->entries[value].units);
        strcat(comment,") ");
      }
    }
    else
    {
      strcpy(comment,"Unknown code table entry");
    }

  }
  else
  {
    strcpy(comment,"Unknown code table entry");
  }

  strcat(comment," (");
  if (table) {
    strcat(comment,table->recomposed_name[0]);
    if (table->recomposed_name[1]!=NULL) {
      strcat(comment," , ");
      strcat(comment,table->recomposed_name[1]);
    }
  }
  strcat(comment,") ");

  grib_dump_long(dumper,a,comment);

}

static int unpack_string (grib_accessor* a, char* buffer, size_t *len)
{
  grib_accessor_codetable* self = (grib_accessor_codetable*)a;
  grib_codetable*          table = NULL;

  size_t size = 1;
  long   value;
  int err = GRIB_SUCCESS;
  char tmp[1024];
  size_t l = 0;

  if( (err = grib_unpack_long(a,&value,&size)) != GRIB_SUCCESS)
    return err;

  if(!self->table) self->table = load_table(self);
  table=self->table;

  if(table && (value >= 0) && (value < table->size) && table->entries[value].abbreviation)
  {
    strcpy(tmp,table->entries[value].abbreviation);
  }
  else
  {

#if 1
    sprintf(tmp,"%d",(int)value);
#else
    return GRIB_DECODING_ERROR;
#endif
  }


  l = strlen(tmp) + 1;

  if(*len < l)
  {
    *len = l;
    return GRIB_BUFFER_TOO_SMALL;
  }

  strcpy(buffer,tmp);
  *len = l;

  return GRIB_SUCCESS;
}

static long value_count(grib_accessor* a)
{
  return 1;
}

static int pack_string(grib_accessor* a, const char* buffer, size_t *len)
{
  grib_accessor_codetable* self = (grib_accessor_codetable*)a;
  grib_codetable*          table ;

  long i;
  size_t size = 1;

  typedef int (*cmpproc)(const char*, const char*);

  cmpproc cmp = a->flags | GRIB_ACCESSOR_FLAG_LOWERCASE ? strcasecmp : strcmp;

  if(!self->table) self->table = load_table(self);
  table=self->table;

  if(!table)
    return GRIB_ENCODING_ERROR;

  if (a->set) {
	  int err=grib_set_string(a->parent->h,a->set,buffer,len);
	  if (err!=0) return err;
  }

  for(i = 0 ; i < table->size; i++)
    if(table->entries[i].abbreviation)
		  if(cmp(table->entries[i].abbreviation,buffer) == 0)
        return grib_pack_long(a,&i,&size);

  if (a->flags & GRIB_ACCESSOR_FLAG_NO_FAIL) {
	  grib_action* act=(grib_action*)(a->creator);
	  if (act->default_value!=NULL) {
		  const char* p = 0;
		  size_t len = 1;
		  long l;
		  int ret=0;
		  double d;
		  char tmp[1024];
		  grib_expression* expression=grib_arguments_get_expression(a->parent->h,act->default_value,0);
		  int type = grib_expression_native_type(a->parent->h,expression);
		  switch(type) {
			  case GRIB_TYPE_DOUBLE:
				  grib_expression_evaluate_double(a->parent->h,expression,&d);
				  grib_pack_double(a,&d,&len);
				  break;

			  case GRIB_TYPE_LONG:
				  grib_expression_evaluate_long(a->parent->h,expression,&l);
				  grib_pack_long(a,&l,&len);
				  break;

			  default:
				  len = sizeof(tmp);
				  p = grib_expression_evaluate_string(a->parent->h,expression,tmp,&len,&ret);
				  if (ret != GRIB_SUCCESS) {
					  grib_context_log(a->parent->h->context,GRIB_LOG_FATAL,
									   "unable to evaluate %s as string",a->name);
					  return ret;
				  }
				  len = strlen(p)+1;
				  pack_string(a,p,&len);
				  break;
		  }
		  return GRIB_SUCCESS;
	  }
	
  }
  return GRIB_ENCODING_ERROR;
}

static int pack_expression(grib_accessor* a, grib_expression *e){
  const char* cval;
  int ret=0;
  long lval=0;
  size_t len = 1;
  char tmp[1024];

  if (strcmp(e->cclass->name,"long")==0) {
    ret=grib_expression_evaluate_long(a->parent->h,e,&lval);
    ret = grib_pack_long(a,&lval,&len);
  } else {
    len = sizeof(tmp);
    cval = grib_expression_evaluate_string(a->parent->h,e,tmp,&len,&ret);
    if (ret!=GRIB_SUCCESS) {
      grib_context_log(a->parent->h->context,GRIB_LOG_ERROR,"grib_accessor_codetable.pack_expression: unable to evaluate string %s to be set in %s\n",grib_expression_get_name(e),a->name);
	  return ret;
    }
    len = strlen(cval) + 1;
    ret = grib_pack_string(a,cval,&len);
  }
  return ret;
}

static void destroy(grib_context* context,grib_accessor* a)
{
	if (a->vvalue != NULL) {
		grib_context_free(context, a->vvalue);
		a->vvalue=NULL;
	}

}

static int  get_native_type(grib_accessor* a){
  int type=GRIB_TYPE_LONG;
  /*printf("---------- %s flags=%ld GRIB_ACCESSOR_FLAG_STRING_TYPE=%d\n",
         a->name,a->flags,GRIB_ACCESSOR_FLAG_STRING_TYPE);*/
  if (a->flags & GRIB_ACCESSOR_FLAG_STRING_TYPE)
    type=GRIB_TYPE_STRING;
  return type;
}

static int    unpack_long   (grib_accessor* a, long* val, size_t *len)
{
  grib_accessor_codetable* self = (grib_accessor_codetable*)a;
  unsigned long rlen = grib_value_count(a);
  unsigned long i = 0;
  long pos = a->offset*8;

  if(!self->table) self->table = load_table(self);

  if(*len < rlen)
  {
    grib_context_log(a->parent->h->context, GRIB_LOG_ERROR, " wrong size (%ld) for %s it contains %d values ",*len, a->name , rlen);
    *len = 0;
    return GRIB_ARRAY_TOO_SMALL;
  }

  if (a->flags & GRIB_ACCESSOR_FLAG_TRANSIENT) {
    *val=a->vvalue->lval;
    *len=1;
    return GRIB_SUCCESS;
  }

  for(i=0; i< rlen;i++){
    val[i] = (long)grib_decode_unsigned_long(a->parent->h->buffer->data , &pos, self->nbytes*8);
  }

  *len = rlen;
  return GRIB_SUCCESS;
}


