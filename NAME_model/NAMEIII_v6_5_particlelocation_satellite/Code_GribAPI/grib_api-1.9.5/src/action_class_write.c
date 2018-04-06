/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/***************************************************************************
 *  Enrico Fucile                                                          *
 ***************************************************************************/
#include "grib_api_internal.h"
/*
   This is used by make_class.pl

   START_CLASS_DEF
   CLASS      = action
   IMPLEMENTS = dump
   IMPLEMENTS = destroy;execute
   MEMBERS    = char *name
   MEMBERS    = int append;
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
static void destroy         (grib_context*,grib_action*);
static int execute(grib_action* a,grib_handle* h);


typedef struct grib_action_write {
    grib_action          act;  
/* Members defined in write */
	char *name;
	int append;
} grib_action_write;


static grib_action_class _grib_action_class_write = {
    0,                              /* super                     */
    "action_class_write",                              /* name                      */
    sizeof(grib_action_write),            /* size                      */
    0,                                   /* inited */
    &init_class,                         /* init_class */
    0,                               /* init                      */
    &destroy,                            /* destroy */

    &dump,                               /* dump                      */
    0,                               /* xref                      */

    0,             /* create_accessor*/

    0,                            /* notify_change */
    0,                            /* reparse */
    &execute,                            /* execute */
    0,                            /* compile */
};

grib_action_class* grib_action_class_write = &_grib_action_class_write;

static void init_class(grib_action_class* c)
{
}
/* END_CLASS_IMP */

extern int errno;

grib_action* grib_action_create_write( grib_context* context, const char* name,int append)
{
  char buf[1024];

  grib_action_write* a =NULL;
  grib_action_class* c   = grib_action_class_write;
  grib_action* act       = (grib_action*)grib_context_malloc_clear_persistent(context,c->size);
  act->op              =  grib_context_strdup_persistent(context,"section");

  act->cclass            = c;
  a                      = (grib_action_write*)act;
  act->context           = context;

  a->name            = grib_context_strdup_persistent(context,name);

  sprintf(buf,"write%p",(void*)a->name);

  act->name      = grib_context_strdup_persistent(context,buf);
  a->append=append;

  return act;
}

static int execute(grib_action* act, grib_handle *h)
{
  int ioerr=0;
  grib_action_write* a = (grib_action_write*) act;
  int err =GRIB_SUCCESS;
  size_t size;
  const void* buffer=NULL;
  char* filename;
  char string[1024]={0,};

  grib_file* of=NULL;

  if ((err=grib_get_message(h,&buffer,&size))!= GRIB_SUCCESS) {
    grib_context_log(act->context,GRIB_LOG_ERROR,"unable to get message\n");
    return err;
  }

  if (strlen(a->name)!=0) {
    err = grib_recompose_name(h,NULL,a->name,string,0);
    filename=string;
  } else {
    filename = act->context->outfilename ? act->context->outfilename : "filter.out";
  }

  if (a->append) of=grib_file_open(filename,"a",&err);
  else of=grib_file_open(filename,"w",&err);

  if (!of || !of->handle) {
    grib_context_log(act->context,GRIB_LOG_ERROR,"unable to open file %s\n",filename);
    return GRIB_IO_PROBLEM;
  }

  if (h->gts_header)
          fwrite(h->gts_header,1,h->gts_header_len,of->handle);

  if(fwrite(buffer,1,size,of->handle) != size) {
      ioerr=errno;
      grib_context_log(act->context,(GRIB_LOG_ERROR)|(GRIB_LOG_PERROR),
                   "Error writing to %s",filename);
      return GRIB_IO_PROBLEM;
  }

  if (h->gts_header) {
        char gts_trailer[4]={'\x0D','\x0D','\x0A','\x03'};
        fwrite(gts_trailer,1,4,of->handle);
  }

  grib_file_close(filename,&err);

  if (err != GRIB_SUCCESS) {
    grib_context_log(act->context,GRIB_LOG_ERROR,"unable to get message\n");
    return err;
  }

  return err;
}


static void dump(grib_action* act, FILE* f, int lvl)
{
}


static void destroy(grib_context* context,grib_action* act)
{
  grib_action_write* a = (grib_action_write*) act;

  grib_context_free_persistent(context, a->name);
  grib_context_free_persistent(context, act->name);
  grib_context_free_persistent(context, act->op);

}
