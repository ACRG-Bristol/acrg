/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 *
 * Description: file pool
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */

#include "grib_api_internal.h"
#define GRIB_MAX_OPENED_FILES 200

static short next_id=0;

GRIB_INLINE static int grib_inline_strcmp(const char* a,const char* b) {
  if (*a != *b) return 1;
  while((*a!=0 && *b!=0) &&  *(a) == *(b) ) {a++;b++;}
  return (*a==0 && *b==0) ? 0 : 1;
}

static grib_file_pool file_pool= {
             0,                    /* grib_context* context;*/
             0,                    /* grib_file* first;*/
             0,                    /* grib_file* current; */
             0,                    /* size_t size;*/
             0,                    /* int number_of_opened_files;*/
  GRIB_MAX_OPENED_FILES            /* int max_opened_files; */
};

void grib_file_pool_clean() {
   grib_file *file,*next;

   if (!file_pool.first) return;

   file=file_pool.first;
   while(file) {
     next=file->next;
     grib_file_delete(file);
     file=next;
   }
}

static void grib_file_pool_change_id() {
   grib_file *file;

   if (!file_pool.first) return;

   file=file_pool.first;
   while(file) {
	 file->id+=1000;
     file=file->next;
   }
}

static grib_file* grib_read_file(grib_context *c,FILE* fh,int *err) {
	short marker=0;
	short id=0;
	grib_file* file;
	*err = grib_read_short(fh,&marker);
	if(!marker) return NULL;

	file=grib_context_malloc_clear(c,sizeof(grib_file));
	file->buffer=0;
	file->name=grib_read_string(c,fh,err);
	if (*err) return NULL;

	*err=grib_read_short(fh,&id);
	file->id=id;
	if (*err) return NULL;

	file->next=grib_read_file(c,fh,err);
	if (*err) return NULL;

	return file;
}

static int grib_write_file(FILE *fh,grib_file* file) {
	int err=0;

	if (!file)
		return grib_write_null_marker(fh);

	err=grib_write_not_null_marker(fh);
	if (err) return err;

	err=grib_write_string(fh,file->name);
	if (err) return err;

	err=grib_write_short(fh,(short)file->id);
	if (err) return err;

	return grib_write_file(fh,file->next);
}

grib_file* grib_file_pool_get_files() {
	return file_pool.first;
}

int grib_file_pool_read(grib_context* c,FILE* fh) {
	int err=0;
	short marker=0;
	grib_file* file;

	if (!c) c=grib_context_get_default();

	err = grib_read_short(fh,&marker);
	if(!marker) {
		grib_context_log(c,GRIB_LOG_ERROR,
		"Unable to find file information in index file\n");
		return GRIB_INVALID_FILE;
	}

	grib_file_pool_change_id();
	file=file_pool.first;

	while (file->next) 
		file=file->next;
	 
	file->next=grib_read_file(c,fh,&err);
	if (err) return err;

	return GRIB_SUCCESS;
}

int grib_file_pool_write(FILE* fh) {
	int err=0;
	if (!file_pool.first)
		return grib_write_null_marker(fh);

	err=grib_write_not_null_marker(fh);
	if (err) return err;

	return grib_write_file(fh,file_pool.first);
}

grib_file* grib_file_open(const char* filename, const char* mode,int* err) {
  grib_file *file=0,*prev=0;
  int same_mode=0;
  int is_new=0;

  if (!file_pool.context) file_pool.context=grib_context_get_default();

  if (file_pool.current && !grib_inline_strcmp(filename,file_pool.current->name)) {
    file=file_pool.current;
  } else {
    file=file_pool.first;
    while (file) {
      if (!grib_inline_strcmp(filename,file->name)) break;
      prev=file;
      file=file->next;
    }
    if (!file) {
      is_new=1;
      file=grib_file_new(file_pool.context,filename,err);
      if (prev) prev->next=file;
      file_pool.current=file;
      if (!prev) file_pool.first=file;
      file_pool.size++;
    }
  }

  if (file->mode) same_mode=grib_inline_strcmp(mode,file->mode) ? 0 : 1;
  if (file->handle && same_mode) {
    *err=0;
    return file;
  }

  if (!same_mode && file->handle) {
    /*printf("========== mode=%s file->mode=%s\n",mode,file->mode);*/
    fclose(file->handle);
  }

  if (!file->handle) {
   /*printf("-- opening file %s %s\n",file->name,mode);*/
   if (!is_new && *mode == 'w')
    file->handle = fopen(file->name,"a");
   else
    file->handle = fopen(file->name,mode);

   file->mode=strdup(mode);
   if (!file->handle) {
     grib_context_log(file->context,GRIB_LOG_PERROR,"grib_file_open: cannot open file %s",file->name);
     *err=GRIB_IO_PROBLEM;
     return NULL;
   }
   if (file_pool.context->io_buffer_size) {
#ifdef POSIX_MEMALIGN
      if (posix_memalign((void**)&(file->buffer),sysconf(_SC_PAGESIZE),file_pool.context->io_buffer_size) ) {
        grib_context_log(file->context,GRIB_LOG_FATAL,"posix_memalign unable to allocate io_buffer\n");
      }
#else
      file->buffer = (void*)malloc(file_pool.context->io_buffer_size);
      if (!file->buffer) {
        grib_context_log(file->context,GRIB_LOG_FATAL,"Unable to allocate io_buffer\n");
      }
#endif
      setvbuf(file->handle,file->buffer,_IOFBF,file_pool.context->io_buffer_size);
   }

   file_pool.number_of_opened_files++;
#if 0
   if (file_pool.number_of_opened_files > file_pool.max_opened_files) {
     grib_file *f=file_pool.first;
     while (f) {
        if (f != file && f->handle) {
          fclose(f->handle);
          break;
        }
     }
   }
#endif
  }
  return file;
}

void grib_file_close(const char* filename,int* err) {
  grib_file* file=NULL;

  if (file_pool.number_of_opened_files > GRIB_MAX_OPENED_FILES) {
    /*printf("++ closing file %s\n",filename);*/
    file=grib_get_file(filename,err);
    fclose(file->handle);
	if (file->buffer) {
		free(file->buffer);
		file->buffer=0;
	}
    file->handle=NULL;
    file_pool.number_of_opened_files--;
  }
}

grib_file* grib_get_file(const char* filename,int* err) {
   grib_file* file=NULL;

   if (file_pool.current->name && !grib_inline_strcmp(filename,file_pool.current->name)) {
      return file_pool.current;
   }

   file=file_pool.first;
   while (file) {
     if (!grib_inline_strcmp(filename,file->name)) break;
     file=file->next;
   }
   if (!file) file=grib_file_new(0,filename,err);

   return file;
}

grib_file* grib_find_file(short id) {
   grib_file* file=NULL;

   if (file_pool.current->name && id==file_pool.current->id) {
      return file_pool.current;
   }

   file=file_pool.first;
   while (file) {
     if (id==file->id) break;
     file=file->next;
   }

   return file;
}

grib_file* grib_file_new(grib_context* c, const char* name, int* err) {
   grib_file* file;

   if (!c) c=grib_context_get_default( );
   file=grib_context_malloc_clear( c,sizeof(grib_file));

   if (!file) {
     grib_context_log(c,GRIB_LOG_ERROR,"grib_file_new: unable to allocate memory");
     *err=GRIB_OUT_OF_MEMORY;
     return NULL;
   }
   file->name=strdup(name);
   file->id=next_id;
   next_id++;
   file->mode=0;
   file->handle=0;
   file->refcount=0;
   file->context=c;
   file->next=0;
   file->buffer=0;
   return file;
}

void grib_file_delete(grib_file* file) {
   if (!file) return;
   if(file->name) free(file->name);
   if (file->mode) free(file->mode);
   if (file->handle) {
	   fclose(file->handle);
   }
   if (file->buffer) {
   	  free(file->buffer);
   }
   grib_context_free(file->context,file);
}

