/**
 * Copyright 2005-2007 ECMWF
 *
 * Licensed under the GNU Lesser General Public License which
 * incorporates the terms and conditions of version 3 of the GNU
 * General Public License.
 * See LICENSE and gpl-3.0.txt for details.
 */

#include "grib_api_internal.h"
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>

#if MANAGE_MEM

#else
static long cnt = 0;
static long cntp = 0;

static void  default_long_lasting_free(const grib_context* c, void* p)   {free(p);cntp--;}
static void* default_long_lasting_malloc(const grib_context* c, size_t size) {cntp++;return malloc(size);}

static void  default_buffer_free(const grib_context* c, void* p)   {free(p);cntp--;}
static void* default_buffer_malloc(const grib_context* c, size_t size) {cntp++;return malloc(size);}
static void* default_buffer_realloc(const grib_context* c, void* p, size_t size) {return realloc(p,size);}

static void  default_free(const grib_context* c, void* p)   {free(p);cnt--;}
static void* default_malloc(const grib_context* c, size_t size) {cnt++;return malloc(size);}
static void* default_realloc(const grib_context* c, void* p, size_t size) {return realloc(p,size);}
#endif

static size_t  default_read(const grib_context* c, void *ptr, size_t size, void *stream) {
	return fread(ptr,  1, size,  stream);
}

static off_t default_tell(const grib_context* c, void *stream) {
	return ftello(stream);
}

static off_t default_seek(const grib_context* c, off_t offset,int whence, void *stream) {
	return fseeko(stream,offset,whence);
}

static int default_feof(const grib_context* c, void *stream) {
	return feof((FILE*)stream);
}
static size_t  default_write(const grib_context* c,const void *ptr, size_t size, void *stream)  {
	return fwrite(ptr,  1, size,  stream);
}

size_t  grib_context_read(const grib_context* c, void *ptr, size_t size, void *stream) {
	if (!c) c=grib_context_get_default();
	return c->read(c,ptr,  size,  stream);
}

off_t  grib_context_tell(const grib_context* c, void *stream) {
	if (!c) c=grib_context_get_default();
	return c->tell(c,stream);
}

int    grib_context_seek(const grib_context* c, off_t offset ,int whence , void *stream) {
	if (!c) c=grib_context_get_default();
	return c->seek(c,offset,whence,stream);
}

int    grib_context_eof(const grib_context* c, void *stream) {
	if (!c) c=grib_context_get_default();
	return c->eof(c,stream);
}
size_t  grib_context_write(const grib_context* c,const void *ptr, size_t size, void *stream)  {
	if (!c) c=grib_context_get_default();
	return c->write(c,ptr,  size,  stream);
}

static void  default_log(const grib_context* c, int level, const char* mess){
	if (!c) c=grib_context_get_default();
	if(level == GRIB_LOG_ERROR)   {
		fprintf(stderr, "GRIB_API ERROR   :  %s\n", mess);
		/*Assert(1==0);*/
	}
	if(level == GRIB_LOG_FATAL)   fprintf(stderr, "GRIB_API ERROR   :  %s\n", mess);
	if(level == GRIB_LOG_DEBUG && c->debug>0)   fprintf(stderr, "GRIB_API DEBUG   :  %s\n", mess);
	if(level == GRIB_LOG_WARNING) fprintf(stderr, "GRIB_API WARNING :  %s\n", mess);
	if(level == GRIB_LOG_INFO)    fprintf(stderr, "GRIB_API INFO    :  %s\n", mess);

	if(level == GRIB_LOG_FATAL) { Assert(0);}

	if(getenv("GRIB_API_FAIL_IF_LOG_MESSAGE"))
	{
		long n = atol(getenv("GRIB_API_FAIL_IF_LOG_MESSAGE"));
		if(n >= 1 && level == GRIB_LOG_ERROR) exit(1);
		if(n >= 2 && level == GRIB_LOG_WARNING) exit(1);
	}

}

static void  default_print(const grib_context* c, void* descriptor, const char* mess)
{
	fprintf(descriptor, "%s", mess);
}

void grib_context_set_print_proc(grib_context* c, grib_print_proc p)
{
	c = c ? c : grib_context_get_default();
	c->print = p;
}

void grib_context_set_debug(grib_context* c, int mode)
{
	c = c ? c : grib_context_get_default();
	c->debug = mode;
}

void grib_context_set_logging_proc(grib_context* c, grib_log_proc p)
{
	c = c ? c : grib_context_get_default();
	c->output_log = p;
}


long grib_get_api_version(){
	return GRIB_API_VERSION;
}

void grib_print_api_version(FILE* out) {
	fprintf(out,"%d.%d.%d",
			GRIB_API_MAJOR_VERSION,
			GRIB_API_MINOR_VERSION,
			GRIB_API_REVISION_VERSION);
}

static  grib_context default_grib_context = {
	0,                            /* inited                     */
	0,                            /* debug                      */
	0,                            /* write_on_fail              */
	0,                            /* no_abort                   */
	0,                            /* io_buffer_size             */
	0,                            /* grib_definition_files_path */
	0,                            /* grib_templates_path        */
	0,                            /* grib_concept_path          */
	0,                            /* grib_reader                */
	0,                            /* user data                  */
	GRIB_REAL_MODE8,              /* real mode for fortran      */

#if MANAGE_MEM
	&grib_transient_free,                /* free_mem                   */
	&grib_transient_malloc,              /* alloc_mem                  */
	&grib_transient_realloc,              /* realloc_mem                  */

	&grib_permanent_free,   /* free_persistant_mem        */
	&grib_permanent_malloc, /* alloc_persistant_mem       */

	&grib_buffer_free,                /* buffer_free_mem                   */
	&grib_buffer_malloc,              /* buffer_alloc_mem                  */
	&grib_buffer_realloc,              /* buffer_realloc_mem                  */

#else

	&default_free,                /* free_mem                   */
	&default_malloc,              /* alloc_mem                  */
	&default_realloc,             /* realloc_mem                */

	&default_long_lasting_free,   /* free_persistant_mem        */
	&default_long_lasting_malloc, /* alloc_persistant_mem       */

	&default_buffer_free,         /* free_buffer_mem            */
	&default_buffer_malloc,        /* alloc_buffer_mem           */
	&default_buffer_realloc,        /* realloc_buffer_mem           */
#endif

	&default_read,                /* file read procedure        */
	&default_write,               /* file write procedure       */
	&default_tell,                /* lfile tell procedure       */
	&default_seek,                /* lfile seek procedure       */
	&default_feof,                /* file feof procedure        */

	&default_log,                 /* logging_procedure          */
	&default_print,               /* print procedure            */
	0,                            /* code tables                */
	0,                            /* files                      */
	0,                            /* multigrib support on       */
	0,                            /* multigrib support          */
	0,                            /* grib_definition_files_dir  */
	0,                            /* handle_file_count          */
	0,                            /* handle_total_count         */
	0,                            /* message_file_offset        */
	0,                            /* no_fail_on_wrong_length    */
	0,                            /* gts_header_on */
	0,                            /* gribex_mode_on */
	0,                            /* keys (grib_trie*) */
	0,                            /* keys_count */
	0,                            /* concepts_index */
	0,                            /* concepts_count */
	{0,},                         /* concepts */
	0,                            /* def_files */
	0,                            /* ieee_packing */
	0,                            /* blacklist */
    0                             /* classes */
#if GRIB_PTHREADS
		,PTHREAD_MUTEX_INITIALIZER     /*mutex*/
#endif
};

#if GRIB_PTHREADS
static pthread_once_t once  = PTHREAD_ONCE_INIT;
static void init() {
	pthread_mutexattr_t attr;
	pthread_mutexattr_init(&attr);
	pthread_mutexattr_settype(&attr,PTHREAD_MUTEX_RECURSIVE);
	pthread_mutex_init(&(default_grib_context.mutex),&attr);
	pthread_mutexattr_destroy(&attr);

}
#endif

grib_context* grib_context_get_default(){
	GRIB_PTHREAD_ONCE(&once,&init);
	GRIB_MUTEX_LOCK(&(default_grib_context.mutex));

	if(!default_grib_context.inited)
	{
		const char * write_on_fail = getenv("GRIB_API_WRITE_ON_FAIL");
		const char * no_abort = getenv("GRIB_API_NO_ABORT");
		const char * debug = getenv("GRIB_API_DEBUG");
		const char *gribex=getenv("GRIB_GRIBEX_MODE_ON");
		const char *ieee_packing=getenv("GRIB_IEEE_PACKING");
		const char *io_buffer_size=getenv("GRIB_API_IO_BUFFER_SIZE");
		default_grib_context.inited = 1;
		default_grib_context.io_buffer_size = io_buffer_size ? atoi(io_buffer_size) : 0;
		default_grib_context.write_on_fail  = write_on_fail ? atoi(write_on_fail) : 0;
		default_grib_context.no_abort  = no_abort ? atoi(no_abort) : 0;
		default_grib_context.debug  = debug ? atoi(debug) : 0;
		default_grib_context.gribex_mode_on=gribex ? atoi(gribex) : 0;
		default_grib_context.ieee_packing=ieee_packing ? atoi(ieee_packing) : 0;
		default_grib_context.grib_templates_path = getenv("GRIB_SAMPLES_PATH");
		if (!default_grib_context.grib_templates_path)
			default_grib_context.grib_templates_path = getenv("GRIB_TEMPLATES_PATH");
#ifdef GRIB_TEMPLATES_PATH
		if(!default_grib_context.grib_templates_path)
			default_grib_context.grib_templates_path = GRIB_TEMPLATES_PATH ;
#endif
		default_grib_context.grib_definition_files_path = getenv("GRIB_DEFINITION_PATH");
#ifdef GRIB_DEFINITION_PATH
		if(!default_grib_context.grib_definition_files_path)
			default_grib_context.grib_definition_files_path = GRIB_DEFINITION_PATH ;
#endif
		default_grib_context.keys_count=0;
		default_grib_context.keys=grib_hash_keys_new(&(default_grib_context),
				&(default_grib_context.keys_count));

		default_grib_context.concepts_index=grib_itrie_new(&(default_grib_context),
				&(default_grib_context.concepts_count));
		default_grib_context.def_files=grib_trie_new(&(default_grib_context));
		default_grib_context.classes=grib_trie_new(&(default_grib_context));
	}


	GRIB_MUTEX_UNLOCK(&(default_grib_context.mutex));

	return &default_grib_context;
}

/* TODO: use parent */

grib_context* grib_context_new(grib_context* parent)
{
	grib_context* c;
#if GRIB_PTHREADS
	pthread_mutexattr_t attr;
#endif

	if (!parent) parent=grib_context_get_default();

	GRIB_MUTEX_LOCK(&(parent->mutex));

	c = (grib_context*)grib_context_malloc_clear_persistent(&default_grib_context,sizeof(grib_context));

	c->inited              = default_grib_context.inited;
	c->debug               = default_grib_context.debug;

	c->real_mode           = default_grib_context.real_mode;

	c->free_mem            = default_grib_context.free_mem;
	c->alloc_mem           = default_grib_context.alloc_mem;

	c->free_persistent_mem = default_grib_context.free_persistent_mem;
	c->alloc_persistent_mem= default_grib_context.alloc_persistent_mem;

	c->read                = default_grib_context.read;
	c->write               = default_grib_context.write;
	c->tell                = default_grib_context.tell;

	c->output_log          = default_grib_context.output_log;
	c->print               = default_grib_context.print    ;
	c->user_data           = default_grib_context.user_data;
	c->def_files           = default_grib_context.def_files;
	
#if GRIB_PTHREADS
	pthread_mutexattr_settype(&attr,PTHREAD_MUTEX_RECURSIVE);
	pthread_mutex_init(&(c->mutex),&attr);
	pthread_mutexattr_destroy(&attr);
#endif

	GRIB_MUTEX_UNLOCK(&(parent->mutex));
	return c;
}

static int init_definition_files_dir(grib_context* c) {
	int err=0;
	char* path=NULL;
	char* p=NULL;
	char* dir=NULL;
	grib_string_list* next=NULL;

	if (!c) c=grib_context_get_default();

	if (c->grib_definition_files_dir) return 0;
	if (!c->grib_definition_files_path) return GRIB_NO_DEFINITIONS;

	path=c->grib_definition_files_path;

	GRIB_MUTEX_LOCK(&(c->mutex));

	p=path;
	while(*p!=':' && *p!='\0') p++;

	if (*p != ':') {
		c->grib_definition_files_dir=(grib_string_list*)grib_context_malloc_clear_persistent(c,sizeof(grib_string_list));
		c->grib_definition_files_dir->value=grib_context_strdup(c,path);
	} else {
		dir=strtok(path,":");

		while (dir != NULL) {
			if (next) {
				next->next=(grib_string_list*)grib_context_malloc_clear_persistent(c,sizeof(grib_string_list));
				next=next->next;
			} else {
				c->grib_definition_files_dir=(grib_string_list*)grib_context_malloc_clear_persistent(c,sizeof(grib_string_list));
				next=c->grib_definition_files_dir;
			}
			next->value=grib_context_strdup(c,dir);
			dir=strtok(NULL,":");
		}
	}

	GRIB_MUTEX_UNLOCK(&(c->mutex));

	return err;
}

char *grib_context_full_path(grib_context* c,const char* basename)
{
	int err=0;
	char full[1024]={0,};
	grib_string_list* dir=NULL;
	grib_string_list* fullpath=0;
	if (!c) c=grib_context_get_default();

	GRIB_MUTEX_LOCK(&(c->mutex));

	if(*basename == '/' || *basename ==  '.') {
		return (char*)basename;
	} else {
		if ((fullpath=grib_trie_get(c->def_files,basename))!=NULL) {
			return fullpath->value;
		}
		if (!c->grib_definition_files_dir)
			err=init_definition_files_dir(c);

		if (err != GRIB_SUCCESS) {
			grib_context_log(c,GRIB_LOG_ERROR,
					"Unable to find definition files directory");
			return NULL;
		}

		dir=c->grib_definition_files_dir;

		while (dir) {
			sprintf(full,"%s/%s",dir->value,basename);
			if (!access(full,F_OK)) {
				fullpath=grib_context_malloc_clear_persistent(c,sizeof(grib_string_list));
				Assert(fullpath);
				fullpath->value=grib_context_strdup(c,full);
				grib_trie_insert(c->def_files,basename,fullpath);
				grib_context_log(c,GRIB_LOG_DEBUG,"Found def file %s",full);
				GRIB_MUTEX_UNLOCK(&(c->mutex));
				return fullpath->value;
			}
			dir=dir->next;
		}

	}

	/* Store missing files so we don't check for them again and again */
	grib_trie_insert(c->def_files,basename,&grib_file_not_found);
	/*grib_context_log(c,GRIB_LOG_ERROR,"Def file \"%s\" not found",basename);*/
	GRIB_MUTEX_UNLOCK(&(c->mutex));
	full[0]=0;
	return NULL;
}

void grib_context_free(const grib_context* c, void* p){
	if (!c) c=grib_context_get_default();
	if(p)  c->free_mem(c,p);
}

void grib_context_free_persistent(const grib_context* c, void* p){
	if (!c) c=grib_context_get_default();
	if(p)  c->free_persistent_mem(c,p);
}

void grib_context_reset(grib_context* c){
	if (!c) c=grib_context_get_default();

	if(c->grib_reader)
	{
		grib_action_file *fr = c->grib_reader->first;
		grib_action_file *fn = fr;
		grib_action* a;

		while(fn){
			fr = fn;
			fn = fn->next;

			a = fr->root;
			while(a)
			{
				grib_action *na = a->next;
				grib_free_action(c, a);
				a = na;
			}
			grib_context_free_persistent(c, fr->filename);
			grib_context_free_persistent(c, fr);
		}
		grib_context_free_persistent(c, c->grib_reader);

	}

	c->grib_reader = NULL;

	if(c->codetable) grib_codetable_delete(c);
	c->codetable = NULL;

	if(c->grib_definition_files_dir)
		grib_context_free(c,c->grib_definition_files_dir);

	if(c->multi_support_on)
		grib_multi_support_reset(c);

}

void grib_context_delete( grib_context* c){
	if (!c) c=grib_context_get_default();
	
	grib_hash_keys_delete( c->keys);
	grib_trie_delete(c->def_files);

	grib_context_reset( c );
	if(c != &default_grib_context)
		grib_context_free_persistent(&default_grib_context,c);
}

void* grib_context_malloc_persistent(const grib_context* c, size_t size){
	void* p =  c->alloc_persistent_mem(c,size);
	if(!p) {
		grib_context_log(c,GRIB_LOG_FATAL,"grib_context_malloc: error allocating %lu bytes",(unsigned long)size);
		Assert(1);
	}
	return p;
}

char* grib_context_strdup_persistent(const grib_context* c,const char* s){
	char *dup = (char*)grib_context_malloc_persistent(c,(strlen(s)*sizeof(char))+1);
	if(dup) strcpy(dup,s);
	return dup;
}

void* grib_context_malloc_clear_persistent(const grib_context* c, size_t size){
	void *p = grib_context_malloc_persistent(c,size);
	if(p) memset(p,0,size);
	return p;
}


void* grib_context_malloc(const grib_context* c, size_t size){
	void* p = NULL;
	if (!c) c=grib_context_get_default();
	if(size == 0) return p;
	else p=c->alloc_mem(c,size);
	if(!p) {
		grib_context_log(c,GRIB_LOG_FATAL,"grib_context_malloc: error allocating %lu bytes",(unsigned long)size);
		Assert(1);
	}
	return p;
}

void* grib_context_realloc(const grib_context* c, void *p,size_t size){

	void* q;
	if (!c) c=grib_context_get_default();
	q=c->realloc_mem(c,p,size);

	if(!q) {
		grib_context_log(c,GRIB_LOG_FATAL,"grib_context_realloc: error allocating %lu bytes",(unsigned long)size);
		exit(1);
	}
	return q;
}

char* grib_context_strdup(const grib_context* c,const char* s){
	char *dup = (char*)grib_context_malloc(c,(strlen(s)*sizeof(char))+1);
	if(dup) strcpy(dup,s);
	return dup;
}

void* grib_context_malloc_clear(const grib_context* c, size_t size){
	void *p = grib_context_malloc(c,size);
	if(p) memset(p,0,size);
	return p;
}

void* grib_context_buffer_malloc(const grib_context* c, size_t size){
	void* p = NULL;
	if (!c) c=grib_context_get_default();
	if(size == 0) return p;
	else p=c->alloc_buffer_mem(c,size);
	if(!p) {
		grib_context_log(c,GRIB_LOG_FATAL,"grib_context_buffer_malloc: error allocating %lu bytes",(unsigned long)size);
		exit(1);
	}
	return p;
}

void grib_context_buffer_free(const grib_context* c, void* p){
	if (!c) c=grib_context_get_default();
	if(p)  c->free_buffer_mem(c,p);
}

void* grib_context_buffer_realloc(const grib_context* c, void *p,size_t size){

	void* q=c->realloc_buffer_mem(c,p,size);

	if(!q) {
		grib_context_log(c,GRIB_LOG_FATAL,"grib_context_buffer_realloc: error allocating %lu bytes",(unsigned long)size);
		exit(1);
	}
	return q;
}

void* grib_context_buffer_malloc_clear(const grib_context* c, size_t size){
	void *p = grib_context_buffer_malloc(c,size);
	if(p) memset(p,0,size);
	return p;
}

void grib_context_set_memory_proc(grib_context* c, grib_malloc_proc m, grib_free_proc f,grib_realloc_proc r)
{
	c->free_mem = f;
	c->alloc_mem = m;
	c->realloc_mem = r;
}

void grib_context_set_persistent_memory_proc(grib_context* c, grib_malloc_proc m, grib_free_proc f)
{
	c->free_persistent_mem = f;
	c->alloc_persistent_mem = m;
}

void grib_context_set_buffer_memory_proc(grib_context* c, grib_malloc_proc m, grib_free_proc f,grib_realloc_proc r)
{
	c->free_buffer_mem = f;
	c->alloc_buffer_mem = m;
	c->realloc_buffer_mem = r;
}


void grib_context_set_data_accessing_proc(grib_context* c, grib_data_read_proc read, grib_data_write_proc write, grib_data_tell_proc             tell)
{
	c->read  = read;
	c->write = write;
	c->tell  = tell;

}

/*              logging procedure                    */
void grib_context_log(const grib_context *c,int level, const char* fmt, ...)
{
	char msg[1024];
	va_list list;

	/* Save some CPU */
	if( level == GRIB_LOG_DEBUG && c->debug<1 )
		return;

	va_start(list,fmt);
	vsprintf(msg, fmt, list);
	va_end(list);

	if(level & GRIB_LOG_PERROR)
	{
		level = level & ~GRIB_LOG_PERROR;

		/* #if HAS_STRERROR */
#if 1
		strcat(msg," (");
		strcat(msg,strerror(errno));
		strcat(msg,")");
#else
		if(errno > 0 && errno < sys_nerr)
		{
			strcat(msg," (");
			strcat(msg,sys_errlist[errno]);
			strcat(msg," )");
		}
#endif
	}


	if(c->output_log)
		c->output_log(c,level,msg);
}

/*              logging procedure                    */
void grib_context_print(const grib_context *c, void* descriptor,const char* fmt, ...)
{
	char msg[1024];
	va_list list;
	va_start(list,fmt);
	vsprintf(msg, fmt, list);
	va_end(list);
	c->print(c,descriptor,msg);
}

