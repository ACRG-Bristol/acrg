/**
 * Copyright 2005-2007 ECMWF
 *
 * Licensed under the GNU Lesser General Public License which
 * incorporates the terms and conditions of version 3 of the GNU
 * General Public License.
 * See LICENSE and gpl-3.0.txt for details.
 */

#include <assert.h>
#include <stdio.h>
#include "grib_api_internal.h"


#define  GRIB 0x47524942
#define  BUDG 0x42554447
#define  TIDE 0x54494445
#define  BUFR 0x42554652

#define GRIB_API_READS_BUFR 1

typedef int   (*readproc)(void*,void*,int,int*);
typedef int   (*seekproc)(void*,off_t);
typedef off_t   (*tellproc)(void*);
typedef void* (*allocproc)(void*,size_t*,int*);


typedef struct reader {
	void *read_data;
	readproc read;

	void *alloc_data;
	allocproc alloc;
	int headers_only;

	seekproc seek;
	tellproc tell;
	off_t offset;

} reader;


static int read_the_rest(reader* r,size_t message_length,unsigned char* tmp, int already_read)
{
	int err            = 0;
	size_t buffer_size = message_length;
	size_t rest=message_length-already_read;

	unsigned char* buffer = (unsigned char*)r->alloc(r->alloc_data,&buffer_size,&err);
	if(err) return err;

	if(buffer_size < message_length)
	{
		if(buffer_size < already_read)
		{
			memcpy(buffer,tmp,buffer_size);
			return GRIB_BUFFER_TOO_SMALL;
		}

		memcpy(buffer,tmp,already_read);
		if((r->read(r->read_data,buffer+already_read,buffer_size-already_read,&err) != buffer_size-already_read) || err)
			return err;

		return GRIB_BUFFER_TOO_SMALL;
	}

	memcpy(buffer,tmp,already_read);

	if((r->read(r->read_data,buffer+already_read,rest,&err) != rest) || err)
			return err;

	if(!r->headers_only && (buffer[message_length-4] != '7' ||
			buffer[message_length-3] != '7' ||
			buffer[message_length-2] != '7' ||
			buffer[message_length-1] != '7')) {

		return GRIB_WRONG_LENGTH;
	}


	return GRIB_SUCCESS;
}

#define CHECK_TMP_SIZE(a) if(sizeof(tmp)<(a)) { fprintf(stderr,"%s:%d sizeof(tmp)<%s %d<%d\n", __FILE__,__LINE__,#a,(int)sizeof(tmp),(int)(a)); return GRIB_INTERNAL_ARRAY_TOO_SMALL; }

#define UINT3(a,b,c) (size_t)((a<<16) + (b<<8) + c);

static int read_GRIB(reader* r)
{
	/* unsigned char tmp[16368];  */
	unsigned char tmp[700368]; /* Should be enough */
	size_t length = 0;
	size_t total_length = 0;
	long edition = 0;
	int  err = 0;
	int i = 0 ,j;
	size_t sec1len = 0;
	size_t sec2len = 0;
	size_t sec3len = 0;
	size_t sec4len = 0;
	unsigned long flags;

	tmp[i++] = 'G';
	tmp[i++] = 'R';
	tmp[i++] = 'I';
	tmp[i++] = 'B';

	r->offset=r->tell(r->read_data)-4;

	if(r->read(r->read_data,&tmp[i],3,&err) != 3 || err)
		return err;

	length= UINT3(tmp[i],tmp[i+1],tmp[i+2]);
	i+=3;

	/* Edition number */
	if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
		return err;

	edition = tmp[i++];

	switch(edition)
	{
		case 1:
			if (r->headers_only) {
				/* Read section 1 length */
				if(r->read(r->read_data,&tmp[i],3,&err) != 3 || err)
					return err;

				sec1len=UINT3(tmp[i],tmp[i+1],tmp[i+2]);
				i+=3;
				/* Read section 1. 3 = length */
				if((r->read(r->read_data,tmp+i,sec1len-3,&err) != sec1len-3) || err)
					return err;
				flags = tmp[15];

				i += sec1len-3;

				CHECK_TMP_SIZE(8+ sec1len +  4 + 3 );

				if(flags & (1<<7)) {
					/* Section 2 */
					if(r->read(r->read_data,&tmp[i],3,&err) != 3 || err)
						return err;

					sec2len=UINT3(tmp[i],tmp[i+1],tmp[i+2]);
					i+=3;
					/* Read section 2 */
					if((r->read(r->read_data,tmp+i,sec2len-3,&err) != sec2len-3) || err)
						return err;
					i += sec2len-3;
				}

				CHECK_TMP_SIZE(8+sec1len +  sec2len + 4 + 3 );

				total_length=length;
				length=8+sec1len + sec2len;

			}
			else if(length & 0x800000)
			{

				/* Large GRIBs */

				/* Read section 1 length */
				for(j=0;j<3;j++)
				{
					if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
						return err;

					sec1len <<= 8;
					sec1len |= tmp[i];
					i++;
				}

				/* table version */
				if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
				/* center */
				if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
				/* process */
				if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
				/* grid */
				if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
				/* flags */
				if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err) return err;
				flags = tmp[i++];

				/* fprintf(stderr," sec1len=%d i=%d flags=%x\n",sec1len,i,flags); */

				CHECK_TMP_SIZE(8+sec1len +  4 + 3 );

				/* Read section 1. 3 = length, 5 = table,center,process,grid,flags */
				if((r->read(r->read_data,tmp+i,sec1len-3-5,&err) != sec1len-3-5) || err)
					return err;

				i += sec1len-3-5;

				if(flags & (1<<7)) {
					/* Section 2 */
					for(j=0;j<3;j++)
					{
						if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
							return err;

						sec2len <<= 8;
						sec2len |= tmp[i];
						i++;
					}
					/* Read section 2 */
					if((r->read(r->read_data,tmp+i,sec2len-3,&err) != sec2len-3) || err)
						return err;
					i += sec2len-3;
				}

				CHECK_TMP_SIZE(sec1len +  sec2len + 4 + 3 );

				if(flags & (1<<6)) {

					/* Section 3 */
					for(j=0;j<3;j++)
					{
						if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
							return err;

						sec3len <<= 8;
						sec3len |= tmp[i];
						i++;
					}

					/* Read section 3 */
					if((r->read(r->read_data,tmp+i,sec3len-3,&err) != sec3len-3) || err)
						return err;
					i += sec3len-3;
				}

				/* fprintf(stderr,"%s sec1len=%d i=%d\n",type,sec1len,i); */

				CHECK_TMP_SIZE(sec1len + sec2len + sec3len + 4 + 3 );


				for(j=0;j<3;j++)
				{
					if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
						return err;

					sec4len <<= 8;
					sec4len |= tmp[i];
					i++;
				}

				if(sec4len < 120)
				{
					/* Special coding */
					length &= 0x7fffff;
					length *= 120;
					length -= sec4len;
					length += 4;
				}
				else
				{
					/* length is already set to the right value */
				}

			}
			break;

		case 2:
			length = 0;

			if(sizeof(long) >= 8) {
				for(j=0;j<8;j++)
				{
					if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
						return err;

					length <<= 8;
					length |= tmp[i];
					i++;
				}
			}
			else
			{
				/* Check if the length fits in a long */
				for(j=0;j<4;j++)
				{
					if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
						return err;

					length <<= 8;
					length |= tmp[i];
					i++;
				}

				if(length)
					return GRIB_MESSAGE_TOO_LARGE; /* Message too large */

				for(j=0;j<4;j++)
				{
					if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
						return err;

					length <<= 8;
					length |= tmp[i];
					i++;
				}
			}
			break;

		default:
			/* fprintf(stderr,"GRIB edition is %d len=%d\n",(int)edition,length); */
			return GRIB_NOT_IMPLEMENTED;
			break;
	}

	Assert(i <= sizeof(tmp));
	err=read_the_rest(r,length,tmp,i);
	if (r->headers_only && edition==1) { 
		 err=r->seek(r->read_data,total_length-length);
	 }

	return err;

}

static int read_TIDE(reader *r,const char* type)
{
	unsigned char tmp[32]; /* Should be enough */
	size_t sec1len = 0;
	size_t sec4len = 0;
	int  err = 0;
	int i = 0, j = 0;

	for(j = 0; j < 4; j++)
	{
		tmp[i] = type[i];
		i++;
	}


	for(j=0;j<3;j++)
	{
		if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
			return err;

		sec1len <<= 8;
		sec1len |= tmp[i];
		i++;
	}

	/* fprintf(stderr,"%s sec1len=%d i=%d\n",type,sec1len,i); */

	CHECK_TMP_SIZE(sec1len + 4 + 3 );


	/* Read sectoin1 */
	if((r->read(r->read_data,tmp+i,sec1len-3,&err) != sec1len-3) || err)
		return err;

	i += sec1len-3;

	for(j=0;j<3;j++)
	{
		if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
			return err;

		sec4len <<= 8;
		sec4len |= tmp[i];
		i++;
	}

	/* fprintf(stderr,"%s sec4len=%d i=%d l=%d\n",type,sec4len,i,4+sec1len+sec4len+4); */

	Assert(i <= sizeof(tmp));
	return read_the_rest(r,4+sec1len+sec4len+4,tmp,i);

}

static int read_BUFR(reader *r)
{
	unsigned char tmp[65536]; /* Should be enough */
	size_t length = 0;
	long edition = 0;
	int  err = 0;
	int i = 0 ,j;

	tmp[i++] = 'B';
	tmp[i++] = 'U';
	tmp[i++] = 'F';
	tmp[i++] = 'R';

	for(j=0;j<3;j++)
	{
		if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
			return err;

		length <<= 8;
		length |= tmp[i];
		i++;
	}

	/* Edition number */
	if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
		return err;

	edition = tmp[i++];

	/* assert(edition != 1); */

	if(edition<2)
	{
		int n;
		size_t sec1len = 0;
		size_t sec2len = 0;
		size_t sec3len = 0;
		size_t sec4len = 0;
		unsigned long flags;

		sec1len = length;

		/* table version */
		if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
		/* center */
		if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
		/* update */
		if(r->read(r->read_data,&tmp[i++],1,&err) != 1 || err) return err;
		/* flags */
		if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err) return err;
		flags = tmp[i++];


		CHECK_TMP_SIZE(sec1len +  4 + 3 );

		/* Read section 1. 3 = length, 5 = table,center,process,flags */

		n = sec1len - 8; /* Just a guess */
		if((r->read(r->read_data,tmp+i,n,&err) != n) || err)
			return err;

		i += n;

		if(flags & (1<<7)) {
			/* Section 2 */
			for(j=0;j<3;j++)
			{
				if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
					return err;

				sec2len <<= 8;
				sec2len |= tmp[i];
				i++;
			}

			/* Read section 2 */
			if((r->read(r->read_data,tmp+i,sec2len-3,&err) != sec2len-3) || err)
				return err;
			i += sec2len-3;
		}

		CHECK_TMP_SIZE(sec1len +  sec2len + 4 + 3 );


		/* Section 3 */
		for(j=0;j<3;j++)
		{
			if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
				return err;

			sec3len <<= 8;
			sec3len |= tmp[i];
			i++;
		}

		/* Read section 3 */
		if((r->read(r->read_data,tmp+i,sec3len-3,&err) != sec3len-3) || err)
			return err;
		i += sec3len-3;


		CHECK_TMP_SIZE( sec1len + sec2len + sec3len + 4 + 3 );


		for(j=0;j<3;j++)
		{
			if(r->read(r->read_data,&tmp[i],1,&err) != 1 || err)
				return err;

			sec4len <<= 8;
			sec4len |= tmp[i];
			i++;
		}

		/* fprintf(stderr," sec1len=%d sec2len=%d sec3len=%d sec4len=%d\n",sec1len, sec2len,sec3len,sec4len); */
		length = 4 + sec1len + sec2len + sec3len + sec4len + 4;
		/* fprintf(stderr,"length = %d i = %d\n",length,i); */
	}

	Assert(i <= sizeof(tmp));
	return read_the_rest(r,length,tmp,i);
}

static int read_any(reader *r,int grib_ok,int bufr_ok)
{
	unsigned char c;
	int err = 0;
	unsigned long magic = 0;

	while(r->read(r->read_data,&c,1,&err) == 1 && err == 0)
	{
		magic <<= 8;
		magic |= c;

		switch(magic & 0xffffffff)
		{

			case GRIB:
				if(grib_ok) 
				{
					err =  read_GRIB(r);
					return err == GRIB_END_OF_FILE ? GRIB_PREMATURE_END_OF_FILE : err; /* Premature EOF */
				}
				break;

			case BUFR:
				if(bufr_ok)
				{
					err =  read_BUFR(r);
					return err == GRIB_END_OF_FILE ? GRIB_PREMATURE_END_OF_FILE : err; /* Premature EOF */
				}
				break;

			case BUDG:
			case TIDE:
				if(grib_ok) 
				{
					err =  read_TIDE(r,magic == TIDE ? "TIDE" : "BUDG");
					return err == GRIB_END_OF_FILE ? GRIB_PREMATURE_END_OF_FILE : err; /* Premature EOF */
				}
				break;

		}
	}

	return err;

}

GRIB_INLINE off_t stdio_tell(void* data) {
	FILE* f = (FILE*)data;
	return ftello(f);
}

GRIB_INLINE int stdio_seek(void* data,off_t len) {
	FILE* f = (FILE*)data;
	int err=0;
	if (fseeko(f,len,SEEK_CUR)) err=GRIB_IO_PROBLEM;
	return err;
}

GRIB_INLINE int stdio_read(void* data,void* buf,int len,int* err)
{
	FILE* f = (FILE*)data;
	int n;
	/* char iobuf[1024*1024]; */

	if (len==0) return 0;

	/* setvbuf(f,iobuf,_IOFBF,sizeof(iobuf)); */
	n   = fread(buf,1,len,f);
	/* fprintf(stderr,"read %d = %x %c\n",1,(int)buf[0],buf[0]); */
	if(n != len) {
		/* fprintf(stderr,"Failed to read %d, only got %d\n",len,n); */
		*err               = GRIB_IO_PROBLEM;
		if(feof(f))   *err = GRIB_END_OF_FILE;
		if(ferror(f)) *err = GRIB_IO_PROBLEM;
	}
	return n;
}


/*================== */
typedef struct user_buffer {
	void*    user_buffer;
	size_t   buffer_size;
	size_t   message_size;
} user_buffer;

static void* user_provider_buffer(void *data,size_t* length,int *err)
{
	user_buffer *u  = (user_buffer*)data;
	u->message_size = *length;
	*length = u->buffer_size;
	return u->user_buffer;
}

static
int _wmo_read_any_from_file(FILE* f,void* buffer,size_t* len,int grib_ok,int bufr_ok)
{
	int         err;
	user_buffer u; 
	reader      r; 

	u.user_buffer  = buffer;
	u.buffer_size  = *len;
	u.message_size = 0;

	r.read_data    = f;
	r.read         = &stdio_read;
	r.seek		   = &stdio_seek;
	r.tell		   = &stdio_tell;
	r.alloc_data   = &u;
	r.alloc        = &user_provider_buffer;
	r.headers_only = 0;


	err            = read_any(&r,grib_ok,bufr_ok);
	*len           = u.message_size;

	return err;
}

int wmo_read_any_from_file(FILE* f,void* buffer,size_t* len)
{
	return _wmo_read_any_from_file(f,buffer,len,1,1);
}

int wmo_read_grib_from_file(FILE* f,void* buffer,size_t* len)
{
	return _wmo_read_any_from_file(f,buffer,len,1,0);
}

int wmo_read_bufr_from_file(FILE* f,void* buffer,size_t* len)
{
	return _wmo_read_any_from_file(f,buffer,len,0,1);
}

/*================== */

typedef struct stream_struct {

	void* stream_data;
	long (*stream_proc)(void*,void* buffer,long len);

} stream_struct;

static off_t stream_tell(void* data)
{
	return 0;
}

static int stream_seek(void* data,off_t len)
{
	return 0;
}
static int stream_read(void* data,void* buffer,int len,int* err)
{
	stream_struct *s = (stream_struct*)data;
	long n =  s->stream_proc(s->stream_data,buffer,len);
	if(n != len) {
		*err = GRIB_END_OF_FILE;
	}
	return n;
}


int wmo_read_any_from_stream(void* stream_data,long (*stream_proc)(void*,void* buffer,long len) ,void* buffer,size_t* len)
{
	int           err;
	stream_struct s;
	user_buffer   u; 
	reader        r; 

	s.stream_data = stream_data;
	s.stream_proc = stream_proc;

	u.user_buffer  = buffer;
	u.buffer_size  = *len;
	u.message_size = 0;

	r.read_data    = &s;
	r.read         = &stream_read;
	r.seek         = &stream_seek;
	r.tell         = &stream_tell;
	r.alloc_data   = &u;
	r.alloc        = &user_provider_buffer;
	r.headers_only = 0;

	err            = read_any(&r,1,1);
	*len           = u.message_size;

	return err;
}

/*================== */

typedef struct alloc_buffer {
	void* buffer;
} alloc_buffer;

static void* allocate_buffer(void *data,size_t* length,int *err)
{
	alloc_buffer *u  = (alloc_buffer*)data;
	u->buffer = malloc(*length);
	if(u->buffer == NULL)
		*err = GRIB_OUT_OF_MEMORY; /* Cannot allocate buffer */
	return u->buffer;
}


static
void *_wmo_read_any_from_file_malloc(FILE* f,int* err,int grib_ok,int bufr_ok)
{
	alloc_buffer u; 
	reader       r; 

	u.buffer       = NULL;

	r.read_data    = f;
	r.read         = &stdio_read;
	r.seek         = &stdio_seek;
	r.tell         = &stdio_tell;
	r.alloc_data   = &u;
	r.alloc        = &allocate_buffer;
	r.headers_only = 0;
  
	*err           = read_any(&r,grib_ok,bufr_ok);

	return u.buffer;
}

void *wmo_read_any_from_file_malloc(FILE* f,int* err)
{
	return _wmo_read_any_from_file_malloc(f,err,1,1);
}

void *wmo_read_grib_from_file_malloc(FILE* f,int* err)
{
	return _wmo_read_any_from_file_malloc(f,err,1,0);
}

void *wmo_read_bufr_from_file_malloc(FILE* f,int* err)
{
	return _wmo_read_any_from_file_malloc(f,err,0,1);
}

/* ======================================= */

typedef struct context_alloc_buffer {
	grib_context*  ctx;
	void*          buffer;
	size_t         length;
} context_alloc_buffer;

static void* context_allocate_buffer(void *data,size_t* length,int *err)
{
	context_alloc_buffer *u  = (context_alloc_buffer*)data;
	u->buffer = grib_context_malloc(u->ctx,*length);
	u->length = *length;

	if(u->buffer == NULL)
		*err = GRIB_OUT_OF_MEMORY; /* Cannot allocate buffer */
	return u->buffer;
}


int grib_read_any_headers_only_from_file_alloc(grib_context* ctx,FILE* f,void **buffer,size_t* length,off_t* offset)
{
	int err;
	context_alloc_buffer u; 
	reader       r; 

	u.buffer       = NULL;
	u.length       = 0;
	u.ctx          = ctx ? ctx : grib_context_get_default();

	r.read_data    = f;
	r.read         = &stdio_read;
	r.seek         = &stdio_seek;
	r.tell         = &stdio_tell;
	r.alloc_data   = &u;
	r.alloc        = &context_allocate_buffer;
	r.headers_only = 1;

	err            = read_any(&r,1,GRIB_API_READS_BUFR);

	*buffer        = u.buffer;
	*length        = u.length;
	*offset		   = r.offset;

	return err;
}

int grib_read_any_from_file_alloc(grib_context* ctx,FILE* f,void **buffer,size_t* length)
{
	int err;
	context_alloc_buffer u; 
	reader       r; 

	u.buffer       = NULL;
	u.length       = 0;
	u.ctx          = ctx ? ctx : grib_context_get_default();

	r.read_data    = f;
	r.read         = &stdio_read;
	r.seek         = &stdio_seek;
	r.tell         = &stdio_tell;
	r.alloc_data   = &u;
	r.alloc        = &context_allocate_buffer;
	r.headers_only = 0;

	err            = read_any(&r,1,GRIB_API_READS_BUFR);

	*buffer        = u.buffer;
	*length        = u.length;

	return err;
}


int grib_read_any_headers_only_from_file(grib_context* ctx,FILE* f,void* buffer,size_t* len)
{
	int         err;
	user_buffer u; 
	reader      r; 

	u.user_buffer  = buffer;
	u.buffer_size  = *len;
	u.message_size = 0;

	r.read_data    = f;
	r.read         = &stdio_read;
	r.seek         = &stdio_seek;
	r.tell         = &stdio_tell;
	r.alloc_data   = &u;
	r.alloc        = &user_provider_buffer;
	r.headers_only = 1;

	err            = read_any(&r,1,GRIB_API_READS_BUFR);

	*len           = u.message_size;

	return err;
}

int grib_read_any_from_file(grib_context* ctx,FILE* f,void* buffer,size_t* len)
{
	int         err;
	user_buffer u; 
	reader      r; 

	u.user_buffer  = buffer;
	u.buffer_size  = *len;
	u.message_size = 0;

	r.read_data    = f;
	r.read         = &stdio_read;
	r.seek		   = &stdio_seek;
	r.tell		   = &stdio_tell;
	r.alloc_data   = &u;
	r.alloc        = &user_provider_buffer;
	r.headers_only = 0;

	err            = read_any(&r,1,GRIB_API_READS_BUFR);

	*len           = u.message_size;

	return err;
}

/* ======================================= */

typedef struct memory_read_data {
	unsigned char  *data;
	size_t          data_len;
} memory_read_data;

static off_t memory_tell(void* data)
{
	return 0;
}

static int memory_seek(void* data,off_t len)
{
	return 0;
}

static int memory_read(void* data,void* buf,int len,int* err)
{
	memory_read_data *m = (memory_read_data*)data;

	if(len == 0)
	{
		*err = GRIB_END_OF_FILE;
		return 0;
	}
	else {
		size_t l = len > m->data_len ? m->data_len : len;
		memcpy(buf,m->data,l);
		m->data_len -= l;
		m->data     += l;
		return l;
	}
}

int grib_read_any_from_memory_alloc(grib_context* ctx,unsigned char** data,size_t* data_length,void **buffer,size_t* length)
{
	int err;
	memory_read_data     m;
	context_alloc_buffer u; 
	reader       r; 

	m.data         = *data;
	m.data_len     = *data_length;

	u.buffer       = NULL;
	u.length       = 0;
	u.ctx          = ctx ? ctx : grib_context_get_default();

	r.read_data    = &m;
	r.read         = &memory_read;
	r.seek		   = &memory_seek;
	r.tell		   = &memory_tell;
	r.alloc_data   = &u;
	r.alloc        = &context_allocate_buffer;
	r.headers_only = 0;

	err            = read_any(&r,1,GRIB_API_READS_BUFR);
	*buffer        = u.buffer;
	*length        = u.length;

	*data_length   = m.data_len;
	*data          = m.data;

	return err;
}


int grib_read_any_from_memory(grib_context* ctx,unsigned char** data,size_t* data_length,void* buffer,size_t* len)
{
	int         err;
	memory_read_data     m;
	user_buffer u; 
	reader      r; 

	m.data         = *data;
	m.data_len     = *data_length;


	u.user_buffer  = buffer;
	u.buffer_size  = *len;
	u.message_size = 0;

	r.read_data    = &m;
	r.read         = &memory_read;
	r.seek		   = &memory_seek;
	r.tell		   = &memory_tell;
	r.alloc_data   = &u;
	r.alloc        = &user_provider_buffer;
	r.headers_only = 0;

	err            = read_any(&r,1,GRIB_API_READS_BUFR);
	*len           = u.message_size;

	*data_length   = m.data_len;
	*data          = m.data;


	return err;
}

int grib_count_in_file(grib_context* c, FILE* f,int* n) {
	grib_handle* h;
	int error=0;
	*n=0;
	while ((h=grib_handle_headers_only_new_from_file(c,f,&error))!=NULL) {
		(*n)++;
		grib_handle_delete(h);
	}
	rewind(f);
	return error==GRIB_END_OF_FILE ? 0 : error ;
}



