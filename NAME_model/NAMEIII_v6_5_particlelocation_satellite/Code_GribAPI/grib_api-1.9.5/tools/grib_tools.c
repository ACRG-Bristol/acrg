/**
 * Copyright 2005-2007 ECMWF
 *
 * Licensed under the GNU Lesser General Public License which
 * incorporates the terms and conditions of version 3 of the GNU
 * General Public License.
 * See LICENSE and gpl-3.0.txt for details.
 */

/*
 * C Implementation: grib_tools
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */

#include "grib_tools.h"
#if HAVE_LIBJASPER
#include "jasper/jasper.h"
#endif

GRIB_INLINE static int grib_inline_strcmp(const char* a,const char* b) {
	if (*a != *b) return 1;
	while((*a!=0 && *b!=0) &&  *(a) == *(b) ) {a++;b++;}
	return (*a==0 && *b==0) ? 0 : 1;
}

static void grib_print_header(grib_runtime_options* options,grib_handle* h);
static void grib_tools_set_print_keys(grib_runtime_options* options,grib_handle* h,const char* ns);
static int grib_tool_with_orderby(grib_runtime_options* options);
static int grib_tool_without_orderby(grib_runtime_options* options);
static int grib_tool_onlyfiles(grib_runtime_options* options);
static int grib_tool_index(grib_runtime_options* options);
static int process(grib_context* c,grib_runtime_options* options,const char* path);
static int scan(grib_context* c,grib_runtime_options* options,const char* dir);

grib_runtime_options options={
		0,         /* verbose       */
		0,         /* fail          */
		0,         /* skip          */
		12,        /* default_print_width */
		0,         /* print_header */
		0,         /* name_space */
		0,         /* print_number */
		1,         /* print_statistics */
        {{0,},},   /* grib_values requested_print_keys[MAX_KEYS] */
		0,         /* requested_print_keys_count */
		{{0,},},   /* grib_values print_keys[MAX_KEYS] */
		0,         /* print_keys_count  */
		0,         /* strict            */
		0,         /* multi_support     */
		0,         /* set_values_count  */
		{{0,},},   /* grib_values set_values[MAX_KEYS] */
		{{0,},},   /* grib_values constraints[MAX_KEYS] */
		0,         /* constraints_count */
		{{0,},},   /* grib_values compare[MAX_KEYS] */
		0,         /* compare_count */
		0,         /* handle_count      */
		0,         /* filter_handle_count */
		0,         /* file_count     */
		0,         /* grib_tools_file infile_extra */
		0,         /* grib_tools_file current_infile */
		0,         /* grib_tools_file infile */
		0,         /*grib_tools_file outfile */
		0,         /* grib_action action */
		0,         /* grib_rule rules */
		0,         /* int dump_flags; */
		0,         /* char* dump_mode; */
		0,         /* repack    */
		0,         /* error    */
		0,          /* gts    */
		0,          /* orderby    */
		0,          /* latlon    */
		{0,},
		{0,},
		{0,},
		{0,},
		{0,},
		4,
		0,
		-1,
		{0,},
		0,       /* index */
		0,       /* index_on */
		0,        /* constant */
		0,         /* dump_filename*/
		0,         /* index */
		0,         /* random */
        0,         /* format */
        0,         /* onlyfiles */
        0,         /* tolerance_count  */
		0,			/* through_index */
        0,         /* index1  */
        0,         /* index2  */
        0,         /* context  */
        0,         /* stop  */
        0,         /* headers_only  */
        {{0,},}    /* grib_values tolerance[MAX_KEYS] */

	};

static int is_index_file(const char* filename) {
	FILE* fh;
	char buf[8]={0,};
	char* str="GRBIDX";
	int ret=0;

	fh=fopen(filename,"r");
	if (!fh) return 0;

	fread(buf,1,1,fh);
	fread(buf,6,1,fh);

	ret=!strcmp(buf,str);

	fclose(fh);

	return ret;
}

static grib_handle* grib_handle_new_from_file_x(grib_context* c,FILE* f,int headers_only,int *err) {
	if (headers_only) 
		return grib_handle_headers_only_new_from_file(c,f,err);

	return grib_handle_new_from_file(c,f,err);
}

int grib_tool(int argc, char **argv)
{
	int ret=0;
	grib_context* c=grib_context_get_default();
	options.context=c;

	if (getenv("DOXYGEN_USAGE") && argc==1 ) usage_doxygen();

	grib_get_runtime_options(argc,argv,&options);

	grib_tool_before_getopt(&options);

	grib_process_runtime_options(c,argc,argv,&options);

	grib_tool_init(&options);
	if (options.dump_filename) {
		dump_file= fopen(options.dump_filename,"w");
		if(!dump_file) {
			perror(options.dump_filename);
			exit(1);
		}
	} else {
		dump_file=stdout;
	}

	if (is_index_file(options.infile->name) && 
		( options.infile_extra && is_index_file(options.infile_extra->name))) { 
			options.through_index=1;
			return grib_tool_index(&options);
		}

	if (options.onlyfiles) 
		ret=grib_tool_onlyfiles(&options);
	else {
		if (options.orderby)
			ret=grib_tool_with_orderby(&options);
		else
			ret=grib_tool_without_orderby(&options);
	}

	if (options.dump_filename) fclose(dump_file);
	return ret;

}

static int grib_tool_with_orderby(grib_runtime_options* options) {
	int err=0;
	grib_failed *failed=NULL,*p=NULL;
	grib_handle* h=NULL;
	grib_context* c=NULL;
	grib_tools_file* infile=NULL;
	char** filenames;
	size_t files_count=0;
	grib_fieldset* set=NULL;
	int i=0;
    c=grib_context_get_default();

	infile=options->infile;
    if(infile) infile->failed=NULL;
    
	files_count=0;
	while(infile) {files_count++;infile=infile->next;}

	filenames=(char**)grib_context_malloc_clear(c,files_count*sizeof(char*));

	infile=options->infile;
	for (i=0;i<files_count;i++) {filenames[i]=infile->name;infile=infile->next;}

	if (grib_options_on("7")) c->no_fail_on_wrong_length=1;

	set=grib_fieldset_new_from_files(0,filenames,files_count,0,0,0,options->orderby,&err);
	if (err) {
		grib_context_log(c,GRIB_LOG_FATAL,"unable to create index for input file %s (%s)",
				filenames[0],grib_get_error_message(err));
		exit(err);
	}

	options->handle_count=0;
	while((h = grib_fieldset_next_handle(set,&err))
			!= NULL || err != GRIB_SUCCESS ) {
		options->handle_count++;
		options->error=err;

		if (!h) {
			fprintf(dump_file,"\t\t*** unable to read message ***\n");
			if (options->fail || err==GRIB_WRONG_LENGTH) GRIB_CHECK_NOLINE(err,0);

			failed=(grib_failed*)grib_context_malloc_clear(c,sizeof(grib_failed));
			failed->count=infile->handle_count;
			failed->error=err;
			failed->next=NULL;

			if (!infile->failed) {
				infile->failed=failed;
			} else {
				p=infile->failed;
				while (p->next) p=p->next;
				p->next=failed;
			}

			continue;
		}

		grib_print_header(options,h);

		grib_skip_check(options,h);

		if (options->skip && options->strict) {
			grib_tool_skip_handle(options,h);
			continue;
		}

		grib_tool_new_handle_action(options,h);

		grib_tool_print_key_values(options,h);

		grib_handle_delete(h);

	}


	grib_tool_finalise_action(options);

	return 0;
}

char iobuf[1024*1024];

static int grib_tool_without_orderby(grib_runtime_options* options) {
	int err=0;
	grib_failed *failed=NULL,*p=NULL;
	grib_handle* h=NULL;
	grib_context* c=NULL;
	grib_tools_file* infile=NULL;

	c=grib_context_get_default();
	options->file_count=0;
	options->handle_count=0;
	options->filter_handle_count=0;
	options->current_infile=options->infile;
	infile=options->infile;
    infile->failed=NULL;

	if (grib_options_on("7")) c->no_fail_on_wrong_length=1;

	while (infile!=NULL && infile->name!=NULL) {

		if (options->print_statistics && options->verbose) fprintf(dump_file,"%s\n",infile->name);
		infile->file = fopen(infile->name,"r");
		if(!infile->file) {
			perror(infile->name);
			exit(1);
		}

		setvbuf(infile->file,iobuf,_IOFBF,sizeof(iobuf));


		options->file_count++;
		infile->handle_count=0;
		infile->filter_handle_count=0;

		grib_tool_new_file_action(options,infile);


		while((h = grib_handle_new_from_file_x(c,infile->file,options->headers_only,&err))
				!= NULL || err != GRIB_SUCCESS ) {
			infile->handle_count++;
			options->handle_count++;
			options->error=err;

			if (!h) {
				fprintf(dump_file,"\t\t*** unreadable message ***\n");
				if (options->fail ) GRIB_CHECK_NOLINE(err,0);

				failed=(grib_failed*)grib_context_malloc_clear(c,sizeof(grib_failed));
				failed->count=infile->handle_count;
				failed->error=err;
				failed->next=NULL;

				if (!infile->failed) {
					infile->failed=failed;
				} else {
					p=infile->failed;
					while (p->next) p=p->next;
					p->next=failed;
				}

				continue;
			}

			grib_print_header(options,h);

			grib_skip_check(options,h);

			if (options->skip && options->strict) {
				grib_tool_skip_handle(options,h);
				continue;
			}

			grib_tool_new_handle_action(options,h);

			grib_print_key_values(options,h);

			grib_handle_delete(h);

		}

		grib_print_file_statistics(options,infile);

		if (infile->file) fclose(infile->file);

		if (infile->handle_count==0) {
			fprintf(dump_file,"no grib messages found in %s\n", infile->name);
			if (options->fail) exit(1);
		}

		infile=infile->next;
		options->current_infile=infile;

	}

	grib_print_full_statistics(options);

	grib_tool_finalise_action(options);

	return 0;
}

static int navigate(grib_field_tree* fields,grib_runtime_options* options) {
	int err=0;
	if (!fields || options->stop) return 0;

	if (fields->field) {
		grib_handle* h=grib_index_get_handle(fields->field,&err);
		if (!options->index2->current) 
			options->index2->current=grib_context_malloc_clear(options->context,sizeof(grib_field_list));
		options->index2->current->field=fields->field;
		if (!h) return err;
		grib_skip_check(options,h);
		if (options->skip && options->strict) {
			grib_tool_skip_handle(options,h);
	    } else {
			grib_tool_new_handle_action(options,h);
			grib_handle_delete(h);
		}
	}

	err=navigate(fields->next_level,options);
	if (err) return err;

	err=navigate(fields->next,options);

	return err;
}


static int grib_tool_index(grib_runtime_options* options) {
	int err=0;
	grib_context* c=NULL;
	char* f1=options->infile->name;
	char* f2=options->infile_extra->name;
	grib_index_key *k1,*k2;
	int found=0;

	c=grib_context_get_default();

	options->index1=grib_index_read(c,f1,&err);
	if (err) 
		grib_context_log(c,(GRIB_LOG_FATAL) | (GRIB_LOG_PERROR) ,
			"unable to read index from %s",f1);
	
	options->index2=grib_index_read(c,f2,&err);
	if (err) 
		grib_context_log(c,(GRIB_LOG_FATAL) | (GRIB_LOG_PERROR) ,
			"unable to read index from %s",f2);

	k1=options->index1->keys;
	while ( k1 ) {
		k2=options->index2->keys;
		found=0;
		while (k2) {
			if ( !strcmp(k1->name,k2->name) ) {
				found=1;
				break;
			}
			k2=k2->next;
		}
		if (!found) {
			printf("Indexes contained in the input files have different keys\n");
			printf("keys in file %s:\n",f1);
			k1=options->index1->keys;
			while (k1) {
				printf("\t%s\n",k1->name);
				k1=k1->next;
			}
			printf("keys in file %s:\n",f2);
			k2=options->index2->keys;
			while (k2) {
				printf("\t%s\n",k2->name);
				k2=k2->next;
			}

			exit(1);

		}

		k1->value[0]=0;
		k1=k1->next;
	}

	k2=options->index2->keys;
	while ( k2 ) {
		k1=options->index1->keys;
		found=0;
		while (k1) {
			if ( !strcmp(k1->name,k2->name) ) {
				found=1;
				break;
			}
			k1=k1->next;
		}
		if (!found) {
			printf("Indexes contained in the input files have different keys\n");
			printf("keys in file %s:\n",f2);
			k2=options->index2->keys;
			while (k2) {
				printf("\t%s\n",k2->name);
				k2=k2->next;
			}
			printf("keys in file %s:\n",f1);
			k1=options->index1->keys;
			while (k1) {
				printf("\t%s\n",k1->name);
				k1=k1->next;
			}

			exit(1);

		}

		k2=k2->next;
	}

	navigate(options->index2->fields,options);

	grib_tool_finalise_action(options);

	return 0;
}

static int scan(grib_context* c,grib_runtime_options* options,const char* dir) {
	struct dirent *s;
	DIR *d;
	int err=0;

	d= opendir(dir);
	if (!d) {
		grib_context_log(c,(GRIB_LOG_ERROR) | (GRIB_LOG_PERROR) , "opendir %s",dir);
		return GRIB_IO_PROBLEM;
	}

	while ((s=readdir(d)) && (err==0)) {
		if(strcmp(s->d_name,".") != 0 && strcmp(s->d_name,"..") != 0) {
			char buf[1024];
			sprintf(buf,"%s/%s",dir,s->d_name);
			process(c,options,buf);
		}
	}
	return 0;
}


static int process(grib_context* c,grib_runtime_options* options,const char* path) {
	struct stat s;
	int err=0;
	int ioerr=0;
	if ( (err = lstat(path,&s)) ) {
		ioerr=errno;
		grib_context_log(c,(GRIB_LOG_ERROR) | (GRIB_LOG_PERROR),"Cannot stat %s",path);
		return GRIB_IO_PROBLEM;
	}

	if (S_ISDIR(s.st_mode) && !S_ISLNK(s.st_mode)) {
		scan(c,options,path);
	} else {
		grib_tool_new_filename_action(options,path);
	}
	return 0;

}

static int grib_tool_onlyfiles(grib_runtime_options* options) {
	grib_context* c=NULL;
	grib_tools_file* infile=NULL;

	c=grib_context_get_default();
	infile=options->infile;

	while (infile!=NULL && infile->name!=NULL) {

		process(c,options,infile->name);

		infile=infile->next;
	}

	grib_tool_finalise_action(options);

	return 0;
}

static void grib_print_header(grib_runtime_options* options,grib_handle* h) {

	size_t strlenkey=0;
	int width;
	if (!options->print_keys || options->handle_count!=1)
		return;

    grib_tools_set_print_keys(options,h,options->name_space);

	if (options->print_keys
			&& options->verbose
			&& options->print_header) {
		int j=0;
		for (j=0;j<options->print_keys_count;j++) {
			strlenkey=strlen(options->print_keys[j].name);
			width= strlenkey < options->default_print_width  ?
				options->default_print_width+2 : strlenkey+2;
            if (options->default_print_width < 0)
				width=strlenkey+1;
			fprintf(dump_file,"%-*s",(int)width,options->print_keys[j].name);
		}
		if (options->latlon) {
			if (options->latlon_mode==4) {
				fprintf(dump_file,"       value1 ");
				fprintf(dump_file," value2 ");
				fprintf(dump_file," value3 ");
				fprintf(dump_file," value4 ");
			} else fprintf(dump_file," value ");
		}
		if (options->index_on) {
			fprintf(dump_file,"        value(%d) ",(int)options->index);
		}
		fprintf(dump_file,"\n");
	}
}


static void grib_tools_set_print_keys(grib_runtime_options* options, grib_handle* h, const char* ns) {
	int i=0;
	grib_keys_iterator* kiter=NULL;

    options->print_keys_count=0;
    
    for (i=0;i<options->requested_print_keys_count;i++) {
      options->print_keys[options->print_keys_count].name=options->requested_print_keys[i].name;
	  if (strlen(options->requested_print_keys[i].name)>options->default_print_width)
		  options->default_print_width=strlen(options->requested_print_keys[i].name);
      options->print_keys[options->print_keys_count].type=options->requested_print_keys[i].type;
      options->print_keys_count++;
    }
    
    if (ns) {
      kiter=grib_keys_iterator_new(h,0,(char*)ns);
      if (!kiter) {
        fprintf(dump_file,"ERROR: Unable to create keys iterator\n");
        exit(1);
      }

      while(grib_keys_iterator_next(kiter))
      {
        const char* name = grib_keys_iterator_get_name(kiter);

        if (options->print_keys_count >= MAX_KEYS ) {
          fprintf(stderr,"ERROR: keys list too long (more than %d keys)\n",
                  options->print_keys_count);
          exit(1);
        }
        options->print_keys[options->print_keys_count].name=strdup(name);
		if (strlen(name)>options->default_print_width)
			options->default_print_width=strlen(name);
        options->print_keys[options->print_keys_count].type=GRIB_TYPE_STRING;
        options->print_keys_count++;
      }

      grib_keys_iterator_delete(kiter);
      if (options->print_keys_count==0 && options->latlon == 0 ) {
        printf("ERROR: name space \"%s\" does not contain any key\n",ns);
        exit(1);
      }

    }

}

static int to_skip(grib_handle* h,grib_values* v,int *err) {
  double dvalue=0;
  int ret=0;
  long lvalue=0;
  char value[MAX_STRING_LEN]={0,};
  size_t len=MAX_STRING_LEN;
  *err=0;
  
  switch (v->type) {
    case GRIB_TYPE_STRING:
      *err=grib_get_string( h,v->name,value,&len);
      ret = v->equal ? grib_inline_strcmp(value,v->string_value) : !grib_inline_strcmp(value,v->string_value);
      break;
    case GRIB_TYPE_DOUBLE:
      *err=grib_get_double( h,v->name,&dvalue);
      ret = v->equal ? (dvalue != v->double_value) : (dvalue == v->double_value);
      break;
    case GRIB_TYPE_LONG:
      *err=grib_get_long( h,v->name,&lvalue);
      ret = v->equal ? (lvalue != v->long_value) : (lvalue == v->long_value);
      break;
    case GRIB_TYPE_MISSING:
      lvalue=grib_is_missing( h,v->name,err);
      ret = (lvalue == v->equal) ? 0 : 1;
      break;
	default:
	  fprintf(dump_file,"invalid type for %s\n",v->name);
	  exit(1);

  }

  return ret;
}

void grib_skip_check(grib_runtime_options* options,grib_handle* h) {
	int i,ret=0;
    grib_values* v=NULL;
    for (i=0;i < options->constraints_count ;i++) {
        v=&(options->constraints[i]);
        if (v->equal) {
          options->skip=1;
          while (v) {
            if (!to_skip(h,v,&ret)) {
              options->skip=0;
              break;
            }
            if (ret != GRIB_SUCCESS && options->fail) {
                grib_context_log(h->context,GRIB_LOG_ERROR,"unable to get \"%s\" (%s)",
                        v->name,grib_get_error_message(ret));
                exit(ret);
            }
            v=v->next;
          }
        } else {
          options->skip=0;
          while (v) {
            if (to_skip(h,v,&ret)) {
              options->skip=1;
              break;
            }
            if (ret != GRIB_SUCCESS && options->fail) {
              grib_context_log(h->context,GRIB_LOG_ERROR,"unable to get \"%s\" (%s)",
                               v->name,grib_get_error_message(ret));
              exit(ret);
            }
            v=v->next;
          }
        }
        if (options->skip==1)
          break;
	}

	if (!options->skip) {
		options->filter_handle_count++;
		if (options->current_infile)
			options->current_infile->filter_handle_count++;
	}
}

void grib_print_key_values(grib_runtime_options* options,grib_handle* h) {
	int i=0;
	int ret=0,width=0;
	size_t strlenvalue=0;
	double dvalue=0;
	long lvalue=0;
	char value[MAX_STRING_LEN];
	char* notfound="not_found";

	if (!options->verbose) return;

	for (i=0;i<options->print_keys_count;i++) {
		size_t len=MAX_STRING_LEN;
		ret=GRIB_SUCCESS;

		if (grib_is_missing(h,options->print_keys[i].name,&ret) && ret==GRIB_SUCCESS)
			sprintf(value,"MISSING");
		else if ( ret == GRIB_SUCCESS ) {
          if (options->print_keys[i].type == GRIB_TYPE_UNDEFINED)
            grib_get_native_type(h,options->print_keys[i].name,&(options->print_keys[i].type));
			switch (options->print_keys[i].type) {
				case GRIB_TYPE_STRING:
					ret=grib_get_string( h,options->print_keys[i].name,value,&len);
					break;
				case GRIB_TYPE_DOUBLE:
					ret=grib_get_double( h,options->print_keys[i].name,&dvalue);
					sprintf(value,options->format,dvalue);
					break;
				case GRIB_TYPE_LONG:
					ret=grib_get_long( h,options->print_keys[i].name,&lvalue);
					sprintf(value,"%d",(int)lvalue);
					break;
				default:
					fprintf(dump_file,"invalid format option for %s\n",options->print_keys[i].name);
					exit(1);
			}
        }
		if (ret != GRIB_SUCCESS) {
			if (options->fail) GRIB_CHECK_NOLINE(ret,options->print_keys[i].name);
			if (ret == GRIB_NOT_FOUND) strcpy(value,notfound);
			else {
				fprintf(dump_file,"%s %s\n",grib_get_error_message(ret),options->print_keys[i].name);
				exit(ret);
			}
		}

        strlenvalue = strlen(value);

        width = strlenvalue < options->default_print_width ?
                    options->default_print_width + 2 :
                    strlenvalue + 2;

        if (options->default_print_width < 0) width = strlenvalue + 1;

        if (options->print_keys_count==i+1 && options->latlon==0) width--;


		fprintf(dump_file,"%-*s",(int)width,value);
	}

	if (options->latlon) {

		if (options->latlon_mode==4){
			int i=0;
            for (i=0;i<4;i++) {
              fprintf(dump_file,options->format,options->values[i]);
              fprintf(dump_file," ");
            }
		} else if (options->latlon_mode==1) {
			sprintf(value,options->format,options->values[options->latlon_idx]);
			strlenvalue = strlen(value);
			width = strlenvalue < options->default_print_width ?
                    options->default_print_width + 2 :
					strlenvalue + 2;
			fprintf(dump_file,"%-*s",(int)width,value);
		}
	}
	if (options->index_on) {
		double v=0;
		/*if (grib_get_double_element(h,"values",options->index,&v) != GRIB_SUCCESS) {*/
		if (1) {
			size_t size;
			double* values;
			int err=0;

			err=grib_get_size(h,"values",&size);
			if (err) {
				sprintf(value,"unknown");
				if (!options->fail) exit(err);
				return;
			}
			values=grib_context_malloc_clear(h->context,size*sizeof(double));
			grib_get_double_array(h,"values",values,&size);
			v=values[options->index];
			grib_context_free(h->context,values);		
		}
		
		sprintf(value,options->format,v);
		strlenvalue = strlen(value);
		width = strlenvalue < options->default_print_width ?
                    options->default_print_width + 2 :
				strlenvalue + 2;
		fprintf(dump_file,"%-*s",(int)width,value);
		
	}
	fprintf(dump_file,"\n");

}


void grib_print_file_statistics(grib_runtime_options* options,grib_tools_file* file) {
	grib_failed* failed=NULL;

	Assert(file);

	failed=file->failed;

	if (!options->print_statistics || !options->verbose) return;

	fprintf(dump_file,"%d of %d grib messages in %s\n\n",
			file->filter_handle_count,
			file->handle_count,
			file->name);
	if (!failed) return;
	/*
	   fprintf(dump_file,"Following bad grib messages found in %s\n",
	   file->name);
	   fprintf(dump_file,"N      Error\n");
	   while (failed){
	   fprintf(dump_file,"%-*d    %s\n",
	   7,failed->count,
	   grib_get_error_message(failed->error));
	   failed=failed->next;
	   }
	   fprintf(dump_file,"\n");
	 */
}

	void grib_print_full_statistics(grib_runtime_options* options) {
		if (options->print_statistics && options->verbose)
			fprintf(dump_file,"%d of %d total grib messages in %d files\n",
					options->filter_handle_count,options->handle_count,options->file_count);
	}

void grib_tools_write_message(grib_runtime_options* options, grib_handle* h) {
	const void *buffer; size_t size;
	grib_file* of=NULL;
	int err=0;
	int ioerr=0;
	char filename[1024]={0,};
	Assert(options->outfile!=NULL && options->outfile->name!=NULL);

	if (options->error == GRIB_WRONG_LENGTH) return;

	if ((err=grib_get_message(h,&buffer,&size))!= GRIB_SUCCESS) {
		grib_context_log(h->context,GRIB_LOG_ERROR,"unable to get binary message\n");
		exit(err);
	}

	err = grib_recompose_name(h,NULL,options->outfile->name,filename,0);

	of=grib_file_open(filename,"w",&err);

	if (!of || !of->handle) {
		ioerr=errno;
		grib_context_log(h->context,(GRIB_LOG_ERROR)|(GRIB_LOG_PERROR),
				"unable to open file %s\n",filename);
		exit(GRIB_IO_PROBLEM);
	}

	if (options->gts && h->gts_header)
		fwrite(h->gts_header,1,h->gts_header_len,of->handle);

	if(fwrite(buffer,1,size,of->handle) != size) {
		ioerr=errno;
		grib_context_log(h->context,(GRIB_LOG_ERROR)|(GRIB_LOG_PERROR),
				"Error writing to %s",filename);
		exit(GRIB_IO_PROBLEM);
	}

	if (options->gts && h->gts_header) {
		char gts_trailer[4]={'\x0D','\x0D','\x0A','\x03'};
		fwrite(gts_trailer,1,4,of->handle);
	}

	grib_file_close(filename,&err);

	if (err != GRIB_SUCCESS) {
		grib_context_log(h->context,GRIB_LOG_ERROR,"unable to write message\n");
		exit(err);
	}

	options->outfile->file=NULL;

#if 0
	if (!options->outfile->file)  {
		options->outfile->file = fopen(options->outfile->name,"w");
		if(!options->outfile->file) {
			perror(options->outfile->name);
			exit(1);
		}
	}
	GRIB_CHECK_NOLINE(grib_get_message(h,&buffer,&size),0);
	if (options->gts && h->gts_header)
		fwrite(h->gts_header,1,h->gts_header_len,options->outfile->file);

	if(fwrite(buffer,1,size,options->outfile->file) != size)
	{
		perror(options->outfile->name);
		exit(1);
	}

	if (options->gts && h->gts_header) {
		char gts_trailer[4]={'\x0D','\x0D','\x0A','\x03'};
		fwrite(gts_trailer,1,4,options->outfile->file);
	}
#endif

}


