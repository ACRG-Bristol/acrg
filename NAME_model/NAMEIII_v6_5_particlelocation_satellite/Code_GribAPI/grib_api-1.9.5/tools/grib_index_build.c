/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: grib_copy
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */

#include "grib_tools.h"

char* grib_tool_description="Build index from a set of files";
char* grib_tool_name="grib_index_build";
char* grib_tool_usage="[options] grib_file grib_file ... ";
grib_index* idx=NULL;
char* keys;
char* default_keys = "mars";

grib_option grib_options[]={
/*  {id, args, help}, on, command_line, value */
    {"f",0,0,0,1,0},
	{"o:",0,0,1,1,"gribidx"},
    {"k:",0,0,0,1,0},
    {"V",0,0,0,1,0},
    {"M",0,0,0,1,0}
};

int grib_options_count=sizeof(grib_options)/sizeof(grib_option);

int main(int argc, char *argv[]) {
  int ret=grib_tool(argc,argv);
  return ret;
}

int grib_tool_before_getopt(grib_runtime_options* options) {
  return 0;
}

int grib_tool_init(grib_runtime_options* options) {

  int ret=0;
  grib_context* c=grib_context_get_default();
  if (grib_options_on("k:")) 
      keys=grib_options_get_option("k:");
  else
	  keys=default_keys;

  options->onlyfiles=1;

  idx=grib_index_new(c,keys,&ret);

  if (!idx || ret)
    grib_context_log(c,GRIB_LOG_FATAL,
        "Unable to create index %s",grib_get_error_message(ret));

  return 0;
}

int grib_tool_new_filename_action(grib_runtime_options* options,const char* file) {
	int ret=0;
	printf("--- %s: processing %s\n",grib_tool_name,file);
	ret=grib_index_add_file(idx,file);
	if (ret) {printf("error: %s\n",grib_get_error_message(ret)); exit(ret);}
   return 0;
}

int grib_tool_new_file_action(grib_runtime_options* options,grib_tools_file* file) {
   return 0;
}

int grib_tool_new_handle_action(grib_runtime_options* options, grib_handle* h) {
  return 0;
}

int grib_tool_skip_handle(grib_runtime_options* options, grib_handle* h) {
  return 0;
}

void grib_tool_print_key_values(grib_runtime_options* options,grib_handle* h) {
}

int grib_tool_finalise_action(grib_runtime_options* options) {
	grib_index_key* keys;
	grib_string_list* values;
	int first;
	grib_index_compress(idx);
	printf("--- %s: keys included in the index file %s:\n",
			grib_tool_name,options->outfile->name);
	printf("--- "); 
	first=1;
	keys=idx->keys;
	while (keys) {
		if (!first) printf(", ");
		printf("%s",keys->name);
		keys=keys->next;
		first=0;
	}
	printf("\n");
	keys=idx->keys;
	while (keys) {
		printf("--- %s = { ",keys->name);
		values=keys->values;
		first=1;
		while (values) {
			if (!first) printf(", ");
			printf("%s",values->value);
			first=0;
			values=values->next;
		}
		printf(" }\n");
		keys=keys->next;
	}
	printf("--- %d messages indexed\n",idx->count);

	if (idx->count)
		grib_index_write(idx,options->outfile->name);
	grib_index_delete(idx);
	return 0;
}
