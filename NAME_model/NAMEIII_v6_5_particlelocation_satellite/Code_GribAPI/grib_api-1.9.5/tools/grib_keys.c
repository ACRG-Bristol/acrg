/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: grib_get
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */

#include "grib_tools.h"
#include "grib_msg_types.h"

grib_option grib_options[]={
/*  {id, args, help}, on, command_line, value*/
    {"T:","type","\n\t\tTo print the keys available in the given grib type."
        "\n\t\tFor a list of the availeble types see -L option.\n"
    ,0,1,0},
    {"F:","file","\n\t\tTo print the keys available in the grib file.\n",0,1,0},
    {"x",0,"Print the extended set of keys.\n",0,1,0},
    {"c",0,"Print only coded keys.\n",0,1,0},
    {"L",0,"List of available types.\n",0,1,0},
    {"t",0,0,0,1,0},
    {"a",0,0,0,1,0}
};

char* grib_tool_description="Lists the keys available for a type of grib (-T option) or \n"
                 "\tin a grib message from a file (-F option).";
char* grib_tool_name="grib_keys";
char* grib_tool_usage="[options]";
char name[1024];

int grib_options_count=sizeof(grib_options)/sizeof(grib_option);

int main(int argc, char *argv[]) { return grib_tool(argc,argv);}

int grib_tool_before_getopt(grib_runtime_options* options) {
  grib_context *c=grib_context_get_default( );
  options->dump_mode=strdup("keys");
  if (grib_options_on("T:")) {
    grib_tools_file* infile=(grib_tools_file*)calloc(1,sizeof(grib_tools_file));
    sprintf(name,"%s",grib_options_get_option("T:"));
    infile->name=grib_external_template_path(c,name);
    infile->next=NULL;
    options->infile=infile;
  } else if (grib_options_on("F:")) {
    grib_tools_file* infile=(grib_tools_file*)calloc(1,sizeof(grib_tools_file));
    sprintf(name,"%s",grib_options_get_option("F:"));
    infile->name=name;
    infile->next=NULL;
    options->infile=infile;
  }
  if (grib_options_on("L")) {
    int i;
    for (i=0;i<types_count;i++) printf("%s\n",types_list[i]);
    exit(0);
  }

  return 0;
}

int grib_tool_init(grib_runtime_options* options) {
  options->print_header=0;
  options->print_statistics=0;
  options->default_print_width=-1;

  options->dump_flags =GRIB_DUMP_FLAG_DUMP_OK;

  if (grib_options_on("x")) {
    options->dump_flags =0;
  } else if (grib_options_on("c")) {
    options->dump_flags =GRIB_DUMP_FLAG_CODED;
  }

  if (grib_options_on("a"))
      options->dump_flags |= GRIB_DUMP_FLAG_ALIASES;

  if (grib_options_on("t"))
      options->dump_flags |= GRIB_DUMP_FLAG_TYPE;

  if (!grib_options_on("T:") && !grib_options_on("F:") && !grib_options_on("A") ) {

  }

  return 0;
}
int grib_tool_new_filename_action(grib_runtime_options* options,const char* file) {
   return 0;
      }


int grib_tool_new_file_action(grib_runtime_options* options,grib_tools_file* file) {
  if (!options->current_infile->name) return 0;

  fprintf(stdout,"=================== %s \n",name);

  return 0;
}

int grib_tool_new_handle_action(grib_runtime_options* options,grib_handle* h) {
  if(grib_options_on("F:")) printf("------ %d --------\n",options->handle_count);
  grib_dump_content(h,stdout,options->dump_mode,options->dump_flags,0);
  return 0;
}

int grib_tool_skip_handle(grib_runtime_options* options, grib_handle* h) {
  grib_handle_delete(h);
  return 0;
}

void grib_tool_print_key_values(grib_runtime_options* options,grib_handle* h) {
  grib_print_key_values(options,h);
}

int grib_tool_finalise_action(grib_runtime_options* options) {
  return 0;
}
