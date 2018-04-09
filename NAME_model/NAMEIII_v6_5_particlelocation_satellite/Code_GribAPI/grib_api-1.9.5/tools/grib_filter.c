/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: grib_filter
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */

#include "grib_tools.h"


grib_option grib_options[]={
/*  {id, args, help}, on, command_line, value */
    {"f",0,0,0,1,0},
    {"f",0,0,1,0,0},
    {"F",0,0,1,0,0},
    {"o:",0,0,1,1,"filter.out"},
    {"q",0,0,1,0,0},
    {"M",0,0,0,1,0},
    {"I",0,0,1,0,0},
    {"V",0,0,0,1,0},
    {"g",0,0,0,1,0},
    {"G",0,0,0,1,0},
    {"7",0,0,0,1,0},
    {"v",0,0,0,1,0}
};
char* grib_tool_description="Apply the rules defined in rules_file to each grib "
   "message\n\tin the grib files provided as arguments.";
char* grib_tool_name="grib_filter";
char* grib_tool_usage="[options] rules_file "
                      "grib_file grib_file ...";
int fail=0;
int grib_options_count=sizeof(grib_options)/sizeof(grib_option);

int main(int argc, char *argv[]) { return grib_tool(argc,argv);}

int grib_tool_before_getopt(grib_runtime_options* options) {
  return 0;
}

int grib_tool_init(grib_runtime_options* options) {

  options->action = grib_action_from_filter(options->infile_extra->name);
  if (!options->action) {
      fprintf(stderr,"%s: error unable to create action\n",options->infile_extra->name);
      exit(1);
  }

  if ( options->outfile && options->outfile->name )
    options->action->context->outfilename=options->outfile->name;

  return 0;
}

int grib_tool_new_filename_action(grib_runtime_options* options,const char* file) {
   return 0;
}  

int grib_tool_new_file_action(grib_runtime_options* options,grib_tools_file* file) {
   return 0;
}

int grib_tool_new_handle_action(grib_runtime_options* options, grib_handle* h) {
  int err=0;
  
  if (options->current_infile->name) {
    size_t len=strlen(options->current_infile->name);
    grib_set_string(h,"file",options->current_infile->name,&len);
  }

  err=grib_handle_apply_action(h,options->action);
  if (err != GRIB_SUCCESS && options->fail && err!=GRIB_NOT_FOUND) {
     printf("ERROR: %s\n",grib_get_error_message(err));
     exit(err);
  }
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
