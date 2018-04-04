/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: grib_debug
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */


#include "grib_tools.h"


grib_option grib_options[]={
/*  {id, args, help}, on, command_line, value */
    {"D",0,0,1,0,0},
    {"M",0,0,0,1,0},
    {"V",0,0,0,1,0},
    {"v",0,0,0,0,0},
    {"7",0,0,0,1,0},
    {"G",0,0,0,1,0},
    {"q",0,0,1,0,0}
};

char* grib_tool_description="Dump the content of a grib file in debug mode.";
char* grib_tool_name="grib_debug";
char* grib_tool_usage="[options] grib_file grib_file ...";

int grib_options_count=sizeof(grib_options)/sizeof(grib_option);

int main(int argc, char *argv[]) { return grib_tool(argc,argv);}

int grib_tool_before_getopt(grib_runtime_options* options) {
  return 0;
}

int grib_tool_init(grib_runtime_options* options) {
  options->dump_flags = GRIB_DUMP_FLAG_VALUES
                 |  GRIB_DUMP_FLAG_READ_ONLY;
  options->dump_mode = "debug";
  return 0;
}

int grib_tool_new_filename_action(grib_runtime_options* options,const char* file) {
   return 0;
}

int grib_tool_new_file_action(grib_runtime_options* options,grib_tools_file* file) {

   return 0;
}

int grib_tool_new_handle_action(grib_runtime_options* options, grib_handle* h) {

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
