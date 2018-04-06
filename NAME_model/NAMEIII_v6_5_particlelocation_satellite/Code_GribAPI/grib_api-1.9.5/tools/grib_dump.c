/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: grib_dump
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */

#include "grib_tools.h"

grib_option grib_options[]={
/*  {id, args, help}, on, command_line, value*/
    {"S",0,0,1,0,0},
    {"O",0,"Octet mode. WMO documentation style dump.\n",0,1,0},
    {"D",0,0,0,1,0},
    {"P:",0,0,0,1,0},
    {"d",0,"Print all data values. Available only in C mode\n",0,1,0},
    {"C",0,0,0,1,0},
    {"t",0,0,0,1,0},
    {"H",0,0,0,1,0},
    {"a",0,0,0,1,0},
    {"w:",0,0,0,1,0},
    {"M",0,0,0,1,0},
    {"7",0,0,0,1,0},
    {"V",0,0,0,1,0},
    {"G",0,0,0,1,0},
    {"q",0,0,1,0,0}
};

char* grib_tool_description="Dump the content of a grib file in different formats.";
char* grib_tool_name="grib_dump";
char* grib_tool_usage="[options] grib_file grib_file ...";

int grib_options_count=sizeof(grib_options)/sizeof(grib_option);

/**
*grib_dump
*Dump the content of a grib file
*
*/
int main(int argc, char *argv[]) { return grib_tool(argc,argv);}

int grib_tool_before_getopt(grib_runtime_options* options) {
  return 0;
}

int grib_tool_init(grib_runtime_options* options) {

  int opt=grib_options_on("C")+grib_options_on("O")+grib_options_on("D");

  options->dump_mode = "default";


  if (opt > 1) {
    printf("%s: simultaneus C/O/D options not allowed\n",grib_tool_name);
    exit(1);
  }

  if (grib_options_on("C")) {
    options->dump_mode = "c_code";
    if (grib_options_on("d"))
      options->dump_flags = 0;
    else
      options->dump_flags = GRIB_DUMP_FLAG_NO_DATA;
  }

  if  (grib_options_on("O")) {
    options->dump_mode = "wmo";
    options->dump_flags = GRIB_DUMP_FLAG_CODED
               | GRIB_DUMP_FLAG_OCTECT
               | GRIB_DUMP_FLAG_VALUES
               | GRIB_DUMP_FLAG_READ_ONLY;
  }

  if (grib_options_on("D")) {
    options->dump_mode = "debug";
    options->dump_flags = GRIB_DUMP_FLAG_VALUES
               |  GRIB_DUMP_FLAG_READ_ONLY;
  }

  if (grib_options_on("a"))
      options->dump_flags |= GRIB_DUMP_FLAG_ALIASES;

  if (grib_options_on("t"))
      options->dump_flags |= GRIB_DUMP_FLAG_TYPE;

  if (grib_options_on("H"))
      options->dump_flags |= GRIB_DUMP_FLAG_HEXADECIMAL;

  if (grib_options_on("d"))
      options->dump_flags |= GRIB_DUMP_FLAG_ALL_DATA;

  return 0;
}

int grib_tool_new_filename_action(grib_runtime_options* options,const char* file) {
   return 0;
}

int grib_tool_new_file_action(grib_runtime_options* options,grib_tools_file* file) {
  char tmp[1024];
  if (!options->current_infile->name) return 0;
  sprintf(tmp,"FILE: %s ",options->current_infile->name);
  if (!grib_options_on("C"))
    fprintf(stdout,"***** %s\n",tmp);
   return 0;
}

int grib_tool_new_handle_action(grib_runtime_options* options, grib_handle* h) {
  long length=0;
  char tmp[1024];
  char identifier[100];
  size_t idlen=100;
  int i;
  GRIB_CHECK_NOLINE(grib_get_long(h,"totalLength",&length),0);

  for (i=0;i<options->print_keys_count;i++)
    grib_set_flag(h,options->print_keys[i].name,GRIB_ACCESSOR_FLAG_DUMP);

  sprintf(tmp,"MESSAGE %d ( length=%ld )",options->handle_count,length);
  if (!grib_options_on("C"))
    fprintf(stdout,"#==============   %-38s   ==============\n",tmp);
  if (!strcmp(options->dump_mode,"default")) {
    GRIB_CHECK_NOLINE(grib_get_string(h,"identifier",identifier,&idlen),0);
    printf("%s {\n",identifier);
  }

  grib_dump_content(h,stdout,options->dump_mode,options->dump_flags,0);

  if (!strcmp(options->dump_mode,"default"))
    printf("}\n");
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
