/**
 * Copyright 2005-2007 ECMWF
 *
 * Licensed under the GNU Lesser General Public License which
 * incorporates the terms and conditions of version 3 of the GNU
 * General Public License.
 * See LICENSE and gpl-3.0.txt for details.
 */

/*
 * C Implementation: grib_set
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */
#include "grib_tools.h"

grib_option grib_options[]={
  /*  {id, args, help}, on, command_line, value*/
  {"s:",0,0,1,1,0},
  {"r",0,0,0,1,0},
  {"d:",0,0,0,1,0},
  /*{"n:","noise percentage","\n\t\tAdd noise to the data values. The noise added is the given percentage of the data value.\n",0,1,0},*/
  {"p:",0,0,1,1,0},
  {"P:",0,0,0,1,0},
  {"w:","key[:{s/d/l}]=value,key[:{s/d/l}]=value,...",
  "\n\t\tWhere clause.\n\t\tSet is only executed for grib messages matching all the "
      "key/value constraints.\n\t\tIf a grib message does not match the constraints it is"
      " copied unchanged\n\t\tto the output_grib_file. This behaviour can be changed "
      "setting the option -S.\n\t\tFor each key a string (key:s) or a double (key:d) or"
      " a long (key:l)\n\t\ttype can be defined. Default type is string.\n",0,1,0},
    {"q",0,0,1,0,0},
    {"7",0,0,0,1,0},
    {"S",0,0,0,1,0},
    {"U",0,0,1,0,0},
    {"V",0,0,0,1,0},
    {"M",0,0,0,1,0},
    {"g",0,0,0,1,0},
    {"G",0,0,0,1,0},
    {"f",0,0,0,1,0},
    {"v",0,0,0,1,0}
};

char* grib_tool_description="Sets key/value pairs in the input grib file and writes"
"\n\teach message to the output_grib_file."
"\n\tIt fails when an error occurs (e.g. key not found).";
char* grib_tool_name="grib_set";
char* grib_tool_usage="[options] grib_file grib_file ... output_grib_file";

int grib_options_count=sizeof(grib_options)/sizeof(grib_option);

int main(int argc, char *argv[]) { return grib_tool(argc,argv);}

int grib_tool_before_getopt(grib_runtime_options* options) {
  return 0;
}

int grib_tool_init(grib_runtime_options* options) {
  if (options->set_values_count == 0 && !options->repack && !options->constant) {
    printf("ERROR: please provide some keys to set through the -s option or use the -r/-d options\n");
    exit(1);
  }
  if (options->verbose) options->print_header=1;
  /*if (grib_options_on("n:")) {
    noise=atof(grib_options_get_option("n:"));
    options->repack=1;
  }*/

  if (grib_options_on("n:") && grib_options_on("d:")) {
    printf("Error: -n and -d options are incompatible. Choose one of the two please.\n");
    exit(1);
  }
  
#if 0
  if (options->outfile && options->outfile->name) {
    options->outfile->file = fopen(options->outfile->name,"w");
    if(!options->outfile->file) {
      perror(options->outfile->name);
      exit(1);
    }
  }
#endif
  return 0;
}

int grib_tool_new_filename_action(grib_runtime_options* options,const char* file) {
   return 0;
   }


int grib_tool_new_file_action(grib_runtime_options* options,grib_tools_file* file) {
  return 0;
}

int grib_tool_new_handle_action(grib_runtime_options* options, grib_handle* h) {

  int i=0;
  int err=0;

  if (!options->skip) {
    double* v=NULL;
    size_t size=0;
    if ( options->repack ) {
      GRIB_CHECK_NOLINE(grib_get_size(h,"values",&size),0);

      v    = (double*)calloc(size,sizeof(double));
      if(!v) {
        fprintf(stderr,"failed to allocate %d bytes\n",(int)(size*sizeof(double)));
        exit(1);
      }

      GRIB_CHECK_NOLINE(grib_get_double_array(h,"values",v,&size),0);
    }

    if (options->set_values_count != 0)
      err=grib_set_values(h,options->set_values,options->set_values_count);

    if ( options->repack ) {

      if (grib_options_on("d:")) {
        for(i = 0; i< size; i++)
             v[i] =  options->constant;
      }
#if 0
      if (grib_options_on("n:")) {
        for(i = 0; i< size; i++)
          v[i] =  options->constant;
      }
#endif

      GRIB_CHECK_NOLINE(grib_set_double_array(h,"values",v,size),0);
      free(v);
    }

    if( err != GRIB_SUCCESS && options->fail) exit(err);
  }

  if (!options->skip || !options->strict) grib_tools_write_message(options,h);

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



