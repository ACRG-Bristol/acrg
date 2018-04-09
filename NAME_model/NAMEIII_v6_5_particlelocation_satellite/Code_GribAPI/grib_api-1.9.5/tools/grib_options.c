/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: grib_options
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */
#include <unistd.h>
#include <string.h>
#include "grib_tools.h"

char* names[]={"parameter","vertical","geography","data","mars","local"};
int names_count=6;

/* id,args,help */
grib_options_help grib_options_help_list[] ={
  {"a",0,"Dump aliases.\n"},
  {"b:","key,key,...",
     "\n\t\tAll the keys in this list are skipped in the comparison. Bit-by-bit compare on.\n"},
  {"B:","order by directive",
     "\n\t\tOrder by. The output will be ordered according the order by directive."
     "\n\t\tOrder by example: \"step asc, centre desc\" (step ascending and centre discending)\n"},
  {"c:","key[:l/d/s/n],key[:l/d/s/n],...",
   "\n\t\tOnly the listed keys or namespaces (:n) are compared. The optional letter after the colon is used "
   "\n\t\tto force the type in the comparison: l->integer, d->float, s->string, n->namespace."
   "\n\t\tSee -a option. Incompatible with -H option.\n"},
  {"d:","value",
   "\n\t\tSet all the data values to \"value\".\n"},
  {"e:","tolerance","\n\t\tOnly values whose difference is more than tolerance are considered different.\n"},
  {"f",0,"Force. Force the execution not to fail on error.\n"},
  {"F:","format","\n\t\tC style format for floating point values.\n"},
  {"g",0,"Copy GTS header. \n"},
  {"G",0,"GRIBEX compatibility mode.\n"},
  {"i:","index",
    "\n\t\tData value corresponding to the given index is printed.\n"},
  {"l:","Latitude,Longitude[,MODE,file]",
   "\n\t\tValue close to the point of a Latitude/Longitude."
   "\n\t\tAllowed values for MODE are:"
   "\n\t\t4 (4 values in the nearest points are printed) Default"
   "\n\t\t1 (the value at the nearest point is printed)"
   "\n\t\tfile (file is used as mask. The closer point with mask value>=0.5 is printed)\n"
       },
  {"n:","namespace",
   "\n\t\tAll the keys belonging to namespace are printed.\n"},
  {"m",0,"Mars keys are printed.\n"},
  {"o:","output_grib_file",
   "\n\t\tOutput grib is written to ouput_grib_file."
   "\n\t\tIf an ouput grib file is required and -o is not used, the"
   "ouput grib is written to filtered.out\n"},
  {"p:","key[:{s/d/l}],key[:{s/d/l}],...",
   "\n\t\tDeclaration of keys to print."
   "\n\t\tFor each key a string (key:s) or a double (key:d) or a long (key:l)"
   "\n\t\ttype can be requested. Default type is string.\n"},
  {"q",0,"Quiet.\n"},
  {"r",0,"Repack data. Sometimes after setting some keys involving properties"
         "\n\t\tof the packing algorithm a repacking of data is needed."
         "\n\t\tThis repacking is performed setting this -r option.\n"},
  {"s:","key[:{s/d/l}]=value,key[:{s/d/l}]=value,...",
   "\n\t\tKey/values to set."
   "\n\t\tFor each key a string (key:s) or a double (key:d) or a long (key:l)"
   "\n\t\ttype can be defined. By default the native type is set.\n"},
  {"t",0,"Print type information.\n"},
  {"w:","key[:{s/d/l}]{=/!=}value,key[:{s/d/l}]{=/!=}value,...",
   "\n\t\tWhere clause."
   "\n\t\tGrib messages are processed only if they match all the"
   " key/value constraints."
   "\n\t\tA valid constraint is of type key=value or key!=value."
   "\n\t\tFor each key a string (key:s) or a double (key:d) or"
   " a long (key:l)\n\t\ttype can be specified. Default type is string.\n"},
  {"v",0,"Verbose.\n"},
  {"7",0,"Does not fail when the message has wrong length\n"},
  {"A:","absolute error\n",
  "\tCompare floating point values using the absolute error as tolerance.\n\t\tDefault is absolute error=0\n"},
  {"C",0,"C code mode. A C code program generating the grib message is dumped.\n"},
  {"D",0,"Debug mode.\n"},
  {"H",0,"Print octet content in hexadecimal format.\n"},
  {"M",0,"Multi-grib support off. Turn off support for multiple fields in single grib message\n"},
  {"O",0,"Octet mode. WMO documentation style dump.\n"},
  {"P:","key[:{s/d/l}],key[:{s/d/l}],...",
   "\n\t\tAs -p adding the declared keys to the default list.\n"},
  {"R:","key1=relative_error1,key2=relative_error2,...\n",
   "\tCompare floating point values using the relative error as tolerance."
"\n\t\tkey1=relative_error will compare key1 using relative_error1."
"\n\t\tall=relative_error will compare all the floating point keys using relative_error. Default all=0.\n"},
  {"S",0,"Strict. Only grib messages matching all the constraints are copied to"
   "\n\t\tthe output file\n"},
  {"V",0,"Version.\n"},
  {"W:","width","\n\t\tMinimum width of each column in output. Default is 10.\n"},
  {"x",0,"Fast parsing option, only headers are loaded.\n"}

};

int grib_options_help_count=sizeof(grib_options_help_list)/sizeof(grib_options_help);

static void set_value(grib_values* value,char* str,int equal);

char* grib_options_get_option(const char* id) {
  int i=0;

  for (i=0; i<grib_options_count; i++) {
    if (!strcmp(id,grib_options[i].id))
      return grib_options[i].value;
  }
  return NULL;
}

int grib_options_command_line(const char* id) {
  int i=0;

  for (i=0; i<grib_options_count; i++) {
    if (!strcmp(id,grib_options[i].id))
      return grib_options[i].command_line;
  }
  return 0;
}

int grib_options_on(const char* id) {
  int i=0;

  for (i=0; i<grib_options_count; i++) {
    if (!strcmp(id,grib_options[i].id))
      return grib_options[i].on;
  }
  return 0;
}

int grib_get_runtime_options(int argc,char** argv,grib_runtime_options* options) {
  int i=0,c=0;
  char* optstr=(char*)calloc(1,2*grib_options_count*sizeof(char));

  if (!optstr) return GRIB_OUT_OF_MEMORY;

  for (i=0;i<grib_options_count;i++)
    if (grib_options[i].command_line)
     strncat(optstr,grib_options[i].id,2);


  while ((c = getopt (argc, argv, optstr)) != -1) {
    if (c == '?') usage(argv[0]);
    i=0;
    while ( i < grib_options_count && grib_options[i].id[0] != c) i++;
    grib_options[i].on=1;
    if ( grib_options[i].id[1] == ':' )
      grib_options[i].value=optarg;
  }
  free(optstr);
  return 0;
}

int grib_process_runtime_options(grib_context* context,int argc,char** argv,grib_runtime_options* options) {
  int i=0,ret=0;
  int has_output=0;int has_input_extra=0,nfiles=0;
  char *karg=NULL,*warg=NULL,*sarg=NULL,*barg=NULL;
  
  if (grib_options_on("V")) {
    printf("\ngrib_api Version ");
    grib_print_api_version(stdout);
    printf("\n\n");
    exit(0);
  }

  if (grib_options_on("B:"))
    options->orderby=grib_options_get_option("B:");

  if (grib_options_on("x")) options->headers_only=1;

  if (grib_options_on("F:"))
    options->format=grib_options_get_option("F:");
  else 
    options->format=strdup("%g");

  if (grib_options_on("i:")) {
    options->index_on=1;
    options->index=atoi(grib_options_get_option("i:"));
  }

  if (grib_options_on("l:"))
    options->latlon=grib_options_get_option("l:");

  has_output=grib_options_on("U");
  has_input_extra=grib_options_on("I");
  options->repack=grib_options_on("r");
  options->gts=grib_options_on("g");

  if (grib_options_on("d:")) {
    options->constant=atof(grib_options_get_option("d:"));
    options->repack=1;
  }

  if (grib_options_on("G")) grib_gribex_mode_on(context);

  nfiles=argc-optind;
  if ( nfiles < (1+has_output+has_input_extra) && !options->infile ) usage(argv[0]);

  if (has_input_extra) {
    options->infile_extra=(grib_tools_file*)calloc(1,sizeof(grib_tools_file));
    options->infile_extra->name=argv[optind];
  }

  if (!options->infile) {
    for (i=optind+has_input_extra;i<argc-has_output;i++) {
      grib_tools_file* p=NULL;
      grib_tools_file* infile=(grib_tools_file*)calloc(1,sizeof(grib_tools_file));
      infile->name=argv[i];
      if (!options->infile) options->infile=infile;
      else {
        p=options->infile;
        while(p->next) p=p->next;
        p->next=infile;
      }
    }
  }

  if (has_output) {
    options->outfile=(grib_tools_file*)calloc(1,sizeof(grib_tools_file));
    options->outfile->name=strdup(argv[argc-1]);
  }

  if (grib_options_on("o:")) {
    options->outfile=(grib_tools_file*)calloc(1,sizeof(grib_tools_file));
    options->outfile->name=grib_options_get_option("o:");
  }

  options->print_number=grib_options_on("N");
  options->print_header=grib_options_on("H");
  options->verbose=grib_options_on("v");

  if (grib_options_on("q") && grib_options_command_line("q"))
    options->verbose=0;

  options->fail=!grib_options_on("f");
  if (grib_options_get_option("W:"))
    options->default_print_width=atoi(grib_options_get_option("W:"));

  if (grib_options_on("n:"))
    options->name_space=grib_options_get_option("n:");
  
  if (grib_options_on("m"))
    options->name_space=strdup("mars");

  
  if (grib_options_on("P:")) karg=grib_options_get_option("P:");
  else if (grib_options_on("p:")) {
    karg=grib_options_get_option("p:");
    options->name_space=NULL;
  }
  
  options->requested_print_keys_count=MAX_KEYS;
  GRIB_CHECK_NOLINE(grib_options_get_keys(argv[0],karg,0,GRIB_TYPE_UNDEFINED,
      options->requested_print_keys,&(options->requested_print_keys_count)),0);


  options->strict=grib_options_on("S");

  if (grib_options_on("M")) grib_multi_support_off(context);
  else grib_multi_support_on(context);

  if (grib_options_on("g")) grib_gts_header_on(context);
  else grib_gts_header_off(context);

  if (grib_options_on("V")) {
    printf("\ngrib_api Version ");
    grib_print_api_version(stdout);
    printf("\n\n");
  }

  if (grib_options_on("s:")) {
    sarg=grib_options_get_option("s:");
    options->set_values_count=MAX_KEYS;
    ret=grib_options_get_keys(argv[0],sarg,1,GRIB_TYPE_UNDEFINED,options->set_values,&(options->set_values_count));
  }

  if (grib_options_on("R:")) {
    sarg=grib_options_get_option("R:");
    options->tolerance_count=MAX_KEYS;
    ret=grib_options_get_keys(argv[0],sarg,1,GRIB_TYPE_DOUBLE,options->tolerance,&(options->tolerance_count));
  }

  if (grib_options_on("b:")) {
    barg=grib_options_get_option("b:");
    options->set_values_count=MAX_KEYS;
    ret=grib_options_get_keys(argv[0],barg,0,GRIB_TYPE_STRING,options->set_values,&(options->set_values_count));
  }

  if (grib_options_on("c:")) {
    sarg=grib_options_get_option("c:");
    options->compare_count=MAX_KEYS;
    ret=grib_options_get_keys(argv[0],sarg,0,GRIB_TYPE_UNDEFINED,options->compare,
                &(options->compare_count));
  }
  if (grib_options_on("e")) {
    for (i=0;i<names_count; i++) {
      options->compare[i+options->compare_count].name=names[i];
      options->compare[i+options->compare_count].type=GRIB_NAMESPACE;
    }
    options->compare_count+=names_count;
  }

  warg=grib_options_get_option("w:");

  options->constraints_count=MAX_KEYS;
  ret=grib_options_get_keys(argv[0],warg,1,GRIB_TYPE_STRING,
      options->constraints,&(options->constraints_count));

  return GRIB_SUCCESS;
}

int grib_options_get_keys(char* progname,char* arg,int values_required,int default_type,grib_values values[],int* count) {
   char* p;
   int i=0;
   if (arg == NULL) {
     *count=0;
     return GRIB_SUCCESS;
   }
   p=strtok(arg,",");
   while (p != NULL) {
     values[i].name=(char*)calloc(1,strlen(p)+1);
     strcpy((char*)values[i].name,p);
     p=strtok(NULL,",");
     i++;
     if (i > *count) return GRIB_ARRAY_TOO_SMALL;
   }
   *count=i;

   for (i=0;i<*count;i++) {
     int equal=1;
     char* value=NULL;
     if (values_required) {
       p=(char*)values[i].name;
       while( *p != '=' && *p != '!' && *p != '\0') p++;
       if (*p == '=') {
         *p='\0';
         p++;
         value=p;
         equal=1;
       } else if (*p == '!' && *(++p) == '=' ) {
         *p='\0';
         *(p-1)='\0';
         p++;
         value=p;
         equal=0;
       } else usage( progname );
     }
     p=(char*)values[i].name;
     while(*p != ':' && *p != '\0') p++;
     if ( *p == ':' ) {
       values[i].type = grib_type_to_int(*(p+1));
       if (*(p+1)=='n') values[i].type = GRIB_NAMESPACE;
       *p='\0';
       p++;
     } else {
       values[i].type = default_type;
     }
     if (values_required) {
      if (strlen(value)==0) {
        printf("%s error: no value provided for key \"%s\"\n",grib_tool_name,values[i].name);
        exit(GRIB_INVALID_ARGUMENT);
      }
      set_value(&values[i],value,equal);
     }
   }
   return GRIB_SUCCESS;
}

static void set_value(grib_values* value,char* str,int equal) {
  char *p=0,*q=0,*s=0;
  char buf[1000]={0,};
  grib_context* c=grib_context_get_default();
  
  value->equal=equal;
  q=str;
  
  while (*q != '/' && *q!=0 ) q++;
  if (*q=='/') {
    s=grib_context_strdup(c,q+1);
    value->next=grib_context_malloc_clear(c,sizeof(grib_values));
    value->next->type=value->type;
    value->next->name=grib_context_strdup(c,value->name);
    set_value(value->next,s,equal);
    grib_context_free(c,s);
  }

  memcpy(buf,str,q-str);
  
  switch (value->type) {
    case GRIB_TYPE_DOUBLE:
      value->double_value = strtod(buf,&p);
      if(*p != 0) value->has_value=1;
      else if (!strcmp(str,"missing") ||
               !strcmp(str,"MISSING") ||
               !strcmp(str,"Missing")  ) {
        value->type=GRIB_TYPE_MISSING;
        value->has_value=1;
      }
    break;
    case GRIB_TYPE_LONG:
      value->long_value = strtol(buf,&p,10);
      if(*p != 0) value->has_value=1;
      else if (!strcmp(buf,"missing") ||
                !strcmp(buf,"MISSING") ||
                !strcmp(buf,"Missing")  ) {
        value->type=GRIB_TYPE_MISSING;
        value->has_value=1;
      }
    break;
    case GRIB_TYPE_STRING:
      if ( !strcmp(buf,"missing") ||
            !strcmp(buf,"MISSING") ||
            !strcmp(buf,"Missing")  ) {
            value->type=GRIB_TYPE_MISSING;
            value->has_value=1;
      } else {
        value->string_value = grib_context_strdup(c,buf);
        value->has_value=1;
      }
      break;
    case GRIB_TYPE_UNDEFINED:
      value->long_value = strtol(buf,&p,10);
      if(*p == 0) {
        value->type = GRIB_TYPE_LONG;
        value->has_value=1;
      } else {
        value->double_value = strtod(buf,&p);
        if(*p == 0) {
            value->type = GRIB_TYPE_DOUBLE;
            value->has_value=1;
        } else if ( !strcmp(buf,"missing") ||
                    !strcmp(buf,"MISSING") ||
                    !strcmp(buf,"Missing")  ) {
            value->type=GRIB_TYPE_MISSING;
            value->has_value=1;
        } else {
          value->string_value = grib_context_strdup(c,buf);
          value->type         = GRIB_TYPE_STRING;
          value->has_value=1;
        }
      }
      break;
   }
}

char* grib_options_get_help(char* id) {
   int i=0;
   char msg[]="ERROR: help not found for option ";
   char* err=(char*)calloc(1,sizeof(msg)+3);
   sprintf(err,"%s%c\n",msg,*id);
   for (i=0; i<grib_options_count;i++) {
     if (!strcmp(id,grib_options[i].id)) {
       if (grib_options[i].help != NULL) return grib_options[i].help;
       else break;
     }
   }
   for (i=0; i<grib_options_help_count;i++) {
     if (!strcmp(id,grib_options_help_list[i].id)) {
        return grib_options_help_list[i].help != NULL ? grib_options_help_list[i].help : err;
     }
   }
   return err;
}

char* grib_options_get_args(char* id) {
   int i=0;
   char empty[]="";
   char msg[]="ERROR: help not found for option -";
   char* err=NULL;
   if (id[1] != ':') return strdup(empty);
   err=(char*)calloc(1,sizeof(msg)+3);
   sprintf(err,"%s%c\n",msg,*id);
   for (i=0; i<grib_options_count;i++) {
     if (!strcmp(id,grib_options[i].id)) {
       if (grib_options[i].args != NULL) return grib_options[i].args;
       else break;
     }
   }
   for (i=0; i<grib_options_help_count;i++) {
     if (!strcmp(id,grib_options_help_list[i].id)) {
        return grib_options_help_list[i].args != NULL ? grib_options_help_list[i].args : err;
     }
   }
   return err;
}

void usage() {
  int i=0;
  printf("\nNAME \t%s\n\n",grib_tool_name);
  printf("DESCRIPTION\n\t%s\n\n",grib_tool_description);
  printf("USAGE \n\t%s %s\n\n",grib_tool_name,grib_tool_usage);
  printf("OPTIONS\n");
  for (i=0;i<grib_options_count;i++) {
    if (grib_options[i].command_line)
      printf("\t-%c %s\t%s",grib_options[i].id[0],
              grib_options_get_args(grib_options[i].id),
              grib_options_get_help(grib_options[i].id));
  }
  printf("\n\n");
  exit(1);
}

void usage_doxygen() {
  int i=0;
  printf("/*!  \\page %s %s\n",grib_tool_name,grib_tool_name);
  printf("\\section DESCRIPTION \n %s\n\n",grib_tool_description);
  printf("\\section USAGE \n %s \n %s\n\n",grib_tool_name,grib_tool_usage);
  printf("\\section OPTIONS \n");
  for (i=0;i<grib_options_count;i++) {
    if (grib_options[i].command_line) {
      printf("-%c %s \\n",
              grib_options[i].id[0],
              grib_options_get_args(grib_options[i].id));
      printf(" %s \\n \\n ",
              grib_options_get_help(grib_options[i].id));
    }
  }
  exit(1);
}

#if 0
void usage_doxygen() {
  int i=0;
  printf("/*!  \\page %s %s\n",grib_tool_name,grib_tool_name);
  printf("\\section DESCRIPTION \n%s\n\n",grib_tool_description);
  printf("\\section USAGE \n%s \n%s\n\n",grib_tool_name,grib_tool_usage);
  printf("\\section OPTIONS\n");
  printf("<table frame=void border=0>\n");
  for (i=0;i<grib_options_count;i++) {
    if (grib_options[i].command_line) {
    printf("<tr>\n");
      printf("<td colspan=2>-%c %s</td>\n",
              grib_options[i].id[0],
              grib_options_get_args(grib_options[i].id));
    printf("</tr><tr>\n");
      printf("<td width=20></td><td>%s</td>",
              grib_options_get_help(grib_options[i].id));
    printf("</tr><tr><td></td></tr>\n");
  }
  }
  printf("</table>\n");
  exit(1);
}
#endif

