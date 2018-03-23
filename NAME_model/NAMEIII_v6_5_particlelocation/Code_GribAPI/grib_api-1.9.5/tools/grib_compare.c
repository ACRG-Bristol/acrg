/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/*
 * C Implementation: grib_compare
 *
 * Author: Enrico Fucile <enrico.fucile@ecmwf.int>
 *
 *
 */

#include "grib_tools.h"


GRIB_INLINE static int grib_inline_strcmp(const char* a,const char* b) {
  if (*a != *b) return 1;
  while((*a!=0 && *b!=0) &&  *(a) == *(b) ) {a++;b++;}
  return (*a==0 && *b==0) ? 0 : 1;
}

GRIB_INLINE static int grib_inline_rstrcmp(const char* a,const char* b) {
  char* p=(char*)a;
  char* q=(char*)b;
  while (*p != 0) p++;
  while (*q != 0) q++;
  q--;p--;
  if (*q != *p) return 1;
  while((p!=a && q!=b) &&  *(p) == *(q) ) {p--;q--;}
  return (q==b) ? 0 : 1;
}

typedef double (*compare_double_proc) (double*,double*,double*);

typedef struct grib_error grib_error;
struct grib_error {
  char* key;
  int count;
  grib_error* next;
};

grib_error* error_summary;

compare_double_proc compare_double;
double global_tolerance=0;
int packingCompare=0;
grib_string_list* blacklist=0;
int compareAbsolute=1;

static int compare_handles(grib_handle* h1,grib_handle* h2,grib_runtime_options* options);
static int compare_values(grib_runtime_options* options,grib_handle* h1,grib_handle *h2,const char *name,int type);
int error=0;
int count=0;
int lastPrint=0;
int force=0;
double maxAbsoluteError = 1e-19;
int onlyListed=1;
int headerMode=0;
int morein1=0;
int morein2=0;
int listFromCommandLine;
int verbose=0;
int tolerance_factor=1;

GRIB_INLINE static double compare_double_absolute(double *a,double *b,double *err) {
  return fabs(*a-*b) > *err ? fabs(*a-*b) : 0;
}

static double compare_double_relative(double *a,double *b,double *err) {
  double relativeError;

  if(fabs(*a) <= maxAbsoluteError || fabs(*b) <= maxAbsoluteError)
    relativeError = fabs(*a-*b);
  else if (fabs(*b) > fabs(*a))
    relativeError = fabs((*a-*b) / *b);
  else
    relativeError = fabs((*a-*b) / *a);
  
  return relativeError > *err ? relativeError : 0;
}

static int blacklisted(const char* name) {
  grib_string_list* b=blacklist;
  while (b) {
    if (!strcmp(name,b->value))
      return 1;
    b=b->next;
  }
  return 0;
}

static double relative_error(double a,double b,double err) {
  double relativeError;
  double maxAbsoluteError = 1e-19;

  if(fabs(a) <= maxAbsoluteError || fabs(b) <= maxAbsoluteError)
    relativeError = fabs(a-b);
  else if (fabs(b) > fabs(a))
    relativeError = fabs((a-b) / b);
  else
    relativeError = fabs((a-b) / a);
  
  return relativeError ;
}


grib_option grib_options[]={
/*  {id, args, help}, on, command_line, value*/
    {"r",0,"Compare files in which the messages are not in the same order. This option is time expensive.\n",0,1,0},
    {"b:",0,0,0,1,0},
    {"e",0,"Edition independent compare. It is used to compare grib edition 1 and 2.\n",0,1,0},
    {"c:",0,0,0,1,0},
    {"S:","start","First field to be processed.\n",0,1,0},
    {"E:","end","Last field to be processed.\n",0,1,0},
    {"a",0,"-c option modifier. The keys listed with the option -c will be added to the list of keys compared without -c.\n"
        ,0,1,0},
    {"H",0,"Compare only message headers. Bit-by-bit compare on. Incompatible with -c option.\n",0,1,0},
    {"R:",0,0,0,1,0},
    {"A:",0,0,0,1,0},
    {"P",0,"Compare data values using the packing error as tolerance.\n",0,1,0},
    {"T:","factor","Compare data values using factor multiplied by the tolerance specified in options -P -R -A.\n",0,1,0},
    {"w:",0,0,0,1,0},
    {"f",0,0,0,1,0},
    {"F",0,0,1,0,0},
    {"q",0,0,1,0,0},
    {"M",0,0,1,0,0},
    {"I",0,0,1,0,0},
    {"V",0,0,0,1,0},
    {"7",0,0,0,1,0},
    {"v",0,0,0,1,0}
};

grib_handle* h1=NULL;
int counter=0;
int start=-1;
int end=-1;

char* grib_tool_description=
         "Compare grib messages contained in two files."
    "\n\tIf some differences are found it fails returning an error code."
    "\n\tFloating point values are compared exactly by default, different tolerance can be defined see -P -A -R."
    "\n\tDefault behaviour: absolute error=0, bit-by-bit compare, same order in files.";

char* grib_tool_name="grib_compare";
char* grib_tool_usage="[options] "
                      "grib_file grib_file";

int grib_options_count=sizeof(grib_options)/sizeof(grib_option);

int main(int argc, char *argv[]) { return grib_tool(argc,argv);}

int grib_tool_before_getopt(grib_runtime_options* options) {
  return 0;
}

int grib_tool_init(grib_runtime_options* options) {
  int ret=0,i;
  int nfiles=1;
  char orderby[]="md5Headers";
  grib_context* context=grib_context_get_default();

  options->strict=1;
  if (grib_options_on("S:")) 
      start=atoi(grib_options_get_option("S:"));

  if (grib_options_on("E:")) 
      end=atoi(grib_options_get_option("E:"));

  if (grib_options_on("f")) force=1;
  else force=0;

  verbose = grib_options_on("v");

  listFromCommandLine=0;
  if (grib_options_on("c:") || grib_options_on("e")) 
    listFromCommandLine=1;

  if (grib_options_on("a")) onlyListed=0;
  else onlyListed=1;

  if (grib_options_on("H")) headerMode=1;
  else headerMode=0;

  if (grib_options_on("H") && grib_options_on("c:")) {
    printf("Error: -H and -c options are incompatible. Choose one of the two please.\n");
    exit(1);
  }
  if (grib_options_on("a") && !grib_options_on("c:")) {
    printf("Error: -a option requires -c option. Please define a list of keys with the -c option.\n");
    exit(1);
  }
    
  if (grib_options_on("b:")) {
    grib_string_list *next=0;
    int i=0;
    blacklist=grib_context_malloc_clear(context,sizeof(grib_string_list));
    blacklist->value=grib_context_strdup(context,options->set_values[0].name);
    next=blacklist;
    for (i=1;i<options->set_values_count;i++) {
      next->next=grib_context_malloc_clear(context,sizeof(grib_string_list));
      next->next->value=grib_context_strdup(context,options->set_values[i].name);
      next=next->next;
    }
    context->blacklist=blacklist;
  }

  if (grib_options_on("r")) {
    char* filename[1];
    filename[0]=options->infile_extra->name;
    options->random=1;
    options->orderby=strdup(orderby);
    options->idx=grib_fieldset_new_from_files(context,filename,
        nfiles,0,0,0,orderby,&ret);
    if (ret) {
      printf("unable to create index for input file %s (%s)",
                       options->infile_extra->name,grib_get_error_message(ret));
      exit(ret);
    }
  } else {
    options->random=0;
    options->infile_extra->file=fopen(options->infile_extra->name,"r");

    if (!options->infile_extra->file) {
      perror(options->infile_extra->name);
      exit(1);
    }
  }

  global_tolerance=0;
  compare_double= &compare_double_absolute;
  if (grib_options_on("R:")) {
    global_tolerance=0;
    for (i=0;i<options->tolerance_count;i++) {
      if (!strcmp((options->tolerance[i]).name,"all")) {
        global_tolerance=(options->tolerance[i]).double_value;
        break;
      }
      if (!strcmp((options->tolerance[i]).name,"global")) {
        global_tolerance=(options->tolerance[i]).double_value;
        break;
      }
    }
    compare_double= &compare_double_relative;
    compareAbsolute=0;
  }
  if (grib_options_on("A:")){
    if (grib_options_on("R:")) {
      maxAbsoluteError = atof(grib_options_get_option("A:"));
    } else {
      compare_double= &compare_double_absolute;
      global_tolerance = atof(grib_options_get_option("A:"));
    }
  }
  if (grib_options_on("P")) {
    packingCompare=1;
    compare_double= &compare_double_absolute;
  }

  if (grib_options_on("T:")) 
    tolerance_factor=atof(grib_options_get_option("T:"));
  
  
  return 0;
}

int grib_tool_new_filename_action(grib_runtime_options* options,const char* file) {
   return 0;
}
int grib_tool_new_file_action(grib_runtime_options* options,grib_tools_file* file) {
   return 0;
}

static void printInfo(grib_handle* h) {
  char shortName[254]={0,};
  char levelType[254]={0,};
  char level[254]={0,};
  char paramId[254]={0,};
  char packingType[254]={0,};
  char gridType[254]={0,};
  char identifier[254]={0,};
  size_t len=254;
  char stepRange[254]={0,};
  if (lastPrint==count) return;

  len=254;
  grib_get_string(h,"shortName",shortName,&len);
  len=254;
  grib_get_string(h,"stepRange",stepRange,&len);
  len=254;
  grib_get_string(h,"levelType",levelType,&len);
  len=254;
  grib_get_string(h,"level",level,&len);
  len=254;
  grib_get_string(h,"paramId",paramId,&len);
  len=254;
  grib_get_string(h,"packingType",packingType,&len);
  len=254;
  grib_get_string(h,"gridType",gridType,&len);
  len=254;
  grib_get_string(h,"identifier",identifier,&len);

  printf("\n-- %s #%d -- shortName=%s paramId=%s stepRange=%s levelType=%s level=%s packingType=%s gridType=%s --\n",
         identifier,count,shortName,paramId,stepRange,levelType,level,packingType,gridType);
  lastPrint=count;

}

static void print_index_key_values(grib_index* index,int counter) {
	grib_index_key* keys=index->keys;
	printf("== %d == ",counter);
	while (keys) {
		printf("%s=%s ",keys->name,keys->value);
		keys=keys->next;
	}
	printf("\n");
}
int grib_tool_new_handle_action(grib_runtime_options* options, grib_handle* h) {
  int err=0;
  count++;

  if (options->through_index) {
	grib_index* idx1=options->index1;
	verbose=0;
	counter++;

	if ( start>0 && counter < start ) return 0;
	if ( end>0 && counter > end ) { 
		options->stop=1;
		return 0;
	}

	grib_index_search_same(idx1,h);
	h1=grib_handle_new_from_index(idx1,&err);
	if (options->verbose) {
		off_t offset=0;
		char* filename=grib_get_field_file(options->index2,&offset);
		printf("file1=\"%s\" ",filename);
		filename=grib_get_field_file(options->index1,&offset);
		printf("file2=\"%s\" \n",filename);
		print_index_key_values(options->index1,counter);
	}

	if (!h1) {
		if (!options->verbose)
			print_index_key_values(idx1,counter);
		printf("====== NOT FOUND in %s\n",options->infile->name);
	}

	if (!h1 || err!= GRIB_SUCCESS ) {
		morein1++;
		if (h1) grib_handle_delete(h1);
		return 0;
    }

    if(compare_handles(h,h1,options)) {
		error++;
		if (!force) exit(1);
    }

	grib_handle_delete(h1);

	return 0;
	
  } else if (options->random) 
    h1 = grib_fieldset_next_handle(options->idx,&err);
  else
    h1=grib_handle_new_from_file(h->context,options->infile_extra->file,&err);
  
  if (!h1 || err!= GRIB_SUCCESS ) {
    morein2++;
    if (h1) grib_handle_delete(h1);
    return 0;
  }

  if(compare_handles(h1,h,options)) {
    error++;
    if (!force) exit(1);
  }

  grib_handle_delete(h1);

  return 0;
}

int grib_tool_skip_handle(grib_runtime_options* options, grib_handle* h) {
  int err=0;
  if (!options->through_index && !options->random)  {
	  h1=grib_handle_new_from_file(h->context,options->infile_extra->file,&err);
  
	  if (!h1 || err!= GRIB_SUCCESS) 
		morein2++;

	  grib_handle_delete(h1);
	  
	  
  }
  
  grib_handle_delete(h);
  count++;
  
  return 0;
}

void grib_tool_print_key_values(grib_runtime_options* options,grib_handle* h) {
  grib_print_key_values(options,h);
}

int grib_tool_finalise_action(grib_runtime_options* options) {
  grib_error* e=error_summary;
  int err=0;
  grib_context* c=grib_context_get_default();
  error+=morein1+morein2;

  /*if (grib_options_on("w:")) return 0;*/
  
  if (error) {
    printf("\n## ERRORS SUMMARY #######\n");
  }
  while ((h1=grib_handle_new_from_file(c,options->infile_extra->file,&err))) {
      morein1++;
      if (h1) grib_handle_delete(h1);
  }
  if (morein1>0) {
    printf("##\n## Different number of messages \n");
    printf("## %d more messages in %s than in %s\n",morein1,
           options->infile_extra->name,options->infile->name);
  }

  if (morein2>0) {
    printf("##\n## Different number of messages \n");
    printf("## %d more messages in %s than in %s\n",morein2,
           options->infile->name,options->infile_extra->name);
  }
  
  if (error) {
    printf("##\n## Summary of different key values \n");
    while (e) {
      printf ("## %s ( %d different )\n",e->key,e->count);
      e=e->next;
    }
    
    printf("##\n## %d different messages out of %d\n\n",error,count);
  }
  if (options->through_index) {
  	grib_index_delete(options->index1);
  	grib_index_delete(options->index2);
  }

  
  if (error !=0) exit(1);
  return 0;
}

static void save_error(grib_context* c,const char* key) {
  grib_error* e=0;
  grib_error* next=0;
  int saved=0;

  if (!error_summary) {
     error_summary=grib_context_malloc_clear(c,sizeof(grib_error));
     error_summary->count=1;
     error_summary->key=grib_context_strdup(c,key);
     return;
  }

  e=error_summary;
  next=e;

  while (next) {
    if (!strcmp(next->key,key)) {
      next->count++;
      saved=1;
      break;
    }
    e=next;
    next=next->next;
  }

  if (!saved) {
    e->next=grib_context_malloc_clear(c,sizeof(grib_error));
    e->next->count=1;
    e->next->key=grib_context_strdup(c,key);
  }
    
  
}

static int compare_values(grib_runtime_options* options,grib_handle* h1,grib_handle *h2,const char *name,int type) {
  size_t len1 = 0;
  size_t len2 = 0;
  int err=0,i=0;
  int err1;
  int err2;
  int type1,type2;
  int countdiff;
  int isangle=0;
  int isMissing1,isMissing2;

  char *sval1 = NULL,*sval2 = NULL;
  unsigned char *uval1 = NULL,*uval2 = NULL;
  double *dval1 = NULL, *dval2 = NULL;
  long *lval1 = NULL, *lval2 = NULL;
  int failed=0;
  double maxdiff=0;
  double packingError1,packingError2;
  double value_tolerance=0;
  grib_context* c=h1->context;


  type1=type;
  type2=type;
  if (verbose) printf("  comparing %s",name);
  
  if( type1==GRIB_TYPE_UNDEFINED && (err = grib_get_native_type(h1,name,&type1)) != GRIB_SUCCESS)
  {
    printInfo(h1);
    printf("Oops... cannot get type of [%s] in 1st field: %s\n",name,grib_get_error_message(err));
    save_error(c,name);
    return err;
  }

  if(type2==GRIB_TYPE_UNDEFINED && (err = grib_get_native_type(h2,name,&type2)) != GRIB_SUCCESS)
  {
    if(err == GRIB_NOT_FOUND)
    {
      printInfo(h1);
      printf("[%s] not found in 2nd field\n",name);
      save_error(c,name);
      return err;
    }
    printInfo(h1);
    printf("Oops... cannot get type of [%s] in 2nd field: %s\n",name,grib_get_error_message(err));
    save_error(c,name);
    return err;
  }

  /*
  if(type1 != type2)
  {
    printInfo(h1);
    printf("Warning, [%s] has different types: 1st field: [%s], 2nd field: [%s]\n",
        name,grib_get_type_name(type1),grib_get_type_name(type2));
    return GRIB_TYPE_MISMATCH; 
  }
  */

  if(type1 == GRIB_TYPE_LABEL)
    return err;

  if(type1 == GRIB_TYPE_SECTION)
    return err;


  if((err = grib_get_size(h1,name,&len1)) != GRIB_SUCCESS)
  {
    printInfo(h1);
    printf("Oops... cannot get size of [%s] in 1st field: %s\n",name,grib_get_error_message(err));
    save_error(c,name);
    return err;
  }

  if((err = grib_get_size(h2,name,&len2)) != GRIB_SUCCESS)
  {
    if(err == GRIB_NOT_FOUND)
    {
      printInfo(h1);
      printf("[%s] not found in 2nd field\n",name);
      save_error(c,name);
      return err;
    }

    printInfo(h1);
    printf("Oops... cannot get size of [%s] in 2nd field: %s\n",name,grib_get_error_message(err));
    save_error(c,name);
    return err;
  }

  /*
  if(len1 != len2 && type1 != GRIB_TYPE_STRING)
  {
    printInfo(h1);
    printf("[%s] has different size: 1st field: %ld, 2nd field: %ld\n",name,(long)len1,(long)len2);
    save_error(c,name);
    return GRIB_COUNT_MISMATCH;
  }
  */

  isMissing1= ( (grib_is_missing(h1,name,&err1)==1) && (err1 == 0) ) ? 1 : 0;
  isMissing2= ( (grib_is_missing(h2,name,&err2)==1) && (err2 == 0) ) ? 1 : 0;

  if ((isMissing1==1) && (isMissing2==1)) {
    if (verbose) printf(" is set to missing in both fields\n");
    return GRIB_SUCCESS;
  }

  if (isMissing1==1) {
    if (verbose) printf(" is set to missing in 1st field\n");
    printInfo(h1);
    printf("%s is set to missing in 1st field is not missing in 2nd field\n",name);
    err1 = GRIB_VALUE_MISMATCH;
    save_error(c,name);
    return GRIB_VALUE_MISMATCH;
  }

  if (isMissing2==1) {
    if (verbose) printf(" is set to missing in 1st field\n");
    printInfo(h1);
    printf("%s is set to missing in 2nd field is not missing in 1st field\n",name);
    err1 = GRIB_VALUE_MISMATCH;
    save_error(c,name);
    return GRIB_VALUE_MISMATCH;
  }

  switch(type1)
  {
    case GRIB_TYPE_STRING:
      if (verbose) printf(" as string\n");
      len1=512;
      len2=512;
      sval1 = grib_context_malloc(h1->context,len1*sizeof(char));
      sval2 = grib_context_malloc(h2->context,len2*sizeof(char));

      if((err1 = grib_get_string(h1,name,sval1,&len1)) != GRIB_SUCCESS)
      {
        printInfo(h1);
        printf("Oops... cannot get string value of [%s] in 1st field: %s\n",
          name,grib_get_error_message(err1));
        save_error(c,name);
      }

      if((err2 = grib_get_string(h2,name,sval2,&len2)) != GRIB_SUCCESS)
      {
        printInfo(h1);
        printf("Oops... cannot get string value of [%s] in 2nd field: %s\n",
          name,grib_get_error_message(err2));
        save_error(c,name);
      }

      if(err1 == GRIB_SUCCESS && err2 == GRIB_SUCCESS)
      {
        if(grib_inline_strcmp(sval1,sval2) != 0)
        {
          printInfo(h1);
          printf("string [%s]: [%s] != [%s]\n",
            name,sval1,sval2);
          err1 = GRIB_VALUE_MISMATCH;
          save_error(c,name);
        }
      }

      grib_context_free(h1->context,sval1);
      grib_context_free(h2->context,sval2);

      if(err1) return err1;
      if(err2) return err2;

      break;

    case GRIB_TYPE_LONG:
      if (verbose) printf(" as long\n");
      
      lval1 = grib_context_malloc(h1->context,len1*sizeof(long));
      lval2 = grib_context_malloc(h2->context,len2*sizeof(long));

      if((err1 = grib_get_long_array(h1,name,lval1,&len1)) != GRIB_SUCCESS)
      {
        printInfo(h1);
        printf("Oops... cannot get long value of [%s] in 1st field: %s\n",
          name,grib_get_error_message(err1));
        save_error(c,name);
      }

      if((err2 = grib_get_long_array(h2,name,lval2,&len2)) != GRIB_SUCCESS)
      {
        printInfo(h1);
        printf("Oops... cannot get long value of [%s] in 2nd field: %s\n",
          name,grib_get_error_message(err2));
        save_error(c,name);
      }

      
      if(err1 == GRIB_SUCCESS && err2 == GRIB_SUCCESS)
      {
        int i;
        countdiff=0;
        for(i = 0; i < len1; i++)
          if(lval1[i] != lval2[i])  countdiff++;

        if (countdiff) {
          printInfo(h1);
          save_error(c,name);
          err1 = GRIB_VALUE_MISMATCH;
          if(len1 == 1)
            printf("long [%s]: [%ld] != [%ld]\n",
                   name,*lval1,*lval2);
          else
            printf("long [%s] %d out of %ld different\n",
                   name,countdiff,(long)len1);
        }
      }


      grib_context_free(h1->context,lval1);
      grib_context_free(h2->context,lval2);

      if(err1) return err1;
      if(err2) return err2;
      break;

    case GRIB_TYPE_DOUBLE:
      if (verbose) printf(" as double");
      dval1 = grib_context_malloc(h1->context,len1*sizeof(double));
      dval2 = grib_context_malloc(h2->context,len2*sizeof(double));

      isangle=0;
      value_tolerance=global_tolerance;
      if (!grib_inline_strcmp(name,"packedValues") || !grib_inline_strcmp(name,"values")
           || !grib_inline_strcmp(name,"codedValues")) {
        packingError1=0;
        packingError2=0;
        err1=grib_get_double(h1,"packingError",&packingError1);
        err2=grib_get_double(h2,"packingError",&packingError2);
        if (packingCompare)
          value_tolerance = packingError1 > packingError2 ? packingError1 : packingError2;
      } else if (!grib_inline_strcmp(name,"unpackedValues") ) {
        packingError1=0;
        packingError2=0;
        err1=grib_get_double(h1,"unpackedError",&packingError1);
        err2=grib_get_double(h2,"unpackedError",&packingError2);
        if (packingCompare)
          value_tolerance = packingError1 > packingError2 ? packingError1 : packingError2;
      } else if ( !grib_inline_rstrcmp(name,"InDegrees")) {
          packingError1=0.0005;
          packingError2=0.0005;
          isangle=1;
          value_tolerance = packingError1 > packingError2 ? packingError1 : packingError2;
      } else if (!grib_inline_strcmp(name,"referenceValue") ) {
        packingError1=0;
        packingError2=0;
        err1=grib_get_double(h1,"referenceValueError",&packingError1);
        err2=grib_get_double(h2,"referenceValueError",&packingError2);
        value_tolerance = packingError1 > packingError2 ? packingError1 : packingError2;
      }

      if (!compareAbsolute) {
         for (i=0;i<options->tolerance_count;i++) {
            if (!strcmp((options->tolerance[i]).name,name)) {
              value_tolerance=(options->tolerance[i]).double_value;
              break;
            }
         }
      }

      if((err1 = grib_get_double_array(h1,name,dval1,&len1)) != GRIB_SUCCESS)
      {
        printInfo(h1);
        printf("Oops... cannot get double value of [%s] in 1st field: %s\n",
          name,grib_get_error_message(err1));
        save_error(c,name);
      }

      if((err2 = grib_get_double_array(h2,name,dval2,&len2)) != GRIB_SUCCESS)
      {
        printInfo(h1);
        printf("Oops... cannot get double value of [%s] in 2nd field: %s\n",
          name,grib_get_error_message(err2));
        save_error(c,name);
      }
      
      if(err1 == GRIB_SUCCESS && err2 == GRIB_SUCCESS && len1!=len2)
	  {
	  	printInfo(h1);
		printf("Different size for \"%s\"  [%ld]  [%ld]\n",name,(long)len1,(long)len2);
		save_error(c,name);
	  }
      if(err1 == GRIB_SUCCESS && err2 == GRIB_SUCCESS && len1==len2)
      {
        int i,imaxdiff;
        double diff;
        double *pv1,*pv2,dnew1,dnew2;
        maxdiff=0;
        imaxdiff=0;
        countdiff=0;
        pv1=dval1;
        pv2=dval2;
        if (isangle) {
          dnew1=*dval1; dnew2=*dval2;
          pv1=&dnew1; pv2=&dnew2;
          if (*dval1 < 0 ) dnew1 += 360.0 ;
          if (*dval2 < 0 ) dnew2 += 360.0 ;
          if (*dval1 > 360 ) dnew1 -= 360.0 ;
          if (*dval2 > 360 ) dnew2 -= 360.0 ;
        }
        value_tolerance*=tolerance_factor;
        if (verbose) printf("  (%d values) tolerance=%g\n",(int)len1,value_tolerance);
        for(i = 0; i < len1; i++) {
          if((diff=compare_double(pv1++,pv2++,&value_tolerance))!=0) {
            failed=1;
            countdiff++;
            if (maxdiff < diff) {maxdiff=diff;imaxdiff=i;}
            err1 = GRIB_VALUE_MISMATCH;
          }
        }

        if (countdiff) {
          printInfo(h1);
          save_error(c,name);
          if (len1>1) {
            printf("double [%s]: %d out of %ld different, ",name,countdiff,(long)len1);
            if (compareAbsolute) printf(" max");
            printf(" absolute diff. = %g,",fabs(dval1[imaxdiff]-dval2[imaxdiff]));
            if (!compareAbsolute) printf(" max");
            printf(" relative diff. = %g",relative_error(dval1[imaxdiff],dval2[imaxdiff],value_tolerance));
            printf("\n\tmax diff. element %d: %.20e %.20e",
                   imaxdiff,dval1[imaxdiff],dval2[imaxdiff]);
            printf("\n\ttolerance=%g",value_tolerance);
            if (packingError2!=0 || packingError1!=0)
              printf(" packingError: [%g] [%g]",packingError1,packingError2);

            
            if (!grib_inline_strcmp(name,"packedValues") || !grib_inline_strcmp(name,"values")
                 || !grib_inline_strcmp(name,"codedValues")) {
              double max1,min1,max2,min2;
              grib_get_double(h1,"max",&max1);
              grib_get_double(h1,"min",&min1);
              grib_get_double(h2,"max",&max2);
              grib_get_double(h2,"min",&min2);
              printf("\n\tvalues max= [%g]  [%g]         min= [%g] [%g]",max1,max2,min1,min2);
            }
            printf("\n");
          } else {
            printf("double [%s]: [%.20e] != [%.20e]\n",
                   name,dval1[0],dval2[0]);
            printf("\tabsolute diff. = %g,",fabs(dval1[0]-dval2[0]));
            printf(" relative diff. = %g\n",relative_error(dval1[0],dval2[0],value_tolerance));
            printf("\ttolerance=%g\n",value_tolerance);
          }
        }
      }

      grib_context_free(h1->context,dval1);
      grib_context_free(h2->context,dval2);

      if(err1) return err1;
      if(err2) return err2;
      break;

    case GRIB_TYPE_BYTES:
      if (verbose) printf(" as bytes\n");
      if (len1==0) len1=512;
      if (len2==0) len2=512;
      uval1 = grib_context_malloc(h1->context,len1*sizeof(unsigned char));
      uval2 = grib_context_malloc(h2->context,len2*sizeof(unsigned char));

      if((err1 = grib_get_bytes(h1,name,uval1,&len1)) != GRIB_SUCCESS)
      {
        printInfo(h1);
        save_error(c,name);
        printf("Oops... cannot get bytes value of [%s] in 1st field: %s\n",
          name,grib_get_error_message(err1));
      }

      if((err2 = grib_get_bytes(h2,name,uval2,&len2)) != GRIB_SUCCESS)
      {
        printInfo(h1);
        save_error(c,name);
        printf("Oops... cannot get bytes value of [%s] in 2nd field: %s\n",
          name,grib_get_error_message(err2));
      }

      if(err1 == GRIB_SUCCESS && err2 == GRIB_SUCCESS)
      {
        if(memcmp(uval1,uval2,len1) != 0)
        {
        int i;
        for(i = 0; i < len1; i++)
          if(uval1[i] != uval2[i])
          {
            printInfo(h1);
            save_error(c,name);
              if(len1 == 1)
                printf("[%s] byte values are different: [%02x] and [%02x]\n",
                  name,uval1[i],uval2[i]);
              else
                printf("[%s] byte value %d of %ld are different: [%02x] and [%02x]\n",
                  name,i,(long)len1,uval1[i],uval2[i]);

            err1 = GRIB_VALUE_MISMATCH;
            break;
          }
          err1 = GRIB_VALUE_MISMATCH;
        }
      }

      grib_context_free(h1->context,uval1);
      grib_context_free(h2->context,uval2);

      if(err1) return err1;
      if(err2) return err2;
      break;

    case GRIB_TYPE_LABEL:
      if (verbose) printf(" as label\n");
      break;

    default:
      if (verbose) printf("\n");
      printInfo(h1);
      save_error(c,name);
      printf("Cannot compare [%s], unsupported type %d\n",name,type1);
      return GRIB_UNABLE_TO_COMPARE_ACCESSORS;
      break;
  }

  return GRIB_SUCCESS;

}


static int compare_handles(grib_handle* h1,grib_handle* h2,grib_runtime_options* options)
{
  int err = 0;
  int i=0;
  const char* name=NULL;
  grib_keys_iterator* iter  = NULL;

  /* mask only if no -c option or headerMode (-H)*/
  if (blacklist && ( !listFromCommandLine || headerMode )) {
    grib_string_list* nextb=blacklist;
    while (nextb) {
      grib_clear(h1,nextb->value);
      grib_clear(h2,nextb->value);
      nextb=nextb->next;
    }
  }

  if (headerMode) {
    const void *msg1=NULL,*msg2=NULL;
    size_t size1=0,size2=0;
    grib_handle *h11, *h22;
    GRIB_CHECK_NOLINE(grib_get_message_headers(h1,&msg1,&size1),0);
    GRIB_CHECK_NOLINE(grib_get_message_headers(h2,&msg2,&size2),0);
    if (size1==size2 && !memcmp(msg1,msg2,size1))
      return 0;

    err=0;
    h11=grib_handle_new_from_partial_message(h1->context,(void*)msg1,size1);
    h22=grib_handle_new_from_partial_message(h1->context,(void*)msg2,size2);

    iter=grib_keys_iterator_new(h11,
                                GRIB_KEYS_ITERATOR_SKIP_COMPUTED,NULL);

    if (!iter) {
      printf("ERROR: unable to get iterator\n");
      exit(1);
    }

    while(grib_keys_iterator_next(iter))
    {
      name=grib_keys_iterator_get_name(iter);
      /*printf("----- comparing %s\n",name);*/

      if (blacklisted(name)) continue;
      if(compare_values(options,h11,h22,name,GRIB_TYPE_UNDEFINED))  err++;
    }

    grib_keys_iterator_delete(iter);
    grib_handle_delete(h11);
    grib_handle_delete(h22);
    return err;
  }

  if ( listFromCommandLine && onlyListed ) {
    for (i=0; i< options->compare_count; i++) {
	  if (blacklisted((char*)options->compare[i].name)) continue;
      if (options->compare[i].type == GRIB_NAMESPACE) {
        iter=grib_keys_iterator_new(h1,0,(char*)options->compare[i].name);
        if (!iter) {
          printf("ERROR: unable to get iterator\n");
          exit(1);
        }
        while(grib_keys_iterator_next(iter))
        {
          name=grib_keys_iterator_get_name(iter);
          /*printf("----- comparing %s\n",name);*/

          if (blacklisted(name)) continue;
          if(compare_values(options,h1,h2,name,GRIB_TYPE_UNDEFINED))  err++;
        }
        grib_keys_iterator_delete(iter);
      } else {
        if( compare_values(options,h1,h2,options->compare[i].name,options->compare[i].type))
          err++;
      }
    }
  } else {
    const void *msg1=NULL,*msg2=NULL;
    size_t size1=0,size2=0;
    GRIB_CHECK_NOLINE(grib_get_message(h1,&msg1,&size1),0);
    GRIB_CHECK_NOLINE(grib_get_message(h2,&msg2,&size2),0);
    if (size1==size2 && !memcmp(msg1,msg2,size1))
       return 0;

    iter=grib_keys_iterator_new(h1,GRIB_KEYS_ITERATOR_SKIP_COMPUTED,NULL);

    if (!iter) {
      printf("ERROR: unable to get iterator\n");
      exit(1);
    }

    while(grib_keys_iterator_next(iter))
    {
      name=grib_keys_iterator_get_name(iter);
      /*printf("----- comparing %s\n",name);*/

      if (blacklisted(name)) continue;
      if(compare_values(options,h1,h2,name,GRIB_TYPE_UNDEFINED))  err++;
    }

    grib_keys_iterator_delete(iter);

    if ( listFromCommandLine ) {
      for (i=0; i< options->compare_count; i++) {
        if (blacklisted(name)) continue;
        if (options->compare[i].type == GRIB_NAMESPACE) {
          iter=grib_keys_iterator_new(h1,0,(char*)options->compare[i].name);
          if (!iter) {
            printf("ERROR: unable to get iterator for %s\n",options->compare[i].name );
            exit(1);
          }
          while(grib_keys_iterator_next(iter))
          {
            name=grib_keys_iterator_get_name(iter);
            /*printf("----- comparing %s\n",name);*/

            if (blacklisted(name)) continue;
            if(compare_values(options,h1,h2,name,GRIB_TYPE_UNDEFINED))  err++;
          }
          grib_keys_iterator_delete(iter);
        } else {
          if( compare_values(options,h1,h2,options->compare[i].name,options->compare[i].type))
            err++;
        }
      }
    } 
    
  }
  return err;
}
