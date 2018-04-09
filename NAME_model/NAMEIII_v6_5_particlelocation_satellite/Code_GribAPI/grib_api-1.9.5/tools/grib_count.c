#include "grib_api_internal.h"

void usage(char* prog) {
  printf("usage: %s infile1 infile2 ... \n",prog);
  exit(1);
}

static int check_file(FILE* in,long *count) {
  void* mesg=NULL;
  size_t size=0;
  int err=0;
  grib_context* c=grib_context_get_default();

  if (!in) return 1;

  while ( (err=grib_read_any_from_file_alloc (c, in,  &mesg , &size))==GRIB_SUCCESS) {
      grib_context_free(c,mesg);
     (*count)++;
  }

  if (err==GRIB_END_OF_FILE) err=GRIB_SUCCESS;
  
  return err;
}

int main(int argc,char* argv[]) {

  FILE* infh;
  char* filename;
  int i;
  int err=0;
  long n=0,nn=0;

  if (argc <2) usage(argv[0]);

  n=0;
  for (i=1;i<argc;i++) {
    filename=argv[i];

	infh=fopen(filename,"r");
	if (!infh) {
	  perror(filename);
	  exit(1);
	}

    nn=0;
    err=check_file(infh,&nn);
	if (err!=0) exit(err);
    n+=nn;

	fclose(infh);
  }
  printf("%ld\n",n);

  return 0;
}
