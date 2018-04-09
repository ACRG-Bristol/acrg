#include "grib_api.h"

void usage(char* prog) {
	printf("usage: %s in1.grib in2.grib out.grib\n",prog);
	exit(1);
}

int main ( int argc, char* argv[]) {
	grib_handle *hfrom,*hto,*h;
	FILE *in;
	char *in_name,*in_name1,*out_name;
	int err=0;
	int what=0;

	if (argc<4) usage(argv[0]);

	in_name=argv[1];
	in_name1=argv[2];
	out_name=argv[3];

	in=fopen(in_name,"r");
	if (!in) {
		perror(in_name);
		exit(1);
	}

	hfrom=grib_handle_new_from_file(0,in,&err);
	GRIB_CHECK(err,0);
	fclose(in);

	in=fopen(in_name1,"r");
	if (!in) {
		perror(in_name1);
		exit(1);
	}

	hto=grib_handle_new_from_file(0,in,&err);
	GRIB_CHECK(err,0);
	fclose(in);

	what=GRIB_SECTION_PRODUCT | GRIB_SECTION_LOCAL;
	h=grib_util_sections_copy(hfrom,hto,what,&err);
	GRIB_CHECK(err,0);

	err=grib_write_message(h,out_name,"w");

	return err;
}
