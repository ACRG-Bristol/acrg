#include "grib_api.h"

void usage(char* prog) {
	printf("usage: %s in.nc\n",prog);
	exit(1);
}

int main(int argc,char* argv[]) {

	char* file;
	int err=0;
	grib_handle* h;
	char identifier[7]={0,};
	size_t len=7;
	grib_context* c=grib_context_get_default();
	 
	if (argc>2) usage(argv[0]);

	file=argv[1];

	h=grib_handle_new_from_nc_file(c,file,&err);
	grib_get_string(h,"identifier",identifier,&len);
	printf("%s\n",identifier);
	GRIB_CHECK(err,0);

	return err;


}
