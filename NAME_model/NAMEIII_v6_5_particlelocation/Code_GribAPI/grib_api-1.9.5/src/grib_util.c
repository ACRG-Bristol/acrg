/**
 * Copyright 2005-2007 ECMWF
 *
 * Licensed under the GNU Lesser General Public License which
 * incorporates the terms and conditions of version 3 of the GNU
 * General Public License.
 * See LICENSE and gpl-3.0.txt for details.
 */

/***************************************************************************
 *   Enrico Fucile
 *                                                                         *
 ***************************************************************************/
#include "grib_api_internal.h"


static void set_total_length(unsigned char* buffer,long *section_length,long *section_offset,int edition,size_t totalLength) {
	long off;
	switch (edition) {
		case 1:
			if(totalLength < 0x800000 ) {
				off=32;
				grib_encode_unsigned_long(buffer, (unsigned long)totalLength ,  &off, 24);
			} else {
				long s4len,t120;
				totalLength -= 4;
				t120  = (totalLength+119)/120;
				s4len  = t120*120 - totalLength;
				totalLength  = 0x800000 | t120;
				off=32;
				grib_encode_unsigned_long(buffer, (unsigned long)totalLength ,  &off, 24);
				off=section_offset[4]*8;
				grib_encode_unsigned_long(buffer, (unsigned long)s4len ,  &off, 24);
			}
			break;
		case 2:
			off=64;
			grib_encode_unsigned_long(buffer, (unsigned long)totalLength ,  &off, 64);
			break;
	}

}

static grib_handle* grib_sections_copy_internal(grib_handle* hfrom,grib_handle* hto,int sections[],int *err) {
	int i;
	size_t totalLength=0;
	unsigned char* buffer;
	unsigned char *p;
	long edition=0;
	long section_length[MAX_NUM_SECTIONS]={0,};
	long section_offset[MAX_NUM_SECTIONS]={0,};
	long off=0;
	grib_handle* h;
	char section_length_str[]="section0Length";
	char section_offset_str[]="offsetSection0";
	long length,offset;

	*err=grib_get_long(hfrom,"edition",&edition);
	if (*err) return NULL;
	
	for (i=0;i<=hfrom->sections_count;i++) {

		if (sections[i]) {h=hfrom;}
		else {h=hto;}

		sprintf(section_length_str,"section%dLength",i);
		if (grib_get_long(h,section_length_str,&length)) continue;
		section_length[i]=length;

		sprintf(section_offset_str,"offsetSection%d",i);
		if (grib_get_long(h,section_offset_str,&offset)) continue;
		section_offset[i]=offset;

		totalLength+=section_length[i];

	}

	buffer=grib_context_malloc_clear(hfrom->context,totalLength*sizeof(char));

	p=buffer;
	off=0;
	for (i=0;i<=hfrom->sections_count;i++) {
		grib_handle* h;
        if (sections[i]) h=hfrom;
		else h=hto;
		p=memcpy(p,h->buffer->data+section_offset[i],section_length[i]);
		section_offset[i]=off;
		off+=section_length[i];
		p+=section_length[i];
   }

   /* copy section 3 present flag*/
   if (edition==1) {
		const void* buffer_to=NULL;
		size_t size_to=0;
		grib_get_message(hto,&buffer_to,&size_to);
		memcpy(buffer+15,((unsigned char*)buffer_to)+15,1);
	}

   set_total_length(buffer,section_length,section_offset,edition,totalLength);

   h=grib_handle_new_from_message(hfrom->context,buffer,totalLength);

   switch (edition) {
   	case 1:
		if (sections[1] && sections[2]) break;
		
		if (sections[1]) {
			long PVPresent;
			grib_get_long(hfrom,"PVPresent",&PVPresent);
			if (PVPresent) {
				double *pv;
				long numberOfVerticalCoordinateValues;
				size_t size=0;

				grib_get_long(hfrom,"numberOfVerticalCoordinateValues",&numberOfVerticalCoordinateValues);
				size=numberOfVerticalCoordinateValues;
				pv=grib_context_malloc_clear(hfrom->context,numberOfVerticalCoordinateValues*sizeof(double));
				grib_get_double_array(hfrom,"pv",pv,&size);
				grib_set_long(h,"PVPresent",1);
				grib_set_double_array(h,"pv",pv,size);

				grib_context_free(hfrom->context,pv);

			} else {
				grib_set_long(h,"PVPresent",0);
			}
		}
		if (sections[2]) {
			long PVPresent;
			grib_get_long(hto,"PVPresent",&PVPresent);
			if (PVPresent) {
				double *pv;
				long numberOfVerticalCoordinateValues;
				size_t size=0;

				grib_get_long(hto,"numberOfVerticalCoordinateValues",&numberOfVerticalCoordinateValues);
				size=numberOfVerticalCoordinateValues;
				pv=grib_context_malloc_clear(hto->context,numberOfVerticalCoordinateValues*sizeof(double));
				grib_get_double_array(hto,"pv",pv,&size);
				grib_set_long(h,"PVPresent",1);
				grib_set_double_array(h,"pv",pv,size);

				grib_context_free(hto->context,pv);

			} else {
				grib_set_long(h,"PVPresent",0);
			}
		}
		break;
	case 2:
		if (sections[1]) {
			long discipline;
			grib_get_long(hfrom,"discipline",&discipline);
			grib_set_long(h,"discipline",discipline);
		}
		break;
   }

   return h;
}

grib_handle* grib_util_sections_copy(grib_handle* hfrom,grib_handle* hto,int what,int *err) {
	long edition_from=0;
	long edition_to=0;
	long localDefinitionNumber=-1;
	int sections_to_copy[MAX_NUM_SECTIONS]={0,};

	*err=grib_get_long(hfrom,"edition",&edition_from);
	if (*err) return NULL;
	*err=grib_get_long(hto,"edition",&edition_to);
	if (*err) return NULL;

	if (edition_to != 1 && edition_to != 2 ) {
		*err=GRIB_NOT_IMPLEMENTED;
		return NULL;
	}

	if (edition_from!=edition_to) {
		*err=GRIB_DIFFERENT_EDITION;
		return NULL;
	}

	if (what & GRIB_SECTION_GRID) {
		switch (edition_from) {
			case 1:
				sections_to_copy[2]=1;
				break;
			case 2:
				sections_to_copy[3]=1;
				break;
		}
	}

	if (what & GRIB_SECTION_DATA) {
		switch (edition_from) {
			case 1:
				sections_to_copy[3]=1;
				sections_to_copy[4]=1;
				break;
			case 2:
				sections_to_copy[5]=1;
				sections_to_copy[6]=1;
				sections_to_copy[7]=1;
				break;
		}
	}

	if (what & GRIB_SECTION_LOCAL) {
		switch (edition_from) {
			case 1:
				sections_to_copy[1]=1;
				break;
			case 2:
				sections_to_copy[2]=1;
				break;
		}
	}

	if (what & GRIB_SECTION_PRODUCT) {
		switch (edition_from) {
			case 1:
				grib_get_long(hfrom,"localDefinitionNumber",&localDefinitionNumber);
				if (localDefinitionNumber==13) {
					sections_to_copy[4]=1;
				}
				sections_to_copy[1]=1;
				break;
			case 2:
				sections_to_copy[1]=1;
				sections_to_copy[4]=1;
				break;
		}
	}

	if (what & GRIB_SECTION_BITMAP) {
		switch (edition_from) {
			case 1:
				sections_to_copy[3]=1;
				break;
			case 2:
				sections_to_copy[6]=1;
				break;
		}
	}

	return grib_sections_copy_internal(hfrom,hto,sections_to_copy,err);
	
}

static grib_trie* init_list(const char* name);
static grib_trie* param_id_list = NULL;
static grib_trie* mars_param_list = NULL;
/* TODO thread safe */
grib_string_list* grib_util_get_param_id(const char* mars_param) {
	if (!mars_param_list && (mars_param_list=init_list("mars_param.table"))==NULL) return NULL;
	return grib_trie_get(mars_param_list,mars_param);
}

grib_string_list* grib_util_get_mars_param(const char* param_id) {
	if (!param_id_list && (param_id_list=init_list("param_id.table"))==NULL) return NULL;
	return grib_trie_get(param_id_list,param_id);
}

static grib_trie* init_list(const char* name) {
	char *full_path=0;
	FILE* fh;
	char s[100];
	char param[100];
	grib_string_list* list=0;
	grib_string_list* next=0;
	grib_trie* trie_list;
	grib_context* c=grib_context_get_default();
	full_path=grib_context_full_path(c,name);

	fh=fopen(full_path,"r");
	if (!fh) {
		grib_context_log(c,GRIB_LOG_PERROR,"unable to read %s",full_path);
		return NULL;
	}

	list=grib_context_malloc_clear(c,sizeof(grib_string_list));
	trie_list=grib_trie_new(c);
	if (fscanf(fh,"%s",param)==EOF) return NULL;
	while (fscanf(fh,"%s",s)!=EOF) {
		if (!strcmp(s,"|")) {
			grib_trie_insert(trie_list, param,list);
			if (fscanf(fh,"%s",param)==EOF) return trie_list;
			list=NULL;
		} else {
			if (!list) {
				list=grib_context_malloc_clear(c,sizeof(grib_string_list));
				list->value=grib_context_strdup(c,s);
			} else {
				next=list;
				while(next->next) next=next->next;
				next->next=grib_context_malloc_clear(c,sizeof(grib_string_list));
				next->next->value=grib_context_strdup(c,s);
			}
		}
	}

	fclose(fh);

	return 0;
}

static void print_values(const grib_util_grid_spec* spec,const double* data_values,size_t data_values_count,const grib_values *values,int count) {
		int i;
		printf("\n");
		printf("XXXXXX => grib_set_values: setting %d values \n",count);

		for(i = 0; i < count ; i++)
		{
			printf("%s = ",values[i].name);
			switch(values[i].type)
			{
				case GRIB_TYPE_LONG: printf(" %ld;",(long)values[i].long_value); break;
				case GRIB_TYPE_DOUBLE: printf(" %g;",values[i].double_value); break;
				case GRIB_TYPE_STRING: printf(" \"%s\";",values[i].string_value); break;
			}
			printf("\n");
		}
		printf("\n");


		if(spec->bitmapPresent) {
			int missing = 0;
			double min = 1e100;
			for(i = 0; i < data_values_count ; i++)
			{
				double d = data_values[i] - spec->missingValue;
				if(d < 0) d = -d;

				if(d < min) {
					min = d;
				}
				
				if(data_values[i] == spec->missingValue)
					missing++;


			}

			printf("MISSING VALUES : %f %d out of %d (min = %f)\n", spec->missingValue, missing, data_values_count, min);

		}
}

#define     ISECTION_2  3000
#define     ISECTION_4  512

grib_handle* grib_util_set_spec(grib_handle* h, 
		const grib_util_grid_spec    *spec,
		const grib_util_packing_spec *packing_spec,
		int		flags,
		const double*                data_values,
		size_t                       data_values_count,
		int* err)
{

#define SET_LONG_VALUE(n,v)   do { Assert(count<1024); values[count].name = n; values[count].type = GRIB_TYPE_LONG;   values[count].long_value = v; count++; } while(0)
#define SET_DOUBLE_VALUE(n,v) do { Assert(count<1024); values[count].name = n; values[count].type = GRIB_TYPE_DOUBLE; values[count].double_value = v; count++; } while(0)
#define SET_STRING_VALUE(n,v) do { Assert(count<1024); values[count].name = n; values[count].type = GRIB_TYPE_STRING;   values[count].string_value = v; count++; } while(0)

#define COPY_SPEC_LONG(x)     do { Assert(count<1024); values[count].name = #x; values[count].type = GRIB_TYPE_LONG;   values[count].long_value = spec->x; count++; } while(0)
#define COPY_SPEC_DOUBLE(x)   do { Assert(count<1024); values[count].name = #x; values[count].type = GRIB_TYPE_DOUBLE; values[count].double_value = spec->x; count++; } while(0)

	grib_values  values[1024];
	size_t       count = 0;
	int i;
	long editionNumber;
	grib_handle* outh = NULL;
	grib_handle* tmp = NULL;
	const char*  grid_type = NULL;
	char name[1024];
	char input_grid_type[100];
	char input_packing_type[100];
	long input_bits_per_value=0;
	long input_decimal_scale_factor=0;
	size_t len=100;
	size_t input_grid_type_len=100;
	long P;

	static grib_util_packing_spec default_packing_spec = {0, };

	int debug = 1;

	if(!packing_spec) {
		packing_spec = &default_packing_spec;
	}

	len=100;
	grib_get_string(h,"packingType",input_packing_type,&len);
	grib_get_long(h,"bitsPerValue",&input_bits_per_value);
	grib_get_long(h,"decimalScaleFactor",&input_decimal_scale_factor);
	printf("XXXXXX => input_packing_type = %s\n",input_packing_type);
	printf("XXXXXX => input_bits_per_value = %ld\n",input_bits_per_value);
	printf("XXXXXX => input_decimal_scale_factor = %ld\n",input_decimal_scale_factor);

	if (flags & GRIB_UTIL_SET_SPEC_FLAGS_ONLY_PACKING) {
		if (packing_spec->packing == GRIB_UTIL_PACKING_USE_PROVIDED) {
			switch (packing_spec->packing_type) {
				case GRIB_UTIL_PACKING_TYPE_SPECTRAL_COMPLEX:
					if (strcmp(input_packing_type,"spectral_complex") && !strcmp(input_packing_type,"spectral_simple"))
						SET_STRING_VALUE("packingType","spectral_complex");
					break;
				case GRIB_UTIL_PACKING_TYPE_SPECTRAL_SIMPLE:
					if (strcmp(input_packing_type,"spectral_simple") && !strcmp(input_packing_type,"spectral_complex"))
						SET_STRING_VALUE("packingType","spectral_simple");
					break;
				case GRIB_UTIL_PACKING_TYPE_GRID_SIMPLE:
					if (strcmp(input_packing_type,"grid_simple") && !strcmp(input_packing_type,"grid_complex"))
						SET_STRING_VALUE("packingType","grid_simple");
					break;
				case GRIB_UTIL_PACKING_TYPE_GRID_COMPLEX:
					if (strcmp(input_packing_type,"grid_complex") && !strcmp(input_packing_type,"grid_simple"))
						SET_STRING_VALUE("packingType","grid_complex");
					break;
				case GRIB_UTIL_PACKING_TYPE_JPEG:
					if (strcmp(input_packing_type,"grid_jpeg") && !strcmp(input_packing_type,"grid_simple"))
						SET_STRING_VALUE("packingType","grid_jpeg");
					break;
				default :
					printf("invalid packing_spec->packing_type = %ld\n",(long)packing_spec->packing_type);
					abort();
			}
		}
		switch(packing_spec->accuracy) {

			case GRIB_UTIL_ACCURACY_SAME_BITS_PER_VALUES_AS_INPUT:
			break;

			case GRIB_UTIL_ACCURACY_USE_PROVIDED_BITS_PER_VALUES:
				if (input_bits_per_value!=packing_spec->bitsPerValue)
					SET_LONG_VALUE("bitsPerValue", packing_spec->bitsPerValue);
				break;

			case GRIB_UTIL_ACCURACY_SAME_DECIMAL_SCALE_FACTOR_AS_INPUT:
				break;

			case GRIB_UTIL_ACCURACY_USE_PROVIDED_DECIMAL_SCALE_FACTOR:
				if (input_decimal_scale_factor!=packing_spec->decimalScaleFactor)
					SET_LONG_VALUE("decimalScaleFactor", packing_spec->decimalScaleFactor);
				break;

			default:
				printf("invalid packing_spec->accuracy = %ld\n",(long)packing_spec->accuracy);
				abort();
		}

		if (packing_spec->deleteLocalDefinition) {
			/* TODO: not working for grib2 */
			SET_LONG_VALUE("deleteLocalDefinition",1);
		}

		/*nothing to be changed*/
		if (count==0) return h;

		if (debug) print_values(spec,data_values,data_values_count,values,count);
		if((*err = grib_set_values(h,values,count)) != 0)
		{
			fprintf(stderr,"GRIB_UTIL_SET_SPEC: Cannot set values  %s\n",grib_get_error_message(*err));

			for(i = 0; i < count; i++)
				if(values[i].error)	
					fprintf(stderr," %s %s\n",values[i].name,grib_get_error_message(values[i].error));
			goto cleanup;
		}
		if (debug) printf("XXXXXX => grib_set_double_array \n");
		{
			int j=0;
			for (j=0;j<20;j++) printf("XXXXXX => %g\n",data_values[j]);
			printf("XXXXXX => data_values_count=%d \n",(int)data_values_count);
		}
		if((*err = grib_set_double_array(h,"values",data_values,data_values_count)) != 0)
		{
			goto cleanup;
		}
		if (debug) printf("XXXXXX <= grib_set_double_array \n");

		return h;
	}

	/* Get edition number from input handle */
	if((*err = grib_get_long(h,"editionNumber",&editionNumber)) != 0)
	{
		return NULL;
	}


	switch(spec->grid_type) {

		case GRIB_UTIL_GRID_SPEC_REGULAR_LL:
			grid_type = "regular_ll";
			break;

		case GRIB_UTIL_GRID_SPEC_ROTATED_LL:
			grid_type = "rotated_ll";
			break;

		case GRIB_UTIL_GRID_SPEC_REGULAR_GG:
			grid_type = "regular_gg";
			break;

		case GRIB_UTIL_GRID_SPEC_ROTATED_GG:
			grid_type = "rotated_gg";
			break;

		case GRIB_UTIL_GRID_SPEC_REDUCED_LL:
			grid_type = "reduced_ll";
			break;

		case GRIB_UTIL_GRID_SPEC_POLAR_STEREOGRAPHIC:
			grid_type = "polar_stereographic";
			break;

		case GRIB_UTIL_GRID_SPEC_REDUCED_GG:
			grid_type = "reduced_gg";
			break;

		case GRIB_UTIL_GRID_SPEC_SH:
			grid_type = "sh";
			break;

		default:
			*err = GRIB_NOT_IMPLEMENTED;
			return NULL;

	}

	SET_STRING_VALUE("gridType", grid_type);

	/* The "pl" is given from the template, but "section_copy" will take care of setting the right headers */

	{
		char levtype[80];
		size_t n = sizeof(levtype);
		Assert(grib_get_string(h,"levelType",levtype,&n) == 0);
		switch (spec->grid_type) {
			case GRIB_UTIL_GRID_SPEC_REDUCED_GG:
				sprintf(name, "%s_pl_%ld_grib%ld", grid_type,spec->N, editionNumber);
				break;
			default :
				sprintf(name, "%s_pl_grib%ld", grid_type, editionNumber);
		}
	}

	/* TODO: recycle tmp handle */
	if (debug) printf("XXXXXX => grib_handle_new_from_template %s\n", name);
	tmp = grib_handle_new_from_samples(NULL, name);
	if(!tmp) {
		*err = GRIB_INVALID_FILE;
		return NULL;
	}
	if (debug) printf("XXXXXX <= grib_handle_new_from_template %s\n", name);

	/* Set grid  */
	switch(spec->grid_type) {
		case GRIB_UTIL_GRID_SPEC_REGULAR_LL:
		case GRIB_UTIL_GRID_SPEC_ROTATED_LL:

			COPY_SPEC_LONG  (bitmapPresent);
			if (spec->missingValue) COPY_SPEC_DOUBLE(missingValue);

			SET_LONG_VALUE  ("ijDirectionIncrementGiven",    1);

			/* default iScansNegatively=0 jScansPositively=0 is ok */
			COPY_SPEC_LONG(iScansNegatively);
			COPY_SPEC_LONG(jScansPositively);

			COPY_SPEC_LONG(Ni);
			COPY_SPEC_LONG(Nj);

			COPY_SPEC_DOUBLE(iDirectionIncrementInDegrees);
			COPY_SPEC_DOUBLE(jDirectionIncrementInDegrees);

			COPY_SPEC_DOUBLE(longitudeOfFirstGridPointInDegrees);
			COPY_SPEC_DOUBLE(longitudeOfLastGridPointInDegrees);

			COPY_SPEC_DOUBLE(latitudeOfFirstGridPointInDegrees);
			COPY_SPEC_DOUBLE(latitudeOfLastGridPointInDegrees);

			if(!fabs((spec->latitudeOfFirstGridPointInDegrees - spec->latitudeOfLastGridPointInDegrees)/spec->jDirectionIncrementInDegrees) == spec->Nj)
			{
				printf("XXXXXXXXXXXXXXXX %g %ld\n",
						(spec->latitudeOfFirstGridPointInDegrees - spec->latitudeOfLastGridPointInDegrees)/spec->jDirectionIncrementInDegrees,
						spec->Nj);
			}

			/*
			   Assert((spec->latitudeOfFirstGridPointInDegrees - spec->latitudeOfLastGridPointInDegrees)/spec->jDirectionIncrementInDegrees == spec->Nj));
			   Assert((spec->longitudeOfLastGridPointInDegrees - spec->longitudeOfFirstGridPointInDegrees)/spec->iDirectionIncrementInDegrees == spec->Ni));
			 */

			break;

		case GRIB_UTIL_GRID_SPEC_REGULAR_GG:
		case GRIB_UTIL_GRID_SPEC_ROTATED_GG:

			COPY_SPEC_LONG  (bitmapPresent);
			if (spec->missingValue) COPY_SPEC_DOUBLE(missingValue);
			SET_LONG_VALUE("ijDirectionIncrementGiven", 1);

			/* TODO: add Assert */

			COPY_SPEC_LONG(Ni);
			COPY_SPEC_DOUBLE(iDirectionIncrementInDegrees);

			COPY_SPEC_LONG(Nj);
			COPY_SPEC_LONG(N);

			/* TODO: Compute here ... */
			COPY_SPEC_DOUBLE(longitudeOfFirstGridPointInDegrees);
			COPY_SPEC_DOUBLE(longitudeOfLastGridPointInDegrees);

			COPY_SPEC_DOUBLE(latitudeOfFirstGridPointInDegrees);
			COPY_SPEC_DOUBLE(latitudeOfLastGridPointInDegrees);
			break;


		case GRIB_UTIL_GRID_SPEC_REDUCED_LL:
			COPY_SPEC_LONG  (bitmapPresent);
			if (spec->missingValue) COPY_SPEC_DOUBLE(missingValue);
			SET_LONG_VALUE("ijDirectionIncrementGiven", 0);

			COPY_SPEC_LONG(Nj);

			COPY_SPEC_DOUBLE(longitudeOfFirstGridPointInDegrees);
			COPY_SPEC_DOUBLE(longitudeOfLastGridPointInDegrees);

			COPY_SPEC_DOUBLE(latitudeOfFirstGridPointInDegrees);
			COPY_SPEC_DOUBLE(latitudeOfLastGridPointInDegrees);
			break;

		case GRIB_UTIL_GRID_SPEC_POLAR_STEREOGRAPHIC:
			COPY_SPEC_LONG  (bitmapPresent);
			if (spec->missingValue) COPY_SPEC_DOUBLE(missingValue);

			COPY_SPEC_DOUBLE(longitudeOfFirstGridPointInDegrees);
			COPY_SPEC_LONG(Ni);
			COPY_SPEC_LONG(Nj);

			/* default iScansNegatively=0 jScansPositively=0 is ok */
			COPY_SPEC_LONG(iScansNegatively);
			COPY_SPEC_LONG(jScansPositively);

			COPY_SPEC_DOUBLE(orientationOfTheGridInDegrees);

			COPY_SPEC_LONG(DxInMetres);
			COPY_SPEC_LONG(DyInMetres);

			break;

		case GRIB_UTIL_GRID_SPEC_REDUCED_GG:
			COPY_SPEC_LONG  (bitmapPresent);
			if (spec->missingValue) COPY_SPEC_DOUBLE(missingValue);
			SET_LONG_VALUE("ijDirectionIncrementGiven", 0);

			/* TODO: add Assert */
			COPY_SPEC_LONG(Nj);
			COPY_SPEC_LONG(N);

			/* TODO: Compute here ... */
			COPY_SPEC_DOUBLE(longitudeOfFirstGridPointInDegrees);
			COPY_SPEC_DOUBLE(longitudeOfLastGridPointInDegrees);

			COPY_SPEC_DOUBLE(latitudeOfFirstGridPointInDegrees);
			COPY_SPEC_DOUBLE(latitudeOfLastGridPointInDegrees);

			break;

		case GRIB_UTIL_GRID_SPEC_SH:
			*err=grib_get_string(h,"gridType",input_grid_type,&input_grid_type_len);

			SET_LONG_VALUE("J", spec->truncation);
			SET_LONG_VALUE("K", spec->truncation);
			SET_LONG_VALUE("M", spec->truncation);

			if(packing_spec->packing_type == GRIB_UTIL_PACKING_TYPE_SPECTRAL_COMPLEX)
			{
				SET_STRING_VALUE("packingType", "spectral_complex");
				SET_LONG_VALUE("JS", 20);
				SET_LONG_VALUE("KS", 20);
				SET_LONG_VALUE("MS", 20);

				if ((!(*err) && strcmp(input_grid_type,"sh")) || packing_spec->computeLaplacianOperator )
					SET_LONG_VALUE("computeLaplacianOperator", 1);
				else {
					SET_LONG_VALUE("computeLaplacianOperator", 0);
					*err=grib_get_long(h,"P",&P);
					SET_LONG_VALUE("P", P);
				}

			}

			break;
	}


	/* Set rotation */

	switch(spec->grid_type) {
		case GRIB_UTIL_GRID_SPEC_ROTATED_LL:
		case GRIB_UTIL_GRID_SPEC_ROTATED_GG:
			COPY_SPEC_LONG(uvRelativeToGrid);
			COPY_SPEC_DOUBLE(latitudeOfSouthernPoleInDegrees);
			COPY_SPEC_DOUBLE(longitudeOfSouthernPoleInDegrees);
			break;
	}

	/* process packing options */
	if(packing_spec->editionNumber) {
		SET_LONG_VALUE("editionNumber", packing_spec->editionNumber);
	}

	if (packing_spec->packing == GRIB_UTIL_PACKING_USE_PROVIDED) {
		switch (packing_spec->packing_type) {
			case GRIB_UTIL_PACKING_TYPE_SPECTRAL_COMPLEX:
				if (strcmp(input_packing_type,"spectral_complex") && !strcmp(input_packing_type,"spectral_simple"))
					SET_STRING_VALUE("packingType","spectral_complex");
				break;
			case GRIB_UTIL_PACKING_TYPE_SPECTRAL_SIMPLE:
				if (strcmp(input_packing_type,"spectral_simple") && !strcmp(input_packing_type,"spectral_complex"))
					SET_STRING_VALUE("packingType","spectral_simple");
				break;
			case GRIB_UTIL_PACKING_TYPE_GRID_SIMPLE:
				if (strcmp(input_packing_type,"grid_simple") && !strcmp(input_packing_type,"grid_complex"))
					SET_STRING_VALUE("packingType","grid_simple");
				break;
			case GRIB_UTIL_PACKING_TYPE_GRID_COMPLEX:
				if (strcmp(input_packing_type,"grid_complex") && !strcmp(input_packing_type,"grid_simple"))
					SET_STRING_VALUE("packingType","grid_complex");
				break;
			case GRIB_UTIL_PACKING_TYPE_JPEG:
				if (strcmp(input_packing_type,"grid_jpeg") && !strcmp(input_packing_type,"grid_simple"))
					SET_STRING_VALUE("packingType","grid_jpeg");
				break;
			default :
				printf("invalid packing_spec->packing_type = %ld\n",(long)packing_spec->packing_type);
				abort();
		}
	}

	switch(packing_spec->accuracy) {

		case GRIB_UTIL_ACCURACY_SAME_BITS_PER_VALUES_AS_INPUT:
			{
				long bitsPerValue = 0;
				Assert(grib_get_long(h, "bitsPerValue", &bitsPerValue) == 0);
				SET_LONG_VALUE("bitsPerValue", bitsPerValue);
			}
			break;

		case GRIB_UTIL_ACCURACY_USE_PROVIDED_BITS_PER_VALUES:
			SET_LONG_VALUE("bitsPerValue", packing_spec->bitsPerValue);
			break;

		case GRIB_UTIL_ACCURACY_SAME_DECIMAL_SCALE_FACTOR_AS_INPUT:
			{
				long decimalScaleFactor = 0;
				Assert(grib_get_long(h, "decimalScaleFactor", &decimalScaleFactor) == 0);
				SET_LONG_VALUE("decimalScaleFactor", decimalScaleFactor);
			}
			break;

		case GRIB_UTIL_ACCURACY_USE_PROVIDED_DECIMAL_SCALE_FACTOR:
			SET_LONG_VALUE("decimalScaleFactor", packing_spec->decimalScaleFactor);
			break;



		default:
			printf("invalid packing_spec->accuracy = %ld\n",(long)packing_spec->accuracy);
			abort();
	}

	if(packing_spec->extra_settings_count) {
		for(i = 0; i < packing_spec->extra_settings_count; i++) {
			Assert(count < 1024);
			values[count++] = packing_spec->extra_settings[i];
		}
	}
	if (packing_spec->deleteLocalDefinition) {
		/* TODO: not working for grib2 */
		SET_LONG_VALUE("deleteLocalDefinition",1);
	}

	/* grib_write_message(h,"input.grib","w"); */
	/* grib_write_message(tmp,"geo.grib","w"); */
	if (debug) printf("XXXXXX => grib_utils_sections_copy %s\n", name);
	if((outh = grib_util_sections_copy(h, tmp, GRIB_SECTION_PRODUCT | GRIB_SECTION_LOCAL,err)) == NULL)
	{
		goto cleanup;
	}
	Assert(*err == 0);
	if (debug) printf("XXXXXX <= grib_utils_sections_copy %s\n", name);

	if (debug) print_values(spec,data_values,data_values_count,values,count);
	if((*err = grib_set_values(outh,values,count)) != 0)
	{
		fprintf(stderr,"SET_GRID_DATA_DESCRIPTION: Cannot set values  %s\n",grib_get_error_message(*err));

		for(i = 0; i < count; i++)
			if(values[i].error)	
				fprintf(stderr," %s %s\n",values[i].name,grib_get_error_message(values[i].error));
		goto cleanup;
	}
	if (spec->pl_size!=0 && spec->grid_type==GRIB_UTIL_GRID_SPEC_REDUCED_GG) {
		grib_set_long_array(outh,"pl",spec->pl,spec->pl_size);
	}

	if (debug) printf("XXXXXX <= grib_set_values \n");

	if (debug) printf("XXXXXX => grib_set_double_array \n");
	if((*err = grib_set_double_array(outh,"values",data_values,data_values_count)) != 0)
	{
		goto cleanup;
	}
	if (debug) printf("XXXXXX <= grib_set_double_array \n");
	fflush(stdout);
	/* grib_write_message(outh,"h.grib","w"); */

	return outh;

cleanup:
	if(outh)  grib_handle_delete(outh);
	return NULL;
}
