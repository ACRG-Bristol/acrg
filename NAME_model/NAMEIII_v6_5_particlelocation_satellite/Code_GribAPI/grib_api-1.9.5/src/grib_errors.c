/* This file is automatically generated by ./errors.pl, do not edit */

#include "grib_api_internal.h"

static const char *errors[] = {
"No error",		/* 0 GRIB_SUCCESS */
"End of resource reached",		/* -1 GRIB_END_OF_FILE */
"Internal error",		/* -2 GRIB_INTERNAL_ERROR */
"Passed buffer is too small",		/* -3 GRIB_BUFFER_TOO_SMALL */
"Function not yet implemented",		/* -4 GRIB_NOT_IMPLEMENTED */
"Missing 7777 at end of message",		/* -5 GRIB_7777_NOT_FOUND */
"Passed array is too small",		/* -6 GRIB_ARRAY_TOO_SMALL */
"File not found",		/* -7 GRIB_FILE_NOT_FOUND */
"Code not found in code table",		/* -8 GRIB_CODE_NOT_FOUND_IN_TABLE */
"Array size mismatch",		/* -9 GRIB_WRONG_ARRAY_SIZE */
"Key/value not found",		/* -10 GRIB_NOT_FOUND */
"Input output problem",		/* -11 GRIB_IO_PROBLEM */
"Message invalid",		/* -12 GRIB_INVALID_MESSAGE */
"Decoding invalid",		/* -13 GRIB_DECODING_ERROR */
"Encoding invalid",		/* -14 GRIB_ENCODING_ERROR */
"Code cannot unpack because of string too small",		/* -15 GRIB_NO_MORE_IN_SET */
"Problem with calculation of geographic attributes",		/* -16 GRIB_GEOCALCULUS_PROBLEM */
"Out of memory",		/* -17 GRIB_OUT_OF_MEMORY */
"Value is read only",		/* -18 GRIB_READ_ONLY */
"Invalid argument",		/* -19 GRIB_INVALID_ARGUMENT */
"Null handle",		/* -20 GRIB_NULL_HANDLE */
"Invalid section number",		/* -21 GRIB_INVALID_SECTION_NUMBER */
"Value cannot be missing",		/* -22 GRIB_VALUE_CANNOT_BE_MISSING */
"Wrong message length",		/* -23 GRIB_WRONG_LENGTH */
"Invalid key type",		/* -24 GRIB_INVALID_TYPE */
"Unable to set step",		/* -25 GRIB_WRONG_STEP */
"Wrong units for step (step must be integer)",		/* -26 GRIB_WRONG_STEP_UNIT */
"Invalid file id",		/* -27 GRIB_INVALID_FILE */
"Invalid grib id",		/* -28 GRIB_INVALID_GRIB */
"Invalid index id",		/* -29 GRIB_INVALID_INDEX */
"Invalid iterator id",		/* -30 GRIB_INVALID_ITERATOR */
"Invalid keys iterator id",		/* -31 GRIB_INVALID_KEYS_ITERATOR */
"Invalid nearest id",		/* -32 GRIB_INVALID_NEAREST */
"Invalid order by",		/* -33 GRIB_INVALID_ORDERBY */
"Missing a key from the fieldset",		/* -34 GRIB_MISSING_KEY */
"The point is out of the grid area",		/* -35 GRIB_OUT_OF_AREA */
"Concept no match",		/* -36 GRIB_CONCEPT_NO_MATCH */
"Definitions files not found",		/* -37 GRIB_NO_DEFINITIONS */
"Wrong type while packing",		/* -38 GRIB_WRONG_TYPE */
"End of resource",		/* -39 GRIB_END */
"Unable to code a field without values",		/* -40 GRIB_NO_VALUES */
"Grid description is wrong or inconsistent",		/* -41 GRIB_WRONG_GRID */
"End of index reached",		/* -42 GRIB_END_OF_INDEX */
"Null index",		/* -43 GRIB_NULL_INDEX */
"End of resource reached when reading message",		/* -44 GRIB_PREMATURE_END_OF_FILE */
"An internal array is too small",		/* -45 GRIB_INTERNAL_ARRAY_TOO_SMALL */
"Message is too large for the current architecture",		/* -46 GRIB_MESSAGE_TOO_LARGE */
"Constant field",		/* -47 GRIB_CONSTANT_FIELD */
"Switch unable to find a matching case",		/* -48 GRIB_SWITCH_NO_MATCH */
"Underflow",		/* -49 GRIB_UNDERFLOW */
"Message malformed",		/* -50 GRIB_MESSAGE_MALFORMED */
"Index is corrupted",		/* -51 GRIB_CORRUPTED_INDEX */
"Invalid number of bits per value",		/* -52 GRIB_INVALID_BPV */
"Edition of two messages is different",		/* -53 GRIB_DIFFERENT_EDITION */
"Value is different",		/* -54 GRIB_VALUE_DIFFERENT */
"Value mismatch",		/* 1 GRIB_VALUE_MISMATCH */
"double values are different",		/* 2 GRIB_DOUBLE_VALUE_MISMATCH */
"long values are different",		/* 3 GRIB_LONG_VALUE_MISMATCH */
"byte values are different",		/* 4 GRIB_BYTE_VALUE_MISMATCH */
"string values are different",		/* 5 GRIB_STRING_VALUE_MISMATCH */
"Offset mismatch",		/* 6 GRIB_OFFSET_MISMATCH */
"Count mismatch",		/* 7 GRIB_COUNT_MISMATCH */
"Name mismatch",		/* 8 GRIB_NAME_MISMATCH */
"Type mismatch",		/* 9 GRIB_TYPE_MISMATCH */
"Type and value mismatch",		/* 10 GRIB_TYPE_AND_VALUE_MISMATCH */
"Unable to compare accessors",		/* 11 GRIB_UNABLE_TO_COMPARE_ACCESSORS */
"Unable to reset iterator",		/* 12 GRIB_UNABLE_TO_RESET_ITERATOR */
"Assertion failure",		/* 13 GRIB_ASSERTION_FAILURE */
};

#define NUMBER(a) sizeof(a)/sizeof(a[0])

const char* grib_get_error_message(int code)
{
  code = -code;
  if(code <0 || code >= NUMBER(errors)) {
    static char mess[80];
    sprintf(mess,"Unknow error %d",code);
    return mess;
    }
  return errors[code];
}

void grib_check(const char* call,const char*  file,int line,int e,const char* msg)
{
	grib_context* c=grib_context_get_default();
    if(e) {
		if (file) {
			fprintf(stderr,"%s at line %d: %s failed: %s",
				file,line, call,grib_get_error_message(e));
			if (msg) fprintf(stderr," (%s)",msg);
			printf("\n");
		} else {
			grib_context_log(c,GRIB_LOG_ERROR,"%s",grib_get_error_message(e));
		}
        exit(e);
    }
}

void grib_fail(const char* expr,const char*  file,int line) {
   grib_context* c=grib_context_get_default();
   fprintf(stderr,"%s at line %d: assertion failure Assert(%s)",file,line,expr);
   if (c->no_abort) return;
   abort();
}

