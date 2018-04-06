/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

#include <stddef.h>

void f_sizeof(void *x,void *y, int *size) {
  *size=((char*)y)-((char*)x);
}
void f_sizeof_(void *x,void *y, int *size) {
  *size=((char*)y)-((char*)x);
}
void f_sizeof__(void *x,void *y, int *size) {
  *size=((char*)y)-((char*)x);
}

void check_double(double *x,double *y,char* ret) {
  *ret = ((char*)y)-((char*)x) == sizeof(*x) ?  't' : 'f';
}
void check_double_(double *x,double *y,char* ret) {check_double(x,y,ret);}
void check_double__(double *x,double *y,char* ret) {check_double(x,y,ret);}

void check_float(float *x,float *y,char* ret) {
  *ret = ((char*)y)-((char*)x) == sizeof(*x) ?  't' : 'f';
}
void check_float_(float *x,float *y,char* ret) { check_float(x,y,ret); }
void check_float__(float *x,float *y,char* ret) { check_float(x,y,ret); }

void check_int(int *x,int *y,char* ret) {
  *ret = ((char*)y)-((char*)x) == sizeof(*x) ?  't' : 'f';
}
void check_int_(int *x,int *y,char* ret) { check_int(x,y,ret); }
void check_int__(int *x,int *y,char* ret) { check_int(x,y,ret); }

void check_long(long *x,long *y,char* ret) {
  *ret = ((char*)y)-((char*)x) == sizeof(*x) ?  't' : 'f';
}
void check_long_(long *x,long *y,char* ret) {check_long(x,y,ret);}
void check_long__(long *x,long *y,char* ret) {check_long(x,y,ret);}

void check_size_t(size_t *x,size_t *y,char* ret) {
  *ret = ((char*)y)-((char*)x) == sizeof(*x) ?  't' : 'f';
}
void check_size_t_(size_t *x,size_t *y,char* ret) {check_size_t(x,y,ret);}
void check_size_t__(size_t *x,size_t *y,char* ret) {check_size_t(x,y,ret);}

