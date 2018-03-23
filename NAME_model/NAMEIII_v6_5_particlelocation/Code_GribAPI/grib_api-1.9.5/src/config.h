/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to 1 if the `closedir' function returns void instead of `int'. */
/* #undef CLOSEDIR_VOID */

/* Grib Api version */
#define GRIB_API_MAIN_VERSION 1.9.5

/* Grib Api Major release */
#define GRIB_API_MAJOR_VERSION 1

/* Grib Api Minor release */
#define GRIB_API_MINOR_VERSION 9

/* Grib Api Revision release */
#define GRIB_API_REVISION_VERSION 5

/* Directory where definition files are */
#define GRIB_DEFINITION_PATH "/home/h03/apdg/NameIII/NameIIILibrary/Code/Version6_5/SharedLibraries_Linux/share/definitions"

/* inline if available */
#define GRIB_INLINE inline

/* 1->pthreads enabled 0->pthreads disabled */
#define GRIB_LINUX_PTHREADS 0

/* memory alignment required */
#define GRIB_MEM_ALIGN 0

/* 1->pthreads enabled 0->pthreads disabled */
#define GRIB_PTHREADS 0

/* Directory where samples are */
#define GRIB_SAMPLES_PATH "/home/h03/apdg/NameIII/NameIIILibrary/Code/Version6_5/SharedLibraries_Linux/share/samples"

/* Directory where templates are */
#define GRIB_TEMPLATES_PATH "/home/h03/apdg/NameIII/NameIIILibrary/Code/Version6_5/SharedLibraries_Linux/share/samples"

/* 1->Timer on 0->Timer off */
#define GRIB_TIMER 0

/* Define to 1 if you have the <assert.h> header file. */
#define HAVE_ASSERT_H 1

/* Define to 1 if you have the `bzero' function. */
#define HAVE_BZERO 1

/* Define to 1 if you have the <ctype.h> header file. */
#define HAVE_CTYPE_H 1

/* Define to 1 if you have the <dirent.h> header file, and it defines `DIR'.
   */
#define HAVE_DIRENT_H 1

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
/* #undef HAVE_DOPRNT */

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if fseeko (and presumably ftello) exists and is declared. */
#define HAVE_FSEEKO 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <inttypes.h> header file. */
/* #undef HAVE_INTTYPES_H */

/* JPEG enabled */
/* #undef HAVE_JPEG */

/* Define if you have EMOS library */
/* #undef HAVE_LIBEMOS */

/* Define if you have JPEG version 2 "Jasper" library */
/* #undef HAVE_LIBJASPER */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define if you have JPEG version 2 "Openjpeg" library */
/* #undef HAVE_LIBOPENJPEG */

/* Define to 1 if you have the png library (-lpng) */
/* #undef HAVE_LIBPNG */

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the <memory.h> header file. */
/* #undef HAVE_MEMORY_H */

/* Define to 1 if you have the <ndir.h> header file, and it defines `DIR'. */
/* #undef HAVE_NDIR_H */

/* Define to 1 if you have the <stdarg.h> header file. */
#define HAVE_STDARG_H 1

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
/* #undef HAVE_STDINT_H */

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
/* #undef HAVE_STRINGS_H */

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/dir.h> header file, and it defines `DIR'.
   */
/* #undef HAVE_SYS_DIR_H */

/* Define to 1 if you have the <sys/ndir.h> header file, and it defines `DIR'.
   */
/* #undef HAVE_SYS_NDIR_H */

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
/* #undef HAVE_SYS_STAT_H */

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
/* #undef HAVE_SYS_TYPES_H */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vprintf' function. */
#define HAVE_VPRINTF 1

/* 1-> ieee big endian float/double 0->no ieee big endian float/double */
#define IEEE_BE 0

/* 1-> ieee little endian float/double 0->no ieee little endian float/double
   */
#define IEEE_LE 1

/* 1-> big endian 0->little endian */
#define IS_BIG_ENDIAN 0

/* memory management */
#define MANAGE_MEM 0

/* 1->OpenMP packing 0->single thread packing */
#define OMP_PACKING 0

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "Software.Services@ecmwf.int"

/* Define to the full name of this package. */
#define PACKAGE_NAME "src/grib_api.h"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "src/grib_api.h  "

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "src-grib_api-h"

/* Define to the version of this package. */
#define PACKAGE_VERSION " "

/* posix_memalign present */
#define POSIX_MEMALIGN 1

/* Define as the return type of signal handlers (`int' or `void'). */
#define RETSIGTYPE void

/* Define to 1 if you have the ANSI C header files. */
/* #undef STDC_HEADERS */

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* vectorised code */
#define VECTOR 0

/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
#define YYTEXT_POINTER 1

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define to 1 to make fseeko visible on some hosts (e.g. glibc 2.2). */
/* #undef _LARGEFILE_SOURCE */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Needs to be undefined on some AIX */
/* #undef _LARGE_FILE_API */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */
