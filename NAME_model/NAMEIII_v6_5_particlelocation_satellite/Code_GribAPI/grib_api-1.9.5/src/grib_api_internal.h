/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/**
    include private headers used for all internal functions of
    grib_api, not seen by the user of the API
  */

#ifndef grib_api_internal_H
#define grib_api_internal_H


#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#ifdef _LARGE_FILES
#undef _LARGE_FILE_API
#endif
#endif

#ifndef GRIB_INLINE
#define GRIB_INLINE
#endif

#if IS_BIG_ENDIAN

#if GRIB_MEM_ALIGN
#define FAST_BIG_ENDIAN 0
#else
#define FAST_BIG_ENDIAN 1
#endif

#endif

#if IEEE_BE
#define IEEE
#else
#if IEEE_LE
#define IEEE
#endif
#endif

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <inttypes.h>

#ifdef  HAVE_STRING_H
#include <string.h>
#else
#include <strings.h>
#endif


#if GRIB_LINUX_PTHREADS
#define PTHREAD_MUTEX_RECURSIVE PTHREAD_MUTEX_RECURSIVE_NP
extern int pthread_mutexattr_settype(pthread_mutexattr_t* attr,int type);
#endif

#if GRIB_PTHREADS

#include <pthread.h>
#define GRIB_PTHREAD_ONCE(a,b) pthread_once(a,b);
#define GRIB_MUTEX_LOCK(a) pthread_mutex_lock(a);
#define GRIB_MUTEX_UNLOCK(a) pthread_mutex_unlock(a);

#else

#define GRIB_PTHREAD_ONCE(a,b)
#define GRIB_MUTEX_LOCK(a)
#define GRIB_MUTEX_UNLOCK(a)

#endif

#ifndef HAVE_FSEEKO
#define fseeko fseek
#define ftello ftell
#endif

#define Assert(a) {if(!(a)) grib_fail(#a,__FILE__,__LINE__);}

#include "grib_api.h"

#define GRIB_UNKNOWN_VALUE   -9999.999
#define GRIB_KEY_UNDEF	"undef"

#define GRIB_HANDLE_BIG_ECMWF_GRIB1    1

#define MAX_FILE_HANDLES_WITH_MULTI 10
#define ACCESSORS_ARRAY_SIZE 2000
#define MAX_NUM_CONCEPTS 2000

#define GRIB_NAMESPACE      10

#define GRIB_MY_BUFFER      0
#define GRIB_USER_BUFFER    1

#define GRIB_REAL_MODE4    4
#define GRIB_REAL_MODE8    8

#define MAX_NUM_SECTIONS  9

#define GRIB_DISPOSABLE_MEMORY      0
#define GRIB_LONG_LASTING_MEMORY    1

#define GRIB_LOG_PERROR          (1<<10)

/* ACCESSOR COMPARE FLAGS */
#define GRIB_COMPARE_NAMES          (1<<0)
#define GRIB_COMPARE_TYPES          (1<<1)

typedef               struct     grib_expression         grib_expression;
typedef               struct     grib_arguments          grib_arguments;

typedef               struct     grib_action_file        grib_action_file;
typedef               struct     grib_action_file_list   grib_action_file_list;
typedef               struct     grib_block_of_accessors grib_block_of_accessors;
typedef               struct     grib_buffer             grib_buffer;
typedef               struct     grib_accessor_class     grib_accessor_class;
typedef               struct     grib_action             grib_action;
typedef               struct     grib_action_class       grib_action_class;
typedef               struct     grib_section            grib_section;
typedef               struct     grib_packer             grib_packer;
typedef               struct     grib_codetable          grib_codetable;

typedef               struct     grib_accessor           grib_accessor;
typedef               struct     grib_iterator_class     grib_iterator_class;
typedef               struct     grib_nearest_class      grib_nearest_class;
typedef               struct     grib_box_class          grib_box_class;
typedef               struct     grib_dumper             grib_dumper;
typedef               struct     grib_dumper_class       grib_dumper_class;
typedef               struct     grib_dependency         grib_dependency;
typedef               struct     string_feed             string_feed;

typedef void           (*nearest_init_class_proc)       (grib_nearest_class*);
typedef int            (*nearest_init_proc)             (grib_nearest* i,grib_handle*,grib_arguments*);

typedef int            (*nearest_find_proc)    (grib_nearest* nearest, grib_handle* h,
                                                  double inlat, double inlon,
                                                  unsigned long flags, double* outlats,
                                                  double* outlons,double *values,
                                                  double* distances, int* indexes,size_t *len);
typedef int            (*nearest_destroy_proc)   (grib_nearest* nearest);

typedef void           (*box_init_class_proc)    (grib_box_class*);
typedef int            (*box_destroy_proc)       (grib_box*);
typedef int            (*box_init_proc)          (grib_box* ,grib_handle*,grib_arguments*);
typedef grib_points*   (*box_get_points_proc)    (grib_box*, double, double, double,double, int*);

typedef void           (*iterator_init_class_proc)       (grib_iterator_class*);
typedef int            (*iterator_init_proc)             (grib_iterator* i,grib_handle*,grib_arguments*);

typedef int            (*iterator_next_proc)             (grib_iterator* i, double *lat, double *lon, double *val);
typedef int            (*iterator_previous_proc)         (grib_iterator* i, double *lat, double *lon, double *val);
typedef int            (*iterator_reset_proc)            (grib_iterator* i);
typedef int            (*iterator_destroy_proc)          (grib_iterator* i);
typedef long           (*iterator_has_next_proc)         (grib_iterator* i);

typedef  int  (*grib_pack_proc)                          (grib_handle* h,const double* in, size_t inlen, void*   out, size_t* outlen);
typedef  int  (*grib_unpack_proc)                        (grib_handle* h,const void*   in, size_t inlen, double* out, size_t* outlen);


typedef  void   (*accessor_destroy_proc)                  (grib_context* ,   grib_accessor* );

typedef  int     (*accessor_unpack_long_proc)             (grib_accessor*, long*,   size_t *len);
typedef  int     (*accessor_unpack_double_proc)           (grib_accessor*, double*, size_t *len);
typedef  int     (*accessor_unpack_double_element_proc)         (grib_accessor*, size_t, double*);
typedef  int     (*accessor_unpack_double_subarray_proc)         (grib_accessor*, double*,size_t , size_t);
typedef  int     (*accessor_unpack_string_proc)           (grib_accessor*, char*,   size_t *len);
typedef  int     (*accessor_unpack_bytes_proc)            (grib_accessor*, unsigned char*, size_t *len);
typedef  int     (*accessor_get_native_type_proc)         (grib_accessor*);
typedef  int     (*accessor_notify_change_proc)           (grib_accessor*,grib_accessor*);
typedef  void    (*accessor_update_size_proc)             (grib_accessor*,size_t);
typedef  size_t  (*accessor_preferred_size_proc)          (grib_accessor*,int);
typedef  void    (*accessor_resize_proc)                  (grib_accessor*,size_t);

typedef grib_accessor*  (*accessor_next_proc)            (grib_accessor*,int);
typedef  grib_section*   (*accessor_sub_section_proc)    (grib_accessor*);


typedef  int   (*accessor_pack_missing_proc)                (grib_accessor*);
typedef  int   (*accessor_pack_is_missing_proc)                (grib_accessor*);
typedef  int   (*accessor_pack_long_proc)                (grib_accessor*, const long*,   size_t *len);
typedef  int   (*accessor_pack_double_proc)              (grib_accessor*, const double*, size_t *len);
typedef  int   (*accessor_pack_string_proc)              (grib_accessor*, const char*,   size_t *len);
typedef  int   (*accessor_pack_bytes_proc)               (grib_accessor*, const unsigned char*, size_t *len);
typedef  int   (*accessor_pack_expression_proc)           (grib_accessor*, grib_expression*);
typedef  int   (*accessor_clear_proc)           (grib_accessor*);

typedef  void  (*accessor_init_class_proc)               (grib_accessor_class*);

typedef  int   (*accessor_compare_proc)                  (grib_accessor*, grib_accessor*);
typedef  long  (*accessor_value_proc)                    (grib_accessor*);
typedef  void  (*accessor_dump_proc)                     (grib_accessor*, grib_dumper*);
typedef  void  (*accessor_init_proc)                     (grib_accessor*, const long len, grib_arguments*);
typedef  void  (*accessor_post_init_proc)                (grib_accessor*);

typedef int (*accessor_nearest_proc) (grib_accessor*, double,double*);

typedef  long  (*grib_binop_long_proc)                        (long,long);
typedef  long  (*grib_unop_long_proc)                         (long);

typedef  double  (*grib_binop_double_proc)                        (double,double);
typedef  double  (*grib_unop_double_proc)                         (double);

typedef  int  (*grib_binop_string_proc)                        (char*,char*);

typedef struct second_order_packed second_order_packed;
typedef  void  grib_expression_visit_proc                (void* udata, grib_expression *e);



struct grib_key_value_list {
  const char* name;
  int         type;
  int         size;
  long*       long_value;
  double*     double_value;
  grib_key_value_list* namespace_value;
  char* string_value;
  int         has_value;
  int         error;
  grib_key_value_list* next;
} ;


struct second_order_packed {
  unsigned long nbits_per_widths;
  unsigned long nbits_per_group_size;
  size_t size_of_group_array;
  size_t packed_byte_count;
  unsigned long *array_of_group_size;
  unsigned long *array_of_group_width;
   long *array_of_group_refs;
} ;

/**
*  an grib_compression
*  Structure supporting the packing and unpacking procedures
*
*  @see  grib_action_create_data
*/
struct grib_packer {
    const char* name;
    grib_pack_proc   pack;     /** <  packing procedure                    */
    grib_unpack_proc unpack;   /** < unpacking procedure                    */
};


/* --------------- */

typedef struct grib_loader grib_loader;
typedef int (*grib_loader_init_accessor_proc)(grib_loader*,grib_accessor*,grib_arguments*);
typedef int (*grib_loader_lookup_long_proc) (grib_context*,grib_loader*,const char* name, long* value);

struct grib_loader {
  void                          *data;
  grib_loader_init_accessor_proc init_accessor;
  grib_loader_lookup_long_proc   lookup_long;
  int                            list_is_resized; /** will be true if we resize a list */
  int  changing_edition;
};

/**
*  an action
*  Structure supporting the creation of accessor, resulting of a statement during a definition file parsing
*
*  @see  grib_action_class
*/
struct grib_action
{
    char                     *name;   /**  name of the definition statement            */
    char                     *op;     /**  operator of the definition statement        */
    char                     *name_space;   /**  namspace of the definition statement  */
    grib_action              *next;   /**  next action in the list                     */
    grib_action_class        *cclass; /**  link to the structure containing a specific behavior */
    grib_context             *context;/**  Context                                     */
    unsigned long            flags;
    char                    *defaultkey; /** name of the key used as default if not found  */
	grib_arguments*         default_value; /** default expression as in .def file */
	char*					set;
    /* If you had something, don't forget to update grib_action_compile */
};

typedef struct grib_accessor_list grib_accessor_list;

struct grib_accessor_list {
	grib_accessor* accessor;
	grib_accessor_list* next;
};

/* compile */

typedef struct grib_compiler {
    int   cnt;
    int   max;
    FILE *out;
    const char *var;
} grib_compiler;



typedef  int  (*action_create_accessors_handle_proc)        (grib_section* p, grib_action* a, grib_loader* h);
typedef  int  (*action_notify_change_proc)                   (grib_action* a, grib_accessor* observer,grib_accessor * observed);

typedef  void  (*grib_dump_proc)                         (grib_action*, FILE*, int );
typedef  void  (*grib_xref_proc)                         (grib_action*, FILE*,const char*);
typedef  void  (*grib_compile_proc)                         (grib_action*, grib_compiler*);
typedef  void  (*action_init_class_proc)            (grib_action_class* a);
typedef  void  (*action_init_proc)                  (grib_action* a);
typedef  void  (*action_destroy_proc)                  (grib_context* context,  grib_action* a);
typedef  grib_action*  (*action_reparse_proc)              (grib_action* a,grib_accessor*,int*);
typedef  int  (*action_execute_proc)              (grib_action* a,grib_handle*);

/**
*  an action_class
*  Structure supporting the specific behavior of an action
*
*  @see  grib_action
*/
struct grib_action_class
{
    grib_action_class          **super; /** < link to a more general behavior                         */
    const char*                name;    /** < name of the behavior class                              */
    size_t                     size;    /** < size in bytes of the structure                          */

  int                        inited;
  action_init_class_proc     init_class;

    action_init_proc           init;
    action_destroy_proc        destroy;
                                        /** < destructor method to realease the memory    */

    grib_dump_proc             dump;    /** < dump method of the action                               */
    grib_xref_proc             xref;    /** < dump method of the action                               */
    action_create_accessors_handle_proc       create_accessor;
                                        /** < method to create the corresponding accessor from a handle*/
    action_notify_change_proc                      notify_change;
                                        /** < method to create the corresponding accessor from a handle*/

  action_reparse_proc              reparse;
  action_execute_proc              execute;

    grib_compile_proc             compile;    /** < compile method of the action                               */
};



/**
*  a buffer
*  Structure containing the datas of a Grib
*/
struct grib_buffer
{
    int              property;   /** < property parameter of buffer         */
    int              validity;   /** < validity parameter of buffer         */
    int              growable;   /** < buffer can be grown         */
    size_t           length;     /** < Buffer length                        */
    size_t           ulength;    /** < length used of the buffer            */
    unsigned char*   data;       /** < the data byte array                  */
};

/**
*  an Accessor
*  Structure supporting each single data unit and allowing its access
*  @see  grib_accessor_class
*/

#define MAX_ACCESSOR_NAMES 20

typedef struct grib_virtual_value grib_virtual_value;

struct grib_virtual_value {
  long     lval;       
  double   dval;  
  char*    cval; 
  int      missing;
  int      length;
  int      type;
};

struct grib_accessor
{
  const char             *name  ;     /** < name of the accessor                       */
  const char*             name_space;  /** < namespace to which the accessor belongs    */
  grib_action            *creator  ;  /** < action that created the accessor           */
  long                   length ;     /** < byte length of the accessor                */
  long                   offset ;     /** < offset of the data in the buffer           */
  grib_section           *parent;     /** < section to which the accessor is attached  */
  grib_accessor          *next  ;     /** < next accessor in list                      */
  grib_accessor          *previous;   /** < next accessor in list                      */
  grib_accessor_class    *cclass;     /** < behavior of the accessor                   */
  unsigned long           flags;      /** < Various flags                              */
  grib_section*           sub_section;

  const char*             all_names[MAX_ACCESSOR_NAMES]  ;   /** < name of the accessor  */
  const char*             all_name_spaces[MAX_ACCESSOR_NAMES]; /** < namespace to which the accessor belongs    */
  int                     dirty;

  grib_accessor          *same;      /** < accessors with the same name */
  long                   loop;      /** < used in lists */
  grib_virtual_value*    vvalue;    /** < virtual value used when transient flag on **/
  const char*            set;

};


#define GRIB_ACCESSOR_FLAG_READ_ONLY        (1<<1)
#define GRIB_ACCESSOR_FLAG_DUMP             (1<<2)
#define GRIB_ACCESSOR_FLAG_EDITION_SPECIFIC (1<<3)
#define GRIB_ACCESSOR_FLAG_CAN_BE_MISSING   (1<<4)
#define GRIB_ACCESSOR_FLAG_HIDDEN           (1<<5)
#define GRIB_ACCESSOR_FLAG_CONSTRAINT       (1<<6)
#define GRIB_ACCESSOR_FLAG_OVERRIDE         (1<<7)
#define GRIB_ACCESSOR_FLAG_NO_COPY          (1<<8)
#define GRIB_ACCESSOR_FLAG_COPY_OK          (1<<9)
#define GRIB_ACCESSOR_FLAG_FUNCTION         (1<<10)
#define GRIB_ACCESSOR_FLAG_DATA             (1<<11)
#define GRIB_ACCESSOR_FLAG_NO_FAIL          (1<<12)
#define GRIB_ACCESSOR_FLAG_TRANSIENT        (1<<13)
#define GRIB_ACCESSOR_FLAG_STRING_TYPE      (1<<14)
#define GRIB_ACCESSOR_FLAG_LONG_TYPE        (1<<15)
#define GRIB_ACCESSOR_FLAG_LOWERCASE        (1<<16)
/* when adding a flag, update grib_compile_flags*/

/**
*  a section accessor
*  Structure supporting hierarchical naming of the accessors
*  @see  grib_accessor
*/
struct grib_section
{
  grib_accessor                *owner;
  grib_handle                  *h;         /** < Handles of all accessors and buffer  */
  grib_accessor                *aclength;  /** < block of the length of the block     */
  grib_block_of_accessors      *block;     /** < block                                */
  grib_action                  *branch;    /** < branch that created the bolck        */
  size_t                       length;
  size_t                       padding;
};



struct grib_iterator_class{
   grib_iterator_class**           super;
   char* name;
   size_t                         size;

   int                            inited;
   iterator_init_class_proc       init_class;

   iterator_init_proc             init;
   iterator_destroy_proc        destroy;

   iterator_next_proc        next;
   iterator_previous_proc    previous;
   iterator_reset_proc       reset;
   iterator_has_next_proc    has_next;

};

struct grib_nearest_class{
   grib_nearest_class**           super;
   char* name;
   size_t                         size;

   int                            inited;
   nearest_init_class_proc       init_class;

   nearest_init_proc             init;
   nearest_destroy_proc          destroy;

   nearest_find_proc             find;

};

struct grib_box_class{
   grib_box_class**           super;
   char* name;
   size_t                     size;
   int                        inited;
   box_init_class_proc        init_class;
   box_init_proc              init;
   box_destroy_proc           destroy;
   box_get_points_proc        get_points;

};

/* --------------- */
/* --------------- */
typedef void (*search_all_callback_proc)(grib_accessor*,void* data);
/* --------------- */



typedef int  (*dumper_init_proc)        (grib_dumper*);
typedef void (*dumper_dump_proc)        (grib_dumper*,grib_accessor*,const char* comment);
typedef void (*dumper_dump_section_proc)(grib_dumper*,grib_accessor*,grib_block_of_accessors* block);
typedef void (*dumper_dump_values_proc) (grib_dumper*,grib_accessor*);
typedef int  (*dumper_destroy_proc)     (grib_dumper*);
typedef void (*dumper_header_proc)      (grib_dumper*,grib_handle*);
typedef void (*dumper_footer_proc)      (grib_dumper*,grib_handle*);
typedef void (*dumper_init_class_proc)  (grib_dumper_class*);

struct grib_dumper {
  FILE*            out;
   unsigned long     option_flags;
   void*             arg;
   int               depth;
   grib_handle       *handle;
   grib_dumper_class *cclass;
 };

struct grib_dumper_class {
   grib_dumper_class**    super;
   char*                    name;
   size_t                   size;
   int                      inited;
   dumper_init_class_proc   init_class;
   dumper_init_proc         init;
   dumper_destroy_proc      destroy;
   dumper_dump_proc         dump_long;
   dumper_dump_proc         dump_double;
   dumper_dump_proc         dump_string;
   dumper_dump_proc         dump_label;
   dumper_dump_proc         dump_bytes;
   dumper_dump_proc         dump_bits;
   dumper_dump_section_proc dump_section;
   dumper_dump_values_proc  dump_values;
   dumper_header_proc       header;
   dumper_footer_proc       footer;
};

struct grib_iterator{
   grib_arguments  *args;                   /**  args of iterator   */
   grib_handle* h;
   long        e;                           /**  current element    */
   size_t     nv;                           /**  number of values   */
   double*  data;                           /**  data values        */
   grib_iterator_class* cclass;
   unsigned long flags;
};

struct grib_nearest{
   grib_arguments              *args;      /**  args of iterator   */
   grib_handle*                h;
   grib_context*               context;
   double*                     values;
   size_t                      values_count;
   grib_nearest_class*         cclass;
   unsigned long               flags;

};

struct grib_box {
   grib_box_class*             cclass;
   grib_context*               context;
   grib_arguments              *args;      
   grib_handle*                h;
   unsigned long               flags;
   grib_points*                points;
};


struct grib_dependency {
   grib_dependency* next;
   grib_accessor*   observed;
   grib_accessor*   observer;
   int              run;
};


struct grib_block_of_accessors
{
    grib_accessor*    first;
    grib_accessor*    last ;
};


typedef struct grib_trie grib_trie;
typedef struct grib_itrie grib_itrie;


struct grib_darray {
  double* v;
  size_t size;
  size_t n;
  size_t incsize;
} ;

struct grib_iarray {
  long* v;
  size_t size;
  size_t n;
  size_t incsize;
} ;


#define MAX_SET_VALUES      10
#define MAX_ACCESSOR_CACHE  100



struct grib_handle
{
    grib_context*           context;       /** < context attached to this handle    */
    grib_buffer*            buffer ;       /** < buffer attached to the handle      */
    grib_section*           root;          /**  the root      section*/
    grib_section*           asserts;       /** the assertion section*/
    grib_section*           rules;         /** the rules     section*/
    grib_dependency*        dependencies;  /** List of dependencies */
    grib_handle*            main;           /** Used during reparsing */
    grib_handle*            kid;           /** Used during reparsing */
    grib_loader*            loader;        /** Used during reparsing */
    int                     values_stack;
    const grib_values*      values[MAX_SET_VALUES];       /** Used when setting multiple values at once */
    size_t                  values_count[MAX_SET_VALUES];  /** Used when setting multiple values at once */
    int                     dont_trigger;  /** Don't notify triggers */
    int                     partial;       /** Not a complete message (just headers) */
    int                     header_mode;   /** Header not jet complete */
    char* gts_header;
    size_t gts_header_len;
    int  use_trie;
    int  trie_invalid;
    grib_accessor* accessors[ACCESSORS_ARRAY_SIZE];
    char* section_offset[MAX_NUM_SECTIONS];
    char* section_length[MAX_NUM_SECTIONS];
    int sections_count;
    off_t offset;
};

struct grib_multi_handle {
  grib_context*           context;       /** < context attached to this handle    */
  grib_buffer*            buffer ;       /** < buffer attached to the handle      */
  size_t                  offset ;
  size_t                  length ;
};


struct grib_accessor_class
{
    grib_accessor_class             **super;
    const char*                     name;
    size_t                          size;

    int                             inited;
    accessor_init_class_proc        init_class;

    accessor_init_proc              init;
    accessor_post_init_proc         post_init;
    accessor_destroy_proc           destroy;

    accessor_dump_proc              dump;
    accessor_value_proc             next_offset;

    accessor_value_proc             value_count;

    accessor_value_proc             byte_count;
    accessor_value_proc             byte_offset;

    accessor_get_native_type_proc   get_native_type;

    accessor_sub_section_proc       sub_section;

    accessor_pack_missing_proc      pack_missing  ;
    accessor_pack_is_missing_proc   is_missing  ;

    accessor_pack_long_proc         pack_long  ;
    accessor_unpack_long_proc       unpack_long  ;

    accessor_pack_double_proc       pack_double;
    accessor_unpack_double_proc     unpack_double;

    accessor_pack_string_proc       pack_string;
    accessor_unpack_string_proc     unpack_string;

    accessor_pack_bytes_proc        pack_bytes;
    accessor_unpack_bytes_proc      unpack_bytes;

    accessor_pack_expression_proc   pack_expression;

    accessor_notify_change_proc     notify_change;
    accessor_update_size_proc       update_size;

    accessor_preferred_size_proc    preferred_size;
    accessor_resize_proc            resize;

    accessor_nearest_proc           nearest_smaller_value;
    accessor_next_proc              next;
    accessor_compare_proc           compare;
    accessor_unpack_double_element_proc     unpack_double_element;
    accessor_unpack_double_subarray_proc    unpack_double_subarray;
    accessor_clear_proc             clear;
};

typedef struct grib_multi_support grib_multi_support;

struct grib_multi_support {
    FILE*                           file;
    size_t                          offset;
    unsigned char*                  message;
    size_t                          message_length;
    unsigned char*                  sections[8];
    unsigned char*                  bitmap_section;
    size_t                          bitmap_section_length;
    size_t                          sections_length[9];
    int                             section_number;
    grib_multi_support*             next;
};

/* Concepts */
typedef struct grib_concept_condition grib_concept_condition;

struct grib_concept_condition {
  grib_concept_condition* next;
  char*               name;
  grib_expression*    expression;
};

typedef struct grib_concept_value_name grib_concept_value_name;
struct grib_concept_value_name {
  grib_concept_value_name*          next;
  char*               name;
} ;

typedef struct grib_concept_value grib_concept_value;

struct grib_concept_value {
  grib_concept_value*          next;
  char*                        name;
  grib_concept_condition*      conditions;
  grib_trie*		           index;
};

/* ----------*/

struct grib_context
{
    int                             inited;
    int                             debug;
    int                             write_on_fail;
    int                             no_abort;
	int 							io_buffer_size;
    char*                           grib_definition_files_path;
    char*                           grib_templates_path;
    char*                           grib_concept_path;

    grib_action_file_list*          grib_reader;
    void*                           user_data;
    int                             real_mode;

    grib_free_proc                  free_mem;
    grib_malloc_proc                alloc_mem;
    grib_realloc_proc               realloc_mem;

    grib_free_proc                  free_persistent_mem;
    grib_malloc_proc                alloc_persistent_mem;

    grib_free_proc                  free_buffer_mem;
    grib_malloc_proc                alloc_buffer_mem;
    grib_realloc_proc               realloc_buffer_mem;

    grib_data_read_proc             read;
    grib_data_write_proc            write;
    grib_data_tell_proc             tell;
    grib_data_seek_proc             seek;
    grib_data_eof_proc              eof;

    grib_log_proc                   output_log;
    grib_print_proc                 print;

    grib_codetable*                 codetable;
    char*                           outfilename;
    int                             multi_support_on;
    grib_multi_support*             multi_support;
    grib_string_list*               grib_definition_files_dir;
    int                             handle_file_count;
    int                             handle_total_count;
    off_t                           message_file_offset;
    int                             no_fail_on_wrong_length;
    int                             gts_header_on;
    int                             gribex_mode_on;
    grib_itrie*                     keys;
    int                             keys_count;
    grib_itrie*                     concepts_index;
    int                             concepts_count;
    grib_concept_value*             concepts[MAX_NUM_CONCEPTS];
    grib_trie*                      def_files;
    
    grib_string_list*                blacklist;
    int                             ieee_packing;
	grib_trie*                      classes;
#if GRIB_PTHREADS
    pthread_mutex_t                 mutex;
#endif

};

/* file_pool */
grib_string_list grib_file_not_found;

typedef struct grib_file grib_file;
typedef struct grib_file_pool grib_file_pool;

struct grib_file {
  grib_context* context;
  char* name;
  FILE* handle;
  char* mode;
  char* buffer;
  long refcount;
  grib_file* next;
  short id;
} ;

struct grib_file_pool {
  grib_context* context;
  grib_file* first;
  grib_file* current;
  size_t size;
  int number_of_opened_files;
  int max_opened_files;
};

/* fieldset */
typedef struct grib_field grib_field;
typedef struct grib_column grib_column;
typedef struct grib_fields grib_fields;
typedef struct grib_int_array grib_int_array;

struct grib_where {
  grib_context* context;
  char* string;
};

struct grib_column {
  grib_context* context;
  int refcount;
  char* name;
  int type;
  size_t size;
  size_t values_array_size;
  long* long_values;
  double* double_values;
  char** string_values;
  int* errors;
} ;

struct grib_order_by {
  char* key;
  int idkey;
  int mode;
  grib_order_by* next;
} ;

#ifdef NEWDB
struct grib_query {
  grib_context* context;
  char* where_string;
  grib_order_by* order_by;
};
#endif

struct grib_field {
  grib_file* file;
  off_t offset;
  long length;
  grib_field* next;
};

struct grib_int_array {
  grib_context* context;
  size_t size;
  int* el;
} ;

#ifndef NEWDB
struct grib_fieldset {
  grib_context* context;
  grib_int_array* filter;
  grib_int_array* order;
  size_t fields_array_size;
  size_t size;
  grib_column* columns;
  size_t columns_size;
  grib_where* where;
  grib_order_by* order_by;
  long current;
  grib_field** fields;
};
#endif

#ifdef NEWDB
/* grib db */
struct grib_db {
  grib_context* context;
  size_t size;
  size_t fields_array_size;
  grib_column* columns;
  size_t columns_size;
  grib_field** fields;
};

struct grib_fieldset {
  grib_context* context;
  grib_db* db;
  grib_int_array* filter;
  grib_int_array* order;
  size_t size;
  grib_query* query;
  long current;
};
#endif


/* index structures */

#define STRING_VALUE_LEN 100

typedef struct grib_field_tree grib_field_tree;

struct grib_field_tree {
  grib_field* field;
  char* value;
  grib_field_tree* next;
  grib_field_tree* next_level;
};

typedef struct grib_index_key grib_index_key;

struct grib_index_key {
  char* name;
  int type;
  char value[STRING_VALUE_LEN];
  grib_string_list* values;
  grib_string_list* current;
  int values_count;
  int count;
  grib_index_key* next;
};

typedef struct grib_field_list grib_field_list;
struct grib_field_list {
  grib_field* field;
  grib_field_list* next;
};


struct grib_index {
  grib_context* context;
  grib_index_key* keys;
  int rewind;
  int orderby;
  grib_index_key* orederby_keys;
  grib_field_tree* fields;
  grib_field_list* fieldset;
  grib_field_list* current;
  grib_file* files;
  int count;
};

/* header compute */
typedef struct grib_math grib_math;

struct grib_math{
  struct grib_math *left;
  struct grib_math *right;
  char        *name;
  int         arity;
};

typedef double (*mathproc)();
typedef int    (*funcproc)(grib_math*,mathproc);

typedef struct func {
  char    *name;
  funcproc addr;
  mathproc proc;
  int      arity;
  char     *info;
} func;

/* action file */
struct grib_action_file
{
    char*             filename   ;
    grib_action*      root       ;
    grib_action_file* next       ;
};

struct grib_action_file_list
{
    grib_action_file * first;
    grib_action_file * last ;
};

#include "grib_expression.h"

/* ----------*/
/* md5 */
typedef unsigned long cvs_uint32;

struct cvs_MD5Context {
  cvs_uint32 buf[4];
  cvs_uint32 bits[2];
  unsigned char in[64];
};
/* --- */

typedef struct grib_rule_entry grib_rule_entry;

struct grib_rule_entry {
  grib_rule_entry *next;
  char            *name;
  grib_expression *value;
};

typedef struct grib_rule grib_rule;

struct grib_rule {
  grib_rule        *next;
  grib_expression  *condition;
  grib_rule_entry  *entries;
};

typedef struct grib_case grib_case;

struct grib_case {
   grib_arguments* values;
   grib_action* action;
   grib_case* next;
};

/* ----------*/

typedef struct code_table_entry {
  char* abbreviation;
  char* title;
  char* units;
} code_table_entry;

struct grib_codetable {
  char*            filename[2];
  char*            recomposed_name[2];
  grib_codetable*  next;
  size_t           size;
  code_table_entry entries[1];
};

#if GRIB_TIMER
typedef struct grib_timer {

    struct timeval start_;
    double timer_;
    int   active_;
    char   *name_;
    int    count_;
    long total_;

    int elapsed_;
    double  cpu_;
    double  total_cpu_;

    char   *statname_;
    grib_context* context;

    struct grib_timer *next_;
} grib_timer;
#else
typedef struct grib_timer {
    char nothing;
} grib_timer;
#endif

typedef struct j2k_encode_helper {

  size_t           buffer_size;

  long             width;
  long             height;
  long             bits_per_value;

  float            compression;

  long             no_values;
  const double    *values;
  double           reference_value;
  double           divisor;
  double           decimal;

  long            jpeg_length;
  unsigned char*  jpeg_buffer;

} j2k_encode_helper;


#include "grib_api_prototypes.h"


#ifdef __cplusplus
}
#endif
#endif
/* This part is automatically generated by ./errors.pl, do not edit */
#ifndef grib_errors_internal_H
#define grib_errors_internal_H
/** Value mismatch */
#define GRIB_VALUE_MISMATCH		1
/** double values are different */
#define GRIB_DOUBLE_VALUE_MISMATCH		2
/** long values are different */
#define GRIB_LONG_VALUE_MISMATCH		3
/** byte values are different */
#define GRIB_BYTE_VALUE_MISMATCH		4
/** string values are different */
#define GRIB_STRING_VALUE_MISMATCH		5
/** Offset mismatch */
#define GRIB_OFFSET_MISMATCH		6
/** Count mismatch */
#define GRIB_COUNT_MISMATCH		7
/** Name mismatch */
#define GRIB_NAME_MISMATCH		8
/** Type mismatch */
#define GRIB_TYPE_MISMATCH		9
/** Type and value mismatch */
#define GRIB_TYPE_AND_VALUE_MISMATCH		10
/** Unable to compare accessors */
#define GRIB_UNABLE_TO_COMPARE_ACCESSORS		11
/** Unable to reset iterator */
#define GRIB_UNABLE_TO_RESET_ITERATOR		12
/** Assertion failure */
#define GRIB_ASSERTION_FAILURE		13
#endif
