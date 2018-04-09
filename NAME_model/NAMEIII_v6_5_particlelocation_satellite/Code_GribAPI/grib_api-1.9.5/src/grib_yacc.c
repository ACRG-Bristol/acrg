/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with grib_yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum grib_yytokentype {
     LOWERCASE = 258,
     IF = 259,
     IF_TRANSIENT = 260,
     ELSE = 261,
     END = 262,
     UNSIGNED = 263,
     TEMPLATE = 264,
     TEMPLATE_NOFAIL = 265,
     TRIGGER = 266,
     ASCII = 267,
     KSEC1EXPVER = 268,
     LABEL = 269,
     LIST = 270,
     WHILE = 271,
     IBMFLOAT = 272,
     SIGNED = 273,
     BYTE = 274,
     CODETABLE = 275,
     COMPLEX_CODETABLE = 276,
     LOOKUP = 277,
     ALIAS = 278,
     UNALIAS = 279,
     META = 280,
     POS = 281,
     INTCONST = 282,
     TRANS = 283,
     FLAGBIT = 284,
     CONCEPT = 285,
     CONCEPT_NOFAIL = 286,
     NIL = 287,
     DUMMY = 288,
     MODIFY = 289,
     READ_ONLY = 290,
     STRING_TYPE = 291,
     LONG_TYPE = 292,
     NO_COPY = 293,
     DUMP = 294,
     NO_FAIL = 295,
     EDITION_SPECIFIC = 296,
     OVERRIDE = 297,
     HIDDEN = 298,
     CAN_BE_MISSING = 299,
     MISSING = 300,
     CONSTRAINT = 301,
     COPY_OK = 302,
     WHEN = 303,
     SET = 304,
     SET_NOFAIL = 305,
     WRITE = 306,
     APPEND = 307,
     PRINT = 308,
     EXPORT = 309,
     REMOVE = 310,
     SKIP = 311,
     PAD = 312,
     SECTION_PADDING = 313,
     PADTO = 314,
     PADTOEVEN = 315,
     PADTOMULTIPLE = 316,
     G1_HALF_BYTE = 317,
     G1_MESSAGE_LENGTH = 318,
     G1_SECTION4_LENGTH = 319,
     SECTION_LENGTH = 320,
     FLAG = 321,
     ITERATOR = 322,
     NEAREST = 323,
     BOX = 324,
     KSEC = 325,
     ASSERT = 326,
     CASE = 327,
     SWITCH = 328,
     DEFAULT = 329,
     EQ = 330,
     NE = 331,
     GE = 332,
     LE = 333,
     BIT = 334,
     BITOFF = 335,
     AND = 336,
     OR = 337,
     NOT = 338,
     IS = 339,
     IDENT = 340,
     STRING = 341,
     INTEGER = 342,
     FLOAT = 343
   };
#endif
/* Tokens.  */
#define LOWERCASE 258
#define IF 259
#define IF_TRANSIENT 260
#define ELSE 261
#define END 262
#define UNSIGNED 263
#define TEMPLATE 264
#define TEMPLATE_NOFAIL 265
#define TRIGGER 266
#define ASCII 267
#define KSEC1EXPVER 268
#define LABEL 269
#define LIST 270
#define WHILE 271
#define IBMFLOAT 272
#define SIGNED 273
#define BYTE 274
#define CODETABLE 275
#define COMPLEX_CODETABLE 276
#define LOOKUP 277
#define ALIAS 278
#define UNALIAS 279
#define META 280
#define POS 281
#define INTCONST 282
#define TRANS 283
#define FLAGBIT 284
#define CONCEPT 285
#define CONCEPT_NOFAIL 286
#define NIL 287
#define DUMMY 288
#define MODIFY 289
#define READ_ONLY 290
#define STRING_TYPE 291
#define LONG_TYPE 292
#define NO_COPY 293
#define DUMP 294
#define NO_FAIL 295
#define EDITION_SPECIFIC 296
#define OVERRIDE 297
#define HIDDEN 298
#define CAN_BE_MISSING 299
#define MISSING 300
#define CONSTRAINT 301
#define COPY_OK 302
#define WHEN 303
#define SET 304
#define SET_NOFAIL 305
#define WRITE 306
#define APPEND 307
#define PRINT 308
#define EXPORT 309
#define REMOVE 310
#define SKIP 311
#define PAD 312
#define SECTION_PADDING 313
#define PADTO 314
#define PADTOEVEN 315
#define PADTOMULTIPLE 316
#define G1_HALF_BYTE 317
#define G1_MESSAGE_LENGTH 318
#define G1_SECTION4_LENGTH 319
#define SECTION_LENGTH 320
#define FLAG 321
#define ITERATOR 322
#define NEAREST 323
#define BOX 324
#define KSEC 325
#define ASSERT 326
#define CASE 327
#define SWITCH 328
#define DEFAULT 329
#define EQ 330
#define NE 331
#define GE 332
#define LE 333
#define BIT 334
#define BITOFF 335
#define AND 336
#define OR 337
#define NOT 338
#define IS 339
#define IDENT 340
#define STRING 341
#define INTEGER 342
#define FLOAT 343




/* Copy the first part of user declarations.  */
#line 11 "griby.y"


#include "grib_api_internal.h"
/* #include "grib_parser.h" */

extern int grib_yylex();
extern int grib_yyerror(const char*);

extern   grib_action*           grib_parser_all_actions;
extern   grib_concept_value*    grib_parser_concept;
extern   grib_context*          grib_parser_context;
extern   grib_rule*             grib_parser_rules;

static grib_concept_value* reverse(grib_concept_value* r);
static grib_concept_value *reverse_concept(grib_concept_value *r,grib_concept_value *s);

/* typedef int (*testp_proc)(long,long); */
/* typedef long (*grib_op_proc)(long,long);   */




/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 33 "griby.y"
{
    char                    *str;
    long                    lval;
    double                  dval;
    grib_darray             *dvalue;
    grib_iarray             *ivalue;
    grib_action             *act;
    grib_arguments          *explist;
    grib_expression         *exp;
    grib_concept_condition  *concept_condition;
    grib_concept_value      *concept_value;
	grib_case               *case_value;
  grib_rule               *rules;
  grib_rule_entry         *rule_entry;
}
/* Line 187 of yacc.c.  */
#line 310 "y.tab.c"
	YYSTYPE;
# define grib_yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 323 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 grib_yytype_uint8;
#else
typedef unsigned char grib_yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 grib_yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char grib_yytype_int8;
#else
typedef short int grib_yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 grib_yytype_uint16;
#else
typedef unsigned short int grib_yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 grib_yytype_int16;
#else
typedef short int grib_yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined grib_yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined grib_yyoverflow || YYERROR_VERBOSE */


#if (! defined grib_yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union grib_yyalloc
{
  grib_yytype_int16 grib_yyss;
  YYSTYPE grib_yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union grib_yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (grib_yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T grib_yyi;				\
	  for (grib_yyi = 0; grib_yyi < (Count); grib_yyi++)	\
	    (To)[grib_yyi] = (From)[grib_yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T grib_yynewbytes;						\
	YYCOPY (&grib_yyptr->Stack, Stack, grib_yysize);				\
	Stack = &grib_yyptr->Stack;						\
	grib_yynewbytes = grib_yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	grib_yyptr += grib_yynewbytes / sizeof (*grib_yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  156
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1247

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  108
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  45
/* YYNRULES -- Number of rules.  */
#define YYNRULES  206
/* YYNRULES -- Number of states.  */
#define YYNSTATES  644

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   343

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? grib_yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const grib_yytype_uint8 grib_yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,   104,     2,     2,
      93,    94,   102,   105,    89,   100,    95,   103,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    99,    90,
     107,    96,   106,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    91,     2,    92,   101,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    97,     2,    98,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const grib_yytype_uint16 grib_yyprhs[] =
{
       0,     0,     3,     5,     7,     9,    11,    13,    14,    16,
      20,    22,    26,    28,    31,    35,    38,    41,    43,    45,
      47,    49,    51,    53,    55,    57,    60,    62,    64,    66,
      70,    72,    80,    91,    99,   107,   115,   126,   134,   142,
     153,   162,   175,   186,   195,   206,   215,   225,   233,   236,
     239,   244,   251,   259,   262,   268,   274,   279,   286,   294,
     297,   303,   312,   321,   325,   331,   337,   343,   349,   353,
     357,   361,   367,   370,   378,   383,   392,   403,   409,   415,
     421,   427,   430,   435,   439,   444,   449,   456,   461,   464,
     466,   469,   471,   474,   480,   482,   490,   502,   510,   522,
     529,   537,   549,   554,   559,   562,   566,   568,   571,   573,
     576,   578,   582,   584,   586,   588,   590,   592,   594,   596,
     598,   600,   602,   604,   606,   608,   610,   619,   627,   635,
     642,   652,   665,   680,   691,   706,   719,   728,   740,   747,
     757,   770,   781,   796,   809,   818,   830,   832,   835,   837,
     840,   845,   856,   866,   874,   880,   886,   892,   898,   900,
     903,   908,   910,   912,   914,   916,   918,   920,   922,   926,
     929,   933,   938,   942,   944,   948,   952,   956,   960,   964,
     966,   970,   974,   976,   980,   984,   988,   992,   996,  1000,
    1004,  1007,  1009,  1013,  1015,  1019,  1021,  1023,  1025,  1027,
    1032,  1035,  1037,  1040,  1042,  1050,  1052
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const grib_yytype_int16 grib_yyrhs[] =
{
     109,     0,    -1,   110,    -1,   131,    -1,   112,    -1,   152,
      -1,     1,    -1,    -1,    88,    -1,   111,    89,    88,    -1,
      87,    -1,   111,    89,    87,    -1,   113,    -1,   113,   112,
      -1,   113,    90,   112,    -1,   113,    90,    -1,   118,    90,
      -1,   119,    -1,   127,    -1,   128,    -1,   129,    -1,   130,
      -1,   120,    -1,   134,    -1,    90,    -1,   114,    90,    -1,
     110,    -1,   116,    -1,   117,    -1,   117,    89,   116,    -1,
     146,    -1,     8,    91,    87,    92,    85,   123,   124,    -1,
       8,    91,    87,    92,    85,    91,   115,    92,   123,   124,
      -1,    12,    91,    87,    92,    85,   123,   124,    -1,    12,
      91,    87,    92,    86,   123,   124,    -1,    19,    91,    87,
      92,    85,   123,   124,    -1,    19,    91,    87,    92,    85,
      91,   115,    92,   123,   124,    -1,    13,    91,    87,    92,
      85,   123,   124,    -1,    18,    91,    87,    92,    85,   123,
     124,    -1,    18,    91,    87,    92,    85,    91,   115,    92,
     123,   124,    -1,    20,    91,    87,    92,    85,   117,   123,
     124,    -1,    20,    91,    87,    92,    85,   117,   123,    49,
      93,    85,    94,   124,    -1,    20,    91,    87,    92,    85,
      93,   115,    94,   123,   124,    -1,    21,    91,    87,    92,
      85,   117,   123,   124,    -1,    21,    91,    87,    92,    85,
      93,   115,    94,   123,   124,    -1,    66,    91,    87,    92,
      85,   117,   123,   124,    -1,    22,    91,    87,    92,    85,
      93,   115,    94,   124,    -1,    29,    85,    93,   115,    94,
     123,   124,    -1,    14,    85,    -1,    14,    86,    -1,    17,
      85,   123,   124,    -1,    17,    85,    95,    85,   123,   124,
      -1,    17,    85,    91,   117,    92,   123,   124,    -1,    26,
      85,    -1,    27,    85,    96,   117,   124,    -1,    28,    85,
      96,   117,   124,    -1,    88,    85,   123,   124,    -1,    88,
      85,    95,    85,   123,   124,    -1,    88,    85,    91,   117,
      92,   123,   124,    -1,    62,    85,    -1,    65,    91,    87,
      92,    85,    -1,    63,    91,    87,    92,    85,    93,   115,
      94,    -1,    64,    91,    87,    92,    85,    93,   115,    94,
      -1,    70,    85,   117,    -1,    57,    85,    93,   115,    94,
      -1,    59,    85,    93,   115,    94,    -1,    60,    85,    93,
     115,    94,    -1,    61,    85,    93,   115,    94,    -1,    58,
      85,   124,    -1,     9,    85,    86,    -1,    10,    85,    86,
      -1,    23,    85,    96,    85,   124,    -1,    24,    85,    -1,
      23,    85,    95,    85,    96,    85,   124,    -1,    24,    85,
      95,    85,    -1,    25,    85,    85,    93,   115,    94,   123,
     124,    -1,    25,    85,    95,    85,    85,    93,   115,    94,
     123,   124,    -1,    67,    85,    93,   115,    94,    -1,    68,
      85,    93,   115,    94,    -1,    69,    85,    93,   115,    94,
      -1,    54,    85,    93,   115,    94,    -1,    55,   115,    -1,
      71,    93,   146,    94,    -1,    34,    85,   124,    -1,    49,
      85,    96,    45,    -1,    49,    85,    96,   146,    -1,    49,
      85,    96,    97,   111,    98,    -1,    50,    85,    96,   146,
      -1,    51,    86,    -1,    51,    -1,    52,    86,    -1,    52,
      -1,    53,    86,    -1,    53,    93,    86,    94,    86,    -1,
      53,    -1,     4,    93,   146,    94,    97,   112,    98,    -1,
       4,    93,   146,    94,    97,   112,    98,     6,    97,   112,
      98,    -1,     5,    93,   146,    94,    97,   112,    98,    -1,
       5,    93,   146,    94,    97,   112,    98,     6,    97,   112,
      98,    -1,    48,    93,   146,    94,   121,   114,    -1,    48,
      93,   146,    94,    97,   122,    98,    -1,    48,    93,   146,
      94,    97,   122,    98,     6,    97,   122,    98,    -1,    49,
      85,    96,   146,    -1,    50,    85,    96,   146,    -1,   121,
     114,    -1,   122,   121,   114,    -1,   110,    -1,    96,   115,
      -1,   110,    -1,    99,   125,    -1,   126,    -1,   125,    89,
     126,    -1,    35,    -1,     3,    -1,    39,    -1,    38,    -1,
      40,    -1,    43,    -1,    41,    -1,    44,    -1,    46,    -1,
      42,    -1,    47,    -1,    28,    -1,    36,    -1,    37,    -1,
      85,    15,    93,   146,    94,    97,   112,    98,    -1,    16,
      93,   146,    94,    97,   112,    98,    -1,    11,    93,   115,
      94,    97,   112,    98,    -1,    30,    85,    97,   131,    98,
     124,    -1,    30,    85,    93,    85,    94,    97,   131,    98,
     124,    -1,    30,    85,    93,    85,    89,    86,    89,    85,
      89,    85,    94,   124,    -1,    30,    85,    93,    85,    89,
      86,    89,    85,    89,    85,    89,    85,    94,   124,    -1,
      30,    85,    93,    85,    89,    86,    89,    85,    94,   124,
      -1,    30,    85,    95,    85,    93,    85,    89,    86,    89,
      85,    89,    85,    94,   124,    -1,    30,    85,    95,    85,
      93,    85,    89,    86,    89,    85,    94,   124,    -1,    30,
      85,    95,    85,    97,   131,    98,   124,    -1,    30,    85,
      95,    85,    93,    85,    94,    97,   131,    98,   124,    -1,
      31,    85,    97,   131,    98,   124,    -1,    31,    85,    93,
      85,    94,    97,   131,    98,   124,    -1,    31,    85,    93,
      85,    89,    86,    89,    85,    89,    85,    94,   124,    -1,
      31,    85,    93,    85,    89,    86,    89,    85,    94,   124,
      -1,    31,    85,    95,    85,    93,    85,    89,    86,    89,
      85,    89,    85,    94,   124,    -1,    31,    85,    95,    85,
      93,    85,    89,    86,    89,    85,    94,   124,    -1,    31,
      85,    95,    85,    97,   131,    98,   124,    -1,    31,    85,
      95,    85,    93,    85,    94,    97,   131,    98,   124,    -1,
     135,    -1,   131,   135,    -1,   133,    -1,   132,   133,    -1,
      72,   116,    99,   112,    -1,    73,    93,   115,    94,    97,
     132,    74,    99,   112,    98,    -1,    73,    93,   115,    94,
      97,   132,    74,    99,    98,    -1,    73,    93,   115,    94,
      97,   132,    98,    -1,    86,    96,    97,   136,    98,    -1,
      85,    96,    97,   136,    98,    -1,    87,    96,    97,   136,
      98,    -1,    88,    96,    97,   136,    98,    -1,   137,    -1,
     137,   136,    -1,    85,    96,   146,    90,    -1,    85,    -1,
      86,    -1,   138,    -1,    87,    -1,    88,    -1,    32,    -1,
      33,    -1,    93,   146,    94,    -1,   100,   139,    -1,    85,
      93,    94,    -1,    85,    93,   115,    94,    -1,   139,   101,
     140,    -1,   139,    -1,   141,   102,   140,    -1,   141,   103,
     140,    -1,   141,   104,   140,    -1,   141,    79,   140,    -1,
     141,    80,   140,    -1,   140,    -1,   142,   105,   141,    -1,
     142,   100,   141,    -1,   141,    -1,   143,   106,   142,    -1,
     143,    75,   142,    -1,   143,   107,   142,    -1,   143,    77,
     142,    -1,   143,    78,   142,    -1,   143,    76,   142,    -1,
     138,    84,   138,    -1,    83,   143,    -1,   142,    -1,   144,
      81,   143,    -1,   143,    -1,   145,    82,   144,    -1,   144,
      -1,   145,    -1,   150,    -1,   151,    -1,    85,    96,   146,
      90,    -1,    56,    90,    -1,   148,    -1,   148,   149,    -1,
     148,    -1,     4,    93,   146,    94,    97,   149,    98,    -1,
     147,    -1,   147,   152,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const grib_yytype_uint16 grib_yyrline[] =
{
       0,   204,   204,   205,   206,   207,   209,   212,   215,   216,
     217,   218,   220,   221,   222,   223,   226,   227,   228,   229,
     230,   231,   232,   233,   237,   238,   241,   242,   245,   246,
     249,   253,   256,   259,   263,   266,   269,   272,   275,   278,
     281,   284,   288,   291,   294,   297,   300,   303,   306,   309,
     312,   315,   318,   321,   324,   327,   330,   333,   336,   339,
     342,   345,   348,   351,   354,   357,   360,   363,   366,   369,
     371,   374,   377,   380,   384,   388,   391,   394,   406,   418,
     430,   433,   436,   439,   442,   443,   444,   446,   449,   450,
     451,   452,   454,   455,   456,   460,   461,   462,   463,   467,
     468,   469,   472,   473,   476,   477,   481,   482,   485,   486,
     489,   490,   493,   494,   495,   496,   497,   498,   499,   500,
     501,   502,   503,   504,   505,   506,   509,   512,   515,   518,
     519,   520,   521,   522,   523,   524,   525,   526,   527,   528,
     529,   530,   531,   532,   533,   534,   538,   539,   542,   543,
     546,   550,   551,   552,   555,   557,   559,   561,   565,   566,
     569,   572,   573,   576,   577,   578,   580,   581,   582,   583,
     584,   585,   589,   590,   593,   594,   595,   596,   597,   598,
     601,   602,   603,   606,   608,   609,   610,   611,   612,   613,
     618,   619,   622,   623,   626,   627,   630,   636,   637,   640,
     641,   644,   645,   648,   652,   655,   656
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const grib_yytname[] =
{
  "$end", "error", "$undefined", "LOWERCASE", "IF", "IF_TRANSIENT",
  "ELSE", "END", "UNSIGNED", "TEMPLATE", "TEMPLATE_NOFAIL", "TRIGGER",
  "ASCII", "KSEC1EXPVER", "LABEL", "LIST", "WHILE", "IBMFLOAT", "SIGNED",
  "BYTE", "CODETABLE", "COMPLEX_CODETABLE", "LOOKUP", "ALIAS", "UNALIAS",
  "META", "POS", "INTCONST", "TRANS", "FLAGBIT", "CONCEPT",
  "CONCEPT_NOFAIL", "NIL", "DUMMY", "MODIFY", "READ_ONLY", "STRING_TYPE",
  "LONG_TYPE", "NO_COPY", "DUMP", "NO_FAIL", "EDITION_SPECIFIC",
  "OVERRIDE", "HIDDEN", "CAN_BE_MISSING", "MISSING", "CONSTRAINT",
  "COPY_OK", "WHEN", "SET", "SET_NOFAIL", "WRITE", "APPEND", "PRINT",
  "EXPORT", "REMOVE", "SKIP", "PAD", "SECTION_PADDING", "PADTO",
  "PADTOEVEN", "PADTOMULTIPLE", "G1_HALF_BYTE", "G1_MESSAGE_LENGTH",
  "G1_SECTION4_LENGTH", "SECTION_LENGTH", "FLAG", "ITERATOR", "NEAREST",
  "BOX", "KSEC", "ASSERT", "CASE", "SWITCH", "DEFAULT", "EQ", "NE", "GE",
  "LE", "BIT", "BITOFF", "AND", "OR", "NOT", "IS", "IDENT", "STRING",
  "INTEGER", "FLOAT", "','", "';'", "'['", "']'", "'('", "')'", "'.'",
  "'='", "'{'", "'}'", "':'", "'-'", "'^'", "'*'", "'/'", "'%'", "'+'",
  "'>'", "'<'", "$accept", "all", "empty", "dvalues", "instructions",
  "instruction", "semi", "argument_list", "arguments", "argument",
  "simple", "if_block", "when_block", "set", "set_list", "default",
  "flags", "flag_list", "flag", "list_block", "while_block",
  "trigger_block", "concept_block", "concept_list", "case_list",
  "case_value", "switch_block", "concept_value", "concept_conditions",
  "concept_condition", "string_or_ident", "atom", "power", "factor",
  "term", "condition", "conjonction", "disjonction", "expression", "rule",
  "rule_entry", "rule_entries", "fact", "conditional_rule", "rules", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const grib_yytype_uint16 grib_yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,    44,
      59,    91,    93,    40,    41,    46,    61,   123,   125,    58,
      45,    94,    42,    47,    37,    43,    62,    60
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const grib_yytype_uint8 grib_yyr1[] =
{
       0,   108,   109,   109,   109,   109,   109,   110,   111,   111,
     111,   111,   112,   112,   112,   112,   113,   113,   113,   113,
     113,   113,   113,   113,   114,   114,   115,   115,   116,   116,
     117,   118,   118,   118,   118,   118,   118,   118,   118,   118,
     118,   118,   118,   118,   118,   118,   118,   118,   118,   118,
     118,   118,   118,   118,   118,   118,   118,   118,   118,   118,
     118,   118,   118,   118,   118,   118,   118,   118,   118,   118,
     118,   118,   118,   118,   118,   118,   118,   118,   118,   118,
     118,   118,   118,   118,   118,   118,   118,   118,   118,   118,
     118,   118,   118,   118,   118,   119,   119,   119,   119,   120,
     120,   120,   121,   121,   122,   122,   123,   123,   124,   124,
     125,   125,   126,   126,   126,   126,   126,   126,   126,   126,
     126,   126,   126,   126,   126,   126,   127,   128,   129,   130,
     130,   130,   130,   130,   130,   130,   130,   130,   130,   130,
     130,   130,   130,   130,   130,   130,   131,   131,   132,   132,
     133,   134,   134,   134,   135,   135,   135,   135,   136,   136,
     137,   138,   138,   139,   139,   139,   139,   139,   139,   139,
     139,   139,   140,   140,   141,   141,   141,   141,   141,   141,
     142,   142,   142,   143,   143,   143,   143,   143,   143,   143,
     143,   143,   144,   144,   145,   145,   146,   147,   147,   148,
     148,   149,   149,   150,   151,   152,   152
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const grib_yytype_uint8 grib_yyr2[] =
{
       0,     2,     1,     1,     1,     1,     1,     0,     1,     3,
       1,     3,     1,     2,     3,     2,     2,     1,     1,     1,
       1,     1,     1,     1,     1,     2,     1,     1,     1,     3,
       1,     7,    10,     7,     7,     7,    10,     7,     7,    10,
       8,    12,    10,     8,    10,     8,     9,     7,     2,     2,
       4,     6,     7,     2,     5,     5,     4,     6,     7,     2,
       5,     8,     8,     3,     5,     5,     5,     5,     3,     3,
       3,     5,     2,     7,     4,     8,    10,     5,     5,     5,
       5,     2,     4,     3,     4,     4,     6,     4,     2,     1,
       2,     1,     2,     5,     1,     7,    11,     7,    11,     6,
       7,    11,     4,     4,     2,     3,     1,     2,     1,     2,
       1,     3,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     8,     7,     7,     6,
       9,    12,    14,    10,    14,    12,     8,    11,     6,     9,
      12,    10,    14,    12,     8,    11,     1,     2,     1,     2,
       4,    10,     9,     7,     5,     5,     5,     5,     1,     2,
       4,     1,     1,     1,     1,     1,     1,     1,     3,     2,
       3,     4,     3,     1,     3,     3,     3,     3,     3,     1,
       3,     3,     1,     3,     3,     3,     3,     3,     3,     3,
       2,     1,     3,     1,     3,     1,     1,     1,     1,     4,
       2,     1,     2,     1,     7,     1,     2
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const grib_yytype_uint8 grib_yydefact[] =
{
       0,     6,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    89,    91,    94,     0,     7,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     2,     4,
      12,     0,    17,    22,    18,    19,    20,    21,     3,    23,
     146,   205,   203,   197,   198,     5,     0,     0,     0,     0,
       0,     7,     0,     0,    48,    49,     0,     7,     0,     0,
       0,     0,     0,     0,    72,     0,    53,     0,     0,     0,
       0,     0,     7,     0,     0,     0,    88,    90,    92,     0,
       0,   166,   167,     0,   161,   162,   164,   165,     0,     0,
      26,    81,    27,    28,   163,   173,   179,   182,   191,   193,
     195,   196,    30,   200,     0,     7,     0,     0,     0,    59,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     7,
       0,     0,     0,     0,     7,     0,     1,     0,     0,     0,
      15,    13,    16,     0,     0,   147,     0,     0,   206,     0,
       0,     0,    69,    70,     0,     0,     0,     0,     0,     0,
       7,   106,     7,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     7,     0,     0,     0,     0,
       0,     0,     0,   108,    83,     0,     0,     0,     0,     7,
     190,     0,     0,   163,   169,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     7,    68,     7,     7,     7,     0,     0,
       0,     0,     7,     7,     7,    63,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     7,     0,     0,    14,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       7,   107,    50,     0,     0,     0,     0,     0,     0,     7,
      74,     7,     0,     7,     7,     0,     0,     0,     0,     0,
       0,     0,   113,   123,   112,   124,   125,   115,   114,   116,
     118,   121,   117,   119,   120,   122,   109,   110,     0,    84,
       0,    85,    87,     0,     0,   170,     0,   168,    29,   161,
     189,   172,   177,   178,   174,   175,   176,   181,   180,   184,
     188,   186,   187,   183,   185,   192,   194,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    82,     0,
       0,     0,     0,   158,   199,     0,     0,     0,     7,    56,
       0,     0,     0,     0,     0,     7,     0,     7,     7,     7,
       0,     7,     7,     7,     7,     0,     0,     0,     0,    71,
       0,     0,    54,    55,     7,     0,     0,     0,     0,     7,
       0,     0,     0,     0,     7,     0,     0,     0,     0,     0,
      10,     8,     0,    93,    80,   171,    64,    65,    66,    67,
       0,     0,    60,     0,    77,    78,    79,     0,     0,     0,
     155,   159,   154,   156,     7,     7,   157,     0,     0,     0,
       0,   201,     0,     0,     7,     7,     0,     7,     7,     7,
       0,     7,    51,     7,     7,     7,     7,     7,     7,     7,
       7,     7,     7,     7,     7,     7,     0,     0,     0,     0,
     129,     0,     0,     0,     0,   138,   111,     0,     0,     0,
       0,    24,    99,     0,    86,     7,     7,     7,     0,     0,
     148,     0,     0,     7,    57,     0,     0,    95,   202,   204,
      97,     0,    31,   128,    33,    34,    37,   127,    52,     0,
      38,     0,    35,     0,    30,     7,     0,     7,     0,    73,
       7,     0,    47,     0,     0,     0,     0,     7,     0,     0,
       0,     0,     7,     0,     0,   104,   100,     0,    25,    11,
       9,     0,     0,     7,     0,     0,   153,   149,     0,   160,
      58,     0,     0,     7,     7,     7,     7,     0,    40,     7,
      43,     7,    75,     7,     0,     7,     0,     0,   136,     0,
       7,     0,     0,   144,   102,   103,     0,   105,    61,    62,
      45,     0,     0,   126,     0,     0,     7,     7,     7,     7,
       0,     7,    46,     7,     0,     7,   130,     0,     0,     0,
       7,   139,     0,     0,     0,   150,   152,     0,     0,     0,
      32,    39,    36,    42,     0,    44,    76,     0,   133,     0,
       7,     0,   141,     0,     7,     0,   151,    96,    98,     7,
       0,     7,     0,     7,   137,     7,     0,     7,   145,   101,
      41,     0,   131,     0,   135,   140,     0,   143,     7,     7,
       7,   132,   134,   142
};

/* YYDEFGOTO[NTERM-NUM].  */
static const grib_yytype_int16 grib_yydefgoto[] =
{
      -1,    57,   203,   402,   430,    60,   472,   121,   122,   123,
      61,    62,    63,   469,   470,   182,   204,   306,   307,    64,
      65,    66,    67,    68,   479,   480,    69,    70,   352,   353,
     124,   125,   126,   127,   128,   129,   130,   131,   132,    71,
     431,   432,    73,    74,    75
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -454
static const grib_yytype_int16 grib_yypact[] =
{
     905,  -454,   -52,   -26,    -5,   -20,    17,    25,    13,    50,
      15,    65,    66,    82,    88,    97,   100,   108,    86,   146,
     150,   153,   156,   162,   169,   186,   191,   201,   203,   217,
     235,   240,   242,   -67,   267,    89,   253,   271,   274,   276,
     277,   278,   282,   269,   286,   294,   295,   283,   284,   288,
     302,   296,   297,    -7,   292,   298,   -60,   392,  -454,  -454,
     990,   305,  -454,  -454,  -454,  -454,  -454,  -454,   263,  -454,
    -454,     2,  -454,  -454,  -454,  -454,    89,    89,   309,   311,
     313,    89,   314,   315,  -454,  -454,    89,   -58,   316,   317,
     318,   319,   320,   102,   330,   -61,  -454,   304,   332,   321,
     210,   226,   310,    89,   333,   334,  -454,  -454,  -454,   325,
     339,  -454,  -454,    89,   340,  -454,  -454,  -454,    89,   213,
    -454,  -454,  -454,   324,   350,   335,  -454,   101,   -46,   -31,
     354,   360,  -454,  -454,   351,   310,   352,   356,   358,  -454,
     372,   375,   377,   378,   373,   374,   376,    89,    89,    89,
     379,    20,   381,   382,   117,   383,  -454,   384,   455,   386,
    1159,  -454,  -454,   389,   393,  -454,   395,   394,  -454,   397,
     398,   401,  -454,  -454,   400,   403,   404,   406,    89,   388,
      89,  -454,   310,   405,   409,   410,   411,   412,   390,   391,
     413,   414,   420,    89,    89,    89,   421,   423,   263,   425,
     429,   263,   380,  -454,  -454,   422,   -15,    89,   427,    89,
     -31,    59,   430,  -454,  -454,    89,    30,   213,   213,   213,
     213,   213,   213,   213,   213,   213,   213,   213,   213,   213,
     213,    89,    89,    89,  -454,    89,    89,    89,   431,   433,
     438,   439,    89,    89,    89,  -454,   434,   440,    89,   447,
     443,   447,   447,    89,   450,   310,   447,    89,  -454,   441,
      89,    89,   442,   444,   451,   445,   172,   452,   446,   448,
     449,  -454,  -454,   459,   461,   462,   463,   464,   458,   310,
    -454,    89,   470,   310,   310,   465,   -33,   124,    11,   -32,
     137,    80,  -454,  -454,  -454,  -454,  -454,  -454,  -454,  -454,
    -454,  -454,  -454,  -454,  -454,  -454,   385,  -454,   -28,  -454,
     207,  -454,  -454,   472,   467,  -454,   468,  -454,  -454,  -454,
    -454,  -454,  -454,  -454,  -454,  -454,  -454,   101,   101,   -46,
     -46,   -46,   -46,   -46,   -46,   -31,   354,   469,   471,   473,
     474,   479,   481,   484,   485,   482,   483,   486,  -454,   476,
     488,   487,   480,   447,  -454,   489,   490,   494,   449,  -454,
     491,   496,   498,  1077,  1159,   -12,  1159,   449,   449,   449,
    1159,   449,   310,    52,    58,   107,   204,   426,   508,  -454,
     500,   502,  -454,  -454,   449,   511,   501,   514,   263,   310,
     515,   503,   517,   263,   310,   380,   518,   519,   275,   516,
    -454,  -454,   -34,  -454,  -454,  -454,  -454,  -454,  -454,  -454,
     512,   520,  -454,    89,  -454,  -454,  -454,   536,   513,    89,
    -454,  -454,  -454,  -454,   449,   310,  -454,   521,   522,    -6,
     523,   -42,   524,   525,    89,   310,   526,   310,   310,   310,
     528,   310,  -454,    89,   310,    89,   310,    89,   449,    89,
     449,    89,   310,   449,    89,   310,   531,   263,    61,    98,
    -454,   538,   263,    68,   141,  -454,  -454,   532,   534,   516,
     -21,  -454,   527,   258,  -454,    89,    89,   449,    89,    21,
    -454,  1159,   542,   310,  -454,  1159,   -42,   603,  -454,  -454,
     605,   541,  -454,  -454,  -454,  -454,  -454,  -454,  -454,   543,
    -454,   544,  -454,   545,   430,   -36,   546,   310,   547,  -454,
     310,   548,  -454,   529,   187,   551,   549,   310,   530,   195,
     552,   557,   310,    89,    89,   527,   606,   516,  -454,  -454,
    -454,   550,   554,   310,   559,   561,  -454,  -454,   558,  -454,
    -454,   564,   565,   449,   449,   449,   449,   556,  -454,   449,
    -454,   310,  -454,   449,   112,   310,   574,   263,  -454,   125,
     310,   576,   263,  -454,  -454,  -454,   569,   527,  -454,  -454,
    -454,  1159,   711,  -454,  1159,  1159,   310,   310,   310,   310,
     582,   310,  -454,   310,   583,   310,  -454,   585,   224,   589,
     310,  -454,   590,   229,   275,  -454,  -454,   578,   579,   584,
    -454,  -454,  -454,  -454,   587,  -454,  -454,   159,  -454,   166,
     310,   591,  -454,   175,   310,   -18,  -454,  -454,  -454,   310,
     598,   310,   601,   310,  -454,   310,   602,   310,  -454,  -454,
    -454,   595,  -454,   596,  -454,  -454,   597,  -454,   310,   310,
     310,  -454,  -454,  -454
};

/* YYPGOTO[NTERM-NUM].  */
static const grib_yytype_int16 grib_yypgoto[] =
{
    -454,  -454,     7,  -454,     0,  -454,  -453,   436,  -212,  -143,
    -454,  -454,  -454,  -307,    49,  -106,    74,  -454,   299,  -454,
    -454,  -454,  -454,  -178,  -454,   219,  -454,   -66,  -241,  -454,
     -92,   573,   113,   131,   111,  -108,   475,  -454,   -37,  -454,
      12,   265,  -454,  -454,   629
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -8
static const grib_yytype_int16 grib_yytable[] =
{
      59,   399,   165,   318,   245,   210,   166,    58,   150,   150,
     355,   356,    72,   547,    36,   360,   525,   111,   112,   108,
     288,   396,   397,   291,   191,   154,   109,   213,   396,   397,
     309,   396,   397,   178,   192,   269,   155,   179,   180,   169,
     170,    76,   120,   167,   225,   226,   227,   228,   255,   177,
     283,   284,   111,   112,   223,   473,   385,   390,    36,   224,
     161,   386,   391,   202,   474,    79,   205,    77,   113,   398,
     114,   115,   116,   117,   567,   229,   230,   526,   118,   434,
     629,   212,   310,    72,   180,   119,    78,   167,   120,   151,
     261,   111,   112,   478,   181,   535,   163,    54,    55,   164,
      84,    85,    80,   113,    82,   114,   115,   116,   117,   389,
     357,   246,   421,   118,   250,   319,   115,   249,    81,   536,
     119,   111,   112,   335,   320,   213,   213,   213,   213,   213,
     213,   213,   213,   213,   213,   213,   213,   213,   213,   111,
     112,    83,   113,   443,   114,   115,   116,   117,   180,   445,
     515,    87,   118,   315,   180,   516,   120,   520,    86,   119,
     258,   181,   521,   527,   372,   163,    54,    55,   164,   311,
     312,    93,   113,    88,   114,   115,   116,   117,   394,    89,
     218,   219,   118,   163,    54,    55,   164,   120,    90,   119,
     113,    91,   114,   115,   116,   117,   517,   188,   189,    92,
     447,   584,   120,   220,   221,   222,   585,   119,   253,   234,
     459,   350,   254,   180,   589,   464,   120,   387,   120,   590,
     361,   388,   165,   362,   250,   165,   163,    54,    55,   164,
     392,    94,   448,   450,   393,    95,   111,   112,    96,   522,
     120,    97,   120,   120,   120,   111,   112,    98,   620,   120,
     120,   120,   425,   621,    99,   622,   272,   367,   368,   435,
     623,   437,   438,   439,   626,   441,   534,   444,   446,   627,
     477,   100,   163,    54,    55,   164,   101,   181,   455,   514,
     163,    54,    55,   164,   519,   555,   102,   113,   120,   114,
     115,   116,   117,   560,   400,   401,   103,   449,   114,   115,
     116,   117,   104,   196,   119,   197,   118,   198,   527,   163,
      54,    55,   164,   119,   163,    54,    55,   164,   483,   199,
     105,   200,   610,   201,   396,   397,   106,   614,   107,   359,
     321,   322,   323,   324,   325,   326,   329,   330,   331,   332,
     333,   334,   505,   133,   507,   529,   530,   510,   163,    54,
      55,   164,   110,   379,   327,   328,   134,   382,   383,   135,
     140,   136,   137,   138,   433,   181,   436,   139,   144,   145,
     440,   533,   181,   146,   181,   181,   181,   141,   181,   588,
     181,   181,   482,   292,   593,   142,   143,   147,   152,   148,
     149,   181,   156,   165,   153,   162,   171,   172,   165,   173,
     193,   175,   176,   183,   184,   185,   186,   187,   293,   202,
     504,   208,   504,   215,   195,   294,   295,   296,   297,   298,
     299,   300,   301,   302,   303,   190,   304,   305,   194,   206,
     207,   181,   209,   211,   216,   231,   217,   576,   577,   578,
     579,   120,   232,   581,   233,   235,   442,   583,   165,   236,
     120,   237,   120,   165,   120,   181,   120,   181,   120,   238,
     181,   120,   239,   460,   240,   241,   242,   243,   465,   244,
     150,   154,   248,   270,   395,   278,   279,   257,   251,   252,
     256,   538,   120,   120,   181,   259,   564,   565,   260,   155,
     261,   262,   263,   264,   265,   266,   267,   273,   280,   484,
     268,   274,   275,   276,   277,   282,   286,   281,   287,   492,
     289,   494,   495,   496,   290,   498,   308,   174,   500,   451,
     502,   313,   165,   341,   317,   342,   509,   165,   348,   512,
     343,   344,   351,   354,   349,   358,   365,   369,   249,   363,
     371,   364,   366,   370,   373,   180,   374,   375,   376,   377,
     181,   181,   181,   181,   378,   381,   181,   540,   403,   384,
     181,   404,   405,   406,   410,   407,   411,   408,   409,   412,
     413,   595,   597,   417,   598,   599,   414,   415,   420,   548,
     416,   550,   418,   419,   552,   247,   424,   422,   423,   426,
     427,   558,   428,   452,   453,   454,   563,   456,   457,   458,
     462,   461,   463,   467,   468,   475,   471,   570,   478,   541,
     481,   542,   566,   476,   554,   559,   271,   528,   485,   486,
     513,   487,   489,   490,   493,   582,   497,   518,   523,   586,
     524,   285,   539,   543,   591,   544,   545,   556,   561,   546,
     549,   551,   553,   615,   568,   314,   557,   316,   569,   580,
     600,   601,   602,   603,   562,   605,   573,   606,   571,   608,
     572,   574,   575,   587,   612,   592,   594,   604,   607,   337,
     609,   338,   339,   340,   611,   613,   616,   617,   345,   346,
     347,   619,   618,   631,   624,   625,   633,   636,   628,   638,
     639,   640,   214,   630,   466,   632,   488,   634,   537,   635,
     168,   637,     0,     0,     0,     0,     0,   336,     0,     0,
       0,     0,   641,   642,   643,   157,     3,   380,     0,     4,
       5,     6,     7,     8,     9,    10,     0,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,     0,     0,    27,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    28,
      29,    30,    31,    32,    33,    34,    35,     0,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,     0,    52,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   158,     0,     0,   159,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   596,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     491,     0,     0,     0,     0,     0,     0,     0,     0,   499,
       0,   501,     0,   503,     0,   506,     0,   508,     0,     0,
     511,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    -7,     1,     0,     0,     2,
       3,   531,   532,     4,     5,     6,     7,     8,     9,    10,
       0,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,     0,     0,    27,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,     0,    52,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      53,    54,    55,    56,   157,     3,     0,     0,     4,     5,
       6,     7,     8,     9,    10,     0,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,     0,     0,    27,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    28,    29,
      30,    31,    32,    33,    34,    35,     0,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      50,    51,     0,    52,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   158,     0,     0,   159,     0,
     160,   157,     3,     0,     0,     4,     5,     6,     7,     8,
       9,    10,     0,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,     0,
       0,    27,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    51,     0,
      52,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   429,   157,     3,   159,     0,     4,     5,     6,
       7,     8,     9,    10,     0,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,     0,     0,    27,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    28,    29,    30,
      31,    32,    33,    34,    35,     0,    37,    38,    39,    40,
      41,    42,    43,    44,    45,    46,    47,    48,    49,    50,
      51,     0,    52,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   158,     0,     0,   159
};

static const grib_yytype_int16 grib_yycheck[] =
{
       0,   308,    68,   215,   147,   113,     4,     0,    15,    15,
     251,   252,     0,    49,    56,   256,   469,    32,    33,    86,
     198,    49,    50,   201,    85,    85,    93,   119,    49,    50,
      45,    49,    50,    91,    95,   178,    96,    95,    96,    76,
      77,    93,    35,    85,    75,    76,    77,    78,   154,    86,
     193,   194,    32,    33,   100,    89,    89,    89,    56,   105,
      60,    94,    94,    99,    98,    85,   103,    93,    83,    97,
      85,    86,    87,    88,   527,   106,   107,    98,    93,    91,
      98,   118,    97,    71,    96,   100,    91,    85,    81,    96,
      96,    32,    33,    72,    87,    74,    85,    86,    87,    88,
      85,    86,    85,    83,    91,    85,    86,    87,    88,    98,
     253,   148,   353,    93,   151,    85,    86,    97,    93,    98,
     100,    32,    33,   231,   216,   217,   218,   219,   220,   221,
     222,   223,   224,   225,   226,   227,   228,   229,   230,    32,
      33,    91,    83,    91,    85,    86,    87,    88,    96,    91,
      89,    85,    93,    94,    96,    94,   149,    89,    93,   100,
     160,   154,    94,   470,   270,    85,    86,    87,    88,   206,
     207,    85,    83,    91,    85,    86,    87,    88,    98,    91,
      79,    80,    93,    85,    86,    87,    88,   180,    91,   100,
      83,    91,    85,    86,    87,    88,    98,    95,    96,    91,
      93,    89,   195,   102,   103,   104,    94,   100,    91,   135,
     388,   248,    95,    96,    89,   393,   209,    93,   211,    94,
     257,    97,   288,   260,   261,   291,    85,    86,    87,    88,
      93,    85,   375,   376,    97,    85,    32,    33,    85,    98,
     233,    85,   235,   236,   237,    32,    33,    85,    89,   242,
     243,   244,   358,    94,    85,    89,   182,    85,    86,   365,
      94,   367,   368,   369,    89,   371,   478,   373,   374,    94,
     413,    85,    85,    86,    87,    88,    85,   270,   384,   457,
      85,    86,    87,    88,   462,    98,    85,    83,   281,    85,
      86,    87,    88,    98,    87,    88,    93,    93,    85,    86,
      87,    88,    85,    93,   100,    95,    93,    97,   615,    85,
      86,    87,    88,   100,    85,    86,    87,    88,   424,    93,
      85,    95,    98,    97,    49,    50,    86,    98,    86,   255,
     217,   218,   219,   220,   221,   222,   225,   226,   227,   228,
     229,   230,   448,    90,   450,    87,    88,   453,    85,    86,
      87,    88,    85,   279,   223,   224,    85,   283,   284,    85,
      91,    85,    85,    85,   364,   358,   366,    85,    85,    85,
     370,   477,   365,    85,   367,   368,   369,    91,   371,   557,
     373,   374,   419,     3,   562,    91,    91,    85,    96,    93,
      93,   384,     0,   459,    96,    90,    87,    86,   464,    86,
      96,    87,    87,    87,    87,    87,    87,    87,    28,    99,
     447,    86,   449,    89,    93,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    95,    46,    47,    96,    96,
      96,   424,    93,    93,    84,    81,   101,   543,   544,   545,
     546,   434,    82,   549,    93,    93,   372,   553,   514,    93,
     443,    93,   445,   519,   447,   448,   449,   450,   451,    87,
     453,   454,    87,   389,    87,    87,    93,    93,   394,    93,
      15,    85,    93,    85,    89,    85,    85,    93,    97,    97,
      97,   481,   475,   476,   477,    96,   523,   524,    93,    96,
      96,    94,    94,    92,    94,    92,    92,    92,    85,   425,
      94,    92,    92,    92,    92,    85,    85,    93,    85,   435,
      85,   437,   438,   439,    85,   441,    94,    81,   444,    93,
     446,    94,   588,    92,    94,    92,   452,   593,    94,   455,
      92,    92,    85,    90,    94,    85,    85,    85,    97,    97,
      92,    97,    97,    97,    85,    96,    85,    85,    85,    85,
     543,   544,   545,   546,    96,    85,   549,   483,    86,    94,
     553,    94,    94,    94,    85,    94,    85,    94,    94,    85,
      85,   571,   572,    97,   574,   575,    94,    94,    98,   505,
      94,   507,    94,    96,   510,   149,    92,    98,    98,    98,
      94,   517,    94,    85,    94,    93,   522,    86,    97,    85,
      97,    86,    85,    85,    85,    93,    90,   533,    72,     6,
      97,     6,     6,    93,    85,    85,   180,    90,    97,    97,
      89,    98,    98,    98,    98,   551,    98,    89,    96,   555,
      96,   195,    90,    92,   560,    92,    92,    86,    86,    94,
      94,    94,    94,   594,    94,   209,    97,   211,    94,    93,
     576,   577,   578,   579,    97,   581,    98,   583,    99,   585,
      99,    97,    97,    89,   590,    89,    97,    85,    85,   233,
      85,   235,   236,   237,    85,    85,    98,    98,   242,   243,
     244,    94,    98,    85,   610,    94,    85,    85,   614,    94,
      94,    94,   119,   619,   395,   621,   431,   623,   479,   625,
      71,   627,    -1,    -1,    -1,    -1,    -1,   232,    -1,    -1,
      -1,    -1,   638,   639,   640,     4,     5,   281,    -1,     8,
       9,    10,    11,    12,    13,    14,    -1,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    -1,    -1,    34,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    48,
      49,    50,    51,    52,    53,    54,    55,    -1,    57,    58,
      59,    60,    61,    62,    63,    64,    65,    66,    67,    68,
      69,    70,    71,    -1,    73,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    85,    -1,    -1,    88,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    98,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     434,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   443,
      -1,   445,    -1,   447,    -1,   449,    -1,   451,    -1,    -1,
     454,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,     0,     1,    -1,    -1,     4,
       5,   475,   476,     8,     9,    10,    11,    12,    13,    14,
      -1,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    -1,    -1,    34,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    -1,    73,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      85,    86,    87,    88,     4,     5,    -1,    -1,     8,     9,
      10,    11,    12,    13,    14,    -1,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    -1,    -1,    34,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    48,    49,
      50,    51,    52,    53,    54,    55,    -1,    57,    58,    59,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    69,
      70,    71,    -1,    73,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    85,    -1,    -1,    88,    -1,
      90,     4,     5,    -1,    -1,     8,     9,    10,    11,    12,
      13,    14,    -1,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    -1,
      -1,    34,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    48,    49,    50,    51,    52,
      53,    54,    55,    56,    57,    58,    59,    60,    61,    62,
      63,    64,    65,    66,    67,    68,    69,    70,    71,    -1,
      73,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    85,     4,     5,    88,    -1,     8,     9,    10,
      11,    12,    13,    14,    -1,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    -1,    -1,    34,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    48,    49,    50,
      51,    52,    53,    54,    55,    -1,    57,    58,    59,    60,
      61,    62,    63,    64,    65,    66,    67,    68,    69,    70,
      71,    -1,    73,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    85,    -1,    -1,    88
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const grib_yytype_uint8 grib_yystos[] =
{
       0,     1,     4,     5,     8,     9,    10,    11,    12,    13,
      14,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    34,    48,    49,
      50,    51,    52,    53,    54,    55,    56,    57,    58,    59,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    69,
      70,    71,    73,    85,    86,    87,    88,   109,   110,   112,
     113,   118,   119,   120,   127,   128,   129,   130,   131,   134,
     135,   147,   148,   150,   151,   152,    93,    93,    91,    85,
      85,    93,    91,    91,    85,    86,    93,    85,    91,    91,
      91,    91,    91,    85,    85,    85,    85,    85,    85,    85,
      85,    85,    85,    93,    85,    85,    86,    86,    86,    93,
      85,    32,    33,    83,    85,    86,    87,    88,    93,   100,
     110,   115,   116,   117,   138,   139,   140,   141,   142,   143,
     144,   145,   146,    90,    85,    85,    85,    85,    85,    85,
      91,    91,    91,    91,    85,    85,    85,    85,    93,    93,
      15,    96,    96,    96,    85,    96,     0,     4,    85,    88,
      90,   112,    90,    85,    88,   135,     4,    85,   152,   146,
     146,    87,    86,    86,   115,    87,    87,   146,    91,    95,
      96,   110,   123,    87,    87,    87,    87,    87,    95,    96,
      95,    85,    95,    96,    96,    93,    93,    95,    97,    93,
      95,    97,    99,   110,   124,   146,    96,    96,    86,    93,
     143,    93,   146,   138,   139,    89,    84,   101,    79,    80,
     102,   103,   104,   100,   105,    75,    76,    77,    78,   106,
     107,    81,    82,    93,   124,    93,    93,    93,    87,    87,
      87,    87,    93,    93,    93,   117,   146,   115,    93,    97,
     146,    97,    97,    91,    95,   123,    97,    93,   112,    96,
      93,    96,    94,    94,    92,    94,    92,    92,    94,   117,
      85,   115,   124,    92,    92,    92,    92,    92,    85,    85,
      85,    93,    85,   117,   117,   115,    85,    85,   131,    85,
      85,   131,     3,    28,    35,    36,    37,    38,    39,    40,
      41,    42,    43,    44,    46,    47,   125,   126,    94,    45,
      97,   146,   146,    94,   115,    94,   115,    94,   116,    85,
     138,   140,   140,   140,   140,   140,   140,   141,   141,   142,
     142,   142,   142,   142,   142,   143,   144,   115,   115,   115,
     115,    92,    92,    92,    92,   115,   115,   115,    94,    94,
     146,    85,   136,   137,    90,   136,   136,   117,    85,   124,
     136,   146,   146,    97,    97,    85,    97,    85,    86,    85,
      97,    92,   123,    85,    85,    85,    85,    85,    96,   124,
     115,    85,   124,   124,    94,    89,    94,    93,    97,    98,
      89,    94,    93,    97,    98,    89,    49,    50,    97,   121,
      87,    88,   111,    86,    94,    94,    94,    94,    94,    94,
      85,    85,    85,    85,    94,    94,    94,    97,    94,    96,
      98,   136,    98,    98,    92,   123,    98,    94,    94,    85,
     112,   148,   149,   112,    91,   123,   112,   123,   123,   123,
     112,   123,   124,    91,   123,    91,   123,    93,   117,    93,
     117,    93,    85,    94,    93,   123,    86,    97,    85,   131,
     124,    86,    97,    85,   131,   124,   126,    85,    85,   121,
     122,    90,   114,    89,    98,    93,    93,   117,    72,   132,
     133,    97,   146,   123,   124,    97,    97,    98,   149,    98,
      98,   115,   124,    98,   124,   124,   124,    98,   124,   115,
     124,   115,   124,   115,   146,   123,   115,   123,   115,   124,
     123,   115,   124,    89,   131,    89,    94,    98,    89,   131,
      89,    94,    98,    96,    96,   114,    98,   121,    90,    87,
      88,   115,   115,   123,   116,    74,    98,   133,   112,    90,
     124,     6,     6,    92,    92,    92,    94,    49,   124,    94,
     124,    94,   124,    94,    85,    98,    86,    97,   124,    85,
      98,    86,    97,   124,   146,   146,     6,   114,    94,    94,
     124,    99,    99,    98,    97,    97,   123,   123,   123,   123,
      93,   123,   124,   123,    89,    94,   124,    89,   131,    89,
      94,   124,    89,   131,    97,   112,    98,   112,   112,   112,
     124,   124,   124,   124,    85,   124,   124,    85,   124,    85,
      98,    85,   124,    85,    98,   122,    98,    98,    98,    94,
      89,    94,    89,    94,   124,    94,    89,    94,   124,    98,
     124,    85,   124,    85,   124,   124,    85,   124,    94,    94,
      94,   124,   124,   124
};

#define grib_yyerrok		(grib_yyerrstatus = 0)
#define grib_yyclearin	(grib_yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto grib_yyacceptlab
#define YYABORT		goto grib_yyabortlab
#define YYERROR		goto grib_yyerrorlab


/* Like YYERROR except do call grib_yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto grib_yyerrlab

#define YYRECOVERING()  (!!grib_yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (grib_yychar == YYEMPTY && grib_yylen == 1)				\
    {								\
      grib_yychar = (Token);						\
      grib_yylval = (Value);						\
      grib_yytoken = YYTRANSLATE (grib_yychar);				\
      YYPOPSTACK (1);						\
      goto grib_yybackup;						\
    }								\
  else								\
    {								\
      grib_yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `grib_yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX grib_yylex (YYLEX_PARAM)
#else
# define YYLEX grib_yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (grib_yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (grib_yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      grib_yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
grib_yy_symbol_value_print (FILE *grib_yyoutput, int grib_yytype, YYSTYPE const * const grib_yyvaluep)
#else
static void
grib_yy_symbol_value_print (grib_yyoutput, grib_yytype, grib_yyvaluep)
    FILE *grib_yyoutput;
    int grib_yytype;
    YYSTYPE const * const grib_yyvaluep;
#endif
{
  if (!grib_yyvaluep)
    return;
# ifdef YYPRINT
  if (grib_yytype < YYNTOKENS)
    YYPRINT (grib_yyoutput, grib_yytoknum[grib_yytype], *grib_yyvaluep);
# else
  YYUSE (grib_yyoutput);
# endif
  switch (grib_yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
grib_yy_symbol_print (FILE *grib_yyoutput, int grib_yytype, YYSTYPE const * const grib_yyvaluep)
#else
static void
grib_yy_symbol_print (grib_yyoutput, grib_yytype, grib_yyvaluep)
    FILE *grib_yyoutput;
    int grib_yytype;
    YYSTYPE const * const grib_yyvaluep;
#endif
{
  if (grib_yytype < YYNTOKENS)
    YYFPRINTF (grib_yyoutput, "token %s (", grib_yytname[grib_yytype]);
  else
    YYFPRINTF (grib_yyoutput, "nterm %s (", grib_yytname[grib_yytype]);

  grib_yy_symbol_value_print (grib_yyoutput, grib_yytype, grib_yyvaluep);
  YYFPRINTF (grib_yyoutput, ")");
}

/*------------------------------------------------------------------.
| grib_yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
grib_yy_stack_print (grib_yytype_int16 *bottom, grib_yytype_int16 *top)
#else
static void
grib_yy_stack_print (bottom, top)
    grib_yytype_int16 *bottom;
    grib_yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (grib_yydebug)							\
    grib_yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
grib_yy_reduce_print (YYSTYPE *grib_yyvsp, int grib_yyrule)
#else
static void
grib_yy_reduce_print (grib_yyvsp, grib_yyrule)
    YYSTYPE *grib_yyvsp;
    int grib_yyrule;
#endif
{
  int grib_yynrhs = grib_yyr2[grib_yyrule];
  int grib_yyi;
  unsigned long int grib_yylno = grib_yyrline[grib_yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     grib_yyrule - 1, grib_yylno);
  /* The symbols being reduced.  */
  for (grib_yyi = 0; grib_yyi < grib_yynrhs; grib_yyi++)
    {
      fprintf (stderr, "   $%d = ", grib_yyi + 1);
      grib_yy_symbol_print (stderr, grib_yyrhs[grib_yyprhs[grib_yyrule] + grib_yyi],
		       &(grib_yyvsp[(grib_yyi + 1) - (grib_yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (grib_yydebug)				\
    grib_yy_reduce_print (grib_yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int grib_yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef grib_yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define grib_yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
grib_yystrlen (const char *grib_yystr)
#else
static YYSIZE_T
grib_yystrlen (grib_yystr)
    const char *grib_yystr;
#endif
{
  YYSIZE_T grib_yylen;
  for (grib_yylen = 0; grib_yystr[grib_yylen]; grib_yylen++)
    continue;
  return grib_yylen;
}
#  endif
# endif

# ifndef grib_yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define grib_yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
grib_yystpcpy (char *grib_yydest, const char *grib_yysrc)
#else
static char *
grib_yystpcpy (grib_yydest, grib_yysrc)
    char *grib_yydest;
    const char *grib_yysrc;
#endif
{
  char *grib_yyd = grib_yydest;
  const char *grib_yys = grib_yysrc;

  while ((*grib_yyd++ = *grib_yys++) != '\0')
    continue;

  return grib_yyd - 1;
}
#  endif
# endif

# ifndef grib_yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for grib_yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from grib_yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
grib_yytnamerr (char *grib_yyres, const char *grib_yystr)
{
  if (*grib_yystr == '"')
    {
      YYSIZE_T grib_yyn = 0;
      char const *grib_yyp = grib_yystr;

      for (;;)
	switch (*++grib_yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++grib_yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (grib_yyres)
	      grib_yyres[grib_yyn] = *grib_yyp;
	    grib_yyn++;
	    break;

	  case '"':
	    if (grib_yyres)
	      grib_yyres[grib_yyn] = '\0';
	    return grib_yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! grib_yyres)
    return grib_yystrlen (grib_yystr);

  return grib_yystpcpy (grib_yyres, grib_yystr) - grib_yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
grib_yysyntax_error (char *grib_yyresult, int grib_yystate, int grib_yychar)
{
  int grib_yyn = grib_yypact[grib_yystate];

  if (! (YYPACT_NINF < grib_yyn && grib_yyn <= YYLAST))
    return 0;
  else
    {
      int grib_yytype = YYTRANSLATE (grib_yychar);
      YYSIZE_T grib_yysize0 = grib_yytnamerr (0, grib_yytname[grib_yytype]);
      YYSIZE_T grib_yysize = grib_yysize0;
      YYSIZE_T grib_yysize1;
      int grib_yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *grib_yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int grib_yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *grib_yyfmt;
      char const *grib_yyf;
      static char const grib_yyunexpected[] = "syntax error, unexpected %s";
      static char const grib_yyexpecting[] = ", expecting %s";
      static char const grib_yyor[] = " or %s";
      char grib_yyformat[sizeof grib_yyunexpected
		    + sizeof grib_yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof grib_yyor - 1))];
      char const *grib_yyprefix = grib_yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int grib_yyxbegin = grib_yyn < 0 ? -grib_yyn : 0;

      /* Stay within bounds of both grib_yycheck and grib_yytname.  */
      int grib_yychecklim = YYLAST - grib_yyn + 1;
      int grib_yyxend = grib_yychecklim < YYNTOKENS ? grib_yychecklim : YYNTOKENS;
      int grib_yycount = 1;

      grib_yyarg[0] = grib_yytname[grib_yytype];
      grib_yyfmt = grib_yystpcpy (grib_yyformat, grib_yyunexpected);

      for (grib_yyx = grib_yyxbegin; grib_yyx < grib_yyxend; ++grib_yyx)
	if (grib_yycheck[grib_yyx + grib_yyn] == grib_yyx && grib_yyx != YYTERROR)
	  {
	    if (grib_yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		grib_yycount = 1;
		grib_yysize = grib_yysize0;
		grib_yyformat[sizeof grib_yyunexpected - 1] = '\0';
		break;
	      }
	    grib_yyarg[grib_yycount++] = grib_yytname[grib_yyx];
	    grib_yysize1 = grib_yysize + grib_yytnamerr (0, grib_yytname[grib_yyx]);
	    grib_yysize_overflow |= (grib_yysize1 < grib_yysize);
	    grib_yysize = grib_yysize1;
	    grib_yyfmt = grib_yystpcpy (grib_yyfmt, grib_yyprefix);
	    grib_yyprefix = grib_yyor;
	  }

      grib_yyf = YY_(grib_yyformat);
      grib_yysize1 = grib_yysize + grib_yystrlen (grib_yyf);
      grib_yysize_overflow |= (grib_yysize1 < grib_yysize);
      grib_yysize = grib_yysize1;

      if (grib_yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (grib_yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *grib_yyp = grib_yyresult;
	  int grib_yyi = 0;
	  while ((*grib_yyp = *grib_yyf) != '\0')
	    {
	      if (*grib_yyp == '%' && grib_yyf[1] == 's' && grib_yyi < grib_yycount)
		{
		  grib_yyp += grib_yytnamerr (grib_yyp, grib_yyarg[grib_yyi++]);
		  grib_yyf += 2;
		}
	      else
		{
		  grib_yyp++;
		  grib_yyf++;
		}
	    }
	}
      return grib_yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
grib_yydestruct (const char *grib_yymsg, int grib_yytype, YYSTYPE *grib_yyvaluep)
#else
static void
grib_yydestruct (grib_yymsg, grib_yytype, grib_yyvaluep)
    const char *grib_yymsg;
    int grib_yytype;
    YYSTYPE *grib_yyvaluep;
#endif
{
  YYUSE (grib_yyvaluep);

  if (!grib_yymsg)
    grib_yymsg = "Deleting";
  YY_SYMBOL_PRINT (grib_yymsg, grib_yytype, grib_yyvaluep, grib_yylocationp);

  switch (grib_yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int grib_yyparse (void *YYPARSE_PARAM);
#else
int grib_yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int grib_yyparse (void);
#else
int grib_yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int grib_yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE grib_yylval;

/* Number of syntax errors so far.  */
int grib_yynerrs;



/*----------.
| grib_yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
grib_yyparse (void *YYPARSE_PARAM)
#else
int
grib_yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
grib_yyparse (void)
#else
int
grib_yyparse ()

#endif
#endif
{
  
  int grib_yystate;
  int grib_yyn;
  int grib_yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int grib_yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int grib_yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char grib_yymsgbuf[128];
  char *grib_yymsg = grib_yymsgbuf;
  YYSIZE_T grib_yymsg_alloc = sizeof grib_yymsgbuf;
#endif

  /* Three stacks and their tools:
     `grib_yyss': related to states,
     `grib_yyvs': related to semantic values,
     `grib_yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow grib_yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  grib_yytype_int16 grib_yyssa[YYINITDEPTH];
  grib_yytype_int16 *grib_yyss = grib_yyssa;
  grib_yytype_int16 *grib_yyssp;

  /* The semantic value stack.  */
  YYSTYPE grib_yyvsa[YYINITDEPTH];
  YYSTYPE *grib_yyvs = grib_yyvsa;
  YYSTYPE *grib_yyvsp;



#define YYPOPSTACK(N)   (grib_yyvsp -= (N), grib_yyssp -= (N))

  YYSIZE_T grib_yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE grib_yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int grib_yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  grib_yystate = 0;
  grib_yyerrstatus = 0;
  grib_yynerrs = 0;
  grib_yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  grib_yyssp = grib_yyss;
  grib_yyvsp = grib_yyvs;

  goto grib_yysetstate;

/*------------------------------------------------------------.
| grib_yynewstate -- Push a new state, which is found in grib_yystate.  |
`------------------------------------------------------------*/
 grib_yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  grib_yyssp++;

 grib_yysetstate:
  *grib_yyssp = grib_yystate;

  if (grib_yyss + grib_yystacksize - 1 <= grib_yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T grib_yysize = grib_yyssp - grib_yyss + 1;

#ifdef grib_yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *grib_yyvs1 = grib_yyvs;
	grib_yytype_int16 *grib_yyss1 = grib_yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if grib_yyoverflow is a macro.  */
	grib_yyoverflow (YY_("memory exhausted"),
		    &grib_yyss1, grib_yysize * sizeof (*grib_yyssp),
		    &grib_yyvs1, grib_yysize * sizeof (*grib_yyvsp),

		    &grib_yystacksize);

	grib_yyss = grib_yyss1;
	grib_yyvs = grib_yyvs1;
      }
#else /* no grib_yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto grib_yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= grib_yystacksize)
	goto grib_yyexhaustedlab;
      grib_yystacksize *= 2;
      if (YYMAXDEPTH < grib_yystacksize)
	grib_yystacksize = YYMAXDEPTH;

      {
	grib_yytype_int16 *grib_yyss1 = grib_yyss;
	union grib_yyalloc *grib_yyptr =
	  (union grib_yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (grib_yystacksize));
	if (! grib_yyptr)
	  goto grib_yyexhaustedlab;
	YYSTACK_RELOCATE (grib_yyss);
	YYSTACK_RELOCATE (grib_yyvs);

#  undef YYSTACK_RELOCATE
	if (grib_yyss1 != grib_yyssa)
	  YYSTACK_FREE (grib_yyss1);
      }
# endif
#endif /* no grib_yyoverflow */

      grib_yyssp = grib_yyss + grib_yysize - 1;
      grib_yyvsp = grib_yyvs + grib_yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) grib_yystacksize));

      if (grib_yyss + grib_yystacksize - 1 <= grib_yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", grib_yystate));

  goto grib_yybackup;

/*-----------.
| grib_yybackup.  |
`-----------*/
grib_yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  grib_yyn = grib_yypact[grib_yystate];
  if (grib_yyn == YYPACT_NINF)
    goto grib_yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (grib_yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      grib_yychar = YYLEX;
    }

  if (grib_yychar <= YYEOF)
    {
      grib_yychar = grib_yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      grib_yytoken = YYTRANSLATE (grib_yychar);
      YY_SYMBOL_PRINT ("Next token is", grib_yytoken, &grib_yylval, &grib_yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  grib_yyn += grib_yytoken;
  if (grib_yyn < 0 || YYLAST < grib_yyn || grib_yycheck[grib_yyn] != grib_yytoken)
    goto grib_yydefault;
  grib_yyn = grib_yytable[grib_yyn];
  if (grib_yyn <= 0)
    {
      if (grib_yyn == 0 || grib_yyn == YYTABLE_NINF)
	goto grib_yyerrlab;
      grib_yyn = -grib_yyn;
      goto grib_yyreduce;
    }

  if (grib_yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (grib_yyerrstatus)
    grib_yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", grib_yytoken, &grib_yylval, &grib_yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (grib_yychar != YYEOF)
    grib_yychar = YYEMPTY;

  grib_yystate = grib_yyn;
  *++grib_yyvsp = grib_yylval;

  goto grib_yynewstate;


/*-----------------------------------------------------------.
| grib_yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
grib_yydefault:
  grib_yyn = grib_yydefact[grib_yystate];
  if (grib_yyn == 0)
    goto grib_yyerrlab;
  goto grib_yyreduce;


/*-----------------------------.
| grib_yyreduce -- Do a reduction.  |
`-----------------------------*/
grib_yyreduce:
  /* grib_yyn is the number of a rule to reduce with.  */
  grib_yylen = grib_yyr2[grib_yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  grib_yyval = grib_yyvsp[1-grib_yylen];


  YY_REDUCE_PRINT (grib_yyn);
  switch (grib_yyn)
    {
        case 2:
#line 204 "griby.y"
    { grib_parser_all_actions = 0;grib_parser_concept=0; grib_parser_rules=0; }
    break;

  case 3:
#line 205 "griby.y"
    { grib_parser_concept     = reverse((grib_yyvsp[(1) - (1)].concept_value)); }
    break;

  case 4:
#line 206 "griby.y"
    { grib_parser_all_actions = (grib_yyvsp[(1) - (1)].act); }
    break;

  case 5:
#line 207 "griby.y"
    { grib_parser_rules       = (grib_yyvsp[(1) - (1)].rules); }
    break;

  case 6:
#line 209 "griby.y"
    { grib_parser_all_actions = 0; grib_parser_concept=0; grib_parser_rules=0; }
    break;

  case 8:
#line 215 "griby.y"
    { (grib_yyval.dvalue)=grib_darray_push(0,0,(grib_yyvsp[(1) - (1)].dval));}
    break;

  case 9:
#line 216 "griby.y"
    { (grib_yyval.dvalue)=grib_darray_push(0,(grib_yyvsp[(1) - (3)].dvalue),(grib_yyvsp[(3) - (3)].dval));}
    break;

  case 10:
#line 217 "griby.y"
    { (grib_yyval.dvalue)=grib_darray_push(0,0,(grib_yyvsp[(1) - (1)].lval));}
    break;

  case 11:
#line 218 "griby.y"
    { (grib_yyval.dvalue)=grib_darray_push(0,(grib_yyvsp[(1) - (3)].dvalue),(grib_yyvsp[(3) - (3)].lval));}
    break;

  case 13:
#line 221 "griby.y"
    { (grib_yyvsp[(1) - (2)].act)->next = (grib_yyvsp[(2) - (2)].act); (grib_yyval.act) = (grib_yyvsp[(1) - (2)].act); }
    break;

  case 14:
#line 222 "griby.y"
    { (grib_yyvsp[(1) - (3)].act)->next = (grib_yyvsp[(3) - (3)].act); (grib_yyval.act) = (grib_yyvsp[(1) - (3)].act); }
    break;

  case 15:
#line 223 "griby.y"
    {  (grib_yyval.act) = (grib_yyvsp[(1) - (2)].act);}
    break;

  case 26:
#line 241 "griby.y"
    { (grib_yyval.explist) = 0; }
    break;

  case 29:
#line 246 "griby.y"
    { (grib_yyvsp[(1) - (3)].explist)->next = (grib_yyvsp[(3) - (3)].explist); (grib_yyval.explist) = (grib_yyvsp[(1) - (3)].explist); }
    break;

  case 30:
#line 249 "griby.y"
    { (grib_yyval.explist) = grib_arguments_new(grib_parser_context,(grib_yyvsp[(1) - (1)].exp),NULL); }
    break;

  case 31:
#line 254 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (7)].str),"unsigned",(grib_yyvsp[(3) - (7)].lval),NULL,(grib_yyvsp[(6) - (7)].explist),(grib_yyvsp[(7) - (7)].lval),NULL,NULL);        free((grib_yyvsp[(5) - (7)].str));  }
    break;

  case 32:
#line 257 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (10)].str),"unsigned",(grib_yyvsp[(3) - (10)].lval),(grib_yyvsp[(7) - (10)].explist),(grib_yyvsp[(9) - (10)].explist),(grib_yyvsp[(10) - (10)].lval),NULL,NULL);        free((grib_yyvsp[(5) - (10)].str));  }
    break;

  case 33:
#line 260 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (7)].str),"ascii",(grib_yyvsp[(3) - (7)].lval),NULL,(grib_yyvsp[(6) - (7)].explist),(grib_yyvsp[(7) - (7)].lval),NULL,NULL);  free((grib_yyvsp[(5) - (7)].str));  }
    break;

  case 34:
#line 264 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (7)].str),"ascii",(grib_yyvsp[(3) - (7)].lval),NULL,(grib_yyvsp[(6) - (7)].explist),(grib_yyvsp[(7) - (7)].lval),NULL,NULL);  free((grib_yyvsp[(5) - (7)].str));  }
    break;

  case 35:
#line 267 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (7)].str),"bytes",(grib_yyvsp[(3) - (7)].lval),NULL,(grib_yyvsp[(6) - (7)].explist),(grib_yyvsp[(7) - (7)].lval),NULL,NULL);      free((grib_yyvsp[(5) - (7)].str));  }
    break;

  case 36:
#line 270 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (10)].str),"bytes",(grib_yyvsp[(3) - (10)].lval),(grib_yyvsp[(7) - (10)].explist),(grib_yyvsp[(9) - (10)].explist),(grib_yyvsp[(10) - (10)].lval),NULL,NULL);      free((grib_yyvsp[(5) - (10)].str));  }
    break;

  case 37:
#line 273 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (7)].str),"ksec1expver",(grib_yyvsp[(3) - (7)].lval),NULL,(grib_yyvsp[(6) - (7)].explist),(grib_yyvsp[(7) - (7)].lval),NULL,NULL);  free((grib_yyvsp[(5) - (7)].str));  }
    break;

  case 38:
#line 276 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (7)].str),"signed",(grib_yyvsp[(3) - (7)].lval),NULL,(grib_yyvsp[(6) - (7)].explist),(grib_yyvsp[(7) - (7)].lval),NULL,NULL);      free((grib_yyvsp[(5) - (7)].str));  }
    break;

  case 39:
#line 279 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (10)].str),"signed",(grib_yyvsp[(3) - (10)].lval),(grib_yyvsp[(7) - (10)].explist),(grib_yyvsp[(9) - (10)].explist),(grib_yyvsp[(10) - (10)].lval),NULL,NULL);      free((grib_yyvsp[(5) - (10)].str));  }
    break;

  case 40:
#line 282 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (8)].str),"codetable",(grib_yyvsp[(3) - (8)].lval), (grib_yyvsp[(6) - (8)].explist),(grib_yyvsp[(7) - (8)].explist),(grib_yyvsp[(8) - (8)].lval),NULL,NULL);    free((grib_yyvsp[(5) - (8)].str)); }
    break;

  case 41:
#line 285 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (12)].str),"codetable",(grib_yyvsp[(3) - (12)].lval), (grib_yyvsp[(6) - (12)].explist),(grib_yyvsp[(7) - (12)].explist),(grib_yyvsp[(12) - (12)].lval),NULL,(grib_yyvsp[(10) - (12)].str));
					free((grib_yyvsp[(5) - (12)].str));free((grib_yyvsp[(10) - (12)].str)); }
    break;

  case 42:
#line 289 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (10)].str),"codetable",(grib_yyvsp[(3) - (10)].lval), (grib_yyvsp[(7) - (10)].explist),(grib_yyvsp[(9) - (10)].explist),(grib_yyvsp[(10) - (10)].lval),NULL,NULL);    free((grib_yyvsp[(5) - (10)].str)); }
    break;

  case 43:
#line 292 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (8)].str),"complex_codetable",(grib_yyvsp[(3) - (8)].lval), (grib_yyvsp[(6) - (8)].explist),(grib_yyvsp[(7) - (8)].explist),(grib_yyvsp[(8) - (8)].lval),NULL,NULL);    free((grib_yyvsp[(5) - (8)].str)); }
    break;

  case 44:
#line 295 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (10)].str),"complex_codetable",(grib_yyvsp[(3) - (10)].lval), (grib_yyvsp[(7) - (10)].explist),(grib_yyvsp[(9) - (10)].explist),(grib_yyvsp[(10) - (10)].lval),NULL,NULL);    free((grib_yyvsp[(5) - (10)].str)); }
    break;

  case 45:
#line 298 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (8)].str),"codeflag",(grib_yyvsp[(3) - (8)].lval), (grib_yyvsp[(6) - (8)].explist),(grib_yyvsp[(7) - (8)].explist),(grib_yyvsp[(8) - (8)].lval),NULL,NULL);  free((grib_yyvsp[(5) - (8)].str)); }
    break;

  case 46:
#line 301 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (9)].str),"lookup",(grib_yyvsp[(3) - (9)].lval),(grib_yyvsp[(7) - (9)].explist),NULL,(grib_yyvsp[(9) - (9)].lval),NULL,NULL); free((grib_yyvsp[(5) - (9)].str)); }
    break;

  case 47:
#line 304 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (7)].str),"bit",0,(grib_yyvsp[(4) - (7)].explist),(grib_yyvsp[(6) - (7)].explist),(grib_yyvsp[(7) - (7)].lval),NULL,NULL); free((grib_yyvsp[(2) - (7)].str)); }
    break;

  case 48:
#line 307 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (2)].str),"label",0,NULL,NULL,0,NULL,NULL);   free((grib_yyvsp[(2) - (2)].str));  }
    break;

  case 49:
#line 310 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (2)].str),"label",0,NULL,NULL,0,NULL,NULL);   free((grib_yyvsp[(2) - (2)].str));  }
    break;

  case 50:
#line 313 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (4)].str),"ibmfloat",4,NULL,(grib_yyvsp[(3) - (4)].explist),(grib_yyvsp[(4) - (4)].lval),NULL,NULL);free((grib_yyvsp[(2) - (4)].str));  }
    break;

  case 51:
#line 316 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(4) - (6)].str),"ibmfloat",4,NULL,(grib_yyvsp[(5) - (6)].explist),(grib_yyvsp[(6) - (6)].lval),(grib_yyvsp[(2) - (6)].str),NULL);free((grib_yyvsp[(4) - (6)].str)); free((grib_yyvsp[(2) - (6)].str)); }
    break;

  case 52:
#line 319 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (7)].str),"ibmfloat",4,(grib_yyvsp[(4) - (7)].explist),(grib_yyvsp[(6) - (7)].explist),(grib_yyvsp[(7) - (7)].lval),NULL,NULL);free((grib_yyvsp[(2) - (7)].str));  }
    break;

  case 53:
#line 322 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (2)].str),"position",0,NULL,NULL,0,NULL,NULL);     free((grib_yyvsp[(2) - (2)].str));  }
    break;

  case 54:
#line 325 "griby.y"
    { (grib_yyval.act) = grib_action_create_variable(grib_parser_context,(grib_yyvsp[(2) - (5)].str),"constant",0,(grib_yyvsp[(4) - (5)].explist),NULL,(grib_yyvsp[(5) - (5)].lval),NULL);free((grib_yyvsp[(2) - (5)].str)); }
    break;

  case 55:
#line 328 "griby.y"
    { (grib_yyval.act) = grib_action_create_variable(grib_parser_context,(grib_yyvsp[(2) - (5)].str),"transient",0,(grib_yyvsp[(4) - (5)].explist),(grib_yyvsp[(4) - (5)].explist),(grib_yyvsp[(5) - (5)].lval),NULL);   free((grib_yyvsp[(2) - (5)].str)); }
    break;

  case 56:
#line 331 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (4)].str),"ieeefloat",4,NULL,(grib_yyvsp[(3) - (4)].explist),(grib_yyvsp[(4) - (4)].lval),NULL,NULL);   free((grib_yyvsp[(2) - (4)].str));  }
    break;

  case 57:
#line 334 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(4) - (6)].str),"ieeefloat",4,NULL,(grib_yyvsp[(5) - (6)].explist),(grib_yyvsp[(6) - (6)].lval),(grib_yyvsp[(2) - (6)].str),NULL);  free((grib_yyvsp[(4) - (6)].str));free((grib_yyvsp[(2) - (6)].str));}
    break;

  case 58:
#line 337 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (7)].str),"ieeefloat",4,(grib_yyvsp[(4) - (7)].explist),(grib_yyvsp[(6) - (7)].explist),(grib_yyvsp[(7) - (7)].lval),NULL,NULL);free((grib_yyvsp[(2) - (7)].str));  }
    break;

  case 59:
#line 340 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (2)].str),"g1_half_byte_codeflag",0,NULL,NULL,0,NULL,NULL);free((grib_yyvsp[(2) - (2)].str));  }
    break;

  case 60:
#line 343 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (5)].str),"section_length",(grib_yyvsp[(3) - (5)].lval),NULL,NULL,0,NULL,NULL);free((grib_yyvsp[(5) - (5)].str));  }
    break;

  case 61:
#line 346 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (8)].str),"g1_message_length",(grib_yyvsp[(3) - (8)].lval),(grib_yyvsp[(7) - (8)].explist),NULL,0,NULL,NULL);free((grib_yyvsp[(5) - (8)].str));  }
    break;

  case 62:
#line 349 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(5) - (8)].str),"g1_section4_length",(grib_yyvsp[(3) - (8)].lval),(grib_yyvsp[(7) - (8)].explist),NULL,0,NULL,NULL);free((grib_yyvsp[(5) - (8)].str));  }
    break;

  case 63:
#line 352 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (3)].str),"ksec",0,(grib_yyvsp[(3) - (3)].explist),NULL,0,NULL,NULL);free((grib_yyvsp[(2) - (3)].str)); }
    break;

  case 64:
#line 355 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (5)].str),"pad",0,(grib_yyvsp[(4) - (5)].explist),0,0,NULL,NULL);   free((grib_yyvsp[(2) - (5)].str)); }
    break;

  case 65:
#line 358 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (5)].str),"padto",0,(grib_yyvsp[(4) - (5)].explist),0,0,NULL,NULL);   free((grib_yyvsp[(2) - (5)].str)); }
    break;

  case 66:
#line 361 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (5)].str),"padtoeven",0,(grib_yyvsp[(4) - (5)].explist),0,0,NULL,NULL);   free((grib_yyvsp[(2) - (5)].str)); }
    break;

  case 67:
#line 364 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (5)].str),"padtomultiple",0,(grib_yyvsp[(4) - (5)].explist),0,0,NULL,NULL);   free((grib_yyvsp[(2) - (5)].str)); }
    break;

  case 68:
#line 367 "griby.y"
    { (grib_yyval.act) = grib_action_create_gen(grib_parser_context,(grib_yyvsp[(2) - (3)].str),"section_padding",0,0,0,(grib_yyvsp[(3) - (3)].lval),NULL,NULL);   free((grib_yyvsp[(2) - (3)].str));  }
    break;

  case 69:
#line 370 "griby.y"
    { (grib_yyval.act) = grib_action_create_template(grib_parser_context,0,(grib_yyvsp[(2) - (3)].str),(grib_yyvsp[(3) - (3)].str)); free((grib_yyvsp[(2) - (3)].str)); free((grib_yyvsp[(3) - (3)].str));}
    break;

  case 70:
#line 372 "griby.y"
    { (grib_yyval.act) = grib_action_create_template(grib_parser_context,1,(grib_yyvsp[(2) - (3)].str),(grib_yyvsp[(3) - (3)].str)); free((grib_yyvsp[(2) - (3)].str)); free((grib_yyvsp[(3) - (3)].str));}
    break;

  case 71:
#line 375 "griby.y"
    { (grib_yyval.act) = grib_action_create_alias(grib_parser_context,(grib_yyvsp[(2) - (5)].str),(grib_yyvsp[(4) - (5)].str),NULL,(grib_yyvsp[(5) - (5)].lval));  free((grib_yyvsp[(2) - (5)].str)); free((grib_yyvsp[(4) - (5)].str)); }
    break;

  case 72:
#line 378 "griby.y"
    { (grib_yyval.act) = grib_action_create_alias(grib_parser_context,(grib_yyvsp[(2) - (2)].str),NULL,NULL,0);  free((grib_yyvsp[(2) - (2)].str)); }
    break;

  case 73:
#line 381 "griby.y"
    {
         (grib_yyval.act) = grib_action_create_alias(grib_parser_context,(grib_yyvsp[(4) - (7)].str),(grib_yyvsp[(6) - (7)].str),(grib_yyvsp[(2) - (7)].str),(grib_yyvsp[(7) - (7)].lval));  free((grib_yyvsp[(2) - (7)].str)); free((grib_yyvsp[(4) - (7)].str)); free((grib_yyvsp[(6) - (7)].str));
    }
    break;

  case 74:
#line 385 "griby.y"
    {
         (grib_yyval.act) = grib_action_create_alias(grib_parser_context,(grib_yyvsp[(4) - (4)].str),NULL,(grib_yyvsp[(2) - (4)].str),0);  free((grib_yyvsp[(2) - (4)].str)); free((grib_yyvsp[(4) - (4)].str)); 
    }
    break;

  case 75:
#line 389 "griby.y"
    { (grib_yyval.act) = grib_action_create_meta(grib_parser_context,(grib_yyvsp[(2) - (8)].str),(grib_yyvsp[(3) - (8)].str),(grib_yyvsp[(5) - (8)].explist),(grib_yyvsp[(7) - (8)].explist),(grib_yyvsp[(8) - (8)].lval),NULL); free((grib_yyvsp[(2) - (8)].str));free((grib_yyvsp[(3) - (8)].str));}
    break;

  case 76:
#line 392 "griby.y"
    { (grib_yyval.act) = grib_action_create_meta(grib_parser_context,(grib_yyvsp[(4) - (10)].str),(grib_yyvsp[(5) - (10)].str),(grib_yyvsp[(7) - (10)].explist),(grib_yyvsp[(9) - (10)].explist),(grib_yyvsp[(10) - (10)].lval),(grib_yyvsp[(2) - (10)].str)); free((grib_yyvsp[(4) - (10)].str));free((grib_yyvsp[(5) - (10)].str));free((grib_yyvsp[(2) - (10)].str));}
    break;

  case 77:
#line 395 "griby.y"
    {
      grib_arguments* a = grib_arguments_new(
        grib_parser_context,
        new_accessor_expression(grib_parser_context,(grib_yyvsp[(2) - (5)].str)),
		NULL
        );
      a->next=(grib_yyvsp[(4) - (5)].explist);
      (grib_yyval.act) = grib_action_create_meta(grib_parser_context,
      "ITERATOR","iterator",a,NULL,
      GRIB_ACCESSOR_FLAG_HIDDEN|GRIB_ACCESSOR_FLAG_READ_ONLY,NULL); free((grib_yyvsp[(2) - (5)].str));
    }
    break;

  case 78:
#line 407 "griby.y"
    {
      grib_arguments* a = grib_arguments_new(
        grib_parser_context,
        new_accessor_expression(grib_parser_context,(grib_yyvsp[(2) - (5)].str)),
		NULL
        );
      a->next=(grib_yyvsp[(4) - (5)].explist);
      (grib_yyval.act) = grib_action_create_meta(grib_parser_context,
      "NEAREST","nearest",a,NULL,
      GRIB_ACCESSOR_FLAG_HIDDEN|GRIB_ACCESSOR_FLAG_READ_ONLY,NULL); free((grib_yyvsp[(2) - (5)].str));
    }
    break;

  case 79:
#line 419 "griby.y"
    {
      grib_arguments* a = grib_arguments_new(
        grib_parser_context,
        new_accessor_expression(grib_parser_context,(grib_yyvsp[(2) - (5)].str)),
		NULL
        );
      a->next=(grib_yyvsp[(4) - (5)].explist);
      (grib_yyval.act) = grib_action_create_meta(grib_parser_context,
      "BOX","box",a,NULL,
      GRIB_ACCESSOR_FLAG_HIDDEN|GRIB_ACCESSOR_FLAG_READ_ONLY,NULL); free((grib_yyvsp[(2) - (5)].str));
    }
    break;

  case 80:
#line 431 "griby.y"
    { (grib_yyval.act) = grib_action_create_put(grib_parser_context,(grib_yyvsp[(2) - (5)].str),(grib_yyvsp[(4) - (5)].explist));free((grib_yyvsp[(2) - (5)].str));}
    break;

  case 81:
#line 434 "griby.y"
    { (grib_yyval.act) = grib_action_create_remove(grib_parser_context,(grib_yyvsp[(2) - (2)].explist));}
    break;

  case 82:
#line 437 "griby.y"
    { (grib_yyval.act) = grib_action_create_assert(grib_parser_context,(grib_yyvsp[(3) - (4)].exp));}
    break;

  case 83:
#line 440 "griby.y"
    { (grib_yyval.act) = grib_action_create_modify(grib_parser_context,(grib_yyvsp[(2) - (3)].str),(grib_yyvsp[(3) - (3)].lval)); free((grib_yyvsp[(2) - (3)].str));}
    break;

  case 84:
#line 442 "griby.y"
    { (grib_yyval.act) = grib_action_create_set_missing(grib_parser_context,(grib_yyvsp[(2) - (4)].str)); free((grib_yyvsp[(2) - (4)].str)); }
    break;

  case 85:
#line 443 "griby.y"
    { (grib_yyval.act) = grib_action_create_set(grib_parser_context,(grib_yyvsp[(2) - (4)].str),(grib_yyvsp[(4) - (4)].exp),0); free((grib_yyvsp[(2) - (4)].str)); }
    break;

  case 86:
#line 444 "griby.y"
    { (grib_yyval.act) = grib_action_create_set_darray(grib_parser_context,(grib_yyvsp[(2) - (6)].str),(grib_yyvsp[(5) - (6)].dvalue)); free((grib_yyvsp[(2) - (6)].str)); }
    break;

  case 87:
#line 446 "griby.y"
    { (grib_yyval.act) = grib_action_create_set(grib_parser_context,(grib_yyvsp[(2) - (4)].str),(grib_yyvsp[(4) - (4)].exp),1); free((grib_yyvsp[(2) - (4)].str)); }
    break;

  case 88:
#line 449 "griby.y"
    { (grib_yyval.act) = grib_action_create_write(grib_parser_context,(grib_yyvsp[(2) - (2)].str),0); free((grib_yyvsp[(2) - (2)].str));}
    break;

  case 89:
#line 450 "griby.y"
    { (grib_yyval.act) = grib_action_create_write(grib_parser_context,"",0); }
    break;

  case 90:
#line 451 "griby.y"
    { (grib_yyval.act) = grib_action_create_write(grib_parser_context,(grib_yyvsp[(2) - (2)].str),1); free((grib_yyvsp[(2) - (2)].str));}
    break;

  case 91:
#line 452 "griby.y"
    { (grib_yyval.act) = grib_action_create_write(grib_parser_context,"",1); }
    break;

  case 92:
#line 454 "griby.y"
    { (grib_yyval.act) = grib_action_create_print(grib_parser_context,(grib_yyvsp[(2) - (2)].str),0); free((grib_yyvsp[(2) - (2)].str)); }
    break;

  case 93:
#line 455 "griby.y"
    { (grib_yyval.act) = grib_action_create_print(grib_parser_context,(grib_yyvsp[(5) - (5)].str),(grib_yyvsp[(3) - (5)].str)); free((grib_yyvsp[(5) - (5)].str)); free((grib_yyvsp[(3) - (5)].str));}
    break;

  case 94:
#line 456 "griby.y"
    { (grib_yyval.act) = grib_action_create_print(grib_parser_context,"",0);  }
    break;

  case 95:
#line 460 "griby.y"
    { (grib_yyval.act) = grib_action_create_if(grib_parser_context,(grib_yyvsp[(3) - (7)].exp),(grib_yyvsp[(6) - (7)].act),0,0); }
    break;

  case 96:
#line 461 "griby.y"
    { (grib_yyval.act) = grib_action_create_if(grib_parser_context,(grib_yyvsp[(3) - (11)].exp),(grib_yyvsp[(6) - (11)].act),(grib_yyvsp[(10) - (11)].act),0); }
    break;

  case 97:
#line 462 "griby.y"
    { (grib_yyval.act) = grib_action_create_if(grib_parser_context,(grib_yyvsp[(3) - (7)].exp),(grib_yyvsp[(6) - (7)].act),0,1); }
    break;

  case 98:
#line 463 "griby.y"
    { (grib_yyval.act) = grib_action_create_if(grib_parser_context,(grib_yyvsp[(3) - (11)].exp),(grib_yyvsp[(6) - (11)].act),(grib_yyvsp[(10) - (11)].act),1); }
    break;

  case 99:
#line 467 "griby.y"
    { (grib_yyval.act) = grib_action_create_when(grib_parser_context,(grib_yyvsp[(3) - (6)].exp),(grib_yyvsp[(5) - (6)].act),NULL); }
    break;

  case 100:
#line 468 "griby.y"
    { (grib_yyval.act) = grib_action_create_when(grib_parser_context,(grib_yyvsp[(3) - (7)].exp),(grib_yyvsp[(6) - (7)].act),NULL); }
    break;

  case 101:
#line 469 "griby.y"
    { (grib_yyval.act) = grib_action_create_when(grib_parser_context,(grib_yyvsp[(3) - (11)].exp),(grib_yyvsp[(6) - (11)].act),(grib_yyvsp[(10) - (11)].act)); }
    break;

  case 102:
#line 472 "griby.y"
    { (grib_yyval.act) = grib_action_create_set(grib_parser_context,(grib_yyvsp[(2) - (4)].str),(grib_yyvsp[(4) - (4)].exp),0); free((grib_yyvsp[(2) - (4)].str)); }
    break;

  case 103:
#line 473 "griby.y"
    { (grib_yyval.act) = grib_action_create_set(grib_parser_context,(grib_yyvsp[(2) - (4)].str),(grib_yyvsp[(4) - (4)].exp),1); free((grib_yyvsp[(2) - (4)].str)); }
    break;

  case 105:
#line 477 "griby.y"
    { (grib_yyvsp[(1) - (3)].act)->next = (grib_yyvsp[(2) - (3)].act); (grib_yyval.act) = (grib_yyvsp[(1) - (3)].act); }
    break;

  case 106:
#line 481 "griby.y"
    { (grib_yyval.explist) = NULL ;}
    break;

  case 107:
#line 482 "griby.y"
    { (grib_yyval.explist) = (grib_yyvsp[(2) - (2)].explist) ;}
    break;

  case 108:
#line 485 "griby.y"
    { (grib_yyval.lval) = 0 ; }
    break;

  case 109:
#line 486 "griby.y"
    { (grib_yyval.lval) = (grib_yyvsp[(2) - (2)].lval); }
    break;

  case 111:
#line 490 "griby.y"
    { (grib_yyval.lval) = (grib_yyvsp[(1) - (3)].lval) | (grib_yyvsp[(3) - (3)].lval); }
    break;

  case 112:
#line 493 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_READ_ONLY; }
    break;

  case 113:
#line 494 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_LOWERCASE; }
    break;

  case 114:
#line 495 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_DUMP; }
    break;

  case 115:
#line 496 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_NO_COPY; }
    break;

  case 116:
#line 497 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_NO_FAIL; }
    break;

  case 117:
#line 498 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_HIDDEN; }
    break;

  case 118:
#line 499 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_EDITION_SPECIFIC; }
    break;

  case 119:
#line 500 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_CAN_BE_MISSING; }
    break;

  case 120:
#line 501 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_CONSTRAINT; }
    break;

  case 121:
#line 502 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_OVERRIDE; }
    break;

  case 122:
#line 503 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_COPY_OK; }
    break;

  case 123:
#line 504 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_TRANSIENT; }
    break;

  case 124:
#line 505 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_STRING_TYPE; }
    break;

  case 125:
#line 506 "griby.y"
    { (grib_yyval.lval) = GRIB_ACCESSOR_FLAG_LONG_TYPE; }
    break;

  case 126:
#line 509 "griby.y"
    { (grib_yyval.act) = grib_action_create_list(grib_parser_context,(grib_yyvsp[(1) - (8)].str),(grib_yyvsp[(4) - (8)].exp),(grib_yyvsp[(7) - (8)].act)); free((grib_yyvsp[(1) - (8)].str)); }
    break;

  case 127:
#line 512 "griby.y"
    { (grib_yyval.act) = grib_action_create_while(grib_parser_context,(grib_yyvsp[(3) - (7)].exp),(grib_yyvsp[(6) - (7)].act));  }
    break;

  case 128:
#line 515 "griby.y"
    { (grib_yyval.act) = grib_action_create_trigger(grib_parser_context,(grib_yyvsp[(3) - (7)].explist),(grib_yyvsp[(6) - (7)].act));  }
    break;

  case 129:
#line 518 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(2) - (6)].str),(grib_yyvsp[(4) - (6)].concept_value),0,0,0,0,0,0,(grib_yyvsp[(6) - (6)].lval),0);  free((grib_yyvsp[(2) - (6)].str)); }
    break;

  case 130:
#line 519 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(2) - (9)].str),(grib_yyvsp[(7) - (9)].concept_value),0,0,(grib_yyvsp[(4) - (9)].str),0,0,0,(grib_yyvsp[(9) - (9)].lval),0);  free((grib_yyvsp[(2) - (9)].str));free((grib_yyvsp[(4) - (9)].str)); }
    break;

  case 131:
#line 520 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(2) - (12)].str),0,(grib_yyvsp[(6) - (12)].str),0,(grib_yyvsp[(4) - (12)].str),(grib_yyvsp[(8) - (12)].str),(grib_yyvsp[(10) - (12)].str),0,(grib_yyvsp[(12) - (12)].lval),0);  free((grib_yyvsp[(2) - (12)].str));free((grib_yyvsp[(6) - (12)].str));free((grib_yyvsp[(4) - (12)].str));free((grib_yyvsp[(8) - (12)].str));free((grib_yyvsp[(10) - (12)].str)); }
    break;

  case 132:
#line 521 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(2) - (14)].str),0,(grib_yyvsp[(6) - (14)].str),0,(grib_yyvsp[(4) - (14)].str),(grib_yyvsp[(8) - (14)].str),(grib_yyvsp[(10) - (14)].str),(grib_yyvsp[(12) - (14)].str),(grib_yyvsp[(14) - (14)].lval),0);  free((grib_yyvsp[(2) - (14)].str));free((grib_yyvsp[(6) - (14)].str));free((grib_yyvsp[(4) - (14)].str));free((grib_yyvsp[(8) - (14)].str));free((grib_yyvsp[(10) - (14)].str));free((grib_yyvsp[(12) - (14)].str)); }
    break;

  case 133:
#line 522 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(2) - (10)].str),0,(grib_yyvsp[(6) - (10)].str),0,(grib_yyvsp[(4) - (10)].str),(grib_yyvsp[(8) - (10)].str),0,0,(grib_yyvsp[(10) - (10)].lval),0);  free((grib_yyvsp[(2) - (10)].str));free((grib_yyvsp[(6) - (10)].str));free((grib_yyvsp[(4) - (10)].str));free((grib_yyvsp[(8) - (10)].str)); }
    break;

  case 134:
#line 523 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(4) - (14)].str),0,(grib_yyvsp[(8) - (14)].str),(grib_yyvsp[(2) - (14)].str),(grib_yyvsp[(6) - (14)].str),(grib_yyvsp[(10) - (14)].str),(grib_yyvsp[(12) - (14)].str),0,(grib_yyvsp[(14) - (14)].lval),0);  free((grib_yyvsp[(4) - (14)].str));free((grib_yyvsp[(8) - (14)].str));free((grib_yyvsp[(6) - (14)].str));free((grib_yyvsp[(10) - (14)].str)); free((grib_yyvsp[(2) - (14)].str));}
    break;

  case 135:
#line 524 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(4) - (12)].str),0,(grib_yyvsp[(8) - (12)].str),(grib_yyvsp[(2) - (12)].str),(grib_yyvsp[(6) - (12)].str),(grib_yyvsp[(10) - (12)].str),0,0,(grib_yyvsp[(12) - (12)].lval),0);  free((grib_yyvsp[(4) - (12)].str));free((grib_yyvsp[(8) - (12)].str));free((grib_yyvsp[(6) - (12)].str));free((grib_yyvsp[(10) - (12)].str)); free((grib_yyvsp[(2) - (12)].str));}
    break;

  case 136:
#line 525 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(4) - (8)].str),(grib_yyvsp[(6) - (8)].concept_value),0,(grib_yyvsp[(2) - (8)].str),0,0,0,0,(grib_yyvsp[(8) - (8)].lval),0);  free((grib_yyvsp[(2) - (8)].str));free((grib_yyvsp[(4) - (8)].str)); }
    break;

  case 137:
#line 526 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(4) - (11)].str),(grib_yyvsp[(9) - (11)].concept_value),0,(grib_yyvsp[(2) - (11)].str),(grib_yyvsp[(6) - (11)].str),0,0,0,(grib_yyvsp[(11) - (11)].lval),0);  free((grib_yyvsp[(2) - (11)].str));free((grib_yyvsp[(4) - (11)].str));free((grib_yyvsp[(6) - (11)].str)); }
    break;

  case 138:
#line 527 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(2) - (6)].str),(grib_yyvsp[(4) - (6)].concept_value),0,0,0,0,0,0,(grib_yyvsp[(6) - (6)].lval),1);  free((grib_yyvsp[(2) - (6)].str)); }
    break;

  case 139:
#line 528 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(2) - (9)].str),(grib_yyvsp[(7) - (9)].concept_value),0,0,(grib_yyvsp[(4) - (9)].str),0,0,0,(grib_yyvsp[(9) - (9)].lval),1);  free((grib_yyvsp[(2) - (9)].str));free((grib_yyvsp[(4) - (9)].str)); }
    break;

  case 140:
#line 529 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(2) - (12)].str),0,(grib_yyvsp[(6) - (12)].str),0,(grib_yyvsp[(4) - (12)].str),(grib_yyvsp[(8) - (12)].str),(grib_yyvsp[(10) - (12)].str),0,(grib_yyvsp[(12) - (12)].lval),1);  free((grib_yyvsp[(2) - (12)].str));free((grib_yyvsp[(6) - (12)].str));free((grib_yyvsp[(4) - (12)].str));free((grib_yyvsp[(8) - (12)].str));free((grib_yyvsp[(10) - (12)].str)); }
    break;

  case 141:
#line 530 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(2) - (10)].str),0,(grib_yyvsp[(6) - (10)].str),0,(grib_yyvsp[(4) - (10)].str),(grib_yyvsp[(8) - (10)].str),0,0,(grib_yyvsp[(10) - (10)].lval),1);  free((grib_yyvsp[(2) - (10)].str));free((grib_yyvsp[(6) - (10)].str));free((grib_yyvsp[(4) - (10)].str));free((grib_yyvsp[(8) - (10)].str)); }
    break;

  case 142:
#line 531 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(4) - (14)].str),0,(grib_yyvsp[(8) - (14)].str),(grib_yyvsp[(2) - (14)].str),(grib_yyvsp[(6) - (14)].str),(grib_yyvsp[(10) - (14)].str),(grib_yyvsp[(12) - (14)].str),0,(grib_yyvsp[(14) - (14)].lval),1);  free((grib_yyvsp[(4) - (14)].str));free((grib_yyvsp[(8) - (14)].str));free((grib_yyvsp[(6) - (14)].str));free((grib_yyvsp[(10) - (14)].str));free((grib_yyvsp[(12) - (14)].str)); free((grib_yyvsp[(2) - (14)].str));}
    break;

  case 143:
#line 532 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(4) - (12)].str),0,(grib_yyvsp[(8) - (12)].str),(grib_yyvsp[(2) - (12)].str),(grib_yyvsp[(6) - (12)].str),(grib_yyvsp[(10) - (12)].str),0,0,(grib_yyvsp[(12) - (12)].lval),1);  free((grib_yyvsp[(4) - (12)].str));free((grib_yyvsp[(8) - (12)].str));free((grib_yyvsp[(6) - (12)].str));free((grib_yyvsp[(10) - (12)].str)); free((grib_yyvsp[(2) - (12)].str));}
    break;

  case 144:
#line 533 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(4) - (8)].str),(grib_yyvsp[(6) - (8)].concept_value),0,(grib_yyvsp[(2) - (8)].str),0,0,0,0,(grib_yyvsp[(8) - (8)].lval),1);  free((grib_yyvsp[(2) - (8)].str));free((grib_yyvsp[(4) - (8)].str)); }
    break;

  case 145:
#line 534 "griby.y"
    { (grib_yyval.act) = grib_action_create_concept(grib_parser_context,(grib_yyvsp[(4) - (11)].str),(grib_yyvsp[(9) - (11)].concept_value),0,(grib_yyvsp[(2) - (11)].str),(grib_yyvsp[(6) - (11)].str),0,0,0,(grib_yyvsp[(11) - (11)].lval),1);  free((grib_yyvsp[(2) - (11)].str));free((grib_yyvsp[(4) - (11)].str));free((grib_yyvsp[(6) - (11)].str)); }
    break;

  case 147:
#line 539 "griby.y"
    { (grib_yyval.concept_value) = (grib_yyvsp[(2) - (2)].concept_value); (grib_yyvsp[(2) - (2)].concept_value)->next = (grib_yyvsp[(1) - (2)].concept_value);   }
    break;

  case 149:
#line 543 "griby.y"
    { (grib_yyval.case_value) = (grib_yyvsp[(2) - (2)].case_value); (grib_yyvsp[(2) - (2)].case_value)->next = (grib_yyvsp[(1) - (2)].case_value);   }
    break;

  case 150:
#line 546 "griby.y"
    { (grib_yyval.case_value) = grib_case_new(grib_parser_context,(grib_yyvsp[(2) - (4)].explist),(grib_yyvsp[(4) - (4)].act));  }
    break;

  case 151:
#line 550 "griby.y"
    { (grib_yyval.act) = grib_action_create_switch(grib_parser_context,(grib_yyvsp[(3) - (10)].explist),(grib_yyvsp[(6) - (10)].case_value),(grib_yyvsp[(9) - (10)].act)); }
    break;

  case 152:
#line 551 "griby.y"
    { (grib_yyval.act) = grib_action_create_switch(grib_parser_context,(grib_yyvsp[(3) - (9)].explist),(grib_yyvsp[(6) - (9)].case_value),grib_action_create_noop(grib_parser_context,"continue")); }
    break;

  case 153:
#line 552 "griby.y"
    { (grib_yyval.act) = grib_action_create_switch(grib_parser_context,(grib_yyvsp[(3) - (7)].explist),(grib_yyvsp[(6) - (7)].case_value),0); }
    break;

  case 154:
#line 555 "griby.y"
    {
	  				(grib_yyval.concept_value) = grib_concept_value_new(grib_parser_context,(grib_yyvsp[(1) - (5)].str),(grib_yyvsp[(4) - (5)].concept_condition)); free((grib_yyvsp[(1) - (5)].str));}
    break;

  case 155:
#line 557 "griby.y"
    {
	  				(grib_yyval.concept_value) = grib_concept_value_new(grib_parser_context,(grib_yyvsp[(1) - (5)].str),(grib_yyvsp[(4) - (5)].concept_condition)); free((grib_yyvsp[(1) - (5)].str));}
    break;

  case 156:
#line 559 "griby.y"
    {
					char buf[80]; sprintf(buf,"%ld",(long)(grib_yyvsp[(1) - (5)].lval)); (grib_yyval.concept_value) = grib_concept_value_new(grib_parser_context,buf,(grib_yyvsp[(4) - (5)].concept_condition));}
    break;

  case 157:
#line 561 "griby.y"
    {
					char buf[80]; sprintf(buf,"%g",(double)(grib_yyvsp[(1) - (5)].dval)); (grib_yyval.concept_value) = grib_concept_value_new(grib_parser_context,buf,(grib_yyvsp[(4) - (5)].concept_condition));}
    break;

  case 159:
#line 566 "griby.y"
    { (grib_yyvsp[(1) - (2)].concept_condition)->next = (grib_yyvsp[(2) - (2)].concept_condition); (grib_yyval.concept_condition) = (grib_yyvsp[(1) - (2)].concept_condition); }
    break;

  case 160:
#line 569 "griby.y"
    { (grib_yyval.concept_condition) = grib_concept_condition_new(grib_parser_context,(grib_yyvsp[(1) - (4)].str),(grib_yyvsp[(3) - (4)].exp)); free((grib_yyvsp[(1) - (4)].str)); }
    break;

  case 161:
#line 572 "griby.y"
    { (grib_yyval.exp) = new_accessor_expression(grib_parser_context,(grib_yyvsp[(1) - (1)].str)); free((grib_yyvsp[(1) - (1)].str)); }
    break;

  case 162:
#line 573 "griby.y"
    { (grib_yyval.exp) = new_string_expression(grib_parser_context,(grib_yyvsp[(1) - (1)].str));  free((grib_yyvsp[(1) - (1)].str)); }
    break;

  case 164:
#line 577 "griby.y"
    { (grib_yyval.exp) = new_long_expression(grib_parser_context,(grib_yyvsp[(1) - (1)].lval));  }
    break;

  case 165:
#line 578 "griby.y"
    { (grib_yyval.exp) = new_double_expression(grib_parser_context,(grib_yyvsp[(1) - (1)].dval));  /* TODO: change to new_float_expression*/}
    break;

  case 166:
#line 580 "griby.y"
    { (grib_yyval.exp) = NULL; }
    break;

  case 167:
#line 581 "griby.y"
    { (grib_yyval.exp) = new_true_expression(grib_parser_context); }
    break;

  case 168:
#line 582 "griby.y"
    { (grib_yyval.exp) = (grib_yyvsp[(2) - (3)].exp); }
    break;

  case 169:
#line 583 "griby.y"
    { (grib_yyval.exp) = new_unop_expression(grib_parser_context,&grib_op_neg,&grib_op_neg_d,(grib_yyvsp[(2) - (2)].exp)); }
    break;

  case 170:
#line 584 "griby.y"
    { (grib_yyval.exp) = new_func_expression(grib_parser_context,(grib_yyvsp[(1) - (3)].str),NULL); free((grib_yyvsp[(1) - (3)].str));}
    break;

  case 171:
#line 585 "griby.y"
    { (grib_yyval.exp) = new_func_expression(grib_parser_context,(grib_yyvsp[(1) - (4)].str),(grib_yyvsp[(3) - (4)].explist)); free((grib_yyvsp[(1) - (4)].str));}
    break;

  case 172:
#line 589 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_pow,NULL,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 174:
#line 593 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_mul,&grib_op_mul_d,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 175:
#line 594 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_div,&grib_op_div_d,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 176:
#line 595 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_modulo,NULL,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 177:
#line 596 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_bit,NULL,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 178:
#line 597 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_bitoff,NULL,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 180:
#line 601 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_add,&grib_op_add_d,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 181:
#line 602 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_sub,&grib_op_sub_d,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 183:
#line 606 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_gt,&grib_op_gt_d,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 184:
#line 608 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_eq,&grib_op_eq_d,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 185:
#line 609 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_lt,&grib_op_lt_d,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 186:
#line 610 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_ge,&grib_op_ge_d,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 187:
#line 611 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_le,&grib_op_le_d,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 188:
#line 612 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_ne,&grib_op_ne_d,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 189:
#line 613 "griby.y"
    { (grib_yyval.exp) = new_string_compare_expression(grib_parser_context,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 190:
#line 618 "griby.y"
    { (grib_yyval.exp) = new_unop_expression(grib_parser_context,&grib_op_not,NULL,(grib_yyvsp[(2) - (2)].exp)); }
    break;

  case 192:
#line 622 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_and,NULL,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp)); }
    break;

  case 194:
#line 626 "griby.y"
    { (grib_yyval.exp) = new_binop_expression(grib_parser_context,&grib_op_or,NULL,(grib_yyvsp[(1) - (3)].exp),(grib_yyvsp[(3) - (3)].exp));}
    break;

  case 199:
#line 640 "griby.y"
    { (grib_yyval.rule_entry) = grib_new_rule_entry(grib_parser_context,(grib_yyvsp[(1) - (4)].str),(grib_yyvsp[(3) - (4)].exp)); free((grib_yyvsp[(1) - (4)].str)); }
    break;

  case 200:
#line 641 "griby.y"
    { (grib_yyval.rule_entry) = grib_new_rule_entry(grib_parser_context,"skip",0);}
    break;

  case 202:
#line 645 "griby.y"
    { (grib_yyvsp[(1) - (2)].rule_entry)->next = (grib_yyvsp[(2) - (2)].rule_entry); (grib_yyval.rule_entry) = (grib_yyvsp[(1) - (2)].rule_entry); }
    break;

  case 203:
#line 648 "griby.y"
    { (grib_yyval.rules) = grib_new_rule(grib_parser_context,NULL,(grib_yyvsp[(1) - (1)].rule_entry)); }
    break;

  case 204:
#line 652 "griby.y"
    { (grib_yyval.rules) = grib_new_rule(grib_parser_context,(grib_yyvsp[(3) - (7)].exp),(grib_yyvsp[(6) - (7)].rule_entry)); }
    break;

  case 206:
#line 656 "griby.y"
    { (grib_yyvsp[(1) - (2)].rules)->next = (grib_yyvsp[(2) - (2)].rules); (grib_yyval.rules) = (grib_yyvsp[(1) - (2)].rules); }
    break;


/* Line 1267 of yacc.c.  */
#line 3086 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", grib_yyr1[grib_yyn], &grib_yyval, &grib_yyloc);

  YYPOPSTACK (grib_yylen);
  grib_yylen = 0;
  YY_STACK_PRINT (grib_yyss, grib_yyssp);

  *++grib_yyvsp = grib_yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  grib_yyn = grib_yyr1[grib_yyn];

  grib_yystate = grib_yypgoto[grib_yyn - YYNTOKENS] + *grib_yyssp;
  if (0 <= grib_yystate && grib_yystate <= YYLAST && grib_yycheck[grib_yystate] == *grib_yyssp)
    grib_yystate = grib_yytable[grib_yystate];
  else
    grib_yystate = grib_yydefgoto[grib_yyn - YYNTOKENS];

  goto grib_yynewstate;


/*------------------------------------.
| grib_yyerrlab -- here on detecting error |
`------------------------------------*/
grib_yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!grib_yyerrstatus)
    {
      ++grib_yynerrs;
#if ! YYERROR_VERBOSE
      grib_yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T grib_yysize = grib_yysyntax_error (0, grib_yystate, grib_yychar);
	if (grib_yymsg_alloc < grib_yysize && grib_yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T grib_yyalloc = 2 * grib_yysize;
	    if (! (grib_yysize <= grib_yyalloc && grib_yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      grib_yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (grib_yymsg != grib_yymsgbuf)
	      YYSTACK_FREE (grib_yymsg);
	    grib_yymsg = (char *) YYSTACK_ALLOC (grib_yyalloc);
	    if (grib_yymsg)
	      grib_yymsg_alloc = grib_yyalloc;
	    else
	      {
		grib_yymsg = grib_yymsgbuf;
		grib_yymsg_alloc = sizeof grib_yymsgbuf;
	      }
	  }

	if (0 < grib_yysize && grib_yysize <= grib_yymsg_alloc)
	  {
	    (void) grib_yysyntax_error (grib_yymsg, grib_yystate, grib_yychar);
	    grib_yyerror (grib_yymsg);
	  }
	else
	  {
	    grib_yyerror (YY_("syntax error"));
	    if (grib_yysize != 0)
	      goto grib_yyexhaustedlab;
	  }
      }
#endif
    }



  if (grib_yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (grib_yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (grib_yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  grib_yydestruct ("Error: discarding",
		      grib_yytoken, &grib_yylval);
	  grib_yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto grib_yyerrlab1;


/*---------------------------------------------------.
| grib_yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
grib_yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label grib_yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto grib_yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (grib_yylen);
  grib_yylen = 0;
  YY_STACK_PRINT (grib_yyss, grib_yyssp);
  grib_yystate = *grib_yyssp;
  goto grib_yyerrlab1;


/*-------------------------------------------------------------.
| grib_yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
grib_yyerrlab1:
  grib_yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      grib_yyn = grib_yypact[grib_yystate];
      if (grib_yyn != YYPACT_NINF)
	{
	  grib_yyn += YYTERROR;
	  if (0 <= grib_yyn && grib_yyn <= YYLAST && grib_yycheck[grib_yyn] == YYTERROR)
	    {
	      grib_yyn = grib_yytable[grib_yyn];
	      if (0 < grib_yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (grib_yyssp == grib_yyss)
	YYABORT;


      grib_yydestruct ("Error: popping",
		  grib_yystos[grib_yystate], grib_yyvsp);
      YYPOPSTACK (1);
      grib_yystate = *grib_yyssp;
      YY_STACK_PRINT (grib_yyss, grib_yyssp);
    }

  if (grib_yyn == YYFINAL)
    YYACCEPT;

  *++grib_yyvsp = grib_yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", grib_yystos[grib_yyn], grib_yyvsp, grib_yylsp);

  grib_yystate = grib_yyn;
  goto grib_yynewstate;


/*-------------------------------------.
| grib_yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
grib_yyacceptlab:
  grib_yyresult = 0;
  goto grib_yyreturn;

/*-----------------------------------.
| grib_yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
grib_yyabortlab:
  grib_yyresult = 1;
  goto grib_yyreturn;

#ifndef grib_yyoverflow
/*-------------------------------------------------.
| grib_yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
grib_yyexhaustedlab:
  grib_yyerror (YY_("memory exhausted"));
  grib_yyresult = 2;
  /* Fall through.  */
#endif

grib_yyreturn:
  if (grib_yychar != YYEOF && grib_yychar != YYEMPTY)
     grib_yydestruct ("Cleanup: discarding lookahead",
		 grib_yytoken, &grib_yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (grib_yylen);
  YY_STACK_PRINT (grib_yyss, grib_yyssp);
  while (grib_yyssp != grib_yyss)
    {
      grib_yydestruct ("Cleanup: popping",
		  grib_yystos[*grib_yyssp], grib_yyvsp);
      YYPOPSTACK (1);
    }
#ifndef grib_yyoverflow
  if (grib_yyss != grib_yyssa)
    YYSTACK_FREE (grib_yyss);
#endif
#if YYERROR_VERBOSE
  if (grib_yymsg != grib_yymsgbuf)
    YYSTACK_FREE (grib_yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (grib_yyresult);
}


#line 660 "griby.y"


static grib_concept_value *reverse_concept(grib_concept_value *r,grib_concept_value *s)
{
    grib_concept_value *v;

    if(r == NULL) return s;

    v         = r->next;
    r->next   = s;
    return reverse_concept(v,r);
}


static grib_concept_value* reverse(grib_concept_value* r)
{
    return reverse_concept(r,NULL);
}




