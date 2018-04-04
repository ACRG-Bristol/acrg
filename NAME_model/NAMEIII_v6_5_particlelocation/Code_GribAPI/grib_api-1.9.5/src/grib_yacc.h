/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

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
/* Line 1489 of yacc.c.  */
#line 241 "y.tab.h"
	YYSTYPE;
# define grib_yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE grib_yylval;

