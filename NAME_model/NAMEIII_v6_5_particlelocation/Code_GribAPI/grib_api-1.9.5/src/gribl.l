/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/


%{
#include "grib_api_internal.h"
#include "grib_yacc.h"
#include <ctype.h>

int yylineno;

/* Keep -Wall quiet */


/*

This is needed for implementing "include", otherwise
flex buffer optimization break the includes.

*/

#define YY_INPUT(buf,result,max_size) \
        { \
        int c = fgetc(yyin); \
        result = (c == EOF) ? YY_NULL : (buf[0] = c, 1); \
        }


/*
#ifdef __GNUC__
#ifdef FLEX_SCANNER
static void yyunput (int c,char *buf_ptr  );
void _grib_ignore_yyunput_unused_error() { yyunput(0,0); }
#endif
#endif
*/

%}


SIGN      [\-\+]
E         [eE]
DIGIT     [0-9]
SIGNED    {SIGN}?{DIGIT}+
NUMB      {DIGIT}+
EXP       {E}{SIGNED}
IDENT     [_A-Za-z]+[_0-9A-Za-z]*
FLOAT1    {SIGNED}+"."{DIGIT}*{EXP}?
FLOAT2    {SIGNED}*"."{DIGIT}+{EXP}?
FLOAT3    {SIGNED}+{EXP}?

%%


"=="           return EQ       ;
">="           return GE       ;
"<="           return LE       ;
"!="           return NE       ;
"<>"           return NE       ;
"bit"          return BIT      ;
"notbit"       return BITOFF   ;

"is"       return IS   ;
"not"       return NOT   ;
"!"         return NOT   ;
"and"       return AND   ;
"&&"        return AND   ;
"or"        return OR   ;
"||"        return OR   ;

"null"        return NIL   ;
"~"           return DUMMY   ;

"lowercase" return LOWERCASE;
"if"        return IF       ;
"_if"        return IF_TRANSIENT       ;
"else"      return ELSE       ;
"unsigned"  return UNSIGNED ;
"ascii"     return ASCII    ;
"byte"     return BYTE    ;
"label"     return LABEL    ;
"list"      return LIST     ;
"while"      return WHILE     ;
"template"  return TEMPLATE ;
"template_nofail"  return TEMPLATE_NOFAIL ;
"trigger"  return TRIGGER ;
"end"       return END      ;
"ibmfloat"  return IBMFLOAT ;
"ieeefloat" return FLOAT ;
"signed"    return SIGNED   ;
"codetable" return CODETABLE;
"complex_codetable" return CODETABLE;
"flags"     return FLAG     ;
"lookup"    return LOOKUP   ;
"meta"      return META     ;
"padtoeven"     return PADTOEVEN    ;
"padto"     return PADTO    ;
"padtomultiple"     return PADTOMULTIPLE    ;
"pad"       return PAD      ;
"section_padding"       return SECTION_PADDING      ;
"alias"     return ALIAS    ;
"unalias"     return UNALIAS    ;
"position"  return POS      ;
"constant"  return INTCONST ;
"transient" return TRANS    ;
"string_type" return STRING_TYPE   ;
"long_type" return LONG_TYPE   ;
"iterator"  return ITERATOR ;
"nearest"  return NEAREST ;
"box"  return BOX ;
"ksec"      return KSEC    ;
"flagbit"      return FLAGBIT    ;
"ksec1expver" return KSEC1EXPVER    ;
"modify"      return MODIFY    ;

"g1_half_byte_codeflag" return G1_HALF_BYTE    ;
"g1_message_length" return G1_MESSAGE_LENGTH    ;
"g1_section4_length" return G1_SECTION4_LENGTH    ;

"export"     return EXPORT;
"remove"     return REMOVE;

"length"          return SECTION_LENGTH     ;

"assert"     return ASSERT    ;

"read_only"           return READ_ONLY;
"no_copy"           return NO_COPY;
"edition_specific"    return EDITION_SPECIFIC;
"dump"                return DUMP;
"no_fail"            return NO_FAIL;
"hidden"              return HIDDEN;
"can_be_missing"      return CAN_BE_MISSING;
"MISSING"      return MISSING;
"constraint"      return CONSTRAINT;
"override"      return OVERRIDE;
"copy_ok"      return COPY_OK;

"set"         return SET;
"set_nofail"         return SET_NOFAIL;
"when"        return WHEN;
"case"     return CASE;
"switch"     return SWITCH;
"default"     return DEFAULT;
"concept"     return CONCEPT;
"concept_nofail"     return CONCEPT_NOFAIL;
"write"       return WRITE;
"append"       return APPEND;
"print"       return PRINT;
"skip"         return SKIP;

"include"       {
          int c,q;
          while((c = input()) && isspace(c) && c != '\n') ;
          q = c; /* the quote */

          yyleng = 0;
          while((c = input()) && c != q && c != '\n')
                {
                    if(c == '\\') yytext[yyleng++] = input();
                    else yytext[yyleng++] =  c;
                }

                yytext[yyleng++] = 0;

          grib_parser_include(yytext);
         }

\"|\'  {
           int c,q = yytext[0];

           yyleng = 0;

           while((c = input()) && c != q && c != '\n')
           {
               if(c == '\\') yytext[yyleng++] = input();
               else yytext[yyleng++] =  c;
            }

            yytext[yyleng++] = 0;
            yylval.str = strdup(yytext);
            return STRING;
        }

`       {
           int c;
           unsigned long val = 0;

           while((c = input()) && c != '`' && c != '\n')
           {
                val <<= 8;
                val |= c;
           }
           yylval.lval = val;
           return INTEGER;
        }


{IDENT}       { yylval.str = strdup(yytext); return IDENT; }
{NUMB}        { yylval.lval = atol((const char *)yytext); return INTEGER; }
{FLOAT1}      { yylval.dval = atof((const char *)yytext); return FLOAT; }
{FLOAT2}      { yylval.dval = atof((const char *)yytext); return FLOAT; }
{FLOAT3}      { yylval.dval = atof((const char *)yytext); return FLOAT; }


\#      {
           int c;
           while((c = input()) && (c != '\n')){}

        yylineno++;
        }
[ \t\r]*  ;
\n     yylineno++;



.       return *yytext;

%%

