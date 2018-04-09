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
/* #include "grib_parser.h" */

extern int yylex();
extern int yyerror(const char*);

extern   grib_action*           grib_parser_all_actions;
extern   grib_concept_value*    grib_parser_concept;
extern   grib_context*          grib_parser_context;
extern   grib_rule*             grib_parser_rules;

static grib_concept_value* reverse(grib_concept_value* r);
static grib_concept_value *reverse_concept(grib_concept_value *r,grib_concept_value *s);

/* typedef int (*testp_proc)(long,long); */
/* typedef long (*grib_op_proc)(long,long);   */


%}

%union {
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
};

%start all

%token LOWERCASE
%token IF
%token IF_TRANSIENT
%token ELSE
%token END
%token UNSIGNED
%token TEMPLATE
%token TEMPLATE_NOFAIL
%token TRIGGER
%token ASCII
%token KSEC1EXPVER
%token LABEL
%token LIST
%token WHILE
%token IBMFLOAT
%token SIGNED
%token BYTE
%token CODETABLE
%token COMPLEX_CODETABLE
%token LOOKUP
%token ALIAS
%token UNALIAS
%token META
%token POS
%token INTCONST
%token TRANS
%token FLAGBIT
%token CONCEPT
%token CONCEPT_NOFAIL
%token NIL
%token DUMMY
	
%token MODIFY

%token READ_ONLY
%token STRING_TYPE
%token LONG_TYPE
%token NO_COPY
%token DUMP
%token NO_FAIL
%token EDITION_SPECIFIC
%token OVERRIDE
%token HIDDEN
%token CAN_BE_MISSING
%token MISSING
%token CONSTRAINT
%token COPY_OK

%token WHEN
%token SET
%token SET_NOFAIL
%token WRITE
%token APPEND
%token PRINT

%token EXPORT
%token REMOVE
%token SKIP

%token PAD
%token SECTION_PADDING
%token PADTO
%token PADTOEVEN
%token PADTOMULTIPLE
%token G1_HALF_BYTE
%token G1_MESSAGE_LENGTH
%token G1_SECTION4_LENGTH
%token SECTION_LENGTH
%token FLAG
%token ITERATOR
%token NEAREST
%token BOX
%token KSEC
%token ASSERT

%token CASE
%token SWITCH
%token DEFAULT

%token EQ
%token NE
%token GE
%token LE
%token BIT
%token BITOFF

%token AND
%token OR
%token NOT

%token IS

%token <str>IDENT
%token <str>STRING

%token <lval>INTEGER
%token <dval>FLOAT

%type <dvalue> dvalues

%type <act>  instructions
%type <act>  instruction
%type <act>  simple
%type <act>  if_block
%type <act>  switch_block
%type <act>  list_block
%type <act>  while_block
%type <act>  when_block
%type <act>  trigger_block
%type <act>  concept_block
%type <act>  set
%type <act>  set_list

%type <explist> argument
%type <explist> arguments
%type <explist> argument_list
%type <explist>  default


%type <exp>     expression
%type <exp>     atom
%type <exp>     condition
%type <exp>     factor
%type <exp>     term
%type <exp>     power
%type <exp>     conjonction
%type <exp>     disjonction
%type <exp>     string_or_ident

%type <lval>     flag
%type <lval>     flags
%type <lval>     flag_list

%type <concept_condition> concept_condition
%type <concept_condition> concept_conditions

%type <concept_value> concept_value
%type <concept_value> concept_list

%type <case_value> case_value
%type <case_value> case_list

%type <rule_entry> rule_entry
%type <rule_entry> rule_entries

%type <rules> rules
%type <rules> rule
%type <rules> fact;
%type <rules> conditional_rule;


%%

all        : empty        { grib_parser_all_actions = 0;grib_parser_concept=0; grib_parser_rules=0; }
           | concept_list { grib_parser_concept     = reverse($1); }
           | instructions { grib_parser_all_actions = $1; }
           | rules        { grib_parser_rules       = $1; }
       /* memory leak here */
       | error        { grib_parser_all_actions = 0; grib_parser_concept=0; grib_parser_rules=0; }
   ;

empty:;


dvalues :  FLOAT  { $$=grib_darray_push(0,0,$1);}
    |  dvalues ',' FLOAT { $$=grib_darray_push(0,$1,$3);}
    |  INTEGER { $$=grib_darray_push(0,0,$1);}
    |  dvalues ',' INTEGER { $$=grib_darray_push(0,$1,$3);}

instructions : instruction
         | instruction instructions { $1->next = $2; $$ = $1; }
         | instruction ';'  instructions { $1->next = $3; $$ = $1; }
         | instruction ';'  {  $$ = $1;}
   ;

instruction: simple  ';'
           | if_block
           | list_block
           | while_block
           | trigger_block
           | concept_block
           | when_block
           | switch_block

   ;

 semi: ';'
    | semi ';'


argument_list: empty       { $$ = 0; }
              | arguments
              ;

arguments    : argument
              | argument ',' arguments { $1->next = $3; $$ = $1; }
              ;

argument     : expression { $$ = grib_arguments_new(grib_parser_context,$1,NULL); }
              ;


simple : UNSIGNED '[' INTEGER ']'   IDENT   default flags
        { $$ = grib_action_create_gen(grib_parser_context,$5,"unsigned",$3,NULL,$6,$7,NULL,NULL);        free($5);  }

    | UNSIGNED '[' INTEGER ']'   IDENT   '[' argument_list ']'   default  flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"unsigned",$3,$7,$9,$10,NULL,NULL);        free($5);  }

    | ASCII   '[' INTEGER ']' IDENT       default flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"ascii",$3,NULL,$6,$7,NULL,NULL);  free($5);  }

  /* Special case for '7777' */
    | ASCII   '[' INTEGER ']' STRING       default flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"ascii",$3,NULL,$6,$7,NULL,NULL);  free($5);  }

    | BYTE '[' INTEGER ']'     IDENT     default     flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"bytes",$3,NULL,$6,$7,NULL,NULL);      free($5);  }

    | BYTE '[' INTEGER ']'     IDENT   '[' argument_list ']'   default  flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"bytes",$3,$7,$9,$10,NULL,NULL);      free($5);  }

    | KSEC1EXPVER '[' INTEGER ']' IDENT     default    flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"ksec1expver",$3,NULL,$6,$7,NULL,NULL);  free($5);  }

    | SIGNED '[' INTEGER ']'     IDENT     default     flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"signed",$3,NULL,$6,$7,NULL,NULL);      free($5);  }

    | SIGNED '[' INTEGER ']'     IDENT   '[' argument_list ']'   default  flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"signed",$3,$7,$9,$10,NULL,NULL);      free($5);  }

    | CODETABLE '[' INTEGER ']' IDENT  argument   default flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"codetable",$3, $6,$7,$8,NULL,NULL);    free($5); }

	| CODETABLE '[' INTEGER ']' IDENT  argument   default SET '(' IDENT ')' flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"codetable",$3, $6,$7,$12,NULL,$10);
					free($5);free($10); }
    
    | CODETABLE '[' INTEGER ']' IDENT  '(' argument_list ')'   default flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"codetable",$3, $7,$9,$10,NULL,NULL);    free($5); }

    | COMPLEX_CODETABLE '[' INTEGER ']' IDENT  argument   default flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"complex_codetable",$3, $6,$7,$8,NULL,NULL);    free($5); }
    
    | COMPLEX_CODETABLE '[' INTEGER ']' IDENT  '(' argument_list ')'   default flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"complex_codetable",$3, $7,$9,$10,NULL,NULL);    free($5); }

    | FLAG  '[' INTEGER ']'      IDENT  argument   default flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"codeflag",$3, $6,$7,$8,NULL,NULL);  free($5); }

    | LOOKUP '[' INTEGER ']'     IDENT  '(' argument_list  ')' flags
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"lookup",$3,$7,NULL,$9,NULL,NULL); free($5); }

    | FLAGBIT  IDENT  '(' argument_list ')'  default flags
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"bit",0,$4,$6,$7,NULL,NULL); free($2); }

    | LABEL       IDENT
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"label",0,NULL,NULL,0,NULL,NULL);   free($2);  }

    | LABEL       STRING
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"label",0,NULL,NULL,0,NULL,NULL);   free($2);  }

    | IBMFLOAT    IDENT      default    flags
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"ibmfloat",4,NULL,$3,$4,NULL,NULL);free($2);  }

    | IBMFLOAT    IDENT '.' IDENT  default    flags
	{ $$ = grib_action_create_gen(grib_parser_context,$4,"ibmfloat",4,NULL,$5,$6,$2,NULL);free($4); free($2); }

    | IBMFLOAT    IDENT '[' argument ']'     default     flags
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"ibmfloat",4,$4,$6,$7,NULL,NULL);free($2);  }

    | POS         IDENT
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"position",0,NULL,NULL,0,NULL,NULL);     free($2);  }

    | INTCONST    IDENT   '=' argument flags
        { $$ = grib_action_create_variable(grib_parser_context,$2,"constant",0,$4,NULL,$5,NULL);free($2); }

    | TRANS       IDENT   '=' argument  flags
        { $$ = grib_action_create_variable(grib_parser_context,$2,"transient",0,$4,$4,$5,NULL);   free($2); }

    | FLOAT       IDENT    default   flags
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"ieeefloat",4,NULL,$3,$4,NULL,NULL);   free($2);  }

    | FLOAT       IDENT '.' IDENT   default   flags
	{ $$ = grib_action_create_gen(grib_parser_context,$4,"ieeefloat",4,NULL,$5,$6,$2,NULL);  free($4);free($2);}

   | FLOAT    IDENT '[' argument ']'     default     flags
   { $$ = grib_action_create_gen(grib_parser_context,$2,"ieeefloat",4,$4,$6,$7,NULL,NULL);free($2);  }

   | G1_HALF_BYTE  IDENT
   { $$ = grib_action_create_gen(grib_parser_context,$2,"g1_half_byte_codeflag",0,NULL,NULL,0,NULL,NULL);free($2);  }

    | SECTION_LENGTH  '[' INTEGER ']'   IDENT
	{ $$ = grib_action_create_gen(grib_parser_context,$5,"section_length",$3,NULL,NULL,0,NULL,NULL);free($5);  }

   | G1_MESSAGE_LENGTH  '[' INTEGER ']'   IDENT '(' argument_list ')'
   { $$ = grib_action_create_gen(grib_parser_context,$5,"g1_message_length",$3,$7,NULL,0,NULL,NULL);free($5);  }

  | G1_SECTION4_LENGTH  '[' INTEGER ']'   IDENT  '(' argument_list ')'
  { $$ = grib_action_create_gen(grib_parser_context,$5,"g1_section4_length",$3,$7,NULL,0,NULL,NULL);free($5);  }

    | KSEC        IDENT  argument
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"ksec",0,$3,NULL,0,NULL,NULL);free($2); }

    | PAD       IDENT  '(' argument_list ')'
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"pad",0,$4,0,0,NULL,NULL);   free($2); }

    | PADTO       IDENT  '(' argument_list ')'
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"padto",0,$4,0,0,NULL,NULL);   free($2); }

    | PADTOEVEN       IDENT  '(' argument_list ')'
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"padtoeven",0,$4,0,0,NULL,NULL);   free($2); }

    | PADTOMULTIPLE       IDENT  '(' argument_list ')'
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"padtomultiple",0,$4,0,0,NULL,NULL);   free($2); }

    | SECTION_PADDING     IDENT  flags
	{ $$ = grib_action_create_gen(grib_parser_context,$2,"section_padding",0,0,0,$3,NULL,NULL);   free($2);  }

    | TEMPLATE    IDENT  STRING
        { $$ = grib_action_create_template(grib_parser_context,0,$2,$3); free($2); free($3);}
    | TEMPLATE_NOFAIL    IDENT  STRING
    { $$ = grib_action_create_template(grib_parser_context,1,$2,$3); free($2); free($3);}

    | ALIAS  IDENT '=' IDENT flags
        { $$ = grib_action_create_alias(grib_parser_context,$2,$4,NULL,$5);  free($2); free($4); }

    | UNALIAS  IDENT 
        { $$ = grib_action_create_alias(grib_parser_context,$2,NULL,NULL,0);  free($2); }

    | ALIAS  IDENT '.' IDENT '=' IDENT flags
        {
         $$ = grib_action_create_alias(grib_parser_context,$4,$6,$2,$7);  free($2); free($4); free($6);
    }
    | UNALIAS  IDENT '.' IDENT
        {
         $$ = grib_action_create_alias(grib_parser_context,$4,NULL,$2,0);  free($2); free($4); 
    }
    | META IDENT  IDENT '(' argument_list ')'  default flags
        { $$ = grib_action_create_meta(grib_parser_context,$2,$3,$5,$7,$8,NULL); free($2);free($3);}

    | META IDENT '.' IDENT IDENT '(' argument_list ')'  default flags
    { $$ = grib_action_create_meta(grib_parser_context,$4,$5,$7,$9,$10,$2); free($4);free($5);free($2);}

    | ITERATOR  IDENT '(' argument_list ')'
        {
      grib_arguments* a = grib_arguments_new(
        grib_parser_context,
        new_accessor_expression(grib_parser_context,$2),
		NULL
        );
      a->next=$4;
      $$ = grib_action_create_meta(grib_parser_context,
      "ITERATOR","iterator",a,NULL,
      GRIB_ACCESSOR_FLAG_HIDDEN|GRIB_ACCESSOR_FLAG_READ_ONLY,NULL); free($2);
    }
    | NEAREST  IDENT '(' argument_list ')'
        {
      grib_arguments* a = grib_arguments_new(
        grib_parser_context,
        new_accessor_expression(grib_parser_context,$2),
		NULL
        );
      a->next=$4;
      $$ = grib_action_create_meta(grib_parser_context,
      "NEAREST","nearest",a,NULL,
      GRIB_ACCESSOR_FLAG_HIDDEN|GRIB_ACCESSOR_FLAG_READ_ONLY,NULL); free($2);
    }
    | BOX  IDENT '(' argument_list ')'
        {
      grib_arguments* a = grib_arguments_new(
        grib_parser_context,
        new_accessor_expression(grib_parser_context,$2),
		NULL
        );
      a->next=$4;
      $$ = grib_action_create_meta(grib_parser_context,
      "BOX","box",a,NULL,
      GRIB_ACCESSOR_FLAG_HIDDEN|GRIB_ACCESSOR_FLAG_READ_ONLY,NULL); free($2);
    }
    | EXPORT IDENT '(' argument_list ')'
       { $$ = grib_action_create_put(grib_parser_context,$2,$4);free($2);}

    | REMOVE argument_list
       { $$ = grib_action_create_remove(grib_parser_context,$2);}

    | ASSERT '(' expression ')'
       { $$ = grib_action_create_assert(grib_parser_context,$3);}

    | MODIFY IDENT flags
       { $$ = grib_action_create_modify(grib_parser_context,$2,$3); free($2);}

  | SET IDENT '=' MISSING { $$ = grib_action_create_set_missing(grib_parser_context,$2); free($2); }
  | SET IDENT '=' expression { $$ = grib_action_create_set(grib_parser_context,$2,$4,0); free($2); }
  | SET IDENT '=' '{' dvalues '}' { $$ = grib_action_create_set_darray(grib_parser_context,$2,$5); free($2); }

  | SET_NOFAIL IDENT '=' expression { $$ = grib_action_create_set(grib_parser_context,$2,$4,1); free($2); }

  
  | WRITE STRING { $$ = grib_action_create_write(grib_parser_context,$2,0); free($2);}
  | WRITE { $$ = grib_action_create_write(grib_parser_context,"",0); }
  | APPEND STRING { $$ = grib_action_create_write(grib_parser_context,$2,1); free($2);}
  | APPEND { $$ = grib_action_create_write(grib_parser_context,"",1); }

  | PRINT STRING { $$ = grib_action_create_print(grib_parser_context,$2,0); free($2); }
  | PRINT '(' STRING ')' STRING { $$ = grib_action_create_print(grib_parser_context,$5,$3); free($5); free($3);}
  | PRINT { $$ = grib_action_create_print(grib_parser_context,"",0);  }
   ;

if_block :
  IF '(' expression ')' '{' instructions '}' { $$ = grib_action_create_if(grib_parser_context,$3,$6,0,0); }
| IF '(' expression ')' '{' instructions '}' ELSE '{' instructions '}'  { $$ = grib_action_create_if(grib_parser_context,$3,$6,$10,0); }
| IF_TRANSIENT '(' expression ')' '{' instructions '}' { $$ = grib_action_create_if(grib_parser_context,$3,$6,0,1); }
| IF_TRANSIENT '(' expression ')' '{' instructions '}' ELSE '{' instructions '}'  { $$ = grib_action_create_if(grib_parser_context,$3,$6,$10,1); }
   ;

when_block :
  WHEN '(' expression ')' set semi   { $$ = grib_action_create_when(grib_parser_context,$3,$5,NULL); }
  | WHEN '(' expression ')' '{' set_list '}'   { $$ = grib_action_create_when(grib_parser_context,$3,$6,NULL); }
  | WHEN '(' expression ')' '{' set_list '}' ELSE '{' set_list '}' { $$ = grib_action_create_when(grib_parser_context,$3,$6,$10); }
  ;

set : SET IDENT '=' expression { $$ = grib_action_create_set(grib_parser_context,$2,$4,0); free($2); }
  | SET_NOFAIL IDENT '=' expression { $$ = grib_action_create_set(grib_parser_context,$2,$4,1); free($2); }
  ;

set_list : set semi
         | set_list set semi { $1->next = $2; $$ = $1; }
         ;


default : empty { $$ = NULL ;}
   | '=' argument_list { $$ = $2 ;}
   ;

flags : empty         { $$ = 0 ; }
      | ':' flag_list { $$ = $2; }
      ;

flag_list  : flag
   | flag_list ',' flag { $$ = $1 | $3; }
   ;

flag: READ_ONLY         { $$ = GRIB_ACCESSOR_FLAG_READ_ONLY; }
    | LOWERCASE            { $$ = GRIB_ACCESSOR_FLAG_LOWERCASE; }
    | DUMP            { $$ = GRIB_ACCESSOR_FLAG_DUMP; }
    | NO_COPY            { $$ = GRIB_ACCESSOR_FLAG_NO_COPY; }
	| NO_FAIL            { $$ = GRIB_ACCESSOR_FLAG_NO_FAIL; }
    | HIDDEN            { $$ = GRIB_ACCESSOR_FLAG_HIDDEN; }
    | EDITION_SPECIFIC  { $$ = GRIB_ACCESSOR_FLAG_EDITION_SPECIFIC; }
    | CAN_BE_MISSING    { $$ = GRIB_ACCESSOR_FLAG_CAN_BE_MISSING; }
    | CONSTRAINT        { $$ = GRIB_ACCESSOR_FLAG_CONSTRAINT; }
    | OVERRIDE           { $$ = GRIB_ACCESSOR_FLAG_OVERRIDE; }
    | COPY_OK           { $$ = GRIB_ACCESSOR_FLAG_COPY_OK; }
    | TRANS         { $$ = GRIB_ACCESSOR_FLAG_TRANSIENT; }
    | STRING_TYPE         { $$ = GRIB_ACCESSOR_FLAG_STRING_TYPE; }
    | LONG_TYPE         { $$ = GRIB_ACCESSOR_FLAG_LONG_TYPE; }
    ;

list_block : IDENT LIST '(' expression ')' '{' instructions '}' { $$ = grib_action_create_list(grib_parser_context,$1,$4,$7); free($1); }
  ;

while_block : WHILE '(' expression ')' '{' instructions '}' { $$ = grib_action_create_while(grib_parser_context,$3,$6);  }
  ;

trigger_block : TRIGGER '(' argument_list ')' '{' instructions '}' { $$ = grib_action_create_trigger(grib_parser_context,$3,$6);  }
   ;

concept_block : CONCEPT IDENT '{' concept_list '}' flags { $$ = grib_action_create_concept(grib_parser_context,$2,$4,0,0,0,0,0,0,$6,0);  free($2); }
   | CONCEPT IDENT '(' IDENT ')' '{' concept_list '}' flags { $$ = grib_action_create_concept(grib_parser_context,$2,$7,0,0,$4,0,0,0,$9,0);  free($2);free($4); }
   | CONCEPT IDENT '(' IDENT ',' STRING ',' IDENT ',' IDENT ')' flags { $$ = grib_action_create_concept(grib_parser_context,$2,0,$6,0,$4,$8,$10,0,$12,0);  free($2);free($6);free($4);free($8);free($10); }
   | CONCEPT IDENT '(' IDENT ',' STRING ',' IDENT ',' IDENT ',' IDENT ')' flags { $$ = grib_action_create_concept(grib_parser_context,$2,0,$6,0,$4,$8,$10,$12,$14,0);  free($2);free($6);free($4);free($8);free($10);free($12); }
   | CONCEPT IDENT '(' IDENT ',' STRING ',' IDENT ')' flags { $$ = grib_action_create_concept(grib_parser_context,$2,0,$6,0,$4,$8,0,0,$10,0);  free($2);free($6);free($4);free($8); }
   | CONCEPT IDENT '.' IDENT '(' IDENT ',' STRING ',' IDENT ',' IDENT ')' flags { $$ = grib_action_create_concept(grib_parser_context,$4,0,$8,$2,$6,$10,$12,0,$14,0);  free($4);free($8);free($6);free($10); free($2);}
   | CONCEPT IDENT '.' IDENT '(' IDENT ',' STRING ',' IDENT ')' flags { $$ = grib_action_create_concept(grib_parser_context,$4,0,$8,$2,$6,$10,0,0,$12,0);  free($4);free($8);free($6);free($10); free($2);}
   | CONCEPT IDENT '.' IDENT '{' concept_list '}' flags { $$ = grib_action_create_concept(grib_parser_context,$4,$6,0,$2,0,0,0,0,$8,0);  free($2);free($4); }
   | CONCEPT IDENT '.' IDENT '(' IDENT ')' '{' concept_list '}' flags { $$ = grib_action_create_concept(grib_parser_context,$4,$9,0,$2,$6,0,0,0,$11,0);  free($2);free($4);free($6); }
   |CONCEPT_NOFAIL IDENT '{' concept_list '}' flags { $$ = grib_action_create_concept(grib_parser_context,$2,$4,0,0,0,0,0,0,$6,1);  free($2); }
   | CONCEPT_NOFAIL IDENT '(' IDENT ')' '{' concept_list '}' flags { $$ = grib_action_create_concept(grib_parser_context,$2,$7,0,0,$4,0,0,0,$9,1);  free($2);free($4); }
   | CONCEPT_NOFAIL IDENT '(' IDENT ',' STRING ',' IDENT ',' IDENT ')' flags { $$ = grib_action_create_concept(grib_parser_context,$2,0,$6,0,$4,$8,$10,0,$12,1);  free($2);free($6);free($4);free($8);free($10); }
   | CONCEPT_NOFAIL IDENT '(' IDENT ',' STRING ',' IDENT ')' flags { $$ = grib_action_create_concept(grib_parser_context,$2,0,$6,0,$4,$8,0,0,$10,1);  free($2);free($6);free($4);free($8); }
   | CONCEPT_NOFAIL IDENT '.' IDENT '(' IDENT ',' STRING ',' IDENT ',' IDENT ')' flags { $$ = grib_action_create_concept(grib_parser_context,$4,0,$8,$2,$6,$10,$12,0,$14,1);  free($4);free($8);free($6);free($10);free($12); free($2);}
   | CONCEPT_NOFAIL IDENT '.' IDENT '(' IDENT ',' STRING ',' IDENT ')' flags { $$ = grib_action_create_concept(grib_parser_context,$4,0,$8,$2,$6,$10,0,0,$12,1);  free($4);free($8);free($6);free($10); free($2);}
   | CONCEPT_NOFAIL IDENT '.' IDENT '{' concept_list '}' flags { $$ = grib_action_create_concept(grib_parser_context,$4,$6,0,$2,0,0,0,0,$8,1);  free($2);free($4); }
   | CONCEPT_NOFAIL IDENT '.' IDENT '(' IDENT ')' '{' concept_list '}' flags { $$ = grib_action_create_concept(grib_parser_context,$4,$9,0,$2,$6,0,0,0,$11,1);  free($2);free($4);free($6); }
   
   ;

concept_list : concept_value
             | concept_list concept_value { $$ = $2; $2->next = $1;   }
       ;

case_list : case_value
             | case_list case_value { $$ = $2; $2->next = $1;   }
       ;

case_value :  CASE arguments ':' instructions  { $$ = grib_case_new(grib_parser_context,$2,$4);  }
        ;

switch_block :
  SWITCH '(' argument_list ')' '{' case_list DEFAULT ':' instructions  '}' { $$ = grib_action_create_switch(grib_parser_context,$3,$6,$9); }
  | SWITCH '(' argument_list ')' '{' case_list DEFAULT ':' '}' { $$ = grib_action_create_switch(grib_parser_context,$3,$6,grib_action_create_noop(grib_parser_context,"continue")); }
  | SWITCH '(' argument_list ')' '{' case_list '}' { $$ = grib_action_create_switch(grib_parser_context,$3,$6,0); }
  ;

  concept_value :  STRING '=' '{' concept_conditions '}' {
	  				$$ = grib_concept_value_new(grib_parser_context,$1,$4); free($1);}
  				| IDENT '=' '{' concept_conditions '}' {
	  				$$ = grib_concept_value_new(grib_parser_context,$1,$4); free($1);}
				| INTEGER '=' '{' concept_conditions '}' {
					char buf[80]; sprintf(buf,"%ld",(long)$1); $$ = grib_concept_value_new(grib_parser_context,buf,$4);}
				| FLOAT '=' '{' concept_conditions '}' {
					char buf[80]; sprintf(buf,"%g",(double)$1); $$ = grib_concept_value_new(grib_parser_context,buf,$4);}
        ;

concept_conditions : concept_condition
                | concept_condition concept_conditions { $1->next = $2; $$ = $1; }
        ;

concept_condition   : IDENT '=' expression ';' { $$ = grib_concept_condition_new(grib_parser_context,$1,$3); free($1); }


string_or_ident :  IDENT   { $$ = new_accessor_expression(grib_parser_context,$1); free($1); }
                | STRING  { $$ = new_string_expression(grib_parser_context,$1);  free($1); }
                ;

atom  : string_or_ident
      | INTEGER { $$ = new_long_expression(grib_parser_context,$1);  }
      | FLOAT { $$ = new_double_expression(grib_parser_context,$1);  /* TODO: change to new_float_expression*/}

    | NIL     { $$ = NULL; }
	| DUMMY     { $$ = new_true_expression(grib_parser_context); }
      | '(' expression ')' { $$ = $2; }
      | '-' atom { $$ = new_unop_expression(grib_parser_context,&grib_op_neg,&grib_op_neg_d,$2); }
    | IDENT '(' ')' { $$ = new_func_expression(grib_parser_context,$1,NULL); free($1);}
    | IDENT '(' argument_list ')' { $$ = new_func_expression(grib_parser_context,$1,$3); free($1);}
      ;


power          : atom '^' power     { $$ = new_binop_expression(grib_parser_context,&grib_op_pow,NULL,$1,$3); }
               | atom
               ;

factor         : factor '*' power    { $$ = new_binop_expression(grib_parser_context,&grib_op_mul,&grib_op_mul_d,$1,$3); }
               | factor '/' power    { $$ = new_binop_expression(grib_parser_context,&grib_op_div,&grib_op_div_d,$1,$3); }
               | factor '%' power    { $$ = new_binop_expression(grib_parser_context,&grib_op_modulo,NULL,$1,$3); }
            | factor BIT  power   { $$ = new_binop_expression(grib_parser_context,&grib_op_bit,NULL,$1,$3); }
            | factor BITOFF power { $$ = new_binop_expression(grib_parser_context,&grib_op_bitoff,NULL,$1,$3); }
               | power
               ;

term           : term '+' factor    { $$ = new_binop_expression(grib_parser_context,&grib_op_add,&grib_op_add_d,$1,$3); }
               | term '-' factor    { $$ = new_binop_expression(grib_parser_context,&grib_op_sub,&grib_op_sub_d,$1,$3); }
               | factor
              ;

condition     : condition '>'     term { $$ = new_binop_expression(grib_parser_context,&grib_op_gt,&grib_op_gt_d,$1,$3); }
             /* | condition '='     term { $$ = new_binop_expression(grib_parser_context,&grib_op_eq,$1,$3); } */
             | condition EQ     term { $$ = new_binop_expression(grib_parser_context,&grib_op_eq,&grib_op_eq_d,$1,$3); }
             | condition '<'     term { $$ = new_binop_expression(grib_parser_context,&grib_op_lt,&grib_op_lt_d,$1,$3); }
             | condition  GE     term { $$ = new_binop_expression(grib_parser_context,&grib_op_ge,&grib_op_ge_d,$1,$3); }
             | condition  LE     term { $$ = new_binop_expression(grib_parser_context,&grib_op_le,&grib_op_le_d,$1,$3); }
             | condition  NE     term { $$ = new_binop_expression(grib_parser_context,&grib_op_ne,&grib_op_ne_d,$1,$3); }
             | string_or_ident IS string_or_ident { $$ = new_string_compare_expression(grib_parser_context,$1,$3); }
/*
             | condition  IN     term { $$ = new_binop_expression(grib_parser_context,grib_op_pow,$1,$3); }
             | condition  MATCH  term { $$ = new_binop_expression(grib_parser_context,grib_op_pow,$1,$3); }
*/
             | NOT condition          { $$ = new_unop_expression(grib_parser_context,&grib_op_not,NULL,$2); }
            | term
             ;

conjonction : conjonction AND condition { $$ = new_binop_expression(grib_parser_context,&grib_op_and,NULL,$1,$3); }
            | condition
            ;

disjonction    : disjonction OR conjonction { $$ = new_binop_expression(grib_parser_context,&grib_op_or,NULL,$1,$3);}
            | conjonction
            ;

expression     : disjonction
            ;


/*  */

rule : fact
     | conditional_rule
   ;

rule_entry : IDENT '=' expression ';'  { $$ = grib_new_rule_entry(grib_parser_context,$1,$3); free($1); }
            | SKIP ';' { $$ = grib_new_rule_entry(grib_parser_context,"skip",0);}
           ;

rule_entries : rule_entry
             | rule_entry rule_entries { $1->next = $2; $$ = $1; }
       ;

fact: rule_entry  { $$ = grib_new_rule(grib_parser_context,NULL,$1); }
    ;


conditional_rule: IF '(' expression ')' '{' rule_entries '}' { $$ = grib_new_rule(grib_parser_context,$3,$6); }
                ;

rules : rule
      | rule rules { $1->next = $2; $$ = $1; }
    ;


%%

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



