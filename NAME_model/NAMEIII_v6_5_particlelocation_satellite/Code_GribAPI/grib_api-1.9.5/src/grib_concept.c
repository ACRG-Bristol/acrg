/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

#include "grib_api_internal.h"


grib_concept_value* grib_concept_value_new(grib_context* c, const char* name,grib_concept_condition* conditions)
{
	grib_concept_value* v = grib_context_malloc_clear_persistent(c,sizeof(grib_concept_value));
	
	v->name = grib_context_strdup_persistent(c,name);
	v->conditions = conditions;
	return v;
}

void grib_concept_value_delete(grib_context* c,grib_concept_value* v)
{
	grib_concept_condition* e = v->conditions;
	while(e)
	{
		grib_concept_condition* n = e->next;
		grib_concept_condition_delete(c,e);
		e = n;
	}
	grib_context_free_persistent(c,v->name);
	grib_context_free_persistent(c,v);
}

grib_concept_condition* grib_concept_condition_new(grib_context* c,const char* name,grib_expression* expression)
{
	grib_concept_condition* v = grib_context_malloc_clear_persistent(c,sizeof(grib_concept_condition));
	v->name = grib_context_strdup_persistent(c,name);
	v->expression = expression;
	return v;
}

void grib_concept_condition_delete(grib_context* c,grib_concept_condition* v)
{
	grib_expression_free(c,v->expression);
	grib_context_free_persistent(c,v->name);
	grib_context_free_persistent(c,v);
}
