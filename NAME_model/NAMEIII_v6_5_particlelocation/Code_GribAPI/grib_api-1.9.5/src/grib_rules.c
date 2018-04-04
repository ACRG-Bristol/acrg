/**
* Copyright 2005-2007 ECMWF
* 
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

#include "grib_api_internal.h"

grib_rule_entry *grib_new_rule_entry(grib_context* c,const char* name,grib_expression* expression)
{
	grib_rule_entry* e = grib_context_malloc_clear_persistent(c,sizeof(grib_rule_entry));
	e->name            = grib_context_strdup_persistent(c,name);
	e->value           = expression;
	return e;
}

static void grib_rule_entry_delete(grib_context* c,grib_rule_entry* e)
{
	if(e)
	{
		grib_context_free_persistent(c,e->name);
		grib_rule_entry_delete(c,e->next);
		grib_context_free_persistent(c,e);

	}
}

grib_rule* grib_new_rule(grib_context* c,grib_expression* condition,grib_rule_entry* entries)
{
	grib_rule* r = grib_context_malloc_clear_persistent(c,sizeof(grib_rule));
	r->condition = condition;
	r->entries   = entries;
	return r;
}


