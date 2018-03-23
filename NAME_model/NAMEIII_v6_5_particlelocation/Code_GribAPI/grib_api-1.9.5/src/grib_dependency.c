/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/***************************************************************************
 *   Jean Baptiste Filippi - 01.11.2005
 *   Enrico Fucile
 *                                                                         *
 ***************************************************************************/
#include "grib_api_internal.h"

static grib_handle* handle_of(grib_accessor* observed)
{
  grib_handle *h = observed->parent->h;
  while(h->main) h = h->main;
  return h;
}

void grib_dependency_add(grib_accessor* observer,grib_accessor* observed)
{
  grib_handle     *h = handle_of(observed);
  grib_dependency *d = h->dependencies;
  grib_dependency *last = 0;

/*
  printf("observe %p %p %s %s\n",(void*)observed,(void*)observer, observed?observed->name:"NULL",
    observer?observer->name:"NULL");
*/

  if(!observer || !observed)
  {
    return;
  }

  /* Assert(h == handle_of(observer)); */

  /* Check if already in list */
  while(d)
  {
    if(d->observer == observer && d->observed == observed)
      return;
    last = d;
    d = d->next;
  }

#if 0
  d = h->dependencies;
  while(d)
  {
    last = d;
    d = d->next;
  }

#endif

  d = (grib_dependency*)grib_context_malloc_clear(h->context,sizeof(grib_dependency));
  Assert(d);

  d->observed = observed;
  d->observer = observer;
  d->next     = 0;

/*
  printf("observe %p %p %s %s\n",(void*)observed,(void*)observer, observed->name,observer->name);
*/
#if 0
  d->next     = h->dependencies;
  h->dependencies = d;
#endif

  if(last)
    last->next = d;
  else
    h->dependencies = d;

}

void grib_dependency_remove_observed(grib_accessor* observed)
{
  grib_handle     *h = handle_of(observed);
  grib_dependency *d = h->dependencies;
  /* printf("%s\n",observed->name); */

  while(d)
  {
	if(d->observed == observed)
    {
      /*  TODO: Notify observer...*/
      d->observed = 0;/*
   printf("grib_dependency_remove_observed %s\n",observed->name); */
    }
    d = d->next;
  }
}

/* TODO: Notification must go from outer blocks to inner block */

int grib_dependency_notify_change(grib_accessor* observed)
{
  grib_handle     *h = handle_of(observed);
  grib_dependency *d = h->dependencies;
  int ret = GRIB_SUCCESS;

  /* Do a two pass mark&sweep, in case some depedancies
  are added while we notify */

  while(d)
  {
    d->run = (d->observed == observed && d->observer != 0);
    d = d->next;
  }

  d = h->dependencies;
  while(d)
  {
    if(d->run)
    {
         /*printf("grib_dependency_notify_change %s %s %p\n",observed->name,d->observer ? d->observer->name : "?", (void*)d->observer);*/
      if( d->observer && (ret = grib_accessor_notify_change(d->observer,observed))
          != GRIB_SUCCESS) return ret;
    }
    d = d->next;
  }
  return ret;
}

void grib_dependency_remove_observer(grib_accessor* observer)
{
  grib_handle     *h = handle_of(observer);
  grib_dependency *d = h->dependencies;

  if (!observer) return;
  
  while(d)
  {
	if(d->observer == observer)
    {
      d->observer = 0;
	}
    d = d->next;
  }
}

void grib_dependency_observe_expression(grib_accessor* observer,grib_expression* e)
{
  grib_expression_add_dependency(e,observer);
}

void grib_dependency_observe_arguments(grib_accessor* observer,grib_arguments* a)
{
  while(a)
  {
    grib_dependency_observe_expression(observer,a->expression);
    a = a->next;
  }
}
