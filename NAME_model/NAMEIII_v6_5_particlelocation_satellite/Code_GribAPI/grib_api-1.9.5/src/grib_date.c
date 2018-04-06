/**
* Copyright 2005-2007 ECMWF
*
* Licensed under the GNU Lesser General Public License which
* incorporates the terms and conditions of version 3 of the GNU
* General Public License.
* See LICENSE and gpl-3.0.txt for details.
*/

/***************************************************************************
 *   Enrico Fucile                                                            *
 *                                                                         *
 ***************************************************************************/
#include "grib_api_internal.h"

int grib_julian_to_datetime(double jd,long *year,long* month,long* day,
                            long *hour,long *minute,long *second)
{
  long z,a,alpha,b,c,d,e;
  double dday,dhour,dminute;
  double f;
  
  jd+=0.5;
  z=(long)jd;
  f=jd-z;

  if (z < 2299161) a=z;
  else {
    alpha=(long)((z-1867216.25)/36524.25);
    a=z+1+alpha-(long)(((double)alpha)/4);
  }
  b=a+1524;
  c=(long)((b-122.1)/365.25);
  d=(long)(365.25*c);
  e=(long)(((double)(b-d))/30.6001);

  dday=b-d-(long)(30.6001*e)+f;
  *day=(long)dday;
  dday-=*day;
  dhour=dday*24;
  *hour=(long)dhour;
  dhour-=*hour;
  dminute=dhour*60;
  *minute=(long)dminute;
  *second=(long)((dminute-*minute)*60);
  
  if (e<14) *month=e-1;
  else *month=e-13;

  if (*month>2) *year=c-4716;
  else *year=c-4715;

  return GRIB_SUCCESS;
}

int grib_datetime_to_julian(long year,long month,long day,
                            long hour,long minute,long second,double* jd)
{
  double a,b,dday;
  long y,m;

  dday=(double)(hour*3600+minute*60+second)/86400.0+day;
  
  if ( month < 3) {
    y=year-1;
    m=month+12;
  }else {
    y=year;
    m=month;
  }
  a=(long)(((double)y)/100);

  if ( y > 1582 ) b=2-a+(long)(a/4);
  else if (y == 1582) {
    if (m > 10) b=2-a+(long)(a/4);
    else if (m == 10) {
      if (day >14) b=2-a+(long)(a/4);
      else b=0;
    }
    else b=0;
  }
  else b=0;

  *jd=(long)(365.25*(y+4716))+ (long)(30.6001*(m+1))+dday+b-1524.5;

  return GRIB_SUCCESS;
}


long grib_julian_to_date(long jdate)
{
  long x,y,d,m,e;
  long day,month,year;

  x = 4 * jdate - 6884477;
  y = (x / 146097) * 100;
  e = x % 146097;
  d = e / 4;

  x = 4 * d + 3;
  y = (x / 1461) + y;
  e = x % 1461;
  d = e / 4 + 1;

  x = 5 * d - 3;
  m = x / 153 + 1;
  e = x % 153;
  d = e / 5 + 1;

  if( m < 11 )
    month = m + 2;
  else
    month = m - 10;


  day = d;
  year = y + m / 11;

  return year * 10000 + month * 100 + day;

}

long grib_date_to_julian(long ddate)
{
  long  m1,y1,a,b,c,d,j1;

  long month,day,year;

  /*Asserts(ddate > 0);*/

  year = ddate / 10000;
  ddate %= 10000;
  month  = ddate / 100;
  ddate %= 100;
  day = ddate;


/*  if (year < 100) year = year + 1900; */

  if (month > 2)
  {
    m1 = month - 3;
    y1 = year;
  }
  else
  {
    m1 = month + 9;
    y1 = year - 1;
  }
  a = 146097*(y1/100)/4;
  d = y1 % 100;
  b = 1461*d/4;
  c = (153*m1+2)/5+day+1721119;
  j1 = a+b+c;

  return(j1);
}

/*
   void basedate_to_verifydate(gribsec1 *s1,request *r)
   {
   int bdate,vdate,cdate,vd,vtime,vstep;

   bdate = (s1->century-1)*1000000 + s1->year * 10000 + s1->month * 100 + s1->day;
   cdate = date_to_julian (bdate);
   vtime = cdate * 24 + s1->p1*units[s1->time_unit];
   vd    = vtime / 24;
   vdate = julian_to_date (vd,mars.y2k);
   vtime = vtime % 24;
   vstep = 0;

   set_value(r,"DATE","%d",vdate);
   set_value(r,"TIME","%02d00",vtime);
   set_value(r,"STEP","%d",vstep);
   }
 */

