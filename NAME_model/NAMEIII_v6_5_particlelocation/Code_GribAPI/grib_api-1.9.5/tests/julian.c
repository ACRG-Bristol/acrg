#include "grib_api.h"

int main(int argc,char* argv[]) {
  long year,month,day,hour,minute,second;
  int i,last;
  double jd;
  long jdl,date;
  double jds[]={2451545.0,2446822.5,2446966.0,2447187.5,2447332.0,2415020.5,2305447.5,
  2305812.5,2205512.5,2026871.8,1356001.0,1355866.5,1355671.4,0.0};

  last=sizeof(jds)/sizeof(double);

  year=1957; 
  month=10;
  day=4;
  hour=19;
  minute=26;
  second=24;

  grib_datetime_to_julian(year,month,day,hour,minute,second,&jd);

  printf("%ld %ld %ld %ld:%ld:%ld -> %f\n",year,month,day,hour,minute,second,jd);

  grib_julian_to_datetime(jd,&year,&month,&day,&hour,&minute,&second);

  printf("%ld %ld %ld %ld:%ld:%ld -> %f\n",year,month,day,hour,minute,second,jd);

  printf("\n");
  for (i=0;i<last;i++) {
    jd=jds[i];
    grib_julian_to_datetime(jd,&year,&month,&day,&hour,&minute,&second);
    printf("+ %ld %ld %ld %ld:%ld:%ld -> %f\n",year,month,day,hour,minute,second,jd);
    grib_datetime_to_julian(year,month,day,hour,minute,second,&jd);
    printf("- %ld %ld %ld %ld:%ld:%ld -> %f\n",year,month,day,hour,minute,second,jd);
	jdl=(long)(jd+0.5);
    date=grib_julian_to_date(jdl);
    printf("+ %ld -> %ld\n",date,jdl);
    jdl=grib_date_to_julian(date);
    printf("- %ld -> %ld\n",date,jdl);
	printf("\n");
  }

  return 0;
}
