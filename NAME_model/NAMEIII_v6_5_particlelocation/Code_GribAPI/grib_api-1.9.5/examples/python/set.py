import traceback
import sys

from gribapi import *
from datetime import date

INPUT='../../data/regular_latlon_surface_constant.grib1'
OUTPUT='out.grib'
VERBOSE=1 # verbose error reporting

def example():
    fin = open(INPUT)
    fout = open(OUTPUT,'w')

    gid = grib_new_from_file(fin)

    dt = date.today()
    today = "%d%02d%02d" % (dt.year,dt.month,dt.day)
    grib_set(gid,'date',int(today))

    grib_set(gid,'centre',80)

    print 'get centre as a integer - centre = %d' % grib_get(gid,'centre',int)
    print 'get centre as a string - centre = %s' % grib_get(gid,'centre',str)
    print 'get date as a string - date = %s' % grib_get(gid,'date',str)

    grib_write(gid,fout)

    grib_release(gid)
    fin.close()
    fout.close()

def main():
    try:
        example()
    except GribInternalError,err:
        if VERBOSE:
            traceback.print_exc(file=sys.stderr)
        else:
            print >>sys.stderr,err.msg

        return 1

if __name__ == "__main__":
    sys.exit(main())
