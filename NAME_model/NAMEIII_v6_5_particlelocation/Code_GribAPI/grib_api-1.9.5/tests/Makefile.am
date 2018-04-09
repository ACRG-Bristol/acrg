
     TESTS = definitions.sh \
			 ieee.sh ls.sh grib1to2.sh \
	         grib2to1.sh badgrib.sh convert.sh filter.sh  \
	         multi.sh budg.sh gridType.sh concept.sh decimalPrecision.sh \
			 bitsPerValue.sh get_fail.sh missing.sh \
			 step.sh set.sh iterator.sh \
			 compare.sh level.sh index.sh \
			 bitmap.sh list.sh \
			 multi_from_message.sh change_scanning.sh \
			 julian.sh statistics.sh tigge.sh tigge_conversions.sh \
			 read_any.sh debug.sh $(JPEG_TEST)

     noinst_PROGRAMS = gauss_sub read_any double_cmp packing pack_unpack \
	                   multi_from_message julian read_index index gribex_perf\
	                   jpeg_perf so_perf png_perf bpv_limit laplacian

	 multi_from_message_SOURCES = multi_from_message.c
	 laplacian_SOURCES = laplacian.c
	 packing_SOURCES = packing.c
	 list_keys_SOURCES = list_keys.c
	 read_any_SOURCES = read_any.c
	 double_cmp_SOURCES = double_cmp.c
	 pack_unpack_SOURCES = pack_unpack.c
	 julian_SOURCES = julian.c
	 index_SOURCES = index.c
	 read_index_SOURCES = read_index.c
	 jpeg_perf_SOURCES = jpeg_perf.c
	 so_perf_SOURCES = so_perf.c
	 png_perf_SOURCES = png_perf.c
	 gribex_perf_SOURCES = gribex_perf.c
	 gribex_perf_SOURCES = gauss_sub.c

     LDADD = $(top_builddir)/src/libgrib_api.a $(EMOS_LIB)

     INCLUDES   = -I$(top_builddir)/src

     EXTRA_DIST = $(TESTS) mf.rules filter_rules include.sh
   
