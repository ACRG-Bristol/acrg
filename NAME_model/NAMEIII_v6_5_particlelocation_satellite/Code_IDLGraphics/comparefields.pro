pro comparefields, $
  datadir1,        $
  datadir2,        $
  filename,        $
  timelabel,       $
  selgrid,         $
  selfield,        $
  sellevel,        $
  selspecies,      $
  selname

;--------------------------------------------------------------
; procedure to compare NAME III output field which are
; stored in the files with the same name and different
; directories. 
;
; Eike Mueller 05/08/2010
;--------------------------------------------------------------

; inputs:
;  datadir1           :  data directory with first data file
;  datadir2           :  data directory with second data file
;  filename           :  name of postscript output file
;  timelabel          :  output time label           }
;  selgrid            :  output grid                 }  Defines name
;  selfield           :  output field                }  of NAME data file
;  sellevel           :  output level                }
;  selspecies         :  output species              }
;  selname            :  name of output requirement  }

;  This script compares the two fields, calculates their absolute and
;  relative difference and computes statistics on the difference fields.
;
;  It generates a postscript file <filename> which shows
;
;  1st page:
;  ---------
;  * the fields from the two directories on a logarithmic scale
;    in the top row
;  * the absolute difference between the fields in the lower left
;    corner (log scale)
;  * the relative difference between the fields in the lower right 
;    corner (log scale).
;  * mean, stdev and skewness of absolute difference distribution 
;    on the second page
;
;  2nd page:
;  ---------
;  * histogram of the distribution of the absolute difference
;
;  3rd page:
;  ---------
;  * cumulative histogram of the relative difference, i.e. how big is 
;    the percentage of gridboxes whose relative difference is smaller
;    than a given threshold?
;
;--------------------------------------------------------------

;-----------------------------------------------------------------------  
; Read data from files

  readfield,         $
    datadir1,        $
    timelabel,       $
    selgrid,         $
    selfield,        $
    sellevel,        $
    selspecies,      $
    lon2d_1,         $
    lat2d_1,         $
    fielddata_1,     $
    selname=selname, $
    namever=3

  readfield,         $
    datadir2,        $
    timelabel,       $
    selgrid,         $
    selfield,        $
    sellevel,        $
    selspecies,      $
    lon2d_2,         $
    lat2d_2,         $
    fielddata_2,     $
    selname=selname, $
    namever=3

;-----------------------------------------------------------------------  
; process data, calculate absolute and relative difference and
; convert to logarithmic scale
 
  fielddata_1 = REVERSE(fielddata_1,2)
  fielddata_2 = REVERSE(fielddata_2,2)

  abs_difference=ABS(fielddata_1-fielddata_2)
  signed_difference=fielddata_1-fielddata_2
  rel_difference=abs_difference

  where_1_GT_2 = WHERE(fielddata_1 GT fielddata_2)
  where_2_GT_1 = WHERE(fielddata_2 GT fielddata_1)
  
  IF (where_1_GT_2 NE [-1]) THEN BEGIN
    rel_difference[where_1_GT_2] = $
      abs_difference[where_1_GT_2]/fielddata_1[where_1_GT_2]
  ENDIF
  IF (where_2_GT_1 NE [-1]) THEN BEGIN
    rel_difference[where_2_GT_1] = $
      abs_difference[where_2_GT_1]/fielddata_2[where_2_GT_1]
  ENDIF 
   
  max_1=MAX(fielddata_1)
  max_2=MAX(fielddata_2)
  max_3=MAX(abs_difference)
  max_4=MAX(rel_difference)
  min_1=MIN(fielddata_1[WHERE(fielddata_1 GT 0)])
  min_2=MIN(fielddata_2[WHERE(fielddata_2 GT 0)])
  min_3=MIN(abs_difference[WHERE(abs_difference GT 0)])
  min_4=MIN(rel_difference[WHERE(rel_difference GT 0)])
  normalisation = MAX([max_1,max_2])
  lowerlimit = MIN([min_1,min_2])
  normalisation = 10.^(FLOOR(ALOG10(normalisation))+2)
  lowerlimit = 10.^(FLOOR(ALOG10(lowerlimit)))
  PRINT, "normalisation = ", normalisation
  PRINT, "lower limit = ", lowerlimit

  a = lowerlimit/normalisation
; DATA 1
  fielddata_1=fielddata_1/normalisation
  fielddata_1[WHERE(fielddata_1>0)]=(ALOG(fielddata_1[WHERE(fielddata_1>0)])-ALOG(a))/(-ALOG(a))

; DATA 2   
  fielddata_2=fielddata_2/normalisation
  fielddata_2[WHERE(fielddata_2>0)]=(ALOG(fielddata_2[WHERE(fielddata_2>0)])-ALOG(a))/(-ALOG(a))
  
; Absolute difference
  abs_difference=abs_difference/normalisation
  abs_difference[WHERE(abs_difference > 0)]=(ALOG(abs_difference[WHERE(abs_difference>0)])-ALOG(a))/(-ALOG(a))
  abs_difference[WHERE(abs_difference < 0)] = 0.0

; Relative difference
  rel_differemce=rel_difference
  rel_difference_save = rel_difference
  a = 0.01
  rel_difference[WHERE(rel_difference > 0)]=(ALOG(rel_difference[WHERE(rel_difference>0)])-ALOG(a))/(-ALOG(a))
  rel_difference[WHERE(rel_difference LE 0)]=0.0

;-----------------------------------------------------------------------  
; open postscript file for output

  SET_PLOT, 'PS'
  DEVICE, COLOR=1,/PORTRAIT, Decomposed=1, Bits_per_pixel=8 ;, /Encapsulated
  LOADCT, 5
  L = 20.0
  DEVICE,XSIZE=L,YSIZE=1.4*L
  DEVICE, FILENAME=filename

  
;-----------------------------------------------------------------------  
; plot: 
;        field_1             (upper left)
;        field_2             (upper right)
;        absolute difference (lower left)
;        relative difference (lower right)


  xwidth = 0.4
  ywidth = 0.4
  
  TV, (1.-fielddata_1)*255., /Centimeters,$
    0.05*L, 0.75*L, XSIZE=0.4*L, YSIZE=0.4*L
  xmin = 0.05
  ymin = 0.75
  ARROW, xmin, ymin/1.4, xmin+xwidth, ymin/1.4, $
    /Data, HSIZE=0.0, /Normalized 
  ARROW, xmin+xwidth, ymin/1.4, xmin+xwidth, (ymin+ywidth)/1.4, $
    /Data, HSIZE=0.0, /Normalized
  ARROW, xmin+xwidth, (ymin+ywidth)/1.4, xmin, (ymin+ywidth)/1.4, $
    /Data, HSIZE=0.0,/Normalized
  ARROW, xmin, (ymin+ywidth)/1.4, xmin, ymin/1.4, $
    /Data, HSIZE=0.0,/Normalized
  
  TV, (1.-fielddata_2)*255., /Centimeters, $
    0.55*L, 0.75*L, XSIZE=0.4*L, YSIZE=0.4*L
  xmin = 0.55
  ymin = 0.75
  ARROW, xmin, ymin/1.4, xmin+xwidth, ymin/1.4, $
    /Data, HSIZE=0.0, /Normalized 
  ARROW, xmin+xwidth, ymin/1.4, xmin+xwidth, (ymin+ywidth)/1.4, $
    /Data, HSIZE=0.0, /Normalized
  ARROW, xmin+xwidth, (ymin+ywidth)/1.4, xmin, (ymin+ywidth)/1.4, $
    /Data, HSIZE=0.0,/Normalized
  ARROW, xmin, (ymin+ywidth)/1.4, xmin, ymin/1.4, $
    /Data, HSIZE=0.0,/Normalized

    
  TV, (1.-abs_difference)*255., /Centimeters,$
    0.05*L, 0.05*L, XSIZE=0.4*L, YSIZE=0.4*L
  xmin = 0.05
  ymin = 0.05
  ARROW, xmin, ymin/1.4, xmin+xwidth, ymin/1.4, $
    /Data, HSIZE=0.0, /Normalized 
  ARROW, xmin+xwidth, ymin/1.4, xmin+xwidth, (ymin+ywidth)/1.4, $
    /Data, HSIZE=0.0, /Normalized
  ARROW, xmin+xwidth, (ymin+ywidth)/1.4, xmin, (ymin+ywidth)/1.4, $
    /Data, HSIZE=0.0,/Normalized
  ARROW, xmin, (ymin+ywidth)/1.4, xmin, ymin/1.4, $
    /Data, HSIZE=0.0,/Normalized

  TV, (1.-rel_difference)*255., /Centimeters, $
    0.55*L, 0.05*L, XSIZE=0.4*L, YSIZE=0.4*L
  xmin = 0.55
  ymin = 0.05
  ARROW, xmin, ymin/1.4, xmin+xwidth, ymin/1.4, $
    /Data, HSIZE=0.0, /Normalized 
  ARROW, xmin+xwidth, ymin/1.4, xmin+xwidth, (ymin+ywidth)/1.4, $
    /Data, HSIZE=0.0, /Normalized
  ARROW, xmin+xwidth, (ymin+ywidth)/1.4, xmin, (ymin+ywidth)/1.4, $
    /Data, HSIZE=0.0,/Normalized
  ARROW, xmin, (ymin+ywidth)/1.4, xmin, ymin/1.4, $
    /Data, HSIZE=0.0,/Normalized

;-----------------------------------------------------------------------  
; colour scale
  colour_bar = DBLARR(256,2)
  colour_bar[*,0] = INDGEN(256)
  colour_bar[*,1] = INDGEN(256)
  TV, 255 - colour_bar, /Centimeters, $
    0.05*L, 0.67*L, XSIZE=0.9*L, YSIZE=0.02*L

;-----------------------------------------------------------------------  
; Labels
  XYOUTS, 0.05, 0.97, "Filename: "+selgrid+"_"+timelabel+".txt"
  XYOUTS, 0.05, 0.95, "Directory 1 [field C^(1), top left]: "+datadir1, CHARSIZE=0.6
  XYOUTS, 0.05, 0.93, "Directory 2 [field C^(2), top right]: "+datadir2, CHARSIZE=0.6
  XYOUTS, 0.05, 0.90, "Species : "+selspecies
  XYOUTS, 0.05, 0.88, "Field   : "+selfield
  XYOUTS, 0.05, 0.86, "Name    : "+selname
  XYOUTS, 0.05, 0.84, "Level   : "+sellevel
  XYOUTS, 0.06, 0.51, "C^(1)[i,j]" 
  XYOUTS, 0.56, 0.51, "C^(2)[i,j]" 
  XYOUTS, 0.06, 0.35, "Absolute difference" 
  XYOUTS, 0.56, 0.35, 'Relative difference (in %)'
  XYOUTS, 0.02, 0.46, lowerlimit, CHARSIZE=0.8, ALIGNMENT=0.0
  XYOUTS, 0.02, 0.47, '1%',CHARSIZE=0.8, ALIGNMENT=0.0
  XYOUTS, 0.5, 0.46, "[logarithmic scale]", ALIGNMENT=0.5
  XYOUTS, 0.98, 0.46, normalisation, ALIGNMENT=1.0, CHARSIZE=0.8
  XYOUTS, 0.98, 0.47, '100%', ALIGNMENT=1.0, CHARSIZE=0.8
  XYOUTS, 0.06, 0.02, "|C^(1)[i,j]-C^(2)[i,j]|" 
  XYOUTS, 0.56, 0.02, "|C^(1)[i,j]-C^(2)[i,j]|/max_{k} {C^(k)[i,j]}"

;-----------------------------------------------------------------------  
; Calculate statistics

  mu = MEAN(signed_difference)
  mustr = "mu="+STRING(mu)
  XYOUTS, 0.1, 0.40, mustr

  sigma = STDDEV(signed_difference)
  sigmastr = "sigma="+STRING(sigma)
  XYOUTS, 0.4, 0.40, sigmastr

  g = SKEWNESS(signed_difference)
  gstr = "skewness="+STRING(g)
  XYOUTS, 0.7, 0.40, gstr
  
;-----------------------------------------------------------------------  
; Plot histogram of absolute difference

  nbins = 101
  delta = MAX([MAX(signed_difference),-MIN(signed_difference)])
  histo = HISTOGRAM(signed_difference, NBINS=nbins,MIN=-delta,max=delta)
  PLOT, (FINDGEN(nbins))/(nbins-1.)*(2*delta)-delta, histo, /NORMAL, POSITION=[0.1,0.1,0.9,0.9], PSYM=10, TITLE='Histogram of C^(1)-C^(2)'
 
;-----------------------------------------------------------------------  
; Plot cumulative histogram of relative difference
  
  histo = HISTOGRAM(rel_difference_save, NBINS=nbins,MIN=0,max=1)
  FOR i=1,nbins-1 DO BEGIN
    histo[i]=histo[i]+histo[i-1]
  ENDFOR
  histo = histo*1./histo[nbins-1]
  PLOT, (FINDGEN(nbins))/(nbins-1.), histo, /NORMAL, POSITION=[0.1,0.1,0.9,0.9], PSYM=10, TITLE='Cumulative histogram of the relative difference'

  DEVICE, /Close

end
