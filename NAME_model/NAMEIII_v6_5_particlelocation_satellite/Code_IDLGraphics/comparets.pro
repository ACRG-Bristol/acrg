pro comparets, $
  filename,    $
  datadir1,    $
  datadir2,    $
  location,    $
  field,       $
  level,       $
  species,     $
  grid

;--------------------------------------------------------------
; procedure to compare NAME III output time series which are
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

;  This script plots the two timeseries and a produces a scatter plot
;  with the first time series on the x-axis and the other on the 
;  y-axis. The correlation coefficient and the mean square error 
;  (normalised) by the number of data points are calculated and printed out.
;  The results are plotted to a postscript file.
;

;-----------------------------------------------------------------------  
; Read files

  tssplit = 'LOCATION'
 
  readts,      $
    datadir1,  $
    tssplit,   $
    location,  $
    field,     $
    level,     $
    species,   $
    grid,      $
    date,      $
    data1,     $
    units,     $
    namever=3
    
  readts,      $
    datadir2,  $
    tssplit,   $
    location,  $
    field,     $
    level,     $
    species,   $
    grid,      $
    date,      $
    data2,     $
    units,     $
    namever=3
    
;-----------------------------------------------------------------------  
; Open postscript file for output

  SET_PLOT, 'PS'
  DEVICE, COLOR=1,/PORTRAIT, Decomposed=0, Bits_per_pixel=8 ;, /Encapsulated
  TEK_COLOR
;  LOADCT, 1
;  L = 20.0
;  DEVICE,XSIZE=L,YSIZE=1.4*L
  DEVICE, FILENAME=filename
  !P.multi=[0,1,2]

;-----------------------------------------------------------------------  
; Plot time series

  plot, data1, thick=3, COLOR=0, POSITION=[0.05,-0.5,0.95,0.8], /NORMAL, $
    TITLE='Time series', $
    XTITLE='timestep', $
    YTITLE='Field value'
  oplot, data1, thick=3, COLOR=4
  oplot, data2, thick=3, COLOR=2
  maxC = MAX([MAX(data1),MAX(data2)])
   
;-----------------------------------------------------------------------  
; Labels

;  XYOUTS, 0.05, 0.97, "Filename: "+selgrid+"_"+timelabel+".txt"
  XYOUTS, 0.05, 1.20, "Directory 1 [TS 1]: "+datadir1, CHARSIZE=0.6, /NORMAL, $
    COLOR=4
  XYOUTS, 0.05, 1.17, "Directory 2 [TS 2]: "+datadir2, CHARSIZE=0.6, /NORMAL, $
    COLOR=2
  XYOUTS, 0.05, 1.10, "Location : "+location, /NORMAL
  XYOUTS, 0.05, 1.05, "Species  : "+species, /NORMAL
  XYOUTS, 0.05, 1.00, "Field    : "+field, /NORMAL
  XYOUTS, 0.05, 0.95, "Level    : "+level, /NORMAL

;-----------------------------------------------------------------------  
; Plot scatter plot and print out statistics

  plot, data1, data2, thick=2, COLOR=0, PSYM=1, $
    POSITION=[0.15,0.3,0.45,0.7], /NORMAL, CHARSIZE=0.5, $
    XTITLE="TS 1", $
    YTITLE="TS 2"
  XYOUTS, 0.05, -0.2, "Pearson Correlation Coefficient r = "+string(CORRELATE(data1, data2)), /NORMAL
  XYOUTS, 0.05, -0.3, "Average Mean Square Error = "+string(SQRT(TOTAL((data1-data2)^2))/N_ELEMENTS(data1)), /NORMAL
  DEVICE, /CLOSE
end
