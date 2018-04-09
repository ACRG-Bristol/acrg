PRO tsstats,moddata,obsdata,mdi=mdi,correlation=correlation,bias=bias,nmb=nmb,$
            mge=mge,nmge=nmge,rmse=rmse,fac2=fac2,$
            threshold=threshold,orss=orss,hitrate=hitrate,$
            falsealarmrate=falsealarmrate
            
;Calculate statistics from timeseries data
; Lucy Davis Aug 2010

;Required inputs:
; moddata - array of model data - missing data < 0
; obsdata - array of observation data - missing data < 0

; Optional inputs:
; threshold required for orss, hitrate and falsealarmrate            
; mdi - missing data indicator - value given to any missing data.
;     - if not given, assumes missing data is negative.

;Outputs
; Correlation - r
;        - -1 <= r <= +1
;        - 1 => +ve correlation (good!)
; Bias - mean bias
;        - Same units as obs
;        - 0 => No bias (good!)
; nmb  - Normalised mean bias
;        - -1 <= nmb <= 1.0
;        - 0 => No bias (good!)
;        - 'Acceptable' if -0.2 < nmb < 0.2
; mge  - Mean gross error
;        - Same units as obs
;        - 0 <= mge
;        - 0 => no error (good!)
; nmge - Normalised mean gross error
;        - 0 <= nmge <= 1
;        - 0 => no error (good!)
; rmse - Root mean square error
;        - 0 <= rmse
;        - 0 => no error (good!)
; fac2 - Fraction of model within factor 2 of observations
;        - 0 <= fac2 <= 1.0
;        - 'Acceptable' if fac2 > 0.5
; orss - Odds ratio skills score
;        - -1 <= ORSS <= +1
;        - -1 => Strong -ve association with obs (bad!)
;        -  0 => Random forecast
;        - +1 => Strong +ve association with obs (good!)
; Hitrate
;      - Odds of a success 
;      - probability that an event is forecast given that is observed.
;        - 0 <= hitrate <= 1
;        - 1 => event always forecast when observed (good!)
;        - 0 => event never forecast when observed (bad!)
; Falsealarmrate
;      - Odds of a failure
;      - probability that an event is forecast given that is was not observed.
;        - 0 <= falsealarmrate <= 1
;        - 1 => event is always forecast when it is not observed (bad!)
;        - 0 => event is never forecast when it is not observed (good!)

;-------------------------------------------------------------------------------

if (n_elements(mdi) ne 0) then begin
  index=where((obsdata ne mdi ) and (moddata ne mdi), n)
endif else begin
  index=where((obsdata ge 0) and (moddata ge 0), n)
endelse

if n gt 0 then begin

  model=float(moddata(index))
  obs=float(obsdata(index))
  ;print,model
  ;print,obs

  ;r
  correlation=CORRELATE(model,obs)
  ;print,'correlation:',correlation

  ;Bias
  total=0.0
  total=total(model-obs)
  bias=total/float(n)  
  ;print,'bias:',bias  
  ;Normalised mean bias
  nmb=total/total(obs)

  ;Mean Gross Error
  total=0.0
  total=total( abs (model-obs) )
  mge=total/float(n)
  ;Normalised mean gross error
  nmge=total/total(obs)

  ;RMSE
  total=total((model-obs)^2)
  rmse=sqrt(total/float(n))  
  ;print,'RMSE:',rmse  
  
  ;FAC2
  factor=model/obs
  fac2indices=where( (factor ge 0.5) and (factor le 2.0) , fac2count)
  fac2=float(fac2count)/float(n)
  ;print,'fac2:',fac2

  if keyword_set(threshold) then begin
    ;print,'threshold: ',threshold
  
    thresindex=where(obs ge threshold and model ge threshold, a)
    thresindex=where(obs lt threshold and model ge threshold, b)
    thresindex=where(obs ge threshold and model lt threshold, c)
    thresindex=where(obs lt threshold and model lt threshold, d)
  
    theta=(float(a)*float(b))/(float(c)*float(d))
    ORSS=(theta-1)/(theta+1)
  
    hitrate= float(a) / (float(a)+float(c))
  
    falsealarmrate= float(b) / (float(b)+float(d))
  
  ;  print,'ORSS: ',ORSS
  ;  print,'hitrate:',hitrate
  ;  print,'false alarm rate:', falsealarmrate
  
  endif else begin
    ORSS=-999
    hitrate=-999  
    falsealarmrate=-999
  endelse  

endif else begin
  correlation=-999
  bias=-999
  nmb=-999
  mge=-999
  nmge=-999
  rmse=-999
  fac2=-999
  ORSS=-999
  hitrate=-999  
  falsealarmrate=-999
endelse

END
