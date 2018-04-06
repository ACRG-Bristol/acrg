pro calc_stats,mod_data,obs_data,start_date,end_date,mod_date,$
  obs_date,r,rmse,fac2,fb,nmse,stderr,meanerr,stats_filename,location,species

;-------------------------------------------------------------------------
;procedure to generate statistics for all locations + species
;rocedure called from within the plot_obs_vs_mod.pro routine


;Version number		Date		Owner			Comment
;0.1			Dec 2003	Mark Harrison


; INPUT
; obs (t): timeseries of observations
; obsnoise : level of observational noise for FAC2 calculation
; model (t,x,y): timeseries of scaled model data at MH per grid point
; obsnoise initially set the noise to zero

; Output generated for each species and locations:
; fac2,fb,r,rmse,nmse,stderr,meanerr

;------------------------------------------------------------------------

 print,'first mod_date',mod_date(0)
 


print,'end_date',end_date
stat_start=start_date
stat_end=end_date
print,'line 27 stat_end',stat_end


n_obs=n_elements(obs_data)
print,'n_obs',n_obs

n_mod=n_elements(mod_data)
print,'n_mod',n_mod



;If the start_date specified in the .com is before the available 
;mod_data, then only do the stats over the time period for 
;which mod_data is available

if (stat_start.julian lt mod_date(0).julian) then begin
  stat_start=mod_date(0)
endif
  
;If the start_date specified in the .com is before the available 
;obs_data, then only do the stats over the time period for 
;which obs_data is available     
 
if (stat_start.julian lt obs_date(0).julian) then begin
  stat_start=obs_date(0)
endif
  
;If the start_date specified in the .com is after the available 
;mod_data, then only do the stats over the time period for 
;which mod_data is available 
 
 print,'last mod_date',mod_date(n_mod-1)
 
if (stat_end.julian gt mod_date(n_mod-1).julian) then begin
  stat_end=mod_date(n_mod-1)
endif

  
;If the start_date specified in the .com is after the available 
;obs_data, then only do the stats over the time period for 
;which obs_data is available 
  
if (stat_end.julian gt obs_date(n_obs-1).julian) then begin
  stat_end=obs_date(n_obs-1)
endif

;If the stat_end calculated above is larger than stat_end calculated above 
;then aviod doing the stats by setting stat_start,mod_date and obs_date 
;to stat_end to set mod_data and obs_data to -9.9
 
 
print,'stat_start', stat_start
print,'stat_end', stat_end  
  
if (stat_start.julian ge stat_end.julian) then begin

  print,'No Stats'
 
endif else begin

dt_print,stat_start
dt_print,stat_end
dt_print,mod_date(0)
dt_print,mod_date(n_elements(mod_date)-1)
dt_print,obs_date(0)
dt_print,obs_date(n_elements(obs_date)-1)




;Acquiring the relevant paired data over the period stat_start till stat_end
;This also avoids data less than 0

ind=where(mod_date.julian ge stat_start.julian, nind)
startindex=0
if nind gt 0 then startindex=ind(0)
ind=where(mod_date.julian ge stat_end.julian, nind)

if nind gt 1 then begin
  endindex=ind(0)-1
  endindex=max([0,endindex])
  print,'df1'
endif else if nind eq 1 then begin
  endindex=ind(0)  
endif else begin
  endindex=n_elements(mod_date)-1
  endindex=max([0,endindex])
  print,'df2'
endelse

print,'MOD:',startindex,endindex,n_elements(mod_date)

temp_date=mod_date(startindex:endindex)
temp_mod=mod_data(startindex:endindex)

startindex_obs=0
ind=where(obs_date.julian ge stat_start.julian, nind)
if nind gt 0 then startindex_obs=ind(0)
ind=where(obs_date.julian ge stat_end.julian, nind)
if nind gt 1 then begin
  endindex_obs=ind(0)-1
endif else if nind eq 1 then begin
  endindex_obs=ind(0)  
endif else begin
  endindex_obs=n_elements(obs_date)-1
endelse

print,'OBS:',startindex_obs,endindex_obs,n_elements(obs_date)
temp_obs=obs_data(startindex_obs:endindex_obs)

info,temp_date,temp_mod,temp_obs

; Remove negative values

nostats=0
ind=where(temp_obs ge 0.0,nind)
if nind gt 0 then begin
  temp_date=temp_date(ind)
  temp_obs=temp_obs(ind)
  temp_mod=temp_mod(ind)
endif else begin
  nostats=1
endelse

ind=where(temp_mod ge 0.0,nind)
if nind gt 0 then begin
  temp_date=temp_date(ind)
  temp_obs=temp_obs(ind)
  temp_mod=temp_mod(ind)
endif else begin
  nostats=1
endelse

if nostats eq 0 then begin
  n_temp=n_elements(temp_date)

  print,'n_temp',n_temp

 ;starting to do the actual statistics
 ;obsnoise initially set the noise to zero

  obsnoise=0
  f_obs=float(n_temp)
  av_obs = AVG(temp_obs)
  err=replicate(0.0,n_temp)
  modelvalue=replicate(0.0,n_temp)

  fac2=0.0
  for t=0,n_temp-1 do begin

    if (temp_obs(t) gt obsnoise) then begin
  
      obs_dble=temp_obs(t)*2.0
      obs_half=temp_obs(t)/2.0
           
    endif else begin
      obs_dble=obsnoise*2.0
      obs_half=0.0
     
    endelse

    if ((temp_mod(t) ge obs_half) and $
        (temp_mod(t) le obs_dble)) then fac2=fac2+1.0    

  endfor

  mean_model = AVG(temp_mod)
  print,'mean_model',mean_model 

  err = (temp_mod-temp_obs)

  if n_elements(err) gt 1 then begin
    stderr=STDEV(err,meanerr)
  endif else begin
    stderr=-9.9
  endelse

  print,'stderr',stderr 

  err = err^2.0

  toterr=total(err)
  print,'toterr',toterr
 
  mse = toterr/f_obs
  print,'mse',mse 

  rmse=mse^0.5
  print,'rmse',rmse 

  nmse=mse/(mean_model*av_obs)
  print,'nmse',nmse 
  
  r=correlate(temp_obs,temp_mod)
  print,'r',r 
 
  fac2=fac2/f_obs
  print,'fac2',fac2

  fb=2*(mean_model-av_obs)/(mean_model+av_obs)
  print,'fb',fb
 
endif else begin
  fac2=-9.9
  fb=-9.9
  r=-9.9
  rmse=-9.9
  nmse=-9.9
  stderr=-9.9
  meanerr=-9.9
  n_temp=-9.9
endelse  


 get_lun,fout
 openw,fout,stats_filename,/append

 printf,fout,location,species,fac2,fb,r,rmse,nmse,stderr,meanerr,$
   Format="(A30,A20,7F10.4)"

 close,fout
 free_lun,fout

; Stats required
endelse

end
