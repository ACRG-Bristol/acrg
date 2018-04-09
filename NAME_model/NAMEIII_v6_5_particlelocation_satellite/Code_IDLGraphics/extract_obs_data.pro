pro extract_obs_data,obs_datadir,year,$
  species,obs_loc,obs_date,obs_data,obs_units,end_date
;-------------------------------------------------------------------------------------------
; procedure to read the cleaned up observed data 

; output:
; obs_date			  : (date array) time series date/time
; obs_data   			  : (string)


; Version number	Date		Owner			Comment
; 0.1			Nov 2003	Mark Harrison

; ------------------------------------------------------------------------------------------------

print,obs_loc
print,year

;defining file path
filename=obs_datadir+obs_loc+year+'.txt'

filename= strcompress(filename,/remove_all)
print,'obs filename'
print,filename

; create a string array into which the header text describing the units 
; for each species is held (obs_header), read the data in
    
 obs_header=strarr(18)
       
 status=dc_read_free(filename,obs_header,/row,nrecs=1,nskip=1)
 
   if(status ne 0)then begin
     print,'Cannot open file:',filename
     missing=1
   endif  

;read the data in now
   
;setting up the array into which the data are to be read 
;species_obs are by default set to floating point values

 carbon_monoxide_units=strarr(1)
 nitric_oxide_units=strarr(1)
 nitrogen_dioxide_units=strarr(1)
 nitrogen_oxides_units=strarr(1)
 ozone_units=strarr(1)
 pm10_units=strarr(1)
 pm2p5_units=strarr(1)
 sulphur_dioxide_units=strarr(1)


; IF N_ELEMENTS(obs_date) EQ 0 THEN obs_date = {!dt}

 obs_date=replicate({!dt},100000)

; read in the observed data calling it *_obs
 
 status=dc_read_fixed(filename,obs_date,obs_date,carbon_monoxide_obs,$
   carbon_monoxide_units,nitric_oxide_obs,nitric_oxide_units,$
   nitrogen_dioxide_obs,nitrogen_dioxide_units,nitrogen_oxides_obs,$
   nitrogen_oxides_units,ozone_obs,ozone_units,$
   pm10_obs,pm10_units,pm2p5_obs,pm2p5_units,sulphur_dioxide_obs,$
   sulphur_dioxide_units,/column,nskip=2,resize=[indgen(18)+1],$
   Format="(100(a25))",dt_template=[5,-1])

 missing=0
 if(status ne 0)then begin
     print,'Cannot open file:',filename
     missing=1
 endif 

; for the species that is requested call the *_obs data the obs_data
; if a species has been requested that is not in the array then print 
; out a message listing the species name

print,'Missing:',missing

if missing eq 0 then begin
  CASE species OF
  
    'CO': begin
    	obs_data=carbon_monoxide_obs
        obs_units=carbon_monoxide_units
        end  
  'NO':	begin
  	obs_data=nitric_oxide_obs    
        obs_units=nitric_oxide_units
        end
  'NO2': begin
  	obs_data=nitrogen_dioxide_obs    
        obs_units=nitrogen_dioxide_units
        end
  'O3':	begin
  	obs_data=ozone_obs      
        obs_units=ozone_units
        end
  'PM10_primary': begin
  	obs_data=pm10_obs 
        obs_units=pm10_units
        end
  'PM10_total':	begin
        obs_data=pm10_obs 
        obs_units=pm10_units
        end
  'SO2': begin
  	obs_data=sulphur_dioxide_obs 
        obs_units=sulphur_dioxide_units
        end
  
  
   
  ELSE: begin
      print,'Species:'+species+' Not Defined, Fill with -9.9'  
      obs_data=[-9.9]
      obs_date=[today()]
      obs_units=[' ']
    end  
  ENDCASE
endif

if missing ne 0 then begin
  print,'Species:'+species+' Not Defined, Fill with -9.9'  
  obs_data=[-9.9]
  obs_date=[today()]
  obs_units=[' ']
endif  

End


