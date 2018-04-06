pro convert_mod_data,mod_datadir,field,level,grid,$
  species,mod_loc,mod_date,mod_data,mod_units 

;-------------------------------------------------------------------------------------------
; procedure to read the modelled data acquired using NAME

; Required information         : temp
;                              : press
;                              : mwt for each species of interest
;			       : air density

; calls the routine            : readts.pro

;output:
;mod_date                      : (date array) time series date/time
;mod_data                      : (date array) time series date/time
;mod_units                     : (date array) time series date/time

;Version number		Date		Owner			Comment
;0.1			Nov 2003	Mark Harrison

; Contains sub-procedures :  gramtoppb, gramtoppm, writeheader


;-------------------------------------------------------------------------------------------
air_den=2.55e19 
temp=283.15
press=1000.0

mwt_co=  12. + 16.
mwt_no=  14. + 16.
mwt_no2= 14. + (2.*16.)
mwt_o3= (3. * 16.)
mwt_so2= 32. + (2. * 16.)


; If the species is that requested in the .com file,read the relevant information
; using the readts routine.Then calls the sub-routine at the base of this file 
; to convert to the relevant units

If species eq 'CO' then begin
   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
     'CO','grid1',mod_date,mod_data,mod_units 
   gramtoppm,mod_data,temp,press,mwt_co,ppm_co
   mod_data=ppm_co
   mod_units='/ ppm'	  

EndIf else If species eq 'NO' then begin	  
   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
    'NO','grid1',mod_date,mod_data,mod_units
   gramtoppb,mod_data,temp,press,mwt_no,ppb_no 
   mod_data=ppb_no  
   mod_units='/ ppb' 

EndIf else If species eq 'NO2' then begin	  
   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
     'NO2','grid1',mod_date,mod_data,mod_units
   gramtoppb,mod_data,temp,press,mwt_no2,ppb_no2
   mod_data=ppb_no2 
   mod_units='/ ppb'   

EndIf else If species eq 'O3' then begin	  
;Values for mod_O3 are calculated as the sum of  O3_PART and O3
;as described in NAME

   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
     'O3','grid1',mod_date,mod_data,mod_units
   
   ;convert from VMR units to ppb
   mod_data=mod_data*air_den*48.0*1e12/(2*6.02e23)
   mod_units='/ ppb'
               
EndIf else If species eq 'PM10_total' then begin	  
;Values for PM10_total are calculated as the sum of primary PM10, NH42SO4,  
;SULPHATE, NH4NO3, NAER as described in NAME
 
   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
     'PM10','grid1',mod_date,mod_PM10,mod_units
     
   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
     'NH42SO4','grid1',mod_date,mod_NH42SO4,mod_units
      
   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
     'SULPHATE','grid1',mod_date,mod_SULPHATE,mod_units
       
   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
     'NH4NO3','grid1',mod_date,mod_NH4NO3,mod_units
       
   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
     'NAER','grid1',mod_date,mod_NAER,mod_units 
       
   mod_data=(mod_PM10+mod_NH42SO4+mod_SULPHATE+mod_NH4NO3+mod_NAER)*1.0e6 
   mod_units='/ ugm-3'

EndIf else If species eq 'PM10_primary' then begin	  
;Values for PM10_primary are taken as primary PM10 as described in NAME
   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
     'PM10','grid1',mod_date,mod_data,mod_units
   mod_data=(mod_data)*1.0e6  
   mod_units='/ ugm-3'

EndIf else If species eq 'SO2' then begin	  		  
   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
     'SULPHUR-DIOXIDE','grid1',mod_date,mod_data,mod_units
   gramtoppb,mod_data,temp,press,mwt_so2,ppb_so2
   mod_data=ppb_so2  
   mod_units='/ ppb'   

EndIf else begin
   readts,mod_datadir,'LOCATION',mod_loc,field,level,$
     species,'grid1',mod_date,mod_data,mod_units
   mod_data=mod_data*1.0e6
   mod_units='/ ug/m3'   

endelse

END

;*************************************************************************
pro gramtoppb,gm3,temp,press,mwt,ppb

; this subroutine converts the input concentration in g/m3 to a
; mixing ratio in ppb

 ppb = (1.E6*gm3)*(22.41*temp*1013.25/(273.15*press))/mwt
 END 

;*************************************************************************
pro gramtoppm,gm3,temp,press,mwt,ppm

; this subroutine converts the input concentration in g/m3 to a
; mixing ratio in ppm

 ppm = (1.E3*gm3)*(22.41*temp*1013.25/(273.15*press))/mwt
 END 

;*************************************************************************
