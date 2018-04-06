pro readaqobs,obsdir,location,date,seldata,namespecies=namespecies,obsspecies=obsspecies,tmpdir=tmpdir

;--------------------------------------------------------------
; procedure to read AQ observation .csv files

; Lucy Davis 09/03/2010

; Arguments

; Required:
;  obsdir   : (string) directory containing obs files
;  location : (string) observation location eg 'ED3'
;  date     : (date array) time series date/time with required dates and times
;             date={datetime,year:0,month:0,day:0,hour:0,minute:0,second:0} 

; Output:
;  seldata  : (fltarr) observation time series data
;              Returns -999.9 where data not avaliable 

; Optional inputs - must have one of these two
;  namespecies : (string) species, in format from NAMEII style fields files
;  obsspecies  : (string) species, matching .csv header

;--------------------------------------------------------------

if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse

;-----------------------------------------------------------------------
; Setup output array

seldata=fltarr(n_elements(date))
seldata[*]=-999.9


;-----------------------------------------------------------------------
; filename

filename=strcompress(obsdir+delim+location,/remove_all)

spawn,'ls '+filename+'*.csv',filelist,exit_status=errcode


if errcode ne 0 then begin

  seldata[*]=-999.9
  
endif else begin  

  for ifile=0,n_elements(filelist)-1 do begin

    ; read title text

    ;Remove commas from certain species in file
    spawn,'head -n5 '+filelist(ifile)+' > '+tmpdir+'/sedtmp.txt'
    speciestoreplace=['1,2,3-trimethylbenzene','1,2,4-trimethylbenzene','1,3,5-trimethylbenzene','1,3-butadiene']
    replacewith=     ['1 2 3-trimethylbenzene','1 2 4-trimethylbenzene','1 3 5-trimethylbenzene','1 3-butadiene']
    FOR i_replace=0,n_elements(speciestoreplace)-1 DO BEGIN
      spawn,'sed -i "s%'+speciestoreplace(i_replace)+'%'+replacewith(i_replace)+'%I" '+tmpdir+'/sedtmp.txt'
    ENDFOR
    
    ;Header
    readdata,tmpdir+'/sedtmp.txt',header,nskip=4,nrows=1
  
    spawn,'rm '+tmpdir+'/sedtmp.txt'

    readdata,filelist(ifile),data,nrecs,nskip=5

    for idate=0,n_elements(date)-1 do begin

      ;Check dates format (yyyy-mm-dd or dd-mm-yyyy)
      datadate=data[0,0]
      datasplit=strsplit(datadate,'-',/extract)
      datasplitlen=strlen(datasplit)
      if datasplitlen[0] eq 4 then begin
        ; year first
        yearplace=0
        monthplace=5
        dayplace=8
      endif else if datasplitlen[2] eq 4 then begin
        ; year last
        yearplace=6
        monthplace=3
        dayplace=0
      endif else begin  
        print,'datasplitlen = ',datasplitlen
      endelse 
     
      dateindex=where(     ( strmid(data[0,*],yearplace,4) eq date(idate).year )  $
                       and ( strmid(data[0,*],monthplace,2) eq date(idate).month ) $
                       and ( strmid(data[0,*],dayplace,2) eq date(idate).day )   $
                       and ( strmid(data[1,*],0,2) eq date(idate).hour )  $
                        ,  datecount )
 
      if (date(idate).hour eq 0) and (idate gt 0) then begin
        ;Correction for observation data reporting hour 24:00:00 when 
        ; model reports 00:00:00
        dateindex=where(    ( strmid(data[0,*],yearplace,4) eq date(idate-1).year )  $
                        and ( strmid(data[0,*],monthplace,2) eq date(idate-1).month ) $
                        and ( strmid(data[0,*],dayplace,2) eq date(idate-1).day )   $
                        and ( strmid(data[1,*],0,2) eq date(idate-1).hour )  $
                       , datecount )
        dateindex=dateindex+1
      endif  
                
      if datecount gt 0 then begin
    
        ;Find species
      
        CASE namespecies OF
          'SULPHUR-DIOXIDE': species='Sulphur dioxide'
          'CO'             : species='Carbon monoxide'   ;1hr avg?
          'NO2'            : species='Nitrogen dioxide'
          'NO'             : species='Nitric oxide'
          'O3'             : species='Ozone'             ;1hr avg?
          'TOTAL-PM10'     : species='PM10 Particulate matter (Hourly)'
          'PM10'           : species='PM10 Particulate matter (Hourly)'
          ELSE             : species='' 
        ENDCASE                       

        if not keyword_set(obsspecies) then obsspecies=species

        if obsspecies eq '' then begin
          print,'species unknown'
          return
        endif          
      
        speciesindex=where( header eq obsspecies, speciescount)
        
        if speciescount eq 0 then begin
          speciesindex=where( header eq '"'+obsspecies+'"', speciescount)
        endif
        
        if speciescount eq 0 then begin
          if ( namespecies eq 'TOTAL-PM10' ) or ( namespecies eq 'PM10' ) then begin
            speciesindex=where( header eq '"PM10 particulate matter (Hourly measured)"', speciescount)
          endif
        endif
        
              
        if speciescount gt 0 then begin 
       
          units=data(speciesindex+2,0)
      
          CASE strmid(units,0,5) OF
            'ugm-3' : factor=1.0e-6
            'mgm-3' : factor=1.0e-3
            ELSE    : factor=1.0e-6
          ENDCASE
      
          concn=data[speciesindex,dateindex]
          IF strlen(concn) GT 0 THEN concn=concn*factor
          IF strlen(concn) EQ 0 THEN concn=-999.9
          seldata[idate]=concn
          
        endif else begin
          seldata[idate]=-999.9 
        endelse  
           
      endif ;datecount>0   
    
    endfor ;idate
  
  endfor  ;ifile

endelse ;file exists


end
