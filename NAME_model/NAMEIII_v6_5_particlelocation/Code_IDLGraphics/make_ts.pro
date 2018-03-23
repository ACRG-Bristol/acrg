pro make_ts,datadir,grid,$
            namever=namever,timeframe=timeframe,selcase=selcase


; read timeseries location file to generate time series files 
; compatible with RIMNET
; 
; assumes time series written with NAME namelist option 'TSSPLIT='LOCATION'
;

if(n_elements(namever) eq 0)then namever=2
if(n_elements(timeframe) eq 0)then timeframe='absolute'
if(n_elements(selcase) eq 0)then selcase='C1'

location=STRARR(100)
if (namever eq 3)then begin
  readdata,datadir+'/RIMNET_ts_locations.txt',location,nrecs, $
           nskip=2,nfields=4
endif else begin
  ;NAME
  readdata,datadir+'/tslocs.txt',location,nrecs, $
           nskip=0,nfields=4
endelse


print,'Locations:'
for i=0,n_elements(location)-1 do begin
 print,i,'   ',location(i)
endfor


tssplit='LOCATION'



; define species
print,'selgrid',grid

readtshead,datadir,modtitle,tshead,selgrid=grid,$
           namever=namever,timeframe=timeframe,selcase=selcase


if (namever eq 3)then begin
  array=tshead(where(tshead(*,2)),2)
  specieslist=array[uniq(array,sort(array))]
endif else begin
  ;NAME
  array=tshead(where(tshead(*,4)),4)
  specieslist=array[uniq(array,sort(array))]
endelse


numspecies=n_elements(specieslist)

print,'Species:'
for i=0,numspecies-1 do begin
  print,i,'   ',specieslist(i)
endfor



; loop through species 

for ispec=0,numspecies-1 do begin

  usespecies=specieslist(ispec)

  filename=usespecies
  if (usespecies eq 'IODINE-131') then begin
    filename='volatile'
  endif
  if (usespecies eq 'XENON-133') then begin
    filename='inert'
  endif
  if (usespecies eq 'CAESIUM-137') then begin
    filename='nonvol'
  endif
  if (usespecies eq 'PLUTONIUM-238') then begin
    filename='actinide'
  endif
  
  
  openw,11,datadir+'/'+filename+'.dat'

  for iloc=0,n_elements(location)-1 do begin

    uselocation=strcompress(location(iloc),/remove_all)
    monid=strcompress(strmid(location(iloc),0,4),/remove_all)
    sitename=strcompress(strmid(location(iloc),5,24),/remove_all)
    
    print,monid
    print,sitename
    print,usespecies


; read time series 

if (namever eq 3)then begin
    readts,datadir,tssplit,uselocation,$
     'Air Concentration',$
     'Boundary layer average',$
      usespecies,grid, $
      aircondates,aircon,airconunits,$
      namever=namever,timeframe=timeframe,selcase=selcase
endif else begin
    readts,datadir,tssplit,uselocation,$
     'Air concentration',$
     'Boundary layer',$
      usespecies,grid, $
      aircondates,aircon,airconunits,$
      namever=namever,timeframe=timeframe,selcase=selcase
endelse

if (namever eq 3)then begin
    readts,datadir,tssplit,uselocation,$
     'Dry deposition',$
     ' ',$
      usespecies,grid,drydepdates,drydep,drydepunits,$
      namever=namever,timeframe=timeframe,selcase=selcase
endif else begin
    readts,datadir,tssplit,uselocation,$
     'Dry deposition',$
     'Boundary layer',$
      usespecies,grid,drydepdates,drydep,drydepunits,$
      namever=namever,timeframe=timeframe,selcase=selcase
endelse
    
if (namever eq 3)then begin
    readts,datadir,tssplit,uselocation,$
     'Wet deposition',$
     ' ',$
      usespecies,grid,wetdepdates,wetdep,wetdepunits,$
      namever=namever,timeframe=timeframe,selcase=selcase
endif else begin
    readts,datadir,tssplit,uselocation,$
     'Wet deposition',$
     'Boundary layer',$
      usespecies,grid,wetdepdates,wetdep,wetdepunits,$
      namever=namever,timeframe=timeframe,selcase=selcase
endelse


; generate accumulated deposition

    drydep_acc=drydep
    wetdep_acc=wetdep
    for ih=1,n_elements(drydep)-1 do begin
      drydep_acc(ih)=drydep_acc(ih-1)+drydep(ih)
      wetdep_acc(ih)=wetdep_acc(ih-1)+wetdep(ih)
    endfor
    

; select hourly data
 
    index=where(aircondates.minute eq 0)
    aircondates=aircondates(index)
    aircon=aircon(index)
    drydep_acc=drydep_acc(index)
    wetdep_acc=wetdep_acc(index)
    
 
; start date
     
    sd=aircondates(0)
    

; first line is start (midnight)

    if(iloc eq 0)then begin
      printf,11,0,':',$
                0,' ',$
		sd.day,'/',$
		sd.month,'/',$
		sd.year,',',$
                format='(i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i4.4,a1)'
    endif
    
    
; start time is midnight so fill in untill start of run
    
    printf,11,' '
    printf,11,sitename,',',monid,',' 
   
    if(sd.hour gt 1)then begin
      for ih=1,sd.hour-1 do begin
        printf,11,ih,':',$
                   0,' ',$
		   sd.day,'/',$
		   sd.month,'/',$
		   sd.year,', ',$
		   0.0,',',$
		   0.0,',',$
		   0.0,',',$
         format='(i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i4.4,a2,e8.2,a1,e8.2,a1,e8.2,a1)'
      endfor
    endif
   
   
    for it=0,n_elements(aircondates)-1 do begin
       if(aircondates(it).minute eq 0)then begin
       
       printf,11,aircondates(it).hour,':',$
                 aircondates(it).minute,' ',$
		 aircondates(it).day,'/',$
		 aircondates(it).month,'/',$
		 aircondates(it).year,', ',$
		 aircon(it),',',$
		 drydep_acc(it),',',$
		 wetdep_acc(it),',',$
         format='(i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i4.4,a2,e8.2,a1,e8.2,a1,e8.2,a1)'
      endif
    endfor  

  endfor 

  printf,11,' '
  close,11
endfor

end


