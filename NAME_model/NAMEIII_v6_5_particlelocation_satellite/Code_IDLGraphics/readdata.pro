PRO readdata,filename,data,nrecs,nskip=nskip,nfields=nfields,$
    nrows=nrows,delim=delim,null=null

;-------------------------------------------------------------------------------
; Procedure to read comma-separated data from file.
; Reads in each row as a single string, then splits it into nfields columns
;
; Lois Huggett 10/02/2009
;------------------------------------------------------------------------------- 

;inputs:
;  filename  : 
;  nskip     : number of lines to be skipped at start of file
;  nfields   : number of columns in each line (add 1 if comma after last datum in line)
;  nrows     : (optional) number of rows to be read in
;  delim     : (optional) delimiter in file, comma as default
;  null      : (optional) =1 (default) if null elements to be retained

;outputs:
;  data      : (string) array of file records
;  nrecs     : (optional) number of records in file - returned if nrows not specified

;------------------------------------------------------------------------------- 

if n_elements(delim) eq 0 then delim=','
if n_elements(null) eq 0 then null=1

;------------------------------------------------------------------------------- 

status = QUERY_ASCII(filename, info)
if(status ne 1) then begin
  print,'2. Cannot open file:',filename
  data=0
  return
endif

restore, 'fieldtemplate.sav'
; fieldtemplate reads 1 string per row and skips 0 rows automatically
template=fieldtemplate

if (n_elements(nrows) eq 0) then begin
  dataread=read_ascii(filename,record_start=nskip,template=template,count=nrecs)
endif else begin
  nrecs=nrows
  dataread=read_ascii(filename,num_records=nrecs,record_start=nskip,template=template)
endelse

if (n_elements(nfields) eq 0) then begin
  temp=strsplit(dataread.field1(0),delim,/Extract,Preserve_Null=null)
  nfields=n_elements(temp)
endif

if (nrecs eq 0) then begin
  data=strarr(nfields,1)
  data(0,0)='1'
  data(1,0)='1'
  
endif else begin
  data=strarr(nfields,nrecs)

  for i=long(0),nrecs-1 do begin
    data(*,i)=strsplit(dataread.field1(i),delim,/Extract,Preserve_Null=null)
  endfor

endelse

END
