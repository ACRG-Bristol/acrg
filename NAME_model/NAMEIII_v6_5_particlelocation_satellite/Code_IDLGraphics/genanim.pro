pro genanim,allfiles,datadir,filetem,gengif=gengif,genjpg=genjpg,$
 genanim=genanim,genpng=genpng,method=method,xpixel=xpixel

;-----------------------------------------------------------------------
; procedure to generate gif animation from postscript files
; DBR March 2002
;
;arguments
; allfiles      : string array  list of postscript files
; datadir       : string        data directory
; filetem       : string        part of filename common to all files
; gengif        : integer (0 or 1) =1 to generate gifs
; genjpg        : integer (0 or 1) =1 to generate jpg
; genanim       : integer (0 or 1) =1 to generate gif animations
; genpng        : integer (0 or 1) =1 to generate png
; method        : string ('alchemy' or 'convert') coomand to 
;                  generate gif/jpg images

;-----------------------------------------------------------------------


; arguments

if(n_elements(gengif) eq 0)then gengif=0
if(n_elements(genanim) eq 0)then gengifanim=1
if(n_elements(xpixel) eq 0)then xpixel=595
if(n_elements(genjpg) eq 0)then genjpg=0
if(n_elements(genpng) eq 0)then genpng=0
if(n_elements(method) eq 0)then method='convert'

print,'Using '+method

if(method eq 'convert')then begin
 density=fix(xpixel/8.2639+0.5)
 cdensity=' -density '+string(density,format='(i3.3)')+$
  'x'+string(density,format='(i3.3)')+' '
endif

if(method eq 'pstoimg')then begin
 density=fix(xpixel/8.2639+0.5)
 cdensity=' -density '+string(density,format='(i3.3)')+' '
endif

if(n_elements(allfiles) eq 1)then begin
  singleps=1 
endif else begin
  singleps=0
endelse


; delay between frames

if(n_elements(allfiles) gt 12)then begin
  delay='10'
endif else begin
  delay='30'
endelse


; generate gifs

if( gengif eq 1 or genanim eq 1)then begin

  for ifile=0,n_elements(allfiles)-1 do begin

     ; First determine number of pages in file
     spawn, 'tail -n 3 '+ allfiles(ifile) + ' | grep Pages: ', PageCount
     pspages = strsplit(PageCount,'%%Pages:',/extract)
     npages = fix(pspages[0])

     for page=1,npages(0) do begin
        tempfile = datadir+'/tempplotfile.ps'
        spawn, 'psselect -p'+strtrim(string(page),1)+' '+allfiles(ifile)+' '+tempfile
     

        case method of
           'convert' : begin

              posdot=strpos(allfiles(ifile),'.',/reverse_search)
              newname=strmid(allfiles(ifile),0,posdot)+string(page-1,format='(i03)')+'.gif'

              
              spawn,' convert -page 595x842! -background white -flatten +adjoin '+ cdensity +'"'+$
                    tempfile +'"'+ ' '+'"'+ newname +'"'
   
           end

           'pstoimg': begin

              posdot=strpos(allfiles(ifile),'.',/reverse_search)
              newname=strmid(allfiles(ifile),0,posdot)+string(page-1,format='(i03)')+'.gif'
              
              spawn, 'pstoimg -quiet -type gif -aaliastext -antialias '+ cdensity +$
                     '"'+ tempfile +'" -out "'+ newname +'"'
 
           end
           
        endcase
     endfor
  endfor
endif


; generate jpgs

if(genjpg eq 1)then begin

  for ifile=0,n_elements(allfiles)-1 do begin

     ; First determine number of pages in file
     spawn, 'tail -n 3 '+ allfiles(ifile) + ' | grep Pages: ', PageCount
     pspages = strsplit(PageCount,'%%Pages:',/extract)
     npages = fix(pspages[0])

     for page=1,npages(0) do begin
        tempfile = datadir+'/tempplotfile.ps'
        spawn, 'psselect -p'+strtrim(string(page),1)+' '+allfiles(ifile)+' '+tempfile

        posdot=strpos(allfiles(ifile),'.',/reverse_search)
        newname=strmid(allfiles(ifile),0,posdot)+string(page-1,format='(i03)')+'.jpg'

        spawn,' convert -page 595x842! +adjoin -background white -flatten -quality 100 '+ $
              cdensity +'"'+ tempfile +'"'+ ' ' +'"'+ newname +'"'

     endfor
  endfor
endif


; generate pngs

if(genpng eq 1)then begin

  for ifile=0,n_elements(allfiles)-1 do begin

     ; First determine number of pages in file
     spawn, 'tail -n 3 '+ allfiles(ifile) + ' | grep Pages: ', PageCount
     pspages = strsplit(PageCount,'%%Pages:',/extract)
     npages = fix(pspages[0])

     for page=1,npages(0) do begin
        tempfile = datadir+'/tempplotfile.ps'
        spawn, 'psselect -p'+strtrim(string(page),1)+' '+allfiles(ifile)+' '+tempfile
  
        case method of
           'convert' : begin

              posdot=strpos(allfiles(ifile),'.',/reverse_search)
              newname=strmid(allfiles(ifile),0,posdot)+string(page-1,format='(i03)')+'.png'
              
              spawn,' convert -page 595x842! +adjoin -background white -flatten '+ cdensity +'"'+$
                    tempfile +'"'+ ' '+'"'+ newname +'"'
    
           end

           'pstoimg': begin

              posdot=strpos(allfiles(ifile),'.',/reverse_search)
              newname=strmid(allfiles(ifile),0,posdot)+string(page-1,format='(i03)')+'.png'
              
              spawn, 'pstoimg -quiet -type gif -aaliastext -antialias '+ cdensity +$
                     '"'+ tempfile +'" -out "'+ newname +'"'
   
           end
        endcase

     endfor
    
  endfor
endif


; generate gif animation


if(genanim eq 1)then begin

   moviefiles=datadir+'/'+filetem+'*.gif'       

   if(method eq 'convert')then begin
      spawn,'convert -loop 0 -delay ' + delay + ' ' +'"'+ moviefiles +'"'+ ' ' + datadir + '/Anim_'+'"'+filetem+'"'+'.gif'
   endif

endif


; remove individual gifs if required

if(gengif eq 0 and genanim eq 1)then begin
   spawn,'rm -f '+datadir+'/'+'"'+filetem+'"'+'*.gif'   
endif 

; Remove temporary postscript file
spawn, 'rm -f '+tempfile

   
end



