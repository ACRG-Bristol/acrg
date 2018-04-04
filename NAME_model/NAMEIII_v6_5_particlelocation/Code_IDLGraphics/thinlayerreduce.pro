pro thinlayerreduce,                              $
    datadir,seltime,selgridin,selgridout,         $
    selfield,selspecies,                          $
    thinlayerlist,thicklayerlist,reductionmethod, $
    factor=factor,namever=namever,                $
    boundarythreshold=boundarythreshold

;--------------------------------------------------------------
;
; procedure to reduce sets of thin vertical layers into
; thick layers. 
; [Works only for input files in NAME II file format]
;
; Eike Mueller, 21/01/2011
;
;--------------------------------------------------------------

; inputs:
;
;  ....filename specifiers
;  datadir            : directory with in/out files
;  seltime            : } used to construct filenames in readfield.pro
;  selgridin          : } in=datadir+'/Fields_'+selgridin+'_'+seltime+'.txt'
;  selgridout         : } out=datadir+'/Fields_'+selgridout+'_'+seltime+'.txt'
;
;  ....parameters specifying the fields to process in the input file
;  selfield           : name of field (e.g. 'Air Concentration')
;  selspecies         : name of species (e.g. 'VOLCANIC_ASH')
;
;  ....list defining layer averages (see discussion below)
;  thinlayerlist      : 2d list of sets of thin layers which are to be
;                       reduced to thick layers
;  thicklayerlist     : list of thick layers, into which the thin layers
;                       are to be reduced
;
;  ....reduction method
;  reductionmethod    : defines how layers are combined. Can be 
;                       'average' or 'peak'
;
;  ....optional parameters
;  factor             : if reductionmethod is 'peak', thick layer 
;                       concentrations are multiplied by this 
;                       multiplicative factor, which defaults to 1.0
;  namever            : NAME file version, only 2 is currently supported, 
;                       defaults to 2 if no value is given
;
;
;  description:
;
;  Each element in thinlayerlist is a list of layers that are to be
;  reduced to a single thick layer (Use '' to denote missing layers). 
;  The number of elements in thicklayerlist has to be equal to the number of 
;  thin layer lists, e.g. one can have:
;
;  thinlayerlist = [ ['From FL000 - FL025','From FL025 - FL050'],
;                    ['From FL050 - FL075','From FL075 - FL100'],
;                    ['From FL100 - FL125',''                  ] ]
;  thicklayerlist = ['From FL000 - FL050',
;                    'From FL050 - FL100',
;                    'From FL100 - FL125']
;
;  for the following reduction:
;
;  'From FL000 - FL025' } 
;  'From FL025 - FL050' } => 'From FL025 - FL050'
;
;  'From FL050 - FL075' }
;  'From FL075 - FL100' } => 'From FL050 - FL100'
;
;  'From FL100 - FL125' } => 'From FL100 - FL125'
;
;  To average over all thin layers in a thick layer set 
;  reductionmethod='average'. Note that all layers are given the same
;  weight, i.e. this only gives sensible output if all thin layers have the
;  same thickness.
;
;  To calculate peaks over thin layers, i.e. find the maximal 
;  concentration value in each thick layer, set reductionmethod='peak'.
;  Optionally peak values can be multiplied with an additional factor
;  given by the keyword 'factor'.
;
;  Example: Consider a file Fields_grid1_201101191800.txt (in directory
;  <datadir>) with the following column headers:
;
; ...             VOLCANIC,                VOLCANIC,                VOLCANIC,                VOLCANIC,  
; ...         VOLCANIC_ASH,            VOLCANIC_ASH,            VOLCANIC_ASH,            VOLCANIC_ASH,            
; ... 006 hr time averaged,    006 hr time averaged,    006 hr time averaged,    006 hr time averaged,    
; ...    Air Concentration,       Air Concentration,       Air Concentration,       Air Concentration,       
; ...                 g/m3,                    g/m3,                    g/m3,                    g/m3,                    
; ...   From FL000 - FL025,      From FL025 - FL050,      From FL050 - FL075,      From FL075 - FL100,      
; ...   0600UTC 20/01/2011,      0600UTC 20/01/2011,      0600UTC 20/01/2011,      0600UTC 20/01/2011,      
; 
;  To calculate peaks over 50FL thick layers, multiply by 
;  the factor 2.0, and write output to file 
;  Fields_grid2_201101191800.txt use:
;
; thinlayerpeaks,,                                            \
;    <datadir>,'201101191800','grid1','grid2',                \
;  'Air Concentration','VOLCANIC_ASH','006 hr time averaged', \
;  [['From FL000 - FL025','From FL025 - FL050'],              \
;   ['From FL050 - FL075','From FL075 - FL100']],             \
;  ['From FL000 - FL050','From FL050 - FL100'],               \
;   'peak',factor=2.0   

;--------------------------------------------------------------

; Parameter checking
; Check for NAME version. Only files in NAME II format can be
; processed
if(n_elements(namever) eq 0) then namever=2
if (namever ne 2) then begin
  print,'Invalid name version: ',namever, $
        ' thinlayerpeaks can only read NAME II files.'
  return
endif 
   
; Check that number of thick layers agrees with number of
; thin layer lists
if ( n_elements(thinlayerlist(0,*)) ne n_elements(thicklayerlist) ) then begin
  print, 'thinlayerlist and thicklayerlist have to have same length'
  return
endif

; Check that layer lists are not empty
if ( n_elements(thinlayerlist) eq 0) then begin
  print, 'layer lists have to be non-empty'
  return
endif

; Number of thick layers in output file
n_thicklayers = n_elements(thicklayerlist)
; Total number of thin layers (including layers that are to 
; be ignored and are denoted by '' in thinlayerlist)
n_totalthinlayers = n_elements(thinlayerlist)
; Number of thin layers
dummy = where(thinlayerlist ne '',n_thinlayers)
; (maximal) number of thin layers per thick layer
n_thinlayers_per_thicklayer = n_elements(thinlayerlist(*,0))

; Check that reductionmethod is either 'peak' or 'average'
if ( (reductionmethod ne 'average') and (reductionmethod ne 'peak') ) then begin
  print, 'reductionmethod has to be "peak" or "average"'
  return
endif

; Only use keyword 'factor' for calculation peak concentrations
if ( (reductionmethod ne 'peak') and (n_elements(factor) ne 0 ) ) then begin
  if (verbose eq 1) then begin
    print, 'keyword "factor" will be ignored.'
  endif
endif

; Construct filename for reading file title information
if(strpos(!version.os,'NT') gt 0)then begin
  pc=1
  delim='\'
endif else begin
  pc=0
  delim='/'
endelse
infilename=datadir+delim+'Fields_'+selgridin+'_'+seltime+'.txt'

; Read title information from input file
readtitle,infilename, $
          title,      $
          head1,      $
          head2,      $
          namever=namever

; Construct field selectors for thin layer columns. These selectors 
; are used by readfield to pick out the columns in the input file
; that will be reduced
;
; * Field (e.g. 'Air Concentration')
selfields_thinlayers=strarr(n_thinlayers)
selfields_thinlayers(*) = selfield
; * Vertical (thin) layers 
sellevels_thinlayers=thinlayerlist(where(thinlayerlist ne ''))
; * Species (e.g. 'VOLCANIC_ASH')
selspeciess_thinlayers=strarr(n_thinlayers)
selspeciess_thinlayers(*) = selspecies

; Read thin layers from input file
; Note that nolongitudeshift is set to 1 (= true) to avoid
; shifting the longitude array in the range (-180...180)
; when reading the field
readfield,datadir,seltime,selgridin,  $
          selfields_thinlayers,       $
          sellevels_thinlayers,       $
          selspeciess_thinlayers,     $
          lon2d,                      $
          lat2d,                      $
          fielddata_thinlayers,       $
          fieldhead_thinlayers,       $
          namever=namever,            $
          nolongitudeshift=1

nx = n_elements(fielddata_thinlayers(0,*,0))
ny = n_elements(fielddata_thinlayers(0,0,*))
          
; Construct new column headers for file with thick layers
fieldhead_thicklayers = strarr(n_thicklayers,6)
fieldhead_thicklayers(*,0) = fieldhead_thinlayers(0,0)
fieldhead_thicklayers(*,1) = selspecies
fieldhead_thicklayers(*,2) = fieldhead_thinlayers(0,2)
fieldhead_thicklayers(*,3) = selfield
fieldhead_thicklayers(*,4) = fieldhead_thinlayers(0,4)
fieldhead_thicklayers(*,5) = thicklayerlist(*)

; reduce data to thick layers

; Initialise thick layer data array
fielddata_thicklayers = fltarr(n_thicklayers,nx,ny)

; Peaks
if (reductionmethod eq 'peak') then begin
  for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
      thinlayerindex = 0
      for i=0,n_thicklayers-1 do begin
        maxa = 0.0
        dummy = where(thinlayerlist(*,i) ne '',nmax)
        for j=0,nmax-1 do begin
          if (fielddata_thinlayers(thinlayerindex,ix,iy) gt maxa) then begin
            maxa = fielddata_thinlayers(thinlayerindex,ix,iy)
          endif
          thinlayerindex = thinlayerindex + 1
        endfor
        fielddata_thicklayers(i,ix,iy) = factor*maxa
      endfor
    endfor
  endfor      
; Averages
endif else begin
  for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
      thinlayerindex = 0
      for i=0,n_thicklayers-1 do begin
        avga = 0.0
        dummy = where(thinlayerlist(*,i) ne '',nmax)
        for j=0,nmax-1 do begin
          avga = avga + fielddata_thinlayers(thinlayerindex,ix,iy)
          thinlayerindex = thinlayerindex + 1
        endfor
        fielddata_thicklayers(i,ix,iy) = avga/nmax
      endfor
    endfor
  endfor
endelse

; Write thick layers to output file
writefile,datadir,seltime,selgridout,  $
          title,                       $
          fieldhead_thicklayers,       $
          fielddata_thicklayers,       $
          lon2d,                       $
          lat2d,                       $
          boundarythreshold=boundarythreshold
           
end
