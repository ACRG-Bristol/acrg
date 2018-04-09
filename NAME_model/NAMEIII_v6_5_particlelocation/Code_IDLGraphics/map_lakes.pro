;+
; NAME:
;  map_lakes
; PURPOSE:
;  Plot many of the larger lakes on a map
; CALLING SEQUENCE:
;  map_lakes
; KEYWORD PARAMETERS:
;      Fill - option to fill the lakes using polyfill
;      All parameters available to either POLYFILL or PLOTS depending on the fill option above
; INPUT PARAMETERS:
;      None.
; EXAMPLE:
;      ; Plot the lakes around the Caspian Sea
;      device, decomposed=1
;      blue = 'ff0000'x & green = '00ff00'x & black = '000000'x & white = 'ffffff'x
;      !P.background = white
;      map_set, limit=[10,20,50,70], /NOERASE
;      map_continents, COLOR=green, /FILL_CONTINENTS, /HIRES
;      map_continents, color=black, thick=2, /hires
;      map_lakes, color=white, /fill
;      map_lakes, color=black, thick=2
; LAKES PLOTTED:
;        Lake Athabasca, Nettilling Lake, Lake Onega, Lake Ladoga, Lake Vanern, Reindeer Lake, Lake Baikal, Cedar Lake,
;        Lake Nipigon, Lake Superior, Caspian Sea, Aral Sea, Lake Balkhash, Lake Michigan, Lake Ontario, Issyk Kul,
;        Great Salt Lake, Lake Urmia, Koko Nor, Lake Chad, Lake Tana, Lake Nicaragua, Lake Turkana, Lake Albert, Lake Victoria,
;        Lake Tanganyika, Lake Nyasa, Lake Titicaca, Lake Winnipeg, Lake Winnipegosis, Lake St. Martin, Lake Manitoba,
;        Great Bear Lake, Great Slave Lake, Lake Huron, Lake Erie, Lake St. Claire
; MODIFICATION HISTORY:
;      See FCM
;-
pro map_lakes, _extra=_extra, fill=fill
  ; Procedure to plot the data of the lakes shapefile
  IF ~KEYWORD_SET(fill) THEN fill = 0
    
  myshape=OBJ_NEW('IDLffShape', FILEPATH('lakes.shp', SUBDIR=['resource', 'maps', 'shape']))
  
  myshape->IDLffShape::GetProperty, N_ENTITIES=num_ent
  
  FOR x=0, (num_ent-1) DO BEGIN
    entity = myshape->GetEntity(x)
    
    IF entity.N_PARTS GT 0 THEN BEGIN
      cuts = [*entity.parts, entity.n_vertices]
      
      FOR index = 0L, entity.N_PARTS-1 DO BEGIN
        IF ~fill then BEGIN
          PLOTS, (*entity.vertices)[0,cuts[index]:cuts[index+1]-1], (*entity.vertices)[1,cuts[index]:cuts[index+1]-1], _extra=_extra
        ENDIF ELSE BEGIN
          PolyFill, (*entity.vertices)[0,cuts[index]:cuts[index+1]-1], (*entity.vertices)[1,cuts[index]:cuts[index+1]-1], _extra=_extra
        ENDELSE
      ENDFOR
      
      myshape->DestroyEntity, entity
    ENDIF
    
  END
  OBJ_DESTROY, myshape
END
