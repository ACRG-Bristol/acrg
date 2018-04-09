; altitude levels in World_Bartholomews.Altitude at the following levels:
; 0 - 200m - index 0
; 200-500m - index 1
; 500-1000m - index 2
; 1000 - 2000m - index 3
;
; fill=1 - polyfill
; fill=0 - plots
;
pro map_altitudes, shapefile, Root_dir=root_dir, Fill=fill, Alt_levels=alt_levels, Alt_colors=alt_colors, _extra=extra

  COMPILE_OPT IDL2
  
  IF ~KEYWORD_SET(fill) THEN fill = 0
  if (KEYWORD_SET(Root_dir)) then begin
    fpath = FILE_SEARCH(root_dir, shapefile )
  endif else begin
    ; use ITT resources
    fpath= FILEPATH(shapefile, SUBDIR=['resource', 'maps', 'shape'])
  endelse

  ; check if file exists
  res=FILE_TEST(fpath, /read)
  if (res ne 1) then message, 'Shape file ', shapefile, ' not readable'
  
  nlevels=N_ELEMENTS(Alt_levels)
  ncolors=N_ELEMENTS(Alt_colors)
  
  alevels=['land 0 - 200m', 'land 200 - 500m', 'land 500 - 1000m', 'land 1000 - 2000m']
  
  if (nlevels lt ncolors) then message, 'Specify the right number of colors to match the selected levels'
  
  myshape=OBJ_NEW('IDLffShape', fpath)
  
  myshape->IDLffShape::GetProperty, N_ENTITIES=num_ent

;
;first scan shapes + record levels so we can draw in bottom-up altitude order
; (also then ignores any levels that don't match)
;  
  levels = BYTARR(num_ent)
  FOR iEnt=0L, (num_ent-1) DO BEGIN
    attr= myshape->GetAttributes(iEnt)
    levels[iEnt] = WHERE(attr.attribute_2 EQ alevels[Alt_levels])
  ENDFOR

  FOR iLevel=0,nlevels-1 DO BEGIN
    p_color=ALT_COLORS[iLevel]
    aiEntsThisLevel = WHERE(levels Eq iLevel)
    nAtThisLevel = N_ELEMENTS(aiEntsThisLevel)
    
  FOR iInThisLevel=0L, (nAtThisLevel-1) DO BEGIN
    iEnt = aiEntsThisLevel[iInThisLevel]
    entity = myshape->GetEntity(iEnt)
    
;    res=WHERE(alevels[alt_levels] eq attr.attribute_2, count)
;    if (count gt 0) then begin
      ; get the color for the level
      IF entity.N_PARTS GT 0 THEN BEGIN
        cuts = [*entity.parts, entity.n_vertices]
        
        FOR index = 0L, entity.N_PARTS-1 DO BEGIN
          IF ~fill then BEGIN
            PLOTS, (*entity.vertices)[0,cuts[index]:cuts[index+1]-1], (*entity.vertices)[1,cuts[index]:cuts[index+1]-1], color=p_color, _extra=_extra
          ENDIF ELSE BEGIN
            x_vals=(*entity.vertices)[0,cuts[index]:cuts[index+1]-1]
            y_vals=(*entity.vertices)[1,cuts[index]:cuts[index+1]-1]

            ;note - must avoid values of exactly 180.0 - this confuses CONVERT_COORD !
            ;check for unexpected range
            IF MAX(x_vals LT -180.0001 OR x_vals GT 180.0001) THEN BEGIN
              PRINT,"x-180 : ",x_vals - 180.0
              PRINT,"x+180 : ",x_vals + 180.0
              MESSAGE,"unexpected longitudes for shape #"+STRTRIM(index,2)              
            ENDIF
            x_vals = x_vals < 179.9999

            res=CONVERT_COORD(x_vals, y_vals, DATA=1, TO_NORMAL=1)
            x_vals=res[0,*]
            y_vals=res[1,*]
            PolyFill, /NORMAL, x_vals, y_vals, color=p_color, _extra=_extra
          ENDELSE
        ENDFOR
        
        myshape->DestroyEntity, entity
      ENDIF
;    ENDIF  ;count
  ENDFOR  ;entity at required level
  ENDFOR  ;levels
  OBJ_DESTROY, myshape
  
end
