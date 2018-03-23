pro map_roads, shapefile, Root_dir=root_dir, Fill=fill, Color=color 
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
  if (res ne 1) then message, 'Shape file ' + shapefile + ' not readable'
  
  ;  myshape=OBJ_NEW('IDLffShape', FILEPATH(shapefile, SUBDIR=['resource', 'maps', 'shape']))
  myshape=OBJ_NEW('IDLffShape', fpath)
  
  myshape->IDLffShape::GetProperty, N_ENTITIES=num_ent
  
  FOR x=0L, (num_ent-1) DO BEGIN
    entity = myshape->GetEntity(x)
    attr= myshape->GetAttributes(x)
    ; attribute_4 is the name
    ; get the color for the level
    IF entity.N_PARTS GT 0 THEN BEGIN
      cuts = [*entity.parts, entity.n_vertices]
      
      FOR index = 0L, entity.N_PARTS-1 DO BEGIN
        IF ~fill then BEGIN
          PLOTS, (*entity.vertices)[0,cuts[index]:cuts[index+1]-1], (*entity.vertices)[1,cuts[index]:cuts[index+1]-1], color=color, _extra=_extra
        ENDIF ELSE BEGIN
        
          PolyFill, (*entity.vertices)[0,cuts[index]:cuts[index+1]-1], (*entity.vertices)[1,cuts[index]:cuts[index+1]-1], color=color, _extra=_extra
        ENDELSE
      ENDFOR
      
      myshape->DestroyEntity, entity
    ENDIF
  ENDFOR
  OBJ_DESTROY, myshape
end
