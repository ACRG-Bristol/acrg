pro map_rivers, shapefile, Root_dir=root_dir, G_name=g_name, Fill=fill, R_Color=r_color,$
 Label=label,L_color=l_color, River_name=river_name, RiverSize=RiverSize
  COMPILE_OPT IDL2
  
  IF ~KEYWORD_SET(fill) THEN fill = 0
  IF ~KEYWORD_SET(RiverSize) then RiverSize='A'
  
  if (KEYWORD_SET(Root_dir)) then begin
    fpath = FILE_SEARCH(root_dir, shapefile )
    ;fpath = FILEPATH(shapefile,root_dir=Root_dir, SUBDIR=['rivers', 'other_rivers', 'secondary_rivers', 'tertiary_rivers'])
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
    ; attribute_1 is the river size
    ; get the color for the level
    if strmatch(attr.attribute_1, 'river ['+Riversize+']') then begin
      IF entity.N_PARTS GT 0 THEN BEGIN
        cuts = [*entity.parts, entity.n_vertices]
        
        FOR index = 0L, entity.N_PARTS-1 DO BEGIN
          IF ~fill then BEGIN
            PLOTS, (*entity.vertices)[0,cuts[index]:cuts[index+1]-1], (*entity.vertices)[1,cuts[index]:cuts[index+1]-1], color=gcolor, _extra=_extra
          ENDIF ELSE BEGIN
            PolyFill, (*entity.vertices)[0,cuts[index]:cuts[index+1]-1], (*entity.vertices)[1,cuts[index]:cuts[index+1]-1], color=gcolor, _extra=_extra
            if(KEYWORD_SET(Label)) then begin
              ; find mid point
              centre_x = total((*entity.vertices)[0,cuts[index]:cuts[index+1]-1])/n_elements((*entity.vertices)[0,cuts[index]:cuts[index+1]-1])
              centre_y= total((*entity.vertices)[1,cuts[index]:cuts[index+1]-1])/n_elements((*entity.vertices)[1,cuts[index]:cuts[index+1]-1])
              XYOUTS, centre_x, centre_y, attr.attribute_4, color=l_color
            endif
          ENDELSE
        ENDFOR ; parts
      endif ; if parts
    endif
    myshape->DestroyEntity, entity
  ENDFOR ; entities
  OBJ_DESTROY, myshape
  
end
