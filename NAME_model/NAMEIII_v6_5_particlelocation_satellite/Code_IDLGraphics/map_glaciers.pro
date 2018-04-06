pro map_glaciers, shapefile, Root_dir=root_dir, G_name=g_name, Fill=fill, $
   G_Color=g_color, Add_Label=Add_label,L_color=l_color, L_orientation=l_orientation, Label_only=label_only, Names=names
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
  myshape=OBJ_NEW('IDLffShape', fpath )
  
  myshape->IDLffShape::GetProperty, N_ENTITIES=num_ent
  
  ; we are just going to output the names
  if ARG_PRESENT(Names) then begin
    ; create max size array
    temp_names=strarr(num_ent)
  endif
  
  FOR x=0L, (num_ent-1) DO BEGIN
    entity = myshape->GetEntity(x)
    attr= myshape->GetAttributes(x)
    
    ; check for names
    if ARG_PRESENT(Names) then begin
      temp_names[x]=strtrim(attr.attribute_4,2)
      CONTINUE
    endif
    ; attribute_4 is the name
    if (KEYWORD_SET(G_NAME)) then begin
      ; if we have been given an array of
      ; names then check if we have a match
      ; skip if we don't otherwise continue
      name=strtrim(attr.attribute_4,2)
      if strlen(name) eq 0 then continue
      ; we have a name
      bfound=where(g_name eq name, name_count)
      if (name_count eq 0) then CONTINUE
    endif
    ; get the color for the level
    IF entity.N_PARTS GT 0 THEN BEGIN
      cuts = [*entity.parts, entity.n_vertices]
      
      FOR index = 0L, entity.N_PARTS-1 DO BEGIN
      x_vals=(*entity.vertices)[0,cuts[index]:cuts[index+1]-1]
      y_vals=(*entity.vertices)[1,cuts[index]:cuts[index+1]-1]
        if (~keyword_set(label_only)) then begin
          IF ~fill then BEGIN
            PLOTS, x_vals, y_vals, color=g_color, _extra=_extra
          ENDIF ELSE BEGIN
          
            PolyFill, x_vals, y_vals, color=g_color, _extra=_extra
          ENDELSE
        endif
        if(KEYWORD_SET(Add_Label)) then begin
              ; find mid point
              centre_x = total(x_vals)/n_elements(x_vals)
              centre_y= total(y_vals)/n_elements(y_vals)
              XYOUTS, centre_x, centre_y, attr.attribute_4, color=l_color
        endif
      ENDFOR
      
      myshape->DestroyEntity, entity
    ENDIF
  ENDFOR
  OBJ_DESTROY, myshape
  ; check for names only
  if ARG_PRESENT(Names) then begin
    ; there will be some elements that don't have a name.
    ; process the temp_names array and then assign to names
    res=where(strlen(temp_names) gt 0, count)
    if (count gt 0) then  begin
      Names=temp_names[res]
    endif else begin
      Names=''
    endelse
  endif
end
