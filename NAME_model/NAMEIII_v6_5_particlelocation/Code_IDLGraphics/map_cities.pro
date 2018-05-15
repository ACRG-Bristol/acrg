pro map_cities, shapefile, Root_dir=root_dir, Fill=fill, Color=color, Country=country, Label=label,$
                L_color=l_color, L_charsize=l_charsize, All=all, limit=limit
  COMPILE_OPT IDL2
  
  min_radius=0.05
  
  IF ~KEYWORD_SET(fill) THEN fill = 0
  
  if KEYWORD_SET(Country) and KEYWORD_SET(all) then message, 'Country and ALL are mutually exclusive. Please use only one'
    if KEYWORD_SET(All) then all_cities=1
  if (~KEYWORD_SET(Country)) then begin
    all_cities=1
    Country=''
  endif
  
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
  myshape=OBJ_NEW('IDLffShape',fpath)
  
  myshape->IDLffShape::GetProperty, N_ENTITIES=num_ent
  
  FOR x=0L, (num_ent-1) DO BEGIN
    entity = myshape->GetEntity(x)
    ; cities are shape_type =1
    ; point is just x, y
    attr= myshape->GetAttributes(x)
    
    if ((attr.attribute_4 eq Country) OR (KEYWORD_SET(All))) then begin
      ; bounds contains the x and y  min and max values
      shape_bounds=entity.bounds
      xmin=shape_bounds[0]
      ymin=SHAPE_BOUNDS[1]
      xmax=SHAPE_BOUNDS[4]
      ymax=SHAPE_BOUNDS[5]
      
      ; calculate the radius
      if (((xmax-xmin) eq 0) or ((ymax-ymin) eq 0)) then begin
        radius=min_radius
      endif else begin
        radius=SQRT( ((xmax-xmin)^2) + ((ymax-ymin)^2) )
        if (radius lt min_radius) then radius=min_radius
      endelse
      
      ; calculate centre point
      centre_x=(xmax - xmin)+xmin
      centre_y=(ymax - ymin)+ymin

      ; don't plot cities outside of plot limits
      if (KEYWORD_SET(limit)) then begin
        if (centre_x gt limit[1] or centre_x lt limit[0]) then continue
        if (centre_y gt limit[3] or centre_y lt limit[2]) then continue
      endif

      ; x and y for the label
      if (KEYWORD_SET(Label)) then begin
         ; take centre and add min_radius
         l_centre_x=centre_x + min_radius
         l_centre_y=centre_y
      endif

      
      if (~FILL) then begin
        PLOTS, CIRCLE(centre_x, centre_y, radius), color=COLOR,_extra=_extra
      endif else begin
        POLYFILL, CIRCLE(centre_x, centre_y, radius), color=COLOR,_extra=_extra
      endelse
      if (KEYWORD_SET(Label)) then begin
        xyouts, l_centre_x, l_centre_y, attr.attribute_0, color=L_color, CHARSIZE=L_charsize

      endif
    endif
    myshape->DestroyEntity, entity
  ENDFOR
  OBJ_DESTROY, myshape
end