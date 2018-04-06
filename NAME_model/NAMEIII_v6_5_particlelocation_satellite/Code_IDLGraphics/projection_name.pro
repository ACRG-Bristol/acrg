pro projection_name, projection=projection, proj_name=proj_name
  ;
  ;Inputs - projection, a number or word indicating the projection
  ;Outputs - word describing projection
  ;-------------------------------------------------------------
  ;Procedure to convert between projection name and number
  ;SJL 22 October 2009
  ;
  ;-------------------------------------------------------------

  vartype=strarr(1)
  help, projection, output=vartype
  
  IF strmatch(vartype,'*STRING*') then begin
  
    proj_name=projection
    
  ENDIF ELSE begin
  
    CASE projection OF
    
      1: proj_name='stereographic'
      2: proj_name='orthographic'
      3: proj_name='lambertconic'
      4: proj_name='lambertazimuthal'
      5: proj_name='gnomic'
      6: proj_name='azimuthalequidistant'
      7: proj_name='satellite'
      8: proj_name='cylindrical'
      9: proj_name='mercator'
      10: proj_name='mollweide'
      11: proj_name='sinusoidal'
      12: proj_name='aitoff'
      13: proj_name='hammeraitoff'
      14: proj_name='albersequalareaconic'
      15: proj_name='transversemercator'
      16: proj_name='millercylindrical'
      17: proj_name='robinson'
      18: proj_name='lambertconicellipsoid'
      19: proj_name='goodeshomolosine'
      else: proj_name='cylindrical'
      
    ENDCASE
    
  ENDELSE
  
END
