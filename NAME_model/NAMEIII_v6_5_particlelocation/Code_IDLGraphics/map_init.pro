PRO map_init, map_struct, _ref_extra=_ref_extra
  PLOT, map_struct.uv_box[[0,2]], $
    map_struct.uv_box[[1,3]],     $
    /NoData, XSTYLE=5, YSTYLE=5,  $
    _extra=_ref_extra
END
