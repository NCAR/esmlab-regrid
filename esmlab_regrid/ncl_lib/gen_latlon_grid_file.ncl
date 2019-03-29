load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  dlon = stringtofloat(getenv("DLON"))
  dlat = stringtofloat(getenv("DLAT"))
  left_lon_corner = stringtofloat(getenv("LEFT_LON_CORNER"))
  grid_out_fname = getenv("GRID_OUT_FNAME")
  grid_type = getenv("GRID_TYPE")

  opt = True
  opt@ForceOverwrite = True
  opt@Debug = True
  opt@LLCorner = (/ -90.0+0.5*dlat,       left_lon_corner+0.5*dlon /)
  opt@URCorner = (/  90.0-0.5*dlat, 360.0+left_lon_corner-0.5*dlon /)
  latlon_to_SCRIP(grid_out_fname, grid_type, opt)

  f = addfile(grid_out_fname,"w")
  nx = 360/dlon
  ny = 180/dlat
  lat = fspan(opt@LLCorner(0),opt@URCorner(0),toint(ny))
  lon = fspan(opt@LLCorner(1),opt@URCorner(1),toint(nx))

  lat!0 = "lat"
  lon!0 = "lon"
  lat&lat = lat
  lon&lon = lon
  f->lat = lat
  f->lon = lon

  exit
end
status_exit(1)
