
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  latlon_file = getenv("LATLON_FILE")
  grid_out_fname = getenv("GRID_OUT_FNAME")

  opt = True
  opt@ForceOverwrite = True
  opt@Debug = True

  f = addfile(latlon_file, "r")
  lat = f->lat
  lon = f->lon

  lat!0 = "lat"
  lon!0 = "lon"
  lat&lat = lat
  lon&lon = lon

  rectilinear_to_SCRIP(grid_out_fname, lat, lon, opt)
  f = addfile(grid_out_fname, "w")
  f->lat = lat
  f->lon = lon

  exit
end
status_exit(1)
