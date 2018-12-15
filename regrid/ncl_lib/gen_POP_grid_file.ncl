;
; gen_POP_grid_files.ncl
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

pi = 3.14159265358979323846d0

begin
;----------------------------------------------------------------------
;-- setup
;----------------------------------------------------------------------

  setfileoption("bin", "ReadByteOrder", "BigEndian")

  horiz_grid_fname = getenv("HORIZ_GRID_FNAME")
  topography_fname = getenv("TOPOGRAPHY_FNAME")
  region_mask_fname = getenv("REGION_MASK_FNAME")

  vert_grid_file = getenv("VERT_GRID_FILE")

  ltripole = getenv("GRID_TYPE") .eq. "tripole"
  km = stringtointeger(getenv("NK"))
  nlat = stringtointeger(getenv("NY"))
  nlon = stringtointeger(getenv("NX"))

  grid_out_fname = getenv("GRID_OUT_FNAME")
  vert_grid_file_out = getenv("VERT_GRID_FILE_OUT")

  print("ltripole = "+ltripole)

  ;
  ; read fields from horiz_grid_fname
  ;
  print("reading horiz_grid from "+horiz_grid_fname)
  if isfilepresent(horiz_grid_fname) then
      grid_file = addfile(horiz_grid_fname,"r")
      ULAT := grid_file->ULAT * (pi/180.)
      ULON := grid_file->ULON * (pi/180.)
      HTN  := grid_file->HTN
      HTE  := grid_file->HTE
      delete(grid_file)
  else
      grid_file_contents := cbinread(horiz_grid_fname, (/ 7, nlat, nlon /), "double")
      ULAT := grid_file_contents(0,:,:)
      ULON := grid_file_contents(1,:,:)
      HTN  := grid_file_contents(2,:,:)
      HTE  := grid_file_contents(3,:,:)
  end if

  ;
  ; construct coordinates of grid centers
  ;
  print("constructing coordinates of grid centers")

  TLAT := new((/ nlat, nlon /), double)
  TLON := new((/ nlat, nlon /), double)

  do j = 1, nlat-1
     jm1 = j-1
     do i = 0, nlon-1
        im1 = mod(i-1+nlon, nlon)

        tmp = cos(ULAT(jm1,im1))
        xsw = cos(ULON(jm1,im1))*tmp
        ysw = sin(ULON(jm1,im1))*tmp
        zsw = sin(ULAT(jm1,im1))

        tmp = cos(ULAT(jm1,i))
        xse = cos(ULON(jm1,i))*tmp
        yse = sin(ULON(jm1,i))*tmp
        zse = sin(ULAT(jm1,i))

        tmp = cos(ULAT(j,im1))
        xnw = cos(ULON(j,im1))*tmp
        ynw = sin(ULON(j,im1))*tmp
        znw = sin(ULAT(j,im1))

        tmp = cos(ULAT(j,i))
        xne = cos(ULON(j,i))*tmp
        yne = sin(ULON(j,i))*tmp
        zne = sin(ULAT(j,i))

        xc = 0.25 * (xsw + xse + xnw + xne)
        yc = 0.25 * (ysw + yse + ynw + yne)
        zc = 0.25 * (zsw + zse + znw + zne)

        r = sqrt(xc*xc + yc*yc + zc*zc)

        TLAT(j,i) = asin(zc/r)
        TLON(j,i) = atan2(yc, xc)
     end do
  end do

  ; generate bottom row vals by extrapolating from 2 rows above
  TLAT(0,:) = TLAT(1,:) - (TLAT(2,:) - TLAT(1,:))
  TLON(0,:) = TLON(1,:) - (TLON(2,:) - TLON(1,:))

  ;
  ; construct coordinates of grid centers
  ;
  print("constructing DXT, DYT, TAREA")

  DXT := new((/ nlat, nlon /), double)
  DXT(1:nlat-1,:) = 0.5 * (HTN(0:nlat-2,:) + HTN(1:nlat-1,:))
  ; generate bottom row vals by extrapolating HTN from 2 rows above
  DXT(0,:) = 0.5 * (2*HTN(0,:)-HTN(1,:) + HTN(0,:))

  DYT := new((/ nlat, nlon /), double)
  DYT(:,1:nlon-1) = 0.5 * (HTE(:,0:nlon-2) + HTE(:,1:nlon-1))
  ; assume zonal periodic bc
  DYT(:,0) = 0.5 * (HTE(:,nlon-1) + HTE(:,0))

  TAREA := DXT * DYT

  ;
  ; construct degree versions of coordinates
  ;

  ULATd := (180.0 / pi) * ULAT
  ULONd := (180.0 / pi) * ULON
  TLATd := (180.0 / pi) * TLAT
  TLONd := (180.0 / pi) * TLON

  ;
  ; construct coordinates of grid corners
  ; storage order of corners is counter-clockwise from NE: NE, NW, SW, SE
  ;
  print("constructing coordinates of grid corners")

  GridCornerLat2d := new((/ nlat, nlon, 4 /), double)
  ; NE corner
  GridCornerLat2d(:,:,0) = ULATd
  ; NW corner (copy from NE corner of column to the left, assume zonal periodic bc)
  GridCornerLat2d(:,1:nlon-1,1) = GridCornerLat2d(:,0:nlon-2,0)
  GridCornerLat2d(:,0,1) = GridCornerLat2d(:,nlon-1,0)
  ; SW corner (copy from NW corner of row below, bottom row is extrapolated from 2 rows above)
  GridCornerLat2d(1:nlat-1,:,2) = GridCornerLat2d(0:nlat-2,:,1)
  GridCornerLat2d(0,:,2) = GridCornerLat2d(1,:,2) - (GridCornerLat2d(2,:,2) - GridCornerLat2d(1,:,2))
  ; SE corner (copy from NE corner of row below, bottom row is extrapolated from 2 rows above)
  GridCornerLat2d(1:nlat-1,:,3) = GridCornerLat2d(0:nlat-2,:,0)
  GridCornerLat2d(0,:,3) = GridCornerLat2d(1,:,3) - (GridCornerLat2d(2,:,3) - GridCornerLat2d(1,:,3))

  GridCornerLon2d := new((/ nlat, nlon, 4 /), double)
  ; NE corner
  GridCornerLon2d(:,:,0) = ULONd
  ; NW corner (copy from NE corner of column to the left, assume zonal periodic bc)
  GridCornerLon2d(:,1:nlon-1,1) = GridCornerLon2d(:,0:nlon-2,0)
  GridCornerLon2d(:,0,1) = GridCornerLon2d(:,nlon-1,0)
  ; SW corner (copy from NW corner of row below, bottom row is extrapolated from 2 rows above
  GridCornerLon2d(1:nlat-1,:,2) = GridCornerLon2d(0:nlat-2,:,1)
  GridCornerLon2d(0,:,2) = GridCornerLon2d(1,:,2) - (GridCornerLon2d(2,:,2) - GridCornerLon2d(1,:,2))
  ; SE corner (copy from NE corner of row below, bottom row is extrapolated from 2 rows above
  GridCornerLon2d(1:nlat-1,:,3) = GridCornerLon2d(0:nlat-2,:,0)
  GridCornerLon2d(0,:,3) = GridCornerLon2d(1,:,3) - (GridCornerLon2d(2,:,3) - GridCornerLon2d(1,:,3))

  ;
  ; flatten corner coordinate arrays
  ;
  GridCornerLat := reshape(GridCornerLat2d, (/ nlat*nlon, 4 /))
  GridCornerLon := reshape(GridCornerLon2d, (/ nlat*nlon, 4 /))

  ;
  ; generate mask from topography file
  ;
  print("reading KMT from "+topography_fname)
  if isfilepresent(topography_fname) then
      file_topo = addfile(topography_fname,"r")
      KMT = file_topo->KMT
      delete(file_topo)
  else
      KMT := cbinread(topography_fname, (/ nlat, nlon /), "integer")
  end if
  GridMask := where(KMT .gt. 0, 1, 0)

  ;
  ; read REGION_MASK, for inclusion in generated grid file
  ;
  print("reading REGION_MASK from "+region_mask_fname)
  if isfilepresent(region_mask_fname) then
      file_rmask = addfile(region_mask_fname,"r")
      REGION_MASK = file_rmask->REGION_MASK
      delete(file_rmask)
  else
      REGION_MASK := cbinread(region_mask_fname, (/ nlat, nlon /), "integer")
  end if

  ;
  ; call curvilinear_to_SCRIP function
  ;

  opt = True
  opt@ForceOverwrite = True
  opt@Debug = True
  opt@GridCornerLat = GridCornerLat
  opt@GridCornerLon = GridCornerLon
  opt@GridMask = GridMask
  print("creating grid file "+grid_out_fname)
  curvilinear_to_SCRIP(grid_out_fname, TLATd, TLONd, opt)
  delete(opt)

  ;
  ; add two-dimensional grid variables to generated file
  ;
  ULATd!0 = "nlat"
  ULATd!1 = "nlon"
  ULATd@units = "degrees_north"

  ULONd!0 = "nlat"
  ULONd!1 = "nlon"
  ULONd@units = "degrees_east"

  TLATd!0 = "nlat"
  TLATd!1 = "nlon"
  TLATd@units = "degrees_north"

  TLONd!0 = "nlat"
  TLONd!1 = "nlon"
  TLONd@units = "degrees_east"

  DXT!0 = "nlat"
  DXT!1 = "nlon"
  DXT@units = "cm"

  DYT!0 = "nlat"
  DYT!1 = "nlon"
  DYT@units = "cm"

  TAREA!0 = "nlat"
  TAREA!1 = "nlon"
  TAREA@units = "cm^2"

  KMT!0 = "nlat"
  KMT!1 = "nlon"

  REGION_MASK!0 = "nlat"
  REGION_MASK!1 = "nlon"

  print("appending 2D grid variables to grid file")
  f = addfile(grid_out_fname, "w")
  f@horiz_grid_fname  = horiz_grid_fname
  f@topography_fname  = topography_fname
  f@region_mask_fname = region_mask_fname
  f->ULATd = ULATd
  f->ULONd = ULONd
  f->TLATd = TLATd
  f->TLONd = TLONd
  f->DXT = DXT
  f->DYT = DYT
  f->TAREA = TAREA
  f->KMT = KMT
  f->REGION_MASK = REGION_MASK
  if (ltripole) then
     f@ltripole = 1
  else
     f@ltripole = 0
  end if
  delete(f)

  exit
end
status_exit(1)
