#! /usr/bin/env python
#BSUB -P NCGD0011
#BSUB -W 12:00
#BSUB -n 1
#BSUB -J gen_POP_grid_files
#BSUB -o logs/regrid.%J
#BSUB -e logs/regrid.%J
#BSUB -q geyser
#BSUB -N
import os
from subprocess import call
import regrid

clobber = False

datestr = regrid.nowstr

#------------------------------------------------
#--- generate POP grids
#------------------------------------------------
if clobber or any([not os.path.exists(regrid.grid_file(g))
                   for g in
                   ['POP_gx3v7','POP_gx1v6','POP_gx1v7','POP_tx1v1','POP_tx0.1v2',
                    'POP_tx0.1v3_62lev']]):

    env = os.environ.copy()
    env['GRID_FILES_OUTDIR'] = regrid.diro['grids']

    call(['ncl','-n','gen_POP_grid_files.ncl'],env=env)

dst_grids = ['POP_gx3v7','POP_gx1v6','POP_gx1v7','POP_tx1v1','POP_tx0.1v2']
src_grids = ['POP_gx1v7','POP_gx1v6']

#------------------------------------------------
#--- generate rectilinear grids
#------------------------------------------------

rectilinear_grids = [{'grid_name': 'T62',
                      'latlon_file' : '/glade/p/cesm/cseg/inputdata/atm/datm7/NYF/nyf.ncep.T62.050923.nc'},
                     {'grid_name': 'f09',
                      'latlon_file' : '/glade/p/work/mclong/grids/f09_f09.nc'}]

print('-'*40)
print('generating rectilinear grids')
for d in rectilinear_grids:
    grid = '_'.join(['rectilinear',d['grid_name']])
    src_grids.append(grid)
    fname = regrid.grid_file(grid)
    if not os.path.exists(fname) or clobber:
        print('grid_name = %s'%grid)
        d.update({'grid_out_fname' : fname})
        ok = regrid.gen_rectilinear_grid_file(**d)

#------------------------------------------------
#--- generate latlon grids
#------------------------------------------------
EorW = lambda x: 'W' if x < 0 else 'E'
lon_str = lambda lon: '%d%s'%(abs(lon),EorW(lon))

latlon_grids = [{'grid_type': '1x1',
                 'dlat': 1.,
                 'dlon': 1.,
                 'left_lon_corner' : -180.},
                {'grid_type': '1x1',
                 'dlat': 1.,
                 'dlon': 1.,
                 'left_lon_corner' : 0.},
                {'grid_type': '1x1',
                 'dlat': 1.,
                 'dlon': 1.,
                 'left_lon_corner' : 20.},
                {'grid_type': '0.25x0.25',
                 'dlat': 0.25,
                 'dlon': 0.25,
                 'left_lon_corner' : -180.}]
print('-'*40)
print('generating latlat grids')
for d in latlon_grids:
    left_lon_corner_str = lon_str(d['left_lon_corner'])
    grid = '_'.join(['latlon',d['grid_type'],left_lon_corner_str])
    src_grids.append(grid)
    fname = regrid.grid_file(grid)
    if not os.path.exists(fname) or clobber:
        print('grid_name = %s'%grid)
        d.update({'grid_out_fname' : fname})
        ok = regrid.gen_latlon_grid_file(**d)

#------------------------------------------------
#--- make weights files
#------------------------------------------------
print('-'*40)
print('make weights files')

#-- list src/dst grids
dst_grids = ['POP_tx0.1v3_62lev','POP_gx1v7']
src_grids = ['latlon_0.25x0.25_180W','POP_gx1v6','POP_gx1v7']

for src_grid in src_grids:
    for dst_grid in dst_grids:
        for interp_method in ['bilinear','conserve']:
            wgtFile = regrid.wgt_file(src_grid,dst_grid,interp_method)
            srcGridFile = regrid.grid_file(src_grid)
            dstGridFile = regrid.grid_file(dst_grid)

            if not os.path.exists(wgtFile) or clobber:
                ok = regrid.gen_weight_file(wgtFile = wgtFile,
                                            srcGridFile = srcGridFile,
                                            dstGridFile = dstGridFile,
                                            InterpMethod = interp_method)

print
