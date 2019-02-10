import os
from subprocess import Popen, PIPE

import logging

import numpy as np
import xarray as xr

from netCDF4 import default_fillvals

from .config import known_grids, dir_root
from . import fill_POP_core

logger = logging.getLogger(__name__)


def fill_ocean_POP(da_in, mask, ltripole=False):
    '''Fill missing values on the POP grid by interative smoothing operation.'''

    non_lateral_dims = da_in.dims[:-2]

    if not non_lateral_dims:
        da_in = fill_ocean_POP_single_layer(da_in, mask, ltripole)

    elif non_lateral_dims == ('time',):
        if mask.dims != ('nlat', 'nlon'):
            raise ValueError('Mask dims do not match data')

        for l in range(da_in.shape[0]):
            da_in.values[l, :, :] = fill_ocean_POP_single_layer(da_in[l, :, :],
                                                                mask[:, :],
                                                                ltripole)

    elif non_lateral_dims == ('z_t',):
        if mask.dims != ('z_t', 'nlat', 'nlon'):
            raise ValueError('Mask dims do not match data')

        for k in range(da_in.shape[0]):
            da_in.values[k, :, :] = fill_ocean_POP_single_layer(da_in[k, :, :],
                                                                mask[k, :, :],
                                                                ltripole)
    elif non_lateral_dims == ('time', 'z_t',):
        if mask.dims != ('z_t', 'nlat', 'nlon'):
            raise ValueError('Mask dims do not match data')

        for l in range(da_in.shape[0]):
            for k in range(da_in.shape[1]):
                da_in.values[l, k, :, :] = fill_ocean_POP_single_layer(
                    da_in[l, :, :], mask[k, :, :],
                    ltripole)
    else:
        raise ValueError(f'Unknown dims: {non_lateral_dims}')

    return da_in


def fill_ocean_POP_single_layer(da_in, mask, ltripole=False):

    tol = 1.0e-4

    fillmask = (np.isnan(da_in) & mask).values
    if not fillmask.any():
        return da_in

    var_pass = da_in.values.astype(np.float32).T
    msv = default_fillvals['f4']
    var_pass[np.isnan(var_pass)] = msv

    add_attrs = {'note': 'fill_ocean_POP applied'}

    fill_POP_core.fill_pop_core(var=var_pass,
                                fillmask=fillmask.T,
                                msv=msv,
                                tol=tol,
                                ltripole=ltripole)

    var_pass[var_pass == msv] = np.nan

    return xr.DataArray(var_pass.T.astype(da_in.dtype),
                        dims=da_in.dims,
                        coords=da_in.coords,
                        attrs=da_in.attrs.update(add_attrs),
                        encoding=da_in.encoding)


def mask_3d_POP(grid):
    '''Construct a 3D mask for the POP grid.'''

    assert 'POP_' in grid, f'mask_3d_POP call on non-POP grid: {grid}'

    scrip_grid_file = gen_grid_file(grid)

    ds = xr.open_dataset(scrip_grid_file)
    ds = xr.merge((ds, open_vertical_grid(grid)))
    KMT = ds.KMT
    z_t = ds.z_t

    nk = len(z_t)
    nj = KMT.shape[0]
    ni = KMT.shape[1]

    MASK = (
        xr.DataArray(np.arange(0, len(z_t)), dims=('z_t')) *
        xr.DataArray(np.ones((nk, nj, ni)), dims=('z_t', 'nlat', 'nlon'),
                     coords={'z_t': z_t})
    )

    MASK = MASK.where(MASK <= KMT - 1)
    MASK.values = np.where(MASK.notnull(), True, False)

    return MASK


def _ncl(ncl_script):
    '''interface to NCL.'''

    cmd = ' && '.join(['module load intel/17.0.1',
                       'module load ncl/6.4.0',
                       'ncl ' + ncl_script])

    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = p.communicate()

    if p.returncode != 0:
        stdout = stdout.decode('UTF-8')
        stderr = stderr.decode('UTF-8')
        if stdout:
            logger.error('NCL error stdout: ' + stdout)
        if stderr:
            logger.error('NCL error stderr: ' + stderr)

        raise OSError('NCL failed.')


def gen_POP_grid_file(grid, grid_out_fname,
                      horiz_grid_fname, topography_fname, region_mask_fname,
                      type, lateral_dims, clobber=False):
    '''Generate POP SCRIP grid files.'''

    dir_grid_files = esmlab.get_options()["gridfile_directory"]

    os.environ['HORIZ_GRID_FNAME'] = horiz_grid_fname
    os.environ['TOPOGRAPHY_FNAME'] = topography_fname
    os.environ['REGION_MASK_FNAME'] = region_mask_fname
    os.environ['GRID_TYPE'] = type

    os.environ['NY'] = '{0:d}'.format(lateral_dims[0])
    os.environ['NX'] = '{0:d}'.format(lateral_dims[1])

    os.environ['GRID_OUT_FNAME'] = grid_out_fname
    os.environ['VERT_GRID_FILE_OUT'] = f'{dir_grid_files}/{grid}_vert.nc'

    ncl_script = os.path.join(os.path.dirname(__file__),
                              'ncl_lib/gen_POP_grid_file.ncl')

    _ncl(ncl_script)


def open_vertical_grid(grid):
    '''Return an `xarray.Dataset` with the vertical grid of "grid".'''

    assert grid in known_grids, f'Unknown grid: {grid}'
    info = known_grids[grid]

    assert 'open_vertical_grid' in info, f'No vertical grid defined for {grid}'
    info = info['open_vertical_grid']

    assert 'vert_grid_file' in info, f'No vertical grid file for {grid}'
    vert_grid_file = info['vert_grid_file']

    if not os.path.isfile(vert_grid_file):
        vert_grid_file = os.path.join(dir_root, vert_grid_file)
        if not os.path.isfile(vert_grid_file):
            raise OSError(f'Missing {vert_grid_file}')

    depth_coord_name = info['depth_coord_name']
    depth_units = info['depth_units']

    tmp = np.loadtxt(vert_grid_file)
    dz = xr.DataArray(tmp[:, 0], dims=(depth_coord_name),
                      attrs={'long_name': 'layer thickness',
                             'units': depth_units})

    depth_edges = np.concatenate(([0.], np.cumsum(dz)))
    depth = xr.DataArray(depth_edges[0:-1] + 0.5 * dz,
                         dims=(depth_coord_name),
                         attrs={'long_name': 'depth',
                                'units': depth_units})
    return xr.Dataset({depth_coord_name: depth, 'dz': dz})


def gen_latlon_grid_file(grid, grid_out_fname,
                         grid_type, dlon, dlat, left_lon_corner,
                         clobber=False):
    '''Generate latlon grid file.'''

    os.environ['DLON'] = '{0:f}'.format(dlon)
    os.environ['DLAT'] = '{0:f}'.format(dlat)
    os.environ['LEFT_LON_CORNER'] = '{0:f}'.format(left_lon_corner)
    os.environ['GRID_TYPE'] = grid_type
    os.environ['GRID_OUT_FNAME'] = grid_out_fname

    ncl_script = os.path.join(os.path.dirname(__file__),
                              'ncl_lib/gen_latlon_grid_file.ncl')

    _ncl(ncl_script)


def gen_rectilinear_grid_file(grid, grid_out_fname,
                              latlon_file,
                              clobber=False):
    '''Generate rectilinear grid file.'''

    os.environ['LATLON_FILE'] = latlon_file
    os.environ['GRID_OUT_FNAME'] = grid_out_fname

    ncl_script = os.path.join(os.path.dirname(__file__),
                              'ncl_lib/gen_rectilinear_grid_file.ncl')

    _ncl(ncl_script)


def gen_grid_file(grid, clobber=False):
    '''Generate a SCRIP grid file for "grid"'''

    dir_grid_files = esmlab.get_options()["gridfile_directory"]

    grid_out_fname = f'{dir_grid_files}/{grid}.nc'
    if os.path.exists(grid_out_fname) and not clobber:
        return grid_out_fname

    assert grid in known_grids, f'Unknown grid: {grid}'

    info = known_grids[grid]
    regrid_method = info['gen_grid_file']['method']
    kwargs = info['gen_grid_file']['kwargs']

    logger.info(f'generating grid file {grid_out_fname}')
    gen_grid_method = globals()[regrid_method]
    gen_grid_method(grid, grid_out_fname, **kwargs, clobber=clobber)

    return grid_out_fname


def gen_all_grids(clobber=False):
    '''Generate all known grids.'''
    for grid in known_grids.keys():
        gen_grid_file(grid, clobber=clobber)
