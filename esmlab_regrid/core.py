import os
from datetime import datetime
from warnings import warn

import ESMF
import esmlab
import numpy as np
import xarray as xr
import xesmf as xe


class _GridRef(object):
    """Get ESMF grid object from grid file."""

    def __init__(self, name, overwrite_existing=False):
        self.name = name
        grid_file_dir = esmlab.config.get('regrid.gridfile-directory')
        print(grid_file_dir)
        self.scrip_grid_file = f'{grid_file_dir}/{self.name}.nc'
        self._gen_grid_file(overwrite_existing=overwrite_existing)
        self._esmf_grid_from_scrip()

    def _gen_grid_file(self, overwrite_existing):
        """Generate a grid file. """
        if os.path.exists(self.scrip_grid_file) and not overwrite_existing:
            return
        raise NotImplementedError('Grid file generation is not implemented')

    def _esmf_grid_from_scrip(self):
        """Generate an ESMF grid object from a SCRIP grid file."""
        if not os.path.exists(self.scrip_grid_file):
            raise FileNotFoundError(f'file not found: {os.path.abspath(self.scrip_grid_file)}')

        self.ds = xr.open_dataset(self.scrip_grid_file)

        self.grid = ESMF.Grid(
            filename=self.scrip_grid_file,
            filetype=ESMF.api.constants.FileFormat.SCRIP,
            add_corner_stagger=True,
        )

        self.shape = self.ds.grid_dims.data
        mask = self.grid.add_item(ESMF.GridItem.MASK)
        mask[:] = self.ds.grid_imask.data.reshape(self.shape[::-1]).T


class Regridder(object):
    """Class to enable regridding between named grids."""

    def __init__(self, name_grid_src, name_grid_dst, method='bilinear', overwrite_existing=False):
        """

        Parameters
        ----------
        name_grid_src : string
                      Name of source grid.
        name_grid_dst : string
                      Name of destination grid.
        method : string, optional
                Regridding method. Options are:

                - 'bilinear'
                - 'conservative'
                - 'patch'
                - 'nearest_s2d'
                - 'nearest_d2s'

        overwrite_existing : bool, optional [Default=False]
                           Overwrite previously generated weight files.
        """

        self.name_grid_src = name_grid_src
        self.name_grid_dst = name_grid_dst
        self.grid_ref_src = _GridRef(name_grid_src)
        self.grid_ref_dst = _GridRef(name_grid_dst)
        self.method = method
        self.weight_file = None  # Will be set in _gen_weights() function
        self.N_src = self.grid_ref_src.shape[0] * self.grid_ref_src.shape[1]
        self.N_dst = self.grid_ref_dst.shape[0] * self.grid_ref_dst.shape[1]

        self._gen_weights(overwrite_existing=overwrite_existing)
        self.A = xe.smm.read_weights(self.weight_file, self.N_src, self.N_dst)

    def _gen_weights(self, overwrite_existing):
        """ Generate regridding weights """

        grid_file_dir = esmlab.config.get('regrid.gridfile-directory')
        weights_dir = f'{grid_file_dir}/weights'
        os.makedirs(weights_dir, exist_ok=True)
        self.weight_file = (
            f'{weights_dir}/{self.name_grid_src}_to_{self.name_grid_dst}_{self.method}'
        )

        if os.path.exists(self.weight_file):
            if overwrite_existing:
                print(f'Removing {os.path.abspath(self.weight_file)}')
                os.remove(os.path.abspath(self.weight_file))

            else:
                return

        print(f'Generating {os.path.abspath(self.weight_file)}')
        regrid = xe.backend.esmf_regrid_build(
            sourcegrid=self.grid_ref_src.grid,
            destgrid=self.grid_ref_dst.grid,
            method=self.method,
            filename=self.weight_file,
        )

        xe.backend.esmf_regrid_finalize(regrid)

    def __call__(
        self,
        data_in,
        renormalize=True,
        apply_mask=None,
        interp_coord=None,
        post_method=None,
        post_method_kwargs=None,
    ):
        """ Perform regridding on an `xarray.DataArray` or `xarray.Dataset`.

        Parameters
        ----------

        data_in : `xr.DataArray` or `xr.Dataset`
              The data to regrid
        renormalize : bool, optional [default=True]
                   Logical flag to trigger renormalization of the remapping
                   weights. This is useful if the remapping weight-file was
                   computed with a different missing value mask than `da_in`.
                   For instance, in regridding 3D ocean data, it is possible
                   to use a mapping file computed at the surface at each
                   successive depth level: setting `renormalize=True` will
                   ensure correct handling of missing values.
        interp_coord : dict, optional
                   Dictionary specifying dimension names an new coordinates;
                   passed to `xarray.DataArray.interp`.
                   New coordinate can be a scalar, array-like or DataArray.
                   If DataArrays are passed as new coordates, their dimensions
                   are used for the broadcasting.
        apply_mask : bool, optional [default=False]
                  Apply a missing-values mask after regridding operations.
        post_method : callable, optional
                 If provided, call this function on DataArray after regridding.
        post_method_kwargs : dict, optional
                  Keyword arguments to pass to `post_method`.

        Returns
        -------

        data_out : `xr.DataArray` or `xr.Dataset`
            The dataarray or dataset regridded to the destination grid.
        """

        if isinstance(data_in, xr.DataArray):
            return self.regrid_dataarray(
                data_in,
                renormalize=renormalize,
                apply_mask=apply_mask,
                interp_coord=interp_coord,
                post_method=post_method,
                post_method_kwargs=post_method_kwargs,
            )

        elif isinstance(data_in, xr.Dataset):
            raise NotImplementedError('Regridding method does support xarray Datasets yet!')

        else:
            raise ValueError(f'Data: {data_in} must be an xarray Dataset or DataArray')

    def regrid_dataarray(
        self,
        da_in,
        renormalize=True,
        interp_coord=None,
        apply_mask=None,
        post_method=None,
        post_method_kwargs=None,
    ):
        # Pull data, dims and coords from incoming DataArray
        data_src = da_in.data
        non_lateral_dims = da_in.dims[:-2]
        copy_coords = {d: da_in.coords[d] for d in non_lateral_dims if d in da_in.coords}

        # If renormalize == True, remap a field of ones
        if renormalize:
            ones_src = np.where(np.isnan(data_src), 0.0, 1.0)
            data_src = np.where(np.isnan(data_src), 0.0, data_src)

        # remap the field
        data_dst = xe.smm.apply_weights(
            self.A, data_src, self.grid_ref_dst.shape[1], self.grid_ref_dst.shape[0]
        )

        # Renormalize to include non-missing data_src
        if renormalize:
            ones_dst = xe.smm.apply_weights(
                self.A, ones_src, self.grid_ref_dst.shape[1], self.grid_ref_dst.shape[0]
            )
            ones_dst = np.where(ones_dst > 0.0, ones_dst, np.nan)
            data_dst = data_dst / ones_dst
            data_dst = np.where(ones_dst > 0.0, data_dst, np.nan)

        # reform into xarray.DataArray
        da_out = xr.DataArray(
            data_dst, name=da_in.name, dims=da_in.dims, attrs=da_in.attrs, coords=copy_coords
        )

        da_out.attrs['regrid_method'] = self.method
        new_history = f'\n{datetime.now()} esmlab.regrid <regrid> '
        if 'history' in da_out.attrs:
            da_out.attrs['history'] += new_history

        else:
            da_out.attrs['history'] = new_history

        # Interpolate coordinates (i.e., vertical)
        # Setup to copy lowest/hightest values  where extrapolation is needed
        if interp_coord:
            try:
                for dim, new_coord in interp_coord.items():
                    if dim in da_out.dims:
                        extrap_values = (da_out.isel(**{dim: 0}), da_out.isel(**{dim: -1}))

                        da_out = da_out.interp(
                            coords={dim: new_coord},
                            method='linear',
                            assume_sorted=True,
                            kwargs={'fill_value': extrap_values},
                        )

            except Exception as exc:
                raise exc

        # Appy a missing-values mask
        if apply_mask is not None:
            if apply_mask.dims != da_in.dims:
                warn(f'Masking {apply_mask.dims}; data have dims: {da_in.dims}')

            da_out = da_out.where(apply_mask)

        # Apply a post_method
        if post_method is not None:
            da_out = post_method(da_out, **post_method_kwargs)

        return da_out
