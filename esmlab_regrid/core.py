import os

import ESMF
import esmlab
import numpy as np
import xarray as xr
import xesmf as xe


class _GridRef(object):
    """Get ESMF grid object from grid file."""

    def __init__(self, name, overwrite_existing):
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
