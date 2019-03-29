import itertools

import esmlab
import numpy as np
import pytest
import xarray as xr
from esmlab.datasets import open_dataset

from esmlab_regrid import regridder
from esmlab_regrid.core import _GridRef

esmlab.config.set({'regrid.gridfile_directory': './tests/esmlab-grid-files/'})


grid_src = ['T62']
grid_dst = ['CAM_f09']
methods = ['bilinear']
data_in = ['ncep_forecast_tseries']
xr_type = [xr.DataArray]
params = list(itertools.product(grid_src, grid_dst, methods, data_in, xr_type))


def test_grid_ref_constructor():
    grid_ref = _GridRef(name='T62', overwrite_existing=False)
    assert isinstance(grid_ref.ds, xr.Dataset)

    with pytest.raises(NotImplementedError):
        _GridRef(name='T78', overwrite_existing=False)


def test_regrid_init():

    R = regridder(
        name_grid_src='T62', name_grid_dst='CAM_f09', method='bilinear', overwrite_existing=True
    )

    assert isinstance(R, regridder)


@pytest.mark.parametrize('grid_src, grid_dst, method, data_in, xr_type', params)
def test_regrid_regrid(grid_src, grid_dst, method, data_in, xr_type):

    R = regridder(
        name_grid_src=grid_src, name_grid_dst=grid_dst, method=method, overwrite_existing=False
    )
    ds = open_dataset(name=data_in)
    dao = R(ds.t_10)
    assert isinstance(dao, xr_type)
    assert dao.data.shape == (1, 192, 288)

    with pytest.raises(NotImplementedError):
        dao = R(ds)

    with pytest.raises(ValueError):
        a = np.ones(10)
        dao = R(a)
