import esmlab
import pytest
import xarray as xr

from esmlab_regrid.core import _GridRef

esmlab.config.set({'regrid.gridfile_directory': './tests/esmlab-grid-files/'})


def test_grid_ref_constructor():
    grid_ref = _GridRef(name='T62', overwrite_existing=False)
    assert isinstance(grid_ref.ds, xr.Dataset)

    with pytest.raises(NotImplementedError):
        _GridRef(name='T78', overwrite_existing=False)
