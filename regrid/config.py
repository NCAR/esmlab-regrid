import os
import logging
import yaml

import esmlab

logger = logging.getLogger(__name__)
logger.setLevel(level=logging.INFO)

dir_root = os.path.dirname(os.path.abspath(__file__))
grid_defitions_file = esmlab.get_options()["grid_defitions_file"]

if not os.path.isfile(grid_defitions_file):
    raise OSError('config file is missing')

with open(grid_defitions_file, 'r') as f:
    regrid_database = yaml.load(f)

known_grids = regrid_database['grids']
