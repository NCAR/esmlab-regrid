import os
import logging
import yaml

logging.basicConfig(level=logging.INFO)

dir_root = os.path.dirname(os.path.abspath(__file__))
config_file = os.path.join(dir_root, 'grids.yml')

# TODO: default to current dir and no defined grids
if not os.path.isfile(config_file):
    raise OSError('config file is missing')

with open(config_file, 'r') as f:
    regrid_database = yaml.load(f)

dir_grid_files = regrid_database['config']['dir_grid_files']
if not os.path.exists(dir_grid_files):
    os.makedirs(dir_grid_files)

dir_weight_files = regrid_database['config']['dir_weight_files']
if not os.path.exists(dir_weight_files):
    os.makedirs(dir_weight_files)

known_grids = regrid_database['grids']
