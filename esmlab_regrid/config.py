import os

import esmlab
import yaml

fn = os.path.join(os.path.dirname(__file__), 'regrid.yaml')
esmlab.config.ensure_file(source=fn, comment=False)

with open(fn) as f:
    defaults = yaml.safe_load(f)

esmlab.config.update(esmlab.config.config, defaults, priority='old')
esmlab.config.refresh()
