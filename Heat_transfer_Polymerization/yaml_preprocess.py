#!/usr/bin/python3

import yaml
import sys
file = sys.argv[1]
with open(file) as f:
    y = yaml.safe_load(f)
    y['num_params']['minlevel'] = 6
    print(yaml.dump(y, default_flow_style=False, sort_keys=False))
