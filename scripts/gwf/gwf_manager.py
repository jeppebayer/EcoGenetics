#!/bin/env python

import os, sys, yaml

# Configuration file location
config_file = os.path.dirname(os.path.realpath(sys.argv[0])) + '/gwf_manager_list.yaml'

# Test if configuration file exits and if not creates one
if not os.path.exists(config_file):
    default_config = dict(projects = None)
    with open(os.path.dirname(os.path.realpath(sys.argv[0])) + '/gwf_manager_list.yaml', 'w') as default_file:
        yaml.dump(default_config, default_file)

# Loads configuration file
config=yaml.safe_load(open(config_file))

print(config['projects'])