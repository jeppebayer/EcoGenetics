#!/bin/env python
import re
import sys
import os

def search_and_replace_in_file(file, search_pattern, new_text):
    with open(file, 'r') as file_obj:
        old_text = file_obj.read()
    new_text = re.sub(search_pattern, new_text, old_text, count=1)
    with open (file, 'w') as file_obj:
        file_obj.write(new_text)

if __name__ == '__main__':
    if (len(sys.argv)-1) < 1:
        exit("Usage: {} <config.yaml>".format(sys.argv[0]))
    if not os.path.exists('workflow.py'):
        exit("Can't locate 'workflow.py' in current directory: {}".format(os.path.realpath(os.path.curdir)))
    if not os.path.exists(sys.argv[1]):
        exit("Can't locate indicated configuration file - {} - in current directory: {}".format(sys.argv[1], os.path.realpath(os.path.curdir)))
    if not os.path.splitext(sys.argv[1])[1] in ['.yml', '.yaml']:
        exit("Configuration file must be in YAML format ('.yml' or '.yaml')")
    search_and_replace_in_file('workflow.py', r'config\=yaml\.safe_load\(open\(\'.*\'\)\)', "config=yaml.safe_load(open('{}'))".format(sys.argv[1]))