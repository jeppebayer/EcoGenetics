import sys
import os
sys.path.insert(0, os.path.realpath('../../workflow_source/'))
from workflow_source import *

gwf = create_vcf_workflow()