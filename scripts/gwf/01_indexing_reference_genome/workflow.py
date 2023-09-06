from gwf import Workflow
from workflow_templates import *
import yaml

# --------------------------------------------------
#                  Configuration
# --------------------------------------------------

config = yaml.safe_load(open('config.yaml'))
REFERENCE_GENOME_PATH=config['reference_genome_path']
ACCOUNT=config['account']

# --------------------------------------------------
#                  Workflow
# --------------------------------------------------
gwf = Workflow(
    defaults={'account': ACCOUNT}
)

# Loops through the reference genome path to find all relevant files
reference_genomes = gwf.glob(REFERENCE_GENOME_PATH)

# Indexes each located reference genome and gives the jobs proper names
gwf.map(index_reference_genome_loop, reference_genomes, name=get_reference_name)
