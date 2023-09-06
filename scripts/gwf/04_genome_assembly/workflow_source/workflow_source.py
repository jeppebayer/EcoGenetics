from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def create_hic_workflow(config_file = glob.glob('*config.yaml')[0]):
    
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------

    config=yaml.safe_load(open(config_file))
    ACCOUNT=config['account']
    SPECIES_NAME=config['species_name']
    HIC_DIRECTORY=config['HiC_directory_path']
    REFERENCE_GENOME=config['reference_genome_path']
    WORK_DIR=config['working_directory']

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )

    hic_align = gwf.target_from_template(
        name = 'HiC_alignment_{}'.format(species_abbreviation(SPECIES_NAME)),
        template = hic_alignment(
            hic_directory = HIC_DIRECTORY,
            reference_genome = REFERENCE_GENOME,
            species_name = SPECIES_NAME,
            work_dir = WORK_DIR
        )
    )

    # test = gwf.target_from_template(
    #     name='test',
    #     template=hic_test()
    # )

    return gwf