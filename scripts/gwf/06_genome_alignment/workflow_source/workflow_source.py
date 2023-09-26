from gwf import Workflow
import glob, yaml, os
from workflow_templates import *

def multi_wga(config_file = glob.glob('*config.y*ml')[0]):
    """
    Template: Creates a multiple sequence alignment of whole genome data.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account'].strip()
    NAME: str = config['alignment_name'].strip()
    SPECIES_LIST: list = config['species_list']
    GENOME_LIST: list = config['genome_list']
    NEWICK: str = config['newick_tree']
    WORK_DIR: str = config['working_directory_path']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    top_dir = '{work_dir}/06_genome_alignment/{name}'.format(work_dir=WORK_DIR, name=NAME)
    os.makedirs(top_dir, exist_ok=True)

    cactus_align = gwf.target_from_template(
        name='{}_cactus_alignmnet'.format(NAME),
        template=run_cactus(
            seqfile=cactus_seqFile(SPECIES_LIST, GENOME_LIST, NEWICK, top_dir),
            file_prefix=NAME,
            output_path=top_dir
        )
    )

    return gwf