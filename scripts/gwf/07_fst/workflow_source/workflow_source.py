from gwf import Workflow
from gwf.workflow import collect
from workflow_templates import *
import glob, yaml, os

def fst(config_file = glob.glob('*config.y*ml')[0]):
    """
    Template: description
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account']
    SPECIES_NAME: str = config['species_name']
    SAMPLE_LIST: list = config['sample_list']
    REFERENCE_GENOME: str = config['reference_genome_path']
    WORK_DIR: str = config['working_directory_path']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    partitions = partition_chrom(parse_fasta=parse_fasta(REFERENCE_GENOME), size=200000)
    
    top_dir = '{work_dir}/07_fst/{species_name}'.format(work_dir=WORK_DIR, species_name=SPECIES_NAME.replace(' ', '_'))
    os.makedirs(top_dir, exist_ok=True)
    
    return gwf