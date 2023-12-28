from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates_assembly import *

def genome_assembly_workflow(config_file = glob.glob('*config.y*ml')[0]):
    
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------

    config=yaml.safe_load(open(config_file))
    ACCOUNT: str=config['account']
    SPECIES_NAME: str=config['species_name']
    HIC_READ1: str=config['HiC_read_1_path']
    HIC_READ2: str=config['HiC_read_2_path']
    DRAFT_GENOME: str=config['draft_genome_path']
    WORK_DIR: str=config['working_directory_path']

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )

    top_dir = f'{WORK_DIR}/04_genome_assembly/{SPECIES_NAME.replace(" ", "_")}'
    os.makedirs(top_dir, exist_ok=True)

    index_reference = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_index_reference',
        template=index_reference_genome_loop(
            reference_genome=DRAFT_GENOME)
    )

    align_hic_data = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_hic_alignment',
        template=hic_align(
            hic_read1=HIC_READ1,
            hic_read2=HIC_READ2,
            draft_genome=DRAFT_GENOME,
            indices=index_reference.outputs['path'],
            output_directory_path=top_dir,
            species_name=SPECIES_NAME
        )
    )

    duplicates = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_mark_duplicates',
        template=mark_duplicates(
            bam_file=align_hic_data.outputs['bam'],
            output_directory_path=top_dir
        )
    )

    scaffolding = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_HiC_scaffolding',
        template=hic_scaffolding(
            draft_genome=DRAFT_GENOME,
            hic_to_draft_bam=duplicates.outputs['markdup'],
            output_directory_path=top_dir,
            species_name=SPECIES_NAME
        )
    )

    return gwf