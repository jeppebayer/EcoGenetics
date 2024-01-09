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

    top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}'
    os.makedirs(top_dir, exist_ok=True)

    index_reference = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_index_reference',
        template=index_reference_genome_loop(
            reference_genome=DRAFT_GENOME)
    )

    align_hic_data = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_HiC_alignment',
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

    no_curation_conversion = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_convert_no_cu',
        template=alignment_conversion_no_manual_curation(
            hic_bin=scaffolding.outputs['bin'],
            scaffolds_final_agp=scaffolding.outputs['final'][0],
            draft_assembly_fai_index=index_reference.outputs['path'][5],
            output_directory_path=top_dir,
            species_name=SPECIES_NAME
        )
    )

    curation_conversion = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_convert_cu',
        template=alignment_conversion_manual_curation(
            hic_bin=scaffolding.outputs['bin'],
            scaffolds_final_agp=scaffolding.outputs['final'][0],
            draft_assembly_fai_index=index_reference.outputs['path'][5],
            output_directory_path=top_dir,
            species_name=SPECIES_NAME
        )
    )

    chromosome_sizes = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_chrom_sizes',
        template=chrom_sizes(
            scaffolds_final_fa=scaffolding.outputs['final'][1]
        )
    )

    matrix_no_curation = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_matrix_no_cu',
        template=contact_matrix_no_manual_curation(
            alignments_sorted=no_curation_conversion.outputs['sorted'],
            chrom_sizes=chromosome_sizes.outputs['sizes'],
            output_directory_path=top_dir,
            species_name=SPECIES_NAME
        )
    )

    # matrix_curation = gwf.target_from_template(
    #     name=f'{species_abbreviation(SPECIES_NAME)}_matrix_cu',
    #     template=contact_matrix_manual_curation(
    #         JBAT_text=curation_conversion.outputs['jbat'][0],
    #         chrom_sizes=chromosome_sizes.outputs['sizes'],
    #         output_directory_path=top_dir,
    #         species_name=SPECIES_NAME
    #     )
    # )

    matrix_curation = gwf.target_from_template(
        name=f'{species_abbreviation(SPECIES_NAME)}_matrix_cu',
        template=contact_matrix_manual_curation(
            JBAT_text=curation_conversion.outputs['jbat'][0],
            JBAT_log=curation_conversion.outputs['jbat'][4],
            output_directory_path=top_dir,
            species_name=SPECIES_NAME
        )
    )

    return gwf