from gwf import Workflow
from gwf.workflow import collect
from workflow_templates import *
import glob, yaml, os

def fst_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Creates :format:`mpileup`file and corresponding :format:`sync` file for a list of :format:`BAM` alignment files.
    Partitions reference genome and creates a small :format:`mpileup` file for each partition and for each of these creates a small :format:`sync` file.
    When files for each partition has been made concatenates all parts into a single :format:`mpileup` file and a single :format:`sync` file.
    
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
    OUTPUT_DIR: str = config['output_directory_path']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    # Partitions reference genome
    partitions = partition_chrom(parse_fasta=parse_fasta(REFERENCE_GENOME), size=200000)
    
    # Directory setup
    top_dir = '{work_dir}/07_fst/{species_name}'.format(work_dir=WORK_DIR, species_name=SPECIES_NAME.replace(' ', '_'))
    os.makedirs(top_dir, exist_ok=True)
    output_dir = '{output_path}/{species_abbr}'.format(output_path=OUTPUT_DIR, species_abbr=species_abbreviation(SPECIES_NAME))
    os.makedirs(output_dir, exist_ok=True)

    mpileup = gwf.map(
        name=name_mpileup,
        template_func=mpileup_parts,
        inputs=partitions,
        extra={'bam_files': SAMPLE_LIST,
               'reference_genome': REFERENCE_GENOME,
               'species_name': SPECIES_NAME,
               'output_directory': top_dir}
    )

    sync = gwf.map(
        name=name_sync,
        template_func=mpileup2sync,
        inputs=collect(mpileup.outputs, ['mpileup'])['mpileups']
    )

    concat_mpileup = gwf.target_from_template(
        name='concat_mpileup',
        template=concat(
            files=collect(mpileup.outputs, ['mpileup'])['mpileups'],
            output_name='{}'.format(species_abbreviation(SPECIES_NAME)),
            output_directory=output_dir
        )
    )

    concat_sync = gwf.target_from_template(
        name='concat_sync',
        template=concat(
            files=collect(sync.outputs, ['sync'])['syncs'],
            output_name='{}'.format(species_abbreviation(SPECIES_NAME)),
            output_directory=output_dir
        )
    )

    return gwf