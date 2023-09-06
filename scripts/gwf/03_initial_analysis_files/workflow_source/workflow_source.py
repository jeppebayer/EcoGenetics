from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def create_vcf_workflow(config_file = glob.glob('*config.yaml')[0]):
    
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------

    config=yaml.safe_load(open(config_file))
    ACCOUNT=config['account']
    SPECIES_NAME=config['species_name']
    SAMPLE_LIST=config['sample_list']
    REFERENCE_GENOME=config['reference_genome_path']
    REPEAT_REGIONS=config['repeat_regions_bed']
    TYPE=config['sequencing_type'].lower()
    WORK_DIR=config['working_directory']
    PLOIDY=config['ploidy']
    BESTN=config['bestn']

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------

    # Creates a version of the species name where space characters have been replaced with underscores
    species_name_nospace = SPECIES_NAME.replace(' ', '_')
    # Creates a list of dictionaries where each dictionary contains the name of a region, start and end position for each segment within a region and its chronological number
    partitions = partition_chrom(parse_fasta(REFERENCE_GENOME))
    # Takes the sample list and converts it into a multi-line string where each line was an entry in the list
    sample_string = ' -b '.join(SAMPLE_LIST)
    # Gets the path to snpEff config file
    snpEff_config = glob.glob(sys.prefix + '/share/snpeff*/snpEff.config')[0]
    # Gets the reference genome version based on the used reference genome
    reference_genome_version = os.path.splitext(os.path.basename(REFERENCE_GENOME))[0]

    # Defines and creates a species specific working directory with a folder for temporary files.
    if len(SAMPLE_LIST) == 1:
        # Pooled samples gets a 'sample folder'
        sample_name = os.path.splitext(os.path.basename(SAMPLE_LIST[0]))
        sample_dir = "/" + sample_name
    else:
        # Individually sequenced samples are grouped in VCF files by species, so no need for 'sample folder'
        sample_name = species_abbreviation(SPECIES_NAME)
        sample_dir = ""
        
    work_dir = WORK_DIR + "/03_initial_analysis_files/" + species_name_nospace + sample_dir
    temp_dir = work_dir + "/temp"
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )

    if TYPE == "individual":
        vcf_per_chr_individual = gwf.map(create_vcf_per_chr_individual,
                partitions,
                name = name_vcf,
                extra = {'reference_genome': REFERENCE_GENOME,
                    'sample_list': sample_string,
                    'repeat_regions': REPEAT_REGIONS,
                    'ploidy': PLOIDY,
                    'bestn': BESTN,
                    'temp_dir': temp_dir,
                    'sample_name': sample_name})
        gwf.target_from_template(
            name = 'VCF_concat_{}'.format(sample_name),
            template = vcf_concat(
                regions = collect(vcf_per_chr_individual.outputs, ['region']),
                sample_name = sample_name,
                temp_dir = temp_dir
            )
        )
    elif TYPE == "pooled":
        vcf_per_chr_pooled = gwf.map(create_vcf_per_chr_pooled,
                partitions,
                name = name_vcf,
                extra = {'reference_genome': REFERENCE_GENOME,
                    'sample_list': sample_string,
                    'repeat_regions': REPEAT_REGIONS,
                    'ploidy': PLOIDY,
                    'bestn': BESTN,
                    'temp_dir': temp_dir,
                    'sample_name': sample_name})
        gwf.target_from_template(
            name = 'VCF_concat_{}'.format(sample_name),
            template = vcf_concat(
                regions = collect(vcf_per_chr_pooled.outputs, ['region']),
                sample_name = sample_name,
                temp_dir = temp_dir
            )
        )

    gwf.target_from_template(
        name = 'snpEff_annotation_{}'.format(sample_name),
        template = snpEff_annotation(
            temp_dir = temp_dir,
            work_dir = work_dir,
            sample_name = sample_name,
            snpEff_config = snpEff_config,
            reference_genome_version = reference_genome_version
        )
    )

    return gwf