from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

# TODO Make create_vcf_workflow more streamlined like pooled_species_vcf_workflow

# Workflow for creating VCF file of from alignments of one to many individually sequenced specimen from a population or pool-sequence specimen.
# The reference genome is used to create a list of all 'chromosomes' and their lengths, which are to partitioned into smaller sequences of 500kbs.
# The partitioned list is used to do variant calling in parallel on thousands of segments to reduce the overall time needed.
# A BED file can be used to avoid repeat regions and further decrease the time needed. By default all non-SNP variants, e.g. indels, are removed.
# Once all VCF for all segments have been created they are merged into one file on which SNPeff is run. 
def create_vcf_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow for creating a single :format:`VCF` file of either pooled or individual sequence data.

    The workflow will automatically any file in the execution directory following the naming convention: *config.yaml
    
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
    REPEAT_REGIONS: str | None = config['repeat_regions_bed']
    TYPE: str = config['sequencing_type'].lower()
    WORK_DIR: str = config['working_directory']
    PLOIDY: int = config['ploidy']
    BESTN: int = config['bestn']

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------

    # Creates a version of the species name where space characters have been replaced with underscores
    species_name_nospace = SPECIES_NAME.replace(' ', '_')
    # Creates a list of dictionaries where each dictionary contains the name of a region, start and end position for each segment within a region and its chronological number
    partitions = partition_chrom(parse_fasta(REFERENCE_GENOME))
    # Takes the sample list and converts it into a string including all samples
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
                regions = collect(vcf_per_chr_individual.outputs, ['region'])['regions'],
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
                regions = collect(vcf_per_chr_pooled.outputs, ['region'])['regions'],
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

def pooled_species_vcf_workflow(config_file = glob.glob('*config.y*ml')[0]):
    """
    Workflow for creating a single :format:`VCF` file containing data on all pooled samples within a species.

    The workflow will automatically any file in the execution directory following the naming convention: *config.yaml
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------

    config = yaml.safe_load(open(config_file))
    ACCOUNT: str = config['account']
    SPECIES_NAME: str = config['species_name']
    SPECIES_DIR: list = config['species_directory_path']
    REFERENCE_GENOME: str = config['reference_genome_path']
    REPEAT_REGIONS:str | None = config['repeat_regions_bed']
    SNPEFF: bool = config['snpeff']
    WORK_DIR: str = config['working_directory']
    PLOIDY: int = config['ploidy']
    BESTN: int = config['bestn']

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------

    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    sample_list = gwf.glob('{}/*/*_filtered.bam'.format(SPECIES_DIR))
    sample_string = ' -b '.join(sample_list)
    partitions = partition_chrom(parse_fasta=parse_fasta(REFERENCE_GENOME), size=5000)

    # Directory setup
    top_dir = '{work_dir}/03_initial_analysis_files/{species_name}/all_populations'.format(work_dir=WORK_DIR, species_name=SPECIES_NAME.replace(' ', '_'))
    os.makedirs(top_dir, exist_ok=True)

    if REPEAT_REGIONS == 'None':
        vcf_parts = gwf.map(
            template_func=vcf_per_chr_pooled_all_no_rep,
            inputs=partitions,
            name=name_vcf,
            extra={'reference_genome': REFERENCE_GENOME,
                'sample_list': sample_string,
                'working_directory': top_dir,
                'ploidy': PLOIDY,
                'bestn': BESTN})
    else:
        vcf_parts = gwf.map(
            template_func=vcf_per_chr_pooled_all_rep,
            inputs=partitions,
            name=name_vcf,
            extra={'reference_genome': REFERENCE_GENOME,
                'sample_list': sample_string,
                'repeat_regions': REPEAT_REGIONS,
                'working_directory': top_dir,
                'ploidy': PLOIDY,
                'bestn': BESTN})
    concatenate = gwf.target_from_template(
        name='VCF_concat_{}'.format(species_abbreviation(SPECIES_NAME)),
        template=concatenate_vcf(
            regions=collect(vcf_parts.outputs, ['region'])['regions'],
            output_directory=top_dir,
            species_name=SPECIES_NAME
        )
    )
    if SNPEFF is True:
        snpeff = gwf.target_from_template(
            name='snpEff_annotation_{}'.format(species_abbreviation(SPECIES_NAME)),
            template=snpEff_ann(
                vcf_file=concatenate.outputs['concat_vcf'],
                reference_genome_version=os.path.splitext(os.path.basename(REFERENCE_GENOME))[0],
                output_directory=top_dir
            )
        )
    return gwf

def vcf_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: description
    
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
    REPEAT_REGIONS: str | None = config['repeat_regions_bed']
    if REPEAT_REGIONS == 'None':
        REPEAT_REGIONS = None
    TYPE: str = config['sequencing_type']
    SNPEFF: bool = config['snpeff']
    if SNPEFF == 'False':
        SNPEFF = False
    elif SNPEFF == 'True':
        SNPEFF = True
    WORK_DIR: str = config['working_directory_path']
    PLOIDY: int = config['ploidy']
    BESTN: int = config['bestn']
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    partitions = partition_chrom_real(parse_fasta=parse_fasta(REFERENCE_GENOME), size=100000)

    # Directory setup
    # Defines and creates a species specific working directory with a folder for temporary files.
    if TYPE == 'pooled':
        # Pooled samples gets a 'sample folder'
        sample_name = os.path.basename(os.path.dirname(SAMPLE_LIST[0]))
        sample_dir = '/{}'.format(sample_name)
        output_directory = os.path.dirname(SAMPLE_LIST[0])
    elif TYPE == 'individual':
        # Individually sequenced samples are grouped in VCF files by species, so no need for 'sample folder'
        sample_name = species_abbreviation(SPECIES_NAME)
        sample_dir = ""
        output_directory = '{}/{}'.format(os.path.dirname(os.path.dirname(SAMPLE_LIST[0])), species_abbreviation(SPECIES_NAME))

    top_dir = '{work_dir}/03_initial_analysis_files/{species_name}{sample}'.format(work_dir=WORK_DIR, species_name=SPECIES_NAME.replace(' ', '_'), sample=sample_dir)
    os.makedirs(top_dir, exist_ok=True)
    
    if TYPE == 'individual':
        vcf_parts = gwf.map(
            template_func=vcf_per_chr_individual,
            inputs=partitions,
            name=vcf_name,
            extra={'reference_genome': REFERENCE_GENOME,
                'sample_list': SAMPLE_LIST,
                'sample_name': sample_name,
                'output_directory': top_dir,
                'repeat_regions': REPEAT_REGIONS,
                'ploidy': PLOIDY,
                'bestn': BESTN})
    elif TYPE == 'pooled':
        vcf_parts = gwf.map(
            template_func=vcf_per_chr_pooled,
            inputs=partitions,
            name=vcf_name,
            extra={'reference_genome': REFERENCE_GENOME,
                'sample_list': SAMPLE_LIST,
                'sample_name': sample_name,
                'output_directory': top_dir,
                'repeat_regions': REPEAT_REGIONS,
                'ploidy': PLOIDY,
                'bestn': BESTN})
        
    concatenate = gwf.target_from_template(
        name='VCF_concat_{}'.format(species_abbreviation(SPECIES_NAME)),
        template=concat_vcf(
            regions=collect(vcf_parts.outputs, ['region'])['regions'],
            output_name=sample_name,
            output_directory=output_directory
        )
    )
    if SNPEFF is True:
        snpeff = gwf.target_from_template(
            name='snpEff_annotation_{}'.format(species_abbreviation(SPECIES_NAME)),
            template=snpEff_ann(
                vcf_file=concatenate.outputs['concat_vcf'],
                reference_genome_version=os.path.splitext(os.path.basename(REFERENCE_GENOME))[0],
                output_directory=output_directory
            )
        )

    return gwf

########################## PoolSNP ##########################

def poolsnp_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Create :format:`VCF`files based on :script:`PoolSNP`
    
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
    sample_names: list = [os.path.basename(os.path.dirname(path)) for path in SAMPLE_LIST]
    REFERENCE_GENOME: str = config['reference_genome_path']
    MPILEUP: str = config['mpileup_path']
    WORKING_DIR: str = config['working_directory_path']
    OUTPUT_DIR: str = config['output_directory_path']
    MAXCOV: float = config['max_cov']
    MINCOV: int = config['min_cov']
    MINCOUNT: int = config['min_count']
    MINFREQ: float = config['min_freq']
    MISSFRAC: float = config['miss_frac']
    BASEQUAL: int = config['bq']
    ALLSITES: bool = config['all_sites']
    if ALLSITES == 'True':
        ALLSITES = 1,
    elif ALLSITES == 'False':
        ALLSITES = 0
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    
    # Partitions reference genome
    partitions = partition_chrom_real(parse_fasta=parse_fasta(REFERENCE_GENOME), size=200000)
    # Creates list of contigs
    contigs = [{'contig': contig['sequence_name']} for contig in parse_fasta(REFERENCE_GENOME)]

    top_dir = '{working_directory}/03_initial_analysis_files/{species_name}'.format(working_directory=WORKING_DIR, species_name=SPECIES_NAME.replace(' ', '_'))
    output_dir = '{output_directory}/{species_abbr}'.format(output_directory=OUTPUT_DIR, species_abbr=species_abbreviation(SPECIES_NAME))

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

    coverage_threshold = gwf.map(
        name=name_cov,
        template_func=max_cov,
        inputs=contigs,
        extra={'mpileup': concat_mpileup.outputs['concat_file'],
               'cutoff': MAXCOV,
               'output_directory': top_dir}
    )

    concat_coverage = gwf.target_from_template(
        name='concatenate_coverage',
        template=concat(
            files=collect(coverage_threshold.outputs, ['cutoff'])['cutoffs'],
            output_name='{species_abbr}-cov-{max_cov}'.format(species_abbr=species_abbreviation(SPECIES_NAME), max_cov=MAXCOV),
            output_directory=top_dir
        )
    )

    run_poolsnp = gwf.target_from_template(
        name='poolsnp',
        template=poolsnp(
            mpileup=MPILEUP,
            max_cov=concat_coverage.outputs['concat_file'],
            sample_list=sample_names,
            reference_genome=REFERENCE_GENOME,
            working_directory=top_dir,
            species_name=SPECIES_NAME,
            output_directory=output_dir,
            min_cov=MINCOV,
            min_count=MINCOUNT,
            min_freq=MINFREQ,
            miss_frac=MISSFRAC,
            bq=BASEQUAL,
            sites=ALLSITES
        )
    )

    return gwf