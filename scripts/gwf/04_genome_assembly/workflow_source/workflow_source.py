from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

# Current workflow still being tested. Based on the Dovetail guide: https://omni-c.readthedocs.io/en/latest/index.html.
# First, Hi-C data is aligned to draft genome, with setting to map mates independently and for supposed optimal results with Omni-C proximity library reads.
# Second, identifies ligation events in the alignment file and records the outer-most (5') aligned base pair and their respective strands, further, each pair is assigned a type of event.
# Third, PCR and possible optical duplicates are removed.
# Fourth, a file containing all Omni-C pairs and a BAM file of the alignment is created. 
# Fifth, a Hi-C contat map is generated using the Omni-C pairs and a chrom.sizes file containing a all chromosomes in the draft genome and their respective lengths.
def hic_processing_workflow(config_file = glob.glob('*config.y*ml')[0]):
    
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------

    config=yaml.safe_load(open(config_file))
    ACCOUNT=config['account']
    SPECIES_NAME=config['species_name']
    HIC_DIRECTORY=config['HiC_directory_path']
    REFERENCE_GENOME=config['reference_genome_path']
    WORK_DIR=config['working_directory_path']

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )

    top_dir = '{work_dir}/04_genome_assembly/{species_name}/HiC'.format(work_dir=WORK_DIR, species_name=SPECIES_NAME.replace(' ', '_'))
    os.makedirs(top_dir, exist_ok=True)
    raw_dir = '{}/raw'.format(top_dir)
    os.makedirs(raw_dir, exist_ok=True)
    align_dir = '{}/aligned'.format(top_dir)
    os.makedirs(align_dir, exist_ok=True)
    read1 = '{}'.format(glob.glob('{}/*_1.fq*'.format(HIC_DIRECTORY))[0])
    read2 = '{}'.format(glob.glob('{}/*_2.fq*'.format(HIC_DIRECTORY))[0])
    if os.path.splitext(read1)[1] == '.gz':
        symlink1 = '{raw_dir}/{abbr}_R1.fastq.gz'.format(raw_dir=raw_dir, abbr=species_abbreviation(SPECIES_NAME))
        symlink2 = '{raw_dir}/{abbr}_R2.fastq.gz'.format(raw_dir=raw_dir, abbr=species_abbreviation(SPECIES_NAME))
    else:
        symlink1 = '{raw_dir}/{abbr}_R1.fastq'.format(raw_dir=raw_dir, abbr=species_abbreviation(SPECIES_NAME))
        symlink2 = '{raw_dir}/{abbr}_R2.fastq'.format(raw_dir=raw_dir, abbr=species_abbreviation(SPECIES_NAME))
    if not os.path.islink(symlink1):
        os.symlink(read1, symlink1)
    if not os.path.islink(symlink2):
        os.symlink(read2, symlink2)


    index_reference = gwf.target_from_template(
        name='{}_index_reference'.format(species_abbreviation(SPECIES_NAME)),
        template=index_reference_genome_loop(
            reference_genome=REFERENCE_GENOME)
    )

    # directories = gwf.target_from_template(
    #     name='{}_directory_setup'.format(species_abbreviation(SPECIES_NAME)),
    #     template=directory_setup(
    #         hic_directory=HIC_DIRECTORY,
    #         working_directory=WORK_DIR,
    #         species_name=SPECIES_NAME
    #     )
    # )

    chrom_size = gwf.target_from_template(
        name='{}_chrom_sizes'.format(species_abbreviation(SPECIES_NAME)),
        template=get_chrom_sizes(
            fasta_file=REFERENCE_GENOME,
            output_directory=raw_dir
        )
    )

    align = gwf.target_from_template(
        name='{}_alignment'.format(species_abbreviation(SPECIES_NAME)),
        template=hic_align(
            read1=symlink1,
            read2=symlink2,
            draft_genome=REFERENCE_GENOME,
            output_directory=align_dir,
            species_name=SPECIES_NAME
        )
    )
    
    ligation = gwf.target_from_template(
        name='{}_ligation_events'.format(species_abbreviation(SPECIES_NAME)),
        template=ligation_events(
            sam_file=align.outputs['sam'],
            chrom_sizes=chrom_size.outputs['chrom_sizes'],
            species_name=SPECIES_NAME
        )
    )

    de_dup = gwf.target_from_template(
        name='{}_de_dup'.format(species_abbreviation(SPECIES_NAME)),
        template=de_duplicate(
            pairsam_file=ligation.outputs['pairsam'])
    )

    split = gwf.target_from_template(
        name='{}_splitting'.format(species_abbreviation(SPECIES_NAME)),
        template=split_pairsam(
            pairsam_file=de_dup.outputs['dedup']
        )
    )

    juicer = gwf.target_from_template(
        name='{}_juicer_tools'.format(species_abbreviation(SPECIES_NAME)),
        template=juicer_hic_matrix(
            pairs_file=split.outputs['pairs'],
            reference_genome=REFERENCE_GENOME,
        )
    )

    return gwf