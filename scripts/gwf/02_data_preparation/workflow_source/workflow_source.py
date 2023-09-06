from gwf import Workflow
import os
import yaml
import glob
from workflow_templates import *

def create_data_prep_workflow(config_file = glob.glob('*config.yaml')[0]):

    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------

    config = yaml.safe_load(open(config_file))
    ACCOUNT = config['account']
    SPECIES_NAME = config['species_name']
    SAMPLE_PATH = config['sample_path']
    REFERENCE_GENOME = config['reference_genome_path']
    WORKDIR_PATH = config['working_directory']
    STATISTICS = config ['statistics']

    # Finds all non-empty directories in 'SAMPLE_PATH' and makes a list
    sample_directories = [sample_directory for sample_directory in glob.glob(SAMPLE_PATH + '/*')
                        if os.path.isdir(sample_directory) and any(os.scandir(sample_directory))]

    # Creates a dictionary where each directory in 'sample_directories' is a key which linked to a dictionary with two 
    # keys 'R1' and 'R2' which are each linked to a list of all read 1 and read 2 file for the sample directory respectively
    # masterlist = {i: {'R1': glob.glob(i + '/*_1.fq.gz'), 'R2': glob.glob(i + '/*_2.fq.gz')} for i in sample_directories}
    masterlist = [
        {'sample_dir': i,
        'sample_name': os.path.basename(i),
        'R1': glob.glob(i + '/*_1.fq.gz'),
        'R2': glob.glob(i + '/*_2.fq.gz'),
        'output_dir': '{work_dir}/{species_name}/{sample_name}'.format(work_dir=WORKDIR_PATH, species_name=SPECIES_NAME.replace(' ', '_'), sample_name=os.path.basename(i))
        } for i in sample_directories
        ]

    def name_func(idx, lst, text: str = ''):
        return '{sample}.{text}_{idx}'.format(sample=lst['sample_name'].replace('-', '_'), idx=idx, text=text)

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------

    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    for idx, lst in enumerate(masterlist, start=1):
        merge_runs = gwf.target_from_template(
            name=name_func(idx, lst, 'merge_runs'),
            template=merging_sequence_runs(
                sample_name=lst['sample_name'],
                R1=lst['R1'],
                R2=lst['R2'],
                output_dir=lst['output_dir'])
        )
        adapterremoval = gwf.target_from_template(
            name=name_func(idx, lst, 'AdapterRemoval'),
            template=AdapterRemoval_paired_end(
                read1=merge_runs.outputs['read1'],
                read2=merge_runs.outputs['read2'],
                sample_name=lst['sample_name'],
                output_dir=lst['output_dir']
            )
        )
        alignment_paired = gwf.target_from_template(
            name=name_func(idx, lst, 'align_paired'),
            template=alignment_contemp(
                sample_name=lst['sample_name'],
                reference_genome=REFERENCE_GENOME,
                read1=adapterremoval.outputs['pairs'][0],
                read2=adapterremoval.outputs['pairs'][1],
                output_dir=lst['output_dir']
            )
        )
        alignment_single = gwf.target_from_template(
            name=name_func(idx, lst, 'align_single'),
            template=alignment_contemp(
                sample_name=lst['sample_name'],
                reference_genome=REFERENCE_GENOME,
                read1=adapterremoval.outputs['single'],
                output_dir=lst['output_dir']
            )
        )
        merge_paired_single = gwf.target_from_template(
            name=name_func(idx, lst, 'merge_alignmnets'),
            template=merge_alignments(
                alignment1=alignment_paired.outputs['alignment'],
                alignment2=alignment_single.outputs['alignment'],
                sample_name=lst['sample_name'],
                output_dir=lst['output_dir']
            )
        )
        mark_duplicates = gwf.target_from_template(
            name=name_func(idx, lst, 'mark_duplicates'),
            template=markduplicates(
                alignment=merge_paired_single.outputs['merged'],
                output_dir=lst['output_dir']
            )
        )   
        unmapped = gwf.target_from_template(
            name=name_func(idx, lst, 'extract_unmapped'),
            template=extract_unmapped_reads(
                alignment=mark_duplicates.outputs['markdup'],
                output_dir=lst['output_dir']
            )
        )
        if STATISTICS:
            pre_filter_stats = gwf.target_from_template(
                name=name_func(idx, lst, 'pre_filter_stats'),
                template=samtools_stats(
                    alignment=mark_duplicates.outputs['markdup'],
                    output_dir=lst['output_dir']
                )
            )
        alignment_filtering = gwf.target_from_template(
            name=name_func(idx, lst, 'filter_alignment'),
            template=filter_alignment(
                alignment=mark_duplicates.outputs['markdup'],
                output_dir=lst['sample_dir']
            )
        )
        if STATISTICS:
            post_filter_stats = gwf.target_from_template(
                name=name_func(idx, lst, 'post_filter_stats'),
                template=samtools_stats(
                    alignment=alignment_filtering.outputs['filtered_alignment'],
                    output_dir=lst['output_dir']
                )
            )
        qm = gwf.target_from_template(
            name=name_func(idx, lst, 'qualimap'),
            template=qualimap(
                alignment=alignment_filtering.outputs['filtered_alignment'],
                output_dir=lst['sample_dir']
            )
        )
    return gwf