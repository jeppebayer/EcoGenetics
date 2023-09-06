from gwf import Workflow, AnonymousTarget
import os

def get_reference_name(idx, target):
    """Function to isolate the directory name for each reference genome"""
    filename = os.path.basename(os.path.dirname(target.inputs['path']))
    return 'Index_{}'.format(filename)

def index_reference_genome_loop(reference_genome):
    """Template for indexing all reference genomes
    in path with 'bwa index' and 'samtools faidx'."""
    inputs = {'path': reference_genome}
    outputs = {'path': ['{}.amb'.format(reference_genome),
                        '{}.ann'.format(reference_genome),
                        '{}.pac'.format(reference_genome),
                        '{}.bwt'.format(reference_genome),
                        '{}.sa'.format(reference_genome),
                        '{}.fai'.format(reference_genome)]}
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '02:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate data_prep
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    bwa index \
        -p {reference_genome} \
        {reference_genome}
    
    samtools faidx \
        -o {reference_genome}.fai \
        {reference_genome}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(reference_genome=reference_genome)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)