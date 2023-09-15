from gwf import AnonymousTarget
import os

def get_reference_name(idx, target):
    """Function to isolate the directory name for each reference genome"""
    filename = os.path.basename(os.path.dirname(target.inputs['path']))
    return 'Index_{}'.format(filename)

def index_reference_genome_loop(reference_genome: str):
    """
    Template: Index reference genomes in path with :script:`bwa index` and :script:`samtools faidx`.
    
    Template I/O::

        inputs = {'path': reference_genome}
        outputs = {'path': [reference_genome.amb,
                            reference_genome.ann,
                            reference_genome.pac,
                            reference_genome.bwt,
                            reference_genome.sa,
                            reference_genome.fai]}
    
    :param str reference_genome:
        Reference or draft genome in `FASTA`format
    """
    inputs = {'path': reference_genome}
    outputs = {'path': ['{}.amb'.format(reference_genome),
                        '{}.ann'.format(reference_genome),
                        '{}.pac'.format(reference_genome),
                        '{}.bwt'.format(reference_genome),
                        '{}.sa'.format(reference_genome),
                        '{}.fai'.format(reference_genome)]}
    protect = outputs['path']
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '02:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate omni_c
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    bwa index \
        -p {reference_genome} \
        {reference_genome}
    
    samtools faidx \
        -o {reference_genome}.fai.prog \
        {reference_genome}
    
    mv {reference_genome}.fai.prog {reference_genome}.fai

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(reference_genome=reference_genome)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, protect=protect, spec=spec)