from gwf import AnonymousTarget
import glob, os

def cactus_seqFile(species: list, genomes: list, newick: str, output_path: str = '.'):
    with open('{}/seqFile.txt'.format(output_path), 'w') as seqfile:
        seqfile.write('{}\n\n'.format(newick))
        for idx, name in enumerate(species):
            seqfile.write('{name}\t{path}\n'.format(name=name, path=genomes[idx]))
    return '{}/seqFile.txt'.format(output_path)

def run_cactus(seqfile: str, file_prefix: str = 'msa', output_path: str = None):
    """
    Template: Uses :script:`cactus` for multiple sequence alignment of whole genomes.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    if output_path is None:
        output_path = os.path.dirname(seqfile)
    inputs = {'seqfile': seqfile}
    outputs = {'hal': '{output_path}/{prefix}.hal'.format(output_path=output_path, prefix=file_prefix)}
    options = {
        'cores': 32,
        'memory': '128g',
        'walltime': '48:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate cactus
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_path}/tmp ] || mkdir -m 775 {output_path}/tmp

    cd /faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/06_genome_alignment/workflow_source/cactus-bin-v2.6.7
    source cactus_env/bin/activate
    cactus \
        {output_path}/tmp \
        {seqfile} \
        {file_name}.prog.hal

    mv {file_name}.prog.hal {outputHal}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(output_path=output_path, seqfile=seqfile, file_name='{}/{}'.format(output_path, file_prefix), outputHal=outputs['hal'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)