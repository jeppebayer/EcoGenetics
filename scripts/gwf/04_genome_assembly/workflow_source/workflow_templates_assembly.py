from gwf import AnonymousTarget
import glob, os

def species_abbreviation(species_name: str) -> str:
    """Creates species abbreviation from species name.
    
    :param str species_name:
        Species name written as *genus* *species*"""
    genus, species = species_name.replace(' ', '_').split('_')
    genus = genus[0].upper() + genus[1:3]
    species = species[0].upper() + species[1:3]
    return genus + species

# Draft genome step 1 - Remove potential leftover adapters
def hifiadapterfilt(pacbio_hifi_reads: str, HiFiAdapterFilt_path: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/04_genome_assembly/workflow_source/HiFiAdapterFilt'):
    """
    Template: Removes remaining adapters from PacBio HiFi reads using :script:`HiFiAdapterFilt`.

    Template I/O::

        inputs = {'pb_hifi': pacbio_hifi_reads}

        outputs = {'filt': *.filt.fastq.gz,
                   'cont': *.contaminant.blastout,
                   'block': *.blocklist,
                   'stats': *.stats}
    
    :param str pacbio_hifi_reads:
        File containing PacBio HiFi reads (.bam | .fastq | .fastq.gz | .fq | .fq.gz).
    :param str HiFiAdapterFilt_path:
        Path to :script:`HiFiAdapterFilt` directory.
    """
    output_directory = f'{os.path.dirname(pacbio_hifi_reads)}/HiFiAdapterFilt'
    inputs = {'pb_hifi': pacbio_hifi_reads}
    if pacbio_hifi_reads.endswith('.gz'):
        prefix = os.path.splitext(os.path.splitext(os.path.basename(pacbio_hifi_reads))[0])[0]
        ext = f'{os.path.splitext(os.path.splitext(os.path.basename(pacbio_hifi_reads))[0])[1]}.gz'
    else:
        prefix = os.path.splitext(os.path.basename(pacbio_hifi_reads))[0]
        ext = os.path.splitext(os.path.basename(pacbio_hifi_reads))[1]
    outputs = {'filt': f'{output_directory}/{prefix}.filt.fastq.gz',
               'cont': f'{output_directory}/{prefix}.contaminant.blastout',
               'block': f'{output_directory}/{prefix}.blocklist',
               'stats': f'{output_directory}/{prefix}.stats'}
    protect = [outputs['filt'], outputs['cont'], outputs['block'], outputs['stats']]
    options = {
        'cores': 10,
        'memory': '100g',
        'walltime': '02:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    export PATH=$PATH:{HiFiAdapterFilt_path}
    export PATH=$PATH:{HiFiAdapterFilt_path}/DB

    [ -d {output_directory} ] || mkdir -p {output_directory}

    ln -s {pacbio_hifi_reads} "$(dirname {pacbio_hifi_reads})"/prog.{ext}

    bash {HiFiAdapterFilt_path}/pbadapterfilt.sh \
        -t {options['cores']} \
        -p prog \
        -o {output_directory}
        
    mv {output_directory}/prog.filt.fastq.gz {outputs['filt']}
    mv {output_directory}/prog.contaminant.blastout {outputs['cont']}
    mv {output_directory}/prog.blocklist {outputs['block']}
    mv {output_directory}/prog.stats {outputs['stats']}
    rm "$(dirname {pacbio_hifi_reads})"/prog.{ext}
        
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

# Draft genome step 2 (N50 and k'mer analysis)
#TODO fix this kmer shas
def kmer_analysis(genome_file: str, output_path: str, k: str =27, cannonical: bool = False, cutoff: int = 1000):
    """
    Template: Count number of k'mers in genome file using :script:`jellyfish count`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    output_path = '{}/kmer_analysis'.format(output_path)
    inputs = {'genome': genome_file}
    file_name = '{}'.format(os.path.basename(inputs['genome']))
    outputs = {'counts': '{output_path}/{file_name}.kmer_count_{k}'.format(output_path=output_path, file_name=file_name, k=k),
               'histo': '{output_path}/{file_name}.kmer_histo_{k}'.format(output_path=output_path, file_name=file_name, k=k),
               'genomescope': ['{output_path}/genomescope2_{k}.{file_name}/{k}_linear_plot.png'.format(output_path=output_path, k=k, file_name=file_name),
                               '{output_path}/genomescope2_{k}.{file_name}/{k}_log_plot.png'.format(output_path=output_path, k=k, file_name=file_name),
                               '{output_path}/genomescope2_{k}.{file_name}/{k}_model.txt'.format(output_path=output_path, k=k, file_name=file_name),
                               '{output_path}/genomescope2_{k}.{file_name}/{k}_progress.txt'.format(output_path=output_path, k=k, file_name=file_name),
                               '{output_path}/genomescope2_{k}.{file_name}/{k}_summary.txt'.format(output_path=output_path, k=k, file_name=file_name),
                               '{output_path}/genomescope2_{k}.{file_name}/{k}_transformed_linear_plot.png'.format(output_path=output_path, k=k, file_name=file_name),
                               '{output_path}/genomescope2_{k}.{file_name}/{k}_transformed_log_plot.png'.format(output_path=output_path, k=k, file_name=file_name)]}
    options = {
        'cores': 32,
        'memory': '480g',
        'walltime': '12:00:00'
    }
    if inputs['genome'].endswith('.gz'):
        input_file = '(zcat {})'.format(inputs['genome'])
    else:
        input_file = inputs['genome']
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {kmer_dir} ] || mkdir -m 775 {kmer_dir}

    jellyfish count \
        -t {cores} \
        -s 8G \
        -m {k} \
        -o {kmer_dir}/{file_name}.prog.kmer_count_{k} \
        {input_file}
    
    mv {kmer_dir}/{file_name}.prog.kmer_count_{k} {count_file}
    
    jellyfish histo \
        -t {cores} \
        -h {max_count} \
        -o {kmer_dir}/{file_name}.prog.kmer_histo_{k} \
        {count_file}

    mv {kmer_dir}/{file_name}.prog.kmer_histo_{k} {histo_file}

    [ -d {kmer_dir}/genomescope2_{k}.{filename} ] || mkdir -m 775 {kmer_dir}/genomescope2_{k}.{filename}

    genomescope2 \
        -i {histo_file} \
        -p 2 \
        -o {kmer_dir}/genomescope2_{k}.{filename} \
        -k {k} \
        -n {k} \

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format()
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def sequence_stats(sequence_file: str, output_directory_path: str,
                   get_lengths: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/bin/get_length.py',
                   n50: str = "/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/bin/n50.py"):
    """
    Template: Creates standard sequence statistics; N50
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'sequence_file': sequence_file}
    outputs = {'stats': f'{output_directory_path}/assembly_stats/{os.path.basename(sequence_file)}.stats'}
    protect = [outputs['stats']]
    options = {
        'cores': 2,
        'memory': '12g',
        'walltime': '03:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory_path}/assembly_stats ] || mkdir -p {output_directory_path}/assembly_stats

    python {get_lengths} {sequence_file} | python {n50} - {output_directory_path}/assembly_stats/{os.path.basename(sequence_file)}.prog.stats
    
    mv {output_directory_path}/assembly_stats/{os.path.basename(sequence_file)}.prog.stats {outputs['stats']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

# Draft genome step 3 (hifiasm and gfa to fasta conversion + new stats)
def hifiasm(hifi_sequence_file: str, output_directory_path: str, species_name: str):
    """
    Template: Create initial assembly of genome using hifiasm.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'hifi': hifi_sequence_file}
    outputs = {'primary': f'{output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.bp.p_ctg.gfa',
               'misc': [f'{output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.a_ctg.gfa',
                        f'{output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.r_utg.gfa',
                        f'{output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.p_ctg.gfa']}
    protect = [outputs['primary']]
    options = {
        'cores': 32,
        'memory': '480g',
        'walltime': '10:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory_path}/hifiasm ] || mkdir -p {output_directory_path}/hifiasm

    hifiasm \
        -o {output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.prog \
        -s 0.1 \
        -l 3 \
        --primary \
        {hifi_sequence_file}
    
    mv {output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.bp.p_ctg.gfa {outputs['primary']}
    mv {output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.a_ctg.gfa {outputs['misc'][0]}
    mv {output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.r_utg.gfa {outputs['misc'][1]}
    mv {output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.p_ctg.gfa {outputs['misc'][2]}

    awk \
        'BEGIN{{FS='\t'}}
        {{if ($0 ~ /^S/)
            {{print ">"$2"\n"$3}}}}' \
        {outputs['primary']} \
    | fold \
    > {os.path.splitext(outputs['primary'])[0]}.fasta

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs,protect=protect, options=options, spec=spec)

def busco():
    """
    Template: Runs BUSCO analysis on genome assembly.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {}
    outputs = {}
    options = {
        'cores': 2,
        'memory': '16g',
        'walltime': '10:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    
    
    mv
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)