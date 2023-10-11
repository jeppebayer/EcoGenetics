from gwf import AnonymousTarget
import glob, os

def species_abbreviation(species_name: str):
    """Creates species abbreviation from species name.
    
    :param str species_name:
        Species name written as *genus* *species*"""
    genus, species = species_name.replace(' ', '_').split('_')
    genus = genus[0].upper() + genus[1:3]
    species = species[0].upper() + species[1:3]
    return genus + species

def parse_fasta(fasta_file: str):
    """Parses `FASTA` file returning all sequence names and lengths paired in a list of dictionaries.
    
    ::
    
        return [{'sequence_name': str, 'sequence_length': int}, ...]
    
    :param str fasta_file:
        Sequence file in `FASTA` format.
    """
    fasta_list = []
    seq_name = None
    length = 0
    with open(fasta_file, 'r') as fasta:
        for entry in fasta:
            entry = entry.strip()
            if entry.startswith(">"):
                if seq_name:
                    fasta_list.append({'sequence_name': seq_name, 'sequence_length': length})
                    length = 0
                entry = entry.split(" ", 1)
                seq_name = entry[0][1:]
            else:
                length += len(entry)
        fasta_list.append({'sequence_name': seq_name, 'sequence_length': length})
    return fasta_list

def partition_chrom(parse_fasta: list, size: int = 500000):
    """
    Partitions `FASTA` file parsed with **parse_fasta**.
    
    Uses the list of dictionaries from **parse_fasta** to creates a list of dictionaries
    containing with partition number, sequence name, start and end postion (0 based).
    By default the partition size is 500kbs.
    
    ::
    
        return [{'num': int, 'region': str, 'start': int, 'end': int}, ...]
    
    :param list parse_fasta:
        List of dictionaries produced by **parse_fasta**.
    :param int size:
        Size of partitions. Default 500kb.
    """
    chrom_partition = []
    num = 1
    for chrom in parse_fasta:
        whole_chunks = chrom['sequence_length'] // size
        partial_chunk = chrom['sequence_length'] - whole_chunks * size
        start = 0
        for chunk in range(whole_chunks):
            end = start + size
            chrom_partition.append({'num': num, 'region': chrom['sequence_name'], 'start': start, 'end': end})
            start = end
            num += 1
        chrom_partition.append({'num': num, 'region': chrom['sequence_name'], 'start': start, 'end': start + partial_chunk})
        num += 1
    return chrom_partition

def make_bed(partitions: list, output_directory: str = '.'):
    """
    Creates :format:`BED` file for each partition from **partition_chrom**.

    Uses the partitions from **partition_chrom** to create a series of :format:`BED` files for each partition.
    The :format:`BED`files will be named: 'num'.bed, where 'num' corresponds to the partitions number.
    Each :format:`BED`files contains one line with three tab-separated columns: 'chromosome', 'start', 'end', with coordinates being 0-based.
    Returns list of all files created.

    ::

        return ['1.bed', '2.bed',... 'n.bed']

    :param list partitions:
        List of dictionaries produced by **partitions_chrom**
    :param str output_directory:
        Directory for output :format:`BED`files.
    """
    output_directory = '{}/tmp'.format(output_directory)
    partition_w_bed = []
    for partition in partitions:
        partition_w_bed.append({'num': partition['num'], 'region': partition['region'], 'start': partition['start'], 'end': partition['end'], 'bed_file': '{dir_path}/{num}.bed'.format(dir_path=output_directory, num=partition['num'])})
    return partition_w_bed

def name_mpileup(idx, target):
    return 'mpileup_{idx}'.format(idx=idx+1)

def name_sync(idx, target):
    return 'sync_{idx}'.format(idx=idx+1)

def bed_files(refrence_genome: str, size: int, output_directory: str = '.', script: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/07_fst/workflow_source/make_bed.py'):
    """
    Template: Makes a series of :format:`BED`files according to chromosome.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    output_directory='{}/tmp'.format(output_directory)
    inputs = {'reference': refrence_genome}
    outputs = {'beds': make_bed(partitions=partition_chrom(parse_fasta=parse_fasta(refrence_genome), size=size), output_directory=output_directory)}
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '00:30:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate data_prep
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_dir} ] || mkdir -p {output_dir}

    python {script} \
        {reference} \
        {size} \
        {output_dir}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(script=script, reference=refrence_genome, size=size, output_dir=output_directory)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mpileup_parts(bam_files: list, reference_genome: str, species_name: str, region: str, num: int, start: int, end: int, output_directory: str = '.'):
    """
    Template: Create :format:`mpileup` files for each partition of reference genome from multiple :format:`BAM`files using :script:`samtools mpileup`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    output_directory = '{}/tmp'.format(output_directory)
    inputs = {'bam_files': bam_files,
              'reference': reference_genome}
    file_name = '{output_path}/{abbr}_{num}_{chrom}'.format(output_path=output_directory, abbr=species_abbreviation(species_name), num=num, chrom=region)
    outputs = {'mpileup': '{file_name}.mpileup'.format(file_name=file_name)}
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '10:00:00'
    }
    bam_string = ' '.join(bam_files)
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate data_prep
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    sleep 5m

    [ -d {output_dir} ] || mkdir -p {output_dir}

    echo -e '{chromosome}\t{start}\t{end}' > {output_dir}/{num}.bed
    
    samtools mpileup \
        --max-depth 0 \
        --fasta-ref {reference} \
        --positions {output_dir}/{num}.bed \
        --min-BQ 0 \
        --region {chromosome} \
        --output {file_name}.prog.mpileup \
        {bam_files}
    
    mv {file_name}.prog.mpileup {mpileup}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(output_dir=output_directory, chromosome=region, start=start, end=end, num=num, reference=reference_genome, file_name=file_name, bam_files=bam_string, mpileup=outputs['mpileup'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mpileup2sync(mpileup_file: str, output_directory: str = None, mpileup2sync: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/07_fst/workflow_source/popoolation2_v1.201/mpileup2sync.pl'):
    """
    Template: Makes a :format:`sync` file for a corresponding :format:`mpileup` file using :script:`popoolation2`'s :script:`mpileup2sync.pl`
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    if output_directory is None:
        output_directory = os.path.dirname(mpileup_file)
    inputs = {'mpileup': mpileup_file}
    file_name = '{output_path}/{filename}'.format(output_path=output_directory, filename=os.path.splitext(os.path.basename(mpileup_file))[0])
    outputs = {'sync': '{}.sync'.format(file_name)}
    options = {
        'cores': 2,
        'memory': '16g',
        'walltime': '06:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate data_prep
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    perl {mpileup2sync} \
        --input {mpileup} \
        --output {file_name}.prog.sync \
        --fastq-type sanger
    
    mv {file_name}.prog.sync {sync}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(mpileup2sync=mpileup2sync, mpileup=mpileup_file, file_name=file_name, sync=outputs['sync'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def concat(files: list, output_name: str, output_directory: str = None):
    """
    Template: Name sorts and concatenates files.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'files': files}
    files.sort
    if output_directory is None:
        output_directory = os.path.dirname(files[0])
    file_name = '{output_path}/{output_name}'.format(output_path=output_directory, output_name=output_name)
    outputs = {'concat_file': '{file_name}{ext}'.format(file_name=file_name, ext=os.path.splitext(files[0])[1])}
    options = {
        'cores': 2,
        'memory': '16g',
        'walltime': '06:00:00'
    }
    protect = outputs['concat_file']
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate data_prep
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    cat \
        {sorted_files} \
        > {file_name}.prog{ext}
    
    mv {file_name}.prog{ext} {concat_file}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(sorted_files=' '.join(files), file_name=file_name, ext=os.path.splitext(files[0])[1], concat_file=outputs['concat_file'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, protect=protect, spec=spec)