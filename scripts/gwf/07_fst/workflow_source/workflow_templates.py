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

def name_mpileup(idx, target):
    return 'mpileup_{idx}'.format(idx=idx+1)

def name_sync(idx, target):
    return 'sync_{idx}'.format(idx=idx+1)

def mpileup_parts(bam_files: list, reference_genome: str, species_name: str, region: str, num: int, start: int, end: int, output_directory: str = '.'):
    """
    Template: Create :format:`mpileup` files for each partition of reference genome from multiple :format:`BAM` files using :script:`samtools mpileup`.
    
    Template I/O::
    
        inputs = {'bam_files': bam_files,
                  'reference': reference_genome}
        outputs = {'mpileup': *.mpileup}
    
    :param list bam_files:
        List of all :format:`BAM` files to be included in :format:`mpileup` file.
    :param str reference_genome:
        Path to genome reference file in `FASTA`format.
    :param str species_name:
        Name of species being worked on.
    :param str region:
        Name of chromosome from **partition_chrom**.
    :param int num:
        Partition number from **partition_chrom**.
    :param int start:
        Start position from **partition_chrom**.
    :param int end:
        End position from **partition_chrom**.
    :param str output_directory:
        Path to desired output directory. Defaults to current directory. Always creates 'tmp' directory in output directory.
    """
    output_directory = '{}/tmp'.format(output_directory)
    file_name = '{output_path}/{abbr}_{num}_{chrom}'.format(output_path=output_directory, abbr=species_abbreviation(species_name), num=num, chrom=region)
    bam_string = ' '.join(bam_files)
    inputs = {'bam_files': bam_files,
              'reference': reference_genome}
    outputs = {'mpileup': '{file_name}.mpileup'.format(file_name=file_name)}
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '10:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate data_prep
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    sleep 2m

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
    
        inputs = {'mpileup': mpileup_file}
        outputs = {'sync': *.sync}
    
    :param str mpileup_file:
        Input :format:`mpileup` file.
    :param str output_directory:
        Desired output directory for :format:`sync` file. Defaults to directory of mpileup_file.
    :param str mpileup2sync:
        Path to :script:`mpile2sync.pl`.
    """
    if output_directory is None:
        output_directory = os.path.dirname(mpileup_file)
    file_name = '{output_path}/{filename}'.format(output_path=output_directory, filename=os.path.splitext(os.path.basename(mpileup_file))[0])
    inputs = {'mpileup': mpileup_file}
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

def concat(files: list, output_name: str, output_directory: str = None, compress: bool = False):
    """
    Template: Name-sorts and concatenates files. Optionally compresses output using :script:`gzip`.
    
    Template I/O::
    
        inputs = {'files': files}
        outputs = {'concat_file': output_name.ext | output_name.ext.gzip}
    
    :param list files:
        List containing files to concatenate.
    :param str output_name:
        Desired name of output file, no extension.
    :param str output_directory:
        Path to output directory. Default uses directory of first file in file_list.
    :param bool compress:
        Bool indicating whether the output file should be compressed or not.
    """
    files.sort
    if output_directory is None:
        output_directory = os.path.dirname(files[0])
    file_name = '{output_path}/{output_name}'.format(output_path=output_directory, output_name=output_name)
    inputs = {'files': files}
    if compress:
        outputs = {'concat_file': '{file_name}{ext}.gzip'.format(file_name=file_name, ext=os.path.splitext(files[0])[1])}
    else:
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
    
    if [ {compress} == 'False' ]; then
        cat \
            {sorted_files} \
            > {file_name}.prog{ext}
        
        mv {file_name}.prog{ext} {concat_file}
    else
        cat \
            {sorted_files} \
        | gzip \
            -c \
            - \
            > {file_name}.prog{ext}.gzip
        
        mv {file_name}.prog{ext}.gzip {concat_file}
    fi

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(compress=compress, sorted_files=' '.join(files), file_name=file_name, ext=os.path.splitext(files[0])[1], concat_file=outputs['concat_file'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, protect=protect, spec=spec)