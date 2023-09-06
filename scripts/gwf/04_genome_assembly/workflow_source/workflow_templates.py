from gwf import AnonymousTarget
import glob, os

def species_abbreviation(species_name: str) -> str:
    """Function for creating species abbreviation from species name."""
    genus, species = species_name.split(maxsplit=1)
    genus = genus[0].upper() + genus[1:3]
    species = species[0].upper() + species[1:3]
    return genus + species

def directory_setup(working_directory: str, species_name: str, hic_directory: str):
    """Template for setting up folder structure Juicer jobs."""
    top_dir = '{work_dir}/04_genome_assembly/{species_name}/HiC'.format(work_dir=working_directory, species_name=species_name.replace(' ', '_'))
    inputs = {'read1': '{}'.format(glob.glob('{}/*_1.fq*'.format(hic_directory))[0]),
              'read2': '{}'.format(glob.glob('{}/*_2.fq*'.format(hic_directory))[0])}
    if os.path.splitext(inputs['read1'])[1] == '.gz':
        outputs = {'symlink1': '{top_dir}/fastq/{abbr}_R1.fastq.gz'.format(top_dir=top_dir, abbr=species_abbreviation(species_name)),
                   'symlink2': '{top_dir}/fastq/{abbr}_R2.fastq.gz'.format(top_dir=top_dir, abbr=species_abbreviation(species_name))}
    else:
        outputs = {'symlink1': '{top_dir}/fastq/{abbr}_R1.fastq'.format(top_dir=top_dir, abbr=species_abbreviation(species_name)),
                   'symlink2': '{top_dir}/fastq/{abbr}_R2.fastq'.format(top_dir=top_dir, abbr=species_abbreviation(species_name))}
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '00:05:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    mkdir -p {top_dir}
    mkdir {top_dir}/fastq # fastqdir
    # mkdir {top_dir}/splits # splitdir
    # mkdir {top_dir}/done_splits # donesplitsdir
    mkdir {top_dir}/aligned # outputdir
    mkdir {top_dir}/HIC_tmp # tmpdir

    ln -s {read1} {symlink1}
    ln -s {read2} {symlink2}
    # cp {symlink1} {top_dir}/splits
    # cp {symlink2} {top_dir}/splits
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(top_dir=top_dir)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def index_reference_genome_loop(reference_genome):
    """Template for indexing all reference genomes in path with 'bwa index' and 'samtools faidx'."""
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

def get_chrom_sizes(fasta_file: str, output_file: str):
    """Template for creating chrom.sizes file from FASTA file.\n
    chrom.sizes contains two columns, one with 'chromosome names' and one with the corresponding length."""
    inputs = {'fasta': fasta_file}
    outputs = {'chrom_sizes': output_file}
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '00:05:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    python {script} {fasta_file} {chrom_sizes}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(script='/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/04_genome_assembly/workflow_source/chrom_sizes.py', fasta_file=inputs['fasta'], chrom_sizes=outputs['chrom_sizes'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def bwa_align(read1: str, read2: str, reference_genome: str):
    """Template for aligning Hi-C reads to draft genome."""
    inputs = {'read1': read1,
              'read2': read2,
              'reference_genome': reference_genome}
    outputs = {'sam': '{read_name}.{read_ext}.sam'.format(read_name=read1.split(sep='_R1', maxsplit=1)[0], read_ext=read1.split(sep='_R1', maxsplit=1)[1])}
    options = {
        'cores': 36,
        'memory': '80g',
        'walltime': '12:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    bwa mem \
        -t {cores} \
        -S \
        -P \
        -5 \
        -M \
        {reference_genome} \
        {read1} \
        {read2} \
        > {sam}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], reference_genome=inputs['reference_genome'], read1=inputs['read1'], read2=inputs['read2'], sam=outputs['sam'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def file_sorting(sam_file: str, output_directory: str, temp_directory: str):
    """Template for running series of Juicer scripts dealing with chimeric reads, fragmentation and sorting."""
    blacklist_script='/home/jepe/software/juicer/scripts/chimeric_blacklist.awk'
    frag_script='/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/04_genome_assembly/workflow_source/frag.awk'
    inputs = {'sam': sam_file,
              'blacklist': blacklist_script,
              'frag': frag_script}
    outputs = {'normal': '{}_norm.txt'.format(os.path.splitext(inputs['sam'])[0]),
               'abnormal': '{}_abnorm.sam'.format(os.path.splitext(inputs['sam'])[0]),
               'unmapped': '{}_unmapped.sam'.format(os.path.splitext(inputs['sam'])[0]),
               'res': '{}_norm.txt.res.txt'.format(os.path.splitext(inputs['sam'])[0]),
               'frag': '{}.frag.txt'.format(os.path.splitext(inputs['sam'])[0]),
               'sort': '{output_dir}/{name}.sort.txt'.format(output_dir=output_directory, name=os.path.splitext(os.path.basename(inputs['sam']))[0])}
    options = {
        'cores': 12,
        'memory': '40g',
        'walltime': '24:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    awk \
        -v fname1={norm_file} \
        -v fname2={abnorm_file} \
        -v fname3={unmapped_file} \
        -f {blacklist_script} \
        {sam_file}
    
    awk \
        -f {frag_script} \
        {norm_file} \
        > {frag_file}

    sort \
        --parallel {cores} \
        -S 35G \
        -T {temp_dir} \
        -k2,2d \
        -k6,6d \
        -k4,4n \
        -k8,8n \
        -k1,1n \
        -k5,5n \
        -k3,3n \
        {frag_file} \
        > {sort_file}
        
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], blacklist_script=blacklist_script, norm_file=outputs['normal'], abnorm_file=outputs['abnormal'], unmapped_file=outputs['unmapped'], frag_script=frag_script, frag_file=outputs['frag'], temp_dir=temp_directory, sort_file=outputs['sort'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def de_duplication():
    """Template to run dedup script to remove duplicates from sorted file."""
    inputs = {}
    outputs = {}
    options = {
        'cores': 1,
        'memory': '12g',
        'walltime': '12:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    awk \
        -v 
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format()
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)











def hic_alignment(hic_directory: str, reference_genome: str, species_name: str, work_dir: str):
    """Template to use juicer for HiC alignment.\n
    [hic_directory] will be used to """
    juicer_dir = '/home/jepe/software/juicer'
    top_dir = '{work_dir}/04_genome_assembly/{species_name}/HiC'.format(work_dir=work_dir, species_name=species_name.replace(' ', '_'))
    inputs = {'read1': '{}'.format(glob.glob('{}/*_1.fq*'.format(hic_directory))[0]),
              'read2': '{}'.format(glob.glob('{}/*_2.fq*'.format(hic_directory))[0]),
              'reference': reference_genome}
    outputs = {'chrom_sizes': '{top_dir}/{species_name}.chrom.sizes'.format(top_dir=top_dir, species_name=species_abbreviation(species_name)),}
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00-00:10:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    if [ -e {top_dir} ]; then
        rm -rf {top_dir}
    fi

    mkdir -p {top_dir}/fastq
    ln -s {read1} '{top_dir}/fastq/{abbr}_R1.fastq'
    ln -s {read2} '{top_dir}/fastq/{abbr}_R2.fastq'

    python {chrom_sizes_script} {reference_genome} {chrom_sizes}

    bash {juicer_dir}/scripts/juicer.sh \
        -A EcoGenetics \
        -d {top_dir} \
        -D {juicer_dir} \
        -p {chrom_sizes} \
        -s none \
        -z {reference_genome} \
        -q short \
        -Q 12:00:00 \
        -l normal \
        -L 24:00:00 \
        -t 36
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(top_dir=top_dir, chrom_sizes=outputs['chrom_sizes'], reference_genome=inputs['reference'], juicer_dir=juicer_dir, read1=inputs['read1'], read2=inputs['read2'], abbr=species_abbreviation(species_name), chrom_sizes_script=CHROM_SIZES)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def assemble_3ddna(fasta_file: str, merged_no_dup: str, rounds: int = 2):
    """Template using 3d-dna to assemble draft assembly using the output of Juicer."""
    inputs = {'fasta': fasta_file,
              'mnd': merged_no_dup}
    outputs = {'logfile': '{}/3ddna.log'.format(os.path.dirname(fasta_file))}
    options = {
        'cores': 16,
        'memory': '138g',
        'walltime': '3-00:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    if [ ! -e {work_dir}/temp ]; then
        mkdir {work_dir}/temp
    fi
    export _JAVA_OPTIONS=-Djava.io.tmpdir={work_dir}/temp

    3d-dna \
        --rounds {r} \
        --editor-coarse-stringency 20 \
        --editor-repeat-covereage 30 \
        --splitter-input-size 500000 \
        --splitter-coarse-resolution 500000 \
        --splitter-coarse-stringency 20 \
        {path_to_input_fasta} \
        {path_to_input_mnd} \
        > {log}.prog
    
    mv {log}.prog {log}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(r=rounds, path_to_input_fasta=inputs['fasta'], path_to_input_mnd=inputs['mnd'], log=outputs['logfile'], work_dir=os.path.dirname(inputs['fasta']))
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def post_review_export(review_file: str, fasta_file: str, merged_no_dup: str, chrom_num: int):
    """Template for running the post review pipeline of 3d-dna, which finalizes assemblies, after review in JuiceBox Assembly Tools."""
    inputs = {'review': review_file,
              'fasta': fasta_file,
              'mnd': merged_no_dup}
    outputs = {'logfile': '{}/post_review.log'.format(os.path.dirname(review_file))}
    options = {
        'cores': 6,
        'memory': '40g',
        'walltime': '03-00:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    if [ ! -e {work_dir}/temp ]; then
        mkdir {work_dir}/temp
    fi

    export _JAVA_OPTIONS=-Djava.io.tmpdir={work_dir}/temp

    bash $(dirname $(dirname $(which 3d-dna)))/share/3d-dna/run-asm-pipeline-post-review.sh \
        --chromosome-map {c} \
        --sort-output \
        --stage finalize \
        --review {review} \
        {path_to_input_fasta} \
        {path_to_input_mnd} \
        > {log}.prog
    
    mv {log}.prog {log}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(c=chrom_num, review=inputs['review'], path_to_input_fasta=inputs['fasta'], path_to_input_mnd=inputs['mnd'], log=outputs['logfile'], work_dir=os.path.dirname(inputs['review']))
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def hifiadapterfilt(pacbio_hifi_reads: str):
    """Template for running HiFiAdapterFilt to remove remaining adapters specifically from PacBio HiFi reads."""
    hifiadapterfilt_directory = '/home/jepe/software/HiFiAdapterFilt'
    output_dir = '{}/HiFiAdapterFilt'.format(os.path.dirname(pacbio_hifi_reads))
    os.makedirs(output_dir, exist_ok=True)
    prefix, ext = os.path.splitext(pacbio_hifi_reads)
    if prefix.endswith(('fastq', 'fq')):
        prefix, ext2 = os.path.splitext(prefix)
        ext = ext2 + ext
    inputs = {'pbhifi': pacbio_hifi_reads}
    outputs = {'filt': '{output_dir}/{prefix}.filt.fastq.gz'.format(output_dir=output_dir, prefix=os.path.basename(prefix)),
               'cont': '{output_dir}/{prefix}.contaminant.blastout'.format(output_dir=output_dir, prefix=os.path.basename(prefix)),
               'block': '{output_dir}/{prefix}.blocklist'.format(output_dir=output_dir, prefix=os.path.basename(prefix)),
               'stats': '{output_dir}/{prefix}.stats'.format(output_dir=output_dir, prefix=os.path.basename(prefix))}
    options = {
        'cores': 10,
        'memory': '100g',
        'walltime': '120'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    export PATH=$PATH:{hifiadapterfilt_dir}
    export PATH=$PATH:{hifiadapterfilt_dir}/DB

    cp {pacbiohifi} "$(dirname({pacbiohifi}))"/prog.{ext}

    bash {hifiadapterfilt_dir}/pbadapterfilt.sh \
        -t {cores} \
        -p prog \
        -o {output_dir}
        
    mv {output_dir}/prog.filt.fastq.gz {filt}
    mv {output_dir}/prog.contaminant.blastout {cont}
    mv {output_dir}/prog.blocklist {block}
    mv {output_dir}/prog.stats {stats}
    rm "$(dirname({pacbiohifi}))"/prog.{ext}
        
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(pacbiohifi=inputs['pbhifi'], hifiadapterfilt_dir=hifiadapterfilt_directory, ext=ext, cores=options['cores'], prefix=prefix, output_dir=output_dir, filt=outputs['filt'], cont=outputs['cont'], block=outputs['block'], stats=outputs['stats'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def kmer_count():
    """Template for counting k'mers in genome sequence data"""
    inputs = {}
    outputs = {}
    options = {
        'cores': 32,
        'memory': '480g',
        'walltime': '240'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate genome_assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    jellyfish count \
        -t {cores} \
        -m {kmer_length} {canonical} \
        -s 8G \
        -o {count_file} \
        <(zcat {sequence_file})
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format()
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)