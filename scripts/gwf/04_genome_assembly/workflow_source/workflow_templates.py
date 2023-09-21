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

# def directory_setup(working_directory: str, species_name: str, hic_directory: str):
#     """Template for setting up folder structure Juicer jobs."""
#     top_dir = '{work_dir}/04_genome_assembly/{species_name}/HiC'.format(work_dir=working_directory, species_name=species_name.replace(' ', '_'))
#     inputs = {'read1': '{}'.format(glob.glob('{}/*_1.fq*'.format(hic_directory))[0]),
#               'read2': '{}'.format(glob.glob('{}/*_2.fq*'.format(hic_directory))[0])}
#     if os.path.splitext(inputs['read1'])[1] == '.gz':
#         outputs = {'symlink1': '{top_dir}/fastq/{abbr}_R1.fastq.gz'.format(top_dir=top_dir, abbr=species_abbreviation(species_name)),
#                    'symlink2': '{top_dir}/fastq/{abbr}_R2.fastq.gz'.format(top_dir=top_dir, abbr=species_abbreviation(species_name))}
#     else:
#         outputs = {'symlink1': '{top_dir}/fastq/{abbr}_R1.fastq'.format(top_dir=top_dir, abbr=species_abbreviation(species_name)),
#                    'symlink2': '{top_dir}/fastq/{abbr}_R2.fastq'.format(top_dir=top_dir, abbr=species_abbreviation(species_name))}
#     options = {
#         'cores': 1,
#         'memory': '8g',
#         'walltime': '00:05:00'
#     }
#     spec = """
#     # Sources environment
#     if [ "$USER" == "jepe" ]; then
#         source /home/"$USER"/.bashrc
#         source activate genome_assembly
#     fi
    
#     echo "START: $(date)"
#     echo "JobID: $SLURM_JOBID"
    
#     mkdir -p {top_dir}
#     mkdir {top_dir}/fastq # fastqdir
#     # mkdir {top_dir}/splits # splitdir
#     # mkdir {top_dir}/done_splits # donesplitsdir
#     mkdir {top_dir}/aligned # outputdir
#     mkdir {top_dir}/HIC_tmp # tmpdir

#     ln -s {read1} {symlink1}
#     ln -s {read2} {symlink2}
#     # cp {symlink1} {top_dir}/splits
#     # cp {symlink2} {top_dir}/splits
    
#     echo "END: $(date)"
#     echo "$(jobinfo "$SLURM_JOBID")"
#     """.format(top_dir=top_dir)
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def index_reference_genome_loop(reference_genome):
#     """Template for indexing all reference genomes in path with 'bwa index' and 'samtools faidx'."""
#     inputs = {'path': reference_genome}
#     outputs = {'path': ['{}.amb'.format(reference_genome),
#                         '{}.ann'.format(reference_genome),
#                         '{}.pac'.format(reference_genome),
#                         '{}.bwt'.format(reference_genome),
#                         '{}.sa'.format(reference_genome),
#                         '{}.fai'.format(reference_genome)]}
#     options = {
#         'cores': 1,
#         'memory': '16g',
#         'walltime': '02:00:00'
#     }
#     spec = """
#     if [ "$USER" == "jepe" ]; then
#         source /home/"$USER"/.bashrc
#         source activate data_prep
#     fi

#     echo "START: $(date)"
#     echo "JobID: $SLURM_JOBID"

#     bwa index \
#         -p {reference_genome} \
#         {reference_genome}
    
#     samtools faidx \
#         -o {reference_genome}.fai \
#         {reference_genome}
    
#     echo "END: $(date)"
#     echo "$(jobinfo "$SLURM_JOBID")"
#     """.format(reference_genome=reference_genome)
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def get_chrom_sizes(fasta_file: str, output_file: str):
#     """Template for creating chrom.sizes file from FASTA file.\n
#     chrom.sizes contains two columns, one with 'chromosome names' and one with the corresponding length."""
#     inputs = {'fasta': fasta_file}
#     outputs = {'chrom_sizes': output_file}
#     options = {
#         'cores': 1,
#         'memory': '8g',
#         'walltime': '00:05:00'
#     }
#     spec = """
#     # Sources environment
#     if [ "$USER" == "jepe" ]; then
#         source /home/"$USER"/.bashrc
#         source activate genome_assembly
#     fi
    
#     echo "START: $(date)"
#     echo "JobID: $SLURM_JOBID"
    
#     python {script} {fasta_file} {chrom_sizes}
    
#     echo "END: $(date)"
#     echo "$(jobinfo "$SLURM_JOBID")"
#     """.format(script='/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/04_genome_assembly/workflow_source/chrom_sizes.py', fasta_file=inputs['fasta'], chrom_sizes=outputs['chrom_sizes'])
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def bwa_align(read1: str, read2: str, reference_genome: str):
#     """Template for aligning Hi-C reads to draft genome."""
#     inputs = {'read1': read1,
#               'read2': read2,
#               'reference_genome': reference_genome}
#     outputs = {'sam': '{read_name}.{read_ext}.sam'.format(read_name=read1.split(sep='_R1', maxsplit=1)[0], read_ext=read1.split(sep='_R1', maxsplit=1)[1])}
#     options = {
#         'cores': 36,
#         'memory': '80g',
#         'walltime': '12:00:00'
#     }
#     spec = """
#     # Sources environment
#     if [ "$USER" == "jepe" ]; then
#         source /home/"$USER"/.bashrc
#         source activate genome_assembly
#     fi
    
#     echo "START: $(date)"
#     echo "JobID: $SLURM_JOBID"
    
#     bwa mem \
#         -t {cores} \
#         -S \
#         -P \
#         -5 \
#         -M \
#         {reference_genome} \
#         {read1} \
#         {read2} \
#         > {sam}
    
#     echo "END: $(date)"
#     echo "$(jobinfo "$SLURM_JOBID")"
#     """.format(cores=options['cores'], reference_genome=inputs['reference_genome'], read1=inputs['read1'], read2=inputs['read2'], sam=outputs['sam'])
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def file_sorting(sam_file: str, output_directory: str, temp_directory: str):
#     """Template for running series of Juicer scripts dealing with chimeric reads, fragmentation and sorting."""
#     blacklist_script='/home/jepe/software/juicer/scripts/chimeric_blacklist.awk'
#     frag_script='/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/04_genome_assembly/workflow_source/frag.awk'
#     inputs = {'sam': sam_file,
#               'blacklist': blacklist_script,
#               'frag': frag_script}
#     outputs = {'normal': '{}_norm.txt'.format(os.path.splitext(inputs['sam'])[0]),
#                'abnormal': '{}_abnorm.sam'.format(os.path.splitext(inputs['sam'])[0]),
#                'unmapped': '{}_unmapped.sam'.format(os.path.splitext(inputs['sam'])[0]),
#                'res': '{}_norm.txt.res.txt'.format(os.path.splitext(inputs['sam'])[0]),
#                'frag': '{}.frag.txt'.format(os.path.splitext(inputs['sam'])[0]),
#                'sort': '{output_dir}/{name}.sort.txt'.format(output_dir=output_directory, name=os.path.splitext(os.path.basename(inputs['sam']))[0])}
#     options = {
#         'cores': 12,
#         'memory': '40g',
#         'walltime': '24:00:00'
#     }
#     spec = """
#     # Sources environment
#     if [ "$USER" == "jepe" ]; then
#         source /home/"$USER"/.bashrc
#         source activate genome_assembly
#     fi
    
#     echo "START: $(date)"
#     echo "JobID: $SLURM_JOBID"
    
#     awk \
#         -v fname1={norm_file} \
#         -v fname2={abnorm_file} \
#         -v fname3={unmapped_file} \
#         -f {blacklist_script} \
#         {sam_file}
    
#     awk \
#         -f {frag_script} \
#         {norm_file} \
#         > {frag_file}

#     sort \
#         --parallel {cores} \
#         -S 35G \
#         -T {temp_dir} \
#         -k2,2d \
#         -k6,6d \
#         -k4,4n \
#         -k8,8n \
#         -k1,1n \
#         -k5,5n \
#         -k3,3n \
#         {frag_file} \
#         > {sort_file}
        
#     echo "END: $(date)"
#     echo "$(jobinfo "$SLURM_JOBID")"
#     """.format(cores=options['cores'], blacklist_script=blacklist_script, norm_file=outputs['normal'], abnorm_file=outputs['abnormal'], unmapped_file=outputs['unmapped'], frag_script=frag_script, frag_file=outputs['frag'], temp_dir=temp_directory, sort_file=outputs['sort'])
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# def de_duplication():
#     """Template to run dedup script to remove duplicates from sorted file."""
#     inputs = {}
#     outputs = {}
#     options = {
#         'cores': 1,
#         'memory': '12g',
#         'walltime': '12:00:00'
#     }
#     spec = """
#     # Sources environment
#     if [ "$USER" == "jepe" ]; then
#         source /home/"$USER"/.bashrc
#         source activate genome_assembly
#     fi
    
#     echo "START: $(date)"
#     echo "JobID: $SLURM_JOBID"
    
#     awk \
#         -v 
    
#     echo "END: $(date)"
#     echo "$(jobinfo "$SLURM_JOBID")"
#     """.format()
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)











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

def get_chrom_sizes(fasta_file: str, output_directory: str = None):
    """
    Template: Create *chrom.sizes* from `FASTA` file.
    
    *chrom.sizes* contains two columns, one with chromosome names and one with the corresponding length.
    If no :param:`output_directory` is provided the directory of :param:`fasta_file` will be used.

    Template I/O::
    
        inputs = {'fasta': fasta_file}
        outputs = {'chrom_sizes': output_directory/chrom.sizes}
    
    :param str fasta_file:
        Genome file in `FASTA` format.
    :param str output_directory:
        Directory to place resulting *chrom.sizes*
    """
    if output_directory is None:
        output_directory = os.path.dirname(fasta_file)
    inputs = {'fasta': fasta_file}
    outputs = {'chrom_sizes': '{}/chrom.sizes'.format(output_directory)}
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '00:05:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate omni_c
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    python {script} {fasta_file} {chrom_sizes}.prog
    
    mv {chrom_sizes}.prog {chrom_sizes}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(script='/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/04_genome_assembly/workflow_source/chrom_sizes.py', fasta_file=inputs['fasta'], chrom_sizes=outputs['chrom_sizes'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def hic_align(read1: str, read2: str, draft_genome: str, species_name: str, output_directory: str = None):
    """
    Template: Align Hi-C data to draft genome using :script:`bwa mem`.

    Runs :script:`bwa mem` with the settings:

        `-S`, (Skip mate rescue).

        `-P`, (Skip pairing) In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.
        
        `-5`, For  split  alignment,  mark the segment with the smallest coordinate as the primary. This option may help some Hi-C pipelines. By default, BWA-MEM marks highest scoring segment as primary.
        
        `-T 0`, Don't output alignment with score lower than 0. This is set as all alignmnets are wanted at this stage.

    Template I/O::

        inputs = {'read1': read1,
                  'read2': read2,
                  'reference_genome': draft_genome}
        
        outputs = {'sam': output_directory/*.sam}

    :param str read1:
        Read pair 1 sequence file.
    :param str read2:
        Read pair 2 sequence file.
    :param str draft_genome:
        Assembled draft genome.
    :param str species_name:
        Name of species in format 'genus species' | 'genus_species'.
    :param str output_directory:
        Directory to place resulting *.sam* file
    """
    if output_directory is None:
        output_directory = os.path.dirname(read1)
    inputs = {'read1': read1,
              'read2': read2,
              'reference_genome': draft_genome}
    outputs = {'sam': '{output_dir}/{read_name}{read_ext}.sam'.format(output_dir=output_directory, read_name=os.path.basename(read1).split(sep='_R1', maxsplit=1)[0], read_ext=read1.split(sep='_R1', maxsplit=1)[1])}
    options = {
        'cores': 36,
        'memory': '80g',
        'walltime': '12:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate omni_c
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    bwa mem \
        -t {cores} \
        -R '@RG\tID:{species_abbr}.HiC\tSM:HiC_to_draft' \
        -S \
        -P \
        -5 \
        -T 0 \
        -o {sam}.prog \
        {reference_genome} \
        {read1} \
        {read2}
    
    mv {sam}.prog {sam}
        
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], species_abbr=species_abbreviation(species_name), reference_genome=inputs['reference_genome'], read1=inputs['read1'], read2=inputs['read2'], sam=outputs['sam'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def ligation_events(sam_file: str, chrom_sizes: str, species_name: str):
    """
    Template: Finds ligation junctions using :script:`pairtools parse`.

    Runs :script:`pairtools parse` with the settings:

        `--min-mapq 40`, The minimal MAPQ score to consider a read as uniquely mapped. Alignments with <40 will be marked as a multi-mapping alignment.

        `--walks-policy 5unique`, Walks is the term used to describe multiple ligations events, resulting three alignments (instead of two) for a read pair. Policy for reporting unrescuable walks (reads containing more than one alignment on one or both sides, that can not be explained by a single ligation between two mappable DNA fragments). `5unique` reports the 5'-most unique alignment on each side, if present.
        
        `--max-inter-align-gap 30`, Read segments that are not covered by any alignment and longer than the specified value are treated as "null" alignments. These null alignments convert otherwise linear alignments into walks, and affect how they get reported as a Hi-C pair.
    
    Template I/O::

        inputs = {'sam': sam_file,
                  'chrom_sizes': chrom_sizes}
        
        outputs = {'pairsam': *.sorted.pairsam}
    
    :param str sam_file:
        Hi-C alignmnet in `SAM` format.
    :param str chrom_sizes:
        File listing all chromosomes in draft genome and their respective lengths.
    :param str species_name:
        Name of species in format 'genus species' | 'genus_species'.
    """
    inputs = {'sam': sam_file,
              'chrom_sizes': chrom_sizes}
    outputs = {'pairsam': '{read_name}.sorted.pairsam'.format(read_name=os.path.splitext(inputs['sam'])[0])}
    options = {
        'cores': 12,
        'memory': '96g',
        'walltime': '12:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate omni_c
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    temp_dir="$(dirname {sam})"/temp
    if [ ! -e "$temp_dir" ]; then
        mkdir -m 775 "$temp_dir"
    fi

    pairtools parse \
        --assembly EG_{species_abbr} \
        --min-mapq 40 \
        --walks-policy 5unique \
        --max-inter-align-gap 30 \
        --nproc-in {cores} \
        --nproc-out {cores} \
        --chroms-path {chrom_sizes} \
        {sam} \
    | pairtools sort \
        --tmpdir "$temp_dir" \
        --nproc {cores} \
        > {sorted_pairsam}.prog
    
    mv {sorted_pairsam}.prog {sorted_pairsam}
        
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], chrom_sizes=inputs['chrom_sizes'], sam=inputs['sam'], species_abbr=species_abbreviation(species_name), sorted_pairsam=outputs['pairsam'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def de_duplicate(pairsam_file: str):
    """
    Template: Remove duplicates from alignment using :script:`pairtools dedup`.
    
    Runs `pairtools dedup` with the settings:

        `--mark-dups`, Duplicate pairs are marked as DD in "pair_type" and as a duplicate in the sam entries.

        `--output-stats FILE`, Output file for duplicate statistics.

    :param str pairsam_file:
        `PAIRSAM` file with `SAM` entries together with the Hi-C pair information.
    """
    inputs = {'pairsam': pairsam_file}
    outputs = {'dedup': '{file_name}.dedup.pairsam'.format(file_name=os.path.splitext(inputs['pairsam'])[0]),
               'stat': '{path}/dedup.stats'.format(path=os.path.dirname(inputs['pairsam']))}
    options = {
        'cores': 12,
        'memory': '96g',
        'walltime': '12:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate omni_c
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    pairtools dedup \
        --nproc-in {cores} \
        --nproc-out {cores} \
        --mark-dups \
        --output-stats {stat_file}.prog \
        --output {dedup_pairsam}.prog \
        {sorted_pairsam}
    
    mv {stat_file}.prog {stat_file}
    mv {dedup_pairsam}.prog {dedup_pairsam}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], stat_file=outputs['stat'], dedup_pairsam=outputs['dedup'], sorted_pairsam=inputs['pairsam'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def split_pairsam(pairsam_file: str):
    """
    Template: Splits `PAIRSAM` file into a `BAM` file and a *.pairs* file using `pairtools split`.
    
    Template I/O::

        inputs = {'pairsam': pairsam_file}

        outputs = {'pairs': *.pairs,
                   'bam': *.bam,
                   'index': *.bam.bai}
    
    :param str pairsam_file:
        **Sorted** `PAIRSAM` file with `SAM` entries together with the Hi-C pair information.
    """
    inputs = {'pairsam': pairsam_file}
    outputs = {'pairs': '{file_name}.pairs'.format(file_name=os.path.splitext(inputs['pairsam'])[0]),
               'bam': '{file_name}.bam'.format(file_name=os.path.splitext(inputs['pairsam'])[0]),
               'index': '{file_name}.bam.bai'.format(file_name=os.path.splitext(inputs['pairsam'])[0])}
    options = {
        'cores': 12,
        'memory': '96g',
        'walltime': '12:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate omni_c
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    temp_dir="$(dirname {pairsam})"/temp
    if [ ! -e "$temp_dir" ]; then
        mkdir -m 775 "$temp_dir"
    fi

    pairtools split \
        --nproc-in {cores} \
        --nproc-out {cores} \
        --output-pairs {pairs}.prog \
        --output-sam {bam}.prog \
        {pairsam}
    
    # samtools sort \
    #     -@ {cores} \
    #     -T "$temp_dir" \
    #     -o {bam}.prog
    
    # samtools index \
    #     -@ {cores} \
    #     -b \
    #     -o {index}.prog \
    #     {bam}
    
    mv {pairs}.prog {pairs}
    mv {bam}.prog {bam}
    mv {index}.prog {index}
        
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], pairsam=inputs['pairsam'], pairs=outputs['pairs'], bam=outputs['bam'], index=outputs['index'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def juicer_hic_matrix(pairs_file: str, chrom_sizes: str, juicer_script: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/04_genome_assembly/workflow_source/juicer_tools.2.20.00.jar'):
    """
    Template:Create Hi-C contact matrix from a *.pairs* file using :script:`juicer_tools pre`.
    
    Template I/O::

        inputs = {'pairs': pairs_file,
                  'chrom_sizes': chrom_sizes}

        outputs = {'hic': *.init_contact_map.hic}
    
    :param str pairs_file:
        *.pairs* file containing Hi-C alignment pairs.
    :param str chrom_sizes:
        File listing all chromosomes in draft genome and their respective lengths.
    :param str juicer_script:
        Path to *juicer_tools.jar*. By default should lead to *juicer_tools.jar* in workflow_source directory.
    """
    inputs = {'pairs': pairs_file,
              'chrom_sizes': chrom_sizes}
    outputs = {'hic': '{path}/{base_name}.init_contact_map.hic'.format(path=os.path.dirname(inputs['pairs']), base_name=os.path.basename(inputs['pairs']).split(sep='.')[0])}
    options = {
        'cores': 12,
        'memory': '100g',
        'walltime': '12:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate omni_c
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    export _JAVA_OPTIONS="-Xms{mem} -Xmx{mem}"

    java \
        -Xms {mem} \
        -Djava.awt.headless=true -jar {juicer_script} \
        -jar {juicer_script} \
            pre \
            --threads {cores} \
            {pairs} \
            {hic_matrix}.prog \
            {chrom_sizes}
    
    mv {hic_matrix}.prog {hic_matrix}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(mem=options['memory'], juicer_script=juicer_script, cores=options['cores'], pairs=inputs['pairs'], hic_matrix=outputs['hic'], chrom_sizes=inputs['chrom_sizes'])
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









# TODO Fix
# def hifiadapterfilt(pacbio_hifi_reads: str, HiFiAdapterFilt_path: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/04_genome_assembly/workflow_source/HiFiAdapterFilt'):
#     """Template for running HiFiAdapterFilt to remove remaining adapters specifically from PacBio HiFi reads."""
#     hifiadapterfilt_directory = '/home/jepe/software/HiFiAdapterFilt'
#     output_dir = '{}/HiFiAdapterFilt'.format(os.path.dirname(pacbio_hifi_reads))
#     inputs = {'pb_hifi': pacbio_hifi_reads}
#     if inputs['pb_hifi'].endswith('.gz'):
#         prefix = os.path.splitext(os.path.splitext(os.path.basename(inputs['pb_hifi']))[0])[0]
#         ext = 
#     else:
#         prefix = 
#     outputs = {'filt': '{output_dir}/{prefix}.filt.fastq.gz'.format(output_dir=output_dir, prefix=os.path.basename(prefix)),
#                'cont': '{output_dir}/{prefix}.contaminant.blastout'.format(output_dir=output_dir, prefix=os.path.basename(prefix)),
#                'block': '{output_dir}/{prefix}.blocklist'.format(output_dir=output_dir, prefix=os.path.basename(prefix)),
#                'stats': '{output_dir}/{prefix}.stats'.format(output_dir=output_dir, prefix=os.path.basename(prefix))}
#     options = {
#         'cores': 10,
#         'memory': '100g',
#         'walltime': '120'
#     }
#     spec = """
#     # Sources environment
#     if [ "$USER" == "jepe" ]; then
#         source /home/"$USER"/.bashrc
#         source activate genome_assembly
#     fi
    
#     echo "START: $(date)"
#     echo "JobID: $SLURM_JOBID"
    
#     export PATH=$PATH:{hifiadapterfilt_dir}
#     export PATH=$PATH:{hifiadapterfilt_dir}/DB

#     cp {pacbiohifi} "$(dirname({pacbiohifi}))"/prog.{ext}

#     bash {hifiadapterfilt_dir}/pbadapterfilt.sh \
#         -t {cores} \
#         -p prog \
#         -o {output_dir}
        
#     mv {output_dir}/prog.filt.fastq.gz {filt}
#     mv {output_dir}/prog.contaminant.blastout {cont}
#     mv {output_dir}/prog.blocklist {block}
#     mv {output_dir}/prog.stats {stats}
#     rm "$(dirname({pacbiohifi}))"/prog.{ext}
        
#     echo "END: $(date)"
#     echo "$(jobinfo "$SLURM_JOBID")"
#     """.format(pacbiohifi=inputs['pbhifi'], hifiadapterfilt_dir=hifiadapterfilt_directory, ext=ext, cores=options['cores'], prefix=prefix, output_dir=output_dir, filt=outputs['filt'], cont=outputs['cont'], block=outputs['block'], stats=outputs['stats'])
#     return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

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