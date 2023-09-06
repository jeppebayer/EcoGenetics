from gwf import AnonymousTarget
import os

def species_abbreviation(species_name: str) -> str:
    """Function for creating species abbreviation from species name."""
    genus, species = species_name.split(maxsplit=1)
    genus = genus[0].upper() + genus[1:3]
    species = species[0].upper() + species[1:3]
    return genus + species

def merging_sequence_runs(sample_name: str, R1: list, R2: list, output_dir: str = "."):
    """Template for merging several paired-end read files from different runs for the same sample."""
    inputs = {'read1': R1,
              'read2': R2}
    outputs = {'read1': '{output_dir}/{sample_name}.R1.concat.fq.gz'.format(output_dir=output_dir, sample_name=sample_name),
               'read2': '{output_dir}/{sample_name}.R2.concat.fq.gz'.format(output_dir=output_dir, sample_name=sample_name)}
    options = {
        'cores': 1,
        'memory': '12g',
        'walltime': '02:00:00'
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

    cat {filelist_R1} > {concat_R1}
    cat {filelist_R2} > {concat_R2}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(filelist_R1=' '.join(inputs['read1']), filelist_R2=' '.join(inputs['read2']), concat_R1=outputs['read1'], concat_R2=outputs['read2'], output_dir=output_dir)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def AdapterRemoval_paired_end(read1: str, read2: str, sample_name: str, output_dir: str = "."):
    """Template for running AdapterRemoval to do initial trimming of raw paired-end sequence files."""
    inputs = {'reads': [read1,
                        read2]}
    outputs = {'pairs': ['{output_dir}/{sample_name}.pair1.truncated'.format(output_dir=output_dir, sample_name=sample_name),
                         '{output_dir}/{sample_name}.pair2.truncated'.format(output_dir=output_dir, sample_name=sample_name)],
                'single': '{output_dir}/{sample_name}.all_collapsed'.format(output_dir=output_dir, sample_name=sample_name),
                'discard': ['{output_dir}/{sample_name}.collapsed'.format(output_dir=output_dir, sample_name=sample_name),
                            '{output_dir}/{sample_name}.collapsed.truncated'.format(output_dir=output_dir, sample_name=sample_name),
                            '{output_dir}/{sample_name}.discarded'.format(output_dir=output_dir, sample_name=sample_name),
                            '{output_dir}/{sample_name}.settings'.format(output_dir=output_dir, sample_name=sample_name),
                            '{output_dir}/{sample_name}.singleton.truncated'.format(output_dir=output_dir, sample_name=sample_name)]}
    options = {
        'cores': 8,
        'memory': '64g',
        'walltime': '03:00:00'
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

    AdapterRemoval \
        --threads {cores} \
        --file1 {read1} \
        --file2 {read2} \
        --adapter1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
        --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
        --minquality 25 \
        --minlength 20 \
        --basename {output_dir}/{sample_name} \
        --trimns \
        --trimqualities \
        --collapse

    cat \
        {output_dir}/{sample_name}.collapsed \
        {output_dir}/{sample_name}.collapsed.truncated \
        > {all_collapsed}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], read1=inputs['reads'][0], read2=inputs['reads'][1], sample_name=sample_name, output_dir=output_dir, all_collapsed=outputs['single'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def alignment_contemp(sample_name: str, reference_genome: str, read1: str, read2: str = None, output_dir: str = "."):
    """Template for aligning contemporary reads"""
    if read2:
        read_type='paired_end'
        inputs = {'reads': [read1,
                            read2]}    
    else:
        read_type='single_end'
        inputs = {'reads': [read1]}
    outputs = {'alignment': '{output_dir}/{sample_name}.{read_type}.aligned.bam'.format(output_dir=output_dir, sample_name=sample_name, read_type=read_type)}
    options = {
        'cores': 8,
        'memory': '64g',
        'walltime': '03:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate data_prep
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_dir}/temp ] || mkdir -p {output_dir}/temp

    bwa mem \
        -t {cores} \
        -R "RG\tID:{sample_name}\tSM:{sample_name} \
        {reference_genome_noext} \
        {reads} \
    | samtools sort \
        -@ (({cores} - 1)) \
        -n \
        -O BAM \
        -T {output_dir}/temp \
        -o {alignment} \
        -
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], sample_name=sample_name, reference_genome_noext=os.path.splitext(reference_genome)[0], reads=' '.join(inputs['reads']), alignment=outputs['alignment'], output_dir=output_dir)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def merge_alignments(alignment1: str, alignment2: str, sample_name: str, output_dir: str = "."):
    """Template for merging alignments"""
    inputs = {'alignmnents': [alignment1, alignment2]}
    outputs = {'merged': '{output_dir}/{sample_name}.aligned.merged.bam'.format(output_dir=output_dir, sample_name=sample_name)}
    options = {
        'cores': 8,
        'memory': '64g',
        'walltime': '03:00:00'
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

    samtools merge \
        -@ (({cores} - 1)) \
        -c \
        -p \
        -n \
        -o {merged_alignment} \
        {alignment1} \
        {alignment2}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], alignment1=inputs['alignmnents'][0], alignment2=inputs['alignmnents'][1], merged_alignment=outputs['merged'], output_dir=output_dir)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def markduplicates(alignment: str, output_dir: str = "."):
    """Template for adding mate-tag and marking duplicates in BAM alignment"""
    inputs = {'alignment': alignment}
    outputs = {'markdup': '{output_dir}/{alignment_name}.markdup.bam'.format(output_dir=output_dir, alignment_name=os.path.basename(os.path.splitext(alignment)[0])),
               'stats': '{output_dir}/stats/{alignment_name}.markdupstats'.format(output_dir=output_dir, alignment_name=os.path.basename(os.path.splitext(alignment)[0])),
               'markdup_index': '{output_dir}/{alignment_name}.markdup.bam.bai'.format(output_dir=output_dir, alignment_name=os.path.basename(os.path.splitext(alignment)[0]))}
    options = {
        'cores': 8,
        'memory': '64g',
        'walltime': '02:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate data_prep
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_dir}/temp ] || mkdir -p {output_dir}/temp
    [ -d {output_dir}/stat ] || mkdir -p {output_dir}/stat

    samtools fixmate \
        -@ (({cores} - 1)) \
        -m \
        -O BAM \
        {alignment} \
        - \
    | samtools sort \
        -@ (({cores} - 1)) \
        -O BAM \
        -T {output_dir}/temp \
        - \
    | samtools markdup \
        -@ (({cores} - 1)) \
        -s \
        -f {stat_file} \
        -T {output_dir}/temp \
        - \
        {markdup_out}
    
    samtools index \
        -@ (({cores} - 1)) \
        -b \
        {markdup_index}    
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], alignment=inputs['alignment'], markdup_out=outputs['markdup'], stat_file=outputs['stats'], markdup_index=outputs['markdup_index'], output_dir=output_dir)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def extract_unmapped_reads(alignment: str, output_dir: str):
    """Template for extracting unmapped reads from alingment in BAM format."""
    inputs = {'alignmnet': alignment}
    outputs = {'unmapped_reads': '{output_dir}/{alignment_name}.unmapped.bam'.format(output_dir=output_dir, alignment_name=os.path.basename(os.path.splitext(alignment)[0])),
               'unmapped_index': '{output_dir}/{alignment_name}.unmapped.bam.bai'.format(output_dir=output_dir, alignment_name=os.path.basename(os.path.splitext(alignment)[0]))}
    options = {
        'cores': 8,
        'memory': '64g',
        'walltime': '02:00:00'
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

    samtools view \
        -@ (({cores} - 1)) \
        -b \
        -f 4 \
        -o {unmapped_reads} \
        {alignment}
    
    samtools index \
        -@ (({cores} - 1)) \
        -b \
        {unmapped_reads} \
        > {unmapped_index}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], alignment=inputs['alignmnet'], unmapped_reads=outputs['unmapped_reads'], unmapped_index=outputs['unmapped_index'], output_dir=output_dir)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def samtools_stats(alignment: str, output_dir: str = "."):
    """Template for running several samtools statistics on BAI indexed alignmnets in BAM format."""
    inputs = {'alignment': alignment}
    outputs = {'stats': ['{output_dir}/stats/{alignment_name}.idxstats'.format(output_dir=output_dir, alignment_name=os.path.basename(os.path.splitext(alignment)[0])),
                         '{output_dir}/stats/{alignment_name}.flagstat'.format(output_dir=output_dir, alignment_name=os.path.basename(os.path.splitext(alignment)[0])),
                         '{output_dir}/stats/{alignment_name}.coverage'.format(output_dir=output_dir, alignment_name=os.path.basename(os.path.splitext(alignment)[0])),
                         '{output_dir}/stats/{alignment_name}.genstats'.format(output_dir=output_dir, alignment_name=os.path.basename(os.path.splitext(alignment)[0]))]}
    options = {
        'cores': 8,
        'memory': '64g',
        'walltime': '04:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate data_prep
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_dir}/stats ] || mkdir -p {output_dir}/stats

    samtools idxstats \
        -@ (({cores} - 1)) \
        {alignment} \
        {idx_stats}
    
    samtools flagstat \
        -@ (({cores} - 1)) \
        {alignment} \
        {flag_stats}
    
    samtools coverage \
        -o {coverage_file} \
        {alignment}
    
    samtools stats \
        -@ (({cores} - 1)) \
        -c 1,1000,1 \
        {alignment} \
        {gen_stats}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], alignment=inputs['alignment'], idx_stats=outputs['stats'][0], flag_stats=outputs['stats'][1], coverage_file=outputs['stats'][2], gen_stats=outputs['stats'][3], output_dir=output_dir)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def filter_alignment(alignment: str, output_dir: str = "."):
    """Template for filtering alignment in BAM format."""
    inputs = {'alignment': alignment}
    outputs = {'filtered_alignment': '{output_dir}/{alignment_basename}.filtered.bam'.format(output_dir=output_dir, alignment_basename=os.path.basename(alignment).split(sep='.')[0]),
               'filtered_index': '{output_dir}/{alignment_basename}.filtered.bam.bai'.format(output_dir=output_dir, alignment_basename=os.path.basename(alignment).split(sep='.')[0])}
    options = {
        'cores': 8,
        'memory': '64g',
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
    
    [ -d {output_dir} ] || mkdir -p {output_dir}

    samtools view \
        -@ (({cores} - 1)) \
        -b \
        -F 3844 \
        -q 20 \
        -o {filtered_alignment} \
        {alignment}

    samtools index \
        -@ (({cores} - 1)) \
        -b \
        {filtered_alignment} \
        {filtered_index}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], alignment=inputs['alignment'], filtered_alignment=outputs['filtered_alignment'], filtered_index=outputs['filtered_index'], output_dir=output_dir)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def qualimap(alignment: str, output_dir: str = "."):
    """Tempalte for running qualimap on alignment in BAM format."""
    os.makedirs('{}/qualimap'.format(output_dir), exist_ok=True)
    inputs = {'alignment': alignment}
    outputs = {'qualimap_file': '{output_dir}/qualimap/{alignment_basename}.qualimap.pdf'.format(output_dir=output_dir, alignment_basename=os.path.basename(alignment).split(sep='.')[0])}
    options = {
        'cores': 8,
        'memory': '160g',
        'walltime': '03:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate data_prep
    fi
    
    export _JAVA_OPTIONS="-Xms160G -Xmx160G"

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_dir}/qualimap ] || mkdir -p {output_dir}/qualimap

    qualimap bamqc \
    -bam {alignment} \
    -outdir {output_dir}/qualimap \
    -outfile {qualimap_file} \
    -outformat PDF \
    --java-mem-size=100G
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(alignment=inputs['alignment'], output_dir=output_dir, qualimap_file=outputs['qualimap_file'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)