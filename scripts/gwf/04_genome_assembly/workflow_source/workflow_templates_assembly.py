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
    outputs = {'fasta': f'{output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.bp.p_ctg.fasta',
               'primary': f'{output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.bp.p_ctg.gfa',
               'misc': [f'{output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.a_ctg.gfa',
                        f'{output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.r_utg.gfa',
                        f'{output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.p_ctg.gfa']}
    protect = [outputs['primary'], outputs['fasta']]
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
        > {output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.bp.p_ctg.prog.fasta

    mv {output_directory_path}/hifiasm/{species_abbreviation(species_name)}.asm.bp.p_ctg.prog.fasta {outputs['fasta']}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs,protect=protect, options=options, spec=spec)

def busco(assembly_file: str, output_directory_path: str, dataset: str = "/faststorage/project/EcoGenetics/BACKUP/database/busco/busco_downloads/lineages/arthropoda_odb10"):
    """
    Template: Runs BUSCO analysis on genome assembly.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'assembly': assembly_file}
    outputs = {'busco': [f'{output_directory_path}/busco_{os.path.splitext(os.path.basename(assembly_file))[0]}/run_{os.path.basename(dataset)}/full_table.tsv',
                         f'{output_directory_path}/busco_{os.path.splitext(os.path.basename(assembly_file))[0]}/run_{os.path.basename(dataset)}/missing_busco_list.tsv',
                         f'{output_directory_path}/busco_{os.path.splitext(os.path.basename(assembly_file))[0]}/run_{os.path.basename(dataset)}/short_summary.json',
                         f'{output_directory_path}/busco_{os.path.splitext(os.path.basename(assembly_file))[0]}/run_{os.path.basename(dataset)}/short_summary.txt',
                         f'{output_directory_path}/busco_{os.path.splitext(os.path.basename(assembly_file))[0]}/short_summary.specific.{os.path.basename(dataset)}.busco_{os.path.splitext(os.path.basename(assembly_file))[0]}.json',
                         f'{output_directory_path}/busco_{os.path.splitext(os.path.basename(assembly_file))[0]}/short_summary.specific.{os.path.basename(dataset)}.busco_{os.path.splitext(os.path.basename(assembly_file))[0]}.json']}
    protect = outputs['busco']
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
    
    busco \
        -f \
        -i {assembly_file} \
        -m genome \
        -o busco_{os.path.splitext(os.path.basename(assembly_file))[0]} \
        --out_path {output_directory_path} \
        -l {dataset}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

# Draft genome step 4 (purge_dup'ing preliminary assembly)
def minimap_raw(assembly_file: str, pacbio_hifi_reads: str, output_directory_path: str, species_name: str, round_number: int = 1):
    """
    Template: Aligns raw PacBio HiFi data to current assembly, producing a :format:`paf` file using minimap2.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'assembly': assembly_file,
              'pacbio': pacbio_hifi_reads}
    outputs = {'paf': f'{output_directory_path}/purge_dups/{round_number:02}/{species_abbreviation(species_name)}.paf.gz',
               'stats': [f'{output_directory_path}/purge_dups/{round_number:02}/PB.stat',
                         f'{output_directory_path}/purge_dups/{round_number:02}/PB.base.cov']}
    options = {
        'cores': 16,
        'memory': '160g',
        'walltime': '01:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory_path}/purge_dups/{round_number:02} ] || mkdir -p {output_directory_path}/purge_dups/{round_number:02}

    minimap2 \
        -x map-hifi \
        -t {options['cores']} \
        {assembly_file} \
        {pacbio_hifi_reads} \
    | gzip \
        -c \
        - \
        > {output_directory_path}/purge_dups/{round_number:02}/{species_abbreviation(species_name)}.prog.paf.gz 
    
    pbcstat \
        -O {output_directory_path}/purge_dups/{round_number:02}/ \
        {output_directory_path}/purge_dups/{round_number:02}/{species_abbreviation(species_name)}.prog.paf.gz

    mv {output_directory_path}/purge_dups/{round_number:02}/{species_abbreviation(species_name)}.prog.paf.gz {outputs['paf']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def split_assembly(assembly_file: str):
    """
    Template: Splits assembly in :format:`fasta` format into equal pieces.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'assembly': assembly_file}
    outputs = {'split': f'{os.path.splitext(assembly_file)[0]}.split.{os.path.splitext(assembly_file)[1]}'}
    options = {
        'cores': 2,
        'memory': '20',
        'walltime': '01:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    split_fa \
        {assembly_file} \
        > {os.path.splitext(assembly_file)[0]}.split.prog.{os.path.splitext(assembly_file)[1]}
    
    mv {os.path.splitext(assembly_file)[0]}.split.prog.{os.path.splitext(assembly_file)[1]} {outputs['split']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def minimap_self(split_assembly_file: str, species_name: str):
    """
    Template: Aligns a split assembly to itself using :script:`minimap2`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'assembly': split_assembly_file}
    outputs = {'paf': f'{os.path.dirname(split_assembly_file)}/{species_abbreviation(species_name)}.self.paf.gz'}
    options = {
        'cores': 16,
        'memory': '160g',
        'walltime': '01:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    minimap2 \
        -x asm5 \
        -t {options['cores']} \
        -DP \
        {split_assembly_file} \
        {split_assembly_file} \
    | gzip \
        -c \
        - \
        > {os.path.dirname(split_assembly_file)}/{species_abbreviation(species_name)}.self.prog.paf.gz
    
    mv {os.path.dirname(split_assembly_file)}/{species_abbreviation(species_name)}.self.prog.paf.gz {outputs['paf']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def purge_duplicates(PB_stat: str, PB_base_cov: str, self_alignment_paf: str, assembly_file: str, species_name: str):
    """
    Template: Purge haplotigs and overlaps using :script:`purge_dups`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'stats': [PB_stat,
                        PB_base_cov],
                'paf': self_alignment_paf,
                'assembly': assembly_file}
    outputs = {'cutoffs': f'{os.path.dirname(PB_stat)}/cutoffs',
               'dups': f'{os.path.dirname(self_alignment_paf)}/dups.bed',
               'purged': f'{os.path.dirname(self_alignment_paf)}/{species_abbreviation(species_name)}.purged.fa'}
    options = {
        'cores': 2,
        'memory': '20g',
        'walltime': '02:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    calcuts \
        {PB_stat} \
        > {os.path.dirname(PB_stat)}/cutoffs.prog

    mv {os.path.dirname(PB_stat)}/cutoffs.prog {outputs['cutoffs']}

    purge_dups \
        -2 \
        -T {outputs['cutoffs']} \
        -c {PB_base_cov} \
        {self_alignment_paf} \
        > {os.path.dirname(self_alignment_paf)}/dups.prog.bed
    
    mv {os.path.dirname(self_alignment_paf)}/dups.prog.bed {outputs['dups']}
    
    get_seqs \
        -p {os.path.dirname(self_alignment_paf)}/{species_abbreviation(species_name)}.prog \
        {outputs['dups']} \
        {assembly_file}

    mv {os.path.dirname(self_alignment_paf)}/{species_abbreviation(species_name)}.prog.purged.fa {outputs['purged']}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Hi-C
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
    outputs = {'path': [f'{reference_genome}.amb',
                        f'{reference_genome}.ann',
                        f'{reference_genome}.pac',
                        f'{reference_genome}.bwt',
                        f'{reference_genome}.sa',
                        f'{reference_genome}.fai']}
    protect = outputs['path']
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '02:00:00'
    }
    spec = f"""
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
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
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, protect=protect, spec=spec)

def hic_align(hic_read1: str, hic_read2: str, draft_genome: str, indices: list, output_directory_path: str, species_name: str):
    """
    Template: Align Hi-C data to draft genome using :script:`bwa mem`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'read1': hic_read1,
              'read2': hic_read2,
              'draft': draft_genome,
              'indices': indices}
    outputs = {'bam': f'{output_directory_path}/HiC/alignment/{species_abbreviation(species_name)}.HiC_to_draft.bam'}
    options = {
        'cores': 36,
        'memory': '288g',
        'walltime': '24:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory_path}/HiC/alignment/tmp ] || mkdir -p {output_directory_path}/HiC/alignment/tmp

    bwa mem \
        -t {options['cores']} \
        -R '@RG\\tID:{species_abbreviation(species_name)}.HiC\\tSM:HiC_to_draft' \
        -S \
        -P \
        -5 \
        -T 0 \
        {draft_genome} \
        {hic_read1} \
        {hic_read2} \
    | samtools sort \
        -@ {options['cores']} \
        -O BAM \
        -T {output_directory_path}/HiC/alignment/tmp \
        -o {output_directory_path}/HiC/alignment/{species_abbreviation(species_name)}.HiC_to_draft.prog.bam
    
    mv {output_directory_path}/HiC/alignment/{species_abbreviation(species_name)}.HiC_to_draft.prog.bam {outputs['bam']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mark_duplicates(bam_file: str, output_directory_path: str):
    """
    Template: Mark duplicates in :format:`BAM` alignment file using PICARD MarkDuplicates.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'bam': bam_file}
    outputs = {'markdup': f'{os.path.splitext(bam_file)[0]}.markdup.bam',
               'metrics': f'{os.path.splitext(bam_file)[0]}.markdup.txt'}
    protect = outputs['metrics']
    options = {
        'cores': 2,
        'memory': '100g',
        'walltime': '24:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory_path}/HiC/alignment/tmp ] || mkdir -p {output_directory_path}/HiC/alignment/tmp

    export _JAVA_OPTIONS="-Xmx100G"

    picard MarkDuplicates \
        --INPUT {bam_file} \
        --OUTPUT {os.path.splitext(bam_file)[0]}.markdup.prog.bam \
        --METRICS_FILE {outputs['metrics']} \
        --TMP_DIR {output_directory_path}/HiC/alignment/tmp
    
    mv {os.path.splitext(bam_file)[0]}.markdup.prog.bam {outputs['markdup']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

# The input BAM could either sorted by read names (samtools sort with -n option) or not.
# The behaviours of the program are slightly different, which might lead to slightly different scaffolding results.
# For a BAM input sorted by read names, with each mapped read pair, a Hi-C link is counted between the middle positions of the read alignments; while for a BAM input sorted by coordinates or unsorted, Hi-C links are counted between the start positions of the read alignments.
# Also, for a BAM input not sorted by read names, the mapping quality filtering is suppressed (-q option).
def hic_scaffolding(draft_genome: str, hic_to_draft_bam: str, output_directory_path: str, species_name: str):
    """
    Template: Genome scaffolding with Hi-C data using :script:`YaHS`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'draft': draft_genome,
              'hic': hic_to_draft_bam}
    outputs = {'final': [f'{output_directory_path}/HiC/YaHS/{species_abbreviation(species_name)}_scaffolds_final.agp',
                         f'{output_directory_path}/HiC/YaHS/{species_abbreviation(species_name)}_scaffolds_final.fa'],
                'bin': f'{output_directory_path}/HiC/YaHS/{species_abbreviation(species_name)}.bin'}
    protect = [outputs['final'][0], outputs['final'][1], outputs['bin']]
    options = {
        'cores': 30,
        'memory': '400g',
        'walltime': '24:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory_path}/HiC/YaHS ] || mkdir -p {output_directory_path}/HiC/YaHS

    yahs \
        -r 1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000 \
        -o {output_directory_path}/HiC/YaHS/{species_abbreviation(species_name)} \
        {draft_genome} \
        {hic_to_draft_bam}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def alignment_conversion_no_manual_curation(hic_bin: str, scaffolds_final_agp: str, draft_assembly_fai_index: str, output_directory_path: str, species_name: str):
    """
    Template: Convert Hi-C alignment file to format required by juicer_tools.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'bin': hic_bin,
              'scaffolds': scaffolds_final_agp,
              'index': draft_assembly_fai_index}
    outputs = {'sorted': f'{output_directory_path}/HiC/YaHS/no_curation/{species_abbreviation(species_name)}.alignments_sorted.txt'}
    protect = outputs['sorted']
    options = {
        'cores': 32,
        'memory': '256g',
        'walltime': '12:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory_path}/HiC/YaHS/tmp ] || mkdir -p {output_directory_path}/HiC/YaHS/tmp
    [ -d {output_directory_path}/HiC/YaHS/no_curation ] || mkdir -p {output_directory_path}/HiC/YaHS/no_curation

    juicer pre \
        {hic_bin} \
        {scaffolds_final_agp} \
        {draft_assembly_fai_index} \
    | sort \
        -k2,2d \
        -k6,6d \
        -T {output_directory_path}/HiC/YaHS/tmp \
        --parallel={options['cores']} \
        -S {options['memory']} \
    | awk \
        'NF' \
        > {output_directory_path}/HiC/YaHS/no_curation/{species_abbreviation(species_name)}.alignments_sorted.prog.txt

    mv {output_directory_path}/HiC/YaHS/no_curation/{species_abbreviation(species_name)}.alignments_sorted.prog.txt {outputs['sorted']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def alignment_conversion_manual_curation(hic_bin: str, scaffolds_final_agp: str, draft_assembly_fai_index: str, output_directory_path: str, species_name: str):
    """
    Template: Convert Hi-C alignment file to format required by juicer_tools.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'bin': hic_bin,
              'scaffolds': scaffolds_final_agp,
              'index': draft_assembly_fai_index}
    outputs = {'jbat': [f'{output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.txt',
                        f'{output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.liftover.agp',
                        f'{output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.assembly',
                        f'{output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.assembly.agp',
                        f'{output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.log']}
    protect = outputs['jbat']
    options = {
        'cores': 2,
        'memory': '30g',
        'walltime': '12:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory_path}/HiC/YaHS/tmp ] || mkdir -p {output_directory_path}/HiC/YaHS/tmp
    [ -d {output_directory_path}/HiC/YaHS/curation ] || mkdir -p {output_directory_path}/HiC/YaHS/curation

    juicer pre \
        -a \
        -o {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog \
        {hic_bin} \
        {scaffolds_final_agp} \
        {draft_assembly_fai_index} \
        > {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.log 2>&1
    
    mv {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.txt {outputs['jbat'][0]}
    mv {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.liftover.agp {outputs['jbat'][1]}
    mv {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.assembly {outputs['jbat'][2]}
    mv {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.assembly.agp {outputs['jbat'][3]}
    mv {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.log {outputs['jbat'][4]}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def chrom_sizes(scaffolds_final_fa: str, script: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/bin/get_length.py'):
    """
    Template: Create file containing two columns. In column 1 the name of each sequence, in column 2 the length of each sequence.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'scaffolds': scaffolds_final_fa}
    outputs = {'sizes': f'{scaffolds_final_fa}.chrom_sizes'}
    protect = outputs['sizes']
    options = {
        'cores': 2,
        'memory': '16g',
        'walltime': '06:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    {script} \
        {scaffolds_final_fa} \
        {scaffolds_final_fa}.prog.chrom_sizes
    
    mv {scaffolds_final_fa}.prog.chrom_sizes {outputs['sizes']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def contact_matrix_no_manual_curation(alignments_sorted: str, chrom_sizes: str, output_directory_path: str, species_name: str, juicer_tools: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/04_genome_assembly/workflow_source/juicer_tools.2.20.00.jar'):
    """
    Template: Generate Hi-C contact matrix using :script:`juicer pre`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'sorted': alignments_sorted,
              'sizes': chrom_sizes}
    outputs = {'hic': f'{output_directory_path}/HiC/YaHS/no_curation/{species_abbreviation(species_name)}.hic'}
    protect = outputs['hic']
    options = {
        'cores': 32,
        'memory': '256g',
        'walltime': '48:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory_path}/HiC/YaHS/tmp ] || mkdir -p {output_directory_path}/HiC/YaHS/tmp
    [ -d {output_directory_path}/HiC/YaHS/no_curation ] || mkdir -p {output_directory_path}/HiC/YaHS/no_curation

    java -Djava.awt.headless=true -jar -Xmx{options['memory']} {juicer_tools} pre \
        -t {output_directory_path}/HiC/YaHS/tmp \
        -j {options['cores']} \
        {alignments_sorted} \
        {output_directory_path}/HiC/YaHS/no_curation/{species_abbreviation(species_name)}.prog.hic \
        {chrom_sizes}

    mv {output_directory_path}/HiC/YaHS/no_curation/{species_abbreviation(species_name)}.prog.hic {outputs['hic']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

# def contact_matrix_manual_curation(JBAT_text: str, chrom_sizes: str, output_directory_path: str, species_name: str, juicer_tools: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/04_genome_assembly/workflow_source/juicer_tools.2.20.00.jar'):
#     """
#     Template: Generate Hi-C contact matrix using :script:`juicer pre`.
    
#     Template I/O::
    
#         inputs = {}
#         outputs = {}
    
#     :param
#     """
#     inputs = {'jbat': JBAT_text,
#               'sizes': chrom_sizes}
#     outputs = {'hic': f'{output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.hic'}
#     protect = outputs['hic']
#     options = {
#         'cores': 32,
#         'memory': '256g',
#         'walltime': '48:00:00'
#     }
#     spec = f"""
#     # Sources environment
#     if [ "$USER" == "jepe" ]; then
#         source /home/"$USER"/.bashrc
#         source activate assembly
#     fi
    
#     echo "START: $(date)"
#     echo "JobID: $SLURM_JOBID"
    
#     [ -d {output_directory_path}/HiC/YaHS/tmp ] || mkdir -p {output_directory_path}/HiC/YaHS/tmp
#     [ -d {output_directory_path}/HiC/YaHS/curation ] || mkdir -p {output_directory_path}/HiC/YaHS/curation

#     java -Djava.awt.headless=true -jar -Xmx{options['memory']} {juicer_tools} pre \
#         -t {output_directory_path}/HiC/YaHS/tmp \
#         -j {options['cores']} \
#         {JBAT_text} \
#         {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.hic \
#         {chrom_sizes}

#     mv {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.hic {outputs['hic']}
    
#     echo "END: $(date)"
#     echo "$(jobinfo "$SLURM_JOBID")"
#     """
#     return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def contact_matrix_manual_curation(JBAT_text: str, JBAT_log: str, output_directory_path: str, species_name: str, juicer_tools: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/04_genome_assembly/workflow_source/juicer_tools.2.20.00.jar'):
    """
    Template: Generate Hi-C contact matrix using :script:`juicer pre`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'jbat': JBAT_text,
              'sizes': JBAT_log}
    outputs = {'hic': f'{output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.hic'}
    protect = outputs['hic']
    options = {
        'cores': 32,
        'memory': '256g',
        'walltime': '48:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory_path}/HiC/YaHS/tmp ] || mkdir -p {output_directory_path}/HiC/YaHS/tmp
    [ -d {output_directory_path}/HiC/YaHS/curation ] || mkdir -p {output_directory_path}/HiC/YaHS/curation

    java -jar -Xmx{options['memory']} {juicer_tools} pre \
        -t {output_directory_path}/HiC/YaHS/tmp \
        -j {options['cores']} \
        {JBAT_text} \
        {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.hic \
        <(cat {JBAT_log}  | grep PRE_C_SIZE | awk '{{print $2" "$3}}')

    mv {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.hic {outputs['hic']}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def post_curation(reviewed_assembly: str, liftover_agp: str, draft_assembly: str, output_directory_path: str, species_name: str):
    """
    Template: Generate :format:`AGP` and :format:`FASTA` files for genome assembly after manual curation with :application:`JuiceBox`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    inputs = {'review': reviewed_assembly,
              'liftover': liftover_agp,
              'draft': draft_assembly}
    outputs = {'final': [f'{output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.FINAL.agp'
                         f'{output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.FINAL.fa']}
    protect = outputs['final']
    options = {
        'cores': 2 ,
        'memory': '30g',
        'walltime': '12:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate assembly
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    juicer post \
        -o {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog \
        {reviewed_assembly} \
        {liftover_agp} \
        {draft_assembly}
    
    mv {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.FINAL.agp {outputs['final'][0]}
    mv {output_directory_path}/HiC/YaHS/curation/{species_abbreviation(species_name)}.JBAT.prog.FINAL.fa {outputs['final'][1]}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)