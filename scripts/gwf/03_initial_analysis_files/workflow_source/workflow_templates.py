from gwf import AnonymousTarget
import os

def species_abbreviation(species_name):
    """Function for creating species abbreviation from species name."""
    genus, species = species_name.split()
    genus = genus[0].upper() + genus[1:3]
    species = species[0].upper() + species[1:3]
    return genus + species

def parse_fasta(fasta_file):
    """Function to parse FASTA file, retrieving all sequence names and lengths and returning them paired in a list of dictionaries."""
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

def partition_chrom(parse_fasta, size = 500000):
    """Function to get partitioning of parsed fasta. Creates a list of dictionaries containing 'num', 'region', 'start', 'end'."""
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
        # num += 1 # Forgot this initially
    return chrom_partition

def name_vcf(idx, target):
    chr = os.path.splitext(os.path.basename(target.outputs['region']))[0].split('_', 1)[1]
    return 'VCF_{idx}_{chr}'.format(chr=chr, idx=idx+1)

def create_vcf_per_chr_pooled(region, num, reference_genome, sample_list, repeat_regions, temp_dir, sample_name, start: int, end: int, ploidy = 100, bestn = 3):
    """Template for creating a VCF file for each 'chromosome' in a pooled sample BAM"""
    check_vcf = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/check_vcf_entry.py'
    inputs = {'region': [reference_genome,
              repeat_regions]}
    outputs = {'region': '{temp_dir}/{num}_{sample_name}_{region}.bcf'.format(region=region, num=num, temp_dir=temp_dir, sample_name=sample_name)}
    options = {
        'cores': 1,
        'memory': '100g',
        'walltime': '18:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    export _JAVA_OPTIONS="-Xms100G -Xmx100G"

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    check=$(freebayes \
        -f {reference_genome} \
        -n {bestn} \
        -p {ploidy} \
        -r {region}:{start}-{end} \
        --min-alternate-fraction 0.05 \
        --min-alternate-count 2 \
        --report-monomorphic \
        -b {sample_list} \
    | {check_vcf})

    if [ "$check" -eq 1 ]; then
        freebayes \
            -f {reference_genome} \
            -n {bestn} \
            -p {ploidy} \
            -r {region}:{start}-{end} \
            --min-alternate-fraction 0 \
            --min-alternate-count 2 \
            --report-monomorphic \
            --pooled-discrete  \
            -b {sample_list} \
        | SnpSift intervals \
            -x \
            {repeat_regions} \
        | bcftools filter \
            --SnpGap 5 \
            -O u \
            - \
        | bcftools filter \
            -e 'TYPE~"del" || TYPE~"ins" || TYPE~"complex"' \
            -O u \
            - \
            > {temp_dir}/{num}_{sample_name}_{region}_prog.bcf
    else
        freebayes \
            -f {reference_genome} \
            -n {bestn} \
            -p {ploidy} \
            -r {region}:{start}-{end} \
            --min-alternate-fraction 0.05 \
            --min-alternate-count 2 \
            --report-monomorphic \
            -b {sample_list} \
        | bcftools view \
        -O u \
        >{temp_dir}/{num}_{sample_name}_{region}_prog.bcf
    fi
    
    mv {temp_dir}/{num}_{sample_name}_{region}_prog.bcf {num}_{sample_name}_{region}.bcf

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(reference_genome=reference_genome, bestn=bestn, ploidy=ploidy, region=region, sample_list=sample_list, repeat_regions=repeat_regions, num=num, temp_dir=temp_dir, sample_name=sample_name, start=start, end=end, check_vcf=check_vcf)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def create_vcf_per_chr_individual(region, num, reference_genome, sample_list, repeat_regions, temp_dir, sample_name, start: int, end: int, ploidy = 2, bestn = 3):
    """Template for creating a VCF file for each 'chromosome' in a individual sample BAM"""
    check_vcf = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/check_vcf_entry.py'
    inputs = {'region': [reference_genome,
              repeat_regions]}
    outputs = {'region': '{temp_dir}/{num}_{sample_name}_{region}.bcf'.format(region=region, num=num, temp_dir=temp_dir, sample_name=sample_name)}
    options = {
        'cores': 1,
        'memory': '100g',
        'walltime': '18:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    export _JAVA_OPTIONS="-Xms100G -Xmx100G"
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    check=$(freebayes \
        -f {reference_genome} \
        -n {bestn} \
        -p {ploidy} \
        -r {region}:{start}-{end} \
        --min-alternate-fraction 0.05 \
        --min-alternate-count 2 \
        --report-monomorphic \
        -b {sample_list} \
    | {check_vcf})

    if [ "$check" -eq 1 ]; then
        freebayes \
            -f {reference_genome} \
            -n {bestn} \
            -p {ploidy} \
            -r {region}:{start}-{end} \
            --min-alternate-fraction 0.05 \
            --min-alternate-count 2 \
            --report-monomorphic \
            -b {sample_list} \
        | SnpSift intervals \
            -x \
            {repeat_regions} \
        | bcftools filter \
            --SnpGap 5 \
            -O u \
            - \
        | bcftools filter \
            -e 'INFO/TYPE~"del" || INFO/TYPE~"ins" || INFO/TYPE~"complex"' \
            -O u \
            - \
            > {temp_dir}/{num}_{sample_name}_{region}_prog.bcf
    else
        freebayes \
            -f {reference_genome} \
            -n {bestn} \
            -p {ploidy} \
            -r {region}:{start}-{end} \
            --min-alternate-fraction 0.05 \
            --min-alternate-count 2 \
            --report-monomorphic \
            -b {sample_list} \
        | bcftools view \
        -O u \
        >{temp_dir}/{num}_{sample_name}_{region}_prog.bcf
    fi
    
    mv {temp_dir}/{num}_{sample_name}_{region}_prog.bcf {temp_dir}/{num}_{sample_name}_{region}.bcf

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(reference_genome=reference_genome, bestn=bestn, ploidy=ploidy, region=region, sample_list=sample_list, repeat_regions=repeat_regions, num=num, temp_dir=temp_dir, sample_name=sample_name, start=start, end=end, check_vcf=check_vcf)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def vcf_concat(regions, sample_name, temp_dir):
    """Template for concatenating VCF parts into one VCF file."""
    inputs = {'regions': regions}
    outputs = {'concatenated_vcf': '{temp_dir}/{sample_name}.vcf'.format(temp_dir=temp_dir, sample_name=sample_name)}
    protect = [outputs['concatenated_vcf']]
    options = {
        'cores': 8,
        'memory': '96g',
        'walltime': '05:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    bcftools concat \
        --threads 8 \
        -o {temp_dir}/{sample_name}_prog.vcf \
        -O v \
        {regions}
    
    mv {temp_dir}/{sample_name}_prog.vcf {temp_dir}/{sample_name}.vcf

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(temp_dir=temp_dir, sample_name=sample_name, regions=' '.join(inputs['regions']))
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def snpEff_annotation(temp_dir, work_dir, sample_name, snpEff_config, reference_genome_version):
    """Template for annotation of VCF file using snpEff."""
    inputs = {'concatenated_vcf': '{temp_dir}/{sample_name}.vcf'.format(temp_dir=temp_dir, sample_name=sample_name)}
    outputs = {'snpeff_vcf': '{work_dir}/{sample_name}.ann.vcf'.format(work_dir=work_dir, sample_name=sample_name),
               'snpeff_stat': '{work_dir}/snpEff_summary.csv'.format(work_dir=work_dir)}
    options = {
        'cores': 2,
        'memory': '50g',
        'walltime': '03:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    snpEff ann \
        -csvStats {work_dir}/snpEff_summary.csv \
        -c {snpEff_config} \
        {reference_genome_version} \
        {temp_dir}/{sample_name}.vcf \
        > {work_dir}/{sample_name}_prog.ann.vcf

    mv {work_dir}/{sample_name}_prog.ann.vcf {work_dir}/{sample_name}.ann.vcf

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(temp_dir=temp_dir, work_dir=work_dir, sample_name=sample_name, reference_genome_version=reference_genome_version, snpEff_config=snpEff_config)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)