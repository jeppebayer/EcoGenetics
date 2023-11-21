from gwf import AnonymousTarget
import os

########################## GENERAL FUNCTIONS ##########################

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
        # TODO uncomment next line
        # num += 1
    return chrom_partition

def partition_chrom_real(parse_fasta: list, size: int = 500000):
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

########################## FREEBAYES ##########################

def name_vcf(idx, target):
    chr = os.path.splitext(os.path.basename(target.outputs['region']))[0].split('_', 1)[1]
    return 'VCF_{idx}_{chr}'.format(chr=chr, idx=idx+1)

def vcf_name(idx, target):
    file = os.path.splitext(os.path.basename(target.outputs['region']))[0].replace('-', '_')
    return 'VCF_{}'.format(file)

def create_vcf_per_chr_pooled(region: str, num: int, reference_genome: str, sample_list: str, repeat_regions: str, temp_dir: str, sample_name: str, start: int, end: int, ploidy = 100, bestn = 3):
    """
    Template: Create VCF file for each 'chromosome' in a pooled alignment using :script:`freebayes`, :script:`SnpSift intervals` and :script:`bcftools filter`.

    Template I/O::
        
        inputs = {'region': [reference_genome, repeat_regions]}

        outputs = {'region': *.bcf}

    """
    check_vcf = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/check_vcf_entry.py'
    inputs = {'region': [reference_genome,
              repeat_regions]}
    outputs = {'region': '{temp_dir}/{num}_{sample_name}_{region}.bcf'.format(region=region, num=num, temp_dir=temp_dir, sample_name=sample_name)}
    options = {
        'cores': 1,
        'memory': '100g',
        'walltime': '36:00:00'
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
    
    mv {temp_dir}/{num}_{sample_name}_{region}_prog.bcf {num}_{sample_name}_{region}.bcf

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(reference_genome=reference_genome, bestn=bestn, ploidy=ploidy, region=region, sample_list=sample_list, repeat_regions=repeat_regions, num=num, temp_dir=temp_dir, sample_name=sample_name, start=start, end=end, check_vcf=check_vcf)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def create_vcf_per_chr_individual(region: str, num: int, reference_genome: str, sample_list: str, repeat_regions: str, temp_dir: str, sample_name: str, start: int, end: int, ploidy = 2, bestn = 3):
    """
    Template: Create VCF file for each 'chromosome' in individual alignments using :script:`freebayes`, :script:`SnpSift intervals` and :script:`bcftools filter`.
    
    Template I/O::

        inputs = {'region': [reference_genome, repeat_regions]}

        outputs = {'region': *.bcf}

    """
    check_vcf = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/check_vcf_entry.py'
    inputs = {'region': [reference_genome,
              repeat_regions]}
    outputs = {'region': '{temp_dir}/{num}_{sample_name}_{region}.bcf'.format(region=region, num=num, temp_dir=temp_dir, sample_name=sample_name)}
    options = {
        'cores': 1,
        'memory': '60g',
        'walltime': '96:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    export _JAVA_OPTIONS="-Xms60G -Xmx60G"
    
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

def vcf_concat(regions: list, sample_name: str, temp_dir: str):
    """
    Template: Concatenates `VCF` parts into a single `VCF` file using :script:`bcftools concat`.
    
    Template I/O::

        inputs = {'regions': regions}

        outputs = {'concatenated_vcf': *.vcf}
    
    :param dict regions:
        List of all `VCF` files to concatenate.
    :param str sample_name:
        Name of current sample
    :param str temp_dir:
        Directory path to keep temporary files.
    """
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
    """
    Template: Annotation of :format:`VCF` file using :script:`snpEff`.
    
    Template I/O::

        inputs = {'concatenated_vcf': *.vcf}

        outputs = {'snpeff_vcf': *.ann.vcf,
                'snpeff_stat': snpEff_summary.csv}
    """
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

def vcf_per_chr_pooled_all_rep(reference_genome: str, sample_list: str, repeat_regions: str, working_directory: str, region: str, num: int, start: int, end: int, ploidy = 100, bestn = 3):
    """
    Template: Create VCF file for each 'chromosome' in a pooled alignment using :script:`freebayes`, :script:`SnpSift intervals` and :script:`bcftools filter`.

    Template I/O::
        
        inputs = {'reference': reference_genome, 
                  'repeats': repeat_regions]}

        outputs = {'region': *.bcf}

    """
    check_vcf = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/check_vcf_entry.py'
    inputs = {'reference': reference_genome,
              'repeats': repeat_regions}
    outputs = {'region': '{work_dir}/tmp/{species_abbr}_{num}_{region}.bcf'.format(work_dir=working_directory, species_abbr=species_abbreviation(os.path.basename(os.path.dirname(working_directory))), num=num, region=region)}
    options = {
        'cores': 1,
        'memory': '80g',
        'walltime': '12:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    export _JAVA_OPTIONS="-Xms{mem}G -Xmx{mem}G"

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    [ -d {work_dir}/tmp ] || mkdir -m 775 {work_dir}/tmp

    check=$(freebayes \
        -f {reference_genome} \
        -n {bestn} \
        -p {ploidy} \
        -r {region}:{start}-{end} \
        --min-alternate-fraction 0.05 \
        --min-alternate-count 2 \
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
            -e 'INFO/TYPE~"del" || INFO/TYPE~"ins" || INFO/TYPE~"complex"' \
            -O u \
            - \
            > {bcf}.prog
    else
        freebayes \
            -f {reference_genome} \
            -n {bestn} \
            -p {ploidy} \
            -r {region}:{start}-{end} \
            --min-alternate-fraction 0.05 \
            --min-alternate-count 2 \
            -b {sample_list} \
        | bcftools view \
        -O u \
        > {bcf}.prog
    fi
    
    mv {bcf}.prog {bcf}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(mem=options['memory'], work_dir=working_directory, check_vcf=check_vcf,reference_genome=inputs['reference'], bestn=bestn, ploidy=ploidy, region=region, sample_list=sample_list, repeat_regions=inputs['repeats'], num=num, start=start, end=end, bcf=outputs['region'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def vcf_per_chr_pooled_all_no_rep(reference_genome: str, sample_list: str, working_directory: str, region: str, num: int, start: int, end: int, ploidy = 100, bestn = 3):
    """
    Template: Create VCF file for each 'chromosome' in a pooled alignment using :script:`freebayes` and :script:`bcftools filter`.

    Template I/O::
        
        inputs = {'reference': reference_genome}

        outputs = {'region': *.bcf}

    """
    check_vcf = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/check_vcf_entry.py'
    inputs = {'reference': reference_genome}
    outputs = {'region': '{work_dir}/tmp/{species_abbr}_{num}_{region}.bcf'.format(work_dir=working_directory, species_abbr=species_abbreviation(os.path.basename(os.path.dirname(working_directory))), num=num, region=region)}
    options = {
        'cores': 1,
        'memory': '80g',
        'walltime': '12:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    export _JAVA_OPTIONS="-Xms{mem}G -Xmx{mem}G"

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    [ -d {work_dir}/tmp ] || mkdir -m 775 {work_dir}/tmp

    check=$(freebayes \
        -f {reference_genome} \
        -n {bestn} \
        -p {ploidy} \
        -r {region}:{start}-{end} \
        --min-alternate-fraction 0 \
        --min-alternate-count 2 \
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
            --pooled-discrete  \
            -b {sample_list} \
        | bcftools filter \
            --SnpGap 5 \
            -O u \
            - \
        | bcftools filter \
            -e 'INFO/TYPE~"del" || INFO/TYPE~"ins" || INFO/TYPE~"complex"' \
            -O u \
            - \
            > {bcf}.prog
    else
        freebayes \
            -f {reference_genome} \
            -n {bestn} \
            -p {ploidy} \
            -r {region}:{start}-{end} \
            --min-alternate-fraction 0 \
            --min-alternate-count 2 \
            -b {sample_list} \
        | bcftools view \
        -O u \
        > {bcf}.prog
    fi
    
    mv {bcf}.prog {bcf}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(mem=options['memory'], work_dir=working_directory, check_vcf=check_vcf,reference_genome=inputs['reference'], bestn=bestn, ploidy=ploidy, region=region, sample_list=sample_list, num=num, start=start, end=end, bcf=outputs['region'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def vcf_per_chr_pooled(reference_genome: str, sample_list: list, sample_name: str, output_directory: str, region: str, num: int, start: int, end: int, repeat_regions: str | None = None, ploidy: int = 100, bestn: int = 3):
    """
    Template: Create VCF file for each 'chromosome' in a pooled alignment using :script:`freebayes`, :script:`SnpSift intervals` and :script:`bcftools filter`.

    Template I/O::
        
        inputs = {'reference': reference_genome,
                  'samples': samples_list
                  'repeats': repeat_regions}

        outputs = {'region': *.bcf}

    """
    check_vcf = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/check_vcf_entry.py'
    file_name = '{output_dir}/tmp/{sample_name}_{num}_{region}'.format(output_dir=output_directory, sample_name=sample_name, num=num, region=region)
    if repeat_regions is None:
        repeat_filter = ''
        inputs = {'reference': reference_genome,
                'samples': sample_list}
    else:
        repeat_filter = '| SnpSift intervals -x {} '.format(repeat_regions)
        inputs = {'reference': reference_genome,
                'samples': sample_list,
                'repeats': repeat_regions}
    outputs = {'region': '{}.bcf'.format(file_name)}
    options = {
        'cores': 1,
        'memory': '80g',
        'walltime': '18:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    export _JAVA_OPTIONS="-Xms80G -Xmx80G"

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    [ -d {output_directory}/tmp ] || mkdir -p {output_directory}/tmp

    check=$(freebayes \
        -f {reference_genome} \
        -n {bestn} \
        -p {ploidy} \
        -r {region}:{start}-{end} \
        --min-alternate-fraction 0 \
        --min-alternate-count 2 \
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
            --pooled-discrete  \
            -b {sample_list} \
        {repeat_filter}| bcftools filter \
            --SnpGap 5 \
            -O u \
            - \
        | bcftools filter \
            -e 'INFO/TYPE~"del" || INFO/TYPE~"ins" || INFO/TYPE~"complex"' \
            -O u \
            - \
            > {file_name}.prog.bcf
    else
        freebayes \
            -f {reference_genome} \
            -n {bestn} \
            -p {ploidy} \
            -r {region}:{start}-{end} \
            --min-alternate-fraction 0 \
            --min-alternate-count 2 \
            -b {sample_list} \
        | bcftools view \
        -O u \
        > {file_name}.prog.bcf
    fi
    
    mv {file_name}.prog.bcf {bcf_file}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(output_directory=output_directory, reference_genome=reference_genome, sample_list=' -b '.join(sample_list), repeat_filter=repeat_filter, file_name=file_name, bcf_file=outputs['region'], region=region, num=num, start=start, end=end, bestn=bestn, ploidy=ploidy, check_vcf=check_vcf)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def vcf_per_chr_individual(reference_genome: str, sample_list: list, sample_name: str, output_directory: str, region: str, num: int, start: int, end: int, repeat_regions: str | None = None, ploidy: int = 2, bestn: int = 3):
    """
    Template: Create VCF file for each 'chromosome' in individual alignments using :script:`freebayes`, :script:`SnpSift intervals` and :script:`bcftools filter`.
    
    Template I/O::

        inputs = {'reference': reference_genome,
                  'samples': samples_list
                  'repeats': repeat_regions}

        outputs = {'region': *.bcf}

    """
    check_vcf = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/check_vcf_entry.py'
    file_name = '{output_dir}/tmp/{sample_name}_{num}_{region}'.format(output_dir=output_directory, sample_name=sample_name, num=num, region=region)
    if repeat_regions is None:
        repeat_filter = ''
        inputs = {'reference': reference_genome,
                'samples': sample_list}
    else:
        repeat_filter = '| SnpSift intervals -x {} '.format(repeat_regions)
        inputs = {'reference': reference_genome,
                'samples': sample_list,
                'repeats': repeat_regions}
    outputs = {'region': '{}.bcf'.format(file_name)}
    options = {
        'cores': 1,
        'memory': '80g',
        'walltime': '18:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    export _JAVA_OPTIONS="-Xms80G -Xmx80G"
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory}/tmp ] || mkdir -p {output_directory}/tmp

    check=$(freebayes \
        -f {reference_genome} \
        -n {bestn} \
        -p {ploidy} \
        -r {region}:{start}-{end} \
        --min-alternate-fraction 0.05 \
        --min-alternate-count 2 \
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
            -b {sample_list} \
        {repeat_filter}| bcftools filter \
            --SnpGap 5 \
            -O u \
            - \
        | bcftools filter \
            -e 'INFO/TYPE~"del" || INFO/TYPE~"ins" || INFO/TYPE~"complex"' \
            -O u \
            - \
            > {file_name}.prog.bcf
    else
        freebayes \
            -f {reference_genome} \
            -n {bestn} \
            -p {ploidy} \
            -r {region}:{start}-{end} \
            --min-alternate-fraction 0.05 \
            --min-alternate-count 2 \
            -b {sample_list} \
        | bcftools view \
        -O u \
        > {file_name}.prog.bcf
    fi
    
    mv {file_name}.prog.bcf {bcf_file}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(output_directory=output_directory, reference_genome=reference_genome, sample_list=' -b '.join(sample_list), repeat_filter=repeat_filter, file_name=file_name, bcf_file=outputs['region'], region=region, num=num, start=start, end=end, bestn=bestn, ploidy=ploidy, check_vcf=check_vcf)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def concatenate_vcf(regions: list, output_directory: str, species_name: str):
    """
    Template: Concatenates `VCF` parts into a single `VCF` file using :script:`bcftools concat`.
    
    Template I/O::

        inputs = {'regions': regions}

        outputs = {'concat_vcf': *.vcf}
    
    :param list regions:
        List of all `VCF` files to concatenate.
    :param str output_directory:
        Directory to output resulting :format:`VCF` file.
    :param str species_name:
        Name of species.
    """
    inputs = {'regions': regions}
    outputs = {'concat_vcf': '{output_dir}/{species_abbr}.vcf'.format(output_dir=output_directory, species_abbr=species_abbreviation(species_name))}
    protect = [outputs['concat_vcf']]
    options = {
        'cores': 8,
        'memory': '96g',
        'walltime': '10:00:00'
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
        -o {vcf}.prog \
        -O v \
        {regions}
    
    mv {vcf}.prog {vcf}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(vcf=outputs['concat_vcf'], regions=' '.join(inputs['regions']))
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def concat_vcf(regions: list, output_name: str, output_directory: str | None = None):
    """
    Template: Concatenates `VCF` parts into a single `VCF` file using :script:`bcftools concat`.
    
    Template I/O::

        inputs = {'regions': regions}

        outputs = {'concat_vcf': *.vcf}
    
    :param dict regions:
        List of all `VCF` files to concatenate.
    :param str sample_name:
        Name of current sample
    :param str temp_dir:
        Directory path to keep temporary files.
    """
    if output_directory is None:
        output_directory = os.path.dirname(regions[0])
    file_name = '{output_dir}/{name}'.format(output_dir=output_directory, name=output_name)
    regions.sort()
    inputs = {'regions': regions}
    outputs = {'concat_vcf': '{}.vcf'.format(file_name)}
    protect = [outputs['concat_vcf']]
    options = {
        'cores': 8,
        'memory': '96g',
        'walltime': '10:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory} ] || mkdir -p {output_directory}

    bcftools concat \
        --threads 8 \
        -o {file_name}.prog.vcf \
        -O v \
        {regions}
    
    mv {file_name}.prog.vcf {concat_vcf}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(output_directory=output_directory, file_name=file_name, concat_vcf=outputs['concat_vcf'], regions=' '.join(inputs['regions']))
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def snpEff_ann(vcf_file: str, reference_genome_version: str, output_directory: str | None = None):
    """
    Template: Annotation of :format:`VCF` file using :script:`snpEff`.
    
    Template I/O::

        inputs = {'vcf': vcf_file}

        outputs = {'snpeff_vcf': output_directory/*.ann.vcf,
                   'snpeff_stat': output_directory/snpEff_summary.csv}
    
    :param str vcf_file:
        :format:`VCF` file to annotate.
    :param str reference_genome_version:
        Name of entry in :script:`snpEff` configuration file.
    :param str output_directory:
        Directory for resulting annotated :format:`VCF` file and statistics file.
    """
    # TODO make snpEff configuration file available in workflow_source directory. Might only need to contain custom entries.
    snpeff_config = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/snpEff/snpEff.config'
    snpeff_data = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/snpEff/data'
    if output_directory is None:
        output_directory = os.path.dirname(vcf_file)
    file_name = '{output_dir}/{filename}'.format(output_dir=output_directory, filename=os.path.splitext(os.path.basename(vcf_file))[0])
    inputs = {'vcf': vcf_file}
    outputs = {'snpeff_vcf': '{}.ann.vcf'.format(file_name),
               'snpeff_stat': '{}/snpEff_summary.csv'.format(output_directory)}
    protect = [outputs['snpeff_stat'], outputs['snpeff_vcf']]
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

    [ -d {output_directory} ] || mkdir -p {output_directory}

    snpEff ann \
        -csvStats {stats} \
        -c {snpEff_config} \
        -dataDir {snpEff_data} \
        {reference_genome_version} \
        {vcf} \
        > {file_name}.prog.ann.vcf

    mv {file_name}.prog.ann.vcf {ann_vcf}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(stats=outputs['snpeff_stat'], snpeff_config=snpeff_config, snpEff_data=snpeff_data, reference_genome_version=reference_genome_version, vcf=vcf_file, file_name=file_name, ann_vcf=outputs['snpeff_vcf'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

########################## PoolSNP ##########################

def name_mpileup(idx, target):
    return 'mpileup_{idx}'.format(idx=idx+1)

def name_sync(idx, target):
    return 'sync_{idx}'.format(idx=idx+1)

def name_cov(idx, target):
    return '{}'.format(os.path.basename(target.outputs['cutoff']))

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
        source activate vcf
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    sleep 2m

    [ -d {output_directory} ] || mkdir -p {output_directory}

    echo -e '{chromosome}\t{start}\t{end}' > {output_directory}/{num}.bed
    
    samtools mpileup \
        --max-depth 0 \
        --fasta-ref {reference} \
        --positions {output_directory}/{num}.bed \
        --min-BQ 0 \
        --region {chromosome} \
        --output {file_name}.prog.mpileup \
        {bam_files}
    
    mv {file_name}.prog.mpileup {mpileup}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(output_directory=output_directory, chromosome=region, start=start, end=end, num=num, reference=reference_genome, file_name=file_name, bam_files=bam_string, mpileup=outputs['mpileup'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mpileup2sync(mpileup_file: str, output_directory: str = None, mpileup2sync: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/popoolation2_v1.201/mpileup2sync.pl'):
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
        source activate vcf
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory} ] || mkdir -p {output_directory}

    perl {mpileup2sync} \
        --input {mpileup} \
        --output {file_name}.prog.sync \
        --fastq-type sanger
    
    mv {file_name}.prog.sync {sync}
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(mpileup2sync=mpileup2sync, mpileup=mpileup_file, file_name=file_name, sync=outputs['sync'], output_directory=output_directory)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def max_cov(mpileup: str, contig: str, cutoff: float, output_directory: str, script: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/PoolSNP/scripts/max-cov.py'):
    """
    Template: Calculates coverage thresholds using :script:`max-cov.py`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    file_name = '{output_directory}/tmp/cov/cutoffs/{contig}'.format(output_directory=output_directory, contig=contig)
    inputs = {'mpileup': mpileup}
    outputs = {'cutoff': '{}.txt'.format(file_name)}
    options = {
        'cores': 1,
        'memory': '10g',
        'walltime': '30:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory}/tmp/cov/cutoffs ] || mkdir -p {output_directory}/tmp/cov/cutoffs

    awk \
        -v contig={contig} \
        '{{if ($1 == contig) {{print $0}}}}' \
        {mpileup} \
    | python {script} \
        --mpileup - \
        --cutoff {cutoff} \
        --contig {contig} \
        --out {file_name}.prog.txt
    
    if [ -f {file_name}.prog.txt ]; then
        mv {file_name}.prog.txt {cutoff_file}
    else
        echo -n "" > {cutoff_file}
    fi
    
    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(output_directory=output_directory, contig=contig, mpileup=mpileup, script=script, cutoff=cutoff, file_name=file_name, cutoff_file=outputs['cutoff'])
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
    files.sort()
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
        source activate vcf
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {output_directory}] || mkdir -p {output_directory}

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
    """.format(compress=compress, sorted_files=' '.join(files), file_name=file_name, ext=os.path.splitext(files[0])[1], concat_file=outputs['concat_file'], output_directory=output_directory)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, protect=protect, spec=spec)

def poolsnp(mpileup: str, max_cov: str, sample_list: list, reference_genome: str, working_directory: str, species_name: str, output_directory: str = None, min_cov: int = 10, min_count: int = 3, min_freq: float = 0.01, miss_frac: float = 0.1, bq: int = 15, sites: int = 1, script: str = '/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/gwf/03_initial_analysis_files/workflow_source/PoolSNP/scripts/PoolSnp.py'):
    """
    Template: Creates :format:`VCF` file using :script:`PoolSnp.py`.
    
    Template I/O::
    
        inputs = {}
        outputs = {}
    
    :param
    """
    if output_directory is None:
        output_directory = working_directory
    inputs = {'mpileup': mpileup,
              'max_cov': max_cov}
    outputs = {'vcf': '{output_directory}/{species_abbr}.mincov2.vcf.gz'.format(output_directory=output_directory, species_abbr=species_abbreviation(species_name))}
    protect = outputs['vcf']
    options = {
        'cores': 20,
        'memory': '160g',
        'walltime': '60:00:00'
    }
    spec = """
    # Sources environment
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi
    
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    [ -d {working_directory}/tmp ] || mkdir -p {working_directory}/tmp

    headerfile={working_directory}/tmp/header.txt
    echo -e "##fileformat=VCFv4.2" > "$headerfile"
    echo -e "##fileDate=$(date +%d'/'%m'/'%y)" >> "$headerfile"
    echo -e "##Source=PoolSnp-1.05" >> "$headerfile"
    echo -e "##Parameters=<ID=MinCov,Number={min_cov},Type=Integer,Description=\"Minimum coverage per sample\">" >> "$headerfile"
    echo -e "##Parameters=<ID=MaxCov,Number={max_cov},Type=Integer,Description=\"Maximum chromosome- and sample-specific maximum coverage; Either a precomputed file or the maximum percentile cutoff, eg. 0.95 to consider only reads within the 95% coverage percentile\">" >> "$headerfile"
    echo -e "##Parameters=<ID=MinCount,Number={min_count},Type=Integer,Description=\"Minimum alternative allele count across all samples pooled\">" >> "$headerfile"
    echo -e "##Parameters=<ID=MinFreq,Number={min_freq},Type=Float,Description=\"Minimum alternative allele frequency across all samples pooled\">" >> "$headerfile"
    echo -e "##Parameters=<ID=MaximumMissingFraction,Number={miss_frac},Type=Float,Description=\"Maximum fraction of samples allowed that are not fullfilling all parameters\">" >> "$headerfile"
    echo -e "##Parameters=<ID=BaseQual,Number={bq},Type=Integer,Description=\"Minimum PHRED scaled base quality\">" >> "$headerfile"
    echo -e "##Reference={reference}" >> "$headerfile"
    echo -e "##INFO=<ID=ADP,Number=1,Type=Integer,Description=\"Average per-sample depth of bases with Phred score >={bq}\">" >> "$headerfile"
    echo -e "##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of samples not called\">" >> "$headerfile"
    echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> "$headerfile"
    echo -e "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Reference Counts\">" >> "$headerfile"
    echo -e "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Alternative Counts\">" >> "$headerfile"
    echo -e "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" >> "$headerfile"
    echo -e "##FORMAT=<ID=FREQ,Number=1,Type=Float,Description=\"Variant allele frequency\">" >> "$headerfile"
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}" >> "$headerfile"

    parallel \
        -k \
        -j {cores} \
        --pipepart \
        --no-notice \
        -a {mpileup} \
    --cat python {script} \
        --mpileup {{}} \
        --min-cov {min_cov} \
        --max-cov {max_cov} \
        --min-freq {min_freq} \
        --miss-frac {miss_frac} \
        --min-count {min_count} \
        --base-quality  {bq} \
        --allsites {sites} \
        >  {working_directory}/tmp/SNPs.prog.txt

    mv {working_directory}/tmp/SNPs.prog.txt {working_directory}/tmp/SNPs.txt

    cat {working_directory}/tmp/header.txt {working_directory}/tmp/SNPs.txt > {VCF}

    rm -f {working_directory}/tmp/SNPs.prog.txt
    rm -f {working_directory}/tmp/SNPs.txt

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(cores=options['cores'], mpileup=mpileup, script=script, working_directory=working_directory, output_directory=output_directory, min_cov=min_cov, max_cov=max_cov, min_count=min_count, min_freq=min_freq, miss_frac=miss_frac, bq=bq, reference=reference_genome, sites=sites, VCF=outputs['vcf'], samples='\t'.join(sample_list))
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)