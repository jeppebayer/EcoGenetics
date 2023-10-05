from gwf import AnonymousTarget
import os

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

def name_vcf(idx, target):
    chr = os.path.splitext(os.path.basename(target.outputs['region']))[0].split('_', 1)[1]
    return 'VCF_{idx}_{chr}'.format(chr=chr, idx=idx+1)

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
    """.format(mem=options['memory'], work_dir=working_directory, check_vcf=check_vcf,reference_genome=inputs['reference'], bestn=bestn, ploidy=ploidy, region=region, sample_list=sample_list, num=num, start=start, end=end, bcf=outputs['region'])
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

def snpEff_ann(vcf_file: str, reference_genome_version: str, output_directory: str):
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
    inputs = {'vcf': vcf_file}
    outputs = {'snpeff_vcf': '{output_dir}/{file_name}.ann.vcf'.format(output_dir=output_directory, file_name=os.path.splitext(os.path.basename(inputs['vcf']))[0]),
               'snpeff_stat': '{output_dir}/snpEff_summary.csv'.format(output_dir=output_directory)}
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
        -csvStats {stats} \
        -c {snpEff_config} \
        -dataDir {snpEff_data} \
        {reference_genome_version} \
        {vcf} \
        > {ann_vcf}.prog

    mv {ann_vcf}.prog {ann_vcf}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """.format(stats=outputs['snpeff_stat'], snpeff_config=snpeff_config, snpEff_data=snpeff_data, reference_genome_version=reference_genome_version, vcf=inputs['vcf'], ann_vcf=outputs['snpeff_vcf'])
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)