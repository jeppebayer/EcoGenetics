from gwf import AnonymousTarget
import os, sys

def snpEff_entry_text(reference_genome, species_name):
    """Function to create entry text for snpEff database"""
    reference_genome_version = os.path.splitext(os.path.basename(reference_genome))[0]
    entry_text = """
# {species_name} genome, version {reference_genome_version}, entry added by Centre for EcoGenetics (AU)
{reference_genome_version}.genome : {species_name}""".format(reference_genome_version=reference_genome_version, species_name=species_name)
    return entry_text

def species_abbreviation(species_name: str):
    """Creates species abbreviation from species name.
    
    :param str species_name:
        Species name written as *genus* *species*"""
    genus, species = species_name.replace(' ', '_').split('_')
    genus = genus[0].upper() + genus[1:3]
    species = species[0].upper() + species[1:3]
    return genus + species

def gff_to_gtf(annotation_name):
    """Tempalte for converting GFF annotation file to GTF 2.2 annotation file."""
    inputs = []
    outputs = ['{}.gtf'.format(annotation_name)]
    protect = ['{}.gtf'.format(annotation_name)]
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '01:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    agat_convert_sp_gff2gtf.pl \
        -gff {annotation_name}.gtf \
        --gtf_version 2.2 \
        -o {annotation_name}.gtf

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID"")"
    """.format(annotation_name=annotation_name)
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def create_snpEff_entry(reference_genome, annotation_file, snpEff_species_dir):
    """Template for creating a new snpEff database entry."""
    inputs = [reference_genome,
              '{}.gtf'.format(os.path.splitext(annotation_file)[0])]
    outputs = ['{snpEff_species_dir}/cds.fa'.format(snpEff_species_dir=snpEff_species_dir),
               '{snpEff_species_dir}/protein.fa'.format(snpEff_species_dir=snpEff_species_dir),
               '{snpEff_species_dir}/genes.gtf'.format(snpEff_species_dir=snpEff_species_dir),
               '{snpEff_species_dir}/sequences.fa'.format(snpEff_species_dir=snpEff_species_dir),
               '{snpEff_species_dir}/snpEffectPredictor.bin'.format(snpEff_species_dir=snpEff_species_dir)
               ]
    options = {
        'cores': 1,
        'memory': '100g',
        'walltime': '01:00:00'
    }
    spec = """
    if [ "$USER" == "jepe" ]; then
        source /home/"$USER"/.bashrc
        source activate vcf
    fi

    export _JAVA_OPTIONS="-Xms100G -Xmx100G"

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    agat_sp_extract_sequences.pl \
        -g {gtf22_file} \
        -f {reference_genome} \
        -t cds \
        -o {snpEff_species_dir}/cds.fa
    
    agat_sp_extract_sequences.pl \
        -g {gtf22_file} \
        -f {reference_genome} \
        -t cds \
        -p \
        -o {snpEff_species_dir}/protein.fa
    
    cp {gtf22_file} {snpEff_species_dir}/genes.gtf
    cp {reference_genome} {snpEff_species_dir}/sequences.fa
    
    snpEff build -gtf22 -v {reference_genome_version}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID"")"
    """.format(reference_genome=reference_genome,
               reference_genome_version=os.path.splitext(os.path.basename(reference_genome))[0],
               gtf22_file='{}.gtf'.format(os.path.splitext(annotation_file)[0]),
               snpEff_species_dir=snpEff_species_dir)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)