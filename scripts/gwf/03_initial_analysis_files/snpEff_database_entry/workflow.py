from gwf import Workflow
import os
import yaml
import glob
import sys
from workflow_templates import species_abbreviation, gff_to_gtf, snpEff_entry_text, create_snpEff_entry

# --------------------------------------------------
#                  Configuration
# --------------------------------------------------

config=yaml.safe_load(open('SteBic_config.yaml'))
ACCOUNT = config['account']
SPECIES_NAME = config['species_name']
REFERENCE_GENOME = config['reference_genome']
ANNOTATION_FILE = config['annotation_file']

# --------------------------------------------------
#                  Workflow
# --------------------------------------------------

gwf = Workflow(
    defaults={'account': ACCOUNT}
)

gwf.target_from_template(
    name='GTF_' + species_abbreviation(SPECIES_NAME),
    template=gff_to_gtf(
        annotation_name=os.path.splitext(ANNOTATION_FILE)[0]
        )
)

# Get path to snpEff directory
snpEff_dir = glob.glob(sys.prefix + '/share/snpeff*')[0]

# Creates backup of snpEff config file if it doesn't exist
if not os.path.exists(snpEff_dir + '/snpEff.config.bak'):
    os.makedirs(snpEff_dir + '/snpEff.config.bak')

reference_genome_version = os.path.splitext(os.path.basename(REFERENCE_GENOME))[0]
snpEff_species_dir = snpEff_dir + '/data/' + reference_genome_version

# Create data folder and genome folder in snpEff directory if it doesn't exist
if not os.path.exists(snpEff_species_dir):
    os.makedirs(snpEff_species_dir)

# Checks if entry text is already in snpEff.config and only adds the entry if it doesn't already exist
entry_text = snpEff_entry_text(REFERENCE_GENOME, SPECIES_NAME)
with open(snpEff_dir + '/snpEff.config', 'r') as database:
    if not entry_text in database.read():
        with open(snpEff_dir + '/snpEff.config', 'a') as database:
            database.write(entry_text)

gwf.target_from_template(
    name='snpEff_' + species_abbreviation(SPECIES_NAME),
    template=create_snpEff_entry(
        reference_genome=REFERENCE_GENOME,
        annotation_file=ANNOTATION_FILE,
        snpEff_species_dir=snpEff_species_dir,
        snpEff_dir=snpEff_dir
        )
)

# Makes file with custom database entries listed for easier access
custom_entries=snpEff_dir + '/custom_entries'
species_entry="{species_name}\t{reference_genome_version}".format(species_name=SPECIES_NAME, reference_genome_version=reference_genome_version)
if not os.path.exists(custom_entries):
    os.makedirs(custom_entries)
    with open(custom_entries, 'a') as ce:
        ce.write("#Organism name\tGenome name")
with open(custom_entries, 'r') as ce:
    if not species_entry in ce.read():
        with open(custom_entries, 'a') as ce:
            ce.write(species_entry)