#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --chdir=/faststorage/project/EcoGenetics/BACKUP/database/uniprot
#SBATCH --mem 20G
#SBATCH --cpus-per-task 2
#SBATCH --time 10:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Reference_Proteomes_2022_05.tar.gz

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate ecogen_primary
fi

# Took 5 days
wget -q -i URLs.txt

# Unpack archive
for file in ./*.tar.gz; do
    tar -zxf "$file" -C .
done

# Unpack protein FASTAs for each kingdom
for file in ./*/*/*.fasta.gz; do
    if [[ ! "$file" == *DNA.fasta.gz ]];then
        if [[ ! "$file" == *additional.fasta.gz ]]; then
            gzip -d "$file"
        fi
    fi
done

for file in ./*/*/*.idmapping.gz; do
    gzip -d "$file"
done

# Concatenate all protein sequences into 'uniprot_ref_proteomes.fasta'
cat ./*/*/*.fasta > uniprot_ref_proteomes.fasta

# Simplyfy sequence IDs
cat uniprot_ref_proteomes.fasta | sed -r 's/(^>sp\|)|(^>tr\|)/>/g' | cut -f1 -d"|" > temp; mv temp uniprot_ref_proteomes.fasta

# Subset mapping files to only contain NCBI TaxID entries
cat ./*/*/*.idmapping | rg -j 2 "NCBI_TaxID" > uniprot_ref_proteomes.taxids

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate blobtools
fi

# Make Diamond DB
diamond makedb \
--threads 2 \
--in uniprot_ref_proteomes.fasta \
-d uniprot_ref_proteomes

exit 0