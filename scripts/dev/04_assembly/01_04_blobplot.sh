#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --time=30
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=2
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

WD="$1"
speciesabbr="$2"
rank="$3"

if [ "$USER" == "jepe" ]; then
    # shellcheck disable=1090
    source /home/"$USER"/.bashrc
    # shellcheck disable=1091
    source activate blobtools
fi

phylum="Mucoromycota,Zoopagomycota,Chordata,Proteobacteria,Mollusca,Rotifera,Chytridiomycota"
order="Entomophthorales,Basidiobolales,Glomerales,Cypriniformes,Zoopagales,Pseudomonadales,Mucorales,Dimargaritales,Ostreida,Spizellomycetales,Mytilida,Anura,other"

color=/faststorage/project/EcoGenetics/people/Jeppe_Bayer/scripts/04_assembly/blobtools/modules/colorfile.txt

if [ "$rank" == "phylum" ]; then
    excl=$phylum
fi

if [ "$rank" == "order" ]; then
    excl=$order
fi

# blobtools plot \
# -i "$WD"/"$speciesabbr"_blobtools_db.blobDB.json \
# -o "$WD"/ \
# -x bestsumorder \
# -r "$rank" \
# -n \
# -p 15 \
# --colours "$color" \
# --exclude "$excl"

blobtools plot \
-i "$WD"/"$speciesabbr"_blobtools_db.blobDB.json \
-o "$WD"/ \
-x bestsumorder \
-r "$rank" \
-n \
-p 1 \
--exclude "other"

exit 0