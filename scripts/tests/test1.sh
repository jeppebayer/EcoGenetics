#!/bin/bash

for ((i=1; i<=5; i++)); do
    bash people/Jeppe_Bayer/scripts/02_data_preparation/02_data_preparation_init.sh \
    -s people/Jeppe_Bayer/X_chromosome_social_spiders/mutation_recombination/S_dumicola/family"$i" \
    -r people/Jeppe_Bayer/X_chromosome_social_spiders/X_Evolution/reference_genomes/S_dumicola/DUM_hifi_hic_scaffolded_trim.fa \
    -d people/Jeppe_Bayer/X_chromosome_social_spiders/mutation_recombination/S_dumicola
done

for ((i=1; i<=5; i++)); do
    bash people/Jeppe_Bayer/scripts/02_data_preparation/02_data_preparation_init.sh \
    -s people/Jeppe_Bayer/X_chromosome_social_spiders/mutation_recombination/S_lineatus/family"$i" \
    -r people/Jeppe_Bayer/X_chromosome_social_spiders/X_Evolution/reference_genomes/S_lineatus/LIN_hifi_hic_scaffolded_trim.fa \
    -d people/Jeppe_Bayer/X_chromosome_social_spiders/mutation_recombination/S_lineatus
done

for ((i=1; i<=5; i++)); do
    bash people/Jeppe_Bayer/scripts/02_data_preparation/02_data_preparation_init.sh \
    -s people/Jeppe_Bayer/X_chromosome_social_spiders/mutation_recombination/S_mimosarum/family"$i" \
    -r people/Jeppe_Bayer/X_chromosome_social_spiders/X_Evolution/reference_genomes/S_mimosarum/MIM_hifi_hic_scaffolded_trim.fa \
    -d people/Jeppe_Bayer/X_chromosome_social_spiders/mutation_recombination/S_mimosarum
done

for ((i=1; i<=5; i++)); do
    bash people/Jeppe_Bayer/scripts/02_data_preparation/02_data_preparation_init.sh \
    -s people/Jeppe_Bayer/X_chromosome_social_spiders/mutation_recombination/S_sarasinorum/family"$i" \
    -r people/Jeppe_Bayer/X_chromosome_social_spiders/X_Evolution/reference_genomes/S_sarasinorum/SARA_hifi_hic_scaffolded_trim.fa \
    -d people/Jeppe_Bayer/X_chromosome_social_spiders/mutation_recombination/S_sarasinorum
done

for ((i=1; i<=5; i++)); do
    bash people/Jeppe_Bayer/scripts/02_data_preparation/02_data_preparation_init.sh \
    -s people/Jeppe_Bayer/X_chromosome_social_spiders/mutation_recombination/S_tentoriicola/family"$i" \
    -r people/Jeppe_Bayer/X_chromosome_social_spiders/X_Evolution/reference_genomes/S_tentoriicola/TENT_hifi_hic_scaffolded_trim.fa \
    -d people/Jeppe_Bayer/X_chromosome_social_spiders/mutation_recombination/S_tentoriicola
done

bash people/Jeppe_Bayer/scripts/02_data_preparation/02_data_preparation_init.sh \
-s people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_dumicola/Betta \
-r people/Jeppe_Bayer/X_chromosome_social_spiders/X_Evolution/reference_genomes/S_dumicola/DUM_hifi_hic_scaffolded_trim.fa \
-d people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_dumicola

bash people/Jeppe_Bayer/scripts/02_data_preparation/02_data_preparation_init.sh \
-s people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_dumicola/Karasburg \
-r people/Jeppe_Bayer/X_chromosome_social_spiders/X_Evolution/reference_genomes/S_dumicola/DUM_hifi_hic_scaffolded_trim.fa \
-d people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_dumicola

bash people/Jeppe_Bayer/scripts/02_data_preparation/02_data_preparation_init.sh \
-s people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_dumicola/Otawi \
-r people/Jeppe_Bayer/X_chromosome_social_spiders/X_Evolution/reference_genomes/S_dumicola/DUM_hifi_hic_scaffolded_trim.fa \
-d people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_dumicola

bash people/Jeppe_Bayer/scripts/02_data_preparation/02_data_preparation_init.sh \
-s people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_mimosarum/Antananarivo \
-r people/Jeppe_Bayer/X_chromosome_social_spiders/X_Evolution/reference_genomes/S_mimosarum/MIM_hifi_hic_scaffolded_trim.fa \
-d people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_mimosarum

bash people/Jeppe_Bayer/scripts/02_data_preparation/02_data_preparation_init.sh \
-s people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_mimosarum/Pon \
-r people/Jeppe_Bayer/X_chromosome_social_spiders/X_Evolution/reference_genomes/S_mimosarum/MIM_hifi_hic_scaffolded_trim.fa \
-d people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_mimosarum

bash people/Jeppe_Bayer/scripts/02_data_preparation/02_data_preparation_init.sh \
-s people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_mimosarum/Wee \
-r people/Jeppe_Bayer/X_chromosome_social_spiders/X_Evolution/reference_genomes/S_mimosarum/MIM_hifi_hic_scaffolded_trim.fa \
-d people/Jeppe_Bayer/X_chromosome_social_spiders/population_data/S_mimosarum

exit 0