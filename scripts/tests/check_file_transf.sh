#!/bin/bash
dir=(/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Isotoma_sp/IsoSp_CA_FSJ-C221
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_villosa/OrcVil_CA_KNJ-C208
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Isotoma_sp/IsoSp_CA_KNJ-C209
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nicoleti/EntNic_CA_ST_6J-C84
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nicoleti/EntNic_DSJ-C190
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nicoleti/EntNic_FUR-C230
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Lepidocyrtus_lignorum/LepLig_FUR-C231
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nicoleti/EntNic_HHJ-C162
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Lepidocyrtus_lignorum/LepLig_HHJ-C215
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nicoleti/EntNic_JHJ-C220
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Isotoma_sp/IsoSp_JHJ-C224
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_cincta/OrcCin_JHJ-C222
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Isotoma_sp/IsoSp_K_FSJ-C210
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nicoleti/EntNic_KØJ-C212
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Isotoma_sp/IsoSp_KØJ-C213
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Isotoma_sp/IsoSp_K_KNJ-C217
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nivalis/EntNiv_MSJ-C93
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Orchesella_villosa/OrcVil_MSJ-C3
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nicoleti/EntNic_ÅRJ-C225
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Pogonognathellus_flavescens/PogFla_ÅRJ-C228
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nicoleti/EntNic_SKJ-C16
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Lepidocyrtus_lignorum/LepLig_SKJ-C102
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Isotoma_sp/IsoSp_SKJ-C103
/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nicoleti/EntNic_SKJ-C216)

for i in "${dir[@]}"; do
    ls "$i" | wc -l
done

exit 0