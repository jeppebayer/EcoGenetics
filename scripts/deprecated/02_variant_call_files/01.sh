#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 6G
#SBATCH --cpus-per-task 8
#SBATCH --time 20:00:00

vcftools \
--bcf people/Jeppe_Bayer/steps/02_variant_call_files/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.calls.bcf \
--freq \
--temp people/Jeppe_Bayer/steps/temp/ \
--out people/Jeppe_Bayer/steps/02_variant_call_files/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642 ; \

vcftools \
--bcf people/Jeppe_Bayer/steps/02_variant_call_files/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642.calls.bcf \
--counts \
--temp people/Jeppe_Bayer/steps/temp/ \
--out people/Jeppe_Bayer/steps/02_variant_call_files/Orchesella_cincta/NYS-F/E100050923_L01_ANIcnqkR228563-642
