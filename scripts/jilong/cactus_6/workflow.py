from gwf import *
from workflow_templates import *
from workflow_targets import *

gwf = Workflow()
ste_dict = build_sp_chrom_dict()
code_list = ["DUM","TENT","SARA","BI","MIM","LIN"]
prepare_cactus_fasta(gwf,ste_dict,code_list)

chrom_map = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/cactus/ste_chr_map.tsv"
tree = "(LIN,(((SARA,BI),(DUM,TENT)),MIM))Ste;"
run_cactus_setting(tree,chrom_map,gwf)

run_cactus(ste_dict,gwf)
