from gwf import *
from workflow_templates import *
def prepare_cactus_fasta(gwf,ste_dic,code_list):
    for ste_chr in ste_dic:
        for sp in code_list:
            seqs_name = ste_dic[ste_chr][sp]
            seqs_list = []
            for seq in seqs_name:
                if sp in ["DUM","TENT","LIN"]:
                    seq_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/{sp}/RMdatabase/combine/scaffolds/masked_{seq}.fa".format(sp=sp,seq=seq)
                    seqs_list.append(seq_path)
                if sp in ["SARA","BI","MIM"]:
                    number = seq.split("_")[-1]
                    seq_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/split/{sp}_hifi_hic_scaffolded_trim.fa.masked.{number}".format(sp=sp,number=number)
                    seqs_list.append(seq_path)
            path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/cactus_6/{ste_chr}".format(ste_chr=ste_chr)
            out = path+"/"+"{ste_chr}_{sp}_soft_masked.fa".format(ste_chr=ste_chr,sp=sp)
            gwf.target_from_template(
                name="merge_fa_{ste}_{sp}".format(ste=ste_chr,sp=sp),
                template = merge_fa(
                    fa_list=seqs_list,
                    out=out,
                    path=path
                    )
                )
def run_cactus_setting(tree,chrom_map,gwf):
    gwf.target_from_template(
                name="cactus_prepare_all",
                template = cactus_setting(
                    tree = tree,
                    chrom_map=chrom_map
                    )
                )
    return LOG_PATH+"/cactus_6_prepare.DONE"

def run_cactus(ste_dic,gwf):
    for ste_chr in ste_dic:
        jobstore="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/cactus_6/{ste_chr}/{ste_chr}_cactus".format(ste_chr=ste_chr)
        seqFile="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/cactus_6/{ste_chr}/{ste_chr}_cactus.txt".format(ste_chr=ste_chr)
        outputHal="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/cactus_6/{ste_chr}/{ste_chr}.hal".format(ste_chr=ste_chr)
        gwf.target_from_template(
            name="cactus_{ste_chr}".format(ste_chr=ste_chr),
            template = cactus(
              jobstore=jobstore,
              seqFile=seqFile,
              outputHal=outputHal,
              chrom_group=ste_chr
              )
        )
