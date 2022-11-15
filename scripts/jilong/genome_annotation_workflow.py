from gwf import *
from workflow_dicts import *
from workflow_templates import *
from workflow_targets import *
gwf = Workflow()
## STAR index
for sp in ["SARA","MIM","BI"]:
    build_STAR_index_sp(sp,gwf)
## STAR align
logs_dict ={}
for sp in ["SARA","MIM","BI"]:
    align_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_align".format(sp=sp)
    index_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_index".format(sp=sp)
    logs = STAR_align_dict(sp_fq_dict[sp],sp,align_path,index_path,gwf)
    gwf.target_from_template(
        name = "STAR_align_{sp}".format(sp=sp),
        template = logs_sum(logs,LOG_PATH+"/align_STAR_{sp}.DONE".format(sp=sp))
        )
    logs_dict[sp]=logs
for sp in ["SARA","MIM","BI"]:
    path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_merged".format(sp=sp)
    filepath = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_align".format(sp=sp)
    filenames = sp_fq_dict[sp].keys() ##change points if need to swtich species
    log = merge_STAR_bam(path,sp,filepath,filenames,logs_dict[sp],gwf)
    path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_merged".format(sp=sp)
    bam_in = "STAR_merged_{sp}".format(sp=sp)
    bam_out = "STAR_merged_{sp}_sort".format(sp=sp)
    log = sort_STAR_bam(path,bam_in,bam_out,sp,gwf) 

## repeat_mask genome
for sp in ["SARA","MIM","BI"]:
    genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/{sp}/{sp}_hifi_hic_scaffolded_trim.fa".format(sp=sp)
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/model".format(sp=sp)
    out = sp
    model_sp(genome,sp,path,out,gwf)
    genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/{sp}/{sp}_hifi_hic_scaffolded_trim.fa".format(sp=sp)
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/mask".format(sp=sp)
    lib1 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/repbase/repbase_arthropoda.fa"
    lib2 = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/model/{sp}-families.fa".format(sp=sp)
    lib = "{sp}_repbase.fa".format(sp=sp)
    repeat_sp(genome,lib1,lib2,lib,path,sp,gwf)

## split masked genome by chrom
sp = "SARA"
numbers = 19
logs = grep_chrom_fa(sp,numbers,gwf)
gwf.target_from_template(
    name = "split_masked_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/Split_masked_chrom_{sp}.DONE".format(sp=sp))
    )
sp = "MIM"
numbers = 19
logs = grep_chrom_fa(sp,numbers,gwf)
gwf.target_from_template(
    name = "split_masked_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/Split_masked_chrom_{sp}.DONE".format(sp=sp))
    )
sp = "BI"
numbers = 18
logs = grep_chrom_fa(sp,numbers,gwf)
gwf.target_from_template(
    name = "split_masked_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/Split_masked_chrom_{sp}.DONE".format(sp=sp))
    )
## braker split masked genome
sp = "SARA"
numbers = 17
logs = braker_sp_chrom(sp,numbers,gwf)
gwf.target_from_template(
    name = "braker_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/BRAKER_{sp}.DONE".format(sp=sp))
    )
sp = "MIM"
numbers = 19
logs = braker_sp_chrom(sp,numbers,gwf)
gwf.target_from_template(
    name = "braker_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/BRAKER_{sp}.DONE".format(sp=sp))
    )
sp = "BI"
numbers = 18
logs = braker_sp_chrom(sp,numbers,gwf)
gwf.target_from_template(
    name = "braker_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/BRAKER_{sp}.DONE".format(sp=sp))
    )
## gff2fa
sp = "SARA"
numbers = chrom_id_dict[sp]
logs = agat_gff2fa(sp,numbers,gwf)
gwf.target_from_template(
    name = "braker_gff2fa_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/gff2fa_{sp}.DONE".format(sp=sp))
    )
## cat all results from braker
sp = "SARA"
numbers = chrom_id_dict[sp]
combine_results_braker_fa(sp,numbers,gwf)
## BUSCO check annotation completeness
sp = "SARA"
busco_sp(sp,gwf)
