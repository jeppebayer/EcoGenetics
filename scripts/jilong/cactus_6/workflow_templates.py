from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/cactus_6/logs"
def logs_sum(logs,sum_log):
    inputs = logs
    outputs = [sum_log]
    options = {
              'cores': 1,
              'memory': '1g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """
    echo "finished" > {sum_log}
    echo date
""".format(sum_log=sum_log)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def build_sp_chrom_dict():
    ste_map = open("/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/cactus/ste_chr_map.tsv")
    ste_dic = {}
    code_list = ["DUM","TENT","SARA","LIN","MIM","BI"]
    for ste_chr in ste_map:
        ste_info = ste_chr.strip("\n").split("\t")
        ste_chr_name = ste_info[0]
        sp_dic = {}
        for sp in code_list:
            sp_dic[sp] = []
        for sp_chr in ste_info[1:]:
            sp_chrs = sp_chr.split(";")
            for sp_scaffold in sp_chrs:
                sp = sp_scaffold.split(".")[0]
                scaffold = sp_scaffold.split(".")[1]
                sp_dic[sp].append(scaffold)
        ste_dic[ste_chr_name]=sp_dic
    return ste_dic
def merge_fa(fa_list,out,path):
    inputs = fa_list
    fa_string=""
    for fa in fa_list:
        fa_string = fa_string+fa+" "
    outputs = [out]
    options = {
        'cores':1,
        'memory':'1g',
        'walltime':'1:00:00',
        'account':'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cat {string} > {out}.tmp
    mv {out}.tmp {out}    
    """.format(path=path,string = fa_string, out = out)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def cactus_setting(tree,chrom_map):
    inputs = [chrom_map]
    outputs = [LOG_PATH+"/cactus_6_prepare.DONE"]
    options = {
        'cores':1,
        'memory':'1g',
        'walltime':'2:00:00',
        'account':'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo jobinfo $SLURM_JOBID
    date
    python cactus_setting.py > {log}.tmp
    mv {log}.tmp {log}    
    """.format(log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def cactus(jobstore,seqFile,outputHal,chrom_group):
    inputs=[LOG_PATH+"/cactus_6_prepare.DONE"]
    outputs=[LOG_PATH+"/cactus_{chrom_group}.DONE".format(chrom_group=chrom_group)]
    options = {
        'cores':32,
        'memory':'128g',
        'walltime':'24:00:00',
        'account':'spider2'
    }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate cactus_align
    echo "start cactus test"
    echo jobinfo $SLURM_JOBID
    date
    cd /home/jilong/software/cactus-bin-v2.0.5
    source cactus_env/bin/activate
    cactus --consCores 32 {jobstore} {seqFile} {outputHal}
    echo "finished cactus"
    date
    jobinfo $SLURM_JOBID
    """.format(jobstore=jobstore,seqFile=seqFile,outputHal=outputHal)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
