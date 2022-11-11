from gwf import *
def bwa_index(fasta,species):
    inputs = [fasta]
    outputs = ["{species}_bwa_index.log".format(species=species)]
    options = {
              'cores': 1,
              'memory': '16g',
              'walltime':'12:00:00',
              'account':"spider2"
    }
    spec ="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    bwa index {fasta} > {log}.tmp 
    mv {log}.tmp {log}
    """.format(fasta=fasta,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def juicer_hic_align(species,fasta):
    inputs = []
    outputs = ["{}_juicer.log".format(species)]
    options = {
              'cores': 1,
              'memory': '1g',
              'walltime':'1:00:00',
              'account':"spider2"
    }
    spec = """
    /home/jilong/software/juicer/scripts/juicer.sh -d /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/{species} -D /home/jilong/software/juicer -p /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/{species}/chrom.sizes -s none -z {fasta} -q short -Q 12:00:00 -l normal -L 24:00:00 -t 36 > {species}_juicer.log
    """.format(species=species,fasta=fasta)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def assemble_3ddna(folder,species,fasta,merged,r):
    inputs = []
    outputs = ["{folder}/{species}_3ddna.log".format(folder=folder,species=species)]
    options = {
              'cores': 16,
              'memory': '200g',
              'walltime':'72:00:00',
              'account':"spider2"
    }
    spec = """
    cd {folder}
    /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/3d_dna/3d-dna/run-asm-pipeline.sh -r {r} --editor-repeat-coverage 30 --editor-coarse-stringency 20 {fasta} {merged} > {folder}/{species}_3ddna.log.tmp
    mv {folder}/{species}_3ddna.log.tmp {folder}/{species}_3ddna.log
    """.format(folder=folder,fasta=fasta,merged=merged,r=r,species=species)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def assemble_3ddna_fast(folder,species,fasta,merged,r):
    inputs = []
    outputs = ["{folder}/{species}_3ddna.log".format(folder=folder,species=species)]
    options = {
              'cores': 16,
              'memory': '80g',
              'walltime':'72:00:00',
              'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo jobinfo $SLURM_JOBID
    echo "start 3d dna"
    date
    mkdir -p {folder}
    cd {folder}
    export _JAVA_OPTIONS=-Djava.io.tmpdir=/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/tmp
    /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/3d_dna/3d-dna/run-asm-pipeline.sh -r {r} --early-exit --editor-repeat-coverage 30 --editor-coarse-stringency 20 --splitter-input-size 500000 --splitter-coarse-resolution 500000 --splitter-coarse-stringency 20 {fasta} {merged} > {folder}/{species}_3ddna.log.tmp
    mv {folder}/{species}_3ddna.log.tmp {folder}/{species}_3ddna.log
    date
    """.format(folder=folder,fasta=fasta,merged=merged,r=r,species=species)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def assemble_3ddna_norm(folder,species,fasta,merged,r):
    inputs = []
    outputs = ["{folder}/{species}_3ddna.log".format(folder=folder,species=species)]
    options = {
              'cores': 16,
              'memory': '200g',
              'walltime':'72:00:00',
              'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo jobinfo $SLURM_JOBID
    echo "start 3d dna"
    date
    cd {folder}
    /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/3d_dna/3d-dna/run-asm-pipeline.sh -r {r} --editor-repeat-coverage 30 --editor-coarse-stringency 20 -s split -f --splitter-input-size 500000 --splitter-coarse-resolution 500000 --splitter-coarse-stringency 20 {fasta} {merged} > {folder}/{species}_3ddna.log.tmp
    mv {folder}/{species}_3ddna.log.tmp {folder}/{species}_3ddna.log
    date
    """.format(folder=folder,fasta=fasta,merged=merged,r=r,species=species)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def assemble_3ddna_stage(folder,species,fasta,merged,r,stage):
    inputs = []
    outputs = ["{folder}/{species}_3ddna_{stage}.log".format(stage=stage,folder=folder,species=species)]
    options = {
              'cores': 24,
              'memory': '100g',
              'walltime':'72:00:00',
              'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo jobinfo $SLURM_JOBID
    echo "start 3d dna"
    date
    mkdir -p {folder}
    cd {folder}
    export _JAVA_OPTIONS=-Djava.io.tmpdir=/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/tmp
    /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/3d_dna/3d-dna/run-asm-pipeline.sh -r {r} --editor-repeat-coverage 30 --editor-coarse-stringency 20 --splitter-input-size 500000 --splitter-coarse-resolution 500000 --splitter-coarse-stringency 20 --fast-start --stage {stage} {fasta} {merged} > {folder}/{species}_3ddna.log.tmp
    mv {folder}/{species}_3ddna_{stage}.log.tmp {folder}/{species}_3ddna_{stage}.log
    date
    """.format(folder=folder,fasta=fasta,merged=merged,r=r,species=species,stage=stage)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def export_fasta_review(folder,species,review,draft,merged):
    inputs = [review,draft,merged]
    outputs = ["{folder}/{species}_export.log".format(folder=folder,species=species)]
    options = {
              'cores': 4,
              'memory': '64g',
              'walltime':'72:00:00',
              'account':"spider2"
    }
    spec = """
    echo jobinfo $SLURM_JOBID
    cd {folder}
    export _JAVA_OPTIONS=-Djava.io.tmpdir=/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/tmp
    /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/3d_dna/3d-dna/run-asm-pipeline-post-review.sh --build-gapped-map --sort-output -s finalize -r {review} {draft} {merged} > {folder}/{species}_export.log.tmp
    mv {folder}/{species}_export.log.tmp {folder}/{species}_export.log
    """.format(folder=folder,review=review,draft=draft,merged=merged,species=species)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)



gwf = Workflow()
## Indexing
species = "DUM"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/reference/DUM_hifi.tmp.hic.hap2.p_ctg.fa"
gwf.target_from_template(
        name="bwa_index_{}".format(species),
        template = bwa_index(fasta,species)
    )
species = "TENT"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/TENT/reference/TENT_hifi.tmp.hic.hap1.p_ctg.fa"
gwf.target_from_template(
        name="bwa_index_{}".format(species),
        template = bwa_index(fasta,species)
    )
species = "LIN"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/LIN/reference/LIN_hifi.tmp.hic.hap2.p_ctg.fa"
gwf.target_from_template(
        name="bwa_index_{}".format(species),
        template = bwa_index(fasta,species)
    )
species = "MIM"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/MIM/reference/MIM_hifi.tmp.hic.hap2.p_ctg.fa"
gwf.target_from_template(
        name="bwa_index_{}".format(species),
        template = bwa_index(fasta,species)
    )
species = "BI"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/BI/reference/BI_hifi.tmp.hic.hap1.p_ctg.fa"
gwf.target_from_template(
        name="bwa_index_{}".format(species),
        template = bwa_index(fasta,species)
    )

## Juicer aligning
species = "SARA"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/SARA/reference/SARA_hifi_hic.asm.tmp.hic.hap1.p_ctg.fa"
gwf.target_from_template(
        name="juicer_{}".format(species),
        template = juicer_hic_align(species,fasta)
    )
species = "DUM"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/reference/DUM_hifi.tmp.hic.hap2.p_ctg.fa"
gwf.target_from_template(
        name="juicer_{}".format(species),
        template = juicer_hic_align(species,fasta)
    )
species = "TENT"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/TENT/reference/TENT_hifi.tmp.hic.hap1.p_ctg.fa"
gwf.target_from_template(
        name="juicer_{}".format(species),
        template = juicer_hic_align(species,fasta)
    )
species = "LIN"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/LIN/reference/LIN_hifi.tmp.hic.hap2.p_ctg.fa"
gwf.target_from_template(
        name="juicer_{}".format(species),
        template = juicer_hic_align(species,fasta)
    )

species = "MIM"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/MIM/reference/MIM_hifi.tmp.hic.hap2.p_ctg.fa"
gwf.target_from_template(
        name="juicer_{}".format(species),
        template = juicer_hic_align(species,fasta)
    )

species = "BI"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/BI/reference/BI_hifi.tmp.hic.hap1.p_ctg.fa"
gwf.target_from_template(
        name="juicer_{}".format(species),
        template = juicer_hic_align(species,fasta)
    )
##assemble 3d dna
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/SARA/3d_dna/SARA_r0"
species = "SARA"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/SARA/reference/SARA_hifi_hic.asm.tmp.hic.hap1.p_ctg.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/SARA/aligned/merged_nodups.txt"
r = 0
gwf.target_from_template(
        name="assemble_r{r}_{species}".format(r=r,species=species),
        template = assemble_3ddna(folder,species,fasta,merged,r)
    )
# DUM
#folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/3d_dna"
#species = "DUM"
#fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/reference/DUM_hifi.tmp.hic.hap2.p_ctg.fa"
#merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/aligned/merged_nodups.txt"
#r = 0
#gwf.target_from_template(
#        name="assemble_r{r}_{species}".format(r=r,species=species),
#        template = assemble_3ddna_fast(folder,species,fasta,merged,r)
#    )
# TENT
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/TENT/3d_dna"
species = "TENT"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/TENT/reference/TENT_hifi.tmp.hic.hap1.p_ctg.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/TENT/aligned/merged_nodups.txt"
r = 0
gwf.target_from_template(
        name="assemble_r{r}_{species}".format(r=r,species=species),
        template = assemble_3ddna_fast(folder,species,fasta,merged,r)
    )
# LIN
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/LIN/3d_dna"
species = "LIN"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/LIN/reference/LIN_hifi.tmp.hic.hap2.p_ctg.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/LIN/aligned/merged_nodups.txt"
r = 0
gwf.target_from_template(
        name="assemble_r{r}_{species}".format(r=r,species=species),
        template = assemble_3ddna_fast(folder,species,fasta,merged,r)
    )
# MIM
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/MIM/3d_dna"
species = "MIM"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/MIM/reference/MIM_hifi.tmp.hic.hap2.p_ctg.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/MIM/aligned/merged_nodups.txt"
r = 0
gwf.target_from_template(
        name="assemble_r{r}_{species}".format(r=r,species=species),
        template = assemble_3ddna_fast(folder,species,fasta,merged,r)
    )


# BI
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/BI/3d_dna"
species = "BI"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/BI/reference/BI_hifi.tmp.hic.hap1.p_ctg.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/BI/aligned/merged_nodups.txt"
r = 0
gwf.target_from_template(
        name="assemble_r{r}_{species}".format(r=r,species=species),
        template = assemble_3ddna_fast(folder,species,fasta,merged,r)
    )
#folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/SARA/3d_dna/SARA_r1"
#r = 1
#gwf.target_from_template(
#        name="assemble_r{r}_{species}".format(r=r,species=species),
#        template = assemble_3ddna(folder,species,fasta,merged,r)
#    )

review = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/SARA/3d_dna/SARA_r0/SARA_hifi_hic.asm.tmp.hic.hap1.p_ctg.0.review.assembly"
draft = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/SARA/reference/SARA_hifi_hic.asm.tmp.hic.hap1.p_ctg.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/SARA/aligned/merged_nodups.txt"
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/SARA"
species = "SARA"
gwf.target_from_template(
        name="export_review_{species}".format(species=species),
        template = export_fasta_review(folder,species,review,draft,merged)
    )


review = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/3d_dna/DUM_hifi.tmp.hic.hap2.p_ctg.0.review.assembly"
draft = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/reference/DUM_hifi.tmp.hic.hap2.p_ctg_wrapped.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/aligned/merged_nodups.txt"
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/DUM"
species = "DUM"
gwf.target_from_template(
        name="export_review_{species}".format(species=species),
        template = export_fasta_review(folder,species,review,draft,merged)
    )

review = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/TENT/3d_dna/TENT_hifi.tmp.hic.hap1.p_ctg.0.review.assembly"
draft = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/TENT/reference/TENT_hifi.tmp.hic.hap1.p_ctg_wrapped.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/TENT/aligned/merged_nodups.txt"
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/TENT"
species = "TENT"
gwf.target_from_template(
        name="export_review_{species}".format(species=species),
        template = export_fasta_review(folder,species,review,draft,merged)
    )

review = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/LIN/3d_dna/LIN_hifi.tmp.hic.hap2.p_ctg.0.review.assembly"
draft = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/LIN/reference/LIN_hifi.tmp.hic.hap2.p_ctg_wrapped.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/LIN/aligned/merged_nodups.txt"
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/LIN"
species = "LIN"
gwf.target_from_template(
        name="export_review_{species}".format(species=species),
        template = export_fasta_review(folder,species,review,draft,merged)
    )

review = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/MIM/3d_dna/MIM_hifi.tmp.hic.hap2.p_ctg.0.review.assembly"
draft = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/MIM/reference/MIM_hifi.tmp.hic.hap2.p_ctg_wrapped.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/MIM/aligned/merged_nodups.txt"
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/MIM"
species = "MIM"
gwf.target_from_template(
        name="export_review_{species}".format(species=species),
        template = export_fasta_review(folder,species,review,draft,merged)
    )

review = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/BI/3d_dna/BI_hifi.tmp.hic.hap1.p_ctg.0.review.assembly"
draft = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/BI/reference/BI_hifi.tmp.hic.hap1.p_ctg_wrapped.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/BI/aligned/merged_nodups.txt"
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/BI"
species = "BI"
gwf.target_from_template(
        name="export_review_{species}".format(species=species),
        template = export_fasta_review(folder,species,review,draft,merged)
    )
