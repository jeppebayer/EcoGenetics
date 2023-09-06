import sys
#tree = sys.argv[1]
#chrom_map = sys.argv[2]
#ste_map = open(chrom_map)
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

for ste_chr in ste_dic:
    tmp_file = open("/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/cactus_6/{ste_chr}/{ste_chr}_cactus.txt".format(ste_chr=ste_chr),"w")
    tmp_file.write("(LIN,(((SARA,BI),(DUM,TENT)),MIM))Ste;\n")
#    tmp_file.write("{tree}\n".format(tree=tree))
    tmp_file.write("\n")
    for sp in code_list:
        tmp_file.write(sp+" "+"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/cactus_6/{ste_chr}/{ste_chr}_{sp}_soft_masked.fa".format(ste_chr=ste_chr,sp=sp)+"\n")
