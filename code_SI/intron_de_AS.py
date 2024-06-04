import argparse
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import pickle
import random as rd
import numpy as np
from scipy import stats
from scipy.stats import permutation_test

from scipy.stats import mannwhitneyu

parse = argparse.ArgumentParser()
parse.add_argument("--Normalized")
args = parse.parse_args()


def statistic(x, y):
    return np.mean(x) - np.mean(y)

with open("dmel-all-r6.48.dico.pkl", "rb") as f_o:
    dico_gtf = pickle.load(f_o)
    


intron_size_file = "gene_introns_exons_span.tsv"
dico_intron = {}
with open(intron_size_file) as f_in:
    f_in.readline()
    for l in f_in:
        spt = l.strip().split()
        dico_intron[dico_gtf[spt[1]]["symbol"]] = { "size": int(spt[-2]),
                                                    "prop": float(spt[-1]) }

def get_intron_all(transcript, strand):
    """return all intron"""
    intron_ = set()
    for tr_id, tr_dico in transcript.items():
        exons = tr_dico.get("exon", [])
        utr5 = tr_dico.get("5UTR", [])
        utr3 = tr_dico.get("3UTR", [])
        N = utr5 + exons + utr3
        if len(N) == 1:
            continue
        if strand == "+":
            N = sorted(N, key=lambda x: x["start"])
            for i,e in enumerate(N[:-1]):
                x = N[i+1]["start"] - e["end"]
                if x > 0:
                    intron_.add((N[i+1]["start"], e["end"]))
        else:
            N = sorted(N, key=lambda x: x["end"], reverse=True)
            for i,e in enumerate(N[:-1]):
                x = e["start"] - N[i+1]["end"] 
                if x > 0:
                    intron_.add((e["start"], N[i+1]["end"] ))

        
    return list(intron_)
        

def get_intron(transcript, strand):
    """return largest intron"""
    largest = 0
    for tr_id, tr_dico in transcript.items():
        exons = tr_dico.get("exon", [])
        utr5 = tr_dico.get("5UTR", [])
        utr3 = tr_dico.get("3UTR", [])
        N = utr5 + exons + utr3
        if len(N) == 1:
            continue
        if strand == "+":
            N = sorted(N, key=lambda x: x["start"])
            for i,e in enumerate(N[:-1]):
                x = N[i+1]["start"] - e["end"]
                #print(x)
                if x > largest:
                    largest = x
        else:
            N = sorted(N, key=lambda x: x["end"], reverse=True)
            for i,e in enumerate(N[:-1]):
                x = e["start"] - N[i+1]["end"] 
                if x > largest:
                    largest = x

        
    return largest


dico = {}
intron = []



for gene_id, dico_gene in dico_gtf.items():
    if dico_gene["biotype"] not in  ["protein_coding_gene"]:
        continue
        
    strand = dico_gene["strand"]
    l = get_intron(dico_gene["transcript"], strand)
    intron.extend(get_intron_all(dico_gene["transcript"], strand))

    dico[gene_id] = l # largest intron all transcript considered for a genes
    
df_AS_DE = pd.read_csv("AS_DE_table.csv", index_col=0)

srpk_as = set(df_AS_DE[df_AS_DE["AS_qval_SRPK"]==1]["FlyBase ID"].values)
u2af38_as = set(df_AS_DE[df_AS_DE["AS_qval_U2af38"]==1]["FlyBase ID"].values)

    
ncount = pd.read_csv(args.Normalized)#"Normalized.count.csv")
fb_expr = ncount[ncount["Control_rep1.norm 	Control_rep2.norm 	Control_rep3.norm".split()].sum(axis=1)  > 30]["FlyBaseID"].values
    
    
df = pd.read_csv("AS_DE_table.csv", index_col=0)
liste_genes = ['kl-2', 'kl-3', 'kl-5', 'CCY', 'ORY', 'PRY', 'WDY', 'Ppr-Y']


gene_with_gap = """FBgn0267428 FBgn0267431 
FBgn0267430 FBgn0267433 FBgn0267489 FBgn0267432 FBgn0001313 
FBgn0046697 FBgn0267449 FBgn0046323 FBgn0267592""".split()

genes_nogap_avg = []
long_intron_gene_name = []
for k, v in dico.items():
    if k in gene_with_gap:
        continue
    if v <= 100:
        genes_nogap_avg.append(k)
    if 50_000 <= v <= 100_000:
        long_intron_gene_name.append(k)

Gene_gap_tokeep = ["FBgn0267431", "FBgn0267430", "FBgn0263112"]
#U2af38
genotype = "U2af38"


f2 = open("S2.BB.txt", "w")

filter_ = set(fb_expr).intersection(u2af38_as)
val_avg_gene = df[(df["FlyBase ID"].isin(genes_nogap_avg)) & (df["FlyBase ID"].isin(filter_))]["log2(U2af38/Control)"]
val_long_gene = df[(df["FlyBase ID"].isin(long_intron_gene_name)) & (df["FlyBase ID"].isin(filter_))]["log2(U2af38/Control)"]
# large_Y = df[df["Gene"].isin(liste_genes)]["log2(U2af38/Control)"]
# pzl = df[df["FlyBase ID"].isin(Gene_gap_tokeep)]["log2(U2af38/Control)"]


plt.scatter([0 + rd.gauss(0, 0.01) for x in range(len(val_avg_gene))], val_avg_gene, color="black", s=1)
plt.scatter([1 + rd.gauss(0, 0.01) for x in range(len(val_long_gene))], val_long_gene, color="black", s=1)
# plt.scatter([2 + rd.gauss(0, 0.02) for x in range(len(pzl))], pzl, color="black", s=1)
# plt.scatter([3 + rd.gauss(0, 0.02) for x in range(len(large_Y))], large_Y, color="black", s=1)
plt.xticks([0,1], labels=["small_intron_gene","large non gap gene"], rotation=45, ha="right")
plt.gca().axhline(0, color="black")
plt.title("U2af38")
plt.savefig("final_figure_4_2_B_U2af38.pdf")
plt.savefig("final_figure_4_2_B_U2af38.png")
plt.show()

f2.write("{}\t{}\t{}\n".format(genotype, "<100pb", "\t".join(list(map(str, val_avg_gene)))))
f2.write("{}\t{}\t{}\n".format(genotype, "50-110kb", "\t".join(list(map(str, val_long_gene)))))
# f2.write("{}\t{}\t{}\n".format(genotype, "autosomal gigantic", "\t".join(list(map(str, pzl)))))
# f2.write("{}\t{}\t{}\n".format(genotype, "Y-linked", "\t".join(list(map(str, large_Y)))))


# res = permutation_test((large_Y, val_avg_gene), statistic, 
#                        n_resamples=1000)#, alternative='less')
# print("U2af38, Y vs small, res: ", res.statistic, res.pvalue)

# res = permutation_test((large_Y, val_long_gene), statistic, 
#                        n_resamples=1000)#, alternative='less')
# print("U2af38, Y vs long, res: ", res.statistic, res.pvalue)

res = permutation_test((val_avg_gene, val_long_gene), statistic, 
                       n_resamples=1000)#, alternative='less')
print("U2af38, small vs long, res: ", res.statistic, res.pvalue)


# U1, p  = mannwhitneyu(large_Y, val_avg_gene)
# print(" MannU U2af38, Y vs small, res: ", U1, p)

# U1, p  = mannwhitneyu(large_Y, val_long_gene)
# print(" MannU U2af38, Y vs long, res: ", U1, p)

U1, p  = mannwhitneyu(val_avg_gene, val_long_gene)
print("MannU U2af38, small vs long, res: ", U1, p)





# SRPK
genotype = "SRPK"

filter_ = set(fb_expr).intersection(srpk_as)

val_avg_gene = df[(df["FlyBase ID"].isin(genes_nogap_avg)) & (df["FlyBase ID"].isin(filter_))]["log2(SRPK/Control)"]
val_long_gene = df[(df["FlyBase ID"].isin(long_intron_gene_name)) & (df["FlyBase ID"].isin(filter_))]["log2(SRPK/Control)"]
# large_Y = df[df["Gene"].isin(liste_genes)]["log2(SRPK/Control)"]
# pzl = df[df["FlyBase ID"].isin(Gene_gap_tokeep)]["log2(SRPK/Control)"]

f2.write("{}\t{}\t{}\n".format(genotype, "<100pb", "\t".join(list(map(str, val_avg_gene)))))
f2.write("{}\t{}\t{}\n".format(genotype, "50-110kb", "\t".join(list(map(str, val_long_gene)))))
# f2.write("{}\t{}\t{}\n".format(genotype, "autosomal gigantic", "\t".join(list(map(str, pzl)))))
# f2.write("{}\t{}\t{}\n".format(genotype, "Y-linked", "\t".join(list(map(str, large_Y)))))

plt.scatter([0 + rd.gauss(0, 0.01) for x in range(len(val_avg_gene))], val_avg_gene, color="black", s=1)
plt.scatter([1 + rd.gauss(0, 0.01) for x in range(len(val_long_gene))], val_long_gene, color="black", s=1)
# plt.scatter([2 + rd.gauss(0, 0.02) for x in range(len(pzl))], pzl, color="black", s=1)
# plt.scatter([3 + rd.gauss(0, 0.02) for x in range(len(large_Y))], large_Y, color="black", s=1)
plt.xticks([0,1], labels=["small_intron_gene","large non gap gene"], rotation=45, ha="right")



f2.close()


# res = permutation_test((large_Y, val_avg_gene), statistic,
#                        n_resamples=1000)#, alternative='less')
# print("SRPK, Y vs small, res: ", res.statistic, res.pvalue)

# res = permutation_test((large_Y, val_long_gene), statistic,
#                        n_resamples=1000)#, alternative='less')
# print("SRPK, Y vs long, res: ", res.statistic, res.pvalue)

res = permutation_test((val_avg_gene, val_long_gene), statistic, 
                       n_resamples=1000)#, alternative='less')
print("SRPK, small vs long, res: ", res.statistic, res.pvalue)



# U1, p  = mannwhitneyu(large_Y, val_avg_gene)
# print(" MannU SRPK, Y vs small, res: ", U1, p)

# U1, p  = mannwhitneyu(large_Y, val_long_gene)
# print(" MannU SRPK, Y vs long, res: ", U1, p)

U1, p  = mannwhitneyu(val_avg_gene, val_long_gene)
print("MannU SRPK, small vs long, res: ", U1, p)



plt.gca().axhline(0, color="black")
plt.savefig("final_figure_4_2_B_SRPK.pdf")
plt.savefig("final_figure_4_2_B_SRPK.png")

plt.title("SRPK")
plt.show()



