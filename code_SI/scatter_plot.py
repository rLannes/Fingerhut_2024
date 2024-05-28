import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import re
import sklearn
import numpy as np
from scipy import stats



df = pd.read_csv("gene_introns_exons_span.tsv", sep="\t")
df.columns = ["chr", "ID", "symbol", "tr", "span", "exon", "score"]

df_AS_DE = pd.read_csv("AS_DE_table.csv", index_col=0)

def get_bins(bins, new_score):
    bin_ = bins[::-1][1:]
    k = []
    for b in bin_:

        try:
            k.append(len([x for x in new_score if x >= b]))
            new_score = [x for x in new_score if x < b]
        except:
            print(new_score)
            raise

    return k[::-1]


dico = {}
gene = df_AS_DE[df_AS_DE["AS_qval_SRPK"] == 1 ]["Gene"].values
dico["SRPK_vs_control"] = df[df["symbol"].isin(gene)][["score", "symbol"]].values

gene = df_AS_DE[df_AS_DE["AS_qval_U2af38"] == 1 ]["Gene"].values
dico["U2af38_vs_control"] = df[df["symbol"].isin(gene)][["score", "symbol"]].values


fig, ax = plt.subplots(1,1, figsize=(8, 8))
(n, bins, patches) = plt.hist(df["score"].values, bins=10, color="#457b9d", alpha=0.7, edgecolor='black')
plt.close()
fig, ax = plt.subplots(1,1, figsize=(8, 8))

for feature, score_labels in dico.items():
    if feature not in ["SRPK_vs_control", "U2af38_vs_control"]:
        continue
    y = []
    seen = set()
    score = score_labels
    for r in score:
        # print(r)
        # for x in r:
        try:
            if r[1] not in seen:
                seen.add(r[1])
                y.append(r[0])
        except:
            print(r)
            raise

    x = range(len(y))
    bin_middle = [bins[i] + ((bins[i + 1] - bins[i]) / 2) for i in range(0, len(bins) - 1, )]
    n_s = get_bins(bins, y)
    plt.scatter(bin_middle, np.array(n_s) / np.array(n), label=feature.split("_")[0])
    print(x, n_s, n)
    print(feature, "spearmanr")
    print( bin_middle, np.array(n_s) / np.array(n))
    res = stats.spearmanr( bin_middle, np.array(n_s) / np.array(n))
    print(res)
    print()
    
plt.gca().set_ylabel('Proportion of genes with AS by each bin')

plt.gca().set_xticks(bins)
plt.gca().set_xticklabels([ "{}%".format(round(x * 100)) for x in list(bins)])#, rotation=45, ha="right")


plt.title("Increase of AS event in gene with high intron content")
plt.legend(bbox_to_anchor=(1.15, 1), loc='upper left')
plt.tight_layout()
plt.savefig("barplot_intron_scatter_AS.1.2.png")
plt.savefig("barplot_intron_scatter_AS.1.2.pdf")
#plt.show()