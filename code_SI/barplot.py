import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from pathlib import Path
import re

df = pd.read_csv("AS_DE_table.csv", index_col=0)


parse = argparse.ArgumentParser()
parse.add_argument("--JUM_out")
args = parse.parse_args()
JUMout = Path(args.JUM_out)

def get_gene_list(file, q_val=0.05):
    res = []
    with open(file) as f_in:
        x = f_in.readline()
        for l in f_in:
            spt = l.strip().split()
            if spt[0:3] == ['Gene', 'AS_event_ID', 'chromosome']:
                continue
            try:
                if float(spt[-2]) < q_val:
                    res.append(spt[0])
            except:
                print(spt, file)
                raise
    return res

def parse_feature(dir, all_genes= None):
    dico = {}

    for f in dir.glob("*simplified.txt"):
        gene = list(set(get_gene_list(f)))
        typ = re.search("output_(.*)_pvalue", str(f))
        if all_genes:
            gene = [x for x in gene if x in all_genes]
        dico[typ.group(1)] = gene

    return dico


# for this plot we need to parse jum output as we look at the type of AS event.
# we use q_val < 0.05


dico = {}
for dir_ in JUMout.glob("*control"):
    name = dir_.parts[-1]
    dico[name] = parse_feature(dir_, list(df["Gene"]))
    
gene_AS_SRPK = []
for k,v in dico['SRPK_vs_control'].items():
    gene_AS_SRPK.extend(v)

gene_AS_U2af38 = []
for k,v in dico['U2af38_vs_control'].items():
    gene_AS_U2af38.extend(v)

order_AS = sorted(list(dico["SRPK_vs_control"].keys()))

order_genotype = ["U2af38_vs_control", "SRPK_vs_control"] #sorted(list(dico.keys()))
fig, ax = plt.subplots(1,1)

bottom = [0] * len(order_genotype)
color = ["#ef476f", "#ffd166", "#06d6a0", "#118ab2", "#073b4c", "#e26d5c"]
cpt = 0
for AS in order_AS:
    x=[]
    for genotype in order_genotype:
        x.append(len(dico[genotype].get(AS, 0)))
    ax.bar(range(len(order_genotype)), x, bottom=bottom, color=color[cpt])
    bottom = [x[i] + v for i, v in enumerate(bottom)]
    cpt += 1
ax.set_xticks(range(len(order_genotype)))
ax.set_xticklabels([x.split("_")[0] for x in order_genotype], rotation=45, ha="right")

custom_lines = [Line2D([0], [0], color=color[i], lw=4) for i in range(6)]
ax.legend(custom_lines, [x.replace("events", "").replace("_", " ") for x in order_AS], bbox_to_anchor=(1.05, 1), loc='upper left')

ax.set_ylabel("number of alternative splicing event")
ax.set_xlabel("genotype")
plt.title("number of Alternative Splicing event by genotype")
plt.tight_layout()
plt.savefig("final_figure_S5_1.pdf")
plt.savefig("final_figure_S5_1.png")

plt.close()



v = {}
for genotype, values in dico.items():
    all_ = []
    for kk, vv in values.items():
        all_.extend(vv)
    all_ = list(all_)
    sum_all = len(all_)
    print(sum_all)
    for AS, val in values.items():
        if genotype not in v:
            v[genotype] = {}
        v[genotype][AS] = len(val) / sum_all

fig, ax = plt.subplots(1, 1)
bottom = [0] * len(order_genotype)
color = ["#ef476f", "#ffd166", "#06d6a0", "#118ab2", "#073b4c", "#e26d5c"]
cpt = 0
for AS in order_AS:
    x=[]
    for genotype in order_genotype:
        x.append(v[genotype][AS])
    ax.bar(range(len(order_genotype)), x, bottom=bottom, color=color[cpt])#, cmap="pastel")#, label=order_AS)
    bottom = [x[i] + v for i, v in enumerate(bottom)]
    cpt += 1
ax.set_xticks(range(len(order_genotype)))
ax.set_xticklabels([x.split("_")[0] for x in order_genotype], rotation=45, ha="right")

custom_lines = [Line2D([0], [0], color=color[i], lw=4) for i in range(6)]
ax.legend(custom_lines, [x.replace("events", "").replace("_", " ") for x in order_AS], bbox_to_anchor=(1.05, 1), loc='upper left')

ax.set_ylabel("proportion of alternative splicing event")
ax.set_xlabel("genotype")
plt.title("Proportion of Alternative Splicing event by genotype")
plt.tight_layout()
plt.savefig("final_figure_S5_2.pdf")
plt.savefig("final_figure_S5_2.png")
plt.close()

for genotype, values in dico.items():
    all_ = []
    for kk, vv in values.items():
        all_.extend(vv)
    all_ = list(set(all_))
    sum_all = len(all_)
    print(genotype, sum_all)

    
plt.bar(x=[1,2], height=[1420, 4569])
plt.xticks([1,2], labels=["SRPK", "U2af38"], fontsize="large")
plt.text(x=0.93, y=1420+100, s="1420", fontsize="large")
plt.text(x=1.93, y=4569+100, s="4569", fontsize="large")
plt.ylim((0, 5000))
plt.tight_layout()
plt.savefig("final_figure_S5_3.pdf")
plt.savefig("final_figure_S5_3.png")


