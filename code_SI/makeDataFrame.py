from pathlib import Path
import pandas as pd
import argparse
import re

parse = argparse.ArgumentParser()
parse.add_argument("--deseq2_out")
parse.add_argument("--JUM_out")
args = parse.parse_args()



# getting percent intron per genes
df_intron = pd.read_csv("gene_introns_exons_span.tsv", sep="\t",
                       names=["chr", "FgID", "Symbol", "transcriptID", "Gene_size", "exon_size", "ratio"])
df_intron = df_intron.iloc[1:]
df_intron["percent_intron"] = df_intron.apply(lambda x : float(x["ratio"])* 100, axis=1)

# DE DESeq2 results
df = pd.read_csv(args.deseq2_out)

df  = df.merge(df_intron, left_on="FlyBase ID", right_on="FgID")
df["Gene"] = [x.strip() for x in df["Gene"]]
df = df[df["percent_intron"] > 0] # at least one intron!

# Reading JUM output fuction
JUMout = Path(args.JUM_out)
def parse_f(f):
    set_ = set()
    with open(f) as f_in:
        for l in f_in:
            if l.startswith("Gene	AS_event_ID"):
                continue
            spt = l.strip().split()
            set_.add(spt[0])
    return set_


# getting the JUM data i.e. which genes are AS, already significant pval 0.05
dico = {}
for d in JUMout.glob("*control"):
    name = d.stem.split("_")[0]
    set_c = set()
    for dd in d.iterdir():
        set_c = set_c.union(parse_f(dd))
    dico[name] = set_c

# we upload the info of AS to the dataframe with genes intron percent
for n in ["SRPK", "U2af38"]:
    df["AS_{}".format(n)] = df.apply(lambda x: 1 if x["Gene"] in dico[n] else 0, axis=1)

    
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

dico2 = {"SRPK": gene_AS_SRPK, "U2af38":gene_AS_U2af38}

for n in ["SRPK", "U2af38"]:
    df["AS_qval_{}".format(n)] = df.apply(lambda x: 1 if x["Gene"] in dico2["{}".format(n)] else 0, axis=1)

df.to_csv("AS_DE_table.csv")