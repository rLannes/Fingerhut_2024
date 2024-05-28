import argparse
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from pathlib import Path
import pickle
import pysam
import random as rd
import traceback


parse = argparse.ArgumentParser()
parse.add_argument("--bamDir")
parse.add_argument("--Normalized")
args = parse.parse_args()
bam_dir = args.bamDir

colorblind_palette = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "black", "grey"]
with open("dmel-all-r6.48.dico.pkl", "rb") as f_o:
    dico_gtf = pickle.load(f_o)


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
            #print(tr_id, N)
            for i,e in enumerate(N[:-1]):
                x = N[i+1]["start"] - e["end"]
                #print(x)
                if x > 0:
                    intron_.add((N[i+1]["start"], e["end"]))
        else:
            N = sorted(N, key=lambda x: x["end"], reverse=True)
            #print(tr_id, N)
            for i,e in enumerate(N[:-1]):
                x = e["start"] - N[i+1]["end"] 
                #print(x)
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
            #print(tr_id, N)
            for i,e in enumerate(N[:-1]):
                x = N[i+1]["start"] - e["end"]
                #print(x)
                if x > largest:
                    largest = x
        else:
            N = sorted(N, key=lambda x: x["end"], reverse=True)
            #print(tr_id, N)
            for i,e in enumerate(N[:-1]):
                x = e["start"] - N[i+1]["end"] 
                #print(x)
                if x > largest:
                    largest = x

        
    return largest

def get_tr(gene, dico_gtf=dico_gtf):
    dico = {
            "FBgn0000114": "FBtr0080344",
            "FBgn0259176" : "FBtr0300521"} 
    
    if gene in dico:
        return dico[gene]
    genes_def = dico_gtf[gene]
    s = 0
    tr = None

    for t, v in genes_def["transcript"].items():
        if len(v["exon"]) > s:
            tr = t
            s = len(v["exon"])
    return tr


def get_total_reads(file):
    with open(file) as f_in:
        for l in f_in:
            l = l.strip()
            if l.startswith("Uniquely mapped reads number"):
                return int(l.split()[-1])
            

def get_cover(bam, chr, start, end):

    pcol = []
    pileup_col = []
    cpt_ = []

    with pysam.AlignmentFile(bam, "rb") as bamfile:
        for pileupcolumn in bamfile.pileup(chr, start, end):
            if start < pileupcolumn.pos < end:
                cpt = 0
                for pileupread in pileupcolumn.pileups:
                    #0
                    if not pileupread.is_del and not pileupread.is_refskip:
                        cpt += 1
                    else:
                        pass
                pcol.append((pileupcolumn.pos, cpt))
                cpt_.append(cpt)
                pileup_col.append(pileupcolumn.pos)

    cpt_ = np.array(cpt_)
    pileup_col = np.array(pileup_col)
    return (pileup_col, cpt_)


class Coverage():

    def __init__(self) -> None:
        self.cpt = []
        self.pileup = []
    
    def sort(self) -> None:
        self.cpt, self.pileup = (list(t) for t in zip(*sorted(zip(self.cpt, self.pileup))))
        self.pileup = np.array(self.pileup)



def to_windows(y, s=100):

    ratio = len(y) / s
    b = [int(i*ratio) for i in range(s)]
    b.append(len(y))
    y_new = [np.mean(y[x : b[i+1]]) for i, x in enumerate(b[:-1])]
    return y_new


def cleave_utr(exon, utr):
    start_u = utr[0]["start"]
    end_u = utr[-1]["end"]
    
    r = []
    if start_u < exon[0]["end"]:
        for e in exon:
            if e["start"] > end_u:
                r.append(e)
            elif e["end"] < end_u:
                continue
            else:
                e["start"] = end_u
                r.append(e)
    else:
        for e in exon:
            if e["end"] < start_u:
                r.append(e)
            elif e["start"] > start_u:
                continue
            else:
                e["end"] = start_u
                r.append(e)
    return r 


dico_read_depth = {"Control" : (79262356 + 81332374) / 10_000_000,\
                    "SRPK" : (84650810 + 88023702) / 10_000_000,\
                    "U2af38": (70629264 + 76610096) / 10_000_000 }


def plot_coverage(genes, bam_dico, color, out, bg_col, genotype, N=None,
                   ratio = 0.3, max_reads=None, check=True, y_lim=None,
                   add_slope=False):
    """
    gene the gene to plot
    bams list of bam file
    corresponding names of the bam
    color for the plot in order of bam
    N smoothing windows
    """
    slope = []
    print(max_reads)

    print()

    cpt = 0
    dico_coverage = {}
    for gene in genes: 
        # get tr
        print(gene)
        transcript = get_tr(gene)
        strand = "fwd" if dico_gtf[gene]["strand"] == "+" else "rev"
        
        # get all exon in + strand order
        exons = sorted(dico_gtf[gene]["transcript"][transcript]["exon"], key = lambda x: x["end"]) 
        utr5 = dico_gtf[gene]["transcript"][transcript].get("5UTR")
        utr3 = dico_gtf[gene]["transcript"][transcript].get("3UTR") 
        print(exons)
        if utr5:
            utr5 = sorted(utr5, key = lambda x: x["end"]) 
            exons = cleave_utr(exons, utr5)
        if utr3:
            utr3 = sorted(utr3, key = lambda x: x["end"])
            exons = cleave_utr(exons, utr3)

        dico_coverage[gene] = []
        #control = []
        control = Coverage()
        for e in exons:
            exon_cov = get_cover(bam=bam_dico[strand]["ctl"], chr="chr"+dico_gtf[gene]["chr"],
                                      start=e["start"], end=e["end"])
            
            control.cpt.extend(exon_cov[0])
            control.pileup.extend(exon_cov[1])
        
        control.pileup = np.array(control.pileup)
        if check and (control.pileup.shape == (0,) or np.mean(control.pileup) < 30 or not np.all(control.pileup)):
            continue

        print(transcript, gene,  dico_gtf[gene]["symbol"])
        this_genotype = Coverage()
        for e in exons:
            exon_cov=get_cover(bam=bam_dico[strand][genotype], chr="chr"+dico_gtf[gene]["chr"],
                                            start=e["start"], end=e["end"])

            this_genotype.cpt.extend(exon_cov[0])
            this_genotype.pileup.extend(exon_cov[1])

        this_genotype.pileup = np.array(this_genotype.pileup)

        print(control.pileup)
        print(this_genotype.pileup)

        if check and this_genotype.pileup.shape == (0,):
            continue

        print(this_genotype.pileup.size, this_genotype.pileup.shape)
        if check and np.mean(this_genotype.pileup) < 30:
            continue
        
        this_genotype.sort()
        control.sort()

        cover_this = this_genotype.pileup / max_reads[genotype]
        cover_control = control.pileup / max_reads["Control"]
        X_ = cover_this / cover_control

        dico_coverage[gene] = list(X_)
        cpt += 1

        if cpt >= 15:
            break
            

    s = 100
    fig , ax = plt.subplots(figsize=(12, 4))
    legend_lines = []
    i = 0
    for v, cov in dico_coverage.items():
        try:    
            cov = dico_coverage[v]
            y = cov 
            if all(i_ != i_ for i_ in y):
                print(y)
                raise AssertionError("{} nan only".format(v))
            # smoothing
            if N:
                y  = np.convolve(y, np.ones(N)/N, mode='valid')
            # reverse if needed and plot
            if strand == "rev":
                y = y[::-1]

            y = to_windows(y)  
            ax.plot(list(range(len(y))), y, color=color[i])
            legend_lines.append(Line2D([0], [0], label=dico_gtf[v]["symbol"], color=color[i]))

            if add_slope:
                y = np.array(y)
                Y_ = y[~np.isnan(y)]
                Y_ = Y_[~np.isinf(Y_)]

                p = np.polyfit(list(range(len(Y_))), Y_, deg=1)
                print(p, len(y), Y_)
                slope.append(p[0])

                y_ = [p[1] + (p[0] * c) for c in range(len(y))]
                plt.plot(range(len(y_)), y_, color=color[i], lw=1)
                y = list(y)
            i += 1

        except Exception as e:
            exc_obj = e
            tb_str = ''.join(traceback.format_exception(None, exc_obj, exc_obj.__traceback__))
            print(v, tb_str, "continue")
            continue
        #plt.fill_between(range(len(y)), y, 0, color=color[i], alpha=.15
        # no need to plot intron in this case
        
        if i > 30:
            break

    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend(legend_lines)
    plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.axhline(y=1, xmin=0, xmax=1, color="black",linestyle="--")

    plt.gca().spines[['right', 'top']].set_visible(False)
    xleft, xright = plt.gca().get_xlim()
    ybottom, ytop = plt.gca().get_ylim()
    #plt.gca().set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)
    if y_lim:
        plt.gca().set_ylim(y_lim)
    plt.margins(x=0, tight=True)
    plt.margins(y=0, tight=True)
    plt.tight_layout()
    plt.gca().set_facecolor(bg_col)
    plt.title(genotype)
    plt.savefig(out + ".png")
    plt.savefig(out + ".pdf")

    return slope


bam_dir = str(bam_dir)
dico_bam = {"fwd":{
    "ctl": bam_dir + "Control_merge_fwd_sorted.bam",
    "U2af38": bam_dir + "U2af38_merge_fwd_sorted.bam",
    "SRPK": bam_dir + "SRPK_merge_fwd_sorted.bam"
}, "rev": {
        "ctl": bam_dir + "Control_merge_rev_sorted.bam",
        "U2af38": bam_dir + "U2af38_merge_rev_sorted.bam",
        "SRPK": bam_dir + "SRPK_merge_rev_sorted.bam"}}
fwd_bam = [bam_dir + "Control_merge_fwd_sorted.bam", bam_dir + "U2af38_merge_fwd_sorted.bam", bam_dir + "SRPK_merge_fwd_sorted.bam"]
rev_bam = [bam_dir + "Control_merge_rev_sorted.bam", bam_dir + "U2af38_merge_rev_sorted.bam", bam_dir + "SRPK_merge_rev_sorted.bam"]


dico_read_depth = {"Control" : (79262356 + 81332374) / 10_000_000,\
                    "SRPK" : (84650810 + 88023702) / 10_000_000,\
                    "U2af38": (70629264 + 76610096) / 10_000_000 }


df = pd.read_csv("AS_DE_table.csv", index_col=0)
ncount = pd.read_csv(args.Normalized)

fb_expr = ncount[ncount["Control_rep1.norm 	Control_rep2.norm 	Control_rep3.norm".split()].sum(axis=1)  > 30]["FlyBaseID"].values
df = df[df["FlyBase ID"].isin(fb_expr)]

dico = {}
intron = []
#
for gene_id, dico_gene in dico_gtf.items():
    if gene_id not in df["FlyBase ID"].values:
        continue
    if dico_gene["biotype"] not in  ["protein_coding_gene"]:
        continue 
        
    strand = dico_gene["strand"]
    l = get_intron(dico_gene["transcript"], strand)
    intron.extend(get_intron_all(dico_gene["transcript"], strand))

    dico[gene_id] = l # largest intron all transcript considered for a genes

gene_with_gap = """FBgn0267428 FBgn0267431 
 FBgn0267430 FBgn0267433 FBgn0267489 FBgn0267432 FBgn0001313 
 FBgn0046697 FBgn0267449 FBgn0046323 FBgn0267592""".split()


dico_slope = {}
jackie_list = """FBgn0283741
FBgn0259176
FBgn0004449
FBgn0000479
FBgn0011817
FBgn0032629""".strip().split()

dico_slope["U2af38"] = {}
dico_slope["U2af38"]["large"] =  plot_coverage(jackie_list, bam_dico=dico_bam, genotype="U2af38",
                color=colorblind_palette,
                out="figure_5_large_U2af38_Final_20_12", bg_col="whitesmoke",
                N=None, ratio = 0.3, max_reads=dico_read_depth,
                check=False, add_slope=True, y_lim=(0, 1.2))


jackie_list = """FBgn0032744
FBgn0034592
FBgn0262592
FBgn0042179
FBgn0036808
FBgn0036512""".strip().split()

dico_slope["U2af38"]["small"] = plot_coverage(jackie_list,
            bam_dico=dico_bam,
            genotype="U2af38",
            color=colorblind_palette,
            out="figure_5_small_U2af38_Final_20_12", bg_col="whitesmoke",
            N=None, ratio = 0.3, max_reads=dico_read_depth,
            check=False, add_slope=True, y_lim=(0, 1.2))


dico_slope["SRPK"] = {}

jackie_list = """FBgn0039395
FBgn0262566
FBgn0262878
FBgn0034882
FBgn0030374
FBgn0036153""".strip().split()


dico_slope["SRPK"]["small"] = plot_coverage(jackie_list,
            bam_dico=dico_bam,
            genotype="SRPK",
            color=colorblind_palette,
            out="figure_5_small_SRPK_Final_20_12", bg_col="whitesmoke",
            N=None, ratio = 0.3, max_reads=dico_read_depth,
            check=False, add_slope=True, y_lim=(0, 1.2))

Jackie_list = ["FBgn0262123",
"FBgn0261260",
"FBgn0004449",
"FBgn0000114",
"FBgn0003435",
"FBgn0267336"]


dico_slope["SRPK"]["large"] = plot_coverage(Jackie_list, bam_dico=dico_bam, genotype="SRPK",
                color=colorblind_palette, #["black", "yellow", "green", "red", "orange", "blue", "grey", "brown", "gold", "pink", "magenta"],
                out="figure_5_large_SRPK_JackieList_12_20", bg_col="whitesmoke",
                N=None, ratio = 0.3, max_reads=dico_read_depth,
                check=False, add_slope=True, y_lim=(0, 1.2))


Gene_gap_tokeep = ["FBgn0267431", "FBgn0267430", "FBgn0263112"]


dico_slope["U2af38"]["Auto w/ gap"] = plot_coverage(Gene_gap_tokeep, bam_dico=dico_bam,
                                             genotype="U2af38",
                color=colorblind_palette,
                out="figure_5_gap_U2af38_Final_20_12",
                  bg_col="whitesmoke", N=None, ratio = 0.3,
                  max_reads=dico_read_depth, check=False, 
                add_slope=True)

dico_slope["SRPK"]["Auto w/ gap"] = plot_coverage(Gene_gap_tokeep, bam_dico=dico_bam, genotype="SRPK",
                color=colorblind_palette,
                out="figure_5_gap_SRPK_Final_20_12",
                  bg_col="whitesmoke", N=None, ratio = 0.3,
                  max_reads=dico_read_depth, check=False, 
                add_slope=True)

Y_genes = ['kl-2', 'kl-3', 'kl-5', 'CCY', 'ORY', 'Ppr-Y']

# symbol -> FbgnSymbol
dico_name = {}
for k, v in dico_gtf.items():
    dico_name[v["symbol"]] = k


dico_slope["U2af38"]["Y-linked"] = plot_coverage([dico_name[x] for x in Y_genes], bam_dico=dico_bam,
             genotype="U2af38",
            color=colorblind_palette,
            out="figure_5_Y_U2af38_Final_20_12",
            bg_col="whitesmoke", N=None, ratio = 0.3,
            max_reads=dico_read_depth,
            check=False, y_lim=(0, 1.2), add_slope=True)

dico_slope["SRPK"]["Y-linked"] = plot_coverage([dico_name[x] for x in Y_genes],
               bam_dico=dico_bam, genotype="SRPK",
            color=colorblind_palette,
            out="figure_5_Y_SRPK_Final_20_12",
            bg_col="whitesmoke",
            N=None, ratio = 0.3, max_reads=dico_read_depth,
            check=False, y_lim=(0, 1.2), add_slope=True)

with open("slope_dico_jackie_fig6.pkl", "wb") as f_o:
    pickle.dump(dico_slope, file=f_o)
