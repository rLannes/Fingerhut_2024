import argparse
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
import numpy as np
from pathlib import Path
import pickle
import pysam



parse = argparse.ArgumentParser()
parse.add_argument("--bamDir")
args = parse.parse_args()
bam_dir = Path(args.bamDir)


with open("dmel-all-r6.48.dico.pkl", "rb") as f_o:
    dico_gtf = pickle.load(f_o)

    
    
def get_tr(gene, dico_gtf=dico_gtf):
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
    #max_reads = get_total_reads(bam.replace("Aligned_rev_sorted.bam", "Log.final.out"))

    with pysam.AlignmentFile(bam, "rb") as bamfile:
        for pileupcolumn in bamfile.pileup(chr, start, end):
            if start <= pileupcolumn.pos < end:
                cpt = 0
                for pileupread in pileupcolumn.pileups:
                    0
                    if not pileupread.is_del and not pileupread.is_refskip:
                        cpt += 1
                    else:
                        pass
                pcol.append((pileupcolumn.pos, cpt))
                #pcol.append((pileupcolumn.pos, cpt/max_reads * 1_000_000))
        #sum_read = sum([x[1] for x in pcol])
        #pcol = [(x[0], x[1] / sum_read) for x in pcol]
    return pcol


def plot_coverage(gene, bams, labels, color, out, bg_col, N=None, ratio = 0.3, max_reads=None):
    """
    gene the gene to plot
    bams list of bam file
    corresponding names of the bam
    color for the plot in order of bam
    N smoothing windows
    """

    # get tr
    transcript = get_tr(gene)
    strand = dico_gtf[gene]["strand"]
    
    # get all exon in + strand order
    #print(gene, transcript)
    exons = sorted(dico_gtf[gene]["transcript"][transcript]["exon"], key = lambda x: x["end"]) 

    dico_coverage = {}
    for i, v in enumerate(bams):
        genotype = labels[i]
        dico_coverage[genotype] = []
        for i, e in enumerate(exons):
            cover = get_cover(bam=v, chr="chr"+dico_gtf[gene]["chr"], start=e["start"], end=e["end"])
            dico_coverage[genotype].extend(cover)
            print(genotype, "exon: {}, start depth : {} end depth: {}".format(i, cover[0], cover[1]) )

    fig , ax = plt.subplots(figsize=(12, 4))
    legend_lines = []
    for i, v in enumerate(labels):

        cov = dico_coverage[v]
        y = [x[1] for x in sorted(cov, key=lambda x: x[0])]
        # smoothing
        if N:
            y  = np.convolve(y, np.ones(N)/N, mode='valid')
        # reverse if needed and plot
        if strand == "-":
            y = y[::-1]
        max_re = max_reads[v]
        y = [x / max_re for x in y]
        ax.plot(list(range(len(y))), y, color=color[i])

        plt.fill_between(range(len(y)), y, 0, color=color[i], alpha=.15)

        legend_lines.append(Line2D([0], [0], label=v, color=color[i]))
        # no need to plot intron in this case
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend(legend_lines)
    plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')
   
# 
    plt.gca().spines[['right', 'top']].set_visible(False)
    xleft, xright = plt.gca().get_xlim()
    ybottom, ytop = plt.gca().get_ylim()
    plt.gca().set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)

    plt.margins(x=0, tight=True)
    plt.margins(y=0, tight=True)
    plt.tight_layout()
    plt.gca().set_facecolor(bg_col)
    plt.savefig(out + ".png")
    plt.savefig(out + ".pdf")


def plot_coverage_overctl(gene, bams, labels, color, out, bg_col,\
                           N=None, ratio = 0.3, max_reads=None):
    """
    gene the gene to plot
    bams list of bam file
    corresponding names of the bam
    color for the plot in order of bam
    N smoothing windows
    """

    # get tr
    transcript = get_tr(gene)
    strand = dico_gtf[gene]["strand"]
    
    # get all exon in + strand order
    exons = sorted(dico_gtf[gene]["transcript"][transcript]["exon"], key = lambda x: x["end"]) 

    dico_coverage = {}
    for i, v in enumerate(bams):
        genotype = labels[i]
        dico_coverage[genotype] = []
        for e in exons:
            dico_coverage[genotype].extend(get_cover(bam=v, chr="chr"+dico_gtf[gene]["chr"], start=e["start"], end=e["end"]))

    fig , ax = plt.subplots(figsize=(12, 4))
    legend_lines = []
    ctrl_depth = max_reads["Control"]
    control = dico_coverage["Control"]
    control = [x[1] for x in sorted(control, key=lambda x: x[0])]
    control = np.array(control)
    control = control / ctrl_depth

    for i, v in enumerate(labels):

        cov = dico_coverage[v]
        y = np.array([x[1] for x in sorted(cov, key=lambda x: x[0])])
        y = y / max_reads[v]

        y = list(y/ control)
        # smoothing
        if N:
            y  = np.convolve(y, np.ones(N)/N, mode='valid')
        # reverse if needed and plot
        if strand == "-":
            y = y[::-1]

        ax.plot(list(range(len(y))), y, color=color[i])
        #plt.fill_between(range(len(y)), y, 0, color=color[i], alpha=.15)
        p = np.polyfit(list(range(len(y))), y, deg=1)
        y_ = [p[1] + (p[0] * c) for c in range(len(y))]
        plt.plot(range(len(y_)), y_, color=color[i], lw=1)

        legend_lines.append(Line2D([0], [0], label=v, color=color[i]))
        # no need to plot intron in this case
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend(legend_lines)
    plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')
   
# 
    plt.gca().spines[['right', 'top']].set_visible(False)
    plt.gca().set_ylim((0, 1.2))
    xleft, xright = plt.gca().get_xlim()
    ybottom, ytop = plt.gca().get_ylim()
    plt.gca().set_aspect(abs((xright - xleft) / (ybottom - ytop)) * ratio)

    plt.margins(x=0, tight=True)
    plt.margins(y=0, tight=True)
    plt.tight_layout()
    plt.gca().set_facecolor(bg_col)
    plt.savefig(out + ".png")
    plt.savefig(out + ".pdf")



bam_dir = str(bam_dir)
fwd_bam = [bam_dir + "/Control_merge_fwd_sorted.bam", bam_dir + "/U2af38_merge_fwd_sorted.bam", bam_dir + "/SRPK_merge_fwd_sorted.bam"]
rev_bam = [bam_dir + "/Control_merge_rev_sorted.bam", bam_dir + "/U2af38_merge_rev_sorted.bam", bam_dir + "/SRPK_merge_rev_sorted.bam"]


# X = ["flagstat_Control_merge_fwd.txt",
# "flagstat_Control_merge_rev.txt",
# "flagstat_SRPK_merge_rev.txt",
# "flagstat_SRPK_merge_fwd.txt",
# "flagstat_U2af38_merge_fwd.txt",
# "flagstat_U2af38_merge_rev.txt"]

# dico_max_read = {}
# for file in X:
#     f_spt = file.split("_")
#     genotype = f_spt[1]
#     if genotype not in dico_max_read:
#         dico_max_read[genotype] = 0
#     #strand = "+" if f_spt[3] == "fwd" else "-"
#     with open(bam_dir + file) as f_in:
#         dico_max_read[genotype] += int(f_in.readline().strip().split()[0])


dico_color = {"Control": "#009E73",
              "U2af38": "#E69F00",  # orange
              "U2af50": "navy",  #### "#56B4E9", # blue
              "SRPK": "#CC79A7"}  # vermillon
dico_read_depth = {"Control" : (79262356 + 81332374) / 10_000_000,\
                    "SRPK" : (84650810 + 88023702) / 10_000_000,\
                    "U2af38": (70629264 + 76610096) / 10_000_000 }

#kl3
g="FBgn0267432"
print("kl-3")
plot_coverage_overctl(g,
            fwd_bam, ["Control", "U2af38", "SRPK"],
            color=["#009E73", "#E69F00", "#CC79A7"],
            out="Jackiefig5_kl3_over_ctl",bg_col="whitesmoke",
              N=25, max_reads = dico_read_depth)

plot_coverage(g,
            fwd_bam, ["Control", "U2af38", "SRPK"],
            color=["#009E73", "#E69F00", "#CC79A7"],
            out="Jackiefig5_kl3",bg_col="whitesmoke",
              N=25, max_reads = dico_read_depth)


#kl2
g="FBgn0001313"
print("kl-2")
plot_coverage_overctl(g,
            fwd_bam, ["Control", "U2af38", "SRPK"],
            color=["#009E73", "#E69F00", "#CC79A7"],
            out="Jackiefig5_kl2_over_ctl",bg_col="whitesmoke",
              N=25, max_reads = dico_read_depth)

plot_coverage(g,
            fwd_bam, ["Control", "U2af38", "SRPK"],
            color=["#009E73", "#E69F00", "#CC79A7"],
            out="Jackiefig5_kl2",bg_col="whitesmoke",
              N=25, max_reads = dico_read_depth)