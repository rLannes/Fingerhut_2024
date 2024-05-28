import argparse
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import pandas as pd



parse = argparse.ArgumentParser()
parse.add_argument("--deseq2_out")
args = parse.parse_args()

# load the data
# Deseq2 Output
file = args.deseq2_out
df = pd.read_csv(file)

df["Gene"] = df.apply(lambda x : x["Gene"].strip(), axis=1)


# Y-linked genes with large intron to put in Red and bigger.
list_y_genes = ["kl-2", "kl-3", "ORY", "Ppr-Y", "kl-5", "CCY", "WDY", "PRY"]

# gene target by RNAi blue and bigger
rnai_ = ["SRPK"]

y =  df.sort_values("padj.2")[["Gene", "log2(SRPK/Control)", "padj.2"]]#y = y[~np.isnan(y)]
Ylinked = y[y["Gene"].isin(list_y_genes)]
rnai = y[y["Gene"].isin(rnai_)]

fig , ax = plt.subplots(figsize=(4*3, 3*2.5))
subdf = df.sort_values("padj.2")[["Gene", "log2(SRPK/Control)", "padj.2"]]#y = y[~np.isnan(y)]
subdf = subdf.dropna()

x = subdf["log2(SRPK/Control)"]
y = subdf["padj.2"]
min_ = 0

for e in y:
    print(e)
    if e > min_:
        min_ = e
        break

Ylinked = subdf[subdf["Gene"].isin(list_y_genes)]
rnai = subdf[subdf["Gene"].isin(rnai_)]

subdf = subdf[~subdf["Gene"].isin(rnai_ + list_y_genes )]
x = subdf["log2(SRPK/Control)"]
y = subdf["padj.2"]   

print(min_)

y = [-math.log10(x) if x > 0 else -math.log10(min_) for x in y ]
ax.scatter(x, y, color="darkgray", alpha = 0.6, zorder=5, s=60, linewidths=0)

ax.scatter(x = rnai["log2(SRPK/Control)"], y = -math.log10(rnai["padj.2"]),
           color="royalblue", alpha = 1, zorder=10, s=130, linewidths=0)

ax.text(rnai["log2(SRPK/Control)"]+0.2, -math.log10(rnai["padj.2"])+5, rnai_[0], zorder=15, fontsize=20)


ax.scatter(x = Ylinked["log2(SRPK/Control)"], y = [-math.log10(x) for x in Ylinked["padj.2"]],
           color="crimson", alpha = 1, zorder=10, s=130, linewidths=0)

for r, i in Ylinked.iterrows():
    ax.text(i["log2(SRPK/Control)"]+0.2, -math.log10(i["padj.2"])+5, i["Gene"], zorder=15, fontsize=20)


#plt.yticks(range(0,250, 50), list(range(0,250, 50)))
ax.yaxis.set_major_locator(MultipleLocator(50))
ax.yaxis.set_major_formatter('{x:.0f}')
ax.yaxis.set_minor_locator(MultipleLocator(25))

ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_major_formatter('{x:.0f}')
ax.xaxis.set_minor_locator(MultipleLocator(2.5))

plt.grid(visible=True, which='major',
         axis='both',
         color="black", alpha=0.8,
         linewidth=1, zorder=1)


plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.title("SRPK", fontsize=28)

plt.margins(x=0, tight=True)
plt.margins(y=0, tight=True)

plt.gca().spines[['right', 'top']].set_visible(False)
xleft, xright = plt.gca().get_xlim()
ybottom, ytop = plt.gca().get_ylim()
plt.xlim((xleft-2, xright+2))
plt.ylim((ybottom-10, ytop+10))
plt.ylabel('$-Log_{10}$ adj.P', fontsize=24)
plt.xlabel('$Log_{2}$ fold change', fontsize=24)
# plt.legend(loc="lower left", ncol = len(ax.lines) )
plt.savefig("SRPK_volcanoplot.pdf")
plt.savefig("SRPK_volcanoplot.png")
plt.show()




rnai_ = ["U2af38"]

y =  df.sort_values("padj.3")[["Gene", "log2(U2af38/Control)", "padj.3"]]
Ylinked = y[y["Gene"].isin(list_y_genes)]
rnai = y[y["Gene"].isin(rnai_)]


fig , ax = plt.subplots(figsize=(4*3, 3*2.5))
subdf = df.sort_values("padj.3")[["Gene", "log2(U2af38/Control)", "padj.3"]]
subdf = subdf.dropna()

x = subdf["log2(U2af38/Control)"]
y = subdf["padj.3"]
min_ = 0

for e in y:
    if e > min_:
        min_ = e
        break

Ylinked = subdf[subdf["Gene"].isin(list_y_genes)]
rnai = subdf[subdf["Gene"].isin(rnai_)]

subdf = subdf[~subdf["Gene"].isin(rnai_ + list_y_genes )]
x = subdf["log2(U2af38/Control)"]
y = subdf["padj.3"]   


y = [-math.log10(x) if x > 0 else -math.log10(min_) for x in y ]
ax.scatter(x, y, color="darkgray", alpha = 0.6, zorder=5, s=60, linewidths=0)

ax.scatter(x = rnai["log2(U2af38/Control)"], y = -math.log10(rnai["padj.3"]),
           color="royalblue", alpha = 1, zorder=10, s=130, linewidths=0)

ax.text(rnai["log2(U2af38/Control)"]+0.2, -math.log10(rnai["padj.3"])+5,
        rnai_[0], zorder=15, fontsize=20)


ax.scatter(x = Ylinked["log2(U2af38/Control)"], y = [-math.log10(x) for x in Ylinked["padj.3"]],
           color="crimson", alpha = 1, zorder=10, s=130, linewidths=0)

for r, i in Ylinked.iterrows():
    ax.text(i["log2(U2af38/Control)"]+0.2, -math.log10(i["padj.3"])+5, i["Gene"], zorder=15, fontsize=20)


#plt.yticks(range(0,250, 50), list(range(0,250, 50)))
ax.yaxis.set_major_locator(MultipleLocator(50))
ax.yaxis.set_major_formatter('{x:.0f}')
ax.yaxis.set_minor_locator(MultipleLocator(25))

ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_major_formatter('{x:.0f}')
ax.xaxis.set_minor_locator(MultipleLocator(2.5))

plt.grid(visible=True, which='major',
         axis='both',
         color="black", alpha=0.8,
         linewidth=1, zorder=1)


plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.title("U2af38", fontsize=28)
plt.margins(x=0, tight=True)
plt.margins(y=0, tight=True)

plt.gca().spines[['right', 'top']].set_visible(False)
xleft, xright = plt.gca().get_xlim()
ybottom, ytop = plt.gca().get_ylim()
plt.xlim((xleft-2, xright+2))
plt.ylim((ybottom-10, ytop+10))
plt.ylabel('$-Log_{10}$ adj.P', fontsize=24)
plt.xlabel('$Log_{2}$ fold change', fontsize=24)
# plt.legend(loc="lower left", ncol = len(ax.lines) )
plt.savefig("U2af38_volcanoplot.pdf")
plt.savefig("U2af38_volcanoplot.png")
plt.show()