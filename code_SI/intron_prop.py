import pickle

# gtf as a dict
with open("dmel-all-r6.48.dico.pkl", "rb") as f_o:
    gtf = pickle.load( f_o)
    
    
def get_longest_transcript(x):
    longest = ""
    size = 0
    for trans, val in x["transcript"].items():
        if not longest:
            longest = val
            size =  val["mRNA"][0]["end"] - val["mRNA"][0]["start"]
            continue
        if val["mRNA"][0]["end"] - val["mRNA"][0]["start"] > size:
            size = val["mRNA"][0]["end"] - val["mRNA"][0]["start"]
            longest = val
    return longest

def get_exon_span(x):
    s = 0
    for ex in x["exon"]:
        s+= ex["end"] - ex["start"]
    return s

def tr_span(x):
    return x["mRNA"][0]["end"] - x["mRNA"][0]["start"]


with open("gene_introns_exons_span.tsv", "w") as f_o:
    f_o.write("\t".join(["chr", "gene_id", "symbol", "tr_id", "gene_size", "intron_size", "intron_prop"]) + "\n")
    for gene_id, sub in gtf.items():
        try:
            if sub["transcript"]:
                longest = get_longest_transcript(sub)
            else:
                continue
        except:
            continue
            
        try:
            exon_s = get_exon_span(longest)
            gene_s = tr_span(longest)
        except:
            raise

        prop = (gene_s - exon_s) / gene_s
        f_o.write("\t".join([sub["chr"], gene_id, sub["symbol"], longest["transcript_id"], str(gene_s), str(gene_s - exon_s), str(prop)]) + "\n")
