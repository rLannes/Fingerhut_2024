import pickle
import argparse


gtf_file = "dmel-all-r6.48.gtf"
biotype = "dmelanogaster_gene_biotypes_fb_2022_05.tsv"

parse = argparse.ArgumentParser()
parse.add_argument("--gtf")
parse.add_argument("--biotype")
args = parse.parse_args()

gtf_file = args.gtf
biotype = args.biotype 

dd_bio = {}


with open(biotype) as f_in:
    for l in f_in:
        spt = l.strip().split()
        dd_bio[spt[0]] = spt[-1]
    
    
def get_attr(string):
    dico = {}
    spt = string.split(";")
    
    for x in spt:
        if x:
            dico[x.split()[0]] =  x.split()[1].replace('"', "")
    return dico

    
dico = {}
with open(gtf_file) as f_in:
    
    for line in f_in:
        spt = line.strip().split("\t")
        
        chr_ = spt[0]
        start = int(spt[3])
        end = int(spt[4])
        strand =  spt[6]
        attr = get_attr(spt[-1])
        gene_id = attr["gene_id"]
        gene_symbol = attr["gene_symbol"]
        type_ = spt[2]
            
        if type_ == "gene":

            dico[gene_id] = {
                "chr": chr_,
                "start": start, 
                "end" : end,
                "strand" : strand,
                "symbol" :  gene_symbol,
                "biotype" : dd_bio.get(gene_id, "unknown"),
            "transcript": {}}
            
        else:

            transcript_id  = attr["transcript_id"]
            transcript_symbol = attr["transcript_symbol"]
            
            if transcript_id not in dico[gene_id]["transcript"]:
                dico[gene_id]["transcript"][transcript_id] = {
                "transcript_symbol" : transcript_symbol,
                "transcript_id" : transcript_id
                }
            
            if type_ not in dico[gene_id]["transcript"][transcript_id]: 
                dico[gene_id]["transcript"][transcript_id][type_] = []
                
            dico[gene_id]["transcript"][transcript_id][type_].append({
                "chr": chr_,
                "start": start, 
                "end" : end,
                "strand" : strand,
            })

with open("dmel-all-r6.48.dico.pkl", "wb") as f_o:
    pickle.dump(dico, f_o)