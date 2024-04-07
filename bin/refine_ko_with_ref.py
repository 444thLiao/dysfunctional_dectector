################
"Manual refine the KEGG annotations from KEGG database"

# Input a file with information from GID to abbrev

###############

import click
from subprocess import check_call
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
from Bio import SeqIO
from Bio.KEGG import REST
import pandas as pd
import numpy as np
from os.path import exists


def batch_iter(iter, batch_size):
    # generating batch according batch_size
    iter = list(iter)
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d : len(iter) + 1])
    return n_iter


kegg_df = ''
link_file = ''
# Genome ID\tinfaa\tabbrev
gid = ''
outdir = ''
strict_mode = True  # If it is on, it means that abbrev and the infaa should be highly similar. The ko not found in abbrev will be masked into nan for infaa.

### static setting




@click.command()
@click.option("-i", "link_file", help="input file linking between KEGG abbrev and Genome ID. such as 'sit\tGCA_0011965.1' ")
@click.option("-o", "outdir", help="Output directory ")
@click.option("-k", "KEGG_df_path", help="Should be a tab-delimeter file, which use Genome ID as index and KO as columns name. The values should be comma-separated locus name. Should be the same as the infaa in link file.")
@click.option("-ss","strict_mode",help="Strict mode: default is OFF",default=False,required=False,is_flag=True,)
def cli(link_file, KEGG_df_path, outdir,strict_mode):
    kegg_df = pd.read_csv(KEGG_df_path,sep='\t',index_col=0)
    gid2info = pd.read_csv(link_file,sep='\t',index_col=0)
    shared_gid = set(gid2info.index).intersection(kegg_df.index)
    print(f"Found {gid2info.shape[0]} in link_file and {len(shared_gid)} shared with KEGG_df.")
    gid2info = gid2info.reindex(shared_gid)

    
def main():
    ref_p = f'{outdir}/{abbrev}.faa'
    refined_count = 0
    ## check abbrev
    try:
        koinfo_str = REST.kegg_link("ko", abbrev).read()
    except :
        raise IOError(f'input wrong abbrev {abbrev}')

    locus2ko = {}
    for row in koinfo_str.strip().split('\n'):
        l,k = row.split('\t')
        locus2ko[l.split(':')[-1]] = k.split(':')[-1]
    ko2locus = defaultdict(list)
    for l,k in locus2ko.items():
        ko2locus[k].append(l)


    gid_ko2locus = {k:v for k,v in kegg_df.loc[gid,:].to_dict().items() if str(v)!='nan'}
    _l2ko = {}
    for ko,locuslist_str in gid_ko2locus.items():
        for locus in locuslist_str.split(','):
            _l2ko[locus] = ko

    cmd = f"makeblastdb -in {infaa} -dbtype prot -out {outdir}/{gid}; "
    check_call([cmd])

    #### get aa seqs

    if not exists(ref_p):
        print(f"Downloading {abbrev} protein sequences.")
        aa = ''
        for each10_locus in batch_iter(ko2locus,10):
            aa_str = REST.kegg_get(each10_locus,'aaseq').read().split('\n')
            aa+=aa_str
        with open(ref_p,'w') as f1:
            f1.write(aa)
        f = SeqIO.parse(ref_p,'fasta')
        seqs = []
        for seq in seqs:
            seq.id = seq.split(':')[-1]
        with open(ref_p,'w') as f1:
            SeqIO.write(seqs,f1,'fasta-2line')

    cmd = f"blastp -in {ref_p} -query {outdir}/{gid} -outfmt 6 -out {outdir}/{gid}_{abbrev}.tab"
    check_call(cmd)

    # parse
    gid2abbrev_df = pd.read_csv(f'{outdir}/{gid}_{abbrev}.tab',sep='\t',header=None)
    # gid2abbrev_df.sort_values(2,ascending=False).groupby(0).head(1)
    l2l = {}
    for _,row in gid2abbrev_df.iterrows():
        if row[2]>=99:
            # over 99% identity 
            l2l[row[0]] = row[1]
            # LOCUS to LOCUS: abbrev ID to infaa ID
            
    infaa_locus2corrected_ko = {}
    for l1,l2 in l2l.items():
        ko1 = locus2ko.get(l1,'')
        ko2 = _l2ko.get(l2,'')
        if ko1!='' and ko1!=ko2:
            infaa_locus2corrected_ko[l2] = ko1
            refined_count+=1
            
    ko2l = defaultdict(list)
    for infaa_locus,ko in infaa_locus2corrected_ko.items():
        ko2l[ko].append(infaa_locus)
    ko2l = {k:','.join(sorted(v)) for k,v in ko2l.items()}
    for ko in ko2l:
        kegg_df.loc[gid,ko] = ko2l[ko]
        
    print(f"In total: {refined_count}/{len(_l2ko)} locus are corrected. ")

    if strict_mode:
        masked_count = 0
        for ko in kegg_df.columns:
            if ko not in ko2l:
                kegg_df.loc[gid,ko] = np.nan
                masked_count+=1
        print(f"In total: {masked_count} KO are masked. ")
##### main

if __name__ == '__main__':
    cli()