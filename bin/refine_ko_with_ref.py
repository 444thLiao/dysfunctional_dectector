################
"Manual refine the KEGG annotations from KEGG database"

# Input a file with information from GID to abbrev
# todo: replace print with logger
# todo: add comments (optional)
###############

import os
from collections import defaultdict
from os.path import basename, dirname, exists, realpath
from subprocess import check_call

import click
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.KEGG import REST
from tqdm import tqdm

from dysfunctional_dectector.src.utilities.tk import batch_iter, output_file


def check_format(link_file):
    link_df = pd.read_csv(link_file,sep='\t',index_col=0)
    assert list(link_df.columns) == ['infaa','abbrev']
    for gid,row in link_df.iterrows():
        infaa = row['infaa']
        if not exists(realpath(infaa)):
            nowinfaa = realpath(dirname(realpath(link_file))+infaa.strip('.'))
        else:
            nowinfaa = realpath(infaa)
        if not exists(nowinfaa):
            raise IOError(f'can not find {infaa}. even using {nowinfaa}')
        link_df.loc[gid,'infaa'] = nowinfaa
    return link_df

def prepare_abbrev2files(abbrev,outdir):
    try:
        koinfo_str = REST.kegg_link("ko", abbrev).read()
    except :
        raise IOError(f'input wrong abbrev {abbrev}')

    raw_locus = []
    locus2ko = {}
    for row in koinfo_str.strip().split('\n'):
        l,k = row.split('\t')
        locus2ko[l.split(':')[-1]] = k.split(':')[-1]
        raw_locus.append(l)
    ref_p = f'{outdir}/{abbrev}.faa'
    #### get aa seqs
    if not exists(ref_p):
        print(f"Downloading {abbrev} protein sequences.")
        aa = []
        for each10_locus in tqdm(batch_iter(raw_locus,10),total=int(len(raw_locus)/10)):
            try:
                aa_str = REST.kegg_get(each10_locus,'aaseq').read().split('\n')
                # maybe lack of aaseq. such as sedi:EBB79_14800 which is a tRNA
            except:
                # print(each10_locus)
                continue
            aa+=aa_str
        aa = [_ for _ in aa if _]
        with open(ref_p,'w') as f1:
            f1.write('\n'.join(aa))
        seqs = list(SeqIO.parse(ref_p,'fasta'))
        for seq in seqs:
            seq.id = seq.id.split(':')[-1]
        with open(ref_p,'w') as f1:
            SeqIO.write(seqs,f1,'fasta-2line')
    return ref_p,locus2ko

def main(kegg_df,gid2info,ofile,outdir,strict_mode):
    for gid,row in gid2info.iterrows():
        infaa = row['infaa']
        abbrev = row['abbrev']
        if str(abbrev)=='nan':
            continue
        cmd = f"makeblastdb -in {infaa} -dbtype prot -out {outdir}/{gid}; "
        if not exists(f"{outdir}/{gid}.phr"):
            check_call([cmd],shell=1)
            
        gid_ko2locus = {k:v for k,v in kegg_df.loc[gid,:].to_dict().items() if str(v)!='nan'}
        _l2ko = {}
        for ko,locuslist_str in gid_ko2locus.items():
            for locus in locuslist_str.split(','):
                _l2ko[locus] = ko
        ref_p,abbrev_locus2ko = prepare_abbrev2files(abbrev,outdir)
        refined_count = 0
        cmd = f"blastp -query {ref_p} -db {outdir}/{gid} -outfmt 6 -qcov_hsp_perc 90 -evalue 1e-3 -out {outdir}/{gid}_{abbrev}.tab"
        print(f"perform Blastp {gid} vs {abbrev}")
        check_call([cmd],shell=1)
        # parse blastp table
        gid2abbrev_df = pd.read_csv(f'{outdir}/{gid}_{abbrev}.tab',sep='\t',header=None)
        l2l = {}
        for _,row in gid2abbrev_df.iterrows():
            if row[2]>=99:
                # over 99% identity, found nearly exactly the same
                l2l[row[0]] = row[1]
                # LOCUS to LOCUS: abbrev ID to infaa ID
            
        infaa_locus2corrected_ko = {}
        for l1,l2 in l2l.items():
            ko1 = abbrev_locus2ko.get(l1,'')
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
            
        print(f"Use {abbrev} to correct {gid} annotations: {refined_count}/{len(_l2ko)} ")

        if strict_mode:
            masked_count = 0
            for ko in kegg_df.columns:
                if ko not in ko2l:
                    kegg_df.loc[gid,ko] = np.nan
                    masked_count+=1
            print(f"Use {abbrev} to masked {masked_count} KO of {gid} annotations. ")
    output_file(ofile,kegg_df)


@click.command()
@click.option("-i", "link_file", help="input file linking between KEGG abbrev and Genome ID. such as 'sit\tGCA_0011965.1' ")
@click.option("-o", "outdir", help="Output directory ")
@click.option("-k", "KEGG_df_path", help="Should be a tab-delimeter file, which use Genome ID as index and KO as columns name. The values should be comma-separated locus name. Should be the same as the infaa in link file.")
@click.option("-ss","strict_mode",help="Strict mode: default is OFF. If it is on, it means that abbrev and the infaa should be highly similar. The ko not found in abbrev will be masked into nan for infaa.",default=False,required=False,is_flag=True,)
def cli(link_file, KEGG_df_path, outdir,strict_mode):
    if not exists(outdir):
        os.makedirs(outdir)
    ofile = f"{outdir}/{basename(KEGG_df_path).rpartition('.')[0]+'_refined.'+basename(KEGG_df_path).rpartition('.')[-1]}"
    kegg_df = pd.read_csv(KEGG_df_path,sep='\t',index_col=0)
    gid2info = check_format(link_file)
    shared_gid = set(gid2info.index).intersection(kegg_df.index)
    print(f"Found {gid2info.shape[0]} in link_file and {len(shared_gid)} shared with KEGG_df.")
    gid2info = gid2info.reindex(shared_gid)
    main(kegg_df,gid2info,ofile,outdir,strict_mode)

# python3 refine_ko_with_ref.py -i ../Example/refine_ko_withref/link_file -o ../Example/refine_ko_withref/output -k ../Example/refine_ko_withref/KEGG_anno.tsv 
# link_file, KEGG_df_path, outdir = "../Example/refine_ko_withref/link_file","../Example/refine_ko_withref/KEGG_anno.tsv","../Example/refine_ko_withref/output"
if __name__ == '__main__':
    cli()