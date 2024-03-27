"""
Using reference-based to refine the annotations

Download or prepare a table based on the KEGG genome database, they contain the full list of the presence/absence of a type strain. 



todo:


"""





from glob import glob
import pandas as pd
from tqdm import tqdm
from collections import defaultdict



def get_refined_patterndf(kegg_df,locus_is_pseudo,l2all,locus2ko,
                          targeted_kos,
                          interpro_threshold=0.8,
                          verbose=1):
    confident_presence = defaultdict(dict)
    gid2ko2intact_l = defaultdict(dict)
    gid2ko2pseudo_l = defaultdict(lambda :defaultdict(list))
    ko2intact_l = defaultdict(list)
    gid2cases = defaultdict(dict)

    #i = related_kosy
    for ko in kegg_df.columns:
        ## For ko not found in this kegg_df
        if ko not in kegg_df.columns:
            for gid in kegg_df.index:
                gid2cases[gid][ko] = 'no KEGG-annotated'
            continue
        ## For ko in kegg_df, search it for each gid.
        cols = kegg_df[ko].to_dict()
        pseudo_l = []
        intact_l = []
        for gid,locus in cols.items():
            if str(locus)=='nan':
                cases = 'no KEGG-annotated'
                gid2cases[gid][ko] = cases
                continue
            all_l = locus.split(',')
            for l in all_l:
                if l in locus_is_pseudo:
                    pseudo_l.append(l)
                else:
                    intact_l.append(l)
            ko2intact_l[ko].extend(intact_l)
            confident_presence[gid][ko] = intact_l
            if all([_ in locus_is_pseudo for _ in all_l]):
                cases = 'no intact'
                gid2ko2pseudo_l[gid][ko] = [_  for _ in all_l
                                            if _  in pseudo_l]
            else:
                cases = 'intact'
                gid2ko2intact_l[gid][ko] = [_  for _ in all_l
                                            if _ not in pseudo_l]
            gid2cases[gid][ko] = cases
    
    ko2gid2status_df = pd.DataFrame.from_dict(gid2cases)
    if verbose:
        i = tqdm(ko2gid2status_df.iterrows(),
                 total=ko2gid2status_df.shape[0])
    else:
        i = ko2gid2status_df.iterrows()
        
    for ko,row in i:
        for gid,v in row.to_dict().items():
            if v != 'no KEGG-annotated' or ko not in ko2intact_l:
                continue
            if len(ko2intact_l[ko]) == 0:
                continue
            i = ko2intact_l[ko][0]
            _a = gid+'_'
            sub_l2all = {k:v for k,v in l2all.items()
                         if k.startswith(_a)}
            found_count = defaultdict(int)
            for db,db_v in sorted(l2all[i].items()):
                found = {k:v
                     for k,v in sub_l2all.items()
                     if v.get(db,'') ==db_v }
                for k,v in found.items():
                    found_count[k]+=1
            num_db = len(l2all[i])
            found_l = [k for k,v in found_count.items()
                       if v >=int(num_db*interpro_threshold)]
            found_l = [l for l in found_l if locus2ko.get(l,'') in ['',ko]]
            # if found locus belongs to other KO, ignore it
            if len(found_l)==0:
                continue
            if all([_ in locus_is_pseudo for _ in found_l]):
                cases = 'no intact'
                gid2ko2pseudo_l[gid][ko].extend(found_l)
            else:
                cases = 'intact'
            ko2gid2status_df.loc[ko,gid] = f"RE({cases})"
    return ko2gid2status_df,gid2ko2pseudo_l,confident_presence


def is_confident_pseudo(this_locus,thresholds):
    """
    determined whether it meet the edge of the contig
    """
    pseudo_region = genome_pos.loc[this_locus,'pseudogenized']
    strand = genome_pos.loc[this_locus,'strand']
    cover_regions = genome_pos.loc[genome_pos['pseudogenized']==pseudo_region,:]
    # this_lengs = []
    # for _,row in cover_regions.iterrows():
    #     this_lengs.append(abs(row['end']-row['start']))
    if strand == -1:
        upstream = cover_regions.index[-1]
    else:
        upstream = cover_regions.index[0]
    len_of_upstream = cover_regions.loc[upstream,'end'] - cover_regions.loc[upstream,'start']
    #print(len_of_upstream)
    low_thres,max_thres = thresholds
    if len_of_upstream<=low_thres or len_of_upstream>=max_thres:
        return True,len_of_upstream
    else:
        return False,len_of_upstream


def assess_confident_pseudo_multi(gid2ko2pseudo_l,confident_presence,
                                  len_threshold=0.6):
    ko2pseudo_l = defaultdict(list)
    ko2intact_l = defaultdict(list)
    # re_ko2intact_l = defaultdict(list)
    for gid,_d in gid2ko2pseudo_l.items():
        for ko,l in _d.items():
            # sub_l = [_ for _ in l if _ not in locus_is_pseudo]
            # re_ko2intact_l[ko].extend(sub_l)
            sub_l = [_ for _ in l if _ in locus_is_pseudo]
            ko2pseudo_l[ko].extend(sub_l)
    for gid,_d in confident_presence.items():
        for ko,l in _d.items():
            ko2intact_l[ko].extend(l)
            
    pseudo2assess_result = {}
    for ko, locus_l in tqdm(ko2pseudo_l.items()):
        candidate_l = ko2intact_l[ko]
        if len(candidate_l)==0:
            continue
        lens = []
        for l in candidate_l:
            s,e = genome_pos.loc[l,['start','end']]
            _l = abs(e-s)
            lens.append(_l)
        low_thres = max(lens)*len_threshold
        max_thres = max(lens)*len_threshold*2
        for locus in locus_l:
            is_conf,ll = is_confident_pseudo(locus,(low_thres,max_thres))
            if is_conf:
                # remaining locus is longer than 50% of the target gene
                pseudo2assess_result[locus] = ('confident pseudo',ko,(low_thres,max_thres),ll)
            else:
                pseudo2assess_result[locus] = ('likely pseudo',ko,(low_thres,max_thres),ll)
    return pseudo2assess_result

v2values = {'intact':1, 
            'no KEGG-annotated':0, 
            'confident pseudo':0.2,
            'RE(intact)':0.5, 
            'RE(no intact)':0.5,
            'no intact':0.5
            }




####### better to add another one together with the one waitting to be refined
gid2ipr = {f.split('/')[-2].split('.')[0]:f
           for f in glob('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/ipr/*.anno/*.faa.tsv')}
kegg_df = pd.read_csv(
    '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_anno_Revised.tsv', sep='\t', index_col=0)
locus2ko = {locus: ko
            for ko, _d in kegg_df.to_dict().items()
            for genome, locus_list in _d.items()
            for locus in str(locus_list).split(',')}

genome_pos = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/76genomes_pos.tsv',sep='\t',index_col=0)
genome_pos.loc[:,'is pseudo'] = ['Yes' if ':' in str(_) else 'No' for _ in genome_pos['pseudogenized']]
genome_pos.loc[:,'genome'] = [_.split('_')[0] for _ in genome_pos['contig']]
locus2contig = genome_pos['contig'].to_dict()
locus_is_pseudo = set(list(genome_pos.index[genome_pos['pseudogenized'].str.contains(':')]))


headers = ['Protein accession',
 'Sequence MD5 digest',
 'Sequence length',
 'Analysis',
 'Signature accession',
 'Signature description',
 'Start location',
 'Stop location',
 'Score',
 'Status',
 'Date',
 'InterPro accession',
 'InterPro annotations',
 'GO annotations',
 'Pathways annotations']

l2all = defaultdict(dict)
ko2others = defaultdict(dict)
ko2ipr = {}
for gid, ofile in tqdm(gid2ipr.items()):
    ko2locus = {k:v for k,v in kegg_df.loc[gid,:].to_dict().items() if str(v)!='nan' and k not in ko2ipr}
    locus2ko = {}
    for ko,l in ko2locus.items():
        for _ in l.split(','):
            locus2ko[_] = ko
    ipr_df = pd.read_csv(ofile,sep='\t',index_col=None,low_memory=False,names=headers)
    ipr_df.pop('Sequence MD5 digest')
    # ipr_df.pop('Pathways annotations)
    subdf = ipr_df.loc[ipr_df['Protein accession'].isin(locus2ko),:]
    for i,row in subdf.iterrows():
        l = row['Protein accession']
        l2all[l][row['Analysis']] = row['Signature accession']
        ko = locus2ko[l]
        if ko in ko2ipr:
            continue
        if row['InterPro accession'] != '-':
            ko2ipr[ko] = row['InterPro accession']
        else:
            ko2others[ko][row['Analysis']] = row['Signature accession']


related_kos = list(kegg_df.columns)
ko2gid2status_df,gid2ko2pseudo_l,confident_presence = get_refined_patterndf(kegg_df,locus_is_pseudo,l2all,locus2ko,
                                                                            related_kos,verbose=1)
pseudo2assess_result = assess_confident_pseudo_multi(gid2ko2pseudo_l,confident_presence,
                                  len_threshold=0.6)
copy_df = ko2gid2status_df.copy()
for locus,(_,ko,_,_) in pseudo2assess_result.items():
    gid = locus.split('_')[0]
    ori = copy_df.loc[ko,gid]
    if ori == 'RE(no intact)':
        copy_df.loc[ko,gid] = 'confident pseudo'
    
bin_df = copy_df.applymap(lambda x:v2values[x])
bin_df = bin_df.reindex(index=related_kos)
bin_df.to_excel('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/functional_annotations/Forcurated_biotin.xlsx')      

gid2kos = defaultdict(list)
for gid,col in copy_df.iteritems():
    for ko,v in col.to_dict().items():
        if v not in ['confident pseudo','no KEGG-annotated']:
            gid2kos[gid].append(ko)
gid2kos = {k:list(set(v)) for k,v in gid2kos.items()}





################
"Manual refine the KEGG annotations from KEGG database"
###############
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas

cmd = f"makeblastdb -in /mnt/ivy/thliao/project/ruegeria_prophage/data/prokka_o/GNM000011965/GNM000011965.faa -dbtype prot -out GNM0000111965; "

import pandas as pd
import numpy as np

kegg_df = pd.read_csv(
    '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_anno_Revised.tsv', sep='\t', index_col=0)
my_faa = '/mnt/ivy/thliao/project/ruegeria_prophage/data/prokka_o/GNM000014065/GNM000014065.faa'
_ko2locus = {k:v for k,v in kegg_df.loc['GNM000014065',:].to_dict().items() if str(v)!='nan'}
_l2ko = {}
for ko,ll in _ko2locus.items():
    for l in ll.split(','):
        _l2ko[l] = ko
result = REST.kegg_link("ko", "sit").read()
locus2ko = {}
for row in result.strip().split('\n'):
    l,k = row.split('\t')
    locus2ko[l.split(':')[-1]] = k.split(':')[-1]
ko2locus = defaultdict(list)
for l,k in locus2ko.items():
    ko2locus[k].append(l)
ko2ref = {}
ref_p = '/mnt/ivy/thliao/project/coral_ruegeria/tmp/GCA_000014065.1.faa'
_df = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/tmp/GNM000014065.tab',sep='\t',header=None)
_df.sort_values(2,ascending=False).groupby(0).head(1)
l2l = {}
for _,row in _df.iterrows():
    if row[2]>=99:
        l2l[row[0]] = row[1]
new_l2ko = {}
same_c = 0
for l1,l2 in l2l.items():
    ko1 = locus2ko.get(l1,'')
    ko2 = _l2ko.get(l2,'')
    # if ko1==ko2 and ko1!='':
    #     same_c+=1
    # elif ko1!='' or ko2!='':
    #     print(l1,ko1,l2,ko2)
    if ko1!='':
        new_l2ko[l2] = ko1
ko2l = defaultdict(list)
for l,ko in new_l2ko.items():
    ko2l[ko].append(l)
ko2l = {k:','.join(sorted(v)) for k,v in ko2l.items()}
for k in ko2l:
    kegg_df.loc['GNM000014065',k] = ko2l[k]
for k in kegg_df.columns:
    if k not in ko2l:
        kegg_df.loc['GNM000014065',k] = np.nan


# 'sil'
result = REST.kegg_link("ko", "sil").read()
locus2ko = {}
for row in result.strip().split('\n'):
    l,k = row.split('\t')
    locus2ko[l.split(':')[-1]] = k.split(':')[-1]
_df = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/tmp/GNM000011965.tab',sep='\t',header=None)
_df.sort_values(2,ascending=False).groupby(0).head(1)
l2l = {}
for _,row in _df.iterrows():
    if row[2]>=99:
        l2l[row[0]] = row[1]
new_l2ko = {}
same_c = 0
for l1,l2 in l2l.items():
    ko1 = locus2ko.get(l1,'')
    ko2 = _l2ko.get(l2,'')
    # if ko1==ko2 and ko1!='':
    #     same_c+=1
    # elif ko1!='' or ko2!='':
    #     print(l1,ko1,l2,ko2)
    if ko1!='':
        new_l2ko[l2] = ko1
ko2l = defaultdict(list)
for l,ko in new_l2ko.items():
    ko2l[ko].append(l)
ko2l = {k:','.join(sorted(v)) for k,v in ko2l.items()}
for k in ko2l:
    kegg_df.loc['GNM000011965',k] = ko2l[k]
for k in kegg_df.columns:
    if k not in ko2l:
        kegg_df.loc['GNM000011965',k] = np.nan

kegg_df.to_csv('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_anno_Revised.tsv', sep='\t', index=1)

