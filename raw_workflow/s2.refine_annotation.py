"""
Using reference-based to refine the annotations

Download or prepare a table based on the KEGG genome database, they contain the full list of the presence/absence of a type strain. 



todo:


"""




from subprocess import check_call
from glob import glob
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
import pandas as pd
import numpy as np
from os.path import exists
from dysfunctional_dectector.src.utilities.tk import output_file
## dynamic input
multiple_input = True  # a switch: then the kegg_output should be a file referred to a dataframe
kegg_output = ''
genome_pos_file = '' # from gbk
pseudogenes_output = ''
ipr_output_file = ''


# prepare a file look like this

# Genome ID\tkegg output file\t ipr output file\t gbk file\t pseudogene output file
# 

## static setting
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
all_kos = []

## OUTPUT 
refined_ko_bindf = ''
refined_ko_infodf = ''



    
def refining_KOmatrix(kegg_df,
                      locus_is_pseudo,
                      locus2other_analysis2ID,
                      locus2ko,
                      interpro_threshold=0.8,
                          verbose=1):
    # noted: locus should be look like {GID}_{number}
    
    confident_presence = defaultdict(dict)
    gid2ko2intact_l = defaultdict(dict)
    gid2ko2pseudo_l = defaultdict(lambda :defaultdict(list))
    ko2intact_l = defaultdict(list)
    gid2cases = defaultdict(dict)

    for ko in all_kos:
        ## For ko not found in this kegg_df
        if ko not in kegg_df.columns:
            for gid in kegg_df.index:
                gid2cases[gid][ko] = 'no KEGG-annotated'
                continue
        ## For ko in kegg_df, search it for each gid.
        ## first step of filtering
        gid2locus_str = kegg_df[ko].to_dict()
        pseudo_l = []
        intact_l = []
        for gid,locus in gid2locus_str.items():
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
        ko2row = tqdm(ko2gid2status_df.iterrows(),
                 total=ko2gid2status_df.shape[0])
    else:
        ko2row = ko2gid2status_df.iterrows()
        
    for ko,row in ko2row:
        for gid,v in row.to_dict().items():
            if v != 'no KEGG-annotated' or ko not in ko2intact_l:
                continue
            if len(ko2intact_l[ko]) == 0:
                continue
            i = ko2intact_l[ko][0]
            _a = gid+'_'
            # extract
            sub_l2all = {k:v 
                         for k,v in locus2other_analysis2ID.items()
                         if k.startswith(_a)}
            found_count = defaultdict(int)
            for db,db_v in sorted(locus2other_analysis2ID[i].items()):
                found = {k:v
                     for k,v in sub_l2all.items()
                     if v.get(db,'') ==db_v }
                for k,v in found.items():
                    found_count[k]+=1
            num_db = len(locus2other_analysis2ID[i])
            found_l = [k for k,v in found_count.items()
                       if v >=int(num_db*interpro_threshold)]
            found_l = [l 
                       for l in found_l 
                       if locus2ko.get(l,'') in ['',ko]]
            # if found locus belongs to other KO, ignore it
            if len(found_l)==0:
                continue
            if all([_ in locus_is_pseudo for _ in found_l]):
                cases = 'no intact'
                gid2ko2pseudo_l[gid][ko].extend(found_l)
            else:
                cases = 'intact'
            ko2gid2status_df.loc[ko,gid] = f"RE({cases})"
            # RE mean retrieved 
    return ko2gid2status_df,gid2ko2pseudo_l,confident_presence


def is_confident_pseudo(this_locus,thresholds,genome_pos):
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


def fromKOtoIPR(kegg_df,gid2ipr):
    # read KO and find it ipr-based signature
    locus2other_analysis2ID = defaultdict(dict)
    ko2others = defaultdict(dict)
    ko2ipr = {}
    for gid, ofile in tqdm(gid2ipr.items(),total=len(gid2ipr)):
        ko2locus = {k:v for k,v in kegg_df.loc[gid,:].to_dict().items() if str(v)!='nan' and k not in ko2ipr}
        locus2ko = {}
        for ko,l in ko2locus.items():
            for _ in l.split(','):
                locus2ko[_] = ko
        ipr_df = pd.read_csv(ofile,
                             sep='\t',index_col=None,low_memory=False,names=headers)
        ipr_df.pop('Sequence MD5 digest')
        # ipr_df.pop('Pathways annotations)
        subdf = ipr_df.loc[ipr_df['Protein accession'].isin(locus2ko),:]
        for i,row in subdf.iterrows():
            l = row['Protein accession']
            locus2other_analysis2ID[l][row['Analysis']] = row['Signature accession']
            ko = locus2ko[l]
            if ko in ko2ipr:
                continue
            if row['InterPro accession'] != '-':
                ko2ipr[ko] = row['InterPro accession']
            else:
                ko2others[ko][row['Analysis']] = row['Signature accession']
    return ko2ipr, ko2others, locus2other_analysis2ID

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
# kegg_df = pd.read_csv(    '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_anno_Revised.tsv', sep='\t', index_col=0)


# genome_pos = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/76genomes_pos.tsv',sep='\t',index_col=0)
# genome_pos.loc[:,'is pseudo'] = ['Yes' if ':' in str(_) else 'No' for _ in genome_pos['pseudogenized']]
# genome_pos.loc[:,'genome'] = [_.split('_')[0] for _ in genome_pos['contig']]



###### main ######
if multiple_input:
    kegg_df = pd.read_csv(kegg_output, sep='\t', index_col=0)
    genome_pos  = pd.read_csv(genome_pos_file,sep='\t',index_col=0)
else:
    sub_kegg_df = pd.read_csv(kegg_output, sep='\t', index_col=0)
    genome_pos  = pd.read_csv(genome_pos_file,sep='\t',index_col=0)
    ## todo: turn a single one into a dataframe
    
    
locus2contig = genome_pos['contig'].to_dict()
locus_is_pseudo = set(list(genome_pos.index[genome_pos['pseudogenized'].str.contains(':')]))
    
locus2ko = {locus: ko
            for ko, _d in kegg_df.to_dict().items()
            for genome, locus_list in _d.items()
            for locus in str(locus_list).split(',')}

ko2ipr, ko2others, locus2other_analysis2ID = fromKOtoIPR(kegg_df,gid2ipr)
ko2gid2status_df,gid2ko2pseudo_l,confident_presence = refining_KOmatrix(kegg_df,locus_is_pseudo,locus2other_analysis2ID,locus2ko,verbose=1)
pseudo2assess_result = assess_confident_pseudo_multi(gid2ko2pseudo_l,confident_presence,len_threshold=0.6)

copy_df = ko2gid2status_df.copy()
for locus,(_,ko,_,_) in pseudo2assess_result.items():
    gid = locus.split('_')[0]
    ori = copy_df.loc[ko,gid]
    if ori == 'RE(no intact)':
        copy_df.loc[ko,gid] = 'confident pseudo'
bin_df = copy_df.applymap(lambda x:v2values[x])

output_file(refined_ko_infodf,copy_df)
output_file(refined_ko_bindf,bin_df)







