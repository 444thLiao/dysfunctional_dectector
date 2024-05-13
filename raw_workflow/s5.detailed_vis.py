

import itertools
import os
from collections import defaultdict
from os.path import *
import pandas as pd
import plotly.graph_objects as go
from Bio import SeqIO
from tqdm import tqdm
import json
from geneblocks import CommonBlocks,DiffBlocks
import matplotlib.pyplot as plt
plt.ioff()
pd.options.display.expand_frame_repr = False
import io
from glob import glob

import matplotlib.pyplot as plt
import plotly.express as px
from dna_features_viewer import GraphicFeature, GraphicRecord

import warnings;warnings.filterwarnings('ignore')

import portion as P
import numpy as np
from IPython.display import display,Image



gid2nano_faa = {
    faa.split("/")[-3]: faa
    for faa in glob(f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/*/09_prokka/*.faa") 
    if '_' not in faa.split('/')[-1]}

gid2nano_fna = {gid:realpath(faa).replace('.faa','.fna') 
                for gid,faa in gid2nano_faa.items()}

contig2seq = {}
for _,f in gid2nano_fna.items():
    contig2seq.update({_.id:_ for _ in SeqIO.parse(f,'fasta')})
    
merged_pseudo = pd.read_csv("/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/pseudofinder/pseudofinderALL_w_annot_MERGEDcontinuous.tsv",sep='\t',index_col=0)    

ko_info = pd.read_csv(
    '/home-user/thliao/db/protein_db/kegg/v20230301/ko_list', sep='\t', index_col=0)
ko2info = {ko: i for ko, i in ko_info['definition'].to_dict().items()}
ko2g = pd.read_csv('/mnt/maple/thliao/data/protein_db/kegg/ko_info.tab',
                   sep='\t', header=None, index_col=0)
ko2g = {ko.split(':')[-1]: str(v).split(';')[0].strip()
        for ko, v in ko2g[1].to_dict().items()}

g2chrom = defaultdict(list)
for row in open('/home-user/thliao/project/coral_ruegeria/nanopore_processing/canu_o/chromosome.txt').read().strip().split('\n'):
    genome,contig = row.split('\t')
    g2chrom[genome].append(contig)
    



def get_fna(s,e,c,f,gbk=False):
    full_seq = contig2seq[c]
    seq = full_seq[s:e]
    with open(f,'w') as f1:
        SeqIO.write([seq],f1,'fasta-2line')
    return seq

def get_records(start,end,contig,func_df,pseudo_df_nano=None,extra_fea=True,
               **kwargs):
    
    _s = func_df.loc[func_df.contig==contig,:]
    sub_df = _s.loc[((_s.start<=start) & (_s.end >=end)) | ((_s.end>start) & (_s.end <end)) | ((_s.start>=start) & (_s.start <=end)),:]
    features=[]
    for i,row in sub_df.iterrows():
        f = get_type(row,0,**kwargs)
        features.append(f)
    if pseudo_df_nano is not None:
        _s = pseudo_df_nano.loc[(pseudo_df_nano['contig']==contig),:]
        sub_df2 = _s.loc[((_s.start<=start) & (_s.end >=end)) | ((_s.end>start) & (_s.end <end)) | ((_s.start>=start) & (_s.start <=end)),:]
        if extra_fea:
            for i,row in sub_df2.iterrows():
                features.append(GraphicFeature(start=row['start'],
                                               end=row['end'],
                                               strand=0,
                                               color='#10d269',label='Pseudo'))
        record = GraphicRecord(sequence_length=end, features=features)
        cropped_record = record.crop((start, end))
        cropped_record.sequence_length = end-start
        return cropped_record,features,start,end,sub_df2,sub_df
    else:
        record = GraphicRecord(sequence_length=end, features=features)
        cropped_record = record.crop((start, end))
        cropped_record.sequence_length = end-start
        return cropped_record,features,start,end,sub_df

def draw_regions(gid,c=None, s=None,strand=None,e=None,ko=None,neighbour_bp=5000,title_size=15,
                 highlight_kos=[],not_show=False,**kwargs):
    
    show_records = []
    if type(c)==list:
        for _c,_s,_e,_gid in zip(c,s,e,gid):
            func_df = g2func_df[_gid]
            original_region,features,s1,e1,pseudo,subfunc = get_records(_s,_e,_c,func_df,merged_pseudo,**kwargs)
            
            p = 'Chr' if _c in g2chrom[_gid] else 'Plasmid'
            title = f"{_gid} {p} ({g2pop[_gid]})"
            show_records.append((original_region,title,subfunc))    
        
    elif c is not None:
        func_df = g2func_df[gid]
        original_region,features,s1,e1,pseudo,subfunc = get_records(s,e,c,func_df,merged_pseudo,**kwargs)
        p = 'Chr' if c in g2chrom[gid] else 'Plasmid'
        title = f"{gid} {p} ({g2pop[gid]})"
        show_records.append((original_region,title,subfunc))    
    
    else:
        func_df = g2func_df[gid]
        r = func_df.loc[func_df['KEGG ID']==ko,:]
        for ridx,row in r.iterrows():
            c1 = row['contig']
            s = row['start']
            e = row['end']
            strand = row['strand']
            s1 = s-neighbour_bp
            e1 = e+neighbour_bp
            original_region,features,s1,e1,pseudo,subfunc = get_records(s1,e1,c1,func_df,merged_pseudo,**kwargs)
            p = 'Chr' if c1 in g2chrom[gid] else 'Plasmid'
            title = f"{gid} {p} ({g2pop[gid]})"
            show_records.append((original_region,title,subfunc))
    ### START DRAWING
    if not_show: 
        return [_[-1] for _ in show_records]   
    if len(show_records)==0:
        return
    if len(show_records)!=1:
        height = 3+ 1*len(show_records)
    else:
        height = 2
    fig, ax_list = plt.subplots(len(show_records), 1, figsize=(10, height)) 
    if len(show_records)==1: 
        ax_list = [ax_list]
    for (original_region,title,subfunc),ax in zip(show_records,ax_list):
        if not highlight_kos:
            for fea in original_region.features:
                if fea.start==s and fea.end==e:
                    fea.linecolor='#ff0000'
                    fea.linewidth=2
        else:
            highlight_kos_pos = list([tuple(sorted(_)) 
                                      for _ in subfunc.loc[subfunc['KEGG ID'].isin(highlight_kos),['start','end']].values])

            for fea in original_region.features:
                pos = sorted([fea.start,fea.end])
                pos = tuple(pos)
                if pos in highlight_kos_pos:
                    fea.linecolor='#ff0000'
                    fea.linewidth=2
        _ = ax.set_title(title,fontdict ={'fontsize':title_size})  
        _ = original_region.plot(ax=ax)
        if strand==-1:
            ax.invert_xaxis()
    plt.tight_layout()
    return [_[-1] for _ in show_records],fig

def get_type(row,base_pos = 0,
             func_col = 'top c (COG)',
             description_col='COG name',
             id_col='COG ID',
             label_col='related genes (COG)'):
    s,e,strand = row['start'],row['end'],row['strand']
    s = s-base_pos
    e = e-base_pos
    strand = 0 if str(strand)=='nan' else strand
    fea = GraphicFeature(start=s, end=e, strand=strand)
    if str(row['locus tag'])=='nan':
        fea.color = '#ffffff'
        fea.label =  'Insert Sequence (IS)'
        return fea
    if 'transposase' in str(row[description_col]).lower():
        fea.color = '#f57c00'
        #fea.label =  'transposes'
    elif str(row[func_col]) in ['nan','NA'] :
        fea.color = '#bdbdbd'
        #fea.label =  'Unknown'
    else:
        fea.color = '#0078d7'
        fea.label = str(row[id_col]) if str(row[label_col]) =='nan' else str(row[label_col])
    return fea

g2func_df = {}
for gid in tqdm(gid2nano_fna):
    if gid in g2func_df:
        continue
    efile = f'/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/summarized_tables/{g2pop[gid]}_{gid}.xlsx'
    if not exists(efile):
        continue
    _df = pd.read_excel(efile)
    g2func_df[gid] = _df
    

genome_pos = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/87genomes_pos.tsv',sep='\t',index_col=0)
genome_pos.loc[:,'is pseudo'] = ['Yes' if ':' in str(_) else 'No' for _ in genome_pos['pseudogenized']]
genome_pos.loc[:,'genome'] = [_.split('_')[0] for _ in genome_pos['contig']]
locus2contig = genome_pos['contig'].to_dict()
locus_is_pseudo = set(list(genome_pos.index[genome_pos['pseudogenized'].str.contains(':')]))
g2pop = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1783Ruegeria_MCs.tsv',sep='\t',index_col=0)
g2pop = g2pop['MC'].to_dict()
d = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_anno_Revised.tsv'
kegg_df = pd.read_csv(d,sep='\t', index_col=0)
kegg_df = kegg_df.loc[:,~kegg_df.isna().all(0).sort_values()]



related_kos = 'K02169;K00647;K09458;K00059;K02372;K00208;K02170;K09789;K19560;K19561;K00652;K00833;K01935;K01012;K19563;K19562'.split(';')
c_l = []
s_l = []
e_l = []
gids = []
for l in ['GB10_03259'#,'GK11_01204'
         ]:
          #'B4_03890']:
#for l in ['AN11_04296', 'AN11_03549', 'AN11_01818']:
    genome = l.split('_')[0]  
    gids.append(genome)
    n = 13000
    c,s,e = genome_pos.loc[l,:].values[:3]
    mid = (s+e)/2
    s = mid - n/2
    e = mid + n/2
    c_l.append(c)
    s_l.append(s)
    e_l.append(e)
    #print(abs(e-s))
_s,f = draw_regions(gids,c=c_l,s=s_l,e=e_l,
                    title_size=25,
                  highlight_kos=related_kos,strand=-1)





############### show detailed changes of pseudo genes
pseudo_idx = 'H9_1:3483585-3484881'
pseudo_contig,(pseudo_start,pseudo_end) = pseudo_idx.split(':')[0],pseudo_idx.split(':')[1].split('-')
pseudo_start,pseudo_end = int(pseudo_start),int(pseudo_end)

match_idx = 'AG8_1:3137592-3138887'
match_contig,(match_start,match_end) = match_idx.split(':')[0],match_idx.split(':')[1].split('-')
match_start,match_end = int(match_start),int(match_end)
strand = -1
s1,e1,c1,s2,e2,c2,seq1,seq2 = get_seqs(pseudo_contig,pseudo_start,pseudo_end,
                                       match_contig,match_start,match_end,
                                       strand=strand,contig2seq=contig2seq,
                                       extend_bp=2000)
func_df1 = g2func_df[pseudo_idx.split('_')[0]]
func_df2 = g2func_df[match_idx.split('_')[0]]
fig, ax_list = plt.subplots(3, 1, figsize=(15, 5))
_ = fig.set_facecolor("white")
original_region,features,s1,e1,sub_df2,sub_df = get_records(s1,e1,c1,func_df1,pseudo_df_nano,label_col='KEGG ID')
repeat_to_region,features,s2,e2,sub_df2,sub_df = get_records(s2,e2,c2,func_df2,pseudo_df_nano,label_col='KEGG ID')
_ = original_region.plot(ax_list[0])
_ = repeat_to_region.plot(ax_list[2])     
if strand == -1:
    _ = ax_list[2].invert_xaxis()
diff_blocks = DiffBlocks.from_sequences(seq1, seq2)
_ = diff_blocks.plot(ax=ax_list[1],separate_axes=False) 
_ = ax_list[1].set_xticks([])