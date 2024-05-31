
from collections import defaultdict
from os.path import *
import pandas as pd
import plotly.graph_objects as go
from Bio import SeqIO
pd.options.display.expand_frame_repr = False
import io
from glob import glob

import matplotlib.pyplot as plt
import plotly.express as px
from dna_features_viewer import GraphicFeature, GraphicRecord

import warnings;warnings.filterwarnings('ignore')
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
import matplotlib.pyplot as plt
plt.ioff()



genome_pos = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/87genomes_pos.tsv',sep='\t',index_col=0)

def transform(fea,featype='IS'):
    if featype == 'IS':
        fea.color = '#ffffff'
        fea.label =  'Insert Sequence (IS)'
    elif featype == 'transposase':
        fea.color = '#f57c00'
        fea.label =  'transposes'
    elif featype == 'COG':
        fea.color = '#0078d7'
    else:
        fea.color = '#bdbdbd'
    return fea

### 
def get_fea(s,e,strand,base_pos = 0,NAME='',defaultcolor='#0078d7',**kwargs):
    s = s-base_pos
    e = e-base_pos
    strand = 0 if str(strand)=='nan' else strand
    fea = GraphicFeature(start=s, end=e, strand=strand)
    fea.label = NAME
    fea.color = defaultcolor
    return fea

def renamed_plottitle(gid,contig):
    p = 'Chr' if contig in g2chrom[gid] else 'Plasmid'
    title = f"{gid} {p} ({g2pop[gid]})"
    return title

def get_records(start,end,contig,pseudo_df_nano=None,
                feaname_gen = None,
                **kwargs):
    contig_pos = genome_pos.loc[genome_pos.contig==contig,:]
    subregion_contig_pos = contig_pos.loc[((contig_pos.start<=start) & (contig_pos.end >=end)) | ((contig_pos.end>start) & (contig_pos.end <end)) | ((contig_pos.start>=start) & (contig_pos.start <=end)),:]

    features=[]
    for locus,row in subregion_contig_pos.iterrows():
        if feaname_gen is None:
            f = get_fea(row['start'],row['end'],row['strand'],NAME=locus)
        else:
            f = get_fea(row['start'],row['end'],row['strand'],NAME=feaname_gen(locus))
        features.append(f)
    if pseudo_df_nano is not None:
        pseudo_contig_pos = pseudo_df_nano.loc[(pseudo_df_nano['contig']==contig),:]
        subregion_pseudo_contig_pos = pseudo_contig_pos.loc[((pseudo_contig_pos.start<=start) & (pseudo_contig_pos.end >=end)) | ((pseudo_contig_pos.end>start) & (pseudo_contig_pos.end <end)) | ((pseudo_contig_pos.start>=start) & (pseudo_contig_pos.start <=end)),:]
        for _,row in subregion_pseudo_contig_pos.iterrows():
            features.append(GraphicFeature(start=row['start'],
                                            end=row['end'],
                                            strand=0,
                                            color='#10d269',label='Pseudo'))
    else:
        subregion_pseudo_contig_pos = None
    record = GraphicRecord(sequence_length=end, features=features)
    cropped_record = record.crop((start, end))
    cropped_record.sequence_length = end-start
    return cropped_record,features,start,end,subregion_pseudo_contig_pos,subregion_contig_pos

def draw_regions(info_l,title_size=15,
                 highlight_locus=[],not_show=False,**kwargs):
    ### parse
    show_records = []
    for _c,_s,_e,_gid in zip(info_l['contig'],info_l['start'],info_l['end'],info_l['genome']):
        original_region,features,s1,e1,pseudo,subfunc = get_records(_s,_e,_c,pseudo_df_nano=merged_pseudo,feaname_gen=kwargs['feaname_gen'])
        title = renamed_plottitle(_gid,_c)
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
    for (original_region,title,posdf),ax,strand in zip(show_records,ax_list,info_l['strand']):
        if not highlight_locus:
            for fea in original_region.features:
                if fea.start==s and fea.end==e:
                    fea.linecolor='#ff0000'
                    fea.linewidth=2
        else:
            highlight_locus_pos = [tuple(_) for _ in genome_pos.loc[highlight_locus,['start','end']].values]
            for fea in original_region.features:
                pos = sorted([fea.start,fea.end])
                pos = tuple(pos)
                if pos in highlight_locus_pos:
                    fea.linecolor='#ff0000'
                    fea.linewidth=2
        _ = ax.set_title(title,fontdict ={'fontsize':title_size})  
        _ = original_region.plot(ax=ax)
        if strand==-1:
            ax.invert_xaxis()
    plt.tight_layout()
    return [_[-1] for _ in show_records],fig


# self used
def parse_funcdf2typeinfo():
    pass


### showing confident absence



def show_two_locus(locus_list,l2cog,highlight_locus):
    """
    two locus pair in locus_list is to determine the direction to extend the region.
    """
    n = 10000
    info_l = {'contig':[],
                'start':[],
                'end':[],
                "strand":[],
                "genome":[]}
    for l in locus_list:
        if type(l)==str:
            genome = l.split('_')[0]  
            c,s,e = genome_pos.loc[l,:].values[:3]
            mid = (s+e)/2
            s = mid - n/2
            e = mid + n/2
            info_l['strand'].append(genome_pos.loc[l,'strand'])
        else:
            info_l['strand'].append(genome_pos.loc[l[0],'strand'])
            genome = l[0].split('_')[0]  
            r1,r2 = genome_pos.loc[l,:].values[:3]
            if r2[1]>=r1[1]:
                c,s,e = r1[:3]
                e +=n
            elif r2[1]<=r1[1]:
                c,s,e = r1[:3]
                s-=n
        info_l['genome'].append(genome)
        info_l['contig'].append(c)
        info_l['start'].append(s)
        info_l['end'].append(e)
        
    draw_regions(info_l,title_size=25,
                 highlight_locus=highlight_locus,feaname_gen=lambda x:l2cog.get(x,x))
def is_continuous(alist):
    vs = []
    for i,v in zip(alist,alist[1:]):
        vs.append( abs(v-i)==1)
    if all(vs):
        return True
    else:
        return False
    

def search_with_neighbouring(ko,gid,gid2status):
    ref = [k for k,v in gid2status.items() if 'intact' == v]
    ref_gid = ref[0]
    if not ref:
        return 
    else:
        ref_locus = kegg_df.loc[ref_gid,ko].split(',')[0]
        around_l = genome_pos.index.get_loc(ref_locus)
        around_l = genome_pos.index[around_l-5:around_l+5]
        d = g2func_df[ref_gid]
        l2cog = {k:v 
                for k,v in dict(zip(d['locus tag'],d['related genes (COG)'])).items()
                if str(v)!='nan' }
        l2cog[ref_locus] = ko
        subcog = '-'.join([l2cog.get(_,'') for _ in around_l if _ in l2cog and _ != ref_locus])

        d = g2func_df[gid]
        _l2cog = {k:v 
                for k,v in dict(zip(d['locus tag'],d['related genes (COG)'])).items()
                if str(v)!='nan' }
        l2cog.update(_l2cog)
        sublocus = [(locus,v) 
                    for locus,v in _l2cog.items() 
                    if v in subcog.split('-')]
        kmer = '-'.join(subcog.split('-')[:3])
        for grp in zip(*[sublocus[i:] for i in range(3)]):
            if not is_continuous([int(_[0].split('_')[-1]) for _ in grp]):
                continue
            if '-'.join([_[1] for _ in grp])==kmer:
                show_two_locus([(around_l[0],around_l[1]),
                                (grp[0][0],grp[1][0])],l2cog,highlight_locus=[ref_locus])
                
            if '-'.join([_[1] for _ in grp[::-1]])==kmer:
                show_two_locus([(around_l[0],around_l[1]),
                                (grp[-1][0],grp[-2][0])],l2cog,highlight_locus=[ref_locus])
                       
#### showing DNA blocks
for ko,col in fullcurated_df.loc[:,substats_df.columns].iteritems():
    gid2status = col.to_dict()
    for gid,s in gid2status.items():
        if gid not in substats_df.index or gid not in g2func_df:
            continue
        if s == 'no KEGG-annotated':
            search_with_neighbouring(ko,gid,gid2status)
        if s in ['no intact','RE(no intact)']:
            ref = [k for k,v in gid2status.items() if 'intact' in v]
            ref_locus = kegg_df.loc[ref[0],ko]
            sub_locus = kegg_df.loc[gid,ko].split(',')[0]
            d = g2func_df[gid]
            l2cog = {k:v 
                     for k,v in dict(zip(d['locus tag'],d['related genes (COG)'])).items()
                     if str(v)!='nan'}
            d = g2func_df[ref[0]]
            l2cog.update({k:v 
                          for k,v in dict(zip(d['locus tag'],d['related genes (COG)'])).items()
                          if str(v)!='nan'})
            
            print(ko,ko2g[ko],gid)
            locus_list = [ref_locus,sub_locus]
            show_two_locus(locus_list,l2cog=l2cog,highlight_locus=[ref_locus])
        if 'pseudo' in s:  
            ## confident pseudogene
            ref = [k for k,v in gid2status.items() if 'intact' in v]
            ref_locus = kegg_df.loc[ref[0],ko]
            sub_locus = kegg_df.loc[gid,ko].split(',')[0]
            d = g2func_df[gid]
            l2cog = {k:v 
                     for k,v in dict(zip(d['locus tag'],d['related genes (COG)'])).items()
                     if str(v)!='nan'}
            d = g2func_df[ref[0]]
            l2cog.update({k:v 
                          for k,v in dict(zip(d['locus tag'],d['related genes (COG)'])).items()
                          if str(v)!='nan'})
            
            print(ko,ko2g[ko],gid)
            locus_list = [ref_locus,sub_locus]
            show_two_locus(locus_list,l2cog=l2cog,highlight_locus=[ref_locus])
            #print(gid)



