
from collections import defaultdict
from os.path import *
import pandas as pd
import plotly.graph_objects as go
from Bio import SeqIO
pd.options.display.expand_frame_repr = False
import io
from glob import glob
from os.path import *
import os
import pandas as pd
from pygenomeviz import GenomeViz
import matplotlib.pyplot as plt
import plotly.express as px
from dna_features_viewer import GraphicFeature, GraphicRecord

import warnings;warnings.filterwarnings('ignore')
from IPython.display import display,Image


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



########### validating two locus via alignment .
from Bio.SeqFeature import SeqFeature, FeatureLocation
def add_pseudo_fea(contig,merged_pseudo):
    if merged_pseudo is None:
        return contig
    else:
        sub_pseudo = merged_pseudo.loc[merged_pseudo.contig==contig.id]
        for i,row in sub_pseudo.iterrows():
            fn = SeqFeature(FeatureLocation(row['start'],
                                        row['end']),
                        type="pseudo",
                        qualifiers=dict(locus_tag=['pseudo'],),
                        strand=0,)
            contig.features.append(fn)
        return contig
    
def get_align(genome_pos,
              ref_gid,sub_gid,
              ref_locus,
              subj_locus,dist_bp = 5000,merged_pseudo=None):

    gbk1 = f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/{ref_gid}/09_prokka/{ref_gid}.gbk"
    gbk2 = f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/canu_o/{sub_gid}/09_prokka/{sub_gid}.gbk"

    ref_contig,ref_start,ref_end = genome_pos.loc[ref_locus,['contig','start','end']]
    ref_start,ref_end = ref_start-dist_bp,ref_end+dist_bp

    subj_contig,subj_start,subj_end = genome_pos.loc[subj_locus,['contig','start','end']]
    subj_start,subj_end = subj_start-dist_bp,subj_end+dist_bp  
    
    
    records = {_.id:_ for _ in SeqIO.parse(gbk1,'genbank')}
    contig = records[ref_contig]
    contig = add_pseudo_fea(contig,merged_pseudo)
    ref_gbk = contig[ref_start:ref_end]
    
    f = [fea for fea in ref_gbk.features if fea.qualifiers.get('locus_tag',[''])[0]==ref_locus and fea.type=='CDS']
    if f[0].strand == -1:
        oid = ref_gbk.id
        ref_gbk = ref_gbk.reverse_complement()
        ref_gbk.id = oid
    with open('./ref.fasta','w') as f1:
        SeqIO.write(ref_gbk,f1,'fasta-2line')        
    records = {_.id:_ for _ in SeqIO.parse(gbk2,'genbank')}
    contig = records[subj_contig]
    contig = add_pseudo_fea(contig,merged_pseudo)
    subj_gbk = contig[subj_start:subj_end]
    
    f = [fea for fea in subj_gbk.features if fea.qualifiers.get('locus_tag',[''])[0]==subj_locus and fea.type=='CDS']
    if f[0].strand == -1:
        oid = subj_gbk.id
        subj_gbk = subj_gbk.reverse_complement()
        subj_gbk.id = oid
    with open('./subj.fasta','w') as f1:
        SeqIO.write(subj_gbk,f1,'fasta-2line')   

    cmd = f"blastn -task blastn -query ref.fasta -subject subj.fasta -outfmt 6 > ./tmp.tab"
    os.system(cmd)
    return ref_gbk,subj_gbk

def manual_check_validation(genome_pos,
                            ref_locus,
                            subj_locus,
                            merged_pseudo=None,
                            dist_bp = 5000,
                            title_text=''):
    # name: O6fna_L7nucl.tab
    ref_gid = ref_locus.split('_')[0]
    subj_gid = subj_locus.split('_')[0]

    highlighted_locus = subj_locus.split(',')
    if ',' in subj_locus:
        subj_locus = subj_locus.split(',')[0]
    ref_gbk,subj_gbk = get_align(genome_pos,ref_gid,subj_gid,ref_locus,subj_locus,dist_bp=dist_bp,merged_pseudo=merged_pseudo)


    f1 = [(f.location.start.real,
        f.location.end.real,
        f.location.strand) for f in ref_gbk.features if f.type=='CDS']
    f2 = [(f.location.start.real,
        f.location.end.real,
        f.location.strand) for f in subj_gbk.features if f.type=='CDS']

    p1 = [(f.location.start.real,
        f.location.end.real,
        f.location.strand) for f in ref_gbk.features if f.type=='pseudo']
    p2 = [(f.location.start.real,
        f.location.end.real,
        f.location.strand) for f in subj_gbk.features if f.type=='pseudo']
        
    genome_list = [
        dict(name=ref_gid, size=len(ref_gbk.seq), features=f1,pseudo=p1),
        dict(name=subj_gid, size=len(subj_gbk.seq), features=f2,pseudo=p2),
    ]
    g_feas = [ ref_gbk,subj_gbk]

    gv = GenomeViz(fig_width=10,link_track_ratio=0.5)
    
    for gidx,genome in enumerate(genome_list):
        name, size, features = genome["name"], genome["size"], genome["features"]
        track = gv.add_feature_track(name, size)
   
        for idx, feature in enumerate(features):
            _gbk = g_feas[gidx]
            feas = [_ for _ in _gbk.features if _.type=='CDS']
            fea = feas[idx]
            name = fea.qualifiers.get('gene',[''])[0]
            if not name:
                name = fea.qualifiers.get('locus_tag',[''])[0]
            start, end, strand = feature
            if fea.qualifiers.get('locus_tag',[''])[0] in [ref_locus]+highlighted_locus:
                track.add_feature(start, end, strand, 
                              label=name,
                            facecolor='red',linewidth=1,arrow_shaft_ratio=1.0,
                            )
            else:
                track.add_feature(start, end, strand, 
                            facecolor='#0078d7',
                            label=name,arrow_shaft_ratio=1.0,
                            )
    if getsize('./tmp.tab')==0:
        pass
    else:
        link_text = pd.read_csv('./tmp.tab',sep='\t',header=None,)
        for _,row in link_text.iterrows():
            name1 = row[0].split("_")[0]
            name2 = row[1].split("_")[0]
            s,e = row[[6,7]]
            _s,_e = row[[8,9]]
            gv.add_link((name1,s,e), (name2,_s,_e),curve=True,
                        v=float(row[2]), vmin=50)
            #gv.set_colorbar(["grey",], vmin=50)
    fig = gv.plotfig()
    track_list = gv.feature_tracks
    for track in track_list:
        name = track.name
        d = [_ for _ in genome_list if _['name']==name][0]
        for pseudo_fea in d['pseudo']:
            start, end, strand = pseudo_fea
            start,end = track.transform_coord(start),track.transform_coord(end)
            x, y = (start, end, end, start), (-1, -1, 1, 1)
            track.ax.fill(x, y, fc="#10d269", alpha=0.5, zorder=-1)         
            
    if title_text:
        ax = fig.get_axes()[0]
        ax.set_title(title_text,y=2,fontdict={'fontsize':19}) 
    #os.system(f"rm ./tmp.tab")
    return fig


