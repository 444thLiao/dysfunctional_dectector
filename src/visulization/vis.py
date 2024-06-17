
from .DNAblocks import draw_regions

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
plt.ioff()
import warnings;warnings.filterwarnings('ignore')
from IPython.display import display,Image

##static_input
merged_pseudo_file = ''
genome_pos_file = ''
merged_pseudo = pd.read_csv(merged_pseudo_file,sep='\t',index_col=0)

genome_pos = pd.read_csv(genome_pos_file,sep='\t',index_col=0)
    
    
    
curated_info = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_curated_log.txt',sep='\t',index_col=0)

############
n = 10000
info_l = {'contig':[],
          'start':[],
          'end':[],
          "strand":[],
          "genome":[]}
ref_l = 'B4_03912'
sub_l = 'BK1_02246'
for l in [ref_l,
          sub_l
         ]:
    genome = l.split('_')[0]  
    c,s,e = genome_pos.loc[l,:].values[:3]
    mid = (s+e)/2
    s = mid - n/2
    e = mid + n/2
    info_l['genome'].append(genome)
    info_l['contig'].append(c)
    info_l['start'].append(s)
    info_l['end'].append(e)
    info_l['strand'].append(genome_pos.loc[l,'strand'])
    
    
_ = draw_regions(info_l,title_size=25,
                 highlight_locus=[ref_l,sub_l])

