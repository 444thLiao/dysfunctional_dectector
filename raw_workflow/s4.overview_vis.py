import sys
sys.path.insert(0,'/home-user/thliao/script/evol_tk')
from Bio import Phylo
import plotly.graph_objs as go
from ete3 import Tree
import io,os
from collections import defaultdict
from tqdm import tqdm
import pandas as pd
import numpy as np
import plotly
from glob import glob
from IPython.display import Image,display
import pandas as pd
from api_tools import tree_vis, read_tree, sort_tree
import  pandas as pd
g2pop = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1623Ruegeria_MCs.tsv',sep='\t',index_col=0)
g2pop = g2pop['MC'].to_dict()

import re
def get_kname(ko):
    if not ko.startswith('K'):
        return ko
    name = kegg_get(ko).read().split('\n')[1].split(' ',1)[-1].strip().split(',')
    return name[-1].strip()


def get_f1(module_based_df,mode='binary',
           no_tree=False,
           title_text=None,pseudo_df=None,
           width=1700,height=1500,yaxisfontsize=25,
           rename_y=lambda x:x,
           rename_x=lambda x:['' for _ in x],
           xtickfontsize = 25,xtickangle=0,):

    df = module_based_df.reindex(index = tre.get_leaf_names())
    vt = tree_vis(tre, leaves2top=False)
    fig = plotly.tools.make_subplots(
        rows=1, cols=6, shared_yaxes=True, horizontal_spacing=0.05 / 3
    )
    
    for gid in df.index:
        if no_tree:
            continue
        fig.add_trace(go.Bar(x=[1],y=[1], marker={'color':pop2color.get(g2pop.get(gid,'Unknown'),'#ffffff'),'line':{'width':0}},name=gid+' '+g2pop.get(gid,''),showlegend=False),1,2)

        fig.add_trace(go.Bar(x=[1],y=[1],marker={'color':cmap_cs.get(gid2cs.get(gid,''),'#ffffff'),'line':{'width':0}},name=gid+';'+gid2cs.get(gid,''),showlegend=False),1,3)

        fig.add_trace(go.Bar(x=[1],y=[1],marker={'color':cmap_l.get(gid2l.get(gid,''),'#ffffff'),
                                                'line':{'width':0}},
                             name=gid+';'+gid2l.get(gid,''),
                            showlegend=False),1,4)

        fig.add_trace(go.Bar(x=[1],y=[1],
                             marker={'color':cmap_c.get(gid2c.get(gid,''),'#ffffff'),
                                                'line':{'width':0}},
                             name=gid+';'+gid2c.get(gid,''),
                            showlegend=False),1,5)

    y = [vt.clade_y_pos[vt.get_clade(c)] for c in df.index]
    if pseudo_df is not None:
        heatmap = go.Heatmap(
            x=df.columns,
            y=[_-0.5 for _ in y],
            z=pseudo_df.values,
            zmin=0,zmax=1,
             colorscale=[[0.0, "#ffffff"],
                         
                         [1,'#da16ff'],],
            #reversescale=True,
            showscale=False,
            xgap=1,ygap=1,)
        fig.append_trace(heatmap, 1, 6)   
    
    if mode=='custom':
        cs = [[0,'#bdbdbd'],
              [0.2,'#0078d7'],
              [0.5, "#FF80AB"],
              [1,'#dd011a']]
        ss = False
        heatmap = go.Heatmap(
            x=df.columns,
            y=[_-0.5 for _ in y],
            z=df.values,
            zmin=0,zmax=1,
            colorscale=cs,
            showscale=ss,
            xgap=4,ygap=4,)
    elif mode=='binary':
        cs = [[0,'#bdbdbd'],
              
              [1,'#dd011a']]
        ss = False
        heatmap = go.Heatmap(
            x=df.columns,
            y=[_-0.5 for _ in y],
            z=df.values,
            colorscale=cs,
            showscale=ss,
        xgap=4,ygap=4,)
    else:
        cs="RdBu"
        ss = True
        colorbar=dict(
            lenmode="pixels",
            yanchor="top",
            y=1.2,
            x=1.1,
            len=130,
            tickfont_size=15,
            #ticksuffix='%',
            #zmin=0,zmax=100,
            ticks="outside",
            #dtick=25
        )
        heatmap = go.Heatmap(
            x=df.columns,
            y=[_-0.5 for _ in y],
            z=df.values,
            colorscale=cs,
            showscale=ss,
            reversescale=True,
            xgap=4,ygap=4,
            colorbar=colorbar
        )
    fig.append_trace(heatmap, 1, 6)
    datas1 = vt.get_plotly_data(
        yscale=1, y_shift=-0.25,fix_length=None)
    if not no_tree:
        fig.add_traces(datas1, rows=[1] * len(datas1), cols=[1] * len(datas1))

    default_xaxis = {"zeroline": False,
                     "showticklabels": False,
                     "showspikes": False,}
    def get_axis(in_dict):
        _d = {}
        _d.update(default_xaxis)
        _d.update(in_dict)
        return _d
    if title_text:
        fig.update_layout(title=title_text,
                          title_x=0.5,
                         )
    fig.update_layout(
        barmode='stack',
        margin_b=0,margin_l=10,
    title_font_size=10,
    width=width,
    height=height,
    paper_bgcolor="#FFFFFF",
    plot_bgcolor="#FFFFFF",
    xaxis=get_axis({"domain": [0, 0.1],"showgrid": False}),    
    xaxis2=get_axis({"domain": [0.1, 0.12]}),
    xaxis3=get_axis({"domain": [0.12, 0.14]}),
    xaxis4=get_axis({"domain": [0.14, 0.16]}),
    xaxis5=get_axis({"domain": [0.16, 0.18]}),
    xaxis6=get_axis({"domain": [0.20, 1],
                     "side": "top",
                     #'tickfont_size':25,
                    "ticktext": rename_x(df.columns),
                    "tickvals": df.columns,
                    "tickfont": {'size': xtickfontsize},     
                     "tickangle":xtickangle,
                     "showspikes":True,
                     'showticklabels':True}),
    yaxis=default_xaxis,
    yaxis2=default_xaxis,
    yaxis3=default_xaxis,
    yaxis4=default_xaxis,
    yaxis5=default_xaxis,    
    yaxis6={
        "zeroline": False,
        "showgrid": False,
        "side": "right",
        "ticktext": rename_y(df.index),
        "tickvals": [_-0.5 for _ in y],
        "tickfont": {'size': yaxisfontsize},
        "showticklabels": True,
        "showspikes": False,
        "range": [0, sorted(vt.clade_y_pos.values())[-1]],
    },)
    display(Image(fig.to_image()))
    return fig
    


m = 'M00165'

for m,name in [('M00844','L-Ala'),
        ('M00019','L-Val'),
        ('M00432','L-Leu'),
        ('M00570','L-Ile'),
        ('M00118','L-Glu'),
        ('M00844','L-Arg'),
        ('M00015','L-Pro'),
        ('M00026','L-His'),
        ('M00023','L-Trp'),
        ('M00040','L-Tyr'),
        ('M00024','L-Phe'),
        ('M00017','L-Met'),
        ('M00021','L-Cys'),
        ('M00020','L-Ser'),
        ('M00016','L-Lys'),
        ('M00013','beta-Ala')
                 ]:
    m = 'M00036'
    if m not in _m:
        continue
    related_kos = []
    for _ in _m[m].values():
        for k in _:
            for _k in re.findall("K\d+",str(k)):
                related_kos.append(_k)
    related_kos = [_ for _ in related_kos if _ not in ['K13541','K00246','K00247']]
    subdf = fullcurated_bindf.fillna(0).reindex(columns=related_kos)
    subdf = subdf.loc[:,~subdf.isna().all(0)]
    print(f"www.genome.jp/module/{m}+"+'+'.join(subdf.columns))
    f = get_f1(subdf,
               no_tree=True,
               rename_x=lambda x:[ko2name.get(_,'').split(',')[0].strip()+ f" ({_})" for _ in x],
               rename_y=lambda x:[[gid2name[_]] if _ in gid2name else cus_rename([_]) for _ in x],
               xtickangle=90,
               #title_text=module_information[target_m][0]+' ' + target_m,
               width=1000,height=700,yaxisfontsize=20,mode='custom')

    f = get_f1(subdf,
               no_tree=True,
               rename_x=lambda x:[ko2ec.get(_,'').split(',')[0].strip().split(';')[0]+ f" ({_})" for _ in x],
               rename_y=lambda x:[[gid2name[_]] if _ in gid2name else cus_rename([_]) for _ in x],
               xtickangle=90,
               #title_text=module_information[target_m][0]+' ' + target_m,
               width=1000,height=700,yaxisfontsize=20,mode='custom')
    break
