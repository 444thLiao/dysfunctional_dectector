#import modules
from dysfunctional_dectector.src.utilities.tk import output_file
from dysfunctional_dectector.src.utilities.tk import  ko2gename
import sys
sys.path.insert(0,'/home-user/thliao/script/evol_tk')
from Bio import Phylo
import plotly.graph_objs as go
from ete3 import Tree
import io,os
from collections import defaultdict
from tqdm import tqdm
import numpy as np
import plotly
from glob import glob
# from IPython.display import Image,display
import pandas as pd
from api_tools import tree_vis, read_tree, sort_tree
import re
from dysfunctional_dectector.src.utilities.logging import *
import click

#static setting
g2pop = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1783Ruegeria_MCs.tsv',sep='\t',index_col=0)
g2pop = g2pop['MC'].to_dict()
mc2genomes = {c: [g for g, _c in g2pop.items() if _c == c]for genome, c in g2pop.items()} 
pop2color = {row.split(',')[2]:row.split(',')[1] for row in open('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/genome2MC_lt3_colorstrip.txt').read().split('\n') if len(row.split(','))==3}
pop2color['MC59.1'] = pop2color['MC59']
pop2color['MC59.2'] = pop2color['MC59']
g2l_file = '/mnt/ivy/thliao/project/coral_ruegeria/data_processing/phylogeny/itol/genome2location_colorstrip.txt'
gid2l = {_.split(',')[0]:_.split(',')[2] for _ in open(g2l_file).read().split('\n') if len(_.split(','))==3 and _.split(',')[1].startswith('#') }
cmap_l = {_.split(',')[2]:_.split(',')[1] for _ in open(g2l_file).read().split('\n') if len(_.split(','))==3 and _.split(',')[1].startswith('#') }
g2c_file = '/mnt/ivy/thliao/project/coral_ruegeria/data_processing/phylogeny/itol/genome2compartment_colorstrip.txt'
gid2c = {_.split(',')[0]:_.split(',')[2] for _ in open(g2c_file).read().split('\n') if len(_.split(','))==3 and _.split(',')[1].startswith('#') }
cmap_c = {_.split(',')[2]:_.split(',')[1] for _ in open(g2c_file).read().split('\n')  if len(_.split(','))==3 and _.split(',')[1].startswith('#') }
g2cs_file = '/mnt/ivy/thliao/project/coral_ruegeria/data_processing/phylogeny/itol/coral_species_colorstrip.txt'
gid2cs = {_.split(',')[0]:_.split(',')[2] for _ in open(g2cs_file).read().split('\n') if len(_.split(','))==3 and _.split(',')[1].startswith('#') }

#check the value of the input file 
def format_check(folder_in):
    all_module = os.listdir(folder_in)
    all_path = [f"{in_folder}/{_}" for _ in all_module]
    for file_in in all_path:
        df_in = pd.read_csv(file_in, sep='\t')
        df_geno = df_in.iloc[:-1,1:]
        if df_in.size != df_in.shape[0]*df_in.shape[1]:
            logger.debug(f"The number of value in the record is wrong")
            exit()
        for index, row in df_geno.iterrows():
            for column, value in row.iteritems():
                if value in {"intact","no KEGG-annotated","confident pseudo","RE(intact)","RE(not intact)","not intact"}:
                    pass
                elif pd.isna(value):
                    pass
                else:
                    logger.debug(f"the value of the records is wrong")
                    exit() 
        for number in df_in.columns.tolist()[1:]:
            if number.startswith("K"):
                pass
            else:
                logger.debug(f"the ko annotation is incorrect")
                exit()
def get_exclusion(in_folder,exclude_collection):
    df_list=[]
    if exclude_collection:
        valid_name = []
        _m = set(exclude_collection)
        all_module = os.listdir(in_folder)
        all_name = list(set([_.rsplit('.', 1)[0].rsplit('_', 1)[0] for _ in all_module]))
        for name in all_name:
            if name in _m:
                pass
            else:
                valid_name.append(name)
            # for m,name in [('M00844','L-Ala'),
            #         ('M00019','L-Val'),
            #         ('M00432','L-Leu'),
            #         ('M00570','L-Ile'),
            #         ('M00118','L-Glu'),
            #         ('M00844','L-Arg'),
            #         ('M00015','L-Pro'),
            #         ('M00026','L-His'),
            #         ('M00023','L-Trp'),
            #         ('M00040','L-Tyr'),
            #         ('M00024','L-Phe'),
            #         ('M00017','L-Met'),
            #         ('M00021','L-Cys'),
            #         ('M00020','L-Ser'),
            #         ('M00016','L-Lys'),
            #         ('M00013','beta-Ala')
            #                 ]:
            #     m = 'M00036'
                # if m not in _m:
                #     continue
        valid_path = [f"{in_folder}/{_}_info.tsv" for _ in valid_name]
        # related_kos = []
        # for _ in _m[m].values():
        #     for k in _:
        #         for _k in re.findall("K\d+",str(k)):
        #             related_kos.append(_k)
        for file in valid_path:
            df = pd.read_csv(file, sep='\t')
            column_name = df.columns.tolist()
            if 'K13541' in column_name:
                df = df.drop('K13541', axis=1)      
            if 'K00246' in column_name:
                df = df.drop('K00246', axis=1) 
            if 'K00247' in column_name:
                df = df.drop('K00247', axis=1) 
            subdf = df.fillna(0)
            subdf = subdf.loc[:,~subdf.isna().all(0)]
            df_list.append(subdf)
        #logger.debug(f"www.genome.jp/module/{m}+"+'+'.join(subdf.columns))
    else:
        all_module = os.listdir(in_folder)
        for module in all_module:
            module_path = f"{in_folder}/{module}"
            df = pd.read_csv(module_path, sep='\t') 
            subdf = df.fillna(0)
            subdf = subdf.loc[:,~subdf.isna().all(0)]
            df_list.append(subdf)
    return df_list
###change the value in the given dataframe
def convert_df(df_in):
    v2values = {'intact':1, 
                'no KEGG-annotated':0, 
                'confident pseudo':0.2,
                'RE(intact)':0.5, 
                'RE(not intact)':0.5,
                'not intact':0.5
                }
    value_df = df_in.replace(v2values)
    return value_df
###get the tree for the graph
          
###get graph
def get_kname(ko):
    if not ko.startswith('K'):
        return ko
    name = kegg_get(ko).read().split('\n')[1].split(' ',1)[-1].strip().split(',')
    return name[-1].strip()

def get_ko2name(in_folder):
    kolist = []
    all_module = os.listdir(in_folder)
    all_path = [f"{in_folder}/{_}" for _ in all_module]
    for file in all_path:
        df = pd.read_csv(file, sep='\t')
        column_name = df.columns.tolist()
        kolist.extend(column_name)
    kolist=list(set(kolist))
    k2n = ko2gename(kolist)
    return k2n

def get_f1(module_based_df,
           mode='binary',
           tree = "",
           title_text=None,pseudo_df=None,
           width=1700,height=1500,yaxisfontsize=25,
           rename_y=lambda x:x,
           rename_x=lambda x:['' for _ in x],
           xtickfontsize = 25,
           xtickangle=0,
           grouping = True):
    if tree == "":
        df = module_based_df[:-1]
    else:
        df = module_based_df[:-1].reindex(index = tree.get_leaf_names())
        vt = tree_vis(tree, leaves2top=False)
    fig = plotly.tools.make_subplots(
        rows=1, cols=6, shared_yaxes=True, horizontal_spacing=0.05 / 3
    )
    
    for gid in df.index:
        if tree == "":
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
    if grouping == True:
        step_list = set(module_based_df[-1].value)
        col_num= df.shape[1]
        block_width = (width*0.8 - (col_num-1) * heatmap.xgap) /col_num
        for i in len(step_list):
            s_coord = module_based_df[module_based_df[-1]==f"step{i}"].idxmin(axis=0).astype(int)
            f_coord = module_based_df[module_based_df[-1]==f"step{i}"].idxmax(axis=0).astype(int)
            s = (s_coord-0.5)*(block_width+heatmap.xgap)
            f = (f_coord+0.5)*(block_width+heatmap.xgap)
            heatmap.add_shape(type="rect",x0=s, y0= 0, x1=f, y1= height, fillcolor="red", opacity=0.3, line={"width":heatmap['layout']['xaxis']['xgap']/4, "dash":"dash", "color":"red"})
        fig.append_trace(heatmap, 1, 6)
    datas1 = vt.get_plotly_data(
        yscale=1, y_shift=-0.25,fix_length=None)
    if tree != "":
        fig.add_traces(datas1, rows=[1] * len(datas1), cols=[1] * len(datas1))

###get the axis of the graph
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
def plot_graph(input_dataframe,tre,op):
    f = get_f1(input_dataframe,
                tree=tre,
                rename_x=lambda x:[ko2name.get(_,'').split(',')[0].strip()+ f" ({_})" for _ in x],
                rename_y=lambda x:[[gid2name[_]] if _ in gid2name else cus_rename([_]) for _ in x],
                xtickangle=90,
                #title_text=module_information[target_m][0]+' ' + target_m,
                width=1000,height=700,yaxisfontsize=20,mode='custom')
    outfile_name = basename(input_dataframe).split(".")[-1]
    output_path =  op + "/" +outfile_name +"_ko2name.html"
    f.write_html(output_path)
    f = get_f1(input_dataframe,
                tree=tre,
                rename_x=lambda x:[ko2ec.get(_,'').split(',')[0].strip().split(';')[0]+ f" ({_})" for _ in x],
                rename_y=lambda x:[[gid2name[_]] if _ in gid2name else cus_rename([_]) for _ in x],
                xtickangle=90,
                #title_text=module_information[target_m][0]+' ' + target_m,
                width=1000,height=700,yaxisfontsize=20,mode='custom')
    outfile_name = basename(input_dataframe).split(".")[-1]
    output_path =  op + "/" + outfile_name +"_ko2ec.html"
    f.write_html(output_path)

@click.command()
@click.option('--folder_input','-ni',type = str, nargs = 1 ,required = True, prompt = "Enter the input folder containing dataframes", help = "Please input the result of s3")
@click.option('--tree','-tr',type = str, nargs = 1 , default = "",prompt = "Tree of the genomes", help = "Please input the corresponding tree file ")
@click.argument('exclusion',type = str)
@click.option('--folder_output', "-o", type=str, nargs=1, required=True, prompt="Enter the output folder name", help="Please input path of folder you want to analyze")
def cli(folder_input,tree,exclusion,folder_output):
    format_check(folder_input)
    processed_df_list = get_exclusion(folder_input,exclusion)
    ko2name = get_ko2name(folder_input)
    for dataframe in processed_df_list:
        ready_df = convert_df(dataframe)
        plot_graph(ready_df,tree,folder_output)


######## RUN
if __name__ == '__main__':
    cli()
