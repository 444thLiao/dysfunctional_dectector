"""

Use refined kegg annotation to calculate the completeness of each module and output some tables for further visualizations.


input: refined annotations
output: 
list of genes/locus/modules needed to be visualized


"""

from glob import glob
import sys
import pandas as pd
from os.path import exists,join,dirname
from bioservices.kegg import KEGG
kegg_api = KEGG()

brite_dir = '/mnt/home-db/pub/protein_db/kegg/v20230301/brite'
brite_file = '/mnt/home-db/pub/protein_db/kegg/v20230301/brite/ko00002.keg'

###

import pickle
from ete3 import Tree
import re
from collections import defaultdict
ko2ec_file = '/mnt/home-db/pub/protein_db/kegg/v20230301/link/ko2ec'
ko2pathway_file = '/home-user/thliao/db/protein_db/kegg/v20230301/link/ko2pathway'
ko2module_file = '/home-user/thliao/db/protein_db/kegg/v20230301/link/ko2module'
def parse_linkfile(infile):
    ko2others = defaultdict(list)
    for row in open(infile).read().strip().split('\n')[1:]:
        rows = row.split('\t')
        ko,others = rows[0],rows[1]
        ko,others = ko.split(':')[1],others.split(':')[1]
        ko2others[ko].append(others)
    ko2others = {ko:';'.join(others) for ko,others in ko2others.items()}

    return ko2others
ko2ec = parse_linkfile(ko2ec_file)
ko2pathway = parse_linkfile(ko2pathway_file)
ko2module = parse_linkfile(ko2module_file)





def get_koid(row):
    return re.findall(r'K\d+',row)

def parse_name(row):
    idx = row[0]
    row = row[1:]
    if get_koid(row):
        return get_koid(row)[0].strip()
    if '<b>' in row or idx in 'AB':
        name = row.split('<b>')[-1].split(' ',1)[-1].replace('</b>','').strip()
        return name
    elif idx == 'C':
        name = row.split(' ',1)[-1].strip()
        return name
    else:
        name = row.split(' ',1)[-1].strip()
        #print(row)
        mid = name.split(' ')[0]
        return mid


def load_pickle(dictionary_location):
    pickle_in = open(dictionary_location,"rb")
    dictionary = pickle.load(pickle_in)
    return dictionary

def read_module2ko():
    data_folder = join(dirname(__file__),'data')
    # Import all modules from dictionaries
    regular_modules = load_pickle(data_folder / "01.KEGG_Regular_Module_Information.pickle")
    bifurcation_modules = load_pickle(data_folder / "02.KEGG_Bifurcating_Module_Information.pickle")
    structural_modules = load_pickle(data_folder / "03.KEGG_Structural_Module_Information.pickle")
    return regular_modules,bifurcation_modules,structural_modules

level2depth = {'A':4,'B':3,'C':2,'D':1}
def parse_brite_file(brite_file):
    basic_tre = Tree()
    hier_infos = open(brite_file).read().strip().split('\n')
    br_name = hier_infos[0].split('\t')[-1]
    # if '+' in br_name:
    #     br_name = suppl_dict.get(br, 'MISSING')

    parent_node = ''
    for row in hier_infos:
        if row[0] not in 'ABCD':
            continue
        name = parse_name(row).strip()
        if not name:continue
        if row[0] == 'A':
            parent_node = basic_tre.add_child(name=name)
            parent_node.add_feature('level',row[0])
            continue
        ## trace back
        diff_num = level2depth[row[0]] - level2depth[parent_node.level] + 1
        if diff_num>=1:
            for _ in range(diff_num):
                parent_node = parent_node.up
        if row[0]!='D':
            parent_node = parent_node.add_child(name = name)
            parent_node.add_feature('level',row[0])
        else:
            _ = parent_node.add_child(name = name)
            _.add_feature('level',row[0])

    return basic_tre

module_brite = parse_brite_file(brite_file)


import numpy as np
def map_pseudo(x):
    if str(x)=='nan':
        return np.nan
    else:
        if all([l in locus_is_pseudo for l in x.split(',')]) :
            return np.nan
        else:
            return ','.join([l for l in x.split(',') if l not in locus_is_pseudo])
gnano_kegg = kegg_df.applymap(map_pseudo)
#####
gid2withpseudo_kos = {}
for gid,row in kegg_df.iterrows():
    gid2withpseudo_kos[gid] = list(row.index[~row.isna()])

gid2kos = {}
for gid,row in gnano_kegg.iterrows():
    gid2kos[gid] = list(row.index[~row.isna()])




def main(odir,gid):
    refined_ko_infodf_path = f'{odir}/{gid}_refined_ko_info.tsv'
    refined_ko_bindf_path = f'{odir}/{gid}_refined_ko_bin.tsv'
    ko_infodf = pd.read_csv(refined_ko_infodf_path,sep='\t',index_col=0)
    subko_df = ko_infodf.loc[ko_infodf[gid]!='no KEGG-annotated',:]

    mentioned_kos = list(subko_df.index)
    related_modules = list(set([m for ko in mentioned_kos for m in ko2module.get(ko,'').split(';')]))



######### RUN
if __name__ == '__main__':
    cli()
