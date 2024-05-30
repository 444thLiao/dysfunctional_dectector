"""
Using reference-based to refine the annotations

Download or prepare a table based on the KEGG genome database, they contain the full list of the presence/absence of a type strain. 



todo:


"""
from glob import glob
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
import pandas as pd
from os.path import exists,realpath
import os
import click
from dysfunctional_dectector.src.utilities.tk import output_file,kegg_api
from dysfunctional_dectector.bin.refine_ko_with_ref import prepare_abbrev2files
from dysfunctional_dectector.src.utilities.logging import *


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

v2values = {'intact':1, 
            'no KEGG-annotated':0, 
            'confident pseudo':0.2,
            'RE(intact)':1, 
            'RE(not intact)':0.5,
            'not intact':0.5
            }

KOFAMSCAN_ko_list = '/mnt/home-db/pub/protein_db/kegg/v20230301/ko_list'
def get_all_kos(KOFAMSCAN_ko_list):    
    all_kos = [_.split('\t')[0] for _ in open(KOFAMSCAN_ko_list).readlines()[1:]]
    return all_kos

def refining_KOmatrix(kegg_df,
                      locus_is_pseudo,
                      locus2other_analysis2ID,
                      locus2ko,
                      interpro_threshold=0.8,
                      target_kos=get_all_kos(KOFAMSCAN_ko_list),
                          verbose=1):
    # noted: locus should be look like {GID}_{number}
    
    confident_presence = defaultdict(dict)
    gid2ko2intact_l = defaultdict(dict)
    gid2ko2pseudo_l = defaultdict(lambda :defaultdict(list))
    ko2intact_l = defaultdict(list)
    gid2cases = defaultdict(dict)
    recase2_l = defaultdict(lambda :defaultdict(list))
    for ko in target_kos:
        ## For ko not found in this kegg_df
        if ko not in kegg_df.columns:
            for gid in kegg_df.index:
                gid2cases[gid][ko] = 'no KEGG-annotated'
                continue
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
                cases = 'not intact'
                gid2ko2pseudo_l[gid][ko] = [_  for _ in all_l
                                            if _  in pseudo_l]
            else:
                cases = 'intact'
                gid2ko2intact_l[gid][ko] = [_  for _ in all_l
                                            if _ not in pseudo_l]
            gid2cases[gid][ko] = cases
    logger.debug(f"Start use Interproscan annotations to refine kegg anntoations...")
    ko2gid2status_df = pd.DataFrame.from_dict(gid2cases)
    if verbose:
        ko2row = tqdm(ko2gid2status_df.iterrows(),
                 total=ko2gid2status_df.shape[0])
    else:
        ko2row = ko2gid2status_df.iterrows()
        
    for ko,row in ko2row:
        if ko not in ko2intact_l or len(ko2intact_l[ko]) == 0:
            continue
        target_gids = set([gid for gid,v in row.to_dict().items() if v == 'no KEGG-annotated'])
        ref_locus = ko2intact_l[ko][0]
        found_count = defaultdict(int)
        for db,db_v in sorted(locus2other_analysis2ID[ref_locus].items()):
            found = {k:v
                    for k,v in locus2other_analysis2ID.items()
                    if v.get(db,'') ==db_v and k.split('_')[0] in target_gids} # be careful to this
            for k,v in found.items():
                found_count[k]+=1
        num_db = len(locus2other_analysis2ID[ref_locus])

        for gid in target_gids:
            found_l = [k for k,v in found_count.items()
                       if v >=int(num_db*interpro_threshold)]
            found_l = [l 
                       for l in found_l 
                       if locus2ko.get(l,'') in ['',ko]]
            # if found locus belongs to other KO, ignore it
            if len(found_l)==0:
                continue
            if all([_ in locus_is_pseudo for _ in found_l]):
                cases = 'not intact'
                gid2ko2pseudo_l[gid][ko].extend(found_l)
            else:
                cases = 'intact'
                recase2_l[gid][ko].extend([_  for _ in found_l 
                                           if _ not in  locus_is_pseudo and _.startswith(f"{gid}_")])
            ko2gid2status_df.loc[ko,gid] = f"RE({cases})"
    return ko2gid2status_df,gid2ko2pseudo_l,confident_presence,recase2_l


def is_confident_pseudo(row,thresholds,genome_pos):
    """
    determined whether it meet the edge of the contig
    """
    pseudo_region = row['pseudogenized']
    strand = row['strand']
    cover_regions = genome_pos.loc[genome_pos['pseudogenized']==pseudo_region,:]
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
                                  locus_is_pseudo,genome_pos,
                                  len_threshold=0.6):
    """
    Need results from multiple genomes to assess the confidence of pseudogenes.

    """
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
    
    logger.debug(f"Start assessing identified pseudogenes ...")
    pseudo2assess_result = {}
    for ko, locus_l in tqdm(ko2pseudo_l.items()):
        candidate_l = ko2intact_l[ko]
        if len(candidate_l)==0:
            for locus in locus_l:
                pseudo2assess_result[locus] = ('likely pseudo',ko,(0,0),0)
            continue
        lens = [abs(e-s) for s,e in genome_pos.loc[candidate_l,['start','end']].values]
        low_thres = max(lens)*len_threshold
        max_thres = max(lens)*len_threshold*2
        pseudogen_regions = genome_pos.loc[locus_l,'pseudogenized']
        sub_genomepos = genome_pos.loc[genome_pos['pseudogenized'].isin(pseudogen_regions)]
        for locus,row in genome_pos.loc[locus_l,:].iterrows():
            is_conf,ll = is_confident_pseudo(row,(low_thres,max_thres),sub_genomepos)
            if is_conf:
                # remaining locus is longer than 50% of the target gene
                pseudo2assess_result[locus] = ('confident pseudo',ko,(low_thres,max_thres),ll)
            else:
                pseudo2assess_result[locus] = ('likely pseudo',ko,(low_thres,max_thres),ll)
    return pseudo2assess_result

def fromKOtoIPR(kegg_df,ipr_df,gid):
    # read KO and find it ipr-based signature
    locus2other_analysis2ID = defaultdict(dict)
    ko2others = defaultdict(dict)
    ko2ipr = {}

    ko2locus = {k:v 
                for k,v in kegg_df.loc[gid,:].to_dict().items() 
                if str(v)!='nan' and k not in ko2ipr}
    locus2ko = {}
    for ko,l in ko2locus.items():
        for _ in l.split(','):
            locus2ko[_] = ko
    if 'Sequence MD5 digest' in ipr_df.columns:
        ipr_df.pop('Sequence MD5 digest')
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

def parse_kofamout(inf):
    l2ko = {}
    for row in open(inf).read().strip().split('\n'):
        rows = row.split('\t')
        l2ko[rows[0]] = ';'.join(sorted(rows[1:]))
    return l2ko

def get_name(name):
    names = kegg_api.lookfor_organism("ruegeria")
    logger.debug(f"Use [{name}] to search against KEGG GENOME database...")
    logger.debug(f"Found {len(names)} items...")
    logger.debug(f"Use the following genomes(top5 at most) as accessory genomes: ")
    logger.debug('\n'.join(names))

    abbrev2name = {n.split(' ')[1]:n[0].split(' ',2)[-1].rsplit(' ',1)[0] 
                   for n in names[:5]}
    return abbrev2name

def main(indir,odir,gid,accessory_name=None):
    # input
    kofamout = f'{indir}/KOFAMSCAN/{gid}.kofamout'
    kofamout_tab = kofamout.replace('.kofamout','_anno.tab')
    ipr_output = f'{indir}/ipr/{gid}/{gid}.faa.tsv'
    genome_pos_file = f'{indir}/pos/{gid}_pos.tsv'
    # output
    os.makedirs(odir,exist_ok=True)
    refined_ko_infodf = f'{odir}/{gid}_refined_ko_info.tsv'
    refined_ko_bindf = f'{odir}/{gid}_refined_ko_bin.tsv'

    if not exists(kofamout_tab):
        l2ko = parse_kofamout(kofamout)
        gid2kegg2locus_info = defaultdict(lambda :defaultdict(list))
        for locus,ko_l in l2ko.items():
            for ko in ko_l.split(';'):
                gid2kegg2locus_info[gid][ko].append(locus)
        gid2kegg2locus_info = {genome:{ko:','.join(list(set(l_list))) 
                                    for ko,l_list in _d.items()} 
                                    for genome,_d in gid2kegg2locus_info.items()}
        kegg_df = pd.DataFrame.from_dict(gid2kegg2locus_info, orient='index')
        kegg_df.to_csv(kofamout_tab,sep='\t',index=1)
    else:
        kegg_df = pd.read_csv(kofamout_tab,sep='\t',index_col=0)



    ipr_df = pd.read_csv(ipr_output,sep='\t',index_col=None,low_memory=False,names=headers)
    if kegg_df.shape[1] == 1:
        logger.debug("Refining annotation using genome itself only, it will not be able to distinguish confident/likely pseudogenes.")

    genome_pos  = pd.read_csv(genome_pos_file,sep='\t',index_col=0)
    locus_is_pseudo = set(list(genome_pos.index[genome_pos['pseudogenized'].str.contains(':')]))

    locus2ko = {locus: ko
                for ko, _d in kegg_df.to_dict().items()
                for genome, locus_list in _d.items()
                for locus in str(locus_list).split(',')}
    ko2ipr, ko2others, locus2other_analysis2ID = fromKOtoIPR(kegg_df,ipr_df,gid)
    
    ko2gid2status_df,gid2ko2pseudo_l,confident_presence,recase2_l = refining_KOmatrix(kegg_df,locus_is_pseudo,locus2other_analysis2ID,locus2ko,verbose=1)

    pseudo2assess_result = assess_confident_pseudo_multi(gid2ko2pseudo_l,confident_presence,locus_is_pseudo,genome_pos,len_threshold=0.6)

    bin_df = copy_df.applymap(lambda x:v2values[x])


    # if accessory_name is not None:
    #     abbrev2name = get_name(accessory_name)
    #     for abbrev,gname in abbrev2name.items():
    #         ref_p,locus2ko = prepare_abbrev2files(abbrev,
    #                                               odir+'/accessory/')
    #         for locus,ko in locus2ko.items():
    #             bin_df.loc[ko,gname] = 1
    #             copy_df.loc[ko,gid] = 'intact'

    output_file(refined_ko_infodf,copy_df)
    output_file(refined_ko_bindf,bin_df)
    logger.debug("Done refining the annotation.")
    logger.debug(f"{gid} result has been output to {realpath(refined_ko_bindf)}")
    
    
# parse args
@click.group()
@click.option('--odir','-o',help="output directory")
@click.pass_context
def cli(ctx,odir):
    ctx.ensure_object(dict)
    ctx.obj['odir'] = odir

@cli.command()
@click.option('--genome','-gid',help="Genome ID")
@click.option('--indir','-i',help="Input directory which is also a output directory of s1.")
@click.option('--addbytext','-add',help="Input a name for searching  accessory genomes and use them to correct the genome you want to annotated.",required=False,default=False)
@click.pass_context
def workflow(ctx,genome,indir,addbytext):
    outputdir = ctx.obj['odir']
    indir = indir
    genome = genome
    accessory_name = addbytext
    main(indir,outputdir,genome,accessory_name)


######### RUN
if __name__ == '__main__':
    cli()






