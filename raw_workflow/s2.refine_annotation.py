"""
Using reference-based to refine the annotations

Download or prepare a table based on the KEGG genome database, they contain the full list of the presence/absence of a type strain.

It is better to be runned with mutliple genomes


"""
from glob import glob
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
import pandas as pd
from os.path import exists,realpath,join,dirname,basename
import os
from Bio import SeqIO
import click
import portion as P
import pandas as pd
from dysfunctional_dectector.src.utilities.tk import output_file,kegg_api,get_gid_from_locus,check,parse_gbk
from dysfunctional_dectector.bin.refine_ko_with_ref import prepare_abbrev2files
from dysfunctional_dectector.src.utilities.logging import logger
import warnings;warnings.filterwarnings('ignore')
from subprocess import CalledProcessError,check_call
from multiprocessing import Process
import multiprocessing as mp
## static setting
mergedpseudo_script = join(dirname(dirname(__file__)),'src','pseudofinder_api','merged_pseudo.py')
ref_based_refine_path = join(dirname(dirname(__file__)),'bin','refine_ko_with_ref.py')



v2values = {'intact':1,
            'no KEGG-annotated':0,
            'confident pseudo':0.2,
            'RE(intact)':1,
            'RE(not intact)':0.5,
            'not intact':0.5
            }
def run(cmd):
    check_call(cmd,shell=1)

KOFAMSCAN_ko_list = '/mnt/home-db/pub/protein_db/kegg/v20230301/ko_list'
def get_all_kos(KOFAMSCAN_ko_list):
    all_kos = [_.split('\t')[0] for _ in open(KOFAMSCAN_ko_list).readlines()[1:]]
    return all_kos

# def fromKOtoIPR(kegg_df,ipr_df):
#     # read KO and find it ipr-based signature
#     locus2other_analysis2ID = defaultdict(dict)
#     ko2others = defaultdict(dict)
#     ko2ipr = {}
#     ko2locus = defaultdict(list)
#     for gid,row in kegg_df.iterrows():
#         for ko,locus_l in row.to_dict().items():
#             for locus in str(locus_l).split(','):
#                 if str(locus_l)!='nan':
#                     ko2locus[ko].append(locus)
#     locus2ko = {}
#     for ko,locus_l in ko2locus.items():
#         for locus in locus_l:
#             locus2ko[locus] = ko
#     if 'Sequence MD5 digest' in ipr_df.columns:
#         ipr_df.pop('Sequence MD5 digest')
#     subdf = ipr_df.loc[ipr_df['Protein accession'].isin(locus2ko),:]
#     subdf.loc[:,'match length'] = subdf['Stop location'] - subdf['Start location']
#     subdf.loc[:,'q cov'] = subdf['match length'].abs()/subdf['Sequence length']*100
#     subdf = subdf.loc[subdf['q cov']>=70,:]
#     for i,row in subdf.iterrows():
#         l = row['Protein accession']
#         locus2other_analysis2ID[l][row['Analysis']] = row['Signature accession']
#         ko = locus2ko[l]
#         if ko in ko2ipr:
#             continue
#         if row['InterPro accession'] != '-':
#             ko2ipr[ko] = row['InterPro accession']
#         else:
#             ko2others[ko][row['Analysis']] = row['Signature accession']
#     return locus2ko,ko2ipr, ko2others, locus2other_analysis2ID

# def get_neighbouring_profile(locus,dis=4):
#     gid,idx = get_gid_from_locus(locus),locus.split('_')[-1]
#     num_idx = len(idx)
#     neighbouring_locus = [f"{gid}_" + str.zfill(str(i),num_idx)
#                           for i in range(int(idx)-dis,int(idx)+dis)]
#     return neighbouring_locus

# def compare_two_db(db1,db2):
#     shared_db = set(db1).intersection(set(db2))
#     if len(shared_db)==0:
#         return 0
#     same_db = [k for k in shared_db if db1[k]==db2[k]]
#     return len(same_db)/len(shared_db)

# def compare_neighbouring_profile(batch_locus1,batch_locus2,l2profile,
#                                  based_db='CDD',return_same=False):
#     cdd_1 = [l2profile.get(l,{}).get(based_db,'') for l in batch_locus1]
#     cdd_2 = [l2profile.get(l,{}).get(based_db,'') for l in batch_locus2]
#     cdd_1 = [_ for _ in cdd_1 if _]
#     cdd_2 = [_ for _ in cdd_2 if _]
#     num_cdd = set(cdd_1).union(set(cdd_2))
#     same_cdd = set(cdd_1).intersection(set(cdd_2))
#     if return_same:
#         return [l for l in batch_locus1 if l2profile.get(l,{}).get(based_db,'') in same_cdd] + [l for l in batch_locus2 if l2profile.get(l,{}).get(based_db,'') in same_cdd]
#     if len(num_cdd)==0:
#         return 0
#     return len(same_cdd)/len(num_cdd)

def parse_alignment_(otab,contig,start,end):
    ## todo: speed it
    if os.path.getsize(otab)==0:
        return []
    try:
        d = pd.read_csv(otab,sep='\t',header=None)
    except:
        logger.debug(f'wrong file {otab}')
        exit()
    found_locus = []
    d = d.loc[d[1]==contig,:]
    d.loc[:,'start'] = d[[8,9]].min(1)
    d.loc[:,'end'] = d[[8,9]].max(1)
    subd = d.loc[~((d.end<=start-12000)|(d.start>=end+12000)),:]
    # smaller offset(8000) likely cause the candidate gene are too few
    for i,row in subd.iterrows():
        s,e = row[['start','end']]
        i = P.closed(s, e)
        any_overlap = P.closed(start,end).overlaps(i)
        if any_overlap:
            intersection_r = P.closed(start,end).intersection(i)
            length = intersection_r.upper-intersection_r.lower+1
            if int(row[2]) >= 60 and int(row[5])<=5 and length >=50:
                # pident and gapopen
                found_locus.append(row[0])
    return found_locus

def get_targetKO_from_alignment(arg):
    subj_ffn,ref_fna,oodir = arg
    sub_gid = subj_ffn.split('/')[-1].split('.')[0]
    gid = ref_fna.split('/')[-1].split('.')[0]
    otab = f"{oodir}/{gid}fna_{sub_gid}ffn.tab"
    cmd = f"blastn -task blastn -query {subj_ffn} -db {oodir}/{gid}fna -outfmt 6 -max_target_seqs 5 -evalue 1e-3 -qcov_hsp_perc 60 -num_threads 4  > {otab} "
    if not exists(otab) or os.path.getsize(otab)==0:
        check_call(cmd,shell=1)
    return otab


def batch_process(indir,gid,ref_locus2info,odir,kolog_file):
    subj_ffn = f"{indir}/link_files/{gid}.ffn"
    found_locus = []
    for ref_locus,info in ref_locus2info.items():
        ref_gid = get_gid_from_locus(ref_locus)
        contig,start,end = [info[_] for _ in ['contig','start','end']]
        #logger.debug(f"{contig} {start} {end}")
        ref_fna = f"{indir}/link_files/{ref_gid}.fna"
        otab = get_targetKO_from_alignment((subj_ffn,ref_fna,odir))
        found_locus = parse_alignment_(otab,contig,start,end)
        #logger.debug(f"found {ko} in {gid} with {ref_gid}: {ref_locus}.......")
        if found_locus:
            # logger.debug(f"found !!! {found_locus} based on {ref_locus}")
            with open(kolog_file,'w') as f1:
                f1.write(f"{';'.join(found_locus)}\t{ref_locus}\n")
            return found_locus,ref_locus
    final_l = found_locus
    if len(final_l)==0:
        with open(kolog_file,'w') as f1:
            f1.write(f"{';'.join(found_locus)}\t\n")    
        return final_l,None

def _batch(arg):
    batch_process(*arg)
    
def refining_KOmatrix(kegg_df,
                      locus_is_pseudo,
                      genome_pos,
                      indir,
                      odir,
                      target_kos=get_all_kos(KOFAMSCAN_ko_list),
                      num_threads=8,
                      verbose=1):
    """
    Determining the status of the presence and absence of a specific gene using the following information:
    1. pseudogene
    2. any other similar gene using CDD

    """
    # noted: locus should be look like {GID}_{number}
    _locus2ko = {}
    confident_presence = defaultdict(dict)
    gid2ko2intact_l = defaultdict(dict)
    gid2ko2pseudo_l = defaultdict(lambda :defaultdict(list))
    ko2intact_l = defaultdict(list)
    gid2cases = defaultdict(dict)
    recase2_l = defaultdict(lambda :defaultdict(list))
    for ko in tqdm(target_kos,desc='Parse raw status: '):
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
                _locus2ko[l] = ko
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
    ko2intact_l = {k:list(set(v)) for k,v in ko2intact_l.items()}

    # locus2length = {}
    # for gid in kegg_df.index:
    #     subj_ffn = f"{indir}/link_files/{gid}.ffn"
    #     locus2length.update({_.id:len(_.seq) for _ in SeqIO.parse(subj_ffn,'fasta')})
    #logger.debug(f"Start use Interproscan annotations to refine kegg anntoations...")
    log_odir = f"{odir}/refine_log/"
    if not exists(log_odir):
        os.makedirs(log_odir)
    for gid in tqdm(kegg_df.index,desc='build blast db: '):
        ref_fna = f"{indir}/link_files/{gid}.fna"
        oodir = odir
        if not exists(f"{oodir}/{gid}fna.nhr"):
            os.system(f"makeblastdb -dbtype nucl -in {ref_fna} -out {oodir}/{gid}fna > /dev/null 2>&1 ; ")
    
    ko2gid2status_df = pd.DataFrame.from_dict(gid2cases)
    if verbose:
        ko2row = tqdm(ko2gid2status_df.iterrows(),
                 total=ko2gid2status_df.shape[0],
                 desc="# of KOs: ")
    else:
        ko2row = ko2gid2status_df.iterrows()
    cmds = []
    ko_idx = 0
    for ko,row in ko2row:
        logger.debug(f"Now is {ko_idx}/{len(ko2row)} KOs......")
        if (ko not in ko2intact_l) or len(ko2intact_l[ko]) == 0:
            continue
        # genomes with no found KO
        target_gids = set([gid for gid,v in row.to_dict().items() if v == 'no KEGG-annotated'])
        logger.debug(f"{len(target_gids)} genomes with no KEGG-annotated {ko}")
        # ref_locus = ko2intact_l[ko][0]
        # ref_neighbouring_locus = get_neighbouring_profile(ref_locus)
        # ref_db_v = locus2other_analysis2ID[ref_locus]
        ref_locus_intact = ko2intact_l[ko]
        ref_locus2info = genome_pos.loc[ref_locus_intact,['contig','start','end']].to_dict(orient='index')        
    #     for sub_gid in target_gids:
    #         subj_ffn = f"{indir}/link_files/{gid}.ffn"
    #         for ref_locus,info in ref_locus2info.items():
    #             ref_gid = get_gid_from_locus(ref_locus)
    #             otab = f"{oodir}/{gid}fna_{sub_gid}ffn.tab"
    #             cmd = f"blastn -task blastn -query {subj_ffn} -db {oodir}/{ref_gid}fna -outfmt 6 -max_target_seqs 5 -evalue 1e-3 -qcov_hsp_perc 60 -num_threads 4  > {otab} "   
    #             if not exists(otab) or os.path.getsize(otab)==0:
    #                 cmds.append(cmd)
    # cmds = list(set(cmds))
    # sbatch_all(cmds,thread_per_tasks=4,prefix_name='blast',batch_size=100)
        ### batch module
        
        args = []
        for gid in target_gids:
            kolog_file = f"{log_odir}/{gid}_{ko}.done"
            if not exists(kolog_file):
                args.append((indir,gid,ref_locus2info,odir,kolog_file))
        logger.debug(f"batch run {len(args)} tasks...")
        with mp.Pool(processes=num_threads) as tp:
            _ = list(tqdm(tp.imap(_batch,args), 
                          total=len(args),
                          desc='# tasks: '))   
        ### batch module
        logger.debug(f"parsing....")
        ko_idx += 1
        for gid in target_gids:
            # genomes without this ko
            kolog_file = f"{log_odir}/{gid}_{ko}.done"
            if not exists(kolog_file):
                final_l,ref_locus = batch_process(indir,gid,ref_locus2info,odir,kolog_file)
            else:
                final_l,ref_locus = open(kolog_file).read().strip('\n').split('\n')[0].split('\t')
                final_l = final_l.split(';')
            final_l = [_ for _ in final_l if _]
            final_l = [_ for _ in final_l if not _locus2ko.get(_) ]
            #print(final_l,gid)
            if len(final_l)==0:   
                continue
            recase2_l[gid][ko].extend([(_,ref_locus)  for _ in final_l ])
            if any([_ in locus_is_pseudo for _ in final_l]):
                cases = 'not intact'
                gid2ko2pseudo_l[gid][ko].extend(final_l)
            else:
                cases = 'intact'
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
    for ko, locus_l in tqdm(ko2pseudo_l.items(),desc="# of KOs with pseudo: "):
        candidate_l = ko2intact_l[ko]
        if len(candidate_l)==0:
            for locus in locus_l:
                pseudo2assess_result[locus] = ('likely pseudo',ko,(0,0),0)
            continue
        lens = [abs(e-s) 
                for s,e in genome_pos.loc[candidate_l,['start','end']].values]
        low_thres = max(lens)*len_threshold
        max_thres = max(lens)*len_threshold*2
        pseudogen_regions = genome_pos.loc[locus_l,'pseudogenized']
        sub_genomepos = genome_pos.loc[genome_pos['pseudogenized'].isin(pseudogen_regions)]
        for locus,row in sub_genomepos.loc[locus_l,:].iterrows():
            is_conf,ll = is_confident_pseudo(row,(low_thres,max_thres),sub_genomepos)
            if is_conf:
                # remaining locus is longer than 60% of the target gene
                pseudo2assess_result[locus] = ('confident pseudo',ko,(low_thres,max_thres),ll)
            else:
                pseudo2assess_result[locus] = ('likely pseudo',ko,(low_thres,max_thres),ll)
    return pseudo2assess_result


def parse_kofamout(inf):
    l2ko = {}
    for row in open(inf).read().strip().split('\n'):
        rows = row.split('\t')
        l2ko[rows[0]] = rows[1]
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

def read_from_kofamout(kofamout_tab,kofamout,gid=None):
    """
    row is genome
    column is kegg ID
    """
    if not exists(kofamout_tab):
        l2ko = parse_kofamout(kofamout)
        gid2kegg2locus_info = defaultdict(lambda :defaultdict(list))
        for locus,ko in l2ko.items():
            #for ko in ko_l.split(';'):
            if gid is None:
                gid = get_gid_from_locus(locus)
            gid2kegg2locus_info[gid][ko].append(locus)
        gid2kegg2locus_info = {genome:{ko:','.join(list(set(l_list)))
                                    for ko,l_list in _d.items()}
                                    for genome,_d in gid2kegg2locus_info.items()}
        kegg_df = pd.DataFrame.from_dict(gid2kegg2locus_info, orient='index')
        kegg_df.to_csv(kofamout_tab,sep='\t',index=1)
    else:
        kegg_df = pd.read_csv(kofamout_tab,sep='\t',index_col=0)
    return kegg_df


def processing_f(genome_pos,kegg_df,indir,odir,num_threads=8,kolist=[]):
    """
    Using alignment based method to reannotated missing genes
    """
    #genome_pos = 
    locus_is_pseudo = set(list(genome_pos.index[genome_pos['pseudogenized'].str.contains(':')]))
    # Likely wrong, sometimes (extremely rare) a locus can be annotated to multiple KOs.
    #locus2ko, ko2ipr, ko2others, locus2other_analysis2ID = fromKOtoIPR(kegg_df,ipr_df)
    if not kolist:
        kolist = get_all_kos(KOFAMSCAN_ko_list)
    ko2gid2status_df,gid2ko2pseudo_l,confident_presence,recase2_l = refining_KOmatrix(kegg_df,
                                                                                      locus_is_pseudo,
                                                                                      genome_pos,
                                                                                      indir,
                                                                                      odir,
                                                                                      target_kos=kolist,
                                                                                      num_threads=num_threads,
                                                                                      verbose=1)

    pseudo2assess_result = assess_confident_pseudo_multi(gid2ko2pseudo_l,confident_presence,locus_is_pseudo,genome_pos,len_threshold=0.6)

    refined_kegg_df = kegg_df.copy()
    for gid,ko2l in recase2_l.items():
        for ko,locus_list in ko2l.items():
            v = refined_kegg_df.loc[gid,ko]
            if str(v)=='nan':
                refined_kegg_df.loc[gid,ko] = ','.join([_[0] for _ in locus_list])

    copy_df = ko2gid2status_df.copy()
    bin_df = copy_df.applymap(lambda x:v2values[x])
    for locus,(status,ko,reason,r2) in pseudo2assess_result.items():
        if status == 'likely pseudo':
            status = 'not intact'
        copy_df.loc[ko,get_gid_from_locus(locus)] = status
    return copy_df,bin_df,refined_kegg_df,recase2_l


def gen_genome_pos(gbk,pseudo_file):
    locus2info = parse_gbk(gbk)
    locus2info = pd.DataFrame.from_dict(locus2info).T
    locus2info.columns = ['contig', 'start', 'end', 'strand']

    pseudo_df_nano = pd.read_csv(pseudo_file,sep='\t',index_col=0)

    contig2interval = defaultdict(list)
    for _, row in pseudo_df_nano.iterrows():
        contig2interval[row["contig"]].append((P.closed(row["start"], row["end"]), _ )  )
    sdf = locus2info
    pseudo = []
    for _, row in sdf.iterrows():
        cand = contig2interval[row["contig"]]
        s, e = row['start'], row['end']
        i = P.closed(s, e)
        any_overlap = [(c[0],c[1]) for c in cand if c[0].overlaps(i)]
        pos = []
        for inter,str_inter in any_overlap:
            intersection_i = i.intersection(inter)
            ratio = abs(intersection_i.upper-intersection_i.lower)/abs(e-s)*100
            if ratio<=10:
                continue
            else:
                pos.append(str_inter)
        if len(pos)!=0:
            pseudo.append(';'.join(pos))
        else:
            pseudo.append('NOT')   
    sdf.loc[:, 'pseudogenized'] = pseudo  
    return sdf

def refine_pseudofinder(odir,dry_run=False):
    ## out of iterations
    genomename = '*'
    pseudo_oname = os.path.join(odir,genomename)
    pseudofinder_finalname = os.path.join(pseudo_oname,f'{genomename}_nr_pseudos.gff')
    ingbk_pattern = os.path.join(pseudo_oname,f'{genomename}.gbk')
    ###
    
    ######### Run pseudofinder merged script
    logger.debug("Try to parse pseudofinder result and merged some splitted pseudogenes.")
    odir = dirname(pseudo_oname) # without genomename
    cmd1 = f"python {mergedpseudo_script} -i '{pseudofinder_finalname}' -o {odir} step1 "
    cmd2 = f"python {mergedpseudo_script} -i '{pseudofinder_finalname}' -o {odir} step2 -gp '{ingbk_pattern}' "
    cmd = ';'.join([cmd1,cmd2])
    logger.debug("run : "+ cmd)
    try:
        psdpros = Process(target=os.system,args=tuple([cmd]))
        psdpros.start()
        psdpros.join()
    except CalledProcessError as err:
        logger.error(f"The pseudofinder fail to finish due to {err}")
        exit()    
        
    for output_gff in glob(pseudofinder_finalname):
        genomename = output_gff.split('/')[-2]
        pseudo_oname = os.path.join(odir,genomename)
        pseudofinder_merged_file = os.path.join(pseudo_oname,f'pseudofinderALL_w_annot_MERGEDcontinuous.tsv')             
        pos_oname = os.path.join(dirname(odir),'pos',genomename+'_pos.tsv')
        if exists(pos_oname):
            continue

        ####
        cmd3 = f"python {mergedpseudo_script} -i {output_gff} -o {odir}/{genomename} -snr {odir}/subnr step3 -u {odir}/merged_for"
        cmd4 = f"python {mergedpseudo_script} -i {output_gff} -o {odir}/{genomename} merge "
        cmd = ';'.join([cmd3,cmd4])
        cmds = check(pseudofinder_merged_file,cmd,'Merged pseudofinder',dry_run=dry_run)
        try:
            if cmds:
                psdpros = Process(target=os.system,args=tuple(cmds))
                psdpros.start()
                psdpros.join()
            else:
                logger.debug("The pseudofinder has been merged")
        except CalledProcessError as err:
            logger.error(f"The pseudofinder fail to finish due to {err}")
            exit()
        ## generate genome_pos
        ingbk = os.path.join(odir,genomename,genomename+'.gbk')
        pos_df = gen_genome_pos(ingbk,pseudofinder_merged_file)
        os.makedirs(dirname(pos_oname),exist_ok=True)
        pos_df.to_csv(pos_oname,sep='\t',index=1)

def merged_multiple(indir,odir,genomelist=[],kolist=[],link_file=None,num_threads=8,refine_pseudo=False):
    logger.debug(f'processing the refinement of pseudofinder results...')
    if refine_pseudo:
        refine_pseudofinder(f'{indir}/pseudofinder/',)
    
    logger.debug(f"reading {indir} for {genomelist} ......")
    # input
    kofamouts = glob(f'{indir}/KOFAMSCAN/*.kofamout')
    kofamout_tabs = [kofamout.replace('.kofamout','_anno.tab')
                    for kofamout in kofamouts]
    ipr_outputs = glob(f'{indir}/ipr/*/*.faa.tsv')
    genome_pos_files = glob(f'{indir}/pos/*_pos.tsv')
    if genomelist:
        kofamouts = [_ for _ in kofamouts
                     if _.split('/')[-1].rsplit('.')[0] in genomelist]
        kofamout_tabs = [kofamout.replace('.kofamout','_anno.tab')
                        for kofamout in kofamouts]
        ipr_outputs = [_ for _ in ipr_outputs
                       if _.split('/')[-2] in genomelist]
        genome_pos_files = [_ for _ in genome_pos_files
                            if _.split('/')[-1].rsplit('_pos')[0] in genomelist]
        num_gs = len(kofamout_tabs)
    # output
    os.makedirs(odir,exist_ok=True)
    ko_infodf = f'{odir}/{num_gs}Genomes_raw_ko_info.tsv'
    pos_path = f'{odir}/{num_gs}Genomes_pos.tsv'
    
    refined_ko_infodf = f'{odir}/{num_gs}Genomes_refined_ko_info.tsv'
    refined_ko_bindf = f'{odir}/{num_gs}Genomes_refined_ko_bin.tsv'
    kegg_new_outpath = f'{odir}/{num_gs}Genomes_refined_kegg.tsv'
    kegg_new_outlog = f'{odir}/{num_gs}Genomes_refined_kegg.log'
    
    logger.debug(f"found {len(kofamouts)} kegg,ipr,pos files......")
    if exists(ko_infodf):
        kegg_df = pd.read_csv(ko_infodf,sep='\t',index_col=0,low_memory=False)
    else:
        kegg_dfs = []
        for kofamout_tab,kofamout in zip(kofamout_tabs,kofamouts):
            kegg_df = read_from_kofamout(kofamout_tab,kofamout)
            kegg_dfs.append(kegg_df)
        kegg_df = pd.concat(kegg_dfs,axis=0)
        kegg_df.to_csv(ko_infodf,sep='\t',index=1)
    
    # ipr_dfs = []
    # for ipr_output in ipr_outputs:
    #     ipr_df = pd.read_csv(ipr_output,sep='\t',index_col=None,low_memory=False,names=headers)
    #     ipr_dfs.append(ipr_df)
    # ipr_df = pd.concat(ipr_dfs,axis=0)
    if exists(pos_path):
        genome_pos = pd.read_csv(pos_path,sep='\t',index_col=0,low_memory=False)
    else:
        gposss = []
        for genome_pos_file in genome_pos_files:
            gpos = pd.read_csv(genome_pos_file,sep='\t',index_col=0)
            gposss.append(gpos)
        genome_pos = pd.concat(gposss,axis=0)
        genome_pos.to_csv(pos_path,sep='\t',index=1)
    #### after reading
    logger.debug(f"Done reading")

    if link_file is not None:
        ofile = f"{odir}/refine/{basename(ko_infodf).rpartition('.')[0]+'_refined.'+basename(ko_infodf).rpartition('.')[-1]}"
        if not exists(ofile):
            cmd = f"python {ref_based_refine_path} -i {link_file} -k {ko_infodf} -o {odir}/refine -ss"
            logger.debug(f"run {cmd}")
            check_call(cmd,shell=1)
        else:
            logger.debug(f"existed {ofile}.... continue")
        kegg_df = pd.read_csv(ofile,sep='\t',index_col=0,low_memory=False)
        logger.debug(f"processing {ofile}")
    oodir = f"{odir}/align"
    os.makedirs(oodir,exist_ok=True)
    if genomelist:
        kegg_df = kegg_df.reindex(genomelist)
    copy_df,bin_df,refined_kegg_df,recase2_l = processing_f(genome_pos,kegg_df,indir,oodir,num_threads=num_threads,kolist=kolist)

    output_file(refined_ko_infodf,copy_df)
    output_file(refined_ko_bindf,bin_df)
    output_file(kegg_new_outpath,refined_kegg_df)
    
    with open(kegg_new_outlog,'w') as f1:
        f1.write(f"genome\tKO\tFound locus\tcomparable intact locus\n")
        for gid,ko2i in recase2_l.items():
            for ko,i in ko2i.items():
                for found_l,ref_l in i:
                    f1.write(f"{gid}\t{ko}\t{found_l}\t{ref_l}\n")
                    
    logger.debug("Done refining the annotation.")
    logger.debug(f"All {len(kofamout_tabs)} genomes result has been output to {realpath(refined_ko_bindf)}")


# parse args
@click.group()
@click.option('--odir','-o',help="output directory")
@click.option("-lf", "--link_file", help="input file linking between KEGG abbrev and Genome ID. such as 'sit\tGCA_0011965.1' ",required=False,default=None)
@click.pass_context
def cli(ctx,odir,link_file):
    ctx.ensure_object(dict)
    ctx.obj['odir'] = odir
    ctx.obj['link_file'] = link_file

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


@cli.command()
@click.option('--indir','-i',help="Input directory which is also a output directory of s1.",required=False,default=None)
@click.option('-gl','--genomelist',help="Input file containing the genome IDs you want to use or a comma-separated genome ids. ",required=False,default='')
@click.option('-nt','--num_threads',type = int,help="number of threads.",required=False,default=8)
@click.option('-kl','--kolist',help="Input file containing the KO IDs you want to use or a comma-separated genome ids. ",required=False,default='')
@click.option("-rp","--refine_pseudo",help="refine_pseudo or not ",default=False,required=False,is_flag=True,)
@click.pass_context
def mlworkflow(ctx,indir,genomelist,num_threads,kolist,refine_pseudo):
    outputdir = ctx.obj['odir']
    link_file = ctx.obj['link_file']
    indir = indir
    if indir is None:
        indir = outputdir+'/s1out'
    genomelist = genomelist
    if exists(genomelist):
        genomelist = open(genomelist).read().strip().split('\n')
    else:
        genomelist = [_ for _ in genomelist.split(',') if _]
    kolist = kolist
    if exists(kolist):
        kolist = open(kolist).read().strip().split('\n')
    else:
        kolist = [_ for _ in kolist.split(',') if _]        
    merged_multiple(indir,outputdir,
                    refine_pseudo=refine_pseudo,
                    genomelist=genomelist,link_file=link_file,num_threads=num_threads,kolist=kolist)

######### RUN
if __name__ == '__main__':
    cli()






