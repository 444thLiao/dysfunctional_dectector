from os.path import exists,dirname
import os
import logging

def check(ofile,cmd,name,dry_run=False,LOGGER=logging.getLogger('dysfunction_logger')):
    executed_cmd=[]
    if exists(ofile):
        LOGGER.debug(f'{name} existed')
    else:
        if not dry_run and not exists(dirname(ofile)):
            os.makedirs(dirname(ofile))
        LOGGER.debug(f'{name} run following command: \n {cmd}')
    if (not dry_run) and (not exists(ofile)):
        executed_cmd.append(cmd)
    elif dry_run:
        LOGGER.debug(f" '{cmd}' would not run.")
    return executed_cmd
        

def parse_gbk(gbk_file):
    info = {}
    records = list(SeqIO.parse(gbk_file, 'genbank'))
    for contig in records:
        all_fea = [_ for _ in contig.features if 'locus_tag' in _.qualifiers]
        for fea in all_fea:
            locus_tag = fea.qualifiers['locus_tag'][0]
            start = fea.location.start.real
            end = fea.location.end.real
            strand = fea.location.strand
            contig_name = contig.id
            info[locus_tag] = [contig_name, start, end, strand]
    return info

def batch_iter(iter, batch_size):
    # generating batch according batch_size
    iter = list(iter)
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d : len(iter) + 1])
    return n_iter

def output_file(opath,df):
    if opath.endswith('.xlsx'):
        df.to_excel(opath)
    elif opath.endswith('.tsv') or opath.endswith('.tab'):
        df.to_csv(opath,sep='\t',index=1)
    elif opath.endswith('.csv'):
        df.to_csv(opath,sep=',',index=1)
    else:
        raise IOError(f'Unknown suffix for {opath}')
        