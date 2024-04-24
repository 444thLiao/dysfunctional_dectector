from os.path import exists,dirname
import os
import logging

def check(ofile,cmd,name,dry_run=True,LOGGER=logging.getLogger('dysfunction_logger')):
    executed_cmd=[]
    if exists(ofile):
        LOGGER.append(f'{name} existed')
    else:
        if not dry_run and not exists(dirname(ofile)):
            os.makedirs(dirname(ofile))
        LOGGER.append(f'{name} run {cmd}')
    if (not dry_run) and (not exists(ofile)):
        #executed_cmd.append(cmd)
        executed_cmd.append(cmd)
    return executed_cmd
        
        
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
        