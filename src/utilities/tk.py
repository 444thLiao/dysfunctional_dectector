from os.path import exists,dirname
import os


def check(ofile,cmd,name,dry_run=True):
    if exists(ofile):
        LOGGER.append(f'{name} existed')
    else:
        if not dry_run and not exists(dirname(ofile)):
            os.makedirs(dirname(ofile))
        LOGGER.append(f'{name} run {cmd}')
    if (not dry_run) and (not exists(ofile)):
        #print(name)
        executed_cmd.append(cmd)
        
        
        
        
        