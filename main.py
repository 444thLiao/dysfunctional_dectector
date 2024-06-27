"""
Apart from the workflow that can annotate a single genome (using faa and gbk files).
It needs to add an api to input the names of accessory genomes and download them together for annotation. 
It should also be able to accept multiple genome files together.


Scenarios used:

1. input a single faa and gbk file. annotate it.
2. input a file containing several faa and gbk files. annotate and merge them.
    2.1 Some files may come from the KEGG database. Use bin/refine_ko_withref.py to refine it.

Optional scenario:
a. Enter a name (e.g. ruegeria). We will download the protein sequences of the top 5 from the kegg database. And we will merge them together and annotate them using all methods. And finally merge all these results. 


"""
import click
import pandas as pd
import os
from os.path import exists,join
import multiprocessing as mp
from subprocess import check_call
import glob
from collections import defaultdict
from tqdm import tqdm
from subprocess import CalledProcessError
from dysfunctional_dectector.src.utilities.tk import check,parse_gbk
from dysfunctional_dectector.src.utilities.logging import logger
from dysfunctional_dectector.src.utilities.genome_download import get_search_result,download_genome
# variable initialization

current_directory = os.path.dirname(os.path.abspath(__file__))
s1_path = f"{current_directory}/raw_workflow/s1.annotation.py"
s2_path = f"{current_directory}/raw_workflow/s2.refine_annotation.py"
s3_path = f"{current_directory}/raw_workflow/s3.detector.py"

###function definition
def get_input_info(infile,fi):
    input_info = defaultdict(dict)
    if exists(infile):
        df = pd.read_csv(infile, sep='\t')
        for index, row in df.iterrows():
            input_info[row['genome']]['faa'] = row['protein file']
            input_info[row['genome']]['gbk'] = row['gbk file']    
        any_abrev = df.loc[~df['abbrev'].isna(),:].shape[0]
        return input_info,any_abrev
    elif type(fi)==tuple and infile == "":
        #for faa_file,gbk_file in fi:
        faa_file = [_ for _ in fi if _.endswith(".faa")][0]
        gbk_file = [_ for _ in fi if not _.endswith(".faa")][0]
        genome_id = os.path.basename(faa_file).replace('.faa','')
        input_info[genome_id]['faa'] = faa_file
        input_info[genome_id]['gbk'] = gbk_file
        return input_info,False
    else:
        logger.error("Please input the metadata file or single faa and gbk file.")
        exit()


def run_command(cmd,module_name,dry_run=False):
    try:
        check_call(cmd,shell=1)
        logger.debug(f"{module_name} module has done")
    except CalledProcessError as err:
        logger.error(f"{module_name} module fails to finish due to {err}")
        exit()

def single_workflow(gbk,faa,out_folder,dry_run,addbytext):
    logger.debug(f"Now s1 module is running for {faa.split('/')[-1]}")
    odir = f"{out_folder}/s1out"
    genomename = gbk.split('/')[-1].rsplit('.')[0]
    kegg_oname = join(odir,'KOFAMSCAN',genomename+'.kofamout')
    pseudo_oname = join(odir,'pseudofinder',genomename)
    pseudofinder_finalname = join(pseudo_oname,f'{genomename}_nr_pseudos.gff')     
    #logger.debug(f"Identifying {kegg_oname} and {pseudofinder_finalname}")   
    # if exists(kegg_oname) and exists(pseudofinder_finalname):
    #     #logger.debug(f"s1 module [{genomename}] has finished.")
    #     return 
    # else:
    cmd = f"python3 {s1_path} -fi {faa} {gbk} -o {out_folder}/s1out "
    run_command(cmd,"s1",dry_run=dry_run)
    logger.debug("s1 module has finished.")
    # gbk_name = os.path.basename(gbk)
    # g_id = os.path.splitext(gbk_name)[0]
    # for filename in os.listdir(f"{out_folder}/s1out/ipr/{g_id}"):
    #     old_path = f"{out_folder}/s1out/ipr/{g_id}/{filename}"
    #     new_name = g_id + "." + filename.split(".")[-2] + "." + filename.split(".")[-1]
    #     new_path = f"{out_folder}/s1out/ipr/{g_id}/{new_name}"
    #     os.rename(old_path,new_path)   

def merged_multiple_workflow(out_folder,num_gs,link_file=None,num_threads=8):
    if link_file is None:
        cmd = f"python3 {s2_path} -o {out_folder}/s2out mlworkflow -i {out_folder}/s1out -nt {num_threads} -rp"  
    else:
        cmd = f"python3 {s2_path} -o {out_folder}/s2out -lf {link_file} mlworkflow -i {out_folder}/s1out -nt {num_threads} -rp"  
    logger.debug("running : " + cmd)
    run_command(cmd,"s2")
    logger.debug("s2 module has finished. Now s3 module is running.")
    cmd = f"python3 {s3_path} -o {out_folder}/s3out/{num_gs}Genomes workflow -gid {num_gs}Genomes -i {out_folder}/s2out"
    run_command(cmd,"s3")
    logger.debug("s3 module has finished.")

def combine_result(s3_out):
    ###combine result of s3
    module_name = []  
    if not os.path.exists(s3_out+"/combination"):
        os.makedirs(s3_out+"/combination") 
    folders = list(os.scandir(s3_out))
    for folder in folders:
        module_folder =list(os.scandir(folder))
        for fd in module_folder:
            for file in os.listdir(fd):
                module_name.append(file.split(".")[0])
    module_name = set(module_name)
    for name in module_name:
        file_list=glob.glob(s3_out + f"/*/module_tables/{name}*")
        info_dict = pd.DataFrame({})
        ###
        data = pd.read_csv(file_list[0], sep='\t')
        column_name = data.columns.tolist()
        data_dict = data.to_dict()
        step_df = pd.DataFrame.from_dict(data_dict).iloc[-1,:].T
        for module in file_list:
            data = pd.read_csv(module, sep='\t')
            data_dict = data.to_dict()
            kegg_df = pd.DataFrame.from_dict(data_dict).iloc[:-1,:]
            info_dict = pd.concat([info_dict, kegg_df])
        info_dict = info_dict.reindex(columns=column_name)
        info_dict = info_dict.append(step_df)
        info_dict.index.name = 'Genome'
        info_dict.to_csv(f"{s3_out}/combination/{name}_combination.tsv",sep='\t',index=0)


def run(x):
    return single_workflow(*x)
# parse args
@click.command()
@click.option('--infile','-i',help="An Metadata file as input. See example file in XXXX",required=False,default=str)
@click.option('-fi',"--file_input",type = str, nargs = 2 ,required = False, help = "Please input faa and gbk file you want to analyze",default=None)
@click.option('-o',"--folder_output",type = str, nargs = 1 ,required = True, help = "Please output path of folder you want to store the analysis results. ")
@click.option("-d","--dry_run",help="Generate command only.",default=False,required=False,is_flag=True,)
@click.option('-at','--addbytext',type = str,help="Input a name for searching  accessory genomes and use them to correct the genome you want to annotated.",required=False,default='')
@click.option('-nt','--num_threads',type = int,help="number of threads.",required=False,default=8)
def cli(infile,file_input,folder_output,dry_run,addbytext,num_threads):
    folder_name = ["s1out","s2out","s3out"]
    for folder in folder_name:
        folder_path = f"{folder_output}/{folder}"
        os.makedirs(folder_path, exist_ok=True)
            
    input,any_abrev = get_input_info(infile,file_input)
    args_collect = []
    with mp.Pool(processes=num_threads) as tp:
        for genome_id,info in tqdm(input.items(),desc='Genomes: '):
            faa_path, gbk_path = info['faa'],info['gbk']
            args_collect.append((gbk_path,faa_path,folder_output,dry_run,addbytext))
        r = list(tqdm(tp.imap(run, args_collect), 
                      total=len(args_collect),desc='# of Running tasks: '))

    if len(input)==1:
        combine_result(f"{folder_output}/s3out")
    else:
        num_gs = len(input)
        if any_abrev:
            merged_multiple_workflow(folder_output,num_gs,infile,num_threads=num_threads)
        else:
            merged_multiple_workflow(folder_output,num_gs,num_threads=num_threads)
    
    # logger.debug("The given file has been analyzed with s1,s2 and s3")
    # if addbytext:    
    #     logger.debug("now start accessory genomes annotation.")
    #     ##start accessory genome annotation
    #     search_result = get_search_result(addbytext)
    #     download_result=download_genome(search_result,folder_output)
    #     logger.debug("The accessory genome has been downloaded.")
    #     genome_input = get_input_info(download_result,"")
    #     for key,value in tqdm(genome_input.items(),desc='Genomes: '):
    #         process = multiprocessing.Process(target=single_workflow,args=(value,key,folder_output,dry_run,addbytext))
    #         process.start()
    #         process.join()
    #     combine_result(f"{folder_output}/s3out")    
    #     logger.debug("The accessory genomes annotation has done.")

if __name__ == '__main__':
    cli()



