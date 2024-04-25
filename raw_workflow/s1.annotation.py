"""
KEGG, ipr, pseudofinder annotations
1.
2.
todo:
1. make it accept and run multiple inputs
2. check gbk and infaa format and their consistence
3. enable logging
"""
from os.path import *
from dysfunctional_dectector.src.utilities.tk import check
#import modules
import logging
import click
from sys import exit
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
from subprocess import CalledProcessError
from multiprocessing import Process
from collections import defaultdict
#enable logging
import logging
logger = logging.getLogger('dysfunction_logger')
logger.setLevel(logging.DEBUG)
hd = logging.FileHandler('debuglog.log')
hd.setLevel(logging.DEBUG)
cons = logging.StreamHandler()
cons.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
hd.setFormatter(formatter)
cons.setFormatter(formatter)
logger.addHandler(hd)
logger.addHandler(cons)
## static setting
KOFAMSCAN_exe = '/mnt/storage3/yfdai/download/kofamscan/bin/exec_annotation'
KOFAMSCAN_profiles = '/mnt/home-db/pub/protein_db/kegg/v20230301/profiles'
KOFAMSCAN_ko_list = '/mnt/home-db/pub/protein_db/kegg/v20230301/ko_list'
pseudofinder_diamond_db = '/home-db/pub/protein_db/nr/v20230523/diamond_index/nr.dmnd'
diamond_path = '/home-db/pub/protein_db/nr/v20230523/diamond_index/diamond'
diamond_db = "/home-db/pub/protein_db/nr/v20230523/diamond_index/nr.dmnd"
pseudofinder_exe = 'python3 /mnt/storage3/yfdai/download/pseudofinder/pseudofinder-master/pseudofinder.py'
ipr_exe = '/mnt/storage3/yfdai/download/interproscan/interproscan-5.67-99.0/interproscan.sh'
num_cpu = 5

suffix_gbk = 'gbk'
suffix_protein = 'faa'


def gbk2faa(gbk):
    records = list(SeqIO.parse(gbk, format="genbank"))
    faa_list = []
    for contig in records:
        prefix = contig.name if contig.name else ""
        idx = 1
        for fea in contig.features:
            if fea.type != "CDS":
                continue
            if prefix:
                seq_name = f"{prefix}_" + "{:05d}".format(idx)
            else:
                seq_name = fea.qualifiers.get("locus_tag", [""])[0]
                if not seq_name:
                    seq_name = fea.qualifiers.get("protein_id", [""])[0]
            _seq = Seq(fea.qualifiers.get("translation",[''])[0])
            new_record = SeqRecord(
                seq=_seq,
                id=seq_name,
                name=fea.qualifiers.get("gene", [""])[0],
                description=fea.qualifiers.get("product",[""])[0],
                )
            faa_list.append(new_record)
            idx += 1
    return faa_list


def processing_IO(file_input,folder_output):
    in_files = defaultdict(dict)
    odir = folder_output
    temp_dir = odir + "/temp_file"
    os.makedirs(temp_dir,exist_ok=True)

    #check the consistence of the input file
    file_input = [realpath(_) for _ in file_input]

    gbk_seq = {}
    faa_seq = defaultdict(list)
    for file in file_input:
        extens = basename(file).split(".")[-1]
        genomename = basename(file).replace('.'+extens,'')
        if extens == suffix_gbk:
            gbk = file
            gbk_seq[genomename] = gbk2faa(gbk)
            if len(gbk_seq[genomename])==0:
                logger.debug(f"Genbank file doesn't have protein sequences.")
            in_files[genomename]['gbk'] = gbk
        elif extens == suffix_protein:
            infaa = file
            records = list(SeqIO.parse(file, "fasta"))
            faa_seq[genomename] = [str(_.seq) for _ in records]
        else:
            logging.error(f"The format of {file_input} is incorrect")
            print(f"The format of {file_input} is incorrect")
            exit()
    if set([str(_.seq) for _ in gbk_seq[genomename]]) == set(faa_seq[genomename]):
        in_files[genomename]['faa'] = infaa
        logger.debug("The faa and gbk file are consistent. Input faa file will be used.")
    else:
        filename = genomename+ "_gbk.faa"
        faafgbk = os.path.join(temp_dir,filename)
        if not exists(os.path.dirname(infaa)):
            os.makedirs(os.path.dirname(infaa))
        with open(faafgbk, "w") as f1:
            SeqIO.write(gbk_seq[genomename], f1, "fasta")
        logger.debug(f"Genbank file have been converted to faa file (see {filename})")
        in_files[genomename]['faa'] = faafgbk
        logger.debug("The protein and gbk file are not consistent. The faa file converted from input gbk file will be used.")

    return in_files,odir,temp_dir


def run_annotate(in_files,odir,temp_dir,dry_run=False):
    for genomename,info in in_files.items():
        ingbk = info['gbk']
        infaa = info['faa']
        ## OUTPUT
        kegg_oname = os.path.join(odir,'KOFAMSCAN',genomename+'.kofamout')
        pseudo_oname = os.path.join(odir,'pseudofinder',genomename)
        finalname = os.path.join(pseudo_oname,f'{genomename}_nr_pseudos.gff')
        ipr_oname = os.path.join(odir,'ipr',genomename)
        ######### Run kofamscan
        ogdir = os.path.join(odir,'KOFAMSCAN')
        os.makedirs(ogdir,exist_ok=True)
        kofams = f"{KOFAMSCAN_exe} -p {KOFAMSCAN_profiles} -k {KOFAMSCAN_ko_list} --tmp-dir {temp_dir}/{genomename}_kofamscan --cpu {num_cpu} -o {kegg_oname} -f mapper-one-line --no-report-unannotated {infaa} "
        kofams += f' && rm -rf {temp_dir}/{genomename}_kofamscan'
        cmds = check(kegg_oname,kofams,'kofamscan',dry_run=dry_run)
        logger.debug("Start KoFamscan")
        try:
            if cmds:
                Kofampros = Process(target=os.system, args=tuple(cmds))
                Kofampros.start()
                Kofampros.join()
            logger.debug("The KofanScan has done")
        except CalledProcessError as err:
            logger.error(f"The KofamScan fail to finish due to {err}")
            exit()
        logger.debug("Files for interproscan and pseudofinder are ready")
        ######### Run interpro
        iprcmd = f"""mkdir -p {ipr_oname} && export LD_LIBRARY_PATH='' && {ipr_exe} -i {infaa} -d {ipr_oname} -cpu {num_cpu} -iprlookup -appl CDD,Pfam"""
        cmds = check(ipr_oname,iprcmd,'interpro',dry_run=dry_run)
        logger.debug("Interproscan start")
        try:
            if cmds:
                iprpros = Process(target=os.system,args=tuple(cmds))
                iprpros.start()
                iprpros.join()
            logger.debug("The interproscan has done")
        except CalledProcessError as err:
            logger.error(f"The interproscan fail to finish due to {err}")
            exit()
        ######### Run pseudofinder
        pseudocmd = f"mkdir -p {pseudo_oname} ; {pseudofinder_exe} annotate -db {diamond_db}  -g {ingbk} -t {num_cpu} -skpdb -op {pseudo_oname}/{genomename}_nr -di -l 0.8 --compliant --diamond_path {diamond_path}"
        cmds = check(finalname,pseudocmd,'pseudofinder',dry_run=dry_run)
        logger.debug("Pseudofinder start")
        try:
            if cmds:
                psdpros = Process(target=os.system,args=tuple(cmds))
                psdpros.start()
                psdpros.join()
            logger.debug("The pseudofinder has done")
        except CalledProcessError as err:
            logger.error(f"The pseudofinder fail to finish due to {err}")
            exit()


# parse args
@click.command()
@click.option('-fi',"--file_input",type = str, nargs = 2 ,required = True, prompt = "Enter the input file name", help = "Please input faa or gbk file you want to analyze")
@click.option('-o',"--folder_output",type = str, nargs = 1 ,required = True, prompt = "Enter the output folder name", help = "Please input path of folder you want to analyze")
@click.option("-d","dry_run",help="Generate command only.",default=False,required=False,is_flag=True,)
def cli(file_input,folder_output,dry_run):
    in_files,odir,temp_dir = processing_IO(file_input,realpath(folder_output))
    run_annotate(in_files,odir,temp_dir,dry_run)


######### RUN
if __name__ == '__main__':
    cli()

