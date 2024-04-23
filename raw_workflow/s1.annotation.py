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
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
from subprocess import CalledProcessError
from multiprocessing import Process
## dynamic input
infaa = ''
odir = './'
gbk = ''
in_file=""
dry_run = True
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
cpu = 5
#命令行工具
@click.command()
@click.option('-fi',"--file_input",type = str, nargs = 2 ,required = True, prompt = "Enter the input file name", help = "Please input faa or gbk file you want to analyze")
@click.option('-o',"--folder_output",type = str, nargs = 1 ,required = True, prompt = "Enter the output folder name", help = "Please input path of folder you want to analyze")
def cli(file_input,folder_output):
    odir = folder_output
    temp_dir = odir + "/temp_file"
    os.makedirs(temp_dir,exist_ok=True)
                
#check the consistence of the input file
    gbk_seq=[]
    faa_seq=[]
    for file in file_input:
        extens = basename(file).split(".")[-1]
        if extens == "gbk":
            gbk = file
            for record in SeqIO.parse(gbk, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        aa_seq = feature.qualifiers.get("translation", [""])[0]
                        gbk_seq.append(aa_seq)
                    else:
                        pass
            gbk_seq.sort(key=len)
        elif extens == "faa":
            infaa = file
            with open(infaa, "r") as file:
                records = list(SeqIO.parse(file, "fasta"))
                for record in records:
                    sequence = record.seq
                    faa_seq.append(sequence)
            faa_seq.sort(key=len)
        else:
            logging.error(f"The format of {file_input} is incorrect")
            print(f"The format of {file_input} is incorrect")
            exit()
    if gbk_seq == faa_seq:
        in_file = infaa
        logger.debug("The faa and gbk file are consistent.Input faa file will be used.")
    else:
        filename = gbk.split("/")[-1].split(".")[0]+ "_gbk.faa"
        faafgbk = os.path.join(temp_dir,filename)
        records = list(SeqIO.parse(gbk, format="genbank"))
        faa_list = []
        for contig in records:
            prefix = contig.name if contig.name else ""
            idx = 1
            for fea in contig.features:
                if prefix:
                    seq_name = f"{prefix}_" + "{:05d}".format(idx)
                else:
                    seq_name = fea.qualifiers.get("locus_tag", [""])[0]
                    if not seq_name:
                        seq_name = fea.qualifiers.get("protein_id", [""])[0]
                if fea.type == "CDS":
                    _seq = Seq(fea.qualifiers.get("translation",[''])[0])
                    new_record = SeqRecord(
                        seq=_seq,
                        id=seq_name,
                        name=fea.qualifiers.get("gene", [""])[0],
                        description=fea.qualifiers.get("product",[""])[0],
                    )
                    faa_list.append(new_record)
                    idx += 1
        if not exists(os.path.dirname(infaa)):
            os.makedirs(os.path.dirname(infaa))
        with open(faafgbk, "w") as f1:
            SeqIO.write(faa_list, f1, "fasta")
        logger.debug("Gbk file have been converted to faa file.")
        in_file = faafgbk
        logger.debug("The aa and gbk file are not consistent.The faa file converted from input gbk file will be used.")
    logger.debug("Check for the consistence has done.")
## OUTPUT 
    genome = basename(in_file).rpartition('.')[0]
    kegg_oname = os.path.join(odir,'KOFAMSCAN',genome+'.kofamout')
    pseudo_oname = os.path.join(odir,'pseudofinder',genome)
    finalname = os.path.join(pseudo_oname,f'{genome}_nr_pseudos.gff')
    ipr_oname = os.path.join(odir,'ipr',genome)
######### Run kofamscan
    ogdir = os.path.join(odir,'KOFAMSCAN')
    os.makedirs(ogdir,exist_ok=True)
    kofams = f"{KOFAMSCAN_exe} -p {KOFAMSCAN_profiles} -k {KOFAMSCAN_ko_list} --tmp-dir .{genome}_tmp -o {kegg_oname} -f mapper-one-line --no-report-unannotated {in_file} "
    kofams += f' && rm -rf .{genome}_tmp'
    check(kegg_oname,kofams,'pseudofinder',dry_run=dry_run)
    logger.debug("Start KoFamscan")
    try:
        Kofampros = Process(target=os.system, args=(kofams,))
        Kofampros.start()
        Kofampros.join()
        logger.debug("The KofanScan has done")
    except CalledProcessError as err:
        logger.error("The KofamScan fail to finish.")
        exit()
######### Select file for interproscan
    file_inipr = f"{temp_dir}/{genome}.faa"
    id2p = []
    with open(kegg_oname, "r") as file:
        for line in file:
            module = []
            line = line.strip()
            module = line.split("\t")
            id2p.append(module[0])
    faa2i = []
    records = list(SeqIO.parse(in_file, format="fasta"))
    for record in records:
        sequence_id = record.id
        if sequence_id not in id2p:
            faa2i.append(record)           
    with open(file_inipr, "w") as i:
        SeqIO.write(faa2i, i,"fasta")
    logger.debug("Files for interproscan and pseudofinder are ready")
######### Run interpro
    iprcmd = f"""mkdir -p {ipr_oname} && export LD_LIBRARY_PATH='' && {ipr_exe} -i {file_inipr} -d {ipr_oname} -cpu 12 -iprlookup -appl CDD,Pfam"""
    check(ipr_oname,iprcmd,'interpro')
    logger.debug("Interproscan start")
    try:
        iprpros = Process(target=os.system,args=(iprcmd,))
        iprpros.start()
        iprpros.join()
        logger.debug("The interproscan has done")
    except CalledProcessError as err:
        logger.error("The interproscan fail to finish.")
        exit()
######### Run pseudofinder   
    pseudocmd = f"mkdir -p {pseudo_oname} ; {pseudofinder_exe} annotate -db {diamond_db}  -g {gbk} -t 12 -skpdb -op {pseudo_oname}/{genome}_nr -di -l 0.8 --compliant --diamond_path {diamond_path}"
    check(finalname,pseudocmd,'pseudofinder',dry_run=dry_run)
    logger.debug("Pseudofinder start")
    try:
        psdpros = Process(target=os.system,args=(pseudocmd,))
        psdpros.start()
        psdpros.join()
        logger.debug("The pseudofinder has done")
    except CalledProcessError as err:
        logger.error("The pseudofinder fail to finish.")
        exit()
######### RUN
if __name__ == '__main__':
    cli()

