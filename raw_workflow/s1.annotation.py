
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


## dynamic input
infaa = ''
odir = './'
gbk = ''
dry_run = True


## static setting
KOFAMSCAN_exe = '/home-user/thliao/software/kofamscan/exec_annotation'
KOFAMSCAN_profiles = '/mnt/home-db/pub/protein_db/kegg/v20230301/profiles'
KOFAMSCAN_ko_list = '/mnt/home-db/pub/protein_db/kegg/v20230301/ko_list'
pseudofinder_diamond_db = '/home-db/pub/protein_db/nr/v20230523/diamond_index/nr.dmnd'
diamond_path = '/home-db/pub/protein_db/nr/v20230523/diamond_index/diamond'
pseudofinder_exe = 'python3 /home-user/thliao/software/pseudofinder/pseudofinder.py'
ipr_exe = '/home-user/thliao/software/interproscan-5.63-95.0/interproscan.sh'

cpu = 5
## OUTPUT 
genome = basename(infaa).rpartition('.')[0]
kegg_oname = join(odir,'KOFAMSCAN',genome+'.kofamout')
pseudo_oname = join(odir,'pseudofinder',genome)
finalname = join({pseudo_oname},f'{genome}_nr_pseudos.gff')
ipr_oname = join(odir,'ipr',genome)

################## DEFINE COMMANDS 

cmd = f"{KOFAMSCAN_exe} -p {KOFAMSCAN_profiles} -k {KOFAMSCAN_ko_list} --tmp-dir .{genome}_tmp -o {kegg_oname} -f mapper-one-line --no-report-unannotated {infaa} "
cmd += f' && rm -rf .{genome}_tmp'
check(kegg_oname,cmd,'pseudofinder',dry_run=dry_run)

cmd = f"mkdir -p {pseudo_oname} ; {pseudofinder_exe} annotate -db  -g {gbk} -t 5 -skpdb -op {pseudo_oname}/{genome}_nr -di -l 0.8 --compliant -dip {diamond_path}"
check(finalname,cmd,'pseudofinder',dry_run=dry_run)

odir = "/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/ipr"
cmd = f"""mkdir {ipr_oname} && export LD_LIBRARY_PATH='' &&  -i {infaa} -d {ipr_oname} -cpu 5 -iprlookup -appl CDD,Pfam"""
check(ipr_oname,cmd,'interpro')

######### RUN


