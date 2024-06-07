###This file contains all the static settings of the dysfunctional_detector.
###It makes sense to run this script before using the programme.
###Make sure the working dictionary is */dysfunctional_detector
import os 

def change_settings(head_of_setting,new_path,input):
    file_path = input
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()
    modified_lines = []
    for line in lines:
        if line.startswith(head_of_setting):
            modified_lines.append(f'{new_path}\n')
        else:
            modified_lines.append(line)
    with open(file_path, "w", encoding="utf-8") as file:
        file.writelines(modified_lines)


working_directory = os.getcwd()
##KOFAMSCAN settings
ke = "KOFAMSCAN_exe = '/mnt/storage3/yfdai/download/kofamscan/bin/exec_annotation'"
kp = "KOFAMSCAN_profiles = '/mnt/home-db/pub/protein_db/kegg/v20230301/profiles'"
kl = "KOFAMSCAN_ko_list = '/mnt/home-db/pub/protein_db/kegg/v20230301/ko_list'"
##Pseudofinder settings
pdb = "pseudofinder_diamond_db = '/home-db/pub/protein_db/nr/v20230523/diamond_index/nr.dmnd'"
dp = "diamond_path = '/home-db/pub/protein_db/nr/v20230523/diamond_index/diamond'"
dd = "diamond_db = '/home-db/pub/protein_db/nr/v20230523/diamond_index/nr.dmnd'"
pe = "pseudofinder_exe = 'python3 /mnt/storage3/yfdai/download/pseudofinder/pseudofinder-master/pseudofinder.py'"
##Interproscan settings
ie = "ipr_exe = '/mnt/storage3/yfdai/download/interproscan/interproscan-5.67-99.0/interproscan1.sh'"
nc = "num_cpu = 5"
mergedpseudo_script = f'{working_directory}/src/pseudofinder_api/merged_pseudo.py'
#S2
KOFAMSCAN_ko_list = '/mnt/home-db/pub/protein_db/kegg/v20230301/ko_list'
#S3
k2e = "ko2ec_file = '/mnt/home-db/pub/protein_db/kegg/v20230301/link/ko2ec'"
k2p = "ko2pathway_file = '/home-user/thliao/db/protein_db/kegg/v20230301/link/ko2pathway'"
k2m = "ko2module_file = '/home-user/thliao/db/protein_db/kegg/v20230301/link/ko2module'"
#script_path
s1_path=f"{working_directory}/raw_workflow/s1.annotation.py"
s2_path=f"{working_directory}/raw_workflow/s2.refine_annotation.py"
s3_path=f"{working_directory}/raw_workflow/s3.detector.py"
print(s1_path,s2_path,s3_path)
#set the config of the whole programme
static_info = [["KOFAMSCAN_exe =",ke,s1_path],
               ["KOFAMSCAN_profiles =",kp,s1_path],
               ["KOFAMSCAN_ko_list =",kl,s1_path],
               ["pseudofinder_diamond_db =",pdb,s1_path],
               ["diamond_path =",dp,s1_path],
               ["diamond_db =",dd,s1_path],
               ["pseudofinder_exe = ",pe,s1_path],
               ["ipr_exe =",ie,s1_path],
               ["num_cpu =",nc,s1_path],
               ["mergedpseudo_script =",f"mergedpseudo_script = {mergedpseudo_script}",s1_path],
               ["KOFAMSCAN_ko_list =",f"KOFAMSCAN_ko_list = {kl}",s2_path],
               ["ko2ec_file =",k2e,s3_path],
               ["ko2pathway_file =",k2p,s3_path],
               ["ko2module_file =",k2m,s3_path]]
for info in static_info:
    change_settings(info[0],info[1],info[2])
    
#run setup.sh to set environment variables
# os.system(f"bash {working_directory}/setup.sh")
print("Environment variables, conda environment and static setting of the script are ready.")
