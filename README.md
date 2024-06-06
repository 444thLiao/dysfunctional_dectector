
# Introduction of dysfunctional_detector:


Dysfunctional_detector is one python tool that annotates pseudogenes and maps them to metabolism pathways. It contains five modules which could be used individually or works together firmly to give the annotation and visualization of pseudogenes. It runs well on a Linux system. One file in faa format and the other in gbk format whose contents are consistent should be used as input. One metadata file having several groups of faa and gbk files could be input as well(See [Example_metadata](https://github.com/444thLiao/dysfunctional_dectector/blob/main/Example/L1_example.tsv "Example_metadata")). Additionally, if you want to get some reference genomes from the KEGG database and analyze them, just use the `--addbytext` parameter followed by the name of interested organism.
***

# Modules:

There exist five modules in the dysfunctional_detector. 
The `s1.annotation.py` module annotates the given sequence with `Kofamscan`, `interproscan` and `pseudofinder` sequentially.  Additionally, we merge the splitted pseudogene in the result of pseudofinder to  better further analysis. The shunt of input data could be viewed in the following graph.
![Shunt_of_s1](https://github.com/444thLiao/dysfunctional_dectector/blob/main/anno.drawio.png "shunt")
The `s2.refine_annotation.py` module refines the output of s1 using interproscan database. It highlights the kegg number of pseudogenized  genes as well. The output tells whether specific genes are present or absent and looks like:

| blank | genome_id |  
| K00XXX | 1 |  
| K00XXX | 0 | 
 

The `s3.detector.py` module prepare the result of s2 for further visualizations. It maps the genes of different genomes(in KO number form) to the same metabolic pathway in one single file. The step infomation of metabolism is included in the output as well. 

The remaining modules, `s4.overview_vis.py`and`s5.detailed_vis.py` are both visualization scripts. `s4.overview_vis.py` makes one heatmap showing the state of genes of all input genomes in one single metabolism pathway. `s5.detailed_vis.py` shows the detailed state of genes within a certain range of specific genome.
***

# Dependence download:

Three softwares should be downloaded to support the dysfunctional_detector.

Kofamscan is one tool annotating gene function by HMMER/HMMSEARCH against KOfam(a customized HMM database of KEGG Orthologs (KOs)). To make it work, you should donwload the [Kofamscan](https://www.genome.jp/ftp/tools/kofam_scan/kofam_scan-1.3.0.tar.gz "Kofamscan")
and its database[Ko_list](https://www.genome.jp/ftp/db/kofam/ko_list.gz "Ko_list")and
[Profile](https://www.genome.jp/ftp/db/kofam/profiles.tar.gz "Profile")

Interproscan is one tool telling which family one protein belongs to and the domains it contains. It has several software requirements: Perl5, Python3 and Java JDK/JRE version 11. You may go to its official [website](https://interproscan-docs.readthedocs.io/en/latest/InstallationRequirements.html "ipr_website") and [download_link](https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.68-100.0/interproscan-5.68-100.0-64-bit.tar.gz "ipr_download") to get it.

Pseudofinder annotates the pseudogenes in the given sequence. You could go to its [installing manuel](https://github.com/filip-husnik/pseudofinder/wiki/2.-Installing-Pseudofinder "Pseudofinder_install") for download instructions.
***

# How to set up:

Before setting up the dysfunctional_detector, make sure you have downloaded the python3 and Anaconda. Meanwhile, concerning system path has set. After that, we could start with the following commands.
```shell  
git clone https://github.com/444thLiao/dysfunctional_dectector
cd dysfunctional_detector
```
Then edit the static settings and system path in the initializaion.py and setup.sh respectively.
```python
python3 initialization.py
```
When all these operations has done, the script should be ready for use.
***

# Test and expected output:

(To be finished when all modules has done.)
***

# Commands:
 

The dysfunctional_detector has five parameters:

`-i`.`--infile`:This parameter should be followed by the metadata file containing the path of several gbk and faa file.

`-fi`.`--file_input`:This parameter should be followed by the gbk and faa file. Gbk file first and faa file next.The filename of gbk file should be the genome_id of the organism.

`-o`.`--folder_output`:This parameter should be followed by the output folder.

`-d`.`--dry_run`:This parameter will make the programme generate command only.

`-at`.`--addbytsxt`:This parameter should be followed by the organism name which you want to download accessory genomes and analyze.

For all the parameters, only`-o`is needed for each run.For `-i` and `-fi`, use one or both is OK. The left parameters are optional.
Following is one example command of the programme:
```python
python3 /mnt/storage3/yfdai/download/script/dysfunctional_dectector/main.py -i /mnt/storage3/yfdai/download/script/dysfunctional_dectector/Example/L1_example.tsv -fi /home-user/thliao/script/dysfunctional_dectector/Example/L1.faa /mnt/storage3/yfdai/download/script/dysfunctional_dectector/L1.gbk -o /mnt/storage3/yfdai/download/script/new_output -add ruegeria
```

***
# References:

## Pseudofinder
Syberg-Olsen MJ*, Graber AI*, Keeling PJ, McCutcheon JP, Husnik F. Pseudofinder: detection of pseudogenes in prokaryotic genomes, Molecular Biology and Evolution 2022, 39(7): msac153, doi: https://doi.org/10.1093/molbev/msac153.
## InterPro
The InterPro protein families and domains database: 20 years on Matthias Blum, Hsin-Yu Chang, Sara Chuguransky, Tiago Grego, Swaathi Kandasaamy, Alex Mitchell, Gift Nuka, Typhaine Paysan-Lafosse, Matloob Qureshi, Shriya Raj, Lorna Richardson, Gustavo A Salazar, Lowri Williams, Peer Bork, Alan Bridge, Julian Gough, Daniel H Haft, Ivica Letunic, Aron Marchler-Bauer, Huaiyu Mi, Darren A Natale, Marco Necci, Christine A Orengo, Arun P Pandurangan, Catherine Rivoire, Christian J A Sigrist, Ian Sillitoe, Narmada Thanki, Paul D Thomas, Silvio C E Tosatto, Cathy H Wu, Alex Bateman, Robert D Finn Nucleic Acids Research (2020), gkaa977, PMID: 33156333
## InterProScan
InterProScan 5: genome-scale protein function classification Philip Jones, David Binns, Hsin-Yu Chang, Matthew Fraser, Weizhong Li, Craig McAnulla, Hamish McWilliam, John Maslen, Alex Mitchell, Gift Nuka, Sebastien Pesseat, Antony F. Quinn, Amaia Sangrador-Vegas, Maxim Scheremetjew, Siew-Yit Yong, Rodrigo Lopez, Sarah Hunter Bioinformatics (2014), PMID: 24451626
## Kofamscan< br > 
Aramaki T., Blanc-Mathieu R., Endo H., Ohkubo K., Kanehisa M., Goto S., Ogata H.
KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold.
Bioinformatics. 2019 Nov 19. pii: btz859. doi: 10.1093/bioinformatics/btz859.
