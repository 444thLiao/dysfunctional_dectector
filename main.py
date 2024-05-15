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





# parse args
@click.command()
@click.option('--infile','-i',help="An Metadata file as input. See example file in XXXX",required=False,default=None)
@click.option('-fi',"--file_input",type = str, nargs = 2 ,required = False, prompt = "Enter the input file name", help = "Please input faa and gbk file you want to analyze")
@click.option('-o',"--folder_output",type = str, nargs = 1 ,required = True, prompt = "Enter the output folder name", help = "Please output path of folder you want to store the analysis results. ")
@click.option("-d","dry_run",help="Generate command only.",default=False,required=False,is_flag=True,)
@click.option('--addbytext','-add',help="Input a name for searching  accessory genomes and use them to correct the genome you want to annotated.",required=False,default=False)
def cli(odir):






