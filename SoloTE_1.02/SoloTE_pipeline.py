#SoloTE v1

import pysam
import sys
import re
import os
import pandas as pd
from datetime import datetime
from multiprocessing import Pool
import subprocess
from pathlib import Path
import shutil

starting_time = datetime.now()

starting_time_formatted = starting_time.strftime("%H:%M:%S")
print("SoloTE started at "+ starting_time_formatted)

#validate PATH requirements
required_software = ['samtools','bedtools']
for software in required_software:
    if shutil.which(software) is None:
        print(software+" is not in your PATH environment variable. Add it to your path before running SoloTE")
        exit()
    else:
        print(software+" found!")


def processPerChrom(chromosome,inputfile,outdir,outbasename):
    count = 0
    print(f'Process: {os.getpid()}')
    
    countsname = outbasename+"_countpercell_"+chromosome+".counts"
    if os.path.exists(countsname):
        print(countsname+" exists. Skipping this file")
        return(countsname)

    os.system("samtools view "+inputfile+" "+chromosome+"| awk 'BEGIN{FS=OFS=\"\t\"}{gn=\"\";gx=\"\";cb=\"\";ub=\"\";for(i=12;i<=NF;i++){if($i~ /CB:/){sub(\".*:\",\"\",$i);cb=$i};if($i~ /UB:/){sub(\".*:\",\"\",$i);ub=$i};if($i~ /GX:/){gx=$i};if($i~ /GN:/){sub(\"GN:Z:\",\"\",$i);gn=$i} }; if(cb!=\"-\" && ub != \"-\"){print gn,cb,ub}}'|sort -t'\t' -k1,1 -k2,2 -k3,3|bedtools groupby -g 1,2 -c 3 -o count_distinct >  "+countsname)

    linenumber = subprocess.run(["wc","-l",countsname],stdout=subprocess.PIPE,text=True)
    if re.match("^0",linenumber.stdout):
        print(linenumber.stdout)
    else:
        return(countsname)


inputfile = sys.argv[1]
headerlines = pysam.view("-H",inputfile).split("\n")
chrnames = []
for headerline in headerlines:
#    if "@SQ" in headerline and "chr" in headerline:
    if "@SQ" in headerline:
        split_sq = headerline.split("\t")
        chrname = re.sub("SN:","",split_sq[1])
        chrnames.append(chrname)
print(chrnames)
co_line = (list(filter(lambda x:'@CO' in x,headerlines)))
co_line = co_line[0]
print(co_line.split("--"))


outsammultnmax_value = (list(filter(lambda x:'outSAMmultNmax 1' in x,co_line.split("--"))))
print(len(outsammultnmax_value))
if len(outsammultnmax_value)!=1:
    print("outSAMmultNmax 1 option missing in BAM file")
    exit()

samtags_from_star = (list(filter(lambda x:'outSAMattributes' in x,co_line.split("--"))))
print(samtags_from_star[0])

if "CB" in samtags_from_star[0] and "UB" in samtags_from_star[0]:
    print("CB and UB tags present in BAM file")
else:
    print("CB and UB tags not present in BAM file")
    exit()

cpus = sys.argv[2]
outbase = sys.argv[3]
TE_bed = sys.argv[4]

outprefix = sys.argv[5]

result_list = []

outdir = outbase+"/"+outprefix+"_SoloTE_temp"

workingdir = os.getcwd()
#finaldir = os.getcwd()+"/"+outbase+"_SoloTE_output"
finaldir = outbase+"/"+outprefix+"_SoloTE_output"

outbase = outprefix

inputbam = os.path.abspath(inputfile)
bam_without_genes=outbase+"_nogenes.bam"
bam_without_genes_newtag=outbase+"_nogenes_newtag.bam"
bam_without_genes_oldtag=outbase+"_nogenes_oldtag.bam"
temp_bam="temp.bam"
te_bam=outbase+"_nogenes_overlappingtes.bam"
te_bamAsbed=outbase+"_nogenes_overlappingtes.bed"
selected_TEs=outbase+"_selectedtes.bed"
gene_bam=outbase+"_genes.bam"
full_bam=outbase+"_full.bam"
full_sortedbam=outbase+"_full_sorted.bam"
final_bam=outbase+"_final.bam"
annotated_te_bam=outbase+"_teannotated.bam"

SoloTE_Home=os.path.dirname(__file__)

Path(outdir).mkdir(parents=True, exist_ok=True)
os.chdir(outdir)

if os.path.exists(bam_without_genes) and os.path.exists(gene_bam):
    print(bam_without_genes+" and "+gene_bam+" exist in output folder. Skipping this step")
else:
    cmd="samtools view --threads "+cpus+" -d GN -U "+bam_without_genes_oldtag+" -O BAM -o "+temp_bam+" "+inputbam
    print(cmd)
    os.system(cmd)   

    cmd="samtools view -d GN:- -o "+bam_without_genes_newtag+" -U "+gene_bam+" -O BAM "+temp_bam
    print(cmd)
    os.system(cmd)

    cmd="samtools cat --threads "+cpus+" -o "+bam_without_genes+" "+bam_without_genes_oldtag+" "+bam_without_genes_newtag
    print(cmd)
    os.system(cmd)

if os.path.exists(te_bam):
    print(te_bam+" exists in output folder. Skipping this step")
else:
    cmd="samtools view --threads "+cpus+" -O BAM -o "+te_bam+" -L "+TE_bed+" "+bam_without_genes
    print(cmd)
    os.system(cmd)
    cmd="samtools index "+te_bam
    print(cmd)
    os.system(cmd)

if os.path.exists(te_bamAsbed):
    print(te_bamAsbed+" exists in output folder. Skipping this step")
else:
    cmd="bedtools bamtobed -i "+te_bam+" -split > "+te_bamAsbed
    print(cmd)
    os.system(cmd)

if os.path.exists(selected_TEs):
    print(selected_TEs+" exists in output folder. Skipping this step")
else:
    cmd="bedtools intersect -a "+TE_bed+" -b "+te_bamAsbed+" -u > "+selected_TEs
    print(cmd)
    os.system(cmd)

if os.path.exists(annotated_te_bam):
    print(annotated_te_bam+" exists in output folder. Skipping this step")
else:
    annotateBAMpath=SoloTE_Home+"/annotateBAM.py"
    cmd="python "+annotateBAMpath+" "+te_bam+" "+selected_TEs+" "+annotated_te_bam
    print(cmd)
    os.system(cmd)

if os.path.exists(full_bam):
    print(full_bam+" exists in output folder. Skipping this step")
else:
    cmd="samtools cat --threads "+cpus+" -o "+full_bam+" "+gene_bam+" "+annotated_te_bam
    print(cmd)
    os.system(cmd)

if os.path.exists(full_sortedbam):
    print(full_sortedbam+" exists in output folder. Skipping this step")
else:
    cmd="samtools sort --threads "+cpus+" -O BAM -o "+full_sortedbam+" "+full_bam
    print(cmd)
    os.system(cmd)

if os.path.exists(final_bam):
    print(final_bam+" exists in output folder. Skipping this step")
else:
    cmd="samtools view -U "+final_bam+" --tag CB:- -@ 8 "+full_sortedbam+" > /dev/null"
    print(cmd)
    os.system(cmd)
    cmd="samtools index "+final_bam
    print(cmd)
    os.system(cmd)


pool = Pool(processes=int(cpus))
#Generate counts
for chromosome in chrnames:
    pool.apply_async(processPerChrom,args=(chromosome,final_bam,outdir,outbase),callback=result_list.append)

pool.close()
pool.join()
print(result_list)

allcountsfile = outbase+"_allcounts.txt"

if os.path.exists(allcountsfile):
    print(allcountsfile+" exists. Will be removed")
    os.remove(allcountsfile)

for countfile in result_list:
    if countfile is not None:
        os.system("cat "+countfile+" >> "+allcountsfile)

genes_tsv_filename = outbase +"_genes.tsv"
barcodes_tsv_filename = outbase+"_barcodes.tsv"

genes_filename1 = outbase + "_genes.txt"
telocus_filename1 = outbase + "_locustes.txt"
tesubf_filename1 = outbase + "_subftes.txt"

os.system("grep -v \"SoloTE\" "+allcountsfile+" > "+genes_filename1)
os.system("grep \"SoloTE\" "+allcountsfile+"|grep \"chr\" > "+telocus_filename1)
os.system("grep \"SoloTE\" "+allcountsfile+"|grep -v \"chr\" > "+tesubf_filename1)


new_outfiles = []

for filename in [genes_filename1,telocus_filename1,tesubf_filename1]:
    outfile = re.sub(".txt","_2.txt",filename)
    new_outfiles.append(outfile)
    os.system("sort -t'\t' -k1,1 -k2,2 -k3,3 "+filename+"|bedtools groupby -g 1,2 -c 3 -o sum > "+outfile)

final_countsfile = outbase + "_allcounts_final.txt"
te_threshold = 3
os.system("cat "+new_outfiles[0]+"> "+final_countsfile)
os.system("cat "+new_outfiles[1]+" "+new_outfiles[2]+"|awk '($NF>="+str(te_threshold)+"){print $0}' >> "+final_countsfile)

os.system("awk '{print $1}' "+final_countsfile+"|sort -u|awk 'BEGIN{FS=OFS=\"\t\"}{print $1,$1}' > "+genes_tsv_filename)
os.system("awk '{print $2}' "+final_countsfile+"|sort -u > "+barcodes_tsv_filename)


genenumber = subprocess.run(["wc","-l",genes_tsv_filename],stdout=subprocess.PIPE,text=True)
barcodenumber = subprocess.run(["wc","-l",barcodes_tsv_filename],stdout=subprocess.PIPE,text=True)
allcounts_number = subprocess.run(["wc","-l",final_countsfile],stdout=subprocess.PIPE,text=True)
marketmatrix_line3 = genenumber.stdout.split(" ")[0]+" "+barcodenumber.stdout.split(" ")[0]+" "+allcounts_number.stdout.split(" ")[0]

output_marketmatrix_filename = outbase+"_allcounts_matrix.mtx"
os.system("echo \"%%MatrixMarket matrix coordinate integer general\n%\" > "+output_marketmatrix_filename)
os.system("echo "+marketmatrix_line3+" >> "+output_marketmatrix_filename)
generateMatrix_Rscript=SoloTE_Home+"/generate_mtx.R"
os.system("Rscript "+generateMatrix_Rscript+" "+barcodes_tsv_filename+" "+genes_tsv_filename+" "+final_countsfile+" "+output_marketmatrix_filename)

print("Creating final results directory")
Path(finaldir).mkdir(parents=True, exist_ok=True)
print(finaldir+" was created")

os.rename(barcodes_tsv_filename,finaldir+"/barcodes.tsv")
os.rename(genes_tsv_filename,finaldir+"/features.tsv")
os.rename(output_marketmatrix_filename,finaldir+"/matrix.mtx")

##Generate TE statistics
outstats_basename = outbase+"_SoloTE.stats"
outstats_filename = workingdir+"/"+outstats_basename

input_te_file = re.sub(".txt","_2.txt",telocus_filename1)
file_table = pd.read_table(input_te_file,header=None,sep="\t")
tes_ls_n = str(len(file_table[0].unique()))
tes_ls_celln = str(len(file_table[1].unique()))
tes_ls_totalumicount = str(file_table[2].sum())


threshold_table = file_table[file_table[2]>=te_threshold]
tes_threshold_ls_n = str(len(threshold_table[0].unique()))
tes_threshold_ls_celln = str(len(threshold_table[1].unique()))
tes_threshold_ls_totalumicount = str(threshold_table[2].sum())

input_te_file = re.sub(".txt","_2.txt",tesubf_filename1)
file_table = pd.read_table(input_te_file,header=None,sep="\t")
tes_subf_n = str(len(file_table[0].unique()))
tes_subf_celln = str(len(file_table[1].unique()))
tes_subf_totalumicount = str(file_table[2].sum())

threshold_table = file_table[file_table[2]>=te_threshold]
tes_threshold_subf_n = str(len(threshold_table[0].unique()))
tes_threshold_subf_celln = str(len(threshold_table[1].unique()))
tes_threshold_subf_totalumicount = str(threshold_table[2].sum())

print("Creating "+outstats_basename+" TE statistics file")
outstats_filehandle = open(outstats_filename,"w")
outstats_filehandle.write("Category\tTotal TEs\tTotal cells\tTotal UMIs\n")
outstats_filehandle.write("Locus-specific TEs\t"+tes_ls_n+"\t"+tes_ls_celln+"\t"+tes_ls_totalumicount+"\n")
outstats_filehandle.write("Locus-specific TEs (UMIs >= "+str(te_threshold)+")\t"+tes_threshold_ls_n+"\t"+tes_threshold_ls_celln+"\t"+tes_threshold_ls_totalumicount+"\n")
outstats_filehandle.write("Subfamily-specific TEs\t"+tes_subf_n+"\t"+tes_subf_celln+"\t"+tes_subf_totalumicount+"\n")
outstats_filehandle.write("Subfamily-specific TEs (UMIs >= "+str(te_threshold)+")\t"+tes_threshold_subf_n+"\t"+tes_threshold_subf_celln+"\t"+tes_threshold_subf_totalumicount+"\n")
outstats_filehandle.close()
print ("Finished creating "+outstats_basename)

print("SoloTE finished with "+inputfile)

finishing_time = datetime.now()
elapsed_time = finishing_time - starting_time
finishing_time_formatted = finishing_time.strftime("%H:%M:%S")
print("SoloTE finished at "+finishing_time_formatted)
print("SoloTE total running time: "+str(elapsed_time))

