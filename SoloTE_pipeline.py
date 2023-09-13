#SoloTE

import pysam
import sys
import re
import os
import pandas as pd
import numpy
from datetime import datetime
from multiprocessing import Pool
import multiprocessing
import subprocess
from pathlib import Path
import shutil


import argparse


solote_version = "1.09"
argparse_object = argparse.ArgumentParser(prog="SoloTE version "+solote_version,description="Analysis of transposable elements in single-cell RNA-Seq data using locus-specific expression")

argparse_object.add_argument("-b","--bam",help="Input BAM file with CB and UB tags",required=True)
argparse_object.add_argument("-t","--threads",help="Number of threads to use during the pipeline",type=int,required=True)
argparse_object.add_argument("-d","--outputdir",help="Directory to store the outputs of SoloTE (if it doesn't exist, it will be created)",required=True)
argparse_object.add_argument("-a","--teannotation",help="TE annotation file in BED format",required=True)
argparse_object.add_argument("-o","--outputprefix",help="Prefix for output files",required=True)
argparse_object.add_argument("--dual",help="Consider reads annotated to genes for calculation of TE expression (default = False, only consider non-genic reads).",action='store_true',required=False)
argparse_object.add_argument("--minoverlap",help="Minimum overlap (in bp) between a read and a TE .",required=False,default=1)

commandargs = argparse_object.parse_args()

inputfile = commandargs.bam
cpus = str(commandargs.threads)
outbase = commandargs.outputdir
TE_bed = commandargs.teannotation
outprefix = commandargs.outputprefix
use_dual_mode = commandargs.dual
min_overlap = commandargs.minoverlap


starting_time = datetime.now()

starting_time_formatted = starting_time.strftime("%H:%M:%S")
print("SoloTE started at "+ starting_time_formatted)

#validate PATH requirements
required_software = ['samtools','bedtools']
for software in required_software:
    if shutil.which(software) is None:
        print("[ERROR] "+software+" is not in your PATH environment variable. Add it to your path before running SoloTE")
        exit()
    else:
        print("[OK] "+software+" found!")


def processPerChrom(chromosome,inputfile,outdir,outbasename):

    countsname = outbasename+"_countpercell_"+chromosome+".counts"
    if os.path.exists(countsname) == True:
        if  os.path.getsize(countsname)>0:
                print(countsname+" exists. Skipping this file")
                return(countsname)
    else:
        print(f'Counts for chromosome {chromosome} are being generated in process: {os.getpid()}')

    temp_countsfile = outbasename+"_countpercell_"+chromosome+".counts.tmp"

    cmd="samtools view "+inputfile+" "+chromosome+"|awk 'BEGIN{FS=OFS=\"\\t\"}{for (tag_index=NF-2;tag_index<=NF;tag_index++){ gsub(\":Z:\",\"\\t\",$tag_index); split($tag_index,split_tag,\"\\t\");  tag_name=split_tag[1]; tag_value=split_tag[2]; tag_dict[tag_name] = tag_value}; print tag_dict[\"GN\"],tag_dict[\"CB\"],tag_dict[\"UB\"] }' > "+temp_countsfile
#    print(cmd)
    os.system(cmd)

    test_table = pd.read_table(temp_countsfile,sep="\t",header=None)
    test_table_count_nunique = test_table.groupby([0,1],as_index=False,).agg('nunique')
    test_table_count_nunique.sort_values([0,1],axis=0).to_csv(countsname,sep="\t",header=None,index=False)

    os.remove(temp_countsfile)

    return(countsname)



result_list = []

results_dir = os.path.abspath(outbase)
outdir = os.path.abspath(outbase)+"/"+outprefix+"_SoloTE_temp"

workingdir = os.getcwd()
finaldir = os.path.abspath(outbase)+"/"+outprefix+"_SoloTE_output"
outbase = outprefix

inputbam = os.path.abspath(inputfile)
TE_bed = os.path.abspath(TE_bed)
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

##Set multiprocessing start method to 'fork'
multiprocessing.set_start_method('fork')

print("SoloTE v"+solote_version+" started!")
print("SoloTE Home directory "+SoloTE_Home)
print("SoloTE executed from "+workingdir)
print("Results will be stored in "+results_dir)
print("Input BAM file: "+inputbam)
print("Input TE BED file: "+TE_bed)

mode = "Standard"
if use_dual_mode:
    mode = "Dual"
    print ("Dual mode enabled. SoloTE will calculate TE expression also considering reads annotated to genes.")  

Path(outdir).mkdir(parents=True, exist_ok=True)
os.chdir(outdir)
print("Currently working in temporary directory: "+outdir)



bam_with_proper_cb_ub_tags = "temp1.bam"

if os.path.exists(te_bam):
    print(te_bam+" exists in output folder. Skipping this step")
else:
#    cmd="samtools view --threads "+cpus+" -O BAM -o "+te_bam+" -L "+TE_bed+" "+bam_without_genes
    if mode == "Dual":
       cmd="samtools view -@ "+cpus+" -O BAM -o "+te_bam+" -L "+TE_bed+" -e '(exists([CB]) && exists([UB]) && [CB]!=\"-\" && [UB]!=\"-\") "+inputbam
    else:
       cmd="samtools view -@ "+cpus+" -O BAM -o "+te_bam+" -L "+TE_bed+" -e '(exists([CB]) && exists([UB]) && [CB]!=\"-\" && [UB]!=\"-\") && (!exists([GN]) || [GN]==\"-\")' "+inputbam
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
    temp_annotated_te_bam = "temp_annotated_te.bam"
    cmd="python "+annotateBAMpath+" "+te_bam+" "+selected_TEs+" "+temp_annotated_te_bam+" "+str(min_overlap)
    print(cmd)
    os.system(cmd)
    sorted_bam=annotated_te_bam+".sorted."
    cmd="samtools sort -@ "+cpus+" -O BAM -o "+annotated_te_bam+" "+temp_annotated_te_bam
    print(cmd)
    os.system(cmd)


if os.path.exists(full_bam):
    print(full_bam+" exists in output folder. Skipping this step")
else:
    cmd="samtools merge --threads "+cpus+" -o - "+inputbam+" "+annotated_te_bam+"|samtools view -@ "+cpus+" -O BAM -o "+final_bam+" -e 'exists([CB]) && exists([UB]) && exists([GN]) && [CB]!=\"-\" && [UB]!=\"-\" && [GN]!=\"-\"' --keep-tag GN,CB,UB"
    print(cmd)
    os.system(cmd)
    cmd="samtools index "+final_bam
    print(cmd)
    os.system(cmd)

full_sortedbam = full_bam

chrnames = []
for chromosome_info in pysam.idxstats(final_bam).split("\n"):
    chromosome_splitinfo = chromosome_info.split("\t")
    if(len(chromosome_splitinfo)>3 and int(chromosome_splitinfo[2])>0):
        chrname = chromosome_splitinfo[0]
        chrnames.append(chrname)

pool = Pool(processes=int(cpus))
#Generate counts
for chromosome in chrnames:
    pool.apply_async(processPerChrom,args=(chromosome,final_bam,outdir,outbase),callback=result_list.append)

pool.close()
pool.join()
#print(result_list)

allcountsfile = outbase+"_allcounts.txt"

if os.path.exists(allcountsfile):
    print(allcountsfile+" exists. Will be removed")
    os.remove(allcountsfile)

for countfile in result_list:
    if countfile is not None:
        os.system("cat "+countfile+" >> "+allcountsfile)



allcounts = pd.read_table(allcountsfile,header=None)
tecounts = allcounts[allcounts[0].str.contains("SoloTE")]
tecounts2 = tecounts[0].str.split("|",expand=True)
tecounts2.loc[tecounts2[4].isnull(),4] = tecounts2.loc[tecounts2[4].isnull(),1]
te_table = allcounts[allcounts[0].str.contains("SoloTE")]
te_annotation = tecounts2[4].str.split(":",expand=True)
te_final_df = pd.concat([te_table,te_annotation],axis=1)
te_final_df.columns=['TE_original_id','barcode','counts','Subfamily','Family','Class']
te_final_df.replace("\?","",inplace=True,regex=True)
locus_tes = te_final_df[te_final_df['TE_original_id'].str.contains("SoloTE\|chr")]

legacy_df = te_final_df.groupby(['TE_original_id','barcode'],as_index=False)['counts'].sum()
locus_df = locus_tes.groupby(['TE_original_id','barcode'],as_index=False)['counts'].sum()
class_df = te_final_df.groupby(['Class','barcode'],as_index=False)['counts'].sum()
family_df = te_final_df.groupby(['Family','barcode'],as_index=False)['counts'].sum()
subfamily_df = te_final_df.groupby(['Subfamily','barcode'],as_index=False)['counts'].sum()
family_df = te_final_df.groupby(['Family','barcode'],as_index=False)['counts'].sum()
family_df.replace('',numpy.nan,inplace=True)
family_df.dropna(subset=['Family'], inplace=True)

genecounts = allcounts[~allcounts[0].str.contains("SoloTE")]
genecounts.columns=['Gene_id','barcode','counts']
gene_df = genecounts.groupby(['Gene_id','barcode'],as_index=False)['counts'].sum()

legacy_df.columns = genecounts.columns
locus_df.columns = genecounts.columns
class_df.columns = genecounts.columns
family_df.columns = genecounts.columns
subfamily_df.columns = genecounts.columns

class_df['Gene_id'] = "SoloTE|"+class_df['Gene_id']
family_df['Gene_id'] = "SoloTE|"+family_df['Gene_id']
subfamily_df['Gene_id'] = "SoloTE|"+subfamily_df['Gene_id']

gene_legacytes = pd.concat([gene_df,legacy_df])
gene_locustes = pd.concat([gene_df,locus_df])
gene_classtes = pd.concat([gene_df,class_df])
gene_familytes = pd.concat([gene_df,family_df])
gene_subfamilytes = pd.concat([gene_df,subfamily_df])

gene_legacytes.to_csv(outbase+"_legacytes.txt",sep="\t",index=False)
gene_locustes.to_csv(outbase+"_locustes.txt",sep="\t",index=False)
gene_classtes.to_csv(outbase+"_classtes.txt",sep="\t",index=False)
gene_familytes.to_csv(outbase+"_familytes.txt",sep="\t",index=False)
gene_subfamilytes.to_csv(outbase+"_subfamilytes.txt",sep="\t",index=False)


print("Creating final results directory")
Path(finaldir).mkdir(parents=True, exist_ok=True)
print(finaldir+" was created")

generateMatrix_Rscript=SoloTE_Home+"/generate_mtx.R"
filelist = [outbase+"_legacytes.txt",outbase+"_locustes.txt",outbase+"_classtes.txt",outbase+"_familytes.txt",outbase+"_subfamilytes.txt"]
for countsfile in filelist:
	mtx_outname=countsfile.replace(".txt","_MATRIX")
	cmd="Rscript "+generateMatrix_Rscript+" "+countsfile+" "+mtx_outname
	print(cmd)
	os.system(cmd)
	os.rename(mtx_outname,finaldir+"/"+mtx_outname)


total_gene_counts = gene_df['counts'].sum()
total_te_counts = legacy_df['counts'].sum()
total_locuste_counts = locus_tes['counts'].sum()
total_subfte_counts = total_te_counts-total_locuste_counts
total_counts = total_gene_counts+total_te_counts
total_gene_percentage = round(total_gene_counts*100/total_counts,3)
total_te_percentage = round(total_te_counts*100/total_counts,3)
total_locuste_percentage = round(total_locuste_counts*100/total_te_counts,3)
total_subfte_percentage = round(total_subfte_counts*100/total_te_counts,3)

print("A total of "+str(total_counts)+" UMIs are in the final matrix.")
print("Of these,")
print("\t"+str(total_gene_counts)+" ("+str(total_gene_percentage)+"%) correspond to genes.")
print("\tand "+str(total_te_counts)+" ("+str(total_te_percentage)+"%) correspond to TEs.")
print("TE detected UMIs are distributed as follows:")
print("\tLocus-specific TEs: "+str(total_locuste_counts)+" UMIs ("+str(total_locuste_percentage)+"%).")
print("\tSubfamily TEs: "+str(total_subfte_counts)+" ("+str(total_subfte_percentage)+"%).")

##Generate TE statistics
outstats_basename = outbase+"_SoloTE.stats"
outstats_filename = workingdir+"/"+outstats_basename
print("Creating "+outstats_basename+" TE statistics file")
outstats_filehandle = open(outstats_filename,"w")
outstats_filehandle.write("A total of "+str(total_counts)+" UMIs are in the final matrix.\n")
outstats_filehandle.write("Of these,\n")
outstats_filehandle.write("\t"+str(total_gene_counts)+" ("+str(total_gene_percentage)+"%) correspond to genes.\n")
outstats_filehandle.write("\tand "+str(total_te_counts)+" ("+str(total_te_percentage)+"%) correspond to TEs.\n")
outstats_filehandle.write("TE detected UMIs are distributed as follows:\n")
outstats_filehandle.write("\tLocus-specific TEs: "+str(total_locuste_counts)+" UMIs ("+str(total_locuste_percentage)+"%).\n")
outstats_filehandle.write("\tSubfamily TEs: "+str(total_subfte_counts)+" ("+str(total_subfte_percentage)+"%).\n")
outstats_filehandle.close()
print ("Finished creating "+outstats_basename)
print("SoloTE finished with "+inputfile)

finishing_time = datetime.now()
elapsed_time = finishing_time - starting_time
finishing_time_formatted = finishing_time.strftime("%H:%M:%S")
print("SoloTE finished at "+finishing_time_formatted)
print("SoloTE total running time: "+str(elapsed_time))




