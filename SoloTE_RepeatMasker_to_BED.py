
import argparse
import pandas
import requests
from colorama import Fore, init
from tqdm import tqdm

init(autoreset=True)

argparse_object = argparse.ArgumentParser(description="SoloTE RepeatMasker to BED")

argparse_object.add_argument("-g","--genome",help="Genome assembly identifier to obtain RepeatMasker annotation (available genomes can be checked with option -l).")
argparse_object.add_argument("-l","--list",help="List available genomes at the UCSC database.",action='store_true',required=False)

commandargs = argparse_object.parse_args()

list_genomes = commandargs.list
genome_assembly = commandargs.genome

mode = "RepeatMasker"

def download(url: str, fname: str, chunk_size=1024):
	resp = requests.get(url, stream=True)
	total = int(resp.headers.get('content-length', 0))
	with open(fname, 'wb') as file, tqdm(
		desc=fname,
		total=total,
		unit='iB',
		unit_scale=True,
		unit_divisor=1024,
	) as bar:
		for data in resp.iter_content(chunk_size=chunk_size):
			size = file.write(data)
			bar.update(size)


if list_genomes == True:
	api_url = "https://api.genome.ucsc.edu/list/ucscGenomes"
	response = requests.get(api_url)
	response = response.json()
	ucscGenomes = response["ucscGenomes"]

	for genome in ucscGenomes:
		current = ucscGenomes[genome]
		genome_description = current['description']
		scientificname = current['scientificName']
		organism = current['organism']
		print(genome+"\t| "+organism+" ["+scientificname+", "+genome_description+"]")
	exit()


api_url = "https://api.genome.ucsc.edu/list/ucscGenomes"
response = requests.get(api_url)
response = response.json()
ucscGenomes = response["ucscGenomes"]

if ucscGenomes.get(genome_assembly) is None:
	print(Fore.RED+"[ERROR] ",end='')
	print("Genome assembly "+genome_assembly+" not available. Check available genomes with option -l.")
	exit()
else:
	print(Fore.GREEN+"[OK] ",end='')
	print("Genome assembly "+genome_assembly+" found.")


if mode == 'UCSC':
	api_url = "http://api.genome.ucsc.edu/getData/track?genome="+genome_assembly+";track=rmsk;maxItemsOutput=1"
	response = requests.get(api_url)
	response = response.json()
	rmsk_url = response["dataDownloadUrl"] 
	print("[LOG] URL to fetch RepeatMasker file: "+rmsk_url)

	rmsk_filename = genome_assembly+'_ucsc_rmsk.txt.gz'
	#r = requests.get(rmsk_url, allow_redirects=True)
	#open(rmsk_filename,'wb').write(r.content)
	print("[LOG] Downloading RepeatMasker file to "+rmsk_filename)
	download(url=rmsk_url,fname=rmsk_filename)


	rmsk_bed_filename = genome_assembly+'_ucsc_rmsk.bed'
	print("[LOG] Beginning conversion of RepeatMasker file to BED format")
	rmsk_table = pandas.read_csv(rmsk_filename,compression="gzip",header=None,sep="\t")
	rmsk_table.columns = ["bin","swScore","milliDiv","milliDel","milliIns","genoName","genoStart","genoEnd","genoLeft","strand","repName","repClass","repFamily","repStart","repEnd","repLeft","ID"]
	#print(rmsk_table)
	rmsk_table['fixed_milliDiv'] = rmsk_table['milliDiv'].astype(str).str.replace("([0-9]$)",".\\1",regex=True)
	rmsk_table['genoStart'] = rmsk_table['genoStart'].astype(str)
	rmsk_table['genoEnd'] = rmsk_table['genoEnd'].astype(str)
	rmsk_table['genoName']+"|"+str(rmsk_table['genoStart'])+"|"+str(rmsk_table['genoEnd'])
	#rmsk_table['te_name'] = rmsk_table['genoName']+"|"+rmsk_table['genoStart']+"|"+rmsk_table['genoEnd']+"|"+rmsk_table['repClass']+":"+rmsk_table['repFamily']+":"+rmsk_table['repName']+"|"+rmsk_table['fixed_milliDiv']+"|"+rmsk_table['strand']
	rmsk_table['te_name'] = rmsk_table['genoName']+"|"+rmsk_table['genoStart']+"|"+rmsk_table['genoEnd']+"|"+rmsk_table['repName']+":"+rmsk_table['repFamily']+":"+rmsk_table['repClass']+"|"+rmsk_table['fixed_milliDiv']+"|"+rmsk_table['strand']
	rmsk_table_bed = rmsk_table[['genoName','genoStart','genoEnd','te_name','fixed_milliDiv','strand']]
	rmsk_table_bed = rmsk_table_bed[rmsk_table_bed['te_name'].str.contains("LINE|SINE|LTR|DNA")]
	rmsk_table_bed.to_csv(rmsk_bed_filename,sep="\t",header=None,index=False)
	print(Fore.GREEN+"[OK] ",end='')
	print("Finished generating "+rmsk_bed_filename)


if mode == "RepeatMasker":
	rmsk_url = "https://hgdownload.soe.ucsc.edu/goldenPath/"+genome_assembly+"/bigZips/"+genome_assembly+".fa.out.gz"

	if genome_assembly == "hs1":
		rmsk_url = "https://hgdownload.soe.ucsc.edu/goldenPath/"+genome_assembly+"/bigZips/"+genome_assembly+".repeatMasker.out.gz"

	print("[LOG] URL to fetch RepeatMasker file: "+rmsk_url)
	rmsk_filename = genome_assembly+'.fa.out.gz'
	print("[LOG] Downloading RepeatMasker file to "+rmsk_filename)
	download(url=rmsk_url,fname=rmsk_filename)

	rmsk_bed_filename = genome_assembly+'_rmsk.bed'
	print("[LOG] Beginning conversion of RepeatMasker file to BED format")
	rmsk_out = pandas.read_csv(rmsk_filename,compression="gzip",skiprows=3,header=None,sep=" ",skipinitialspace=True)
	rmsk_out.columns = ['SW_score','percDiv','percDel','percIns','querySeq','queryStart','queryEnd','queryLeft','strand','matchingRepeat','repeatClass_Family','repeatBegin','repeatStart','repeatEnd','ID']
	repeatinfo = rmsk_out['repeatClass_Family'].str.split("/",expand=True)
	rmsk_out['repeatID'] = rmsk_out['matchingRepeat']+":"+repeatinfo[1]+":"+repeatinfo[0]
	rmsk_out['percDiv'] = rmsk_out['percDiv'].astype(str)
	rmsk_out['queryStart'] = rmsk_out['queryStart'].astype(str)
	rmsk_out['queryEnd'] = rmsk_out['queryEnd'].astype(str)
	rmsk_out['fixedStrand'] = rmsk_out['strand'].str.replace("C","-")
	rmsk_out['te_name'] = rmsk_out['querySeq']+"|"+rmsk_out['queryStart']+"|"+rmsk_out['queryEnd']+"|"+rmsk_out['repeatID']+"|"+rmsk_out['percDiv']+"|"+rmsk_out['fixedStrand']
	rmsk_out_bed = rmsk_out[['querySeq','queryStart','queryEnd','te_name','percDiv','fixedStrand']]
	rmsk_out_bed = rmsk_out_bed[rmsk_out['repeatClass_Family'].str.contains("LINE|SINE|LTR|DNA|RC")]
	rmsk_out_bed = rmsk_out_bed[~rmsk_out_bed['querySeq'].str.contains("chrna|_fix|_random|_alt|chrUn")]
	rmsk_out_bed.to_csv(rmsk_bed_filename,sep="\t",header=None,index=False)
	print(Fore.GREEN+"[OK] ",end='')
	print("Finished generating "+rmsk_bed_filename)


os.remove(rmsk_filename)






