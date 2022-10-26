# SoloTE

SoloTE README

## 0. DEPENDENCIES
SoloTE requires the following tools to be installed, and available in your PATH environment variable:
- Samtools v1.11 or higher (http://www.htslib.org/download/)
- BEDTools v2.29.2 or higher (https://github.com/arq5x/bedtools2/releases)
- R v4 or higher (https://www.r-project.org/)

Also, Python3.9.5 or higher should be available in your computer, along with the following modules:
- Pysam (https://pysam.readthedocs.io/en/latest/installation.html)


## 1. SETTING UP NECESSARY FILES

RepeatMasker files (\*rm.out) for several genomes, can be found at the UCSC webpage. A conversion utility is conveniently packaged with SoloTE, in order to tranform the RepeatMasker out file to the BED format required by SoloTE.
It can be run like this:
```
convertRMOut_to_SoloTEinput.sh RepeatMaskerOutfile NameForTEfile
```
where
- RepeatMaskerOutfile: out file obtained from UCSC or with RepeatMasker
- NameForTEfile: Name for the new, properly formatted, TE annotation file in BED format (required for the pipeline)

If you have your own RepeatMasker file and/or a file corresponding to Transposable Elements obtained from another tool, make sure to adapt it to the following format for SpatialTE:
```
sequenceName	startPosition	endPosition	sequenceName|startPosition|endPosition|TE_Subfamily:TE_Family:TE_Class|strand	score(optional)	.
```

So, column 4, the ID, is a concatenation of the locus of the TE and its identifiers at the Subfamily, Family and Class level. The file should look like this:
```
chr1	3000001	3002128	chr1|3000001|3002128|L1_Mus3:L1:LINE|-	12955	-
chr1	3003153	3003994	chr1|3003153|3003994|L1Md_F:L1:LINE|-	1216	-
chr1	3003994	3004054	chr1|3003994|3004054|L1_Mus3:L1:LINE|-	234	-
chr1	3004041	3004206	chr1|3004041|3004206|L1_Rod:L1:LINE|+	3685	+
```


## 2. RUNNING SOLOTE

Once everything is set up, you can run the SoloTE script with the configuration file:
```
python SoloTE_pipeline.py BAMfile Threads OutputDir TEannotationBEDfile OutputPrefix
```
where
- BAMfile: relative or full path to the CellRanger / STAR BAM file containing the CB and UB tags
- Threads: Number of threads to use during the pipeline 
- OutputDir: Name for output directory
- TEannotationBEDfile: Name of TE annotation file in BED format generated with convertRMOut_to_SoloTEinput.sh
- OutputPrefix: Prefix for output files

