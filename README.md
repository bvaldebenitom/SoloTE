# SoloTE
[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs42003--022--04020--5-yellow)](https://www.nature.com/articles/s42003-022-04020-5)

SoloTE README

## 0. DEPENDENCIES
SoloTE requires the following tools to be installed, and available in your PATH environment variable:
- Samtools v1.16 or higher (http://www.htslib.org/download/)
- BEDTools v2.29.2 or higher (https://github.com/arq5x/bedtools2/releases)
- R v4 or higher (https://www.r-project.org/)

Also, Python3.9.5 or higher should be available in your computer, along with the following modules:
- Pysam (https://pysam.readthedocs.io/en/latest/installation.html)
- Pandas v1.5.0 or higher (https://pandas.pydata.org/)


## 1. SETTING UP NECESSARY FILES

RepeatMasker files (\*rm.out) for several genomes, can be found at the UCSC webpage. The helper utility, `SoloTE_RepeatMasker_to_BED.py`, is packed with SoloTE, and it streamlines the download of RepeatMasker file from UCSC, and conversion to BED format. It can be run like this:
```
python SoloTE_RepeatMasker_to_BED.py -g GenomeVersion
```
where
- GenomeVersion: genome version identifier. For example, `hg38` for human, and `mm10` for mouse.

Additionally, the utility can be called with the `-l` option, and it will list identifiers that can be supplied as the `-g` parameter, of all available genomes at UCSC (first 5 lines shown next):
```
python SoloTE_RepeatMasker_to_BED.py -l
ailMel1	| Panda [Ailuropoda melanoleuca, Dec. 2009 (BGI-Shenzhen 1.0/ailMel1)]
allMis1	| American alligator [Alligator mississippiensis, Aug. 2012 (allMis0.2/allMis1)]
anoCar1	| Lizard [Anolis carolinensis, Feb. 2007 (Broad/anoCar1)]
anoCar2	| Lizard [Anolis carolinensis, May 2010 (Broad AnoCar2.0/anoCar2)]
anoGam1	| A. gambiae [Anopheles gambiae, Feb. 2003 (IAGEC MOZ2/anoGam1)]
```

If you have your own RepeatMasker file and/or a file corresponding to Transposable Elements obtained from another tool, make sure to adapt it to the following format for SoloTE:
```
sequenceName	startPosition	endPosition	sequenceName|startPosition|endPosition|TE_Subfamily:TE_Family:TE_Class|strand	score(optional)	.
```

So, column 4, the ID, is a concatenation of the locus of the TE and its identifiers at the Subfamily, Family and Class level. The file should look like this:
```
chr1	11505	11675	chr1|11505|11675|L1MC5a:L1:LINE|25.1|-	25.1	-
chr1	11678	11780	chr1|11678|11780|MER5B:hAT-Charlie:DNA|29.4|-	29.4	-
chr1	15265	15355	chr1|15265|15355|MIR3:MIR:SINE|23.0|-	23.0	-
chr1	18907	19048	chr1|18907|19048|L2a:L2:LINE|33.8|+	33.8	+
chr1	19972	20405	chr1|19972|20405|L3:CR1:LINE|31.2|+	31.2	+
```


## 2. RUNNING SOLOTE

Once everything is set up, you can run the SoloTE script:
```
python SoloTE_pipeline.py --threads NumberOfThreads --bam BAMfile --teannotation BEDfile --outputprefix Prefix --outputdir OutputDirectory
```

where

`--threads`: Number of threads to use

`--bam`: BAM file

`--teannotation`: TE annotation in BED format

`--outputprefix`: Output files prefix

`--outputdir`: Output directory (if it doesn't exist, it will be created)

