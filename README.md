# README #

SPLICIFY is the proteogenomic pipeline for differential splice variant identification.

This pipeline was developed by M A Komor, please cite.

**Version 1.0**

## Introduction

SPLICIFY consists of 2 steps - a genomic and a proteomic part. 

Within the genomic part of the proteogenomic pipeline (step 1) RNA-seq analysis is performed, including quality and adapter trimming with Trimmomatic, reads mapping with STAR, differential splicing analysis with rMATS and post-processing steps to transform
the splice variants into potential protein variant sequence database (FASTA), that
can be further used with MaxQuant - a search engine to identify MS/MS spectra.

Within the proteomic part of the proteogenomic pipeline (step2) down-stream analysis of MaxQuant output is performed with the use of the results from step 1. Variant peptides are extracted and quantified. The step 2 of the pipeline produces a final table with both RNA and protein isoform information.


## Software requirements

Required software:

* Linux
* python 2.7.x
* R version 3.2.x or higher
* bedtools (v2.23.0 of later)
* 3rd party software requirements:
	* Trimmomatic
		* java (here: version 1.7.0_17)
	* STAR
		* x86-64 compatible processors
		* 64 bit Linux or Mac OS X
		* 30GB of RAM for human genome
	* rMATS
		* numPy, sciPy
		* samtools (v1.2 or later)
		
Trimmomatic, STAR and rMATS are included in the proteogenomic pipeline package. The installation of the software is required once the user chooses to use a different version of these tools.

Recommended  software:

* transeq from EMBOSS package
* Windows
* MaxQuant

In case of the use of other search engine for mass spectra identification, the output should be processed by the user so that it fits the step 2 of the pipeline (resembles MaxQuant output). See section *Using a different search engine* for details.

## Quick start
For a quick start the following parameters should be adjusted in the configuration file.

To run **proteogenomics_step1.py**:

* `Input files`:
	* `seq` - `paired` or `single` for paired-end or single-end sequencing
	* `Group1` - RNA-seq *fastq.gz files from condition 1, please see *Configuration file* for details 
	* `Group2` - RNA-seq *fastq.gz files from condition 2, please see *Configuration file* for details 
* `Trimmomatic`:
	* `crop`, `minlen` - the length to which the reads should be trimmed, cannot be longer than the original read length in `*fastq.gz`
	* `threads` - cannot be more than capacity of your machine
* `STAR`:
	* `runThreaN` - cannot be more than capacity of your machine

To run **proteogenomics_step2.py**:

* `MaxQuant output`:
	* `evidence` - path to evidence.txt file, if you cannot run MaxQuant, please see *Using a different search engine* for details 
	* `peptides` - path to peptides.txt file, if you cannot run MaxQuant, please see *Using a different search engine* for details  
* `Extract`:
	* `threads` - cannot be more than capacity of your machine
* `Quantitative Analysis`
	* `group1_column_nr` - numbers of columns with intensity for each sample from condition1, please see *Configuration file* for details 
	* `group2_column_nr` - numbers of columns with intensity for each sample from condition2, please see *Configuration file* for details 

## How to run ?

###1: Genomic part of the pipeline:

	python proteogenomics_step1.py -c src/config.ini
or

	python proteogenomics_step1.py --config src/config.ini

Output defined in configuration file, e.g.:

	[matsToFasta]
	outputPepFasta = data/rmats/output.pep.fasta

 

### 2: MaxQuant
 Use step1 output (`output.pep.fasta`) and  human protein database as databases used in MaxQuant for protein identification.

### 3: Proteomic part of the pipeline:

	python proteogenomics_step2.py -c src/config.ini
or

	python proteogenomics_step2.py --config src/config.ini

### 4: Output 

Output is defined in configuration file, e.g.:

	[Extract]
	output_prefix = data/prot/

In the output directory there will be the following files, please see *Output explanation* for details:

* `variantPeptides.txt` - all peptides specifically mapping to the splice variants identified in the genomics part, together with the RNA information
* `variantPeptidesQA.txt` - file as `variantPeptides.txt` including columns with the results of quantitative analysis
* `variantPeptidesCanonical(QA).txt`, `variantPeptidesNonCanonical(QA).txt` - `variantPeptides(QA).txt` divided into 2 files based on the information if the peptides is in the canonical protein database
* `variantPeptidesInclusion(QA).txt`, `variantPeptidesExclusion(QA).txt` - `variantPeptides(QA).txt` divided into 2 files based on the information if the peptides is inclusion or exclusion variant specific

## Other options

By default full analysis is performed. In case the user wants to skip parts of the analysis, (e.g. no quantitative analysis due to lack of replicates in the proteomics or if the processes was interrupted, resuming from the last step,) it can be done with the use of the additional arguments.

### Step 1

	usage: proteogenomics_step1.py [-h] [-t [TRIM]] [-s [STAR]] [-r [RMATS]]
	                               [-m [MATSTOFASTA]] -c CONFIG
	
	Run step 1 of the proteogenomic pipeline. By default the full analysis will be
	performed. If you want to skip some steps use False in the commands: 
	-t, --trim, -s, --star, -r, --rmats, -m, --matstofasta, 
	but make sure the config file has paths to intermediate outputs and the file names 
	are correct.
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -t [TRIM], --trim [TRIM]
	                        run trimmomatic, True/False (default: True)
	  -s [STAR], --star [STAR]
	                        run STAR, True/False (default: True)
	  -r [RMATS], --rmats [RMATS]
	                        run rMATS, True/False (default: True)
	  -m [MATSTOFASTA], --matstofasta [MATSTOFASTA]
	                        run matsToFasta, True/False (default: True)
	
	required argument:
	  -c CONFIG, --config CONFIG
	                        configuration file, required (e.g. :
	                        src/config.ini)

### Step 2

	usage: proteogenomics_step2.py [-h] [-e [EXTRACT]] [-g [GETRNA]]
	                               [-q [QUANTIFY]] -c CONFIG
	
	Run step 2 of the proteogenomic pipeline. By default the full analysis will be
	performed. If you want to skip some steps use False in the commands: -e,
	--extract, -g, --getrna, -q, --quantify, but make sure the config file has paths
	to intermediate outputs and the file names are correct. In case of no
	replicates, it is recommended to skip the quantitative analysis
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -e [EXTRACT], --extract [EXTRACT]
	                        run ExtractVariantPeptides, True/False (default: True)
	  -g [GETRNA], --getrna [GETRNA]
	                        run getRNAInformation, True/False (default: True)
	  -q [QUANTIFY], --quantify [QUANTIFY]
	                        run quantitativeAnalysis , True/False (default: True)
	
	required argument:
	  -c CONFIG, --config CONFIG
	                        configuration file, required (e.g.:
	                        src/config.ini)
	

## Configuration file

Configuration file is the required input for both steps of the proteogenomic pipeline. It contains all the parameters necessary to run the pipeline. Sample configuration files for paired-end (`src/configPE.ini`) and single-end (`src/configSE.ini`) RNA-seq input are in `src` directory.

Configuration file is divided in 8 sections, all sections except for `[Quantitative Analysis]` are required for the pipeline. As replicates are needed for the quantitative analysis on the peptide level, the user can skip this step.

### Explanation of the sections:

#### Input files ####
This section contains the paths to `*fastq.gz` files divided in 2 groups, representing 2 conditions for which differential splicing analysis will be performed, and the information if paired or single-end RNA sequencing was performed.

**Section keys:**

* `seq` - parameter can only take `single` or `paired` values
* `Group1` - list of comma separated `*.fastq.gz` files from **condition 1**, for paired-end data `*.fastq.gz` with forward (`R1`) and reverse (`R2`) reads should be colon separated. The files per sample should be located next to each other in the string as in *Example section for paired-end data* below.
* `Group2` - list of comma separated `*.fastq.gz` files from **condition 2**, formatting should be done as in Group1

**Example section for paired-end data**

    [INPUT FILES]

	seq = paired
	Group1 = data/fastq/sample1_R1.fq.gz:data/fastq/sample1_R2.fq.gz, data/fastq/sample2_R1.fq.gz:data/fastq/sample2_R2.fq.gz
	Group2 = data/fastq/sample3_R1.fq.gz:data/fastq/sample3_R2.fq.gz, data/fastq/sample4_R1.fq.gz:data/fastq/sample4_R2.fq.gz

**Example section for single-end data**

    [INPUT FILES]

	seq = single
	Group1 = data/fastq/sample1.fastq.gz, data/fastq/sample2.fastq.gz
	Group2 = data/fastq/sample3.fastq.gz, data/fastq/sample4.fastq.gz

#### Trimmomatic ####
This section contains parameters and input files necessary for running Trimmomatic, a tool for quality and adapter trimming. The reads need to be trimmed to the length set by the user, due to the requirements of the later tools in the pipeline. 

**Section keys:**

* `trimmomatic_jar` - path to trimmmatic jar. The program is already included in the pipeline package, so it does not have to be edited. The user can update this parameter with the path to different installation of Trimmomatic.
* `threads` - number of threads used for trimming
* `illuminaclip` - path to adapter sequenced. The adapter file is included in the pipeline package for TruSeq Illumina adapters, this parameter does not have to be edited.
* `minlen` , `crop` - read length after trimming, shorter reads will be discarded, longer reads will be cropped to this length
* `output` - path to the output
* For the following parameters please see the Trimmomatic manual:
	* `phred`, `leading`, `trailing`, `avgqual`, `slidingwindow`

Additional parameters will be discarded

**Example section for Trimmomatic for 2x125bp paired-end Illumina RNA-seq data**

	[TRIMMOMATIC]

	trimmomatic_jar = src/trimmomatic-0.33.jar
	threads = 8
	phred = 33
	illuminaclip = data/adapters/TruSeq3-PE.fa:2:30:10
	leading = 20
	trailing = 20
	crop = 120
	avgqual = 20
	slidingwindow = 4:20
	minlen = 120
	output = data/trimmed/

#### STAR ####
This section contains parameters and input files for reads mapping with STAR.

**Section keys**

* `pathToStar` - path to STAR tool. The program is already included in the proteogenomic pipeline package, the user can update this path with a different version of STAR.
* `genomeDir` - path to genome index generated with STAR. In the proteogenomic pipeline package a STAR genome index to genome built hg19 is already included. If the user wishes to use different STAR genome index, the path should be edited.
* `outFileNamePrefix` - path to the output
* `runThreadN` - number of threads used for mapping
* For the following parameters please see the STAR manual:
	* `outSAMtype`, `outSAMattributes`

**Example section for STAR for Illumina RNA-seq data**

	[STAR]

	pathToStar = src/STAR
	genomeDir = data/genome/star_index
	outSAMtype = BAM SortedByCoordinate
	readFilesCommand = zcat
	outFileNamePrefix = data/mapped/
	runThreadN = 8
	outSAMattributes = All

#### rMATS ####
This section contains the parameters and input files for differential splicing analysis with rMATS.

**Section keys**

* `pathTorMATS` - path to rMATS tool. In the proteogenomic pipeline package the rMATS tool is already included. To perform the analysis with another version of rMATS, the user should edit this parameter.
* `gtf` - path to gtf file. In the proteogenomic pipeline package a gtf file for genome built hg19 with HGNC gene symbols is included.
* `o` - path to output directory
* For the following parameters please see the rMATS documentation:
	* `analysis` - can only take `P` or `U` values, corresponding to paired or unpaired analysis

**Examples section for rMATS**

	[RMATS]

	pathTorMATS = src/rMATS.3.2.5
	gtf = data/genome/refGene-hg19.gtf
	analysis = U
	o = data/rmats

#### matsToFasta ####

matsToFasta is a three step approach to obtain fasta sequences from rMATS differential splice variants.

1. **rMATS output to bed file** - rMATS output is FDR filtered and transformed into a bed file, bed file with the splice variants can be viewed in IGV
2. **Bed file to fasta file** - bed file is transformed into a fasta file with nucleotide sequences of the splice variants
3. **3 frame translation** - nucleotide sequences are translated in 3 frames into a database of potential protein splice variants, a fasta file is produced that has to be taken along for mass spectra identification together with human protein database.

**Section keys**

* `comparisonName` - name of the conditions/experiments compared in the differential analysis, should be string without any special characters (in particular `$,;-%` should not be used)
* `input` - accepts only values: `ReadsOnTargetAndJunctionCounts` or `JunctionCountsOnly`, indicating which rMATS output files should be taken along for further analysis. Please see rMATS documentation for more details.
* `fdr` - FDR threshold for rMATS output, 0.05 is recommended. Please see rMATS documentation for more details. 
* `outputBed` - path to output bed file including splice variants.
* `ifMXE` - accepts values : `yes` or `no` indicating if mutually exclusive events should be taken along in the analysis. Due to the high false positives rate for these events in rMATS analysis, the user can exclude them from the pipeline.
* `Genome_Fasta` - path to genome fasta file. In the proteogenomic pipeline package there is a fasta file included for genome built hg19. The user can adjust this parameter to a local fasta file. Genome fasta file should be indexed, to index the genome run `samtools index *.fa`, where `*.fa` is the genome fasta file.
* `outputFasta` - path to output fasta file with nucleotide sequences of the splice variants
* `outputPepFasta` - path to output fasta file with amino acid sequences of the splice variants
* `pathToTranseq` - path to transeq tool from the EMBOSS package, `FALSE` indicates that transeq is installed globally and does not need a path. If transeq is not installed, a R script will be run to perform the same analysis. Please see `src/threeFrame.R` for details 

**Example section for matsToFasta**

	[matsToFasta]

	comparisonName = group1vsgroup2
	input = ReadsOnTargetAndJunctionCounts
	fdr = 0.05
	outputBed = data/rmats/output.bed
	ifMXE = yes
	Genome_Fasta = data/genome/hg19.fa
	outputFasta = data/rmats/output.fasta
	outputPepFasta = data/rmats/output.pep.fasta
	pathToTranseq = FALSE

#### MaxQuant output ####
This section contains the MaxQuant output files.

**Section keys**

* `evidence` - path to evidence.txt file from MaxQuant output. 
* `peptides` - path to peptides.txt file from MaxQuant output.

If a different search engine is used, the mass spectra identification file should be adjusted. For details, please see section *Using different search engines*.
	
	[MaxQuant output]
	evidence = data/prot/evidence.txt
	peptides = data/prot/peptides.txt

#### Extract ####
This section contains parameters needed to extract variant peptides from the MaxQuant output.

**Section keys**

* `threads` - number of threads to be used
* `canonical` - path to the database with human canonical protein sequences. In the proteogenomic pipeline package the Swissprot canonical database is already included.
* `output_prefix` - path to the output, the file names are generated by the tool.

**Example section for Extract**

	[Extract]
	threads = 8
	canonical = data/prot/uniprot-canonical-swissprot.fasta
	output_prefix = data/prot/

#### Quantitative Analysis ####
This section contains the parameters needed for differential peptide expression analysis.

**Section keys**

* `group1_column_nr` - comma separated column numbers with peptide intensities per sample in the peptides.txt file, for all the samples in **group 1/condition 1** from the differential splicing analysis on RNA level, the columns are counted from 1.
* `group2_column_nr`- comma separated column numbers with peptide intensities per sample in the peptides.txt file, for all the samples in **group 2/condition 2** from the differential splicing analysis on RNA level, the columns are counted from 1.
* `imputation` - accepts values `yes` or `no` indicating if the imputation algorithm should be performed for the peptides with missing intensity values

**Example section for Quantitative Analysis**


	[Quantitative Analysis]
	group1_column_nr = 66, 67, 68
	group2_column_nr = 63, 64, 65
	imputation = yes

## Using a different search engine
It is possible to run the pipeline with different search engines than MaxQuant, however, in that case the full analysis is not supported. The output of the search engine has to be adjusted by the user to resemble evidence and peptides file from MaxQuant. 

The required files should be tab-delimited and have these columns:

### evidence.txt

This is a file with each identification in separate row. 

* Sequence - column with peptide sequences
* Experiment - column with sample names where peptide was identified

Sequence    | Experiment 
----------- | -----------
AAPPLPR     | sample1    
AAPPLPR     | sample2    
AAPPLPR     | sample3    
ACPGLTYHR   | sample2
ACPGLTYHR   | sample3        

In evidence.txt file :

	Sequence	Experiment 
	AAPPLPR		sample1    
	AAPPLPR		sample2
	AAPPLPR		sample3        
	ACPGLTYHR	sample2
	ACPGLTYHR	sample3    

### peptides.txt

This is a file where each row is unique for a peptide sequence. The same peptide sequence cannot be in multiple rows.

* `Sequence` - peptide sequence
* `IntensitySampleX` - raw intensity value for this peptide in sample X

Sequence    | IntensitySample1 | IntensitySample2 | IntensitySample3 
----------- | ---------------- | ---------------- | ---------------- 
AAPPLPR     | 31773000         | 59181000         | 26111000         
ACPGLTYHR   | 0                | 17978000         | 10142000  

In peptides.txt file:

	Sequence	IntensitySample1	IntensitySample2	IntensitySample3 
	AAPPLPR		31773000			59181000			26111000         
	ACPGLTYHR	0					17978000			10142000  

In this case, if `group1 = samples1, sample2` and `group2 = sample3`, values passed to the configuration file are:

	[Quantitative Analysis]
	group1_column_nr = 2, 3
	group2_column_nr = 4

## Output explanation

### variantPeptides(QA).txt

This is the main output of the proteogenomic pipeline. 

Each row represents a peptide with an isoform it's mapping to. If a peptide maps to multiple isoforms, it will be in multiple rows.

Columns explanation:

* `Sequence` - peptide sequence
* `ID` - rMATS splice variant ID where the peptide maps
* `Event` - type of splicing event: se - skipped exon, a3ss - alternative 3' splice site, a5ss - alternative 5' splice site, ri- retained intron, mxe - mutually exclusive exons
* `Gene` - gene name
* `Incl_Excl` - if inclusion or exclusion variant of the isoform. If MXE, inclusion means inclusion of exon1 and exclusion of exon2, exclusion means exclusion of exon1 and inclusion of exon2
* `coordinates` - coordinates of the isoform in the format: chromosome;strand;exonStart-exonEnd;...;exonStart-exonEnd;_frameOfTranslation
* `FDR` - false discovery rate of splice variant calculated my rMATS
* `InclLevel1` - mean of Inclusion Levels of all the samples in Group1, inclusion levels per sample are calculated by rMATS
* `ExclLevel1` - mean of Exclusion Levels of all the samples in Group1, `ExclLevel1 = 1 - InclLevel1`
* `InclLevel2` - mean of Inclusion Levels of all the samples in Group2, inclusion levels per sample are calculated by rMATS
* `ExclLevel2` - mean of Exclusion Levels of all the samples in Group2, `ExclLevel2 = 1 - InclLevel2`
* `IncLevelDifference` - difference in inclusion levels between Group1 and Group2, `IncLevelDifference = InclLevel1 - InclLevel2`
* `Sample` - sample names in which the peptide was identified
* `Aberrant` - TRUE - non-canonical peptide, FALSE - canonical peptide (in SwissProtDB)
* `SecondVariant` - TRUE - peptides for both variants identified, FALSE - peptides identified only for this variant
* `Peptide` - Information if the peptide spans exon-exon junction (split peptide), spans exon-intron junction (spanning peptide) or maps on the spliced region (on target)
* `NumberOfisoforms` - number of isoforms to which the peptide maps
* `Left` - peptide length on the left side of the junction
* `Right` - peptide length on the right side of the junction
* `Location` - peptide coordinates
* `IntensitySampleX` - normalized log10 transformed intensity of the peptide in SampleX, this value was taken along to limma analysis. If `imputation = no`, NA's for missing values.
* `limma.logFC` - limma log2 fold change from the differential peptide expression between the groups
* `limma.p.value` - limma p-value from the differential peptide expression between the groups
* `limma.adj.p.value` - p-value adjusted with Benjamini–Hochberg correction
* `NumberOfSamples` - number of samples in which peptide was identified and intensity was available, useful if `imputation = yes` included to extract intensities from the imputed values for missing values.




## Contribution guidelines

* Author
	* Gosia Komor
* Code review
	* Thang Pham 
* Running tests
	* Linda Bosch

