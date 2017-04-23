##Bia's notes - final project

Using Galaxy to do RNA-Seq Assemble  
Data from the study is at [this NCBI GEO page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17274)  
There is a total of 36 files (human, chimp, rhesus/male1, male2, female1, female2/rep1, rep2)  
The group decided to divide these 36 files among us (9 files each) to start the assemble  
All files were distributed evenly to the group members, so each one of us did the sama steps below: 

###To collect the data  (BM, JL, EC, XF)
1. Galaxy>Get Data > EBI SRA  
2. Searched for "GSE17274"  
3. Clicked on "Submission(Read/Analysis) - SRA010277  
4. Clicked on designated files, getting them in FASTQ format directly into Galaxy (FASTQ files (Galaxy))

##FASTQC analysis  (BM, JL, EC, XF)
Since the study does not give information on how they processed the raw reads, we decided to follow the standard quality control steps on the raw data. We used Galaxy due to its friendly-user interface and also because some members of our groups already had some experience with it.  
1. Used FASTQC tool on Galaxy  
2. NGS:QC and manipulation > FASTQC  
3. Selected all files to be set as input to FASTQC  
4. Downloaded all results to github repository and analyzed them

##Trimming reads  (BM, JL, EC, XF)
1. Since the files are compacted (fastq.gz), we need to unzip in order to load into the trimming tool
2. Clicked on "edit attributes" of each file
2. Clicked on tab "Convert Format"
3. Selected "Convert fastq.gz to fastq" and clickec on "Convert"
4. The Trimmomatic tool on Galaxy accepts only "fastqsanger" files as input
5. Clicked on "edit attributes" > Datatype
6. Change datatype to "fastqsanger"
7. Run Trimmomatic  
**Parameters**  
- Perform initial ILLUMINACLIP step? - **YES**
- Select Trimmomatic operation to perform - **SLIDINGWINDOW**
- Number of bases to average across - **4**
- Average quality required - **20**

##Alignment of reads  (BM)
The study uses the [MAQ](http://maq.sourceforge.net/maq-man.shtml) software to map and align the RNA-Seq reads to the reference genomes.  

### Download the genomes (JL)
Download the genomes, following the information available from the paper and [Supplementary Methods](http://genome.cshlp.org/content/20/2/180/suppl/DC1). 



####Convert fastq to bfq (BM)
MAQ requires the files to be in a special format to run its functions, so the first step was to convert the fasta (genomes) and fastq (reads) files into bfa and bfq respectively.  
I used condo to run this part of the analysis, given that condo already had the software installed and allowed me to run several jobs at the same time, making the process faster than expected. 
-Used default parameters, just like the paper  

`$maq fastq2bfq input_file.fastq output_file.bfq`   

-Did separate runs for human, chimp and rhesus. Scripts are on github repository, under names **human_fq2bfq.sbatch**, **chimp_fq2bfq.sbatch** and **rhesus_fq2bfq.sbatch**

####Convert fasta to bfa  (BM)
-Used default parameters, just like the paper

`$ maq fasta2bfa input_file.fasta output_file.bfa`  

-Did everything in one script. Script is on github repository under name **XXXX** (will upload soon)  

#### Align the reads  (BM)
-`match` command aligns reads on bfq format to the reference in bfa format  
-The reads used on the study are single-end, therefore, we just put one file as input file  

`$ maq match output_file.map reference.bfa input_file.bfq`  

-Did separate runs for human, chimp and rhesus. Scripts are on github repository, under names **aligment_chimp.sbatch**, **alignment_human.sbatch** and **alignment_rhesus.sbatch**  

#### Converting to mapview format  (BM)
-Following the paper steps, we converted the alignments format to the mapview format using the `mapview` command  

`$ maq mapview output_file.map > output_file.txt`  

-The commands were run directly on the command line, but all the commands are on scripts on the github repository, under the names **mapview_human**, **mapview_chimp**, **mapview_rhesus**  

## Defining orthologous exons  (BM, JL)
-Downloaded the genomic coordinates for the human genome [Ensembl database](http://www.ensembl.org/Homo_sapiens/Info/Index) (release 88 instead of release 50).  
- Since the file had exons, transcripts, genes, etc, we are going to select the lines that only corresponds to exons. We are going to do that by using `grep`  

`$ grep "exons" file.gtf > annotations`  

-Uploaded the file to Galaxy > Fetch Alignments/Sequences > Extract Genomic DNA  
-We used the human genome from UCSC (hg38) and standard parameters  
-We got sequence data of 194008 exons. The paper got for 520023 exons. That difference can be due to the usage of a more fragmented genome that the paper research group used compared to ours, given that we used a more recent version of the human genome for this analysis.  



## Differential expression analysis (EC, XF)

### (EC)

### Redo analysis in Figure 2 of the paper (XF)
To check if our differential gene expression analysis can reproduce the results of conserved sexually dimorphic gene expression patterns as described in the Figure 2 of the paper, normalized read counts of four genes (RBM4, PAQR3, PGM2, and KDM5C) were analyzed and visualized using R. The R script `GeneExpressionPattern_Fig2.R` was included in the **script** folder.

Briefly, reads counts of the four genes were extracted from the input file `primate_norm_read_count.csv`. Data were cleaned with `dplyr` and visualized with `ggplot2`. The output figures were included in the **figures** folder.
