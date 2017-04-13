##Bia's notes - final project

Using Galaxy to do RNA-Seq Assemble  
Data from the study is at [this NCBI GEO page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17274)  
There is a total of 36 files (human, chimp, rhesus/male1, male2, female1, female2/rep1, rep2)  
The group decided to divide these 36 files among us (9 files each) to start the assemble  
I am responsible for files 28-36, which corresponds to:  
-GSM432625	Rhesus female 2 rep2
-GSM432626	Rhesus female 3 rep1
-GSM432627	Rhesus female 3 rep2
-GSM432628	Rhesus male 1 rep1
-GSM432629	Rhesus male 1 rep2
-GSM432630	Rhesus male 2 rep1
-GSM432631	Rhesus male 2 rep2
-GSM432632	Rhesus male 3 rep1
-GSM432633	Rhesus male 3 rep2     

###To collect the data  
1. Galaxy>Get Data > EBI SRA  
2. Searched for "GSE17274"  
3. Clicked on "Submission(Read/Analysis) - SRA010277  
4. Clicked on designated files, getting them in FASTQ format directly into Galaxy (FASTQ files (Galaxy))

##FASTQC analysis  
1. Used FASTQC tool on Galaxy
2. NGS:QC and manipulation > FASTQC
3. Selected all files to be set as input to FASTQC
4. Downloaded all results to github repository and analyzed them

##Trimming reads  
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