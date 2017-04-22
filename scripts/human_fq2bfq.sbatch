#!/bin/bash
#Submit this script with: sbatch thefilename
#SBATCH -t 5:00:00   # walltime
#SBATCH -e ./fastq2bfq_human.error
#SBATCH -o ./fastq2bfq_human.out
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 16   # total number of processor cores in this job
#SBATCH -J "MAQ_fastq2bfqhuman"   # job name
#SBATCH --mail-user=biakemi@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load intel
module load allinea
ulimit -s unlimited
module load maq/0.7.1 

maq fastq2bfq SRR032116.fastq human_f1r1.bfq
maq fastq2bfq SRR032117.fastq human_f1r2.bfq
maq fastq2bfq SRR032118.fastq human_f2r1.bfq
maq fastq2bfq SRR032119.fastq human_f2r2.bfq
maq fastq2bfq SRR032120.fastq human_f3r1.bfq
maq fastq2bfq SRR032121.fastq human_f3r2.bfq
maq fastq2bfq SRR032122.fastq human_m1r1.bfq
maq fastq2bfq SRR032123.fastq human_m1r2.bfq
maq fastq2bfq SRR032124.fastq human_m2r1.bfq
maq fastq2bfq SRR032125.fastq human_m2r2.bfq
maq fastq2bfq SRR032126.fastq human_m3r1.bfq
maq fastq2bfq SRR032127.fastq human_m3r2.bfq
