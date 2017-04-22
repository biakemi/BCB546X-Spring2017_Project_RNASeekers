#!/bin/bash
#Submit this script with: sbatch thefilename
#SBATCH -t 5:00:00   # walltime
#SBATCH -e ./fastq2bfq_rhesus.error
#SBATCH -o ./fastq2bfq_rhesus.out
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 16   # total number of processor cores in this job
#SBATCH -J "MAQ_fastq2bfqrhesus"   # job name
#SBATCH --mail-user=biakemi@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load intel
module load allinea
ulimit -s unlimited
module load maq/0.7.1

maq fastq2bfq SRR032140.fastq rhesus_f1r1.bfq
maq fastq2bfq SRR032141.fastq rhesus_f1r2.bfq
maq fastq2bfq SRR032142.fastq rhesus_f2r1.bfq
maq fastq2bfq SRR032143.fastq rhesus_f2r2.bfq
maq fastq2bfq SRR032144.fastq rhesus_f3r1.bfq
maq fastq2bfq SRR032145.fastq rhesus_f3r2.bfq
maq fastq2bfq SRR032146.fastq rhesus_m1r1.bfq
maq fastq2bfq SRR032147.fastq rhesus_m1r2.bfq
maq fastq2bfq SRR032148.fastq rhesus_m2r1.bfq
maq fastq2bfq SRR032149.fastq rhesus_m2r2.bfq
maq fastq2bfq SRR032150.fastq rhesus_m3r1.bfq
maq fastq2bfq SRR032151.fastq rhesus_m3r2.bfq
