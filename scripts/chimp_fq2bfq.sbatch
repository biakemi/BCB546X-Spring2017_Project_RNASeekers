#!/bin/bash
#Submit this script with: sbatch thefilename
#SBATCH -t 5:00:00   # walltime
#SBATCH -e ./fastq2bfq_chimp.error
#SBATCH -o ./fastq2bfq_chimp.out
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 16   # total number of processor cores in this job
#SBATCH -J "MAQ_fastq2bfqchimp"   # job name
#SBATCH --mail-user=biakemi@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load intel
module load allinea
ulimit -s unlimited
module load maq/0.7.1

maq fastq2bfq SRR032128.fastq chimp_f1r1.bfq
maq fastq2bfq SRR032129.fastq chimp_f1r2.bfq
maq fastq2bfq SRR032130.fastq chimp_f2r1.bfq
maq fastq2bfq SRR032131.fastq chimp_f2r2.bfq
maq fastq2bfq SRR032132.fastq chimp_f3r1.bfq
maq fastq2bfq SRR032133.fastq chimp_f3r2.bfq
maq fastq2bfq SRR032134.fastq chimp_m1r1.bfq
maq fastq2bfq SRR032135.fastq chimp_m1r2.bfq
maq fastq2bfq SRR032136.fastq chimp_m2r1.bfq
maq fastq2bfq SRR032137.fastq chimp_m2r2.bfq
maq fastq2bfq SRR032138.fastq chimp_m3r1.bfq
maq fastq2bfq SRR032139.fastq chimp_m3r2.bfq
