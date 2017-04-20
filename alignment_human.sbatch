#!/bin/bash
#Submit this script with: sbatch thefilename
#SBATCH -t 8:00:00   # walltime
#SBATCH -e ./alignment_human.error
#SBATCH -o ./alignment_human.out
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 16   # total number of processor cores in this job
#SBATCH -J "alignment_human"   # job name
#SBATCH --mail-user=biakemi@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load intel
module load allinea
ulimit -s unlimited
module load maq/0.7.1


maq match humanf1r1_alignment.map human_ref_genome.bfa human_f1r1.bfq
maq match humanf1r2_alignment.map human_ref_genome.bfa human_f1r2.bfq 
maq match humanf2r1_alignment.map human_ref_genome.bfa human_f2r1.bfq
maq match humanf2r2_alignment.map human_ref_genome.bfa human_f2r2.bfq
maq match humanf3r1_alignment.map human_ref_genome.bfa human_f3r1.bfq
maq match humanf3r2_alignment.map human_ref_genome.bfa human_f3r2.bfq
maq match humanm1r1_alignment.map human_ref_genome.bfa human_m1r1.bfq  
maq match humanm1r2_alginment.map human_ref_genome.bfa human_m1r2.bfq
maq match humanm2r1_alignment.map human_ref_genome.bfa human_m2r1.bfq 
maq match humanm2r2_alignment.map human_ref_genome.bfa human_m2r2.bfq
maq match humanm3r1_alignment.map human_ref_genome.bfa human_m3r1.bfq
maq match humanm3r2_alignment.map human_ref_genome.bfa human_m3r2.bfq
