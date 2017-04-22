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


maq match rhesusf1r1_alignment.map rhesus_ref_genome.bfa rhesus_f1r1.bfq
maq match rhesusf1r2_alignment.map rhesus_ref_genome.bfa rhesus_f1r2.bfq
maq match rhesusf2r1_alignment.map rhesus_ref_genome.bfa rhesus_f2r1.bfq
maq match rhesusf2r2_alignment.map rhesus_ref_genome.bfa rhesus_f2r2.bfq
maq match rhesusf3r1_alignment.map rhesus_ref_genome.bfa rhesus_f3r1.bfq
maq match rhesusf3r2_alignment.map rhesus_ref_genome.bfa rhesus_f3r2.bfq
maq match rhesusm1r1_alignment.map rhesus_ref_genome.bfa rhesus_m1r1.bfq
maq match rhesusm1r2_alginment.map rhesus_ref_genome.bfa rhesus_m1r2.bfq
maq match rhesusm2r1_alignment.map rhesus_ref_genome.bfa rhesus_m2r1.bfq
maq match rhesusm2r2_alignment.map rhesus_ref_genome.bfa rhesus_m2r2.bfq
maq match rhesusm3r1_alignment.map rhesus_ref_genome.bfa rhesus_m3r1.bfq
maq match rhesusm3r2_alignment.map rhesus_ref_genome.bfa rhesus_m3r2.bfq
~                                                                      
