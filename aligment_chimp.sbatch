#!/bin/bash
#Submit this script with: sbatch thefilename
#SBATCH -t 8:00:00   # walltime
#SBATCH -e ./alignment_chimp.error
#SBATCH -o ./alignment_chimp.out
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 16   # total number of processor cores in this job
#SBATCH -J "alignment_chimp"   # job name
#SBATCH --mail-user=biakemi@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load intel
module load allinea
ulimit -s unlimited
module load maq/0.7.1

maq match chimpf1r1_alignment.map chimp_ref_genome.bfa chimp_f1r1.bfq
maq match chimpf1r2_alignment.map chimp_ref_genome.bfa chimp_f1r2.bfq 
maq match chimpf2r1_alignment.map chimp_ref_genome.bfa chimp_f2r1.bfq
maq match chimpf2r2_alignment.map chimp_ref_genome.bfa chimp_f2r2.bfq
maq match chimpf3r1_alignment.map chimp_ref_genome.bfa chimp_f3r1.bfq
maq match chimpf3r2_alignment.map chimp_ref_genome.bfa chimp_f3r2.bfq
maq match chimpm1r1_alignment.map chimp_ref_genome.bfa chimp_m1r1.bfq  
maq match chimpm1r2_alginment.map chimp_ref_genome.bfa chimp_m1r2.bfq
maq match chimpm2r1_alignment.map chimp_ref_genome.bfa chimp_m2r1.bfq 
maq match chimpm2r2_alignment.map chimp_ref_genome.bfa chimp_m2r2.bfq
maq match chimpm3r1_alignment.map chimp_ref_genome.bfa chimp_m3r1.bfq
maq match chimpm3r2_alignment.map chimp_ref_genome.bfa chimp_m3r2.bfq
