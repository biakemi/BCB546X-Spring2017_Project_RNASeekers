#!/bin/bash
#Submit this script with: sbatch thefilename
#SBATCH -t 1:00:00   # walltime
#SBATCH -e ./genome_conversion.error
#SBATCH -o ./genome_conversion.out
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 16   # total number of processor cores in this job
#SBATCH -J "genome_conversion"   # job name
#SBATCH --mail-user=biakemi@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load intel
module load allinea
ulimit -s unlimited
module load maq/0.7.1

maq fasta2bfa Galaxy52-[Homo_sapiens.GRCh38.dna.primary.assembly.fa].fasta human_ref_genome.bfa
maq fasta2bfa Galaxy53-[panTro4.fa].fasta chimp_ref_genome.bfa
maq fasta2bfa Galaxy56-[rheMac8.fa].fasta rhesus_ref_genome.bfa
