#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

module use ~/modulefiles.local

module load mlst

species=cglabrata
ref_genome=CBS138_ASM254v2

while IFS= read -r line; do
    file="$line"
    python pubmlst_rest.py "$file" "$species" "$ref_genome"
    sleep 5
done < fasta_subsets.txt
