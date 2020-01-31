#!/bin/bash
#SBATCH --job-name=enum_large_fas
#SBATCH --account=def-jpviroid
#SBATCH --mem-per-cpu=5000M
#SBATCH -t 0-00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

#Activate virtal environment of python
source /project/6003961/USERID/bin/virtenv/bin/activate

python /project/6003961/USERID/Scripts/fasta_enumerator.py --fasta-file /home/USERID/scratch/G4prediction/Data/$1/$2 --nt-limit 5000000 --large --output /home/USERID/scratch/G4prediction/Data/$1/SplitFile/ -e

