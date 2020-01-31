#!/bin/bash

#SBATCH --job-name=G4RNAScreener
#SBATCH --mem=30G
#SBATCH --array=0
#SBATCH --account=def-jpviroid
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

source /project/6003961/USERID/bin/virtenv/bin/activate

/project/6003961/USERID/bin/g4rna_screener/screen.py /home/USERID/scratch/G4prediction/Data/$1/SplitFile/$2 -a /project/6003961/USERID/bin/g4rna_screener/G4RNA_2016-11-07.pkl -w 60 -s 10 -c description cGcC G4H G4NN sequence start end -e > /home/USERID/scratch/G4prediction/Data/$1/CSVFile/$2

