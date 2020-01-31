#!/bin/bash

#SBATCH --job-name=pG4Annotation
#SBATCH --mem=30G
#SBATCH --array=0
#SBATCH --account=def-jpviroid
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

source /project/6003961/USERID/bin/virtenv/bin/activate

python ../../Scripts/G4Conserve/G4Annotation.py -chr $1
