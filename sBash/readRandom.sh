#!/bin/bash

#SBATCH --mem=8G
#SBATCH --account=def-jpviroid
#SBATCH --time=0:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

source /project/6003961/USERID/bin/virtenv/bin/activate

python /project/6003961/USERID/Scripts/G4Conserve/readRandom.py -chr $1
