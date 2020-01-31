#!/bin/bash

#SBATCH --mem=8G
#SBATCH --account=def-jpviroid
#SBATCH --time=0:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

source /project/USERID/vana2406/bin/virtenv/bin/activate

python /project/USERID/vana2406/Scripts/G4Conserve/generateRandom.py -chr $1
