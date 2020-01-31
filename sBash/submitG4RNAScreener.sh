#!/bin/bash

#SBATCH --job-name=submitG4RNAScreener
#SBATCH --account=def-jpviroid
#SBATCH --time=168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user anais.vannutelli@usherbrooke.ca

startTime=`date +%s`

# Get list of MSA files in Gblocks repertory to compute Best ML RAxML tree on them
list_chr=$(ls /home/USERID/scratch/G4prediction/Data/ | grep chr)

for chr in $list_chr;
do
   list_file=$(ls /home/USERID/scratch/G4prediction/Data/$chr/SplitFile/)
   for file in $list_file;
   do
      jobs_nb=$(squeue -u USERID | wc -l)
      if [ $jobs_nb -le 999 ];
      then
         echo "sbatch /project/6003961/USERID/Scripts/G4Screener.sh $chr $file"
         sbatch /project/6003961/USERID/Scripts/G4Screener.sh $chr $file
         sleep 0.1
      else
         sleep 42
      fi
   done
done
echo $startTime
