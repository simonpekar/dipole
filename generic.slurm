#!/bin/sh

# Connect through SSH
# ssh -Y laura.domine@hopper.calcul.labos.polytechnique.fr
# Send a file to hopper
# scp utils.cpp laura.domine@hopper.calcul.labos.polytechnique.fr:~
# Get a file from hopper
# scp laura.domine@hopper.calcul.labos.polytechnique.fr:~/displit/tree.root .
# Launch interactive run
# srun --ntasks=1 ./main 1

#SBATCH --job-name=test
#SBATCH --time=120
#SBATCH --account=medium
#SBATCH --mail-user=simon.pekar@polytechnique.edu
#SBATCH --ntasks=1
#SBATCH --partition=interactive 
 
module add root
 
./main