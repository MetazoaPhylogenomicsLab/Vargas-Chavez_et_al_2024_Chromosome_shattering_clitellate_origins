{bash}
#!/bin/sh

#SBATCH -N 1
#SBATCH -C clk
#SBATCH --mem=2G
#SBATCH -o OMA_expanded_db.%A_%a.out
#SBATCH -e OMA_expanded_db.%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user="mail@.com"


#Running db single process
echo "***************" Running OMA first stage "***************"
OMA -c
echo "***************" db process OMA done "***************"