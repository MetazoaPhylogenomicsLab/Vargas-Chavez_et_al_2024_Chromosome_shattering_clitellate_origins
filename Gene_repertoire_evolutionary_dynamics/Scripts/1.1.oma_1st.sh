{bash}
#!/bin/sh

#SBATCH -N 1
#SBATCH -C clk
#SBATCH --mem=2G
#SBATCH --array=1-60
#SBATCH -o OMA_1st.%A_%a.out
#SBATCH -e OMA_1st.%A_%a.err
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user="mail@.com"

export NR_PROCESSES=60

#Running all-versus-all 1st stage
echo "**************" Running OMA first stage "************"
OMA -s
echo "***************" First stage OMA done "***************"
