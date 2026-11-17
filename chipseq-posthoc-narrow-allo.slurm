#!/bin/bash
#SBATCH --job-name=posthoc_chipseq_narrow_allo
#SBATCH --partition=bigmem
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=posthoc_chip_narrow_allo.out
#SBATCH --error=posthoc_chip_narrow_allo.err
#SBATCH --mail-type=END
#SBATCH --mail-user=bradlem9@ccf.org

RUN_CONFIG=config/config_human_narrow_allo.yaml
FINAL_RES_DIR=/home/bradlem9/isilon/cGAS/chipseq_crossrun/crossrun_out_human_narrow_allo

cd /home/bradlem9/bradlem9/chipseq_crossrun
bash run_all.sh $RUN_CONFIG

echo "Syncing results from $(basename $FINAL_RES_DIR) to Isilon ($(dirname $FINAL_RES_DIR))"
#rsync -qavI --remove-source-files $(basename $FINAL_RES_DIR) $(dirname $FINAL_RES_DIR)/
rsync -qavI $(basename $FINAL_RES_DIR) $(dirname $FINAL_RES_DIR)/
# Delete directory structure after transfer
#find $(basename $FINAL_RES_DIR) -depth -type d -empty -delete
cp $RUN_CONFIG $FINAL_RES_DIR
