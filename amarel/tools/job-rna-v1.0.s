#!/bin/bash

#SBATCH --job-name=072681	# xx month xx day x year x job.
#SBATCH --ntasks=1		# set defaut 1
#SBATCH --cpus-per-task=2	# set to the number of fastq files. max 32 for konglab (partition p_kongt_1 or node hal0114)
#SBATCH --mem=16GB		# defaut 64. max 192
#SBATCH --time=2:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
#SBATCH --partition=p_kongt_1		# mask out with ## if to use the "main" partition.




module use /projects/community/modulefiles
module load nextflow 
export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
export NXF_OPTS='-Xms1g -Xmx4g'
export NF_Work_Dir="/scratch/${USER}/NFWorkDir/${PWD}/work"


### 1
export IS_SINGLE="false"  # "true" for single end  and "false" for paired end

## Note: for nextflow follow the * and {} combination for paired end files
## may use only * for single end files
## MUST in quotes to prevent bash expansion.

### 2
export IN_FILES="/scratch/rw409/oarc/pten/link/*_R{1,2}_001.fastq.gz" 






mkdir -p $NF_Work_Dir
date
SECONDS=0


srun nextflow run rna-seq.nf --SingleEnd=$IS_SINGLE --reads=$IN_FILES -w $NF_Work_Dir -with-trace -with-report ${SLURM_JOB_PARTITION}_${SLURM_JOB_NODELIST}_${SLURM_JOB_ID}_RNA-seq-nf.html \
-with-timeline ${SLURM_JOB_PARTITION}_${SLURM_JOB_NODELIST}_${SLURM_JOB_ID}_RNA-seq-nf-timeline.html -resume 

date


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
touch "${SLURM_JOB_PARTITION}-${SLURM_JOB_NODELIST}-${SLURM_JOB_ID}-$(($duration / 3600 ))h_$(($(($duration % 3600)) / 60 ))m_$(($duration % 60))s.time"


# Comments by Wu
# Environment variables
# These are the most basic; there are many more.  By default SLURM changes to the directory from which the job was submitted, so the SLURM_SUBMIT_DIR environment variable is usually not needed.
# SLURM_JOB_ID
# SLURM_SUBMIT_DIR
# SLURM_JOB_PARTITION
# SLURM_JOB_NODELIST
