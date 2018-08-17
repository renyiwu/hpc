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

echo "CML: $0 $@ "

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -s|--se)
    IS_SINGLE="true"
    shift # past argument
    ;;
    -p|--pe)
    IS_SINGLE="false"
    shift # past argument
    ;;
    -i|--input)
    IN_FILES="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done


if [ -z "$IS_SINGLE" ] || [ -z "$IN_FILES" ]; then
 echo '
Input error!
Usage: sbatch job-xxx.s -[sp] -i <string>

-s, --se
	Single end reads.
-p, --pe
	Paired end reads.
-i, --input <string>
	<String> must be quoted to prevent bash expansion.
	eg: "/path/to/fastq/*_foo_R{1,2}_bar.fastq.gz"  This format
	is required for PE but may use "*.gz" for SE reads.
	
'
exit 1
fi

#

export $IS_SINGLE
export $IN_FILES



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
