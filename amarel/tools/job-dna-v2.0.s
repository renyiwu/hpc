#!/bin/bash

#SBATCH --job-name=072681	# xx month xx day x year x job.
#SBATCH --ntasks=1		# set defaut 1
#SBATCH --cpus-per-task=16	# set to the number of fastq files. max 32 for konglab (partition p_kongt_1 or node hal0114)
#SBATCH --mem=64GB		# defaut 64. max 192
#SBATCH --time=24:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
#SBATCH --partition=p_kongt_1		# mask out with ## if to use the "main" partition.


module use /projects/community/modulefiles
module load nextflow 
export NXF_OPTS='-Xms1g -Xmx4g'
export NF_Work_Dir="/scratch/${USER}/NFWorkDir/${PWD}/work"
mkdir -p $NF_Work_Dir



### 1
IS_SINGLE_DEF=""  # "true" for single end  and "false" for paired end

## Note: for nextflow follow the * and {} combination for paired end files
## may use only * for single end files
## MUST in quotes to prevent bash expansion.

### 2
IN_FILES_DEF="" # Eg, "/path/to/*_R{1,2}_001.fastq.gz" 

echo "CML: $0 $@ "

### 3 Capture arguments from CML

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -s|--se)
    IS_SINGLE_CML="true"
    shift # past argument
    ;;
    -p|--pe)
    IS_SINGLE_CML="false"
    shift # past argument
    ;;
    -i|--input)
    IN_FILES_CML="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done

# Get setting for endingness
if [ "$IS_SINGLE_CML" ]
 then
  export IS_SINGLE="$IS_SINGLE_CML"
 else
  export IS_SINGLE="$IS_SINGLE_DEF"
fi

# Get setting for fatstq locations.
if [ "$IN_FILES_CML" ]
 then
  export IN_FILES="$IN_FILES_CML"
 else
  export IN_FILES="$IN_FILES_DEF"
fi


if [ -z "$IS_SINGLE" ] || [ -z "$IN_FILES" ]; then
 echo '
Input error!
Usage: sbatch job-xxx.s -s|p -i <string>

-s, --se
	Single end reads.
-p, --pe
	Paired end reads.
-i, --input <string>
	<String> must be quoted to prevent bash expansion.
	eg: "/path/to/fastq/*_foo_R{1,2}_bar.fastq.gz"  This format
	is required for PE but may use "*.gz" for SE reads.
	
Or, define the following pamameters in the job file.
IS_SINGLE_DEF=""	  # "true" for single end  and "false" for paired end
IN_FILES_DEF=""		 # Eg, "/path/to/*_R{1,2}_001.fastq.gz" 
'
exit 1
fi

#

mkdir -p $NF_Work_Dir
date
SECONDS=0

nextflow run methyl-seq.nf -w $NF_Work_Dir --SingleEnd="$IS_SINGLE" --reads="$IN_FILES" -with-trace -with-report DNA-methyl.html  -with-timeline DNA-methyl-timeline.html -resume

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



# Delete intermediate files.
if  [ $?  -eq  0 ];
then
    echo " Pipeline completed. Removing WorkDir files"
    rm -rf $NF_Work_Dir
else
    echo "Pipeline not completed." 
fi

