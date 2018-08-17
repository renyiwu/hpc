#!/bin/bash
#SBATCH --job-name=m-seq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=185GB
#SBATCH --time=240:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)

#
# Use konglab node
#
#SBATCH --partition=p_kongt_1

module use /projects/community/modulefiles
module load nextflow 
export NXF_OPTS='-Xms1g -Xmx4g'
export NF_Work_Dir="/scratch/${USER}/NFWorkDir/${PWD}/work"
mkdir -p $NF_Work_Dir

#
## Version 1.x EXPECTS USER TO DEFINE --SingleEnd AND --reads IN THIS JOB FILE.
#
nextflow run methyl-seq.nf -w $NF_Work_Dir --SingleEnd="false" --reads="/path/to/*_R{1,2}_001.fastq.gz" -with-trace -with-report DNA-methyl.html  -with-timeline DNA-methyl-timeline.html -resume 
if  [ $?  -eq  0 ];
then
    echo " Pipeline completed. Removing WorkDir files"
    rm -rf $NF_Work_Dir
else
    echo "Pipeline not completed." 
fi

