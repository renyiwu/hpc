#!/usr/bin/bash
echo -e "\n##1. pwd ##"
pwd

echo -e "\n##2. sinfo ##"
sinfo

echo -e  "\n##3. squeue -u rw409 ##"
squeue -u rw409

echo -e "\n##4. sacct ## use scancel <job_id> to cancel/terminate a job"
sacct

echo -e "\n##5. scontrol show node hal0114##"
scontrol show node hal0114

echo -e "\n##6. scontrol show partition p_kongt_1 ##"
scontrol show partition p_kongt_1

date
pwd
