#!/bin/bash -l

# do not edit below this line

number_of_replicas=$(grep 'number_of_replicas' input_file.py); number_of_replicas=${number_of_replicas##*number_of_replicas*=}; number_of_replicas=${number_of_replicas%% \#*}
number_of_subjobs=$(grep 'number_of_subjobs' input_file.py); number_of_subjobs=${number_of_subjobs##*number_of_subjobs*=}; number_of_subjobs=${number_of_subjobs%% \#*}
jobName=$(grep 'jobName' input_file.py); jobName=${jobName##*jobName*=}; jobName=${jobName%% \#*}
partitionName=$(grep 'partitionName' input_file.py); partitionName=${partitionName##*partitionName*=}; partitionName=${partitionName%% \#*}; partitionName=$( echo $partitionName | tr -d '"')
number_of_GPUs_per_node=$(grep 'number_of_GPUs_per_node' input_file.py); number_of_GPUs_per_node=${number_of_GPUs_per_node##*number_of_GPUs_per_node*=}; number_of_GPUs_per_node=${number_of_GPUs_per_node%% \#*}

first_subjob=0
numberOfNodes=`expr $number_of_replicas / $number_of_GPUs_per_node`
swarmNumber_padded=`printf %04d $swarmNumber`

# make sure number_of_replicas is a multiple of number_of_GPUs_per_node
if [ $(($number_of_replicas % $number_of_GPUs_per_node)) != 0 ]; then exit 1; fi

mkdir -p energies propagate_runs simulations submission_logs lambdas 

for (( subjob=0; subjob<$number_of_subjobs; subjob++ )); do
    jobSchedulerOutput=0

    if [ $first_subjob -eq 0 ]; then
        jobSchedulerOutput="$(sbatch -o ./submission_logs/${jobName}.%a.out -J ${jobName} -N ${numberOfNodes} -p $partitionName --mem=80G --gres=gpu:${number_of_GPUs_per_node} -t 0-02:00:00 ./submit-replicas.sh ${number_of_replicas} 0)"

    else
        jobSchedulerOutput="$(sbatch --depend=afterok:${job_scheduler_number} -o ./submission_logs/${jobName}.%a.out  -J ${jobName} -N ${numberOfNodes} -p $partitionName --mem=80G --gres=gpu:${number_of_GPUs_per_node} -t 0-02:00:00 ./submit-replicas.sh ${number_of_replicas})"
    fi

    job_scheduler_number=${jobSchedulerOutput//[!0-9]}
    let first_subjob=1
done

exit
