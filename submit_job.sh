#!/bin/bash -l

number_of_replicas=8      # Must be a multiple of # of GPUs per node
number_of_subjobs=4 
subjobs_before_exchange=2 # Set to 0 if no exchanges desired; if not 0, then must be >= 2
jobName="example"         # No spaces
partitionName=edison      # Slurm partition to run job on
number_of_GPUs_per_node=8 # Must be >=2 if running exchanges


# do not edit below this line


first_subjob=0
numberOfNodes=`expr $number_of_replicas / $number_of_GPUs_per_node`
swarmNumber_padded=`printf %04d $swarmNumber`

# make sure number_of_replicas is a multiple of number_of_GPUs_per_node
if [ $(($number_of_replicas % $number_of_GPUs_per_node)) != 0 ]; then exit 1; fi

mkdir -p energies propagate_runs simulations submission_logs lambdas 

for (( subjob=0; subjob<$number_of_subjobs; subjob++ )); do
    jobSchedulerOutput=0

    if [ $first_subjob -eq 0 ]; then
        #jobSchedulerOutput="$(sbatch -J ${jobName} -N ${numberOfNodes} -p $partitionName --gres=gpu:32g:${number_of_GPUs_per_node} -C cuda-mode-exclusive -t 0-02:00:00 ./submit_Exchange-min-replica.sh ${number_of_replicas} 0)"
        jobSchedulerOutput="$(sbatch -J ${jobName} -N ${numberOfNodes} -p $partitionName --mem=160G --gres=gpu:${number_of_GPUs_per_node} -t 0-02:00:00 ./submit_Exchange-min-replica.sh ${number_of_replicas} 0)"

    else
        if [ $subjobs_before_exchange != 0 ] && [ $(($subjob % $subjobs_before_exchange)) == 0 ]; then
            # Perform exchange
            #jobSchedulerOutput="$(sbatch --depend=afterok:${job_scheduler_number} -J ${jobName} -N 1 -p $partitionName --gres=gpu:32g:${number_of_GPUs_per_node} -C cuda-mode-exclusive -t 0-02:00:00 ./submit_Exchange-min-replica.sh ${number_of_replicas} 1)"
            jobSchedulerOutput="$(sbatch --depend=afterok:${job_scheduler_number} -J ${jobName} -N 1 -p $partitionName --mem=160G--gres=gpu:${number_of_GPUs_per_node} -t 0-02:00:00 ./submit_Exchange-min-replica.sh ${number_of_replicas} 1)"
        else 
            #jobSchedulerOutput="$(sbatch --depend=afterok:${job_scheduler_number} -J ${jobName} -N ${numberOfNodes} -p $partitionName --gres=gpu:32g:${number_of_GPUs_per_node} -C cuda-mode-exclusive -t 0-02:00:00 ./submit_Exchange-min-replica.sh ${number_of_replicas} 0)"
            jobSchedulerOutput="$(sbatch --depend=afterok:${job_scheduler_number} -J ${jobName} -N ${numberOfNodes} -p $partitionName --mem=160G --gres=gpu:${number_of_GPUs_per_node} -t 0-02:00:00 ./submit_Exchange-min-replica.sh ${number_of_replicas} 0)"
        fi
    fi

    job_scheduler_number=${jobSchedulerOutput//[!0-9]}
    let first_subjob=1
done

exit
