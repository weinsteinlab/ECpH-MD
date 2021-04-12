#!/bin/sh

number_of_replicas=6 # must be a multiple of 6
number_of_subjobs=5 
subjobs_before_exchange=2 # set to 0 if no exchanges desired
jobName="exchange_test" # no spaces
partitionName=dcs            #Slurm partition to run job on

# do not edit below this line

first_subjob=0
numberOfNodes=`expr $number_of_replicas / 6`
swarmNumber_padded=`printf %04d $swarmNumber`

# make sure number_of_replicas is a multiple of 6
if [ $(($number_of_replicas % 6)) != 0 ]; then exit 1; fi

mkdir -p energies propagate_runs simulations submission_logs lambdas 

for (( subjob=0; subjob<$number_of_subjobs; subjob++ )); do
    if [ $first_subjob -eq 0 ]; then
        jobSchedulerOutput="$(sbatch -J ./submission_logs/${jobName} -N ${numberOfNodes} -p $partitionName --gres=gpu:32g:6 -C cuda-mode-exclusive -t 0-02:00:00 ./submit_Exchange-min-replica.sh ${number_of_replicas} 0)"

    else
        if [ $subjobs_before_exchange != 0 ] && [ $(($subjob % $subjobs_before_exchange)) == 0 ]; then
            # Perform exchange
            jobSchedulerOutput="$(sbatch --depend=afterok:${job_scheduler_number} -J ./submission_logs/${jobName} -N 1 -p $partitionName --gres=gpu:32g:1 -C cuda-mode-exclusive -t 0-02:00:00 ./submit_Exchange-min-replica.sh ${number_of_replicas} 1)"
        else 
            jobSchedulerOutput="$(sbatch --depend=afterok:${job_scheduler_number} -J ./submission_logs/${jobName} -N ${numberOfNodes} -p $partitionName --gres=gpu:32g:6 -C cuda-mode-exclusive -t 0-02:00:00 ./submit_Exchange-min-replica.sh ${number_of_replicas} 0)"
        fi
    fi

    job_scheduler_number=${jobSchedulerOutput//[!0-9]}
    let first_subjob=1
done

exit
