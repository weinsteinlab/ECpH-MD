#!/bin/sh

number_of_replicas=12
subjobs_per_iteration=1
iterations_per_subjob=2
number_of_subjobs=1 
jobName="dat" # no spaces
partitionName=dcs            #Slurm partition to run job on

# do not edit below this line


# make sure each iteration is able to completely finish
if [ $((number_of_subjobs % subjobs_per_iteration)) != 0 ]; then exit 1; fi

first_subjob=0
numberOfNodes=`expr $number_of_replicas / 6`
swarmNumber_padded=`printf %04d $swarmNumber`

mkdir -p erf_files energies progress propagate_runs simulations submission_logs lambdas 

for (( subjob=1; subjob<=$number_of_subjobs; subjob++ ))
do
  if [ $first_subjob -eq 0 ]; then
    jobSchedulerOutput="$(sbatch -J ./submission_logs/${jobName} -N ${numberOfNodes} -p $partitionName --gres=gpu:32g:6 -C cuda-mode-exclusive -t 0-02:00:00 ./submit_Exchange-min-replica.sh ${number_of_replicas} ${subjobs_per_iteration} ${iterations_per_subjob})"
  else
    jobSchedulerOutput="$(sbatch --depend=afterok:${job_scheduler_number} -J ./submission_logs/${jobName} -N ${numberOfNodes} -p $partitionName --gres=gpu:32g:6 -C cuda-mode-exclusive -t 0-02:00:00 ./submit_Exchange-min-replica.sh ${number_of_replicas} ${subjobs_per_iteration} ${iterations_per_subjob})"
  fi

  job_scheduler_number=$(echo $jobSchedulerOutput | awk '{print $2}' | sed -e 's/<//' | sed -e 's/>//')
  let first_subjob=1
done

exit
