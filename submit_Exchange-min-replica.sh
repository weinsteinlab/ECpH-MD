#!/bin/bash -l

source ~/.bashrc

module load gcc/8.1.0/1
module load cuda/10.1

number_of_replicas=$1
subjobs_per_iteration=$2
iterations_per_subjob=$3

iteration_number=0
subjob_number=0

rm -rf ./progress/* num_finished_jobs.txt
unformatted_output_name=$(grep 'output_name' input_file.py)
eval "${unformatted_output_name// /}" # sets the variable output_name in this scope

echo $number_of_replicas
number_of_log_files=$(ls -1 ./simulations/${output_name}*.log 2>/dev/null | wc -l)


# set subjob_number and iteration_number if previous runs exist
if [ $number_of_log_files != 0 ]; then
    # check that each replica's prior job ran
    if [ $(($number_of_log_files % $number_of_replicas)) != 0 ]; then exit 1; fi
    
    number_of_FINISH=$(grep "FINISH" ./propagate_runs/propagate_runs*.log | wc -l)

    # make sure each job completed correctly
    if [ $number_of_FINISH != $number_of_log_files ]; then exit 1 ; fi

    individual_log_file=$(ls -1 ./propagate_runs/propagate_runs_pH_*.log | head -n1)
    finish_in_individual_log=$(grep "FINISH" $individual_log_file | wc -l) 

    subjob_number=$(($finish_in_individual_log % subjobs_per_iteration))
    iteration_number=$((number_of_log_files/(subjobs_per_iteration * number_of_replicas)))
fi
    
conda activate openmm_7.4.0

python3 -u Exchange-min-replica.py $iteration_number $subjob_number $subjobs_per_iteration $iterations_per_subjob >> overall_job.log
