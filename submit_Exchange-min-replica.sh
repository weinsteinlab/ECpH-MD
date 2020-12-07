#!/bin/bash -l

source ~/.bashrc

module load gcc/8.1.0/1
module load cuda/10.1
conda activate openmm_7.4.0

number_of_replicas=$1
subjobs_per_iteration=$2
iterations_per_subjob=$3

iteration_number=0
subjob_number=0
node_number=1

CWD=`pwd`

pH_low=$(grep 'pH_low' input_file.py); eval "${pH_low// /}"
pH_step=$(grep 'pH_step' input_file.py); eval "${pH_step// /}"
pH_high=$(grep 'pH_high' input_file.py); eval "${pH_high// /}"; pH_high=$(echo $pH_high - $pH_step | bc)
pH_seq=($(seq $pH_low $pH_step $pH_high))
replicas_per_pH=$(grep 'replicas_per_pH' input_file.py); eval "${replicas_per_pH// /}"

MD_nsteps_lambdas=$(grep 'MD_nsteps_lambdas' input_file.py); eval "${MD_nsteps_lambdas// /}"
MD_nsteps_replicas=$(grep 'MD_nsteps_replicas' input_file.py); eval "${MD_nsteps_replicas// /}"

MD_nsteps_lambdas=$(grep 'MD_nsteps_lambdas' input_file.py); eval "${MD_nsteps_lambdas// /}"
MD_nsteps_replicas=$(grep 'MD_nsteps_replicas' input_file.py); eval "${MD_nsteps_replicas// /}"

for ((i=0; i < $iterations_per_subjob; i++)); do
    rm -rf ./progress/* num_finished_jobs.txt
    unformatted_output_name=$(grep 'output_name' input_file.py)
    eval "${unformatted_output_name// /}" # sets the variable output_name in this scope

    number_of_log_files=$(find . -name "*${output_name}*.log" | wc -l)

    # set subjob_number and iteration_number if previous runs exist
    if [ $number_of_log_files != 0 ]; then
        # check that each previous replica job ran
        if [ $(($number_of_log_files % $number_of_replicas)) != 0 ]; then exit 1; fi

        number_of_FINISH=$(grep "FINISH" ./propagate_runs/propagate_runs*.log | wc -l)

        # make sure each job completed correctly
        if [ $number_of_FINISH != $number_of_log_files ]; then exit 1; fi

        individual_log_file=$(ls -1 ./propagate_runs/propagate_runs_pH_*.log | head -n1)
        finish_in_individual_log=$(grep "FINISH" $individual_log_file | wc -l)

        subjob_number=$(($finish_in_individual_log % subjobs_per_iteration))
        iteration_number=$((number_of_log_files/(subjobs_per_iteration * number_of_replicas)))
    fi

    if [ $iteration_number -eq 0 ] && [ $subjob_number -eq 0 ]; then
        echo "generating lambda list"
        python3 -u createLambdaList.py 
    fi

    # propagate_replicas 
    for ((j=0; j < $number_of_replicas; j++)); do
        echo "pH:${pH_seq[j]} iteration_number:${iteration_number} subjobs_per_iteration:${subjobs_per_iteration} MD_nsteps_replica:${MD_nsteps_replicas} MD part 1"
        srun -N1 --gres=gpu:32g:1 --mem=50G python3 -u run_replica.py ${pH_seq[j]} ${iteration_number} ${subjob_number} $subjobs_per_iteration $MD_nsteps_replicas 1 >> ${CWD}/propagate_runs/propagate_runs_pH_${pH_seq[j]}.log & 
       sleep 5
    done

    wait  

    #mix_lambdas
    #jsrun --smpiargs=none  python3 -u mix_lambdas.py ${iteration_number} ${subjob_number} 1
done
