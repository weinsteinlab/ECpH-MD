#!/bin/bash -l

source ~/.bashrc

module load gcc/8.1.0/1
module load cuda/10.1
conda activate openmm_7.4.0

number_of_replicas=$1

CWD=`pwd`

pH_low=$(grep 'pH_low' input_file.py); eval "${pH_low// /}"
pH_step=$(grep 'pH_step' input_file.py); eval "${pH_step// /}"
pH_high=$(grep 'pH_high' input_file.py); eval "${pH_high// /}"; pH_high=$(echo $pH_high - $pH_step | bc)
pH_seq=($(seq $pH_low $pH_step $pH_high))
replicas_per_pH=$(grep 'replicas_per_pH' input_file.py); eval "${replicas_per_pH// /}"
MD_nsteps_replicas=$(grep 'MD_nsteps_replicas' input_file.py); eval "${MD_nsteps_replicas// /}"

unformatted_output_name=$(grep 'output_name' input_file.py)
eval "${unformatted_output_name// /}" # sets the variable output_name in this scope

number_of_log_files=$(find . -name "*${output_name}*.log" | wc -l)
subjob_number=0

# set subjob_number if previous runs exist
if [ $number_of_log_files != 0 ]; then
    # check that each previous replica job ran
    if [ $(($number_of_log_files % $number_of_replicas)) != 0 ]; then exit 1; fi

    number_of_FINISH=$(grep "FINISH" ./propagate_runs/propagate_runs*.log | wc -l)

    # make sure each job completed correctly
    if [ $number_of_FINISH != $number_of_log_files ]; then exit 1; fi

    subjob_number=$(($number_of_log_files / $number_of_replicas))
fi

if [ $subjob_number -eq 0 ]; then
    echo "generating lambda list"
    python3 -u createLambdaList.py 
fi

# propagate_replicas 
for ((j=0; j < $number_of_replicas; j++)); do
    echo "pH:${pH_seq[j]} subjob_number:${subjob_number}"
    srun -N1 --gres=gpu:32g:1 --mem=50G python3 -u run_replica.py ${pH_seq[j]} ${subjob_number} >> ${CWD}/propagate_runs/propagate_runs_pH_${pH_seq[j]}.log & 
   sleep 5
done

wait  

exit
