#!/bin/bash

RESULTS=~/working/results
SCRIPTS=~/working/scripts
LOGGING=~/working/results/batch/log

# require two input arguments
threads=$1
label=$2
pdb_dir=$3

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
	echo "ERROR: requires three arguments (threads, label, pdb_dir)"
	echo "Example: bash batch_dock.sh 50 validation_1 ~/pdb/setA/"
	exit 1
fi

# don't write to a directory which already exists
results_dir=$RESULTS/batch/$label
if [ -d $results_dir ]; then
	echo "Directory $results_dir already exists, specify a new label"
	exit 1
fi

pam_start=0
pam_end=255
range=$(($pam_end-$pam_start+1))
step=$(($range/$threads))
remainder=$(($range % $threads))
pids=()

# amount of remainder to distribute to tasks
remainder_plus=$(($remainder/$threads+1))

# set so pam_first is initially pam_start
pam_last=$(($pam_start-1))

# run docking scripts for ranges in parallel
echo "Launching scripts"
i=0
while [ $i -lt $threads ]; do
	if [ $remainder -lt $(($i+1)) ]; then
		remainder_plus=0
	fi
	pam_first=$(($pam_last+1))
	pam_last=$(($pam_first+$step+$remainder_plus-1))
	log_name="$LOGGING/$(($label))_pf_$(($pam_first))_pl_$(($pam_last))_thread_$(($i))" # name for logging stdout and stderr

	nohup python $SCRIPTS/dock_variants.py -s=$pam_first -e=$pam_last --output_dir=$results_dir --pdb_dir=$pdb_dir > $(($log_name))_out.txt 2> $(($log_name))_err.txt &
	pids+=($!)
	i=$(($i+1))
done

sleep 0.5

for pid in ${pids[*]}; do
	wait $pid
done

sleep 0.5

# run CSV script
echo "Writing to CSV"
python $SCRIPTS/results_csv.py $results_dir

echo "Done!"
