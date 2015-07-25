#!/bin/bash

# TO DO
# add mandatory input args specifying label and threads
threads=6
#label=
pam_start=14
pam_end=63
range=$(($pam_end-$pam_start+1))
results_dir="results/batch/$label/"

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
	python dock_variants.py -s=$pam_first -e=$pam_last & --output_dir=$results_dir &
	pids+=($!)	
	i=$(($i+1))
	sleep 1
done

for pid in ${pids[*]}; do
	wait $pid
done

# run CSV script
echo "Writing to CSV"
python results_csv.py $results_dir

echo "Done!"