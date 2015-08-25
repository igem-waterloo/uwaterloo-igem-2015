#!/bin/bash
RESULTS=~/working/results
SCRIPTS=~/working/scripts
LOGGING=~/working/results/batch/logs

usage() {
cat << EOF

usage $0 options

This script runs dock_variants using with logging. Intended for submitting jobs
on the cluster.

Example: batch_dock.sh -t 10 -l wt_64 -d ~/working/pdb/wt_64 -p C_B

OPTIONS:
  -h  Show this message
  -t  Number of threads (Mandatory)
  -l  Label for output txt files (Mandatory)
  -d  Directory containing PDB files (Mandatory)
  -f  String passed along with --setup_foldtree to dock_variants.py (optional)
  -p  String passed along with --set_partners to dock_variants.py (optional)
  -s  Passes string "--pam64" to dock_variants.py (optional)
  -a  Passes string "--alt_tool" to dock_variants.py (optional)
EOF
}
threads=
label=
pdb_dir=
set_partners=
setup_foldtree=
pam64=
pam_tool=
pam_start=0
pam_end=255

while getopts “ht:l:d:f:p:sa” OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        t)
            threads=$OPTARG
            ;;
        l)
            label=$OPTARG
            ;;
        d)
            pdb_dir=$OPTARG
            ;;
        f) 
            setup_foldtree=$OPTARG
            ;;
        p)
            set_partners=$OPTARG
            ;;
        s)
            pam64=True
            ;;
        a)
            pam_tool=True
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

# check mandatory arguments have values
if [ -z $threads ] || [ -z $label ]  || [ -z $pdb_dir ]; then
    echo "ERROR: requires three arguments (threads, label, pdb_dir)"
    usage
    exit 1
fi

# parse optional arguments for dock_variants.py and csv_results.py
opt_dock_variants_args=""
opt_csv_results_args=""
if [ ! -z $setup_foldtree ]; then
    echo "Adding '--setup_foldtree $setup_foldtree' argument to dock_variants"
    opt_dock_variants_args+=" --setup_foldtree $setup_foldtree"
fi
if [ ! -z $set_partners ]; then
    echo "Adding '--set-partners $set_partners' argument to dock_variants"
    opt_dock_variants_args+=" --set_partners $set_partners"
fi
if [ ! -z $pam64 ]; then
    echo "Passing '--pam64' argument to dock_variants"
    opt_dock_variants_args+=" --pam64"
    pam_end=63
fi
if [ ! -z $pam_tool ]; then
    echo "Adding '--alt_tool' argument to dock_variants"
    echo "Adding '--alt_tool' argument to csv_results"
    opt_dock_variants_args+=" --alt_tool"
    opt_csv_results_args=" --alt_tool"
fi

# don't write to a directory that already exists
results_dir=$RESULTS/batch/$label
if [ -d $results_dir ]; then
    echo "Directory $results_dir already exists, specify a new label"
    exit 1
fi

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
    log_name="$LOGGING/$label"
    log_name+="_pf_$pam_first"
    log_name+="_pl_$pam_last"
    log_name+="_thread_$i" # name for logging stdout and stderr

    out_name=$log_name
    out_name+="_out.txt"

    err_name=$log_name
    err_name+="_err.txt"

    nohup nice -n 10 python $SCRIPTS/dock_variants.py -s=$pam_first -e=$pam_last --output_dir=$results_dir --pdb_dir=$pdb_dir $opt_dock_variants_args > $out_name 2> $err_name &
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
python $SCRIPTS/csv_results.py --results_dir $results_dir $opt_csv_results_args

echo "Done!"
