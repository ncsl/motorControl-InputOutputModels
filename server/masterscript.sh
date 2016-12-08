printf "run sleep (Enter 1, or 0)? "
read RUNSLEEP

if [[ "$RUNSLEEP" -eq 1 ]]; then
	# runs the sleep function on all faulty nodes 
	qsub -l walltime=24:00:00,nodes=node054 run_b_sleep.sh
	qsub -l walltime=24:00:00,nodes=node215 run_b_sleep.sh
	qsub -l walltime=24:00:00,nodes=node232 run_b_sleep.sh
fi

for i in {1..40}; do
	echo $i;
	# run a pbs batch job. Make sure there are no spaces in between the parameters passed
	qsub -v index=$index run_job.pbs
done