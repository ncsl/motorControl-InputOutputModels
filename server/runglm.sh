#!/bin/bash -l

source /etc/profile.modules
module load matlab/matlab2013a

# 0. Load in input parameters
proc="$1"
index="$2"

## 01: Set parameters for matlab to run, and check if matlab is on path
matlab_jvm="matlab -nojvm -nodesktop -nosplash -r"
[[ ! -z "`which matlab`" ]] || \
	{ 
		echo "MATLAB not found on the PATH; please add to path."; 
		exit 1;
	}

echo "Running glm computation."
echo $proc
echo $index
matlab -logfile /home/ali/motorControl-InputOutputModels/server/_log/job$1.txt -nojvm -nodisplay -nosplash -r "index='$index'; \
	runppm(index);\
	exit"