#!/bin/bash
# Submit multiple matlab jobs at individual timesteps 
for it in `seq 2 1 755`; 
do
   JOB=`bsub << EOJ
#!/bin/bash
#==============================================================================
#  Submit matlab batch job for pegasus2
#==============================================================================
#BSUB -n 1
#BSUB -q general
#BSUB -R "rusage[mem=20000]"
#BSUB -o matlabrun.stdout.%J
#BSUB -e matlabrun.stderr.%J
#BSUB -J getQVTForCAPE.m 
#BSUB -W 00:12
# End of options

# go to the matlab archive to run the matlab code
cd /nethome/gchen/scripts/matlab/CAMSAM_exp

# start matlab code 
/share/opt/MATLAB/R2016a/bin/matlab -nodesktop -nosplash -nodisplay  << M_PROG
getQVTForCAPE(${it});
M_PROG

EOJ
`
   echo "JobID = ${JOB} for timesteps ${it} on `date`"
   sleep 0 # pause at every few seconds to allow non-overlapping
#   memory usage on the same node for multiple tasks
done
exit
