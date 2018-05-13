#!/bin/bash
# Submit multiple matlab jobs at individual timesteps it for
# ke_spectrum.m to evaluate.
for it in `seq 41 10 41`; 
do
   fname="/nethome/gchen/var_PC1_${it}_withLand_test2.mat"
   #fname="/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/var_PC1_${it}.mat"
   if [ ! -f $fname ]; # if file non-exists, then run the code and generate the file
   then
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
#BSUB -J getBigQPCVar 
#BSUB -W 00:12
# End of options

# go to the matlab archive to run the matlab code
cd /nethome/gchen/scripts/matlab/CAMSAM_exp

# start matlab code 
/share/opt/MATLAB/R2016a/bin/matlab -nodesktop -nosplash -nodisplay  << M_PROG
filename='${fname}'
getBigQPCVar(${it},filename);
M_PROG

EOJ
`
   echo "JobID = ${JOB} for timesteps ${it} on `date`"
   fi
   sleep 0 # pause at every few seconds to allow non-overlapping
#   memory usage on the same node for multiple tasks
done
exit
