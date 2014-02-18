#!/bin/bash
#$ -cwd
#$ -l h_vmem=1G,s_rt=18:0:0
#$ -t 1-3
#$ -m eas
#
# The lines above are instructions to SGE.
# 
# Just copy this sample script to a directory where you have write
# permission (e. g., a sub-directory of your home-directory).
# Then edit it to suit your needs.
#
# Start batch jobs via:
# qsub ./myjob.sh

NGtype="HNG"

WORKDIR=./NAR_longpattern_HNG155_${NGtype}/job.${SGE_TASK_ID}

# create a local directory for this job only
if [ ! -e $WORKDIR ]; then
    mkdir -p $WORKDIR 
else
    # this should never happen
    echo "Clean up /data/$USER/ directory on $HOSTNAME"\
        | mail -s "SGE-error" $USER
    exit 1
fi

# create a local copy of the program and start the job
OLDDIR=`pwd`

#---------COPY NECESSARY FILES OVER TO THE WORK DIRECTORY-------
cp  GA_GC_NTF.in $WORKDIR/
cat rngSEED.in > $WORKDIR/rngSEED.in
cp  GA_GC_NTF.x  $WORKDIR/
cp  mu_v_irho_folder/mu_v_irho*   $WORKDIR/
#-----------GO TO THE WORK DIRECTORY AND EXECUTE ----------------
cd $WORKDIR
./GA_GC_NTF.x ./  ${SGE_TASK_ID} ${NGtype}  > GA_GC_NTF.stdout 2> GA_GC_NTF.stderr

# copy all output files back to your home directory
# and clean up


rm mu_v_irho?NG*
#------get rid of the interpolation files first, they're just clutter

# if [ ! -e $OLDDIR/job_${JOB_ID} ]; then
# mkdir $OLDDIR/job_${JOB_ID}
# fi
# cp -R $WORKDIR/ $OLDDIR/job_${JOB_ID} && rm -r $WORKDIR
# cp $OLDDIR/job_${JOB_ID}/avg_occ_Ll-* $OLDDIR/job_${JOB_ID}/avg_result_task_$SGE_TASK_ID
