#---carpark_jobsubtail.sh : BELOW HERE IS DIRECTLY COPIED INTO THE SH FILE FROM THE JOBSUB_TAIL.SH TEMPLATE.------
# ---last updated on  Sun Apr 27 14:22:25 CEST 2014  by  bren  at location  , bren-Desktop

#  changes from  Sun Apr 27 14:22:25 CEST 2014 : fixed the modulo submission criteria (-1 from the taskID modulo num_eps)

#----------------------------------------------------------------------

TASKIDm1=$((${SGE_TASK_ID}-1));

muN_i=$(( ${TASKIDm1}/${num_eps} +1))
eps_i=$(( ${TASKIDm1}%${num_eps} +1))

#----------------------------------------------------------------------

muN=$(head -$muN_i  ${muN_array_file} | tail -1)
eps=$(head -$eps_i  ${eps_array_file} | tail -1)

#----------------------------------------------------------------------
# echo " executing job " ${JOB_ID} " with paramaters : " $NGtype $muN $eps >> job_${JOB_ID}.log



D_before=$(date)
echo " executing job " ${JOB_ID}.${SGE_TASK_ID} " with paramaters : " $NGtype $muN $eps  " at " $D_before >> job_${JOB_ID}.${SGE_TASK_ID}.log

 ./GA_GC_NTF.x $NGtype ${SGE_TASK_ID} $muN $eps

D_after=$(date)

echo " Job finished executable at " $D_after ", Now terminating. " >> job_${JOB_ID}.${SGE_TASK_ID}.log

# copy all output files back to your home directory
# and clean up
# mkdir $OLDDIR/job_${JOB_ID}
# cp -R $WORKDIR/ $OLDDIR/job_${JOB_ID} && rm -r $WORKDIR

