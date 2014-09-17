#--- GA_metajob.sh -a script that sets up the GA simulations and then submits them.
# ---last updated on  Sun Apr 27 14:57:27 CEST 2014  by  ga79moz  at location  TUM , lagrange2

#---------------------------------------------------------------------------------------------------------------------------

NGtype="HNG" #----- command line.
#--------------------------------

kS_N=0.0
kA_N=1.0
kS_TF=0.0
kA_TF=1.0

Llim="100"
t0="0.0000000001"
tf="800"
t_trans=100.0
dt_obs=1.0
dtau_plot=0.1
max_tcorr=1; 

a="2"  

mu_TF=0.0

krm_b=0
krm_val="0.0"

should_plot_snapshots=0 
Nplots2makeshort=0
Nplots2makelong=0

should_plot_kymo=0
Nplots2make_kymo=0


BZcond="boltzmann_on_add"
BZalpha=1.0

parity_check="88885888"
num_trials=100

output_subfolder="muN_scan"

#run_time="72:59:0"
run_time="1:59:0"
		# ------------------pathout="./"
		# ----- we don't specify path_out here. path_out="./"

job_sub_script=GA_Hard_Dimers_muscan.sh
WORKDIR=./GA_output_on_space/job_GA_hard_dimers_muscan_a-${a}_L-${Llim}_N-${num_trials}_output-${NGtype}/

#----------------------------------------------------------------------------------------

eps_array_file="/home/t30/ger/ga79moz/grad_research_phd/project/Nucl/scan_space_eps_muN/array_eps_a"${a}".txt"
muN_array_file="/home/t30/ger/ga79moz/grad_research_phd/project/Nucl/scan_space_eps_muN/array_muN_a"${a}".txt"

#---- FOR THE THERMO-SCANS WE USE A SPECIFIC FILE SET
# eps_array_file="/home/t30/ger/ga79moz/grad_research_phd/project/Nucl/scan_space_eps_muN/array_eps_a"${a}"_THERMO.txt"
# muN_array_file="/home/t30/ger/ga79moz/grad_research_phd/project/Nucl/scan_space_eps_muN/array_muN_a"${a}"_THERMO.txt"



if [ -f ${eps_array_file} ];
then
   echo "reading epsilon array from file  ${eps_array_file} exists."
else
   echo "ERROR: File ${eps_array_file} does not exist. Exiting"
   exit 
fi


if [ -f ${muN_array_file} ];
then
   echo "reading epsilon array from file  ${muN_array_file} exists."
else
   echo "ERROR: File ${muN_array_file} does not exist. Exiting"
   exit
fi


num_eps=$(cat $eps_array_file |wc -l)
num_mus=$(cat $muN_array_file |wc -l)

num_tasks=$(( ${num_eps}*${num_mus} ))

echo "num_eps=" ${num_eps} ", num_mus=" ${num_mus}
echo "preparing to submit num_tasks=" ${num_tasks}

#-----------------SET UP THE C++ INPUT FILE ---------------------------

echo $kS_N  $kA_N  $kS_TF $kA_TF  >  GA_GC_NTF.in
echo $Llim                        >> GA_GC_NTF.in

echo $t0 $tf  $t_trans $dt_obs $dtau_plot $max_tcorr >> GA_GC_NTF.in
echo $a                >> GA_GC_NTF.in
echo $mu_TF            >> GA_GC_NTF.in

echo ""                >>  GA_GC_NTF.in

echo $krm_b  $krm_val  >> GA_GC_NTF.in

echo $should_plot_snapshots  $Nplots2makeshort $Nplots2makelong >> GA_GC_NTF.in
echo $should_plot_kymo   $Nplots2make_kymo                      >> GA_GC_NTF.in

echo $BZcond     $BZalpha                                       >> GA_GC_NTF.in

echo ""                                >>  GA_GC_NTF.in

echo  $output_subfolder                >>  GA_GC_NTF.in
echo  $parity_check                    >>  GA_GC_NTF.in
echo  $num_trials                      >>  GA_GC_NTF.in
# ------echo  $pathout                         >>  GA_GC_NTF.in

#---------------------------------------------------------------------

# create a local directory for this job only
if [ ! -e $WORKDIR ]; then
    mkdir -p $WORKDIR
else
    # this should never happen
#   echo "Clean up /data/$USER/ directory on $HOSTNAME"
#        | mail -s "SGE-error" $USER
    exit 1
fi

# create a local copy of the program and start the job
OLDDIR=`pwd`

#---------COPY NECESSARY FILES OVER TO THE WORK DIRECTORY-------
cp  GA_GC_NTF.x                   ${WORKDIR}/
cp  GA_GC_NTF.in                  ${WORKDIR}/
echo $RANDOM  >                   ${WORKDIR}/rngSEED.in

#-----------------  NOW BUILD THE SCRIPT FILE  -----------------

echo "#!/bin/bash"                      >   ${WORKDIR}${job_sub_script}
echo "#$ -S /bin/sh"                    >>  ${WORKDIR}${job_sub_script}

echo "#$ -cwd"                          >>  ${WORKDIR}${job_sub_script}
# echo "#$ -m eas"                       >>  ${WORKDIR}${job_sub_script}
echo "#$ -l h_vmem=500M,s_rt=${run_time}"   >>  ${WORKDIR}${job_sub_script}
echo "#$ -t 1-"${num_tasks}             >>  ${WORKDIR}${job_sub_script}
# echo "#$ -q lagrange"                   >>  ${WORKDIR}${job_sub_script}

echo ""                                 >>  ${WORKDIR}${job_sub_script}
echo "#-----------------------------"   >>  ${WORKDIR}${job_sub_script}
echo ""                                 >>  ${WORKDIR}${job_sub_script}
echo ""                                 >>  ${WORKDIR}${job_sub_script}

echo "NGtype=\""${NGtype}"\""           >>  ${WORKDIR}${job_sub_script}
echo "eps_array_file="${eps_array_file} >>  ${WORKDIR}${job_sub_script}
echo "muN_array_file="${muN_array_file} >>  ${WORKDIR}${job_sub_script}
echo "num_tasks="${num_tasks}           >>  ${WORKDIR}${job_sub_script}

echo  "num_eps="${num_eps}              >>  ${WORKDIR}${job_sub_script}
echo  "num_mus="${num_mus}              >>  ${WORKDIR}${job_sub_script}

cat GA_jobsub_tail.sh                   >>  ${WORKDIR}${job_sub_script}

#-----------GO TO THE WORK DIRECTORY AND EXECUTE ----------------
cd $WORKDIR

# qsub ${job_sub_script}

#--------------------------------------------

