#!/bin/bash

## VERSION DATE: 10-08-2011

# USAGE:
#   bash WIMT_batch.sh Subj101
#      will process a single subject with the subject number "Subj101"
#   bash WIMT_batch.sh
#      will process all the sujbects in your "ANTS normalized anatomy" directory.

# USER DEFINED VARIABLES

# Location of your ANTS normalized anatomy subjects directory:
wimtBaseDir='/mindhive/gablab/GATES/Analysis/ANTS/normbrains'

# Where do you want to write the output to?
wimtOutDir='/mindhive/gablab/GATES/Analysis/ANTS/norm_Nback_letters'

# Location of your pipeline level 1 output directory (top level, with individual subjects)
subjBaseDir='/mindhive/gablab/GATES/Analysis/Nback_letters/l1output/'

# Template for the directory structure containing your con images and spmT images:
# e.g., for con and spmT images in "$subjBaseDir/subj100/reg_cons/" use: 'reg_cons'
# NB: These images *MUST* be coregistered to the subject anatomy you input to ANTS_batch.sh.
conInTempDir='reg_cons'

# If yourcoregistered con / spmT images have a postfix, what is it?
# (e.g., for "con_0001_warped.img" the postfix is '_warped'
#        for "con_0001.img" the postfix is '')
postfix='_out_warped'

# What version of Freesurfer do you want to use?
fsversion='5.1.0'

# Do you want to overwrite existing output files?
# 0 = "no", 1 = "always", 2 = "only if normalized image is older than input image"
overwrite=1



###############################################################################
# DON'T CHANGE ANYTHING BELOW UNLESS YOU ARE WILLING TO ACCEPT THE CONSEQUENCES:

# ADD Freesurfer AND Torque TO YOUR ENVIRONMENT
source /software/common/bin/fss $fsversion
source /etc/profile.d/modules.sh

# ADD ANTS TO YOUR PATH
export PATH=/software/ANTS-dev:$PATH
export ANTSPATH=/software/ANTS-dev/

# CONFIRM TORQUE IS RUNNING
currentUser=`whoami`
module list &> /tmp/${currentUser}.ants.torquecheck.tmp
checkForTorque=`cat /tmp/${currentUser}.ants.torquecheck.tmp | grep -c torque`
if [ $checkForTorque -lt 1 ]
then
	echo
	echo "Cannot run WIMT normalization process."
	echo "This script requires the ability to submit jobs to the PBS Torque queue."
	echo "To activate the torque queue, input:"
	echo -e "\t\e[1;39mmodule add torque\e[0m"
	echo "from a terminal."
	echo
	exit
fi

if [ $# -lt 1 ]
then
	subjlist=`ls -F ${wimtBaseDir} | grep \/ | sed 's/\///'`
else
	subjlist=`echo $@`
fi

echo "Preparing to normalize the following subjects:"
echo ${subjlist}


for subj in  $subjlist
do
	conBaseDir="$subjBaseDir"/"$subj"/"$conInTempDir"
	outDir="$wimtOutDir"/"$subj"/

	if [ ! -d $conBaseDir ]
	then
		echo -e "\n\nNOTICE: Cannot find pipeline output directory!"
		echo "Subject is: ${subj}"
		echo "Skipping..."
		echo `date` >> wimt.error.log
		echo -e "Cannot find pipeline output directory.\nSubject is: ${subj}\n" >> wimt.error.log
		continue
	fi

	if [ ! -d $outDir ]
	then
		cmd="mkdir -p ${outDir}"
		#echo $cmd
		eval $cmd
	fi

	numContrasts=`ls ${conBaseDir}/*.nii* | grep -c con_`
	#echo $numContrasts

	for contrast in `seq 1 ${numContrasts}`
	do

		for imgType in con #spmT var
		do
			cvtInImg="${conBaseDir}/${imgType}_`printf %04d $contrast`${postfix}.nii.gz"
			wimtOutImg="${outDir}/${imgType}_`printf %04d $contrast`${postfix}_ants.nii"

			if [ ! -e $cvtInImg ]
			then
				echo -e "\n\nNOTICE: Cannot find spm con image!"
				echo "Subject is: ${subj}"
				echo "Target file was: ${cvtInImg}"
				echo "Skipping..."
				echo -e "Cannot find spm con image.\nSubject is: ${subj}\nTarget file was: ${cvtInImg}\n" >> wimt.error.log
				continue
			fi

			if [ ! -e $wimtBaseDir/$subj/ants_deformed.nii.gz ] || [ ! -e $wimtBaseDir/$subj/ants_Warp.nii.gz ] || [ ! -e $wimtBaseDir/$subj/ants_Affine.txt ]
			then
				echo -e "\n\nNOTICE: A necessary ANTS output file was not found!"
				echo "${wimtBaseDir}/${subj}"
				echo "WIMT script requires deformed.nii.gz, Warp.nii.gz, and Affine.txt"
				echo "Subject is: ${subj}"
				echo "Skipping..."
				echo -e "Cannot find ANTS normalized anatomy.\nSubject is: ${subj}\n" >> wimt.error.log
				continue
			fi

		### CHECK TO SEE IF OUTPUT FILES ALREADY EXIST AND/OR RUN:

			if [ ! -e $wimtOutImg ] || [ $overwrite -ne 0 ]
			then
				if [ -e $wimtOutImg ] && [ $cvtInImg -ot $wimtOutImg ] && [ $overwrite -eq 2 ]
				then
					echo -e "Output file ${wimtOutImg} is older than input file ${cvtInImg} ... skipping.\n(Specify overwrite=1 in the script to force overwrite)."
				else
					echo -ne "Subject:  ${subj}\tContrast:  ${contrast}\t"
				#	echo "${cvtInImg} -> ${wimtOutImg}"
					cmd="WarpImageMultiTransform 3 \
						${cvtInImg} \
						${wimtOutImg} \
						-R ${wimtBaseDir}/${subj}/ants_deformed.nii.gz \
						${wimtBaseDir}/${subj}/ants_Warp.nii.gz \
						${wimtBaseDir}/${subj}/ants_Affine.txt"
					echo "$cmd" > /tmp/${currentUser}.wimt.${subj}.qsub
					#eval "$cmd"
					qsub -V -e /dev/null -o /dev/null -N WIMT.${subj} /tmp/${currentUser}.wimt.${subj}.qsub
					read -t 0.25  # wait a bit to avoid overwhelming the queue with job submissions
					rm -f /tmp/${currentUser}.wimt.${subj}.qsub # clean up your mess
				fi
			else
				echo -e "File:\n${wimtOutImg}\nalready exists ... skipping.\n(Specify overwrite=1 in the script to force overwrite)."
			fi
		done #con/spmt
	done #contrast num
done #subj
