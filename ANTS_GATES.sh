#!/bin/bash

# USAGE:
#   bash ANTS_batch.sh Subj101
#      will process a single subject with the subject number "Subj101"
#   bash ANTS_batch.sh
#      will process all the sujbects in your Freesurfer subjects directory.

# USER DEFINED VARIABLES

# Location of your target template:
# (e.g., a standard brain anatomy like the MNI template):
targAnat='/software/fsl/fsl-4.1.6/data/standard/MNI152_T1_1mm_brain.nii.gz'

# Location of your output (target) directory:
antsOutDir='/mindhive/gablab/GATES/Analysis/ANTS/normbrains'

# Location of your Freesurfer subjects directory:
subjBaseDir='/mindhive/gablab/GATES/data'

# What version of Freesurfer do you want to use?
fsversion='5.1.0'

# Do you want to overwrite existing output files?
# 0 = "no", 1 = "always", 2 = "only if normalized image is older than input image"
overwrite=2




###############################################################################
# DON'T CHANGE ANYTHING BELOW UNLESS YOU ARE WILLING TO ACCEPT THE CONSEQUENCES:

source /software/common/bin/fss $fsversion
source /etc/profile.d/modules.sh
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/local/lib
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/software/local/lib"
# ADD ANTS TO YOUR PATH
export PATH=/software/ANTS-dev/bin:$PATH
export ANTSPATH=/software/ANTS-dev/bin/

#export PATH=/software/ANTS-dev:$PATH
#export ANTSPATH=/software/ANTS-dev/

if [ ! -d $antsOutDir ]
then
	mkdir -p ${antsOutDir}
fi

if [ -e "$targAnat" ]
then
	templateName=`echo ${targAnat} | awk -F / '{print $NF'}`
	ln -s "$targAnat" "$antsOutDir"/"$templateName" 2>/dev/null
else
	echo "User specified target anatomy does not exist! Exiting..."
	exit
fi

if [ $# -lt 1 ]
then
	subjlist=`ls -F ${subjBaseDir} | grep \/ | sed 's/\///'`
else
	subjlist=`echo $@`
fi

echo "Preparing to normalize the following subjects:"
echo ${subjlist}
#echo ${subjectList}

for subj in $subjlist
#for subj in $subjectList
do
	if [ ! -d $antsOutDir/$subj ]
	then
		cmd="mkdir ${antsOutDir}/${subj}"
		echo $cmd
		eval $cmd
	fi

	cvtInImg="${subjBaseDir}/${subj}/mri/brain.mgz"
	cvtOutImg="${antsOutDir}/${subj}/brain.nii"

	if [ ! -e $cvtInImg ]
	then
		echo -e "\n\nNOTICE: Cannot find original subject anatomy image!"
		echo "Subject is: ${subj}"
		echo "Target file was: ${cvtInImg}"
		echo "Skipping..."
		echo `date` >> ants.error.log
		echo -e "Cannot find original subject anatomy image.\nSubject is: ${subj}\nTarget file was: ${cvtInImg}\n" >> ants.error.log
		continue
	fi

	if [ ! -e $cvtOutImg ] || [ $overwrite == 1 ]
	then
		cmd="mri_convert ${cvtInImg} ${cvtOutImg}"
		echo $cmd
		eval $cmd
	else
		echo `date` >> ants.error.log
		echo -e "File:\n${cvtOutImg}\nalready exists! ... skipping.\n(Specify overwrite=1 in the script to force overwrite)."
		echo -e "File:\n${cvtOutImg}\nalready exists! ... skipping.\n(Specify overwrite=1 in the script to force overwrite)." >> ants.error.log
	fi

	if [ ! -e $antsOutDir/$subj/ants_deformed.nii.gz ] || [ $overwrite == 1 ]
	then
		cmd="bash ${ANTSPATH}/antsIntroduction.sh \
				-d 3 \
				-r ${antsOutDir}/${templateName} \
				-i ${antsOutDir}/${subj}/brain.nii \
				-o ${antsOutDir}/${subj}/ants_"
		echo $cmd
		eval $cmd
	else
		echo `date` >> ants.error.log
		echo -e "File:\n${antsOutDir}/${subj}/ants_deformed.nii\nalready exists! ... skipping.\n(Specify overwrite=1 in the script to force overwrite)."
		echo -e "File:\n${antsOutDir}/${subj}/ants_deformed.nii\nalready exists! ... skipping.\n(Specify overwrite=1 in the script to force overwrite)." >> ants.error.log
	fi

	inverseName=`echo ${templateName} | awk -F . '{print $1'}`

	if [ -e ${antsOutDir}/${inverseName}_InverseWarp.nii.gz ]
	then
		cmd="mv ${antsOutDir}/${inverseName}_InverseWarp.nii.gz ${antsOutDir}/${subj}/${inverseName}_InverseWarp.nii.gz"
		echo $cmd
		eval $cmd
	fi

done
