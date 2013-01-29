#!/bin/bash
stats_basedir='/mindhive/gablab/GATES/Analysis/ROI/datasink_052312'
 


headerdone=false
hed="Subject,"

for header in 'Index' 'SegId' 'NVoxels' 'Volume_mm3' 'StructName' 'Mean' 'StdDev' 'Min' 'Max' 'Range' ;do
	hed=${hed}$(cat ${stats_basedir}/func_seg_stats/_subject_id_306/_func_segstats0/summary.stats | sed -n '/#/!p' | awk ' { print func_"'$header'_"$5 } ' | tr '\n' ',')
	hed=${hed}$(cat ${stats_basedir}/struct_seg_stats/_subject_id_306/_struct_segstats0/summary.stats | sed -n '/#/!p' | awk ' { print struct_"'$header'_"$5 } ' | tr '\n' ',')
done

echo $hed > ${stats_basedir}/segstats_052312.csv

for subj in 306 307 308 309 310 311 312 313 314 315 317 318 319;do

	echo "working on subj ${subj}..."		
	for contrast in 0 1 2 3 4 5 6 7; do
		declare -i i
		i=1
		subjbody="${subj}-WM-con${contrast},"
		for header in 'Index' 'SegId' 'NVoxels' 'Volume_mm3' 'StructName' 'Mean' 'StdDev' 'Min' 'Max' 'Range' ;do
			subjbody=${subjbody}$(cat ${stats_basedir}/func_seg_stats/_subject_id_${subj}/_func_segstats${contrast}/summary.stats | sed -n '/#/!p' | awk ' { print $'$i' } ' | tr '\n' ',')		
			subjbody=${subjbody}$(cat ${stats_basedir}/struct_seg_stats/_subject_id_${subj}/_struct_segstats${contrast}/summary.stats | sed -n '/#/!p' | awk ' { print $'$i' } ' | tr '\n' ',')
		
			i=$i+1
		done
	echo $subjbody >> ${stats_basedir}/segstats_052312.csv
	done

done
