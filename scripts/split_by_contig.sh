#!/bin/bash
ref=$1
taf=$2
out_dir=$3

#echo ${ref}
#echo ${taf}

#rm ${out_dir}/${ref}.done.txt

echo "creating tafs of each chr"
for n in {1..22}
do
	echo "${ref} ${n}"
	taffy view -c --region ${ref}.chr${n} -i $taf > ${out_dir}/${ref}.${n}.taf.gz &
done
wait

echo "indexing tafs"
for n in {1..22}
do
	taffy index -i ${out_dir}/${ref}.${n}.taf.gz &
done
wait

#if this works we write this
echo "done" > ${out_dir}/${ref}.done.txt

