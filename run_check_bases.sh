#!/bin/bash
echo "check bases" > check_bases.txt
for n in {1..22}; do
	echo ${n} >> check_bases.txt 	
	python3 scripts/check_bases.py working/hg38_chr${n}.fa epo/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${n}.fa >> check_bases.txt 
done
