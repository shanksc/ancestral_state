#!/bin/bash
for n in {1..22}; do
	srun --time=00:20:00 python3 scripts/check_windows.py working_rep/hg38_chr${n}.fa epo/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${n}.fa check_bases/${n}_window > check_bases/${n}_window.txt & 
done
wait
