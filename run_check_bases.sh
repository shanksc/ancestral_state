#!/bin/bash
for n in {1..22}; do
	srun python3 scripts/check_bases.py working_single_copy/hg38_chr${n}.fa epo/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${n}.fa ${n}.npy > check_bases/${n}.txt 
done
