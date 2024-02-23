#!/bin/bash

sum_A=0
sum_C=0
sum_G=0
sum_T=0
sum_periods=0

for file in ${1}/hg38_chr*_stats.txt; do
    counts=$(grep -oE 'A: [0-9]+ C: [0-9]+ G: [0-9]+ T: [0-9]+ Periods: [0-9]+' "$file")

    count_A=$(echo "$counts" | awk '{print $2}')
    count_C=$(echo "$counts" | awk '{print $4}')
    count_G=$(echo "$counts" | awk '{print $6}')
    count_T=$(echo "$counts" | awk '{print $8}')
    count_periods=$(echo "$counts" | awk '{print $10}')

    sum_A=$((sum_A + count_A))
    sum_C=$((sum_C + count_C))
    sum_G=$((sum_G + count_G))
    sum_T=$((sum_T + count_T))
    sum_periods=$((sum_periods + count_periods))
done

echo "Sum of A: $sum_A"
echo "Sum of C: $sum_C"
echo "Sum of G: $sum_G"
echo "Sum of T: $sum_T"
echo "Sum of Periods: $sum_periods"

