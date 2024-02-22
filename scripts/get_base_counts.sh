#!/bin/bash
awk '/^>/ {next} {A+=gsub("A", ""); C+=gsub("C", ""); G+=gsub("G", ""); T+=gsub("T", ""); period_count+=gsub("\\.", "")} END {print "A:", A, "C:", C, "G:", G, "T:", T, "Periods:", period_count}' $1
