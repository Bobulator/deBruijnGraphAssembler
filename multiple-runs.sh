#!/bin/sh

# $1: File to analyze
# $2: Smallest kmer size
# $3: Largest kmer size
# $4: Coverage filter limit
# $5: Weighted edge filter limit

for ((i=$2; i<$3+1; i++)); do
	for ((j=0; j<$4+1; j++)); do
		for ((k=0; k<$5+1; k++)); do
			echo contigGenerator.py $1 $i $j $k
   			python contigGenerator.py $1 $i $j $k
		done
	done
done