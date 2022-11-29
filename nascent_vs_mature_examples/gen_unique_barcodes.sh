#!/bin/bash
# Generate 10XV3 R1 read sequences with random barcodes and GGGGGGGG as UMI
# Usage: gen_unique_barcodes.sh <number_of_reads>
array=("A" "T" "C" "G")
end=$(($1/4))
for (( x=0; x<$end; x++ ))
do
	echo ""
	shh=15
	str=""
	for i in {1..16}
	do
		xx=$((($x >> (2*$shh)) & 3))
		chh="${array[$xx]}"
		str=$str$chh
		shh=$(($shh-1))
	done
	echo "$str"GGGGGGGG
	echo ""
	echo ""
done
