#!/bin/bash
# Generate 10XV3 R1 read sequences with random barcodes and GGGGGGGGGGGG as UMI
# Usage: gen_unique_barcodes.sh <number_of_reads>
# (If a second argument is supplied, a random 12-bp UMI will also be generated rather than GGGGGGGGGGGG)
array=("A" "T" "C" "G")
end=$(($1/4))
for (( x=0; x<$end; x++ ))
do
        echo "@read"$x
	shh=15
	str=""
	for i in {1..16}
	do
		xx=$((($x >> (2*$shh)) & 3))
		chh="${array[$xx]}"
		str=$str$chh
		shh=$(($shh-1))
	done
        if [ -z "$2" ]
        then
                echo "$str"GGGGGGGGGGGG
        else
                echo "$str"${str:4:16}
        fi
        echo "+"
        echo "KKKKKKKKKKKKKKKKKKKKKKKKKKKK"
done
