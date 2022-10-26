#!/bin/bash

RepeatMaskerOutfile=$1
outputfile=$2

rm $outputfile

awk 'BEGIN{OFS="\t"}(NR>=4){strand="+"; if($9=="C"){ strand="-" }; split($11,arr,"/");id=arr[2]":"arr[1];id=$5"|"$6"|"$7"|"$10":"id"|"strand; if(arr[1] ~ /DNA|LINE|LTR|RC|SINE/) { print $5,$6,$7,id,$1,strand } }' $RepeatMaskerOutfile > $outputfile
