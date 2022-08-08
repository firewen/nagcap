#!/bin/bash

evid=$1
dep1=$2
dep2=$3
ddep=$4

if [ -e tmp ]; then
	rm tmp
fi

for ((evdp=dep1;evdp<=dep2; evdp=evdp+ddep))
do
	m=$(cat $evid/output_$evdp.txt | awk '{if($1=="M")print $3,$4,$5,$6,$7,$8}')
	exp=$(cat $evid/output_$evdp.txt | awk '{if($1=="exp") print $3+7}')
	err=$(cat $evid/output_$evdp.txt | awk '{if($1=="ERR") print $3}')

	echo $evdp $err 0 $m $exp 0 0 >> tmp
done

range=$(cat tmp	| awk '{print $1,$2}' | minmax -C | awk '{print $1-1"/"$2*(1+0.1)"/"$3*(1-0.1)"/"$4*(1+0.1)}')

gmt begin $evid/depth png
gmt basemap -JX3i -R$range -Bxaf+l"Depth (km)" -Byaf+l"Misfit" -BWSne
cat tmp | gmt meca -Sm0.3 -Gblack
gmt end
	
