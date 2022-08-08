#!/bin/bash

evid=$1
evdp=$2

gmt set FONT_TITLE 10
gmt begin syn_ob png,pdf
	nst=$(cat $evid/sdc1.txt | wc -l)
	echo $nst
gmt subplot begin ${nst}x6 -Fs2c/1.5c  -B
	for ((ist=1;ist<=$nst;ist++)) 
	do
		one=$(awk -v ist=$ist '{if(NR==ist)print $1,$2,$3}' $evid/sdc1.txt)
		stnm=$(echo $one | awk '{print $1}')
		dist=$(echo $one | awk '{printf("%6.2f", sqrt($2*$2+$3*$3))}')
		az=$(echo $one | awk '{PI=3.1415926;az=atan2($2,$3)*180/PI;if($2<0)az=360+az;printf("%6.2f", az)}')
		echo $stnm $dist $az

		twin=$(cat $evid/figd/$stnm.s? | gmt info -C -I1 | awk '{print $2}')

		figid=$(awk -v ist=$ist 'BEGIN{print (ist-1)*6+0}')
#		echo $figid
		gmt subplot set $figid
		echo 0 1.5 $stnm | gmt text -JX? -R-5/5/-6/6 -F+f8p -N
		echo 0 -1.5 $dist","$az | gmt text -F+f6p -N

		figid=$(awk -v ist=$ist 'BEGIN{print (ist-1)*6+1}')
		gmt subplot set $figid
		awk '{print $1,$2}' $evid/figd/$stnm.pr | gmt plot -JX? -R0/$twin/-1/1 -Wblue
		awk '{print $1,$3}' $evid/figd/$stnm.pr | gmt plot -Wred
	
		cor=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1f", $4)}')	
		tf=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%4.1f", $5)}')	
		echo -5 0.2 $cor | gmt text -F+f6p+jML -N
		echo -5 -0.2 $tf | gmt text -F+f6p+jML -N
		x=$(awk -v twin=$twin 'BEGIN{print twin*2/3}')
		obmax=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1e", $6)}')	
		echo $x 0.8 $obmax | gmt text -F+f6p,1,blue -N
		synmax=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1e", $7)}')	
		echo $x -0.8 $synmax | gmt text -F+f6p,1,red -N

		if [ $ist -eq 1 ]
		then
			gmt basemap -JX? -R0/$twin/-1/1 -B+t"Pr"
		fi

		figid=$(awk -v ist=$ist 'BEGIN{print (ist-1)*6+2}')
		gmt subplot set $figid
		awk '{print $1,$2}' $evid/figd/$stnm.pz | gmt plot -JX? -R0/$twin/-1/1 -Wblue
		awk '{print $1,$3}' $evid/figd/$stnm.pz | gmt plot -Wred

		cor=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1f", $8)}')	
		tf=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%4.1f", $9)}')	
		echo -5 0.2 $cor | gmt text -F+f6p+jML -N
		echo -5 -0.2 $tf | gmt text -F+f6p+jML -N
		x=$(awk -v twin=$twin 'BEGIN{print twin*2/3}')
		obmax=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1e", $10)}')	
		echo $x 0.8 $obmax | gmt text -F+f6p,1,blue -N
		synmax=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1e", $11)}')	
		echo $x -0.8 $synmax | gmt text -F+f6p,1,red -N
		figid=$(awk -v ist=$ist 'BEGIN{print (ist-1)*6+3}')

		if [ $ist -eq 1 ]
		then
			gmt basemap -JX? -R0/$twin/-1/1 -B+t"Pz"
		fi

		gmt subplot set $figid
		awk '{print $1,$2}' $evid/figd/$stnm.sr | gmt plot -JX? -R0/$twin/-1/1 -Wblue
		awk '{print $1,$3}' $evid/figd/$stnm.sr | gmt plot -Wred

		cor=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1f", $12)}')	
		tf=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%4.1f", $13)}')	
		echo 0 0.2 $cor | gmt text -F+f6p+jML -N
		echo 0 -0.2 $tf | gmt text -F+f6p+jML -N
		x=$(awk -v twin=$twin 'BEGIN{print twin*2/3}')
		obmax=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1e", $14)}')	
		echo $x 0.8 $obmax | gmt text -F+f6p,1,blue -N
		synmax=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1e", $15)}')	
		echo $x -0.8 $synmax | gmt text -F+f6p,1,red -N
		figid=$(awk -v ist=$ist 'BEGIN{print (ist-1)*6+4}') 

		if [ $ist -eq 1 ]
		then
			gmt basemap -JX? -R0/$twin/-1/1 -B+t"Sr"
		fi

		gmt subplot set $figid
		awk '{print $1,$2}' $evid/figd/$stnm.st | gmt plot -JX? -R0/$twin/-1/1 -Wblue
		awk '{print $1,$3}' $evid/figd/$stnm.st | gmt plot -Wred

		cor=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1f", $16)}')	
		tf=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%4.1f", $17)}')	
		echo 0 0.2 $cor | gmt text -F+f6p+jML -N
		echo 0 -0.2 $tf | gmt text -F+f6p+jML -N
		x=$(awk -v twin=$twin 'BEGIN{print twin*2/3}')
		obmax=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1e", $18)}')	
		echo $x 0.8 $obmax | gmt text -F+f6p,1,blue -N
		synmax=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1e", $19)}')	
		echo $x -0.8 $synmax | gmt text -F+f6p,1,red -N
		figid=$(awk -v ist=$ist 'BEGIN{print (ist-1)*6+5}') 

		if [ $ist -eq 1 ]
		then
			gmt basemap -JX? -R0/$twin/-1/1 -B+t"St"
		fi

		gmt subplot set $figid
		awk '{print $1,$2}' $evid/figd/$stnm.sz | gmt plot -JX? -R0/$twin/-1/1 -Wblue
		awk '{print $1,$3}' $evid/figd/$stnm.sz | gmt plot -Wred

		cor=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1f", $20)}')	
		tf=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%4.1f", $21)}')	
		echo 0 0.2 $cor | gmt text -F+f6p+jML -N
		echo 0 -0.2 $tf | gmt text -F+f6p+jML -N
		x=$(awk -v twin=$twin 'BEGIN{print twin*2/3}')
		obmax=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1e", $22)}')	
		echo $x 0.8 $obmax | gmt text -F+f6p,1,blue -N
		synmax=$(grep $stnm $evid/figd/plotd.txt | awk '{printf("%5.1e", $23)}')	
		echo $x -0.8 $synmax | gmt text -F+f6p,1,red -N

		if [ $ist -eq 1 ]
		then
			gmt basemap -JX? -R0/$twin/-1/1 -B+t"Sz"
		fi

#		echo $ist $nst
		if [ $ist -eq $nst ] 
		then
			dx=$(awk -v x=$twin 'BEGIN{print x/2}')
		 	gmt basemap -JX? -R0/$twin/-1/1 -BS -Bxa${dx}
		fi
	done
gmt subplot end

gmt subplot begin 1x1 -Fs8c/2c -B -Yh+1c
	gmt subplot set 0
	gmt basemap -JX? -R0/10/-1/1	
	gmt plot -W1p,blue << EOF
4 0.4
5 0.4
EOF
	gmt plot -W1p,red << EOF
4 -0.4
5 -0.4
EOF
	echo 6 0.4 "obs." | gmt text -F+f10p,blue
	echo 6 -0.4 "syn." | gmt text -F+f10p,red
	m=$(cat $evid/output_$evdp.txt | awk '{if($1=="M")print $3,$4,$5,$6,$7,$8}')
	exp=$(cat $evid/output_$evdp.txt | awk '{if($1=="exp") print $3+7}')
	echo 2 0 0 $m $exp | gmt meca -Sm1
gmt subplot end

gmt end
