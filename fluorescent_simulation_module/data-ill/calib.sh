#! bash

w=`expr 263 - 1`;
h=`expr 337 - 1`;
mw=60;
step=40;

DEST="clip";

mkdir $DEST;

for (( x=0; x <=$w; x=`expr $x + $step` ))
do
    for (( y=0; y <=$h; y=`expr $y + $step` ))
    do
        xx=$x;
        yy=$y;
        if [ `expr $xx + $mw` -gt $w ]
        then
            xx=`expr $w - $mw`;
        fi
        if [ `expr $yy + $mw` -gt $h ]
        then
            yy=`expr $h - $mw`;
        fi
        mkdir $DEST"/"$xx"-"$yy"-60-60"
        ./calib -i cellG_GAUSS -o $DEST"/"$xx"-"$yy"-60-60" -s -1 -r $xx,$yy,60,60
    done
done
