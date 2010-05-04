#!/bin/bash

if test "x$3" == "x"; then
    echo usage : `basename $0` NCPUS MESH PX
    echo       example:
    echo         `basename $0` 16 32x32 P1Bubble
    exit
fi

NCPUS=$1
MESH=$2
PX=$3

sed "s/MESH/$MESH/g" runTest.Base.sh > runTest.$MESH.${PX}P1.$NCPUS.sh
sed "s/NCPUS/$NCPUS/g" -i runTest.$MESH.${PX}P1.$NCPUS.sh
sed "s/PX/$PX/g" -i runTest.$MESH.${PX}P1.$NCPUS.sh

chmod a+r runTest.$MESH.${PX}P1.$NCPUS.sh
chmod u+x runTest.$MESH.${PX}P1.$NCPUS.sh
ls -l runTest.$MESH.${PX}P1.$NCPUS.sh