#!/bin/bash

# bash script to merge 2 Xdmf files coming from the hdf5 exporter.
# works on Linux platforms

if test $1& test $2& test $3
then
tac $1|sed "1,5 d"|tac > $3;
sed "1,8 d" $2 >> $3;
else
echo 'usage: ./mergeXdmf firstIn.xmf secondIn.xmf out.xmf'
fi
