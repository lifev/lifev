#!/bin/bash

#PBS -l mppwidth=NCPUS
#PBS -l mppnppn=4
#PBS -l walltime=01:00:00

#PBS -A hpx1vp2a

#PBS -V
#PBS -e runTest.MESH.PXP1.NCPUS.err
#PBS -o runTest.MESH.PXP1.NCPUS.out

### Declares  whether  the  job  is rerunnable.  
#PBS -r n

SIZE=MESH

. /usr/local/packages/deisa/bash.bashrc.cray
module switch PrgEnv-pgi PrgEnv-gnu

module load hdf5-parallel
module load acml/4.2.0

# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

# Set the number of MPI tasks
export NPROC=`qstat -f $PBS_JOBID | awk '/mppwidth/ {print $3}'`

# Set the number of MPI tasks per node
export NTASK=`qstat -f $PBS_JOBID | awk '/mppnppn/  {print $3}'`


## configuring my test:


EXEC_DIR="Ethiersteinman/cube${SIZE}.PXP1/${NPROC}"
mkdir -p $EXEC_DIR
cd       $EXEC_DIR

# setting the data file and mesh:
if ! test -a ../data; then
    cp ../../../data ../
    sed "s/cube4x4.mesh/cube$SIZE.mesh/g" -i ../data
    sed "s/P1/PX/g" -i ../data
fi

test -a data || ln -s  ../data .
test -a Mesh || ln -s  ../Mesh .

grep mesh_dir data ../data ../../../data

date

echo running on NPROC=$NPROC and NTASK$NTASK

aprun -n $NPROC -N $NTASK $PBS_O_WORKDIR/test_ethiersteinman

date
