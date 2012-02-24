#!/bin/bash

./test_partition -f data_partition
mpirun -np 4 ./test_cylinder -f data_cylinder

