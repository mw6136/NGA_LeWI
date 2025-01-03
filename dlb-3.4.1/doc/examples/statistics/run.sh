#!/usr/bin/env bash

# variables detected at configure time
MPIEXEC=

APP="./mpi_ompss_pils"
ARGS="--loads 100,250 --grain 1.0 --iterations 100000 --task-duration 50"

export NX_ARGS="--enable-dlb"
export LB_POLICY="No"
export LB_STATISTICS=1

$MPIEXEC -np 2 $APP $ARGS
