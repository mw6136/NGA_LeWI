#!/usr/bin/env bash

# variables detected at configure time
DLB_HOME=/usr/local
MPIEXEC=

# variables to be modified by the user
TRACE=0
DLB=0
LEWI_KEEP_ONE_CPU=0
MPIEXEC_BIND_FLAG=""

APP="./mpi_ompss_pils"
ARGS="--loads 1000,2500 --grain 0.1 --iterations 50 --task-duration 500 --verbose"

if [[ $DLB == 1 ]] ; then
    if [[ $TRACE == 1 ]] ; then
        TRACENAME="pils_dlb.prv"
        PRELOAD="$EXTRAE_HOME/lib/libnanosmpitrace.so:$DLB_HOME/lib/libdlb_mpi_instr.so"
        export EXTRAE_CONFIG_FILE="extrae.xml"
        export NX_ARGS+=" --instrumentation=extrae"
    else
        PRELOAD="$DLB_HOME/lib/libdlb_mpi.so"
    fi
    export NX_ARGS+=" --enable-block --enable-dlb --force-tie-master"
    export DLB_ARGS+=" --lewi"

    # Advanced: use --lewi-keep-one-cpu
    if [[ $LEWI_KEEP_ONE_CPU == 1 ]] ; then
        if [[ $TRACE == 1 ]] ; then
            TRACENAME="pils_dlb_lewi_keep_one_cpu.prv"
        fi
        export DLB_ARGS+=" --lewi-keep-one-cpu"
    fi
else
    if [[ $TRACE == 1 ]] ; then
        TRACENAME="pils.prv"
        PRELOAD="$EXTRAE_HOME/lib/libnanosmpitrace.so"
        export EXTRAE_CONFIG_FILE="extrae.xml"
        export NX_ARGS+=" --instrumentation=extrae"
    fi
    export DLB_ARGS+=" --no-lewi"
fi

$MPIEXEC -np 2 $MPIEXEC_BIND_FLAG env LD_PRELOAD="$PRELOAD" $APP $ARGS

if [[ $TRACE == 1 ]] ; then
    $EXTRAE_HOME/bin/mpi2prv -f TRACE.mpits -no-keep-mpits -o "$TRACENAME"
    rm -f TRACE.spawn
fi
