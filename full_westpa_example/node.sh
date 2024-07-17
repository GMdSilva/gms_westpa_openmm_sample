#!/bin/bash -l

set -x

cd $1; shift
source env.sh
export WEST_JOBID=$1; shift
export SLURM_NODENAME=$1; shift
echo "starting WEST client processes on: "; hostname
echo "current directory is $PWD"
echo "environment is: "
env | sort

$WEST_ROOT/bin/w_run "$@" &> west-$SLURM_NODENAME-node.log

echo "Shutting down."
