#!/bin/bash
#SBATCH --job-name=westpa2
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# #SBATCH --constraint=geforce3090
#SBATCH -p gpu --gres=gpu:8
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
set -x
cd $SLURM_SUBMIT_DIR
source env.sh || exit 1

env | sort

cd $WEST_SIM_ROOT


SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info-$SLURM_JOBID.json
# start server
$WEST_ROOT/bin/w_run --work-manager=zmq --n-workers=0 --zmq-mode=master --zmq-write-host-info=$SERVER_INFO --zmq-comm-mode=tcp &> west-$SLURM_JOBID.log &

for ((n=0; n<60; n++)); do
    if [ -e $SERVER_INFO ] ; then
        echo "== server info file $SERVER_INFO =="
        cat $SERVER_INFO
        break
    fi
    sleep 1
done

# exit if host info file doesn't appear in one minute
if ! [ -e $SERVER_INFO ] ; then
    echo 'server failed to start'
    exit 1
fi

for node in $(scontrol show hostname $SLURM_NODELIST); do
   ssh -o StrictHostKeyChecking=no $node $PWD/node.sh $SLURM_SUBMIT_DIR $SLURM_JOBID $node --verbose \
                     --work-manager=zmq --n-workers=8 --zmq-mode=node \
                     --zmq-read-host-info=$SERVER_INFO \
                     &> west-$JOB_ID-$machine.log &
done

wait
