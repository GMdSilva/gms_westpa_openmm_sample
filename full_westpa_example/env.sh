#!/bin/bash

source /gpfs/home/gmontei2/.bashrc
export WEST_ROOT=/oscar/home/gmontei2/data/gmontei2/.conda/envs/westpa-2.0

export WM_ZMQ_MASTER_HEARTBEAT=100
export WM_ZMQ_WORKER_HEARTBEAT=100
export WM_ZMQ_TIMEOUT_FACTOR=200
export WEST_SIM_ROOT="$PWD"
export SIM_NAME=$(basename $WEST_SIM_ROOT)
export BASH=/bin/bash
export PERL=/usr/bin/perl
export ZSH=/bin/zsh
export IFCONFIG=/bin/ifconfig
export CUT=/usr/bin/cut
export TR=/usr/bin/tr
export LN=/bin/ln
export CP=/bin/cp
export RM=/bin/rm
export SED=/bin/sed
export CAT=/bin/cat
export HEAD=/bin/head
export TAR=/bin/tar
export AWK=/usr/bin/awk
export PASTE=/usr/bin/paste
export GREP=/bin/grep
export SORT=/usr/bin/sort
export UNIQ=/usr/bin/uniq
export HEAD=/usr/bin/head
export MKDIR=/bin/mkdir
export ECHO=/bin/echo
export DATE=/bin/date
