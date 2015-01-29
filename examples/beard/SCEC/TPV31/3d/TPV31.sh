#!/bin/bash
#PBS -j oe
#PBS -o TPV31.out
#PBS -N TPV31
# PBS -l hostlist=compute-2-1+compute-2-11+compute-2-13+compute-2-15+compute-2-17+compute-2-19+compute-2-3+compute-2-5+compute-2-7+compute-2-9
# PBS -l hostlist=compute-3-1+compute-3-11+compute-3-13+compute-3-15+compute-3-17+compute-3-19+compute-3-21+compute-3-23+compute-3-25+compute-3-27+compute-3-29+compute-3-3+compute-3-31+compute-3-5+compute-3-7+compute-3-9
# PBS -l hostlist=compute-7-1+compute-7-11+compute-7-13+compute-7-15+compute-7-17+compute-7-19+compute-7-21+compute-7-23+compute-7-25+compute-7-27+compute-7-29+compute-7-3+compute-7-31+compute-7-33+compute-7-35+compute-7-37+compute-7-39+compute-7-5+compute-7-7+compute-7-9
# PBS -l hostlist=compute-7-1+compute-7-11+compute-7-13+compute-7-15+compute-7-17+compute-7-19+compute-7-21+compute-7-23+compute-7-25+compute-7-27+compute-7-29+compute-7-3+compute-7-31+compute-7-33+compute-7-35+compute-7-37+compute-7-39+compute-7-7+compute-7-9
#PBS -l procs=512
#PBS -l pmem=2gb
#PBS -l naccesspolicy=singlejob
#PBS -m abe
#PBS -M jekozdon@nps.edu
# PBS -q beards
#PBS -l walltime=24:00:00

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

NCPU=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`

echo ------------------------------------------------------
echo ' This job is allocated on '${NCPU}' cpu(s)'
echo 'Job is running on node(s): '
cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: number of nodes is $NNODES
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

source /etc/profile
module load compile/intel/13.0.0 mpi/openmpi
root_dir=/smallwork/jekozdon/data/SCEC_CVWS/TPV31/3d

echo mkdir -p $root_dir/data
mkdir -p $root_dir/data

echo cd $root_dir
cd $root_dir

mpirun ~/bfam_build/compute_3/beard_3d /home/jekozdon/codes/bfam/examples/beard/SCEC/TPV31/3d/TPV31_base.lua
