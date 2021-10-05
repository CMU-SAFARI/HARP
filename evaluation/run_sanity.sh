#!/usr/bin/env bash
help () {
    echo "usage: ./script <harp_executable> <einsim_executable> <output_directory> <int: num_cores (optional - nprocs by default)>";
}
[ "$#" -eq 3 ] || [ "$#" -eq 4 ] || { help && exit -1; }

# https://stackoverflow.com/questions/59895/how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itsel?page=1&tab=votes#tab-top
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

HARP_EXEC=${1}
EINSIM_EXEC=${2}
OUT_DIR=${3}
NCPUS=${4:-$(nproc --all)}

# make necessary directories
if [ -d ${OUT_DIR} ]; then
    echo "[ERROR]: requested output directory already exists: \"${OUT_DIR}\", containing "$(ls ${OUT_DIR}/ | wc -l)" files/directories..."
    read -n 1 -s -r -p "Press [Enter] to overwrite or ^C to exit..."
    echo "OK"
    rm -rf ${OUT_DIR}/*
fi
JSON_DIR=${OUT_DIR}/ecc_code_database
EVAL_DIR=${OUT_DIR}/profiler_evaluation
PROB_DIR=${OUT_DIR}/postcorrection_probabilities
FIG_DIR=${OUT_DIR}/figures
mkdir -p ${OUT_DIR} ${JSON_DIR} ${FIG_DIR}

########################################
###### SANITY-CHECK CONFIGURATION ######
########################################
#
# {NCODES=1, NWORDS=1, NTASKS=16} Takes approximately 30 seconds using i7-7700HQ @ 2.80GHz
# 208 KB of output (12 KB/file)
# NOTE: runtime is exponentially related to K
#
if true; then
    ANALYSIS=evaluations
    K=16
    CSEED=0
    WSEED=0
    NCODES=1
    NWORDS=1
    NTASKS=4
    ${SCRIPT_DIR}/analysis_wrapper.sh ${HARP_EXEC} ${EINSIM_EXEC} ${EVAL_DIR} ${JSON_DIR} ${FIG_DIR} ${ANALYSIS} ${K} ${CSEED} ${WSEED} ${NCODES} ${NWORDS} ${NTASKS} ${NCPUS}
fi

#
# {NCODES=10, NWORDS=100, NTASKS=16} Takes approximately 30 seconds using i7-7700HQ @ 2.80GHz
# 76 MB of output (4.8 MB/file)
# NOTE: output file size scales linearly with both NCODES and NWORDS
#
if false; then
    ANALYSIS=probabilities
    K=64
    CSEED=0
    WSEED=0
    NCODES=10
    NWORDS=10
    NTASKS=16
    ${SCRIPT_DIR}/analysis_wrapper.sh ${HARP_EXEC} ${EINSIM_EXEC} ${PROB_DIR} ${JSON_DIR} ${FIG_DIR} ${ANALYSIS} ${K} ${CSEED} ${WSEED} ${NCODES} ${NWORDS} ${NTASKS} ${NCPUS}
fi

########################################
######## EVALUTION CONFIGURATION #######
########################################
#
# Post-correction error probabilites
# Yields 74 GB total output data (4.6 GB/file)
# ~90 minutes/task (1 CPU-day total) @ AMD EPYC 7742 and Intel(R) Xeon(R) Gold 5118 
#
# K=64
# CSEED=0
# WSEED=0
# NCODES=100
# NWORDS=10000
# NTASKS=16
# run_probabilities_analysis ${K} ${CSEED} ${WSEED} ${NCODES} ${NWORDS} ${NTASKS}

#
# Profiler evaluation analysis
# C++ Monte-Carlo simulation (i.e., data-generation)
#   Yields 3.4 GB total output data (~13 MB/file)
#   ~20 days/task (14 CPU-years total) @ AMD EPYC 7742 and Intel(R) Xeon(R) Gold 5118 
# Python analysis
#   Elapsed time: 3:09:20 (h:mm:ss)
#   Maximum resident memory: 38.13 GB
#
# K=${K}
# CSEED=0
# WSEED=0
# NCODES=10
# NWORDS=100
# NTASKS=256
# run_profiler_evaluation_analysis ${K} ${CSEED} ${WSEED} ${NCODES} ${NWORDS} ${NTASKS}