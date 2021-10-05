#!/usr/bin/env bash
help () {
    echo "usage: ./script <harp_executable> <einsim_executable> <output_directory> <json_output_directory> <figures_output_directory> <analysis_to_perform> <int: dataword_length> <int: base random seed for ECC code generation> <int: base random seed for ECC word generation> <int: number of ECC codes to simulate> <int: number of ECC words per code to simulate> <int: number of tasks with this config to spawn> <int: num_cores (optional - nprocs by default)>";
}
[ "$#" -eq 13 ] || { echo "[ERROR] invalid number of arguments given" && help && exit -1; }

# https://stackoverflow.com/questions/59895/how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itsel?page=1&tab=votes#tab-top
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

HARP_EXEC=${1}
EINSIM_EXEC=${2}
OUT_DIR=${3}
JSON_DIR=${4}
FIG_DIR=${5}
ANALYSIS=${6}
K=${7}
CSEED=${8}
WSEED=${9}
NCODES=${10}
NWORDS=${11}
NTASKS=${12}
N_CPUS=${13}

###################################################
###### postcorrection probabilities analysis ######
###################################################
run_probabilities_analysis() {
    local ANALYSIS=probabilities
    local K=${1}
    local CSEED_MIN=${2}
    local WSEED_MIN=${3}
    local NCODES=${4}
    local NWORDS=${5}
    local NTASKS=${6}
    local CSTRIDE=$(( ${NCODES} ))
    local CMAX=$(( ${CSEED_MIN}+(${NTASKS}-1)*${CSTRIDE} ))
    local WSTRIDE=$(( ${NCODES} ))
    local WMAX=$(( ${WSEED_MIN}+(${NTASKS}-1)*${WSTRIDE} ))
    ${SCRIPT_DIR}/analysis_worker.sh ${HARP_EXEC} ${EINSIM_EXEC} ${OUT_DIR} ${JSON_DIR} ${N_CPUS} ${ANALYSIS} ${K} ${CSEED_MIN} ${CSTRIDE} ${CMAX} ${WSEED_MIN} ${WSTRIDE} ${WMAX} ${NCODES} ${NWORDS}
    if [ $? -eq 0 ]; then
        mkdir -p ${FIG_DIR}
        python3 ${SCRIPT_DIR}/../script/figure_4-parse_postcorrection_probabilities_data.py ${OUT_DIR} -o ${FIG_DIR}/
    fi
}

##########################################
###### profiler evaluation analysis ######
##########################################
run_profiler_evaluation_analysis() {
    local ANALYSIS=evaluations
    local K=${1}
    local CSEED_MIN=${2}
    local WSEED_MIN=${3}
    local NCODES=${4}
    local NWORDS=${5}
    local NTASKS=${6}
    local CSTRIDE=$(( ${NCODES} ))
    local CMAX=$(( ${CSEED_MIN}+(${NTASKS}-1)*${CSTRIDE} ))
    local WSTRIDE=$(( ${NCODES} ))
    local WMAX=$(( ${WSEED_MIN}+(${NTASKS}-1)*${WSTRIDE} ))
    ${SCRIPT_DIR}/analysis_worker.sh ${HARP_EXEC} ${EINSIM_EXEC} ${OUT_DIR} ${JSON_DIR} ${N_CPUS} ${ANALYSIS} ${K} ${CSEED_MIN} ${CSTRIDE} ${CMAX} ${WSEED_MIN} ${WSTRIDE} ${WMAX} ${NCODES} ${NWORDS}
    if [ $? -eq 0 ]; then
        mkdir -p ${FIG_DIR}
        python3 ${SCRIPT_DIR}/../script/figures_6to10-parse_evaluation_data.py ${JSON_DIR} ${OUT_DIR} -o ${FIG_DIR}/
    fi
}

if [ "${ANALYSIS}" = "evaluations" ]; then
    run_profiler_evaluation_analysis ${K} ${CSEED} ${WSEED} ${NCODES} ${NWORDS} ${NTASKS}
elif [ "${ANALYSIS}" = "probabilities" ]; then
    run_probabilities_analysis ${K} ${CSEED} ${WSEED} ${NCODES} ${NWORDS} ${NTASKS}
else
    echo "[ERROR] invalid analysis requested: ${ANALYSIS}"
    exit -1
fi