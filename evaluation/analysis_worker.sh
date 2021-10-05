#!/usr/bin/env bash
help () {
    echo "usage: ./script <harp_executable> <einsim_executable> <output_directory> <json_output_dir> <int: num_cores> <string: analysis> <comma-separated-list: int: K (dataword length)> <int: CSEED-min (inclusive)> <int: CSEED stride> <int: CSEED-max (inclusive)> <int: WSEED-min (inclusive)> <int: WSEED stride> <int: WSEED-max (inclusive)> <int: num codes to simulate> <int: num words to simulate>";
}
[ "$#" -eq 15 ] || { echo "[ERROR] invalid arguments supplied: $@" && help && exit -1; }

HARP_EXEC=${1}
EINSIM_EXEC=${2}
OUTDIR=${3}
JSON_DIR=${4}
NCORES=${5}
ANALYSIS=${6}
KLIST=$(echo "${7}" | tr ',' ' ')
CSEED_MIN=${8}
CSEED_STRIDE=${9}
CSEED_MAX=${10}
WSEED_MIN=${11}
WSEED_STRIDE=${12}
WSEED_MAX=${13}
NCODES=${14}
NWORDS=${15}
# NCORES=$(nproc --all)
[ -f ${HARP_EXEC} ] || { echo "[ERROR]: invalid harp executable: \"${HARP_EXEC}\"" && exit -1; }
[ -f ${EINSIM_EXEC} ] || { echo "[ERROR]: invalid einsim executable: \"${EINSIM_EXEC}\"" && exit -1; }
[ ! -f ${JSON_DIR} ] || { echo "[ERROR]: invalid JSON output directory - detected file: \"${JSON_DIR}\"" && exit -1; }
[ "${NCORES}" -ge 0 ] || { echo "[ERROR]: invalid core count specified: \"${NCORES}\"" && exit -1; }
if [ -d ${OUTDIR} ]; then
    echo "[ERROR]: requested output directory already exists: \"${OUTDIR}\", containing "$(ls ${OUTDIR}/ | wc -l)" files/directories..."
    read -n 1 -s -r -p "Press [Enter] to overwrite or ^C to exit..."
    echo "OK"
    rm -rf ${OUTDIR}/*
fi
for K in ${KLIST}; do
    [ "${K}" -ge 0 ] || { echo "[ERROR]: invalid K value specified: \"${K}\"" && exit -1; }
done
[ "${CSEED_MIN}" -ge 0 ] || { echo "[ERROR]: invalid CSEED_min specified: \"${RMIN}\"" && exit -1; }
[ "${CSEED_STRIDE}" -ge 0 ] || { echo "[ERROR]: invalid CSEED_stride specified: \"${RSTRIDE}\"" && exit -1; }
[ "${CSEED_MAX}" -ge 0 ] || { echo "[ERROR]: invalid CSEED_max specified: \"${RMAX}\"" && exit -1; }
[ "${WSEED_MIN}" -ge 0 ] || { echo "[ERROR]: invalid WSEED_min specified: \"${RMIN}\"" && exit -1; }
[ "${WSEED_STRIDE}" -ge 0 ] || { echo "[ERROR]: invalid WSEED_stride specified: \"${RSTRIDE}\"" && exit -1; }
[ "${WSEED_MAX}" -ge 0 ] || { echo "[ERROR]: invalid WSEED_max specified: \"${RMAX}\"" && exit -1; }
[ "${NCODES}" -ge 0 ] || { echo "[ERROR]: invalid n_codes specified: \"${NCODES}\"" && exit -1; }
[ "${NWORDS}" -ge 0 ] || { echo "[ERROR]: invalid n_words specified: \"${NWORDS}\"" && exit -1; }

# prepare output directories
mkdir -p ${OUTDIR}
[ -d ${JSON_DIR} ] || mkdir -p ${JSON_DIR}

# kill all child processes on exit
reap_children () {
    if [ -n "$(jobs -p)" ]; then 
        kill $(jobs -p)
    fi
}
trap 'echo Reaping child processes and exiting... && reap_children' EXIT

# launch runs
echo "[INFO] launching runs for K in {${KLIST}} with CSEED in {seq ${CSEED_MIN} ${CSEED_STRIDE} ${CSEED_MAX}} and  WSEED in {seq ${WSEED_MIN} ${WSEED_STRIDE} ${WSEED_MAX}}"
RUNNING=0
for K in ${KLIST}
do
    CSEED_LIST=($(seq ${CSEED_MIN} ${CSEED_STRIDE} ${CSEED_MAX}))
    WSEED_LIST=($(seq ${WSEED_MIN} ${WSEED_STRIDE} ${WSEED_MAX}))
    if [ "${#CSEED_LIST[@]}" -ne "${#WSEED_LIST[@]}" ]; then
        echo "[ERROR] mismatch between CSEED and WSEED list"
        exit -1
    fi

    for seed_idx in "${!CSEED_LIST[@]}"
    do
        if (( RUNNING >= NCORES )); then
            wait -n
            echo "[INFO] job complete"
            ((RUNNING--))
        fi

        CSEED=${CSEED_LIST[seed_idx]}
        WSEED=${WSEED_LIST[seed_idx]}
        OUTFILE="${OUTDIR}/${ANALYSIS}_k${K}_r${CSEED}_s${WSEED}_c${NCODES}_w${NWORDS}.log"
        if [[ -f "${OUTFILE}" ]]; then
            continue
        fi
        
        CMD="./${HARP_EXEC} ${EINSIM_EXEC} ${ANALYSIS} -j ${JSON_DIR} -k ${K} -r ${CSEED} -s ${WSEED} -c ${NCODES} -w ${NWORDS}"
        echo "[INFO] launching job ${OUTFILE}"
        if ! command -v /usr/bin/time &> /dev/null
        then
            ${CMD} > ${OUTFILE} 2>&1 &
        else
            /usr/bin/time -v ${CMD} > ${OUTFILE} 2>&1 &
        fi
        ((RUNNING++))
    done
done

while (( RUNNING > 0 ))
do
    echo "[INFO] waiting for ${RUNNING} jobs"    
    wait -n
    ((RUNNING--))
done
