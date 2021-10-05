#!/usr/bin/env bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
EINSIM_DIR=einsim
cd ${SCRIPT_DIR}/${EINSIM_DIR}
make $@
cd -