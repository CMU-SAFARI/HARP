#!/usr/bin/env bash
pushd .
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
Z3_DIR=z3
cd ${SCRIPT_DIR}/${Z3_DIR}
python scripts/mk_make.py --prefix=${SCRIPT_DIR}/${Z3_DIR}
cd build
make $@
make install
popd