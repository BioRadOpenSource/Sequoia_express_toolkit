#!/bin/bash

# enable debugging
[[ "$NXF_DEBUG_ENTRY" ]] && set -x

# wrap cli args with single quote to avoid wildcard expansion
cli=''; for x in "$@"; do cli+="'$x' "; done
cd /opt/biorad/
bash -c "nextflow -log /work/nextflow.log run main.nf $cli"
