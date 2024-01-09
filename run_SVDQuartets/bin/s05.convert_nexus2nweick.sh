#/usr/bin/env bash

Usage="$0 <nexus_file>"

if [ $# -ne 1 ]; then
    echo "$Usage"
    exit 1
fi

source activate bioconvert

bioconvert nexus2newick "$1" "$1.nwk"

cat "$1.nwk" | sed 's/)[0-9]\+/)/g'> "$1.noSupportValueDsuit.nwk"

conda deactivate 
