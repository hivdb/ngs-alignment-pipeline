#! /bin/bash

set -e

UNMAPPED_DIR=$1/unmapped

for file in `ls $UNMAPPED_DIR/*.unmapped.fas`; do
    out="${file%%.sam}.blastresult.txt"
    entrypoints/blastfilter $file refs/hxb2.fasta 50 $out 
    echo $out
done
