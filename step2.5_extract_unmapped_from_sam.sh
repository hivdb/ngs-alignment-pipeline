#! /bin/bash

set -e

SAMDIR=$1

OUTDIR=$1/unmapped

mkdir -p $OUTDIR

for file in `ls $SAMDIR/*.sam`; do
    name=`basename ${file%%.sam}`
    samtools view -S -f4 $file | awk '{print ">" $1 "/" $2 "\n" $10}' > $OUTDIR/$name.unmapped.fas
    echo $OUTDIR/$name.unmapped.fas
done
