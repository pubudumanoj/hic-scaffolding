#!/bin/bash -ue
echo 
mkdir temp
mkdir merged
inputs=`printf I= | sed 's\ \ I=\'`
java -Xmx20G -Djava.io.tmpdir=temp/ -jar $PICARD_HOME/picard.jar MergeSamFiles $inputs OUTPUT=merged/Afraterculus.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT
