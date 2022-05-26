#!/bin/bash -ue
mkdir temp
mkdir merged
inputs=`echo  | sed 's\ \ -I \g'`
echo $inputs
java -Xmx40g -Djava.io.tmpdir=temp/ -jar $PICARD_HOME/picard.jar MergeSamFiles -I $inputs -OUTPUT merged/Afraterculus.bam -USE_THREADING TRUE -ASSUME_SORTED TRUE -VALIDATION_STRINGENCY LENIENT
