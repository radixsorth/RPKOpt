#!/bin/bash

# creates a summary from the reactor study

cd Direct-vs-adjoint

TARGET=summary.txt

echo -e "Summary:\n" > $TARGET
for f in *.log
do
    head -n 1 $f >> $TARGET
    tail -n 3 $f | head -n 2 >> $TARGET
    echo >> $TARGET
done
