#!/bin/bash
for i in {1001..1020}
do
  reroot -l -b -q skimTree.C'('$i')'
done
