#!/bin/bash

for j in {1..289}
do
	echo $j
	./scripts/randomForest_eachTaxon.py $j 100
done
