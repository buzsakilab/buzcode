#!/bin/bash
# syntax : ./makeDataset xx Path
# where :
# xx is the rat number
# Path is the parent path of the Rat Data (without '/' at the end)
# it will look in directories Path/Ratxx/xx*
# it creates datasetRatx in the current directory

datasetFile=$1'/datasetRat'$2
echo $datasetFile

if [ -r $datasetFile ]
	then
	rm $datasetFile
fi

path=$1/Rat$2/
echo $path
for i in $path$2*
	
	do	
	echo ${i/$path} >> $datasetFile
	
done