#!/bin/bash
for i in {1..100}
do
	echo $i
	./Main.exe -NSamples 1 -s 5020 -NQ 1 -QMin 1 -QMax 1 -y 2 -ID $i -doCent 1
done
