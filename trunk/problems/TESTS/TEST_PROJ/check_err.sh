#!/bin/bash

rm -f file.txt

# select proper column from error dump file
cut -d';' -s -f6 "files/$1" > file.txt

# set tolerance
eps=1.00000E-11

# compare error to set tolerance
while read line; do
	flag=$(echo "${line} < ${eps}" | bc -l )
	if [[ "$flag" == 0 ]]
       	then
		echo "Fail!"
		echo $line
		rm -f file.txt
		break
	fi
done < file.txt

if [[ "$flag" != 0 ]]
then
	echo "Pass!"
fi

rm -f file.txt
