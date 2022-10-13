#!/bin/bash
#echo $1
#head -1 $1
#awk -F, '{print $1}' $1
awk -F, '{print $1}' $1 | sort -k1,1 | uniq -c | sort -n

