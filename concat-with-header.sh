#!/bin/bash

if [ "$1" == "-h" -o $# -eq 0 ]; then
		echo "Usage: `basename $0` *files*"
		exit 0
fi

FLS=( $@ )

awk 'FNR==1 && NR!=1{next;}{print}' ${FLS[@]}
