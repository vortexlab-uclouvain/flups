#!/bin/bash
source $2

make -C $1 install -j4

FILE1=libflups_a2a.a
if test -f "$1/lib/$FILE1"; then
    echo "$FILE1 exists."
    exit 0
else
    echo "$FILE1 does not exist."
    exit 1
fi