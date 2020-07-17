#!/bin/bash

# input fasta file from submission
ASSEMBLY=$1

# make folder
mkdir chunks

# remove any chunks and list from previous iteration
rm chr.list
rm chunks/*

# genome headers
grep ">" $ASSEMBLY | sed 's_>__' | shuf | tee -a chr.list


# slit into 80 chunks
split -d -n l/80 chr.list chunks/genomechunk.

cd chunks/
rename genomechunk.0 genomechunk. genomechunk.0*
