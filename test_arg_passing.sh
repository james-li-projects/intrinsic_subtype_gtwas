#!/bin/bash

# Check if the required argument (output subdirectory) is provided
if [ "$#" -ne 1 ]; then
    echo "Approach for running GWAS is missing. Please include approach as a command line argument."
fi

# Define directories
PASSED_ARG="$1" 
echo $PASSED_ARG
