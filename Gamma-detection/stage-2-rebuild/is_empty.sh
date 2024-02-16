#!/bin/bash

mkdir -p output

for ((i = 0; i <= 1999; i++)); do
    file="output/${i}_ans.txt"

    if [ -f "$file" ]; then
        if [ ! -s "$file" ]; then
            echo "File $file exists and is empty."
        fi
    else
        echo "File $file does not exist."
    fi
done
