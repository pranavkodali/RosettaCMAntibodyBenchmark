#!/bin/bash
# run with `bash sc_parser.bash <file>`
file="$1"
new_file="${file%.*}.txt"
incr=0
while read -r line
do
  if ((incr >= 2)); then
    cut="${line:214}.pdb"
    echo "$cut" >> "$new_file"
  fi
  ((++incr))
done < "$file"
