#!/bin/bash

path=/project/davidwcr_264/Projects/KGP/KGLIHFHLP/found.24_wTopoff.release.no_fq.txt

while read -r path; do
  if [ -L "$path" ] && [ ! -e "$path" ]; then
    echo "Dead symlink: $path"
  fi
done < deadfilelist.txt
