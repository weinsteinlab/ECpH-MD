#!/bin/bash -l

#### Generate patch statements needed for psfgen use
#### USE: ./generate_residue_patches.sh myPDB.pdb segmentName

#### NOTE: script does NOT generate patch statements for disulfide bridges

rm ./patches.txt
grep "CA  ASP" $1 | tr -s ' ' ',' | cut -d ',' -f 6 | xargs -n1 echo "patch ASPP ${2}:" | sed 's/\(.*\) /\1/' >> ./patches.txt
grep "CA  GLU" $1 | tr -s ' ' ',' | cut -d ',' -f 6 | xargs -n1 echo "patch GLUP ${2}:" | sed 's/\(.*\) /\1/' >> ./patches.txt
