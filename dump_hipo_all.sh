#!/usr/bin/env bash
# To run this script:
# 1. Open a terminal at /w/hallb-scshelf2102/clas12/nlbucuru/CNDHitFinder_Optimization
# 1. chmod +x dump_hipo_all.sh
# 2. bash
# 3. ./dump_hipo_all.sh or bash dump_hipo_all.sh

module use /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/modulefiles/ module load clas12
module load coatjava/13.5.2

# --- Output file ---
out="/w/hallb-scshelf2102/clas12/nlbucuru/CNDHitFinder_Optimization/hipo_dump_all.txt"
: > "$out"

# --- Base path ---
base="/work/clas12b/users/lixu/singleParticle"

# --- Test folders ---
folders=("testOSG" "testCJ0" "testCJ1")

# --- Loop ---
for name in "${folders[@]}"; do
  echo "========== $name ==========" >> "$out"
  for i in {0..19}; do
    f="$base/$name/output$i/dst-$i.hipo"
    if [[ -f "$f" ]]; then
      echo "==================== $f ====================" >> "$out"
      hipo-utils -info "$f" >> "$out"
      echo >> "$out"
    fi
  done
done

echo "Wrote: $out"

