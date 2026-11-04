#!/usr/bin/env bash
# Download the reference data for the KlebTrace pipeline into databases/ (about 0.9 GB):
# the AMRFinderPlus and PlasmidFinder Apptainer images, the AMRFinderPlus database
# 2026-05-15.1 and the PlasmidFinder database at commit 81c11f4, the versions used in Galaxy.
# Bakta's database (39.4 GB) is not downloaded; the pipeline reuses Bakta's Galaxy output.
#
# Usage: bash scripts/download_databases.sh   (from the repository root)
set -euo pipefail
mkdir -p databases/containers databases/amrfinderplus
cd databases

export APPTAINER_CACHEDIR=$PWD/containers/cache
for image in ncbi-amrfinderplus:4.2.7--hf69ffd2_0 plasmidfinder:2.1.6--py314hdfd78af_2; do
  sif=containers/quay.io-biocontainers-${image/:/-}.img
  [ -s "$sif" ] || apptainer pull --name "$sif" "docker://quay.io/biocontainers/$image"
done
rm -rf containers/cache

# AMRFinderPlus database: fetch each file in the version folder, then index it
v=2026-05-15.1
url=https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.2/$v/
mkdir -p amrfinderplus/$v
for f in $(curl -s "$url" | grep -o 'href="[^"/?][^"]*"' | sed 's/href="//; s/"$//' | grep -v '^http'); do
  curl -s -f -o "amrfinderplus/$v/$f" "$url$f"
done
apptainer exec containers/quay.io-biocontainers-ncbi-amrfinderplus-4.2.7--hf69ffd2_0.img \
  amrfinder_index "amrfinderplus/$v"

# PlasmidFinder database at the commit used in Galaxy
[ -d plasmidfinder/.git ] || git clone -q https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git plasmidfinder
git -C plasmidfinder checkout -q 81c11f4

du -sh containers amrfinderplus plasmidfinder
