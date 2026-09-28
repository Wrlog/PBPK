#!/usr/bin/env bash
# Assemble each app in apps/ into a self-contained directory for the
# browser build: app.R and its R/ and data/ files, plus the shared R code
# and the model files. reference/ (mrgsolve) is deliberately left out.
#
# Usage: bash tools/assemble_apps.sh [output_dir]   (default: build)
set -euo pipefail
out="${1:-build}"
rm -rf "${out:?}"
mkdir -p "$out"
for app in apps/*/; do
  name=$(basename "$app")
  mkdir -p "$out/$name/R" "$out/$name/models"
  cp -r "$app". "$out/$name/"
  cp shared/*.R "$out/$name/R/"
  cp models/*.cpp "$out/$name/models/"
  echo "assembled $name"
done
