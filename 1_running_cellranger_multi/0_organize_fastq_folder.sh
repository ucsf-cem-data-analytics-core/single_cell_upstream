#!/usr/bin/env bash

# Usage:
#   ./0_organize_fastq_folder.sh [-n] /path/to/fastq_folder
#   -n = dry run (no files moved)

set -euo pipefail

DRY_RUN=false

while getopts ":n" opt; do
  case "$opt" in
    n)
      DRY_RUN=true
      ;;
    *)
      echo "Usage: $0 [-n] /path/to/fastq_folder"
      exit 1
      ;;
  esac
done

shift $((OPTIND - 1))

FASTQ_DIR="${1:-.}"
cd "$FASTQ_DIR"

echo "Running in: $FASTQ_DIR"
if $DRY_RUN; then
  echo "DRY RUN MODE (no files will be moved)"
fi

# Only inspect files in the top-level directory, not inside sample folders.
while IFS= read -r -d '' file; do
  base_name="${file#./}"

  # Extract sample name (everything before _L00)
  sample="${base_name%%_L00*}"

  # Skip anything that does not match the expected pattern
  if [[ "$sample" == "$base_name" ]]; then
    echo "Skipping unexpected filename: $base_name" >&2
    continue
  fi

  if $DRY_RUN; then
    echo "[DRY RUN] mkdir -p \"$sample\""
    echo "[DRY RUN] mv \"$base_name\" \"$sample/\""
  else
    mkdir -p -- "$sample"
    mv -- "$base_name" "$sample/"
  fi
done < <(
  find . -maxdepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) -print0
)

echo "Done."