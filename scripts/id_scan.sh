#!/bin/bash
# Pre-commit scan: fails if a study identifier, GUID, family/sample ID or credential pattern is found in tracked files.
# Usage: bash scripts/id_scan.sh   (from the repository root)
set -u
PAT='(MONT|UCLA|PENN|TORO)[0-9]{3,4}|NDAR[A-Z0-9]{8}|SP0[0-9]{6}|SF0[0-9]{5}|REACH[0-9]{4,}|card_[0-9]+|AKIA[0-9A-Z]{16}|aws_secret_access_key|BEGIN (RSA|OPENSSH) PRIVATE KEY'
files=$(git ls-files 2>/dev/null || find . -type f -not -path './.git/*')
hits=$(echo "$files" | grep -v '^scripts/id_scan.sh$' | xargs -d '\n' grep -nE "$PAT" 2>/dev/null)
if [ -n "$hits" ]; then echo "ID/credential patterns found:"; echo "$hits"; exit 1; fi
data=$(echo "$files" | grep -E '\.(csv|tsv|parquet|xlsx|avro|edf|rds|pptx|pdf)$' | grep -v '^headband_sleep_architecture/tests/')
if [ -n "$data" ]; then echo "Data or report files tracked:"; echo "$data"; exit 1; fi
echo "id_scan: clean"
