#!/usr/bin/env bash

root="$(pwd)"

cd "$(dirname "$0")" || exit

gbk="$root/$1"
fna="$root/$2"

tag="${gbk/.gbk}"

./glimmerize.py "$gbk"
build-icm -r "$tag.icm" < "$tag.training.fna"
elph "$tag.upstream.fna" LEN=6 | ./get-motif-counts.awk > "$tag.motif"

exit
