#!/usr/bin/env bash

glimmer3 \
	-l -o 50 -g 110 -t 30 -b "$tag1.motif" -P "$(cat "$tag1.startuse.csv")" \
	"$fna" "$tag1.icm" "$tag2"


#./predict_to_gff.py "$tag2.predict"
#
#
#bedtools getfasta -fi "$fna" -fo "$tag2.fna" -bed "$tag2.gff"
