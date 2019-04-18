#!/usr/bin/env python3

import signal
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from pathlib import Path

from Bio import SeqIO

ext_to_fmt = {
	"fna": "fasta",
	"gbk": "gb"
}


def cslice(s, i, j):
	if i < 0:
		return s[max(len(s) - abs(i), j):] + s[:j]
	elif j > len(s):
		return s[i:] + s[:min(i, j - len(s))]
	else:
		return s[i:j]


def parse_cds(path):
	for record in SeqIO.parse(str(path), "genbank"):
		for feature in record.features:
			if feature.type == "CDS":
				yield record, feature, feature.location


def signal_handler(sig, frame):
	print("You pressed Ctrl+C!")
	sys.exit(0)


def parse_argv(argv):
	parser = ArgumentParser(
		description="glimmerizer",
		formatter_class=ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		"path",
		type=Path
	)
	parser.add_argument(
		"-nups", "--nups", type=int, default=25
	)
	parser.add_argument(
		"-trim", "--trim", type=int, default=3
	)
	parser.add_argument(
		"-key-id", "--key-id", default="protein_id"
	)
	parser.add_argument(
		"-key-description", "--key-description", default="product"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	key_id = args.key_id
	key_desc = args.key_description

	codons = OrderedDict((("ATG", 0), ("GTG", 0), ("TTG", 0)))

	training = []
	upstream = []

	for idx, ele in enumerate(parse_cds(args.path), start=1):
		record, feature, location = ele
		trecord = feature.extract(record)[:-args.trim]
		codon = str(trecord.seq[:3])
		if codon in codons:

			codons[codon] += 1

			start, end, strand = location.start, location.end, " +-"[location.strand]

			fid = feature.qualifiers.get(key_id, [f"cds-{idx}"])[0]
			trecord.id = f"{fid}|{start + 1}-{end - args.trim}|{strand}"
			trecord.description = feature.qualifiers.get(key_desc, [""])[0]

			urecord = cslice(record, start - args.nups, start)
			urecord.id = f"{fid}|{start + 1 - args.nups}-{start}|{strand}"
			urecord.description = trecord.description

			training.append(trecord)
			upstream.append(urecord)

	SeqIO.write(training, args.path.with_suffix(".training.fna"), "fasta")
	SeqIO.write(upstream, args.path.with_suffix(".upstream.fna"), "fasta")
	with args.path.with_suffix(".startuse.csv").open("w") as file:
		total = sum(codons.values())
		print(*(f"{val / total:0.3f}" for val in codons.values()), sep=",", file=file)

	return 0


if __name__ == "__main__":
	signal.signal(signal.SIGINT, signal_handler)
	sys.exit(main(sys.argv))
