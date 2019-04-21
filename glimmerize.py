#!/usr/bin/env python3

import os
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import OrderedDict
from pathlib import Path
from signal import signal, SIGPIPE, SIG_DFL

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
	for record in SeqIO.parse(path, "genbank"):
		for feature in record.features:
			if feature.type == "CDS":
				yield record, feature, feature.location


def parse_argv(argv):
	parser = ArgumentParser(
		description="glimmerizer",
		formatter_class=ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		"file",
		type=FileType()
	)
	parser.add_argument(
		"-prefix", "--prefix"
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

	args.prefix = args.prefix or args.file.name

	return args


def main(argv):
	args = parse_argv(argv[1:])

	path = Path(args.prefix)
	os.makedirs(path.parent, exist_ok=True)

	key_id = args.key_id
	key_desc = args.key_description

	codons = OrderedDict((("ATG", 0), ("GTG", 0), ("TTG", 0)))

	training = []
	upstream = []

	with args.file as file:
		for idx, ele in enumerate(parse_cds(file), start=1):
			record, feature, location = ele

			try:
				trecord = feature.extract(record)[:-args.trim]
				codon = str(trecord.seq[:3])
			except ValueError:
				codon = ""

			if codon in codons:

				codons[codon] += 1

				start, end, strand = location.start, location.end, " +-"[location.strand]

				fid = feature.qualifiers.get(key_id, [f"cds-{idx}"])[0]
				trecord.id = f"{fid}|{start + 1}-{end - args.trim}|{strand}"
				trecord.description = feature.qualifiers.get(key_desc, [""])[0]

				training.append(trecord)

				uidx = start - args.nups
				topo = record.annotations.get("topology", "linear")
				urecord = record[uidx:start] if topo == "linear" else cslice(record, uidx, start)

				if len(urecord) > 0:
					urecord.id = f"{fid}|{start + 1 - args.nups}-{start}|{strand}"
					urecord.description = trecord.description
					upstream.append(urecord)

	SeqIO.write(training, path.with_suffix(".training.fna"), "fasta")
	SeqIO.write(upstream, path.with_suffix(".upstream.fna"), "fasta")
	with path.with_suffix(".startuse.csv").open("w") as file:
		total = sum(codons.values())
		print(*(f"{val / total:0.3f}" for val in codons.values()), sep=",", file=file)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
