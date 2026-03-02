import argparse
import os
import sys

from grimoire.genome import Reader

parser = argparse.ArgumentParser()
parser.add_argument('fasta')
parser.add_argument('gff3')
arg = parser.parse_args()

seqs = {}
genome = Reader(fasta=arg.fasta, gff=arg.gff3)
for chrom in genome:
	for gene in chrom.ftable.build_genes():
		for tx in gene.transcripts():
			seq = tx.seq_str()
			if seq not in seqs: seqs[seq] = []
			seqs[seq].append(tx.id)

for seq, ids in seqs.items():
	print('>', '|'.join(ids), sep='')
	for i in range(0, len(seq), 80):
		print(seq[i:i+80])