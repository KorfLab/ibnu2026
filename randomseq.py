import argparse
import random
import sys

import korflab

parser = argparse.ArgumentParser()
parser.add_argument('fasta')
parser.add_argument('reads', type=int)
parser.add_argument('readlen', type=int)
parser.add_argument('--subrate', type=float, default=0.1,
	help='[%(default).3f]')
parser.add_argument('--indrate', type=float, default=0.1,
	help='[%(default).3f]')
arg = parser.parse_args()

for defline, seq in korflab.readfasta(arg.fasta):
	chrom = defline.split()[0]
	for i in range(arg.reads):
		pos = random.randint(0, len(seq) - arg.readlen)
		rseq = []
		rdef = f'>{chrom}:{pos}-{pos+arg.readlen}'
		for j in range(arg.readlen):
			if random.random() < arg.subrate:
				rseq.append(random.choice('acgt'))
			elif random.random() < arg.indrate:
				if random.random() < 0.5: # insert
					rseq.append(seq[pos+j])
					rseq.append(random.choice('ACGT'))
				# else is deletion, which is not to add anything
			else:
				rseq.append(seq[pos+j])
		rseq = ''.join(rseq)
		if random.random() < 0.5:
			rseq = korflab.anti(rseq)
			rdef += '-'
		else:
			rdef += '+'
		print(rdef, ''.join(rseq))