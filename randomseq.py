import argparse
import random
import sys

import korflab

parser = argparse.ArgumentParser()
parser.add_argument('chromosomes', type=int)
parser.add_argument('length', type=int)
parser.add_argument('--seed', type=int, default=2)
arg = parser.parse_args()

random.seed(arg.seed)

for i in range(arg.chromosomes):
	seq = ''.join(random.choices('ACGT', k=arg.length))
	print(f'>chr_{i}')
	for j in range(0, len(seq), 80):
		print(seq[j:j+80])
	