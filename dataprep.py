import argparse
import glob
import os
import sys
from grimoire.genome import Reader


parser = argparse.ArgumentParser()
parser.add_argument('fasta')
parser.add_argument('gff3')
parser.add_argument('name')
parser.add_argument('--pct', type=float, default=0.90,
	help='pct identity cutoff for uniqueness [%(default).2f]')
parser.add_argument('--build', default='build',
	help='build directory [%(default)s]')
arg = parser.parse_args()


if os.getenv('CONDA_DEFAULT_ENV') != 'setup':
	sys.exit('activate setup environment first')
os.system(f'mkdir -p {arg.build}')

# part 1: run haman
if not os.path.exists(f'{arg.build}/{arg.name}'):
	print('running haman', file=sys.stderr)
	os.system(f'haman {arg.fasta} {arg.gff3} pcg {arg.build}/{arg.name} --plus --issuesok')

# part 2: get canonical/longest tx from each gene
txfiles = []
for fa in glob.glob(f'{arg.build}/{arg.name}/*.fa'):
	gff = fa[:-2] + 'gff3'
	genome = Reader(fasta=fa, gff=gff)
	chrom = next(genome)
	gene = chrom.ftable.build_genes()[0]
	
	if not gene.is_coding(): continue
	txs = gene.transcripts()
	longest_tx = txs[0]
	if len(txs) > 1:
		max_idx = 0
		max_len = txs[0].end - txs[0].beg
		for i, tx in enumerate(txs):
			if tx.end - tx.beg > max_len:
				max_len = tx.end - tx.beg
				max_idx = i
		longest_tx = txs[max_idx]
	
	txf = fa[:-2] + 'tx'
	txfiles.append(txf)
	with open(txf, 'w') as fp:
		print('>', gene.id, sep='', file=fp)
		print(longest_tx.tx_str(), file=fp)

# part 3: blast all-vs-all
db = f'{arg.build}/{arg.name}.txs'
out = f'{arg.build}/{arg.name}.blastn'
if not os.path.exists(out):
	print('running blast', file=sys.stderr)
	os.system(f'cat {arg.build}/{arg.name}/*.tx > {db}')
	os.system(f'formatdb -p F -i {db}')
	p = '-e 1e-10 -F F -m 8'
	os.system(f'blastall -p blastn -d {db} -i {db} {p} > {out}')

# part 4: remove paralogs
keep = set()
kill = set()
with open(out) as fp:
	for line in fp:
		f = line.split()
		keep.add(f[0])
		if f[0] == f[1]: continue # self
		if float(f[2]) > arg.pct:
			kill.add(f[1])
			kill.add(f[2])
print('removed ', len(kill), ' paralogs', file=sys.stderr)

# part 5: create final nuggets
mini_fa = f'{arg.name}-genome.fa'
mini_tx = f'{arg.name}-mRNA.fa'
mini_gf = f'{arg.name}-genome.gff3'
os.system('rm -f {mini_fa} {mini_tx} {mini_gf}')

with open(mini_fa, 'w') as mf, open(mini_tx, 'w') as mt, open(mini_gf, 'w') as gf:
	for fa in glob.glob(f'{arg.build}/{arg.name}/*.fa'):
		gff = fa[:-2] + 'gff3'
		txf = fa[:-2] + 'tx'
		os.system(f'cat {fa} >> {mini_fa}')
		os.system(f'cat {txf} >> {mini_tx}')
		os.system(f'cat {gff} >> {mini_gf}')