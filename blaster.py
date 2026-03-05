import argparse
import glob
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('mRNA_file', help='*.fa.gz')
parser.add_argument('reads_dir', help='*fastq.gz')
parser.add_argument('--build', default='build',
	help='build directory [%(default)s]')
parser.add_argument('--cpus', type=int, default=4, help='[%(default)i]')
parser.add_argument('--testing', action='store_true')
arg = parser.parse_args()

if not os.path.exists(f'{arg.build}/db.fa.nsq'):
	os.system(f'mkdir -p {arg.build}/blast')
	os.system(f'gunzip -c {arg.mRNA_file} > {arg.build}/db.fa')
	os.system(f'formatdb -p F -i {arg.build}/db.fa')

params = ' '.join((
	'-p blastn',
	f'-d {arg.build}/db.fa',    # database file is temporary
	f'-i {arg.build}/temp.fa',  # query file is temporary
	f'-a {arg.cpus}',           # depends on computer
	'-r 1 -q -1',               # scoring system: +1, -1
	'-G 2 -E 1',                # gap system: -2 to open, -1 to extend
	'-e 1e-10',                 # filter off poor matches
	'-m 8')                     # use tabular output
)

for path in glob.glob(f'{arg.reads_dir}/*'):
	filename = path.split('/')[-1].split('.')[0]
	output = f'{arg.build}/blast/{filename}.tsv'
	if os.path.exists(output):
		print('done already, skipping', filename)
		continue
	print('processing', filename)
	os.system(f'python3 fastq2fasta.py {path} > {arg.build}/temp.fa')
	os.system(f'blastall {params} > {output}')
	if arg.testing: sys.exit('stopped after one alignment... testing')
