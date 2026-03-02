import argparse
import glob
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('blast_db')
parser.add_argument('reads_dir')
parser.add_argument('build_dir')
parser.add_argument('--cpus', type=int, default=8, help='[%(default)i]')
arg = parser.parse_args()

params = ' '.join((
	'-r 1 -q -1',  # scoring system: +1, -1
	'-G 2 -E 1',   # gap system: -2 to open, -1 to exten
	'-e 1e-10',    # filter off poor matches
	'-m 8')        # use tabular output
)

for file in glob.glob(f'{arg.reads_dir}/*'):
	os.system(f'python3 fastq2fasta.py {file} > {arg.build_dir}/temp.fa')
	os.system(f'blastall -p blastn -d {arg.blast_db} -i {arg.build_dir}/temp.fa {params} -a {arg.cpus} > {arg.build_dir}/temp.blast')
	sys.exit('testing')

