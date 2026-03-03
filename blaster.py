import argparse
import glob
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('mRNA_file', help='.fa.gz')
parser.add_argument('reads_dir', help='fastq files')
parser.add_argument('--build', default='build',
	help='build directory [%(default)s]')
parser.add_argument('--cpus', type=int, default=8, help='[%(default)i]')
arg = parser.parse_args()

if not os.path.exists(f'{arg.build}/db.fa.nsq'):
	os.system(f'gunzip -c {arg.mRNA_file} > build/db.fa')
	os.system(f'formatdb -p F -i build/db.fa')

params = ' '.join((
	'-r 1 -q -1',  # scoring system: +1, -1
	'-G 2 -E 1',   # gap system: -2 to open, -1 to exten
	'-e 1e-10',    # filter off poor matches
	'-m 8')        # use tabular output
)

for file in glob.glob(f'{arg.reads_dir}/*'):
	os.system(f'python3 fastq2fasta.py {file} > {arg.build}/temp.fa')
	os.system(f'blastall -p blastn -d build/db.fa -i {arg.build}/temp.fa {params} -a {arg.cpus} > {arg.build}/temp.blast')
	sys.exit('testing')
