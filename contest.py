import argparse
import os
import platform
import random
import re
import sys

import korflab

def run(cli):
	pid = os.getpid()
	tempfile = f'{DIR}/{pid}.temp'
	print(cli, file=sys.stderr)
	if os.system(f'{TIMER} {cli} 2> {tempfile}') != 0: sys.exit(f'FAILED: {cli}')
	runstats = {}
	with open(tempfile) as fp:
		for line in fp:
			if re.match(r'\s+User', line):
				runstats['utime'] = float(line.split()[-1])
			elif re.match(r'\s+System', line):
				runstats['stime'] = float(line.split()[-1])
			elif re.match(r'\s+Maximum', line):
				runstats['rsize'] = int(line.split()[-1])
	os.remove(tempfile)
	return runstats

def simulate_reads(infile, outfile, rlen, xcov, rate):
	rate /= 100
	with open(outfile, 'w') as fp:
		for defline, seq in korflab.readfasta(infile):
			chrom = defline.split()[0]
			rnum = int(len(seq) * xcov / rlen)
			for i in range(rnum):
				pos = random.randint(0, len(seq) - rlen)
				if 'NN' in seq[pos:pos+rlen]: continue
				rseq = []
				rdef = f'>{chrom}:{pos+1}-{pos+rlen}'
				for j in range(rlen):
					if random.random() < rate:
						rseq.append(random.choice('acgt'))
					elif random.random() < rate:
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
				print(rdef, ''.join(rseq), sep='\n', file=fp)

def test_bbmap(gfile, rfile, cpus):
	run(f'/usr/bin/time -v conda run -n bbmap bbmap.sh in={rfile} ref={gfile} nodisk=t threads={cpus} out={rfile}.bbmap > {rfile}.stdout 2> {rfile}.stderr')


def test_blast(gfile, rfile, cpus):
	r1 = run(f'conda run -n blast-legacy formatdb -p F -i {gfile}')
	r2 = run(f'conda run -n blast-legacy blastall -p blastn -d {gfile} -i {rfile} -a {cpus} -r 1 -q -1 -G 2 -E 1 -e 1e-10 -m 8 -o {rfile}.blastn')
	seen = set()
	aligned = 0
	total = 0
	with open(f'{rfile}.blastn') as fp:
		for line in fp:
			f = line.split()
			m = re.search(r'^(\w+):(\d+)\-(\d+)(.)$', f[0])
			if m:
				c1 = m.group(1)
				b1 = int(m.group(2))
				e1 = int(m.group(3))
				s1 = m.group(4)
			else:
				sys.exit('blast parsing error')
			if (c1, b1, e1, s1) in seen: continue
			c2 = f[1]
			b2 = int(f[8])
			e2 = int(f[9])
			if s1 == '-': b2, e2 = e2, b2
			tlen = e1 - b1 + 1
			alen = e2 - b2 + 1
			seen.add( (c1, b1, e1, s1) )
			aligned += alen
			total += tlen
	cpu = r1['utime'] + r1['stime'] + r2['utime'] + r2['stime']
	mem = max(r1['rsize'], r2['rsize'])
	return aligned / total, cpu, mem


def test_bowtie2(gfile, rfile, cpus):
	pass
def test_bwa(gfile, rfile, cpus):
	pass
def test_gmap(gfile, rfile, cpus):
	pass
def test_hisat2(gfile, rfile, cpus):
	pass
def test_minimap2(gfile, rfile, cpus):
	pass
def test_pblat(gfile, rfile, cpus):
	pass
def test_segemehl(gfile, rfile, cpus):
	pass
def test_star(gfile, rfile, cpus):
	pass
def test_subread(gfile, rfile, cpus):
	pass

parser = argparse.ArgumentParser(description='compare alignment programs')
parser.add_argument('genome', help='genome file *.fa.gz')
parser.add_argument('--rlen', type=int, default=500,
	help='read length [%(default)i]')
parser.add_argument('--x', type=float, default=0.1,
	help='x-fold coverage of genome [%(default).1f]')
parser.add_argument('--cpus', type=int, default=4, help='[%(default)i]')
parser.add_argument('--build', default='build', help='[%(default)s]')
parser.add_argument('--seed', type=int, default=1)
parser.add_argument('--testing', action='store_true')
arg = parser.parse_args()

## set up GLOBALS
DIR = f'{arg.build}/contest'
if not os.path.exists(f'{DIR}/genome.fa'):
	os.system(f'mkdir -p {DIR}')
	os.system(f'gunzip -c {arg.genome} > {DIR}/genome.fa')

match platform.system():
	case 'Linux':  TIMER = '/usr/bin/time -v'
	case 'Darwin': TIMER = '/usr/bin/time -l'
	case _: sys.exit('error, OS not supported')

## main loop
random.seed(arg.seed)
for err_rate in (1, 10):
	gfile = f'{DIR}/genome.fa'
	rfile = f'{DIR}/reads.{arg.seed}.{err_rate}.fa'
	simulate_reads(gfile, rfile, arg.rlen, arg.x, err_rate)
	#test_bbmap(gfile, rfile, arg.cpus)
	print(err_rate, 'blast', test_blast(gfile, rfile, arg.cpus), sep='\t')
	#test_bowtie2(gfile, rfile, arg.cpus)
	#test_bwa(gfile, rfile, arg.cpus)
	#test_gmap(gfile, rfile, arg.cpus)
	#test_hisat2(gfile, rfile, arg.cpus)
	#test_minimap2(gfile, rfile, arg.cpus)
	#test_pblat(gfile, rfile, arg.cpus)
	#test_segemehl(gfile, rfile, arg.cpus)
	#test_star(gfile, rfile, arg.cpus)
	#test_subread(gfile, rfile, arg.cpus)
