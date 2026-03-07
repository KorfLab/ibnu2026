import argparse
import os
import platform
import random
import re
import subprocess
import sys

import korflab

def run(cli):
	pid = os.getpid()
	tempfile = f'{DIR}/temp.{pid}'
	print(cli, file=sys.stderr)
	args = f'{TIMER} {cli}'.split()
	
	with open(tempfile, 'w') as fp:
		result = subprocess.run(args, stderr=fp)
	if result.stdout is not None: sys.exit('why is there stdout?')
	if result.returncode != 0: sys.exit(f'FAILED, see {tempfile}')

	runstats = {}
	with open(tempfile) as fp:
		if TIMER == '/usr/bin/time -v': # Linux
			for line in fp:
				if re.match(r'\s+User', line):
					runstats['utime'] = float(line.split()[-1])
				elif re.match(r'\s+System', line):
					runstats['stime'] = float(line.split()[-1])
				elif re.match(r'\s+Maximum', line):
					runstats['rsize'] = int(line.split()[-1])
		elif TIMER == '/usr/bin/time -l': # Darwin
			for line in fp:
				m = re.match(r'\s+(\S+) real\s+(\S+) user\s+(\S+) sys', line)
				if m:
					runstats['etime'] = float(m.group(1))
					runstats['utime'] = float(m.group(2))
					runstats['stime'] = float(m.group(3))
				m = re.match(r'\s+(\d+)\s+maximum resident', line)
				if m:
					runstats['rsize'] = int(m.group(1))
	os.remove(tempfile)
	return runstats

def simulate_reads(infile, outfile, rlen, xcov, rate):
	rate /= 100
	truth = []
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
					truth.append( (chrom, pos+1, pos+rlen, '-') )
				else:
					rdef += '+'
					truth.append( (chrom, pos+1, pos+rlen, '+') )
				print(rdef, ''.join(rseq), sep='\n', file=fp)
	return truth

def overlap(c1, b1, e1, c2, b2, e2):
	if c1 != c2: return False
	if e1 < b2 or e2 < b1: return False
	return True

# the question of how much is _supposed_ to align is not addressed...
def proc_alignments(aligns, truth):

	check = {}
	for chrom, beg, end, strand in truth:
		tag = f'{chrom}:{beg}-{end}{strand}'
		check[tag] = None

	for query, c2, b2, e2 in aligns:
		if query not in check: sys.exit('bad programmer')
		if check[query] is not None: continue # use first/best alignment found
		check[query] = (c2, int(b2), int(e2))

	data = []
	for tag in check:
		m = re.match(r'(\w+):(\d+)\-(\d+)', tag)
		c1 = m.group(1)
		b1 = int(m.group(2))
		e1 = int(m.group(3))
		
		if check[tag] is None:
			data.append(0)
			continue # missed entirely
		
		c2, b2, e2 = check[tag]
		if not overlap(c1, b1, e1, c2, b2, e2):
			data.append(0)
			continue # paralogs...

		num = min(e1, e2) - max(b1, b2) +1
		den = max(e1, e2) - min(b1, b2) +1
		pct = num / den
		data.append(pct)
	
	return data

def sam2alignments(file):
	alignments = []
	for sam in korflab.readsam(file):
		chrom = sam.chrom.split()[0]
		if '_part_2' in sam.qname: continue # bbmap weirdness
		alignments.append( (sam.qname, chrom, sam.beg, sam.end) )
	return alignments

#####################
## WORKING TESTERS ##
#####################

def test_bbmap(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n bbmap bbmap.sh in={rfile} ref={gfile} nodisk=t threads={cpus} out={rfile}.bbmap')
	raw = proc_alignments(sam2alignments(f'{rfile}.bbmap'), truth)
	return {
		'cpu': r1['utime'] + r1['stime'] ,
		'mem': r1['rsize'],
		'raw': raw,
	}

def test_blast(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n blast-legacy formatdb -p F -i {gfile}')
	r2 = run(f'conda run -n blast-legacy blastall -p blastn -d {gfile} -i {rfile} -a {cpus} -r 1 -q -1 -G 2 -E 1 -e 1e-10 -m 8 -o {rfile}.blastn')
	alignments = []
	with open(f'{rfile}.blastn') as fp:
		for line in fp:
			f = line.split()
			query = f[0]
			chrom = f[1]
			beg = int(f[8])
			end = int(f[9])
			if beg > end: beg, end = end, beg
			alignments.append( (query, chrom, beg, end) )
	raw = proc_alignments(alignments, truth)
	return {
		'cpu': r1['utime'] + r1['stime'] + r2['utime'] + r2['stime'],
		'mem': max(r1['rsize'], r2['rsize']),
		'raw': raw,
	}

def test_bwa(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n bwa bwa index {gfile}')
	r2 = run(f'conda run -n bwa bwa mem -a -t {cpus} {gfile} {rfile} > {rfile}.bwa')
	raw = proc_alignments(sam2alignments(f'{rfile}.bwa'), truth)
	return {
		'cpu': r1['utime'] + r1['stime'] + r2['utime'] + r2['stime'],
		'mem': max(r1['rsize'], r2['rsize']),
		'raw': raw,
	}

def test_hisat2(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n hisat2 hisat2-build -f {gfile} {gfile} --quiet')
	r2 = run(f'conda run -n hisat2 hisat2 -x {gfile} -U {rfile} -f -p {cpus} --quiet > {rfile}.hisat2')
	raw = proc_alignments(sam2alignments(f'{rfile}.hisat2'), truth)
	return {
		'cpu': r1['utime'] + r1['stime'] + r2['utime'] + r2['stime'],
		'mem': max(r1['rsize'], r2['rsize']),
		'raw': raw,
	}

def test_minimap2(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n minimap2 minimap2 -a {gfile} {rfile} -t {cpus} > {rfile}.minimap2')
	raw = proc_alignments(sam2alignments(f'{rfile}.minimap2'), truth)
	return {
		'cpu': r1['utime'] + r1['stime'] ,
		'mem': r1['rsize'],
		'raw': raw,
	}

def test_segemehl(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n segemehl segemehl.x -x {gfile}.idx -d {gfile}')
	r2 = run(f'conda run -n segemehl segemehl.x -i {gfile}.idx -d {gfile} -q {rfile} -t {cpus} -o {rfile}.segemehl')
	raw = proc_alignments(sam2alignments(f'{rfile}.segemehl'), truth)
	return {
		'cpu': r1['utime'] + r1['stime'] + r2['utime'] + r2['stime'],
		'mem': max(r1['rsize'], r2['rsize']),
		'raw': raw,
	}

def test_subread(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n subread subread-buildindex -o {gfile} {gfile}')
	r2 = run(f'conda run -n subread subread-align -i {gfile} -r {rfile} -t 0 --SAMoutput --multiMapping -B 5 -T {cpus} -o {rfile}.subread')
	raw = proc_alignments(sam2alignments(f'{rfile}.subread'), truth)
	return {
		'cpu': r1['utime'] + r1['stime'] + r2['utime'] + r2['stime'],
		'mem': max(r1['rsize'], r2['rsize']),
		'raw': raw,
	}

#########################
## PROBLEMATIC TESTERS ##
#########################

def test_bowtie2(gfile, rfile, cpus):
	# requires fastq
	pass

def test_gmap(gfile, rfile, cpus):
	r1 = run(f'conda run -n gmap gmap_build -d {gfile}-gmap -D . {gfile}')
	r2 = run(f'conda run -n gmap {rfile} -d {gfile}-gmap -D . -f samse -t {cpus} > {rfile}.gmap')
	cpu = r1['utime'] + r1['stime'] + r2['utime'] + r2['stime']
	mem = max(r1['rsize'], r2['rsize'])
	return align_stats(sam2alignments(f'{rfile}.gmap')), cpu, mem

def test_pblat(gfile, rfile, cpus):
	# fails on MacOS with Bus error when threads > 1
	r1 = run(f'conda run -n pblat pblat {gfile} {rfile} {rfile}.pblat -threads={cpus}')
	cpu = r1['utime'] + r1['stime']
	mem = r1['rsize']
	return align_stats(sam2alignments(f'{rfile}.pblat')), cpu, mem

def test_star(gfile, rfile, cpus):
	# fails on MacOS, unable to run cat
	r1 = run(f'conda run -n star STAR --runMode genomeGenerate --genomeDir {gfile}-star --genomeFastaFiles {gfile} --genomeSAindexNbases 8')
	r2 = run(f'conda run -n star STAR --genomeDir {gfile}-star --readFilesIn {rfile} --readFilesCommand cat --outFileNamePrefix {rfile}.star --runThreadN {cpus}')
	os.rename(f'{rfile}.starAligned.out.sam', f'{rfile}.star')
	cpu = r1['utime'] + r1['stime'] + r2['utime'] + r2['stime']
	mem = max(r1['rsize'], r2['rsize'])
	return align_stats(sam2alignments(f'{rfile}.star')), cpu, mem


#########
## CLI ##
#########

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

## set up build directory
DIR = f'{arg.build}/contest'
if not os.path.exists(f'{DIR}/genome.fa'):
	os.system(f'mkdir -p {DIR}')
	os.system(f'gunzip -c {arg.genome} > {DIR}/genome.fa')
DIR = os.path.abspath(DIR)
print(f'Working directory: {DIR}', file=sys.stderr)

## set up TIMER
match platform.system():
	case 'Linux':  TIMER = '/usr/bin/time -v'
	case 'Darwin': TIMER = '/usr/bin/time -l'
	case _: sys.exit('error, OS not recognized')

## main loop
random.seed(arg.seed)
tests = (
	('bbm', test_bbmap),
	('bst', test_blast),
#	('bwa', test_bwa),
#	('ht2', test_hisat2),
#	('mm2', test_minimap2),
#	('seg', test_segemehl),
#	('sub', test_subread),

#	('bt2', test_bowtie2),
#	('gmp', test_gmap),
#	('pbt', test_pblat),
#	('str', test_star),
)
for err_rate in range(3):
	gfile = f'{DIR}/genome.fa'
	rfile = f'{DIR}/reads.{arg.seed}.{err_rate}.fa'
	truth = simulate_reads(gfile, rfile, arg.rlen, arg.x, err_rate)
	detail = {}
	for prog, tester in tests:
		result = tester(gfile, rfile, arg.cpus, truth)
		detail[prog] = result['raw']
		print(err_rate, prog, result['cpu'], result['mem'], sep='\t')

	for i in range(len(truth)):
		for prog in detail:
			print(i, prog, detail[prog][i])
		print()
