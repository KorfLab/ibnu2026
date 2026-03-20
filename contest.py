import argparse
import os
import platform
import random
import re
import statistics
import subprocess
import sys

import korflab

def chrom_lengths(file):
	chrlen = {}
	for defline, seq in korflab.readfasta(file):
		chrid = defline.split()[0]
		chrlen[chrid] = len(seq)
	return chrlen

def run(cli):
	match platform.system():
		case 'Linux':  time = '/usr/bin/time -v'
		case 'Darwin': time = '/usr/bin/time -l'
		case _: sys.exit('error, OS not recognized')

	tempfile = f'{DIR}/temp.{os.getpid()}' # saving for debugging
	if os.system(f'{time} {cli} 2> {tempfile}') != 0:
		sys.exit(f'**FAILED**\n{cli}\nSee {tempfile}')

	runstats = {}
	with open(tempfile) as fp:
		if time == '/usr/bin/time -v': # Linux
			for line in fp:
				if re.match(r'\s+User', line):
					runstats['utime'] = float(line.split()[-1])
				elif re.match(r'\s+System', line):
					runstats['stime'] = float(line.split()[-1])
				elif re.match(r'\s+Maximum', line):
					runstats['rsize'] = int(line.split()[-1])
		elif time == '/usr/bin/time -l': # Darwin
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
	truth = {}
	uid = 0
	with open(outfile, 'w') as fp:
		for defline, seq in korflab.readfasta(infile):
			chrom = defline.split()[0]
			rnum = int(len(seq) * xcov / rlen)
			for i in range(rnum):
				uid += 1
				pos = random.randint(0, len(seq) - rlen)
				if 'NN' in seq[pos:pos+rlen]: continue
				rseq = []
				rid = f'r{uid}'
				rdef = f'>{rid}#{chrom}:{pos+1}-{pos+rlen}'
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
					truth[rid] = (chrom, pos+1, pos+rlen, '-')
				else:
					rdef += '+'
					truth[rid] = (chrom, pos+1, pos+rlen, '+')
				print(rdef, ''.join(rseq), sep='\n', file=fp)
	return truth

def overlap(c1, b1, e1, c2, b2, e2):
	if c1 != c2: return False
	if e1 < b2 or e2 < b1: return False
	return True

def proc_alignments(aligns, truths):
	d = {}
	for defline, chrom, beg, end in aligns:
		rid, query = defline.split('#')
		if rid not in d: d[rid] = []
		d[rid].append( (chrom, beg, end) )
	
	missed = 0
	paralog = 0
	found = []
	for rid in truths:
		c1, b1, e1, s1 = truths[rid]
		if rid not in d:
			missed += 1
			found.append(0)
			continue
		best_match = 0
		best_align = None
		for align in d[rid]:
			c2, b2, e2 = align
			if not overlap(c1, b1, e1, c2, b2, e2): continue
			num = min(e1, e2) - max(b1, b2) +1
			den = max(e1, e2) - min(b1, b2) +1
			pct = num / den
			if pct > best_match:
				best_match = pct
				best_align = align
		if best_align is None:
			paralog += 1
			found.append(0)
		else:
			found.append(best_match)
	
	return statistics.mean(found), missed/len(truths), paralog/len(truths)

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
	c, m, p = proc_alignments(sam2alignments(f'{rfile}.bbmap'), truth)
	return c, m, p, r1['utime'] + r1['stime'], r1['rsize']

def test_blast(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n blast-legacy formatdb -p F -i {gfile}')
	r2 = run(f'conda run -n blast-legacy blastall -p blastn -d {gfile} -i {rfile} -a {cpus} -r 1 -q -1 -G 2 -E 1 -F F -e 1e-10 -m 8 -o {rfile}.blastn')
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
	c, m, p = proc_alignments(alignments, truth)
	return c, m, p, r2['utime'] + r2['stime'], r2['rsize']

def test_bwa(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n bwa bwa index {gfile}')
	r2 = run(f'conda run -n bwa bwa mem -a -t {cpus} -o {rfile}.bwa {gfile} {rfile}')
	c, m, p = proc_alignments(sam2alignments(f'{rfile}.bwa'), truth)
	return c, m, p, r2['utime'] + r2['stime'],  r2['rsize']

def test_hisat2(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n hisat2 hisat2-build -f {gfile} {gfile} --quiet')
	r2 = run(f'conda run -n hisat2 hisat2 -x {gfile} -U {rfile} -f -p {cpus} --quiet > {rfile}.hisat2')
	c, m, p = proc_alignments(sam2alignments(f'{rfile}.hisat2'), truth)
	return c, m, p, r2['utime'] + r2['stime'], r2['rsize']

def test_minimap2(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n minimap2 minimap2 -a {gfile} {rfile} -t {cpus} -o {rfile}.minimap2')
	c, m, p = proc_alignments(sam2alignments(f'{rfile}.minimap2'), truth)
	return c, m, p, r1['utime'] + r1['stime'], r1['rsize']

def test_segemehl(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n segemehl segemehl.x -x {gfile}.idx -d {gfile}')
	r2 = run(f'conda run -n segemehl segemehl.x -i {gfile}.idx -d {gfile} -q {rfile} -t {cpus} -o {rfile}.segemehl')
	c, m, p = proc_alignments(sam2alignments(f'{rfile}.segemehl'), truth)
	return c, m, p, r2['utime'] + r2['stime'], r2['rsize']

def test_subread(gfile, rfile, cpus, truth):
	r1 = run(f'conda run -n subread subread-buildindex -o {gfile} {gfile}')
	r2 = run(f'conda run -n subread subread-align -i {gfile} -r {rfile} -t 0 --SAMoutput --multiMapping -B 5 -T {cpus} -o {rfile}.subread')
	c, m, p = proc_alignments(sam2alignments(f'{rfile}.subread'), truth)
	return c, m, p, r2['utime'] + r2['stime'], r2['rsize']

#########################
## PROBLEMATIC TESTERS ##
#########################

def test_gmap(gfile, rfile, cpus, truth):
	# fails with error: Permission denied... in a temporary directory
	if platform.system() == 'Darwin':
		return {'cpu': 0, 'mem': 0,'raw': [0] * len(truth)}

	r1 = run(f'conda run -n gmap gmap_build --genomedb genome-gmap --genomedir {DIR} {gfile}')
	r2 = run(f'conda run -n gmap {rfile} --genomedb genome-gmap --genomedir {DIR} -f samse -t {cpus} > {rfile}.gmap')
	c, m, p = proc_alignments(sam2alignments(f'{rfile}.gmap'), truth)
	return c, m, p, r2['utime'] + r2['stime'], r2['rsize']

def test_pblat(gfile, rfile, cpus, truth):
	# fails on MacOS with Bus error when threads > 1
	if platform.system() == 'Darwin': cpus = 1

	r1 = run(f'conda run -n pblat pblat {gfile} {rfile} {rfile}.pblat -threads={cpus} > {rfile}.stdout')
	alignments = []
	with open(f'{rfile}.pblat') as fp:
		for _ in range(5): line = next(fp)
		for line in fp:
			f = line.split()
			query = f[9]
			chrom = f[13]
			beg = int(f[15])
			end = int(f[16])
			alignments.append( (query, chrom, beg, end) )

	c, m, p = proc_alignments(alignments, truth)
	return c, m, p, r1['utime'] + r1['stime'], r1['rsize']

def test_star(gfile, rfile, cpus, truth):
	# fails on MacOS, unable to run cat (permission issue?)
	if platform.system() == 'Darwin':
		return {'cpu': 0, 'mem': 0,'raw': [0] * len(truth)}

	r1 = run(f'conda run -n star STAR --runMode genomeGenerate --genomeDir {gfile}-star --genomeFastaFiles {gfile} --genomeSAindexNbases 8 > /dev/null')
	r2 = run(f'conda run -n star STAR --genomeDir {gfile}-star --readFilesIn {rfile} --readFilesCommand cat --outFileNamePrefix {rfile}.star --runThreadN {cpus} > /dev/null')
	os.rename(f'{rfile}.starAligned.out.sam', f'{rfile}.star')
	c, m, p = proc_alignments(sam2alignments(f'{rfile}.star'), truth)
	return c, m, p, r2['utime'] + r2['stime'], r2['rsize']

#########
## CLI ##
#########

parser = argparse.ArgumentParser(description='compare alignment programs')
parser.add_argument('genome', help='genome file *.fa.gz')
parser.add_argument('name', help='directory name <buiild>/contest/<name>')
parser.add_argument('--rlen', type=int, default=500,
	help='read length [%(default)i]')
parser.add_argument('--x', type=float, default=1.0,
	help='x-fold coverage of genome [%(default).1f]')
parser.add_argument('--cpus', type=int, default=4, help='[%(default)i]')
parser.add_argument('--build', default='build', help='[%(default)s]')
parser.add_argument('--seed', type=int, default=1)
arg = parser.parse_args()

## set up build directory
DIR = f'{arg.build}/contest/{arg.name}'
if not os.path.exists(f'{DIR}/genome.fa'):
	os.system(f'mkdir -p {DIR}')
	os.system(f'gunzip -c {arg.genome} > {DIR}/genome.fa')
DIR = os.path.abspath(DIR)
print(f'Working directory: {DIR}', file=sys.stderr)

## programs and testers
random.seed(arg.seed)
tests = (
	('bbmap', test_bbmap),
	('blast', test_blast),
	('bwa', test_bwa),
	('hisat2', test_hisat2),
#	('gmap', test_gmap),
	('minimap2', test_minimap2),
	('pblat', test_pblat),
	('segemehl', test_segemehl),
#	('star', test_star),
	('subread', test_subread),
)

## Experiment 1: increasing error rate
cov_graph = {}
mis_graph = {}
for err in range(21):
	if err not in cov_graph: cov_graph[err] = {}
	if err not in mis_graph: mis_graph[err] = {}
	gfile = f'{DIR}/genome.fa'
	rfile = f'{DIR}/reads.{arg.seed}.{err}.fa'
	truth = simulate_reads(gfile, rfile, arg.rlen, arg.x, err)
	for prog, tester in tests:
		cov, mis, par, cpu, mem = tester(gfile, rfile, arg.cpus, truth)
		print(f'{prog}\t{cov:.3f}\t{mis:.3f}\t{par:.3f}\t{mem/1e6:.1f}\t{cpu:.1f}')
		cov_graph[err][prog] = cov
		mis_graph[err][prog] = mis

for graph, name in zip((cov_graph, mis_graph), ('coverage', 'missed')):
	with open(f'{arg.name}.{name}.tsv', 'w') as fp:
		print('err', end='\t', file=fp)
		for prog in graph[0]: print(prog, end='\t', file=fp)
		print(file=fp)

		for err in graph:
			print(err, end='\t', file=fp)
			for prog in graph[err]:
				print(f'{graph[err][prog]:.3f}', end='\t', file=fp)
			print(file=fp)

"""

## Experiment 2: badread simulated reads
badfq = f'{DIR}/reads.badread.fq'
badfa = f'{DIR}/reads.badread.fa'
os.system(f'conda run -n badread badread simulate --reference {arg.genome} --quantity {arg.x}x --length {arg.rlen},1 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 --start_adapter_seq "" --end_adapter_seq "" --seed {arg.seed} > {badfq} 2> /dev/null')
uid = 0
truth = []
chrlen = chrom_lengths(arg.genome)
with open(badfq) as ifp, open(badfa, 'w') as ofp:
	while True:
		try: header = next(ifp)
		except: break
		seq = next(ifp)
		plus = next(ifp)
		qual = next(ifp)
		f = header.split()
		uid += 1
		chrom, strand, beg_end = f[1].split(',')
		beg, end = beg_end.split('-')
		beg = int(beg)
		end = int(end)
		strand = strand[0]
		if strand == '-':
			beg = chrlen[chrom] - beg
			end = chrlen[chrom] - end
		print(f'>r{uid}#{chrom}:{beg}-{end}{strand}', file=ofp)
		print(seq, file=ofp, end='')
		truth.append( (chrom, beg, end, strand) )

for prog, tester in tests:
	result = tester(f'{DIR}/genome.fa', badfa, arg.cpus, truth)
	cpu = f'{result["cpu"]:.1f}'
	mem = f'{result["mem"]/1e6:.1f}'
	print(prog, mem, cpu, sep='\t')
	cov = statistics.mean(result['raw'])
	mis = len([x for x in result['raw'] if x == 0]) / len(truth)
	print(cpu, mem, cov, mis)
	#	cov_graph[err][prog] = cov
	#	mis_graph[err][prog] = mis

"""
