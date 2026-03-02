import argparse
import korflab

parser = argparse.ArgumentParser()
parser.add_argument('file', help='fastq')
arg = parser.parse_args()

fp = korflab.getfp(arg.file)
while True:
	try: header = next(fp)
	except: break
	if header == '': break
	seq = next(fp)
	plus = next(fp)
	qual = next(fp)
	uid = '>' + header[1:].split()[0]
	print(uid, seq, sep='\n', end='')
