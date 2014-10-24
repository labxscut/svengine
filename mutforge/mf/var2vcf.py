#!/usr/bin/env python
import csv,time,argparse 
import pysam
import pybedtools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
#properly import muforge
try:
		import mf.mutforge as mf #in debugging
except ImportError:
		import mutforge as mf #in production

def run():
	parser = argparse.ArgumentParser(description="converting mutforge varfile to vcffile, \
                                                ")
	parser.add_argument('-t', '--filetype', dest='filetype', default='bs',
			                        help="which format to convert to ")
	parser.add_argument('-p', '--prefix', dest='prefix', default='example',
			                        help="in/out prefix")
	parser.add_argument('-r', '--reffile', dest='reffile', default='example.fna',
			                        help="reffile")
	parser.add_argument('-v', '--varfile', dest='varfile', default='example.var',
			                        help="varfile")
	parser.add_argument('-f', '--plansize', dest='plansize', default=100000,
			                        help="palnsize + ligation*2 is the grain size for genome planning in mutforge, (default=100000bp)")
	args = parser.parse_args()
	main(args)

def main(args):
	ft = args.filetype; pf = args.prefix; rf = args.reffile; vf = args.varfile; ps = args.plansize
	if ft == 'bs':
		vcffile = csv.writer(open(pf+".bs",'w'),delimiter='\t') 
		for line in mf.CommentedFile(open(vf,'rU')):
			#var_ch, var_st, var_ed, var_tp, var_sz, var_nt ...
			var_ch = line[4]
			var_st = line[5]
			var_sz = int(line[7])
			var_tp = line[1].split('_')[0]
			if var_tp in ("TRA","ITRA"):
				var_tp = "DEL" if line[6] == 'True' else "INS"
			if var_tp in ("DINS","FINS","INS"):
				var_tp = 'INS'; info=mf.readins(line[9]); var_sz=Seq(pysam.Fastafile(info[0]).fetch(info[1],info[2],info[2]+info[3]),generic_dna)
				if info[5]:
					var_sz = var_sz.reverse_complement();
				var_sz = Seq(str(var_sz)*info[4], generic_dna)
				var_ed = int(var_st)+ps
			if var_tp in ("DUP"): #IDUP is currently converted wrong
				var_ed = int(var_st)+var_sz
				info = mf.readins(line[9])
				var_sz = info[4] 
			if var_tp in ("INV"):
				var_ed = int(var_st)+var_sz
			if var_tp in ("DEL"):
				var_ed = int(var_st)+var_sz
				var_sz = "1.0"
			vcffile.writerow( [var_ch, var_st, var_ed, var_tp, var_sz] + line )
	else:
		raise "filetype not implemented yet" % ft

if __name__ == '__main__':
	starttime=time.time()
	run()
