#!/usr/bin/env python

import re, os, sys, random, csv, subprocess, traceback, argparse, tempfile, itertools, operator, collections, ConfigParser, json, imp
from itertools import izip
from collections import Counter
from multiprocessing import Pool
import pysam, pybedtools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

io_quote = "this is for end user"
code_quote = 'this is for program'
maxrefnum = 24 		  #for debugging and human, should set to a big number afterwards

class CommentedFile:
	def __init__(self, f, commentstring="#"):
		self.f = csv.reader(f,delimiter='\t')
		self.commentstring = commentstring
	def next(self):
		line = self.f.next()
		while line[0].startswith(self.commentstring):
			line = self.f.next()
		return line
	def __iter__(self):
		return self

class Bunch(object):
  def __init__(self, adict):
    self.__dict__.update(adict)

#paradigm to make dictionary from local variables
#for i in ('apple', 'banana', 'carrot'):
#    fruitdict[i] = locals()[i]
#paradigm to generate random sequence
#''.join([random.choice('AGTC') for x in range(1200)])
#paradigm to get fasta file size
#pysam.Fastafile('example.fna').lengths
#fas2bed(fasfile)

def CountUnique(mylist):
	h = {}
	for x in mylist: h[x] = h.get(x, 0) + 1
	return h

def isint(s):
	try: int(s); return True;
	except ValueError: return False;

def locatemin(a):
	smallest = min(a)
	return smallest, [index for index, element in enumerate(a) if smallest == element]

def demultiplex(m):
	return [i for i, x in enumerate([int(x) for x in bin(m)[2:][::-1]]) if x == 1]

def main(args):

	print "successfully debugging loading"

	runmode=0 #1 fasta+meta, 2 bam+meta, 3 fasta+var, 4 bam+var
	if args.inbamfilenames!=None and args.metafile!=None:
		runmode=2;
	elif args.inbamfilenames!=None and args.varfile!=None:
		runmode=4;
	elif args.hapfiles!=None and args.metafile!=None:
		runmode=1;
	elif args.hapfiles!=None and args.varfile!=None:
		runmode=3;

	if runmode==0:
		sys.stderr.write("inccorrect input specification, please use on of following: \n \
													(1) fasta+meta, (2) bam+meta, (3) fasta+var, (4) bam+var")

	hapfiles = args.hapfiles.split(':') 		#this is required input
	gapfile = args.gapfile 							  #this is required input
	parfile = args.parfile								#parfile controls WGS
	reffile = args.reffile 							  #this is required input
	nploid = args.nploid 									#this input has default
	plansize = args.plansize							#this input has default
	nligation = args.nligation						#this input has default
	nprocs = args.nprocs
	parlist = ConfigParser.ConfigParser()
	parlist.readfp(parfile)
	oprefix = 'example' if not args.oprefix else args.oprefix
	

	if len(hapfiles) == 1: #reuse hapfiles for haplotype
		hapseq=[pysam.Fastafile(hapfiles[0])]*nploid
	elif len(hapfiles)==nploid:
		hapseq=[pysam.Fastafile(hapfiles[i]) for i in range(0,nploid)]
		for ref in freseq:
			assert ref.lengths == hapseq[0].lengths, "hapfile seqs must be of the same length cross haplos"
	else:
		raise Exception("incorrect number of hapfiles provided")

	try:
		gaptab=pybedtools.BedTool(gapfile.read(),from_string=True)
	except pybedtools.cbedtools.BedToolsFileError:
		raise Exception("incorrect gapfile provided, must in BED format")

	dic={} #common dictionary interface for passing parameters
	for i in ('hapseq','gaptab','reffile','plansize','nligation','nploid','oprefix','nprocs'):
		dic[i] = locals()[i]

	if runmode==1 or runmode==2:
		metafile = args.metafile
		metatab = CommentedFile(metafile)
		dic['metatab']=metatab
		varlist = makevar(Bunch(dic))
		dic['varlist']=varlist
		var2file(Bunch(dic))

	if runmode==3 or runmode==4:
		dic['varfile'] = args.varfile
		varlist = file2var(Bunch(dic))
  	# FIXME: file2var needs to be tested?

	if runmode==1 or runmode==3:
		dic['parlist'] = parlist
		fbam = var2wgs(Bunch(dic)) # run wgsim + bwa + samtools to generate bamfiles
		print fbam

	if runmode==2 or runmode==4:
		inbamfilenames = args.inbamfilenames.split(':')
		var2bs(Bunch(dic))

	return None #done

def var2bs(bun):
	numinbamfiles = len(inbamfilenames)
	oprefix = args.oprefix.split(':') if args.oprefix!=None else args.oprefix.split(':')
	return None

def fas2bed(fasfile):
	fasseq = pysam.Fastafile(fasfile)
	bed = [ fasseq.references[i]+"\t0\t"+str(fasseq.lengths[i]) for i in range(0,min(fasseq.nreferences,maxrefnum)) ]
	return pybedtools.BedTool('\n'.join(bed),from_string=True)

def fas2size(fasfile):
	sizefile = tempfile.NamedTemporaryFile('w',delete=False)
	fasseq = pysam.Fastafile(fasfile)
	size = [ fasseq.references[i]+'\t'+str(fasseq.lengths[i]) for i in range(0,min(fasseq.nreferences,maxrefnum)) ]
	sizefile.write('\n'.join(size))
	sizefile.close()
	return sizefile.name

def chkvar(varbed):
	if varbed.merge().count() == varbed.count(): return True;
	else: return False

def ppbam(bun,bamfiles):
	
	return None

def var2wgs(bun):
	hapseq=bun.hapseq
	nprocs=bun.nprocs
	oprefix=bun.oprefix
	parlist=bun.parlist
	nligation=bun.nligation
	nlibrary=bun.parlist.getint('wgsim', 'nlibrary')
	libnames=json.loads(parlist.get('wgsim', 'libnames'))
	varlist=bun.varlist
	#parfile in python ConfigParser format
	#LibraryNumber = 3
	#Have to use very complex merging
	
	#1, construct continuous contig sets with ligations (nploid)*(variant/non-variant)
	#hapvar{0: bed_features}
	hapvar = collections.OrderedDict(); hapbed = collections.OrderedDict();
	#print varlist.values()
	for var in varlist.values():
		#print var
		#print [var[4],var[5],var[5]+var[7],var[0],var[1],var[2],var[6],var[8]]
		#CHR ST	ED	VID	MID	FREQ	DEL	INS; adding one base size to ensure flanking regions are separate
		hapvar[var[3]] = hapvar.get(var[3],[]) + [[var[4],var[5],var[5]+var[7]+2*nligation+1,var[0],var[1],var[2],var[7],var[9]]]
	#print hapvar
	for hap in hapvar.keys():
		hapvar[hap] = pybedtools.BedTool('\n'.join([ '\t'.join([str(c) for c in x]) for x in hapvar[hap] ]), from_string=True)
#		if hap==0:
#			print "varing region: for hap", hap
#			print hapvar[hap]
		assert chkvar(hapvar[hap]) == True, "varfile contains self overlaps, which is not allowed"
	for hap in hapvar.keys():
		sizefile = fas2size(hapseq[hap].filename)												 #create genome size file
		hapbed[hap] = hapvar[hap].complement(g=sizefile) #create complementary blocks
#		if hap==0:
#			print "complement region: for hap", hap
#			print hapbed[hap].truncate_to_chrom({"22":(0,51304567)})
		hapflank = hapbed[hap].flank(g=sizefile,b=nligation) #create ligated complementary
#		if hap==0:
#			print "flank region: for complement", hap
#			print hapflank.truncate_to_chrom({"22":(0,51304567)})
		hapbed[hap] = hapflank.cat(hapbed[hap]) #connect complementary+ligation
#		if hap==0:
#			print "all dummy region: ", hap
#			print hapbed[hap].truncate_to_chrom({"22":(0,51304567)})
		hapbed[hap] = pybedtools.BedTool('\n'.join([ f[0]+'\t'+str(f[1])+'\t'+str(f[2])+'\tNone\tNone\t0\t0\tNone' for f in hapbed[hap]]),from_string=True)
#		if hap==0:
#			print "fillin region: for hap", hap
#			print hapbed[hap]
		#VID=None, MID=None, FREQ=0, DEL=0, INS=None

	libbams = [ [] ] * nlibrary #bamfiles to be merged: [ [hap1lib1, hap1lib2, ...], [hap2lib1, hap2lib2, ...] ]
	pool = Pool(processes=int(nprocs))
	for hap in hapvar.keys():
		#print hapvar[hap]
		for bed in hapvar[hap]:
			bedbam = makebam(bun,hap,bed)
			#bam = pool.apply_async(makebam, [bun,var,bedline])
			for lib in range(0,nlibrary):
				libbams[lib] += bedbam[lib]

	fbam = [ None ] * nlibrary
	for lib in range(0,nlibrary):
		fbam[lib] = mergebam(libbams[lib],'.'.join([oprefix,libnames[lib]])) #post processing: librarywise merging bams
	
	return fbam

def makebam(bun,hap,bed): # makebam from a pair of seqs of representing the original and sved haplotype
	nligation=bun.nligation
	hapfas=bun.hapseq[hap]
	seq0=SeqRecord(Seq(hapfas.fetch(bed[0],int(bed[1])+1,int(bed[2])), generic_dna),id=bed[0], name="", description="")
	seq1=seq0
	#CHR ST	ED	VID	MID	FREQ	DEL	INS
	if int(bed[6])!=0: #work on deletion
		seq1.seq=seq1.seq[:nligation]+seq1.seq[nligation+int(bed[6]):];
	if bed[7]!="None": #work on insertion
		#print bed[7]
		info = readins(bed[7])
		#print info
		seqi = SeqRecord(Seq(pysam.Fastafile(info[0]).fetch(info[1],info[2],info[2]+info[3]), generic_dna),id=info[0], name="", description="")
		#['example.fasta', '22', 25202000, 1000, 2, False]
		if info[5]:
			seqi.seq = seqi.seq.reverse_complement();
		seqt=seqi
		for i in range(1,info[4]):
			seqi.seq += seqt.seq
		seq1.seq=seq1.seq[:nligation]+seqi.seq+seq1.seq[nligation:]
	parlist=bun.parlist
	nlibrary=parlist.getint('wgsim', 'nlibrary')
	libnames=json.loads(parlist.get('wgsim', 'libnames'))
	freq=float(bed[5])

	bams = [ None ] * nlibrary
	for lib in range(0,nlibrary):
		freq0 = 1-freq
		bam0 = runwgs(bun,seq0,lib,freq0)
		bam1 = runwgs(bun,seq1,lib,freq)
		bams[lib] = mergebam([bam0,bam1],tempfile.NamedTemporaryFile('w',delete=False).name)

	return bams

def runcmd(cmd):
	try:
		subprocess.Popen(cmd,shell=True)
	except OSError:
		raise Exception("cmd %s didn't run successfully" % cmd)

def runwgs(bun,seq,lib,freq): #create bam file for 1 library given seq and freq
	libnames=json.loads(bun.parlist.get('wgsim', 'libnames'))
	coverage=json.loads(bun.parlist.get('wgsim', 'coverage'))
	isize=json.loads(bun.parlist.get('wgsim', 'isize'))
	rl1=json.loads(bun.parlist.get('wgsim', 'read1'))
	rl2=json.loads(bun.parlist.get('wgsim', 'read2'))
	libn=libnames[lib]
	nligation=bun.nligation
	reffile=bun.reffile
	nrdp = (len(seq)*coverage[libn]*freq)/2+1
	nmode = len(isize[libn])
	mode = [ None ] * nmode
	fa=tempfile.NamedTemporaryFile('w',delete=False).name
	SeqIO.write(seq,fa,'fasta')
	if freq == 0 or nrdp == 0:
		return None
	else:
		for m in range(0,nmode):
			mden=isize[libn][m][0]; mnread = int(nrdp * mden)+1
			misize=isize[libn][m][1]; misizesd=isize[libn][m][2]
			fq1=tempfile.NamedTemporaryFile('w',delete=False).name
			fq2=tempfile.NamedTemporaryFile('w',delete=False).name
			mbam=tempfile.NamedTemporaryFile('w',delete=False).name
			#print (misize,misizesd,mnread,rl1[libn],rl2[libn],fa,fq1,fq2)
			wgs_cmd="wgsim -h -d %d -s %d -N %d -1 %d -2 %d %s %s %s" % (misize,misizesd,mnread,rl1[libn],rl2[libn],fa,fq1,fq2)
			#print wgs_cmd
			runcmd(wgs_cmd) #fastq file
			#print (fq1,fq2,nligation,len(seq)-nligation)
			fq3,fq4=filterfq(fq1,fq2,nligation,len(seq)-nligation) #half ligation region reads
			bwasw_cmd="bwa bwasw %s %s %s | samtools view -uS - | samtools sort - %s" % (reffile.name,fq3,fq4,mbam)
			#subprocess.Popen(bwasw_cmd)  #bam file
			runcmd(bwasw_cmd)
			#print bwasw_cmd
			#TODO: for tmpfile in (fq1,fq2,fq3,fq4): os.remove(tmpfile)
			mode[m]=mbam+'.bam'
	libbam = mergebam(mode,tempfile.NamedTemporaryFile('w',delete=False).name)
	#TODO: for tmpfile in mode + [fa]: os.remove(tmpfile)
	return libbam

def filterfq(fq1,fq2,leftcut,rightcut): #half ligation region reads
	fq3seq = []; fq4seq=[]; thincnt=0
	for rfq1,rfq2 in zip(SeqIO.parse(open(fq1,'rU'),'fastq'),SeqIO.parse(open(fq2,'rU'),'fastq')):
		assert rfq1.id[:-1] == rfq2.id[:-1], "unmatching fastq1 and fastq2, error in simulation"
		#@fins.seq1_461_964_1:0:0_1:0:0_0/1
		rdpl = int(rfq1.id.split('_')[1]); rdpr = int(rfq1.id.split('_')[2])
		#print rdpl, rdpr
		if (rdpl < leftcut and rdpr < leftcut) or (rdpl > rightcut and rdpr > rightcut):
			#TODO: currently thin when both in ligation region, future may need check not thin if at chrom start/end
			thincnt+=1
			if thincnt % 2 == 1: #odd counts not thinned, thin exactly 1/2
				fq3seq.append(rfq1); fq4seq.append(rfq2)
		else:
			fq3seq.append(rfq1); fq4seq.append(rfq2)
	fq3=tempfile.NamedTemporaryFile('w',delete=False).name
	fq4=tempfile.NamedTemporaryFile('w',delete=False).name
	SeqIO.write(fq3seq,fq3,'fastq')
	SeqIO.write(fq4seq,fq4,'fastq')
	return fq3,fq4 #NOTE: these are freeed in upper level

def mergebam(bams,oprefix):
	bams = [x for x in bams if x is not None] #remove empty components
	nbams = len(bams); fbam=oprefix;
	#samtools merge can take up to 4092 files in one merge
	for i in range(0,nbams/4000+1):
		if i == nbams/4000:
			mbam=fbam
		else:
			os.remove(mbam)
			mbam=tempfile.NamedTemporaryFile('w',delete=False).name
		merge_cmd = "samtools merge %s %s" % (mbam, " ".join(bams[i*4000:min((i+1)*4000,nbams)]))
		print merge_cmd
		runcmd(merge_cmd)
		#subprocess.Popen(bwasw_cmd)  #bam file
	return fbam+".bam"

def file2var(bun):
	ifile=csv.reader(CommentedFile(bun.varfile)); hapseq=bun.hapseq; nligation=bun.nligation
	varlist=OrderedDict()
	for row in ifile: valist[row[0]]=row[0];
	for var in varlist.keys():
		varlist[var][2] = float(varlist[var][2])
		varlist[var][3] = int(varlist[var][3])
		varlist[var][5] = int(varlist[var][5])-nligation
		varlist[var][6] = bool(varlist[var][6])
		varlist[var][7] = int(varlist[var][7])
		varlist[var][8] = bool(varlist[var][8])
		varlist[var][9] = readins(varlist[var][9])
	return varlist

def readins(cell):
	if cell == 'None': return None;
	tmp = cell.split(','); chr = tmp[1].split(':')[0]; copy = int(tmp[2])
	start = int(tmp[1].split(':')[1].split('-')[0])-1; span = int(tmp[1].split(':')[1].split('-')[1])-start
	return [ tmp[0], chr, start, span, copy, True if tmp[3] == 'r' else False ]

def formatins(cell):
	#cell = ['example.fasta', '22', 25202000, 1000, 2, False] #0-based index cell[1] and 1-based cell[2]
	#print cell
	if not cell: return 'None';
	revcomp = 'r' if cell[5] else 'f'
	coord = cell[1]+':'+str(cell[2]+1)+'-'+str(cell[2]+cell[3])
	return ','.join([cell[0],coord,str(cell[4]),revcomp])

def var2file(bun):
	varlist=bun.varlist; hapseq=bun.hapseq; nligation=bun.nligation; ofile=csv.writer(open(bun.oprefix+'.var','w'),delimiter='\t')
	ofile.writerow(['#VID','MID','VARFREQ','VARHAP','VARCHR','VARPOS','VARDEL','VARDELSIZE','VARINS','VARINSSEQ(HAP/SEQFILE,CHR:START-END,COPY,REVCOMP)'])
	for var in varlist.keys():
		varlist[var][9]=formatins(varlist[var][9])
	for var in varlist.keys():
		varlist[var][5]=varlist[var][5]+nligation
	ofile.writerows(varlist.values())
	return None

def makevar(bun):
	metatab = bun.metatab
	gaptab = bun.gaptab
	hapseq = bun.hapseq
	nploid = bun.nploid
	plansize = bun.plansize
	nligation = bun.nligation
	trunks = maketrunks(bun)   													#trunks is in bed format
	trunks,stretches = blocktrunks(trunks,gaptab)  #block trunks that not usable
	
	metacnt={}; varcnt=1; varlist=collections.OrderedDict()
	
	rows=[]
	for row in metatab:
		rows.append(row)
	rowix=range(0,len(rows))
	#random.shuffle(rowix) #randomized between metatypes
	
	for r in rowix:
		sizevec=[]; copyvec=[]; donorvec=[]; recptvec=[]; donorfreq=[]; recptfreq=[]
		row = rows[r]; metatype = row[0]; cnt = int(row[1])
		
		if metatype != 'FINS':
			sizedist = expanddist(row[2],cnt)
		else:
			#print "foreignseq=", row[2]
			foreignseq = pysam.Fastafile(row[2])
			if(foreignseq.nreferences!=cnt):
				raise Exception("number of foreign sequences provided not equal to specified in FINS!")
			else:
				sizedist = CountUnique(foreignseq.lengths)

		copydist = expanddist(row[3],cnt)
		donordist = expanddist(row[4],cnt)
		recptdist = expanddist(row[5],cnt)
		donorfdist = expanddist(row[6],cnt)
		recptfdist = expanddist(row[7],cnt)
		for key in sizedist.keys():
			sizevec.extend([key]*sizedist[key])
		for key in copydist.keys():
			copyvec.extend([key]*copydist[key])
		for key in donordist.keys():
			donorvec.extend([key]*donordist[key])
		for key in recptdist.keys():
			recptvec.extend([key]*recptdist[key])
		for key in donorfdist.keys():
			donorfreq.extend([key]*donorfdist[key])
		for key in recptfdist.keys():
			recptfreq.extend([key]*recptfdist[key])
		### randomized within metatypes ###
		random.shuffle(copyvec)
		random.shuffle(donorvec)
		random.shuffle(recptvec)
		random.shuffle(donorfreq)
		random.shuffle(recptfreq)
		assert(len(sizevec)==cnt)
		assert(len(copyvec)==cnt)
		assert(len(donorvec)==cnt)
		assert(len(recptvec)==cnt)
		assert(len(donorfreq)==cnt)
		assert(len(recptfreq)==cnt)

		def maketid(stretches,varstretch):
			#print stretches
			#print varstretch
			tids = [t for t, x in enumerate(stretches) if x >= varstretch] #always random pick among min qualifying tids
			if tids == []: raise Exception("cannot allocate genome region for metafile, please reduce plansize!");
			mtids = locatemin([stretches[i] for i in tids])
			#print mtids
			tid = random.sample([tids[i] for i in mtids[1]],1)[0]   #locate the tid
			#print tid
			for i in xrange(tid,tid+varstretch): stretches[i]=0;						#mark head of stretch used
			if mtids[0]-varstretch>0: stretches[tid+varstretch]=mtids[0]-varstretch; #mark remaining stretch by its head
			return tid, stretches
			
		for i in xrange(0,cnt): #this could be parallelized? yes if we preallocate trunks
			metacnt[metatype] = metacnt.get(metatype,0) + 1
			mid = metatype + str(metacnt[metatype])
			#print metatype
			donor=demultiplex(donorvec[i]) #this is a list of donor haplos
			recpt=demultiplex(recptvec[i]) #this is a list of recpt haplos
			dpick = random.sample(donor,1)[0] #have to randomly choose a donor sequence
			vartag = 'VAR'
			
			if metatype == 'DEL':    #this should be only operations to donors
				dtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
				for d in donor:
					varid = vartag+str(varcnt)
					varlist[varid]=[varid,mid,donorfreq[i],d,trunks[dtid].chrom,trunks[dtid].start,True,sizevec[i],False,None]
					varcnt+=1
			elif metatype == 'TRP' or metatype == 'ITRP':  #one deletion at donor followed by a insertion at recepient
				rc = True if metatype == 'ITRP' else False
				dtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
				for d in donor:
					varid = vartag+str(varcnt)
					varlist[varid]=[varid,mid,donorfreq[i],d,trunks[dtid].chrom,trunks[dtid].start,True,sizevec[i],False,None]
					varcnt+=1
				rtid,stretches=maketid(stretches,1)
				for r in recpt:    #inserting sequence: [ filename, chr, start, size, copy, rc=True/False ]
					varid = vartag+str(varcnt)
					varlist[varid]=[varid,mid,recptfreq[i],r,trunks[rtid].chrom,trunks[rtid].start,False,0,True,\
						[hapseq[dpick].filename,trunks[dtid].chrom,trunks[dtid].start,sizevec[i],copyvec[i],rc]]
					varcnt+=1
			elif metatype == 'INV' or metatype == 'IDUP':  #one deletion at donor followed by insertion of reverse
				dtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
				for d in donor:
					varid = vartag+str(varcnt)
					varlist[varid]=[varid,mid,donorfreq[i],d,trunks[dtid].chrom,trunks[dtid].start,True,sizevec[i],True,
						[hapseq[dpick].filename,trunks[dtid].chrom,trunks[dtid].start,sizevec[i],copyvec[i],True]]
					varcnt+=1
			elif metatype == 'DINS':
				dtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
				rtid,stretches=maketid(stretches,1)
				for r in recpt:    #inserting sequence: [ filename, chr, start, size, copy, rc=True/False ]
					varid = vartag+str(varcnt)
					varlist[varid]=[varid,mid,recptfreq[i],r,trunks[rtid].chrom,trunks[rtid].start,False,0,True,\
						[hapseq[dpick].filename,trunks[dtid].chrom,trunks[dtid].start,sizevec[i],copyvec[i],False]]
					varcnt+=1
			elif metatype == 'DUP':
				rtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
				for r in recpt:    #inserting sequence: [ filename, chr, start, size, copy, rc=True/False ]
					varid = vartag+str(varcnt)
					varlist[varid]=[varid,mid,recptfreq[i],r,trunks[rtid].chrom,trunks[rtid].start,False,0,True,\
						[hapseq[r].filename,trunks[rtid].chrom,trunks[rtid].start,sizevec[i],copyvec[i]-1,False]]
					varcnt+=1
			elif metatype == 'FINS':
				#print donor, recpt, dpick
				rtid,stretches=maketid(stretches,1)
				#print foreignseq.filename
				#print foreignseq.references[i]
				for r in recpt:    #inserting sequence: [ filename, chr, start, size, copy, rc=True/False ]
					varid = vartag+str(varcnt)
					#print foreignseq[i]
					varlist[varid]=[varid,mid,recptfreq[i],r,trunks[rtid].chrom,trunks[rtid].start,False,0,True,\
						[foreignseq.filename,foreignseq.references[i],0,foreignseq.lengths[i],copyvec[i],False]]
					varcnt+=1
			else:
				raise Exception("unknown metatype %s" % metatype)

	return varlist #synonym to varlist

def blocktrunks(trunks,gaptab):
	#print(trunks)
	#print(gaptab)
	remains = trunks.intersect(gaptab,v=True) #only keeping non overlapping trunks
	#FIXME: sometimes bedtool error will rise if temporary files were deleted, how to fix in safeway?
	#print(remains) #this has to be loaded to memory
	cnt=1; stretches = [0] * len(remains)
	for j in xrange(1,len(remains)):
		if remains[j].chrom==remains[j-1].chrom and remains[j].start == remains[j-1].end:
			cnt+=1 #continue counting from j-1 to j
		else:    #stop counting, j-1 in stretch starts from , j is not in
			stretches[j-cnt]=cnt
			for k in xrange(j-cnt+1,j): stretches[k] = 0
			cnt=1  #start new counting
		if j == len(remains)-1 : #this is forced ending, we need to trace back (cnt-1)
			stretches[j-cnt+1]=cnt
			for k in xrange(j-cnt+2,j): stretches[k] = 0
	return remains,stretches

def maketrunks(bun):
	#one trunk will be [ sv_pos-buffer, sv_pos+maxsv+buffer ]
	#trunks are stored by [ {1:(chr1,sv_pos),2:(chr2,sv_pos),...}, {1:(chr1,sv_pos),2:(chr2,sv_pos),...}, ... ]
	#all converted to bed format
	hapseq=bun.hapseq
	plansize=bun.plansize
	nligation=bun.nligation
	trunks = [] # trunk[i].start = sv[i].pos - nligation; trunk[i+1] = trunk[i].start + plansize + 2*ligation
	for i in xrange(0,min(hapseq[0].nreferences,maxrefnum)):
		pos=xrange(0,hapseq[0].lengths[i]-plansize-2*nligation,plansize+2*nligation)
		trunks.extend([ '%s\t%s\t%s' % (hapseq[0].references[i], str(p), str(p+plansize+2*nligation)) for p in pos ])
	trunks = pybedtools.BedTool("\n".join([ x for x in trunks ]), from_string=True) #convert to BED
	return trunks

def expanddist(diststring, cnt):
	dist=diststring.split('_')
	uniques = {}
	new_cnt = 0
	if dist[0] == 'fix':
		#print dist[1:len(dist)]
		values = [ int(x) if isint(x) else float(x) for x in dist[1:len(dist)] ]
		uniques = CountUnique(values)
		sorted_uniques = sorted(uniques.items(), key=operator.itemgetter(1)) #increasing order sorting of values
		nvalue = len(values)
		nunique = len(uniques)
		amp = cnt/nvalue
		rem = cnt-amp*nvalue
		#print(sorted_uniques)
		for item in sorted_uniques:
			key = item[0]
			if rem>0:
				delta = uniques[key] if rem-uniques[key]>=0 else rem
				uniques[key] = uniques[key] * amp + delta
				rem = rem - delta
			else:
				delta = 0
				uniques[key] = uniques[key] * amp
			new_cnt = new_cnt + uniques[key]
		assert new_cnt == cnt
		#print(uniques)

	elif dist[0] == 'unif':
		raise Exception("can't expand specified distribution!")
	elif dist[0] == 'exp':
		raise Exception("can't expand specified distribution!")
	elif dist[0] == 'norm':
		raise Exception("can't expand specified distribution!")
	else:
		raise Exception("can't expand specified distribution!")
	return uniques

def run():
	parser = argparse.ArgumentParser(description="spike-in mutations to fasta/bam files, \
                                                  for bam files first assemble haplo-types \
                                                  for fasta files we segment fasta files \
                                                  to confirm if the configuration is viable")
	parser.add_argument('hapfiles', metavar='hapfiles', #example.fasta
                        help="inputfile of hap file(s), required")
	parser.add_argument('gapfile', metavar='gapfile', type=argparse.FileType('rU'), #example.bed
                        help="inputfile of gap or other banned regions, required")
	parser.add_argument('parfile', metavar='parfile', type=argparse.FileType('rU'), #example.par
                        help="inputfile of library-wise parameters, required"),
	parser.add_argument('reffile', metavar='reffile', type=argparse.FileType('rU'), #example.fasta
                        help="inputfile of ref file(s), required"),
	parser.add_argument('-m', '--metafile', dest='metafile', type=argparse.FileType('rU'), #example.meta
                        help="inputfile of meta specification")
	parser.add_argument('-v', '--varfile', dest='varfile', type=argparse.FileType('rU'),  #example.var/op
                        help="inputfile of variant specification")
	parser.add_argument('-b', '--bamfile', dest='inbamfilenames', #in1.bam,in2.bam,in3.bam
                        help="inputfile of library-wise bam file(s)")
	parser.add_argument('-o', '--oprefix', dest='oprefix', default=None, #output1,output2,output3,...
                        help="prefix for output files")
	parser.add_argument('-y', '--nploid', dest='nploid', default=2,
                        help="genome ploidity (default=2)")
	parser.add_argument('-f', '--plansize', dest='plansize', default=10000,
                        help="palnsize + ligation*2 is the grain size for genome planning, (default=10000bp)")
	parser.add_argument('-g', '--nligation', dest='nligation', default=2000,
                        help="genome trunk ligation (default=2000bp)")
	parser.add_argument('-n', '--nprocs', dest='nprocs', default=1,
                        help="split into multiple processes (default=1)")
	parser.add_argument('--layout', action='store_true', default=False, help="dry run to layout")
	args = parser.parse_args()
	main(args)

if __name__ == '__main__':
	run()
