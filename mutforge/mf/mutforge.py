#!/usr/bin/env python

import copy, time, re, os, shutil, sys, random, csv, subprocess, traceback, argparse, tempfile, itertools, operator, collections, ConfigParser, json, imp
from itertools import izip
from collections import Counter
from multiprocessing import Pool
from functools import partial
import pysam, pybedtools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

io_quote = "this is for end user"
code_quote = 'this is for program'
maxrefnum = 24 		  #for debugging and human, should set to a big number afterwards
burnin = 100000		#burnin region for each chromosome

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

def runcmd(cmd):
	try:
		subprocess.call(cmd,shell=True,stderr=open('./stderr.txt','a'),stdout=open('./stdout.txt','a'))
		#subprocess.call(cmd,shell=True)
	except OSError:
		raise Exception("cmd %s didn't run successfully" % cmd)

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

def test_map(a,b):
	print a,b

def main(args):

	nprocs = args.nprocs
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
	else:
		raise Exception("must provide meta/var file and bam/hap files!")

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
	parlist = ConfigParser.ConfigParser()
	parlist.readfp(parfile)
	metafile = args.metafile
	varfile = args.varfile
	oname = metafile.name.rstrip(".meta") if metafile else varfile.name.rstrip(".var")
	oprefix = oname if not args.oprefix else args.oprefix
	

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

	print "entering runmode", runmode

	if runmode==1 or runmode==2:
		metatab = CommentedFile(metafile)
		dic['metatab']=metatab
		varlist = makevar(Bunch(dic))
		print "done makevar"
		dic['varlist']=varlist
		
		print "before", varlist
		var2file(Bunch(dic))
		print "after", varlist
		print "done var2file"


	if runmode==3 or runmode==4:
		dic['varfile'] = varfile
		varlist = file2var(Bunch(dic))
  	# FIXME: file2var needs to be tested?

	var2bed(Bunch(dic))
	print "done var2bed"
	
	#quit()

	if runmode==1 or runmode==3:
		dic['parlist'] = parlist
		fbam = var2wgs(Bunch(dic)) # run wgsim + bwa + samtools to generate bamfiles
		print fbam
		print "done var2wgs"

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
	nploid=bun.nploid
	nligation=bun.nligation
	nlibrary=bun.parlist.getint('wgsim', 'nlibrary')
	libnames=json.loads(parlist.get('wgsim', 'libnames'))
	varlist=bun.varlist
	reffile = bun.reffile.name
	parlist = bun.parlist
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
		hapvar[var[3]] = hapvar.get(var[3],[]) + [[var[4],var[5],var[5]+var[7]+2*nligation+1,var[0],var[1],var[2],var[7],formatins(var[9])]]
	#print hapvar
	for hap in hapvar.keys():
		hapvar[hap] = pybedtools.BedTool('\n'.join([ '\t'.join([str(c) for c in x]) for x in hapvar[hap] ]), from_string=True)
		assert chkvar(hapvar[hap]) == True, "varfile contains self overlaps, which is not allowed"
	for hap in hapvar.keys():
		sizefile = fas2size(hapseq[hap].filename)												 #create genome size file
		if hap==0:
			print "varing region: for hap", hap
			#print hapbed[hap].truncate_to_chrom({"22":(0,51304567)})
			print hapvar[hap]
		hapbed[hap] = hapvar[hap].complement(g=sizefile) #create complementary blocks
		if hap==0:
			print "complement region: for hap", hap
			#print hapbed[hap].truncate_to_chrom({"22":(0,51304567)})
			print hapbed[hap]
		hapbed[hap] = hapbed[hap].slop(g=sizefile,b=nligation)
		if hap==0:
			print "all dummy regions: for hap", hap
			#print hapbed[hap].truncate_to_chrom({"22":(0,51304567)})
			print hapbed[hap]
		hapbed[hap] = pybedtools.BedTool('\n'.join([ f[0]+'\t'+str(f[1])+'\t'+str(f[2])+'\tNone\tNone\t0\t0\tNone' for f in hapbed[hap]]),from_string=True)
		hapvar[hap] = hapvar[hap].cat(hapbed[hap],postmerge=False)
		if hap==0:
			print "all generating regions: for hap", hap
			print hapvar[hap]
		#VID=None, MID=None, FREQ=0, DEL=0, INS=None

	hapiter=[]
	print len(hapvar[0])
	print hapvar[0][0]

	for hap in hapvar.keys():
		hapfas=bun.hapseq[hap].filename
		for feat in hapvar[hap]:
			hapiter.append([hap,hapfas,nligation,parlist,reffile,nploid,feat])
	for it in hapiter:
		print it

	pool = Pool(processes=int(nprocs))
	bedbams=pool.map(fakebam,hapiter,1)

	libbams = {} #bamfiles to be merged: [ [hap1lib1, hap1lib2, ...], [hap2lib1, hap2lib2, ...] ]
	libcnt0 = 0; libcnt1 = 0
	for i in xrange(0,len(bedbams)):
		for lib in xrange(0,nlibrary):
			print bedbams[i][1][lib];
			if lib == 0: 
				libcnt0 += 1
			else: 
				libcnt1 += 1;
			if libbams.get(lib,None) == None:
				libbams[lib]=[bedbams[i][1][lib]]
			else:
				libbams[lib].append(bedbams[i][1][lib]) 

	print "len bedbams:", len(bedbams)
	for lib in range(0,nlibrary):
		print "len libbams:", len(libbams[lib])
	print  libcnt0,libcnt1
	#print bedbams
	#print libbams

	fbam = [None] * nlibrary
	for lib in range(0,nlibrary):
		fbam[lib] = mergebam(libbams[lib],'.'.join([oprefix,libnames[lib]])) #post processing: librarywise merging bams
	
	return fbam

def fakebam(iter): # makebam from a pair of seqs of representing the original and sved haplotype
	hap=iter[0]; hapfas=iter[1]; nligation=iter[2]; parlist=iter[3]; reffile=iter[4]; nploid=iter[5]; bed=iter[6];
	seq0=SeqRecord(Seq(pysam.Fastafile(hapfas).fetch(bed[0],int(bed[1])+1,int(bed[2])), generic_dna),id=bed[0], name="", description="")
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
	nlibrary=parlist.getint('wgsim', 'nlibrary')
	libnames=json.loads(parlist.get('wgsim', 'libnames'))
	freq=float(bed[5])

	bams = [None] * nlibrary
	for lib in range(0,nlibrary):
		freq0 = 1-freq
		print "-----"
		bam0 = runwgs(parlist,nligation,reffile,nploid,seq0,lib,freq0)
		bam1 = runwgs(parlist,nligation,reffile,nploid,seq1,lib,freq)
		print "-----"
		bams[lib] = mergebam([bam0,bam1],tempfile.NamedTemporaryFile('w',delete=False).name)

	return [hap,bams]

def runwgs(parlist,nligation,reffile,nploid,seq,lib,freq): #create bam file for 1 library given seq and freq
	libnames=json.loads(parlist.get('wgsim', 'libnames'))
	coverage=json.loads(parlist.get('wgsim', 'coverage'))
	isize=json.loads(parlist.get('wgsim', 'isize'))
	rl1=json.loads(parlist.get('wgsim', 'read1'))
	rl2=json.loads(parlist.get('wgsim', 'read2'))
	libn=libnames[lib]
	nrdp = len(seq)*coverage[libn]/nploid*freq/(rl1[libn]+rl2[libn])
	nmode = len(isize[libn])
	mode = [ None ] * nmode
	fa=tempfile.NamedTemporaryFile('w',delete=False).name
	print "run wgs for:", libn, fa, coverage[libn], freq, nrdp
	SeqIO.write(seq,fa,'fasta')
	if freq == 0 or nrdp == 0:
		return None
	else:
		for m in range(0,nmode):
			mden=isize[libn][m][0]; mnread = int(nrdp * mden)
			if mnread>0:
				misize=isize[libn][m][1]; misizesd=isize[libn][m][2]
				fq1=tempfile.NamedTemporaryFile('w',delete=False).name
				fq2=tempfile.NamedTemporaryFile('w',delete=False).name
				mbam=tempfile.NamedTemporaryFile('w',delete=False).name
				#print (misize,misizesd,mnread,rl1[libn],rl2[libn],fa,fq1,fq2)
				wgs_cmd="wgsim -h -d %d -s %d -N %d -1 %d -2 %d %s %s %s > /dev/null 2>&1" % (misize,misizesd,mnread,rl1[libn],rl2[libn],fa,fq1,fq2)
				print wgs_cmd
				runcmd(wgs_cmd) #fastq file
				fq3,fq4=filterfq(fq1,fq2,nligation,len(seq)-nligation) #half ligation region reads
				#print (fq1,fq2,fq3,fq4,nligation,len(seq)-nligation)
				bwasw_cmd="bwa bwasw %s %s %s | samtools view -uS - | samtools sort - %s > /dev/null 2>&1" % (reffile,fq3,fq4,mbam)
				print bwasw_cmd
				runcmd(bwasw_cmd)
				#TODO: for tmpfile in (fq1,fq2,fq3,fq4): os.remove(tmpfile)
				mode[m]=mbam+'.bam'
			else:
				mode[m]=None
	libbam = mergebam(mode,tempfile.NamedTemporaryFile('w',delete=False).name)
	#TODO: for tmpfile in mode + [fa]: os.remove(tmpfile)
	return libbam

def filterfq(fq1,fq2,leftcut,rightcut): #half ligation region reads
	fq3seq = []; fq4seq=[]; thincnt=0
	for rfq1,rfq2 in zip(SeqIO.parse(open(fq1,'rU'),'fastq'),SeqIO.parse(open(fq2,'rU'),'fastq')):
		assert rfq1.id[:-1] == rfq2.id[:-1], "unmatching fastq1 and fastq2, error in simulation"
		#@fins.seq1_461_964_1:0:0_1:0:0_0/1
		rdpl = int(rfq1.id.split('_')[1]); rdpr = int(rfq1.id.split('_')[2])
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
	nbams = len(bams); fbam=oprefix+".bam";
	#print nbams
	#samtools merge can take up to 4092 files in one merge
	if nbams==1:
		shutil.move(bams[0],fbam) #if merge only one file, just move it to new name
	else:
		for i in range(0,nbams/4000+1):
			if i == 0: #first iteration
				cbam = ""
				mbam = tempfile.NamedTemporaryFile('w',delete=False).name
			if i == nbams/4000: #last iteration
				if cbam != "": cbam = mbam
				mbam = fbam
			else:
				cbam = mbam
				mbam = tempfile.NamedTemporaryFile('w',delete=False).name
			merge_cmd = "samtools merge -f %s %s %s > /dev/null 2>&1" % (mbam, cbam, " ".join([bam for bam in bams[i*4000:min((i+1)*4000,nbams)]]))
			print merge_cmd
			runcmd(merge_cmd)
	#FIXME: who do we clean up this mess of temporary files here?
	return fbam

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
	#IMPORTANT: OrderedDict and List are all mutable
	varcopy=copy.deepcopy(bun.varlist); nligation=bun.nligation; ofile=csv.writer(open(bun.oprefix+'.var','w'),delimiter='\t')
	ofile.writerow(['#VID','MID','VARFREQ','VARHAP','VARCHR','VARPOS','VARDEL','VARDELSIZE','VARINS','VARINSSEQ(HAP/SEQFILE,CHR:START-END,COPY,REVCOMP)'])
	for var in varcopy.keys():
		varcopy[var][9]=formatins(varcopy[var][9])
	for var in varcopy.keys():
		varcopy[var][5]=varcopy[var][5]+nligation #actual site
	ofile.writerows(varcopy.values())
	return None

def var2bed(bun):
	#IMPORTANT: OrderedDict and List are all mutable
	varcopy=copy.deepcopy(bun.varlist); nligation=bun.nligation; ofile=csv.writer(open(bun.oprefix+'.bed','w'),delimiter='\t')
	ofile.writerow(['#CHROM','START','END','VID','MID','VARFREQ','VARHAP','VARCHR','VARPOS','VARDEL','VARDELSIZE','VARINS','VARINSSEQ(HAP/SEQFILE,CHR:START-END,COPY,REVCOMP)'])
	for var in varcopy.keys():
		varcopy[var][9]=formatins(varcopy[var][9])
	for var in varcopy.keys():
		varcopy[var][5]=varcopy[var][5]+nligation #actual site
		ofile.writerow([varcopy[var][4],varcopy[var][5],varcopy[var][5]+varcopy[var][7]]+varcopy[var])
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
						[hapseq[dpick].filename,trunks[dtid].chrom,trunks[dtid].start+nligation,sizevec[i],copyvec[i],rc]]
					varcnt+=1
			elif metatype == 'INV' or metatype == 'IDUP':  #one deletion at donor followed by insertion of reverse
				dtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
				for d in donor:
					varid = vartag+str(varcnt)
					varlist[varid]=[varid,mid,donorfreq[i],d,trunks[dtid].chrom,trunks[dtid].start,True,sizevec[i],True,
						[hapseq[dpick].filename,trunks[dtid].chrom,trunks[dtid].start+nligation,sizevec[i],copyvec[i],True]]
					varcnt+=1
			elif metatype == 'DINS':
				dtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
				rtid,stretches=maketid(stretches,1)
				for r in recpt:    #inserting sequence: [ filename, chr, start, size, copy, rc=True/False ]
					varid = vartag+str(varcnt)
					varlist[varid]=[varid,mid,recptfreq[i],r,trunks[rtid].chrom,trunks[rtid].start,False,0,True,\
						[hapseq[dpick].filename,trunks[dtid].chrom,trunks[dtid].start+nligation,sizevec[i],copyvec[i],False]]
					varcnt+=1
			elif metatype == 'DUP':
				rtid,stretches=maketid(stretches,sizevec[i]/(plansize+2*nligation)+1)
				for r in recpt:    #inserting sequence: [ filename, chr, start, size, copy, rc=True/False ]
					varid = vartag+str(varcnt)
					varlist[varid]=[varid,mid,recptfreq[i],r,trunks[rtid].chrom,trunks[rtid].start,True,sizevec[i],True,\
						[hapseq[r].filename,trunks[rtid].chrom,trunks[rtid].start+nligation,sizevec[i],copyvec[i],False]]
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
#	for gap in gaptab:
#		print gap
#	for trunk in trunks:
#		print trunk
	remains = trunks.intersect(gaptab,v=True) #only keeping non overlapping trunks
#	for remain in remains:
#		print remain
	#FIXME: sometimes bedtool error will rise if temporary files were deleted, how to fix in safeway?
	nr = len(remains)
	cnt=1; stretches = [0] * nr; j=0;
	for feat in remains:
		if j==0:
			j = j+1; pfeat=feat; continue
		if feat.chrom == pfeat.chrom and feat.start == pfeat.end:
			cnt+=1 #continue counting from j-1 to j
		else:    #stop counting, j-1 in stretch starts from , j is not in
			stretches[j-cnt]=cnt
			for k in xrange(j-cnt+1,j): stretches[k] = 0
			cnt=1  #start new counting
		if j == nr-1 : #this is forced ending, we need to trace back (cnt-1)
			stretches[j-cnt+1]=cnt
			for k in xrange(j-cnt+2,j): stretches[k] = 0
		j=j+1; pfeat=feat
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
		pos=xrange(burnin,hapseq[0].lengths[i]-plansize-2*nligation-burnin,plansize+2*nligation)
		trunks.extend([ '%s\t%s\t%s' % (hapseq[0].references[i], str(p), str(p+plansize+2*nligation)) for p in pos ])
	trunks = pybedtools.BedTool("\n".join([ x for x in trunks ]), from_string=True) #convert to BED
	return trunks

def maketid(stretches,varstretch):
	#print stretches
	#print varstretch
	tids = [t for t, x in enumerate(stretches) if x >= varstretch] #always random pick among min qualifying tids
	if tids == []: raise Exception("cannot allocate genome region for metafile, please reduce plansize!");
	if varstretch == 1: #insertion type sequence, could from any place
		mtids = tids
	else:
		mtids = locatemin([stretches[i] for i in tids])
	#print mtids
	tid = random.sample([tids[i] for i in mtids[1]],1)[0]   #locate the tid
	#print tid
	for i in xrange(tid,tid+varstretch): stretches[i]=0;						#mark head of stretch used
	if mtids[0]-varstretch>0: stretches[tid+varstretch]=mtids[0]-varstretch; #mark remaining stretch by its head
	return tid, stretches

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
	starttime=time.time()
	run()
	print "total runtime: ", time.time()-starttime, " seconds"

#def var2vcf(bun):
#	varlist=bun.varlist; ofile=[ csv.writer(open(bun.oprefix+'.'+str('.vcf','w'),delimiter='\t')
#	vcf_meta = call_meta + sample_meta + contig_meta
#	ofile.writerow(vcf_meta)
#	ofile.writerow(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"])
#	metalist={}
#	for var in varlist.keys():
#		if not metalist.get(varlist[var][1],False):
#			metalist[varlist[var][1]].append(varlist[var])
#		else:
#			metalist[varlist[var][1]]=varlist[var]
#	for mid in metalist.keys():
#		sv_type=''.join([i for i in mid if not i.isdigit()])
#		if len(metalist[mid])
	#NOTE: only when reffile == hapfiles vcf file will make sense
	#NOTE: need to go over by MID and report correspondingly

#meta_keys=["fileformat","fileDate","source","reference","phasing"]
#meta_values=["None","None","None","None","partial"]
#info_values=[
#'<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
#'<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
#'<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
#'<ID=AA,Number=1,Type=String,Description="Ancestral Allele">',
#'<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">',
#'<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">',
#'<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
#'<ID=BKPTID,Number=.,Type=String,Description="ID of the assembled alternate allele in the assembly file">',
#'<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">',
#'<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">',
#'<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
#'<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">',
#'<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">',
#'<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
#'<ID=SOMATIC,Number=0,Type=Flag,Description="SOMATIC structural variation">',
#'<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">',
#'<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">',
#'<ID=METRANS,Number=.,Type=String,Description="Mobile element transduction info of the form CHR,START,END,POLARITY">',
#'<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
#'<ID=DGVID,Number=1,Type=String,Description="ID of this element in Database of Genomic Variation">',
#'<ID=DBVARID,Number=1,Type=String,Description="ID of this element in DBVAR">',
#'<ID=DBRIPID,Number=1,Type=String,Description="ID of this element in DBRIP">',
#'<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">',
#'<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">',
#'<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">',
#'<ID=CILEN,Number=2,Type=Integer,Description="Confidence interval around the length of the inserted material between breakends">',
#'<ID=DPADJ,Number=.,Type=Integer,Description="Read Depth of adjacency">',
#'<ID=CN,Number=1,Type=Integer,Description="Copy number of segment containing breakend">',
#'<ID=CNADJ,Number=.,Type=Integer,Description="Copy number of adjacency">',
#'<ID=CICN,Number=2,Type=Integer,Description="Confidence interval around copy number for the segment">',
#'<ID=CICNADJ,Number=.,Type=Integer,Description="Confidence interval around copy number for the adjacency">'
#]
#format_values=[
#'<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">',
#'<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">',
#'<ID=CNL,Number=.,Type=Float,Description="Copy number genotype likelihood for imprecise events">',
#'<ID=NQ,Number=1,Type=Integer,Description="Phred style probability score that the variant is novel with respect to the genome ancestor">',
#'<ID=HAP,Number=1,Type=Integer,Description="Unique haplotype identifier">',
#'<ID=AHAP,Number=1,Type=Integer,Description="Unique identifier of ancestral haplotype">',
#'<ID=GT,Number=1,Type=String,Description="Genotype">',
#'<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
#'<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
#'<ID=AD,Number=1,Type=Integer,Description="Read Supporting ALT">',#TCGA
#'<ID=BQ,Number=1,Type=Integer,Description="Base Quality ofr Read Supporting ALT">',#TCGA
#'<ID=SS,Number=1,Type=Integer,Description="Somatic Status 0=wt,1=gm,3=so,4=LOH,5=PTM/unknown">',#TCGA
#'<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">',
#'<ID=GL,Number=G,Type=Integer,Description="Genotype Likelihood">',
#'<ID=lW,Number=1,Type=Integer,Description="Max Coverage Likelihood Ratio">',
#'<ID=lWc,Number=1,Type=Integer,Description="Cumulative Max Coverage Likelihood Ratio">',
#'<ID=lCd,Number=1,Type=Integer,Description="Deletion Insert Size Likelihood Ratio">',
#'<ID=lCi,Number=1,Type=Integer,Description="Insertion Insert Size Likelihood Ratio">',
#'<ID=lDl,Number=1,Type=Integer,Description="Left/Plus Anchored Hang Read Likelihood Ratio">',
#'<ID=lDr,Number=1,Type=Integer,Description="Right/Minus Anchored Hang Read Likelihood Ratio">',
#'<ID=RGl,Number=1,Type=Integer,Description="Left Limit of Voted Region">',
#'<ID=RGr,Number=1,Type=Integer,Description="Right Limit of Voted Region">',
#'<ID=BPl,Number=1,Type=Integer,Description="Softclipped Reads Left Break Point">',
#'<ID=BPr,Number=1,Type=Integer,Description="Softclipped Reads Right Break Point">'
#]
#alt_values=[
#'<ID=DEL,Description="Deletion">',
#'<ID=DEL:ME:ALU,Description="Deletion of ALU element">',
#'<ID=DEL:ME:L1,Description="Deletion of L1 element">',
#'<ID=DUP,Description="Duplication">',
#'<ID=DUP:TANDEM,Description="Tandem Duplication">',
#'<ID=INS,Description="Insertion of novel sequence">',
#'<ID=INS:ME:ALU,Description="Insertion of ALU element">',
#'<ID=INS:ME:L1,Description="Insertion of L1 element">',
#'<ID=INV,Description="Inversion">',
#'<ID=CNV,Description="Copy number variable region">'
#]
#filter_values=[
#'<ID=q10,Description="Quality below 10">',
#'<ID=s50,Description="Less than 50% of samples have data">'
#]
#info_keys= ["INFO"] * len(info_values)
#format_keys= ["FORMAT"] * len(format_values)
#alt_keys= ["ALT"] * len(alt_values)
#filter_keys= ["FILTER"] * len(filter_values)
#sv_keys= info_keys + format_keys + alt_keys + filter_keys
#sv_values= info_values + format_values + alt_values + filter_values
#
#def sample_meta(id, genomes, mixture, desc):
#  template=['##SAMPLE=<ID=%s,Genomes=%s,Mixture=%s,Description="%s">'] * len(id)
#  return "\n".join( [ template[i] % (id[i], genomes[i], mixture[i], desc[i]) for i in range(0,len(id)) ] )
#
#def contig_meta(id, length, md5, species):
#  template=['##contig=<ID=%s,length=%s,md5=%s,species="%s">'] * len(id)
#	return "\n".join( [ template[i] % (id[i], length[i], md5[i], species[i]) for i in range(0,len(id)) ] )
#
#def call_meta(keys, values, program="mutforge"):
#  keys.append(["fileformat","fileDate","source"])
#	values.append(["VCFv4.1",time.ctime(),program)])
#  return "\n".join([ '##'+key+'='+value for key,value in izip(keys,values) ])
