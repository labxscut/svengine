#!/usr/bin/env python
#TODO: add FINSERIONS and DINSERTIONS

import copy, time, re, os, shutil, sys, random, csv, subprocess, traceback, argparse, tempfile, itertools, operator, collections, ConfigParser, json, imp, StringIO
from multiprocessing import Pool
import pysam
import pybedtools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import pybedtools.featurefuncs
#properly import muforge
try:
	import mf.mutforge as mf #in debugging
except ImportError:
	import mutforge as mf #in production

def lumpy2var(bp,ref,nligation):
	idx=1; delidx=1; insidx=1; traidx=1; invidx=1; dupidx=1; varlist=collections.OrderedDict()
	for line in bp:
		if line[10] == "TYPE:DELETION":
			mid="DEL_"+str(delidx); var="VAR_"+str(idx); idx+=1; delidx+=1
			varlist[var]=[None] * 10
			varlist[var][0] = var	#VID
			varlist[var][1] = mid	#MID
			varlist[var][2] = 1.0 #VAR FREQ
			varlist[var][3] = 0		#VAR HAP
			varlist[var][4] = line[0]		#VAR CHR
			varlist[var][5] = int(line[1])-nligation #VAR POS
			varlist[var][6] = True	#VAR DEL
			varlist[var][7] = int(line[4])-int(line[1])+1 #VAR DEL SIZE
			varlist[var][8] = False	#VAR INS
			varlist[var][9] = None	#VAR INS SIZE
		if line[10] == "TYPE:INVERSION":
			mid="INV_"+str(invidx); var="VAR_"+str(idx); idx+=1; invidx+=1
			varlist[var]=[None] * 10
			varlist[var][0] = var	#VID
			varlist[var][1] = mid	#MID
			varlist[var][2] = 1.0 #VAR FREQ
			varlist[var][3] = 0		#VAR HAP
			varlist[var][4] = line[0]		#VAR CHR
			varlist[var][5] = int(line[1])-nligation	#VAR POS
			varlist[var][6] = True	#VAR DEL
			varlist[var][7] = int(line[4])-int(line[1])+1 #VAR DEL SIZE
			varlist[var][8] = True	#VAR INS
			varlist[var][9] = (ref.name,	line[0], int(line[1]), int(line[4])-int(line[1])+1, 1, True)	#VAR INS SIZE
		if line[10] == "TYPE:DUPLICATION": #something weird with lumpy's output here, start is bigger than end
			mid="DUP_"+str(dupidx); var="VAR_"+str(idx); idx+=1; dupidx+=1
			varlist[var]=[None] * 10
			varlist[var][0] = var	#VID
			varlist[var][1] = mid	#MID
			varlist[var][2] = 1.0 #VAR FREQ
			varlist[var][3] = 0		#VAR HAP
			varlist[var][4] = line[0]		#VAR CHR
			varlist[var][5] = int(line[4])-nligation	#VAR POS
			varlist[var][6] = True	#VAR DEL
			varlist[var][7] = int(line[1])-int(line[4])+1 #VAR DEL SIZE
			varlist[var][8] = True	#VAR INS
			varlist[var][9] = (ref.name,	line[0], int(line[4]), int(line[1])-int(line[4])+1, 2, False)	#VAR INS SIZE
		if line[10] == "TYPE:TRANSLOCATION":
			line2 = bp.next() #this can be ingnored
			if int(line[2])-int(line[1])==1:#where insertion happens
				delst=int(line[4]); deled=int(line[5]); delchr=line[3]; insst=int(line[1]); insted=int(line[2]); inschr=line[0]
			else:
				delst=int(line[1]); deled=int(line[2]); delchr=line[0]; insst=int(line[4]); insted=int(line[5]); inschr=line[3]
			#do insertion
			mid="TRA_"+str(traidx); var="VAR_"+str(idx); idx+=1; 
			varlist[var]=[None] * 10
			varlist[var][0] = var	#VID
			varlist[var][1] = mid	#MID
			varlist[var][2] = 1.0 #VAR FREQ
			varlist[var][3] = 0		#VAR HAP
			varlist[var][4] = inschr		#VAR CHR
			varlist[var][5] = insst-nligation	#VAR POS
			varlist[var][6] = False	#VAR DEL
			varlist[var][7] = 0 #VAR DEL SIZE
			varlist[var][8] = True	#VAR INS
			varlist[var][9] = (ref.name,	delchr, delst, deled-delst+1, 1, False)	#VAR INS SIZE
			mid="TRA_"+str(traidx); var="VAR_"+str(idx); idx+=1; traidx+=1
			#do deletion
			varlist[var]=[None] * 10
			varlist[var][0] = var	#VID
			varlist[var][1] = mid	#MID
			varlist[var][2] = 1.0 #VAR FREQ
			varlist[var][3] = 0		#VAR HAP
			varlist[var][4] = delchr		#VAR CHR
			varlist[var][5] = delst-nligation	#VAR POS
			varlist[var][6] = True	#VAR DEL
			varlist[var][7] = deled-delst+1 #VAR DEL SIZE
			varlist[var][8] = False	#VAR INS
			varlist[var][9] = None	#VAR INS SIZE
	return varlist

def var2seq(brpt,reffile,nligation,oprefix): #assume varbed is checked for nonoverlapping
	segst=0; ref=pysam.Fastafile(reffile.name); seq=[]
	for i in xrange(0,ref.nreferences):
		rname=ref.references[i]; rlen=ref.lengths[i]
		rseq=Seq(ref.fetch(rname,0,rlen), generic_dna)
		if brpt.get(rname,None) == None: #non var at all
			seg = [ rseq ]
		else:
			brs = collections.OrderedDict(sorted(brpt[rname].items(), key=lambda t: t[0]))
			seg = []; nbrs=len(brs) #member of seg are Seqs
			for j in xrange(0,nbrs):
				#VID MID VARFRQ VARHAP VARCHR VARPOS VARDEL VARDELSZ VARINS VARINSSZ
				var=brs.values()[j]
				if j == (nbrs-1): #last var
					seged=rlen
				else:
					seged=brs.values()[j+1][5]
				if j == 0: #first seg appended without var
					seg.append(rseq[segst:var[5]]); segst=var[5] 
				seq1 = rseq[segst:seged] #seg if connected is the new seq 
				if var[6]: #work on deletion
					seq1 = seq1[:nligation]+seq1[nligation+var[7]:]
					segst = seged
				if var[8]: #work on insertion
					#['example.fasta', '22', 25202000, 1000, 2, False]
					info = var[9]
					seqi = Seq(pysam.Fastafile(info[0]).fetch(info[1],info[2],info[2]+info[3]), generic_dna) #st=0-based, ed=0-based+1
					if info[5]:
						seqi = seqi.reverse_complement();
					seqi = Seq(str(seqi)*info[4], generic_dna)
					seq1 = seq1[:nligation]+seqi+seq1[nligation:]
					segst = seged
				seg.append(seq1)
		#print "old:",len(rseq)
		#print "new:",len(SeqRecord( Seq(''.join([str(s) for s in seg]), generic_dna), id=rname, name=rname, description=oprefix+"."+rname ) )
		seq.append( SeqRecord( Seq(''.join([str(s) for s in seg]), generic_dna), id=rname, name=rname, description=oprefix+"."+rname ) )
	SeqIO.write(seq,oprefix+".fnv",'fasta')
	return seq

def main(args):
	nprocs = args.nprocs
	varfile = args.varfile
	metafile = args.metafile
	reffile = args.reffile
	gapfile = args.gapfile
	parfile = args.parfile
	parlist = ConfigParser.ConfigParser()
	parlist.readfp(parfile)
	plansize = args.plansize
	edgesize = 100000
	layout = args.layout
	runmode=0 #1 meta, 2 var, 3 var+meta
	if varfile!=None and metafile!=None:
		runmode=3;
		oname=varfile.name.rsplit(".")[0]	
	elif varfile!=None and metafile==None:
		runmode=2;
		oname=varfile.name.rsplit(".")[0]	
	elif metafile!=None and varfile==None:
		runmode=1;
		oname=metafile.name.rsplit(".")[0]	
	else:
		raise Exception("must provide meta and/or var files!")
	oprefix = oname if not args.oprefix else args.oprefix
	nligation = 0; nploid = 1; 
	freq1 = [ float(f) for f in args.mixfreq.split(",") ]
	for f in freq1:
		assert f>=0 and f<=1, "mixfreq must be a float between 0 and 1" 
	try:
		sizefile = mf.fas2size(reffile.name) #create genome size file
		gaptab=pybedtools.BedTool(gapfile.read(),from_string=True)
		gaptab=gaptab.slop(g=sizefile, b=edgesize)
	except pybedtools.cbedtools.BedToolsFileError:
		raise Exception("incorrect gapfile provided, must in BED format")
	#create a string meta file:
	#FINS	2500	adeno.fna,fix_...	fix_1	fix_1	fix_1	fix_1	fix_1	#FINS=COPY*INS, foreign insertion
	#DINS	2500	fix_...	fix_1	fix_1	fix_1	fix_1	fix_1	#DINS=COPY*INS, domestic insertion
	#inssize="fix_100_200_300_400_500_600_700_800_900_1000_1000_2000_3000_4000_5000_6000_7000_8000_9000_10000"
	#metastr='\n'.join([ "\t".join(['FINS',str(2),'adeno.fna,'+inssize,'fix_1','fix_1','fix_1','fix_1','fix_1' ]),
	#					"\t".join(['DINS',str(2),'fix_100_200','fix_1','fix_1','fix_1','fix_1','fix_1' ]) ])
	varlist = collections.OrderedDict()
	dic={'nligation':nligation,'varlist':varlist,'oprefix':oprefix,'reffile':reffile,'nploid':nploid,'plansize':plansize,'varcnt':0,'edgesize':edgesize,'gaptab':gaptab}
	print "runmode=", runmode

	if runmode > 1: #have var input
		varintype = varfile.name.split(".")[-1]
		if varintype == 'lumpy':
			bp=mf.CommentedFile(varfile)
			#bp=csv.reader(varfile,delimiter="\t")
			#somehow we need to set ligation setting 
			#assume that every true var pos out of program is true position
			#every true var pos inside the program is position+nligation
			varlist=lumpy2var(bp,reffile,nligation)
		else:
			raise Exception("varintype %s hasn't been implemented" % varintype)
		dic['varcnt']=len(varlist) #this is needed for continuously counting
		dic['varlist']=varlist #this is needed for continuously counting
		#d.iloc[0:1,0:1]; deoverlapping
		vartab=mf.var2bed(mf.Bunch(dic))
		ovlptab=vartab.intersect(vartab,c=True).to_dataframe(names=mf.varbed_names+["ovlap"])
		print "input vars before deovlapping", len(ovlptab)
		assert len(ovlptab) == len(varlist), "assert the same set of vars in ovlptab and varlist"
		ovlpkey=[ varlist.keys()[j] for j in range(0,len(ovlptab)) if int(ovlptab.iloc[j,len(mf.varbed_names)]) > 1 ]
		for key in ovlpkey:
			del varlist[key] #exclude these overlapping
		dic['varlist']=varlist #this is needed for vartab
		vartab=mf.var2bed(mf.Bunch(dic))
		print "input vars after deovlapping", len(vartab)
		assert mf.chkvar(vartab) == True, "==Warn: varfile contains self overlaps, which is not allowed"
		gaptab=vartab.cat(gaptab.each(pybedtools.featurefuncs.extend_fields,n=vartab.field_count()).saveas())
		dic['gaptab']=gaptab

	if runmode%2 == 1:	
		dic['hapseq']=[pysam.Fastafile(reffile.name)];
		metatab = mf.CommentedFile(metafile)
		varlist,varcnt = mf.makevar(mf.Bunch(dic),metatab)
		dic['varlist'] = varlist
		dic['varcnt'] = varcnt
		vartab = mf.var2bed(mf.Bunch(dic))
		if len(vartab)>0:
			print "all output vars", len(varlist)
			ovlptab=vartab.intersect(vartab,c=True).to_dataframe(names=mf.varbed_names+["ovlap"])
			print "input vars before deovlapping", len(ovlptab)
			assert len(ovlptab) == len(varlist), "assert the same set of vars in ovlptab and varlist"
			ovlpkey=[ varlist.keys()[j] for j in range(0,len(ovlptab)) if int(ovlptab.iloc[j,len(mf.varbed_names)]) > 1 ]
			for key in ovlpkey:
				del varlist[key] #exclude these overlapping
			dic['varlist']=varlist #this is needed for vartab
			vartab=mf.var2bed(mf.Bunch(dic))
			print "input vars after deovlapping", len(vartab)
			assert mf.chkvar(vartab) == True, "==Warn: varfile still contains self overlaps, which is not allowed"

	mf.var2file(mf.Bunch(dic))
	#read in bedpe to generate varlist
	#now we need to check whether the  sv has been correctly inserted
	#we also need to check whether the position is properly selected 
	brpt=collections.OrderedDict()
	for key in varlist.keys():
		var = varlist[key]
		if brpt.get(var[4],None) == None:
			brpt[var[4]] = collections.OrderedDict({var[5]:var})
		else:
			brpt[var[4]][var[5]]=var  #brpt[chr][pos]
	ref = pysam.Fastafile(reffile.name)
	seq0 = [ SeqRecord(Seq(ref.fetch(ref.references[i],0,ref.lengths[i]) ),id=ref.references[i],name=ref.references[i],description=ref.references[i]) for i in xrange(0,ref.nreferences) ]
	seq1 = var2seq(brpt,reffile,nligation,oprefix)
	nlibrary=parlist.getint('xwgsim', 'nlibrary')
	bams = [None] * nlibrary

	if layout:
		print "done layout and fas file generation"
		quit()
	
	pariter=[]
	for lib in range(0,nlibrary):
		for fi in range(0,len(freq1)):
			pariter.append([parlist,nligation,reffile.name,nploid,seq0,lib,freq1[fi],seq0,seq1,lib,oprefix])

	pool = Pool(processes=int(nprocs))
	bedbams=[]
	bedbams=pool.map(makebam,pariter,1)
	
	return None

def makebam(parit):		
	print "processing lib %s freq %s regions" % (str(parit[5]),str(parit[6]))
	parlist=parit[0]; nligation=parit[1]; rname=parit[2]; nploid=parit[3]; seq0=parit[4]; lib=parit[5]; fq1=parit[6]; fq0=1-fq1; seq0=parit[7]; seq1=parit[8]; lib=parit[9]; oprefix=parit[10]
	libnames=json.loads(parlist.get('xwgsim', 'libnames'))
	bam0 = mf.runwgs(parlist,nligation,rname,nploid,seq0,lib,fq0)
	bam1 = mf.runwgs(parlist,nligation,rname,nploid,seq1,lib,fq1)
	#print [bam0,bam1]
	bam = mf.mergebam([bam0,bam1],oprefix+"."+libnames[lib])
	print "done lib %s freq %s regions" % (str(parit[5]),str(parit[6]))
	return bam

def run():
	parser = argparse.ArgumentParser(description="spike-in mutations to fasta files, \
                                                  we insert mutations into a provided reference \
                                                  this script is haplotype insensitive and ligation less")
	parser.add_argument('gapfile', metavar='gapfile', type=argparse.FileType('rU'), #example.bed
                        help="inputfile of gap or other banned regions, required")
	parser.add_argument('parfile', metavar='parfile', type=argparse.FileType('rU'), #example.par
                        help="inputfile of library-wise parameters, required"),
	parser.add_argument('reffile', metavar='reffile', type=argparse.FileType('rU'), #example.fasta
                        help="inputfile of ref file(s), required"),
	parser.add_argument('-m', '--metafile', dest='metafile', type=argparse.FileType('rU'), #example.meta
                        help="inputfile of meta specification")  # variants to be randomly added
	parser.add_argument('-v', '--varfile', dest='varfile', type=argparse.FileType('rU'),  #example.var/op
                        help="inputfile of variant specification") # prespecified variants, bedpe format or others
	parser.add_argument('-o', '--oprefix', dest='oprefix', default=None, #output1,output2,output3,...
                        help="prefix for output files")
	parser.add_argument('-f', '--plansize', dest='plansize', type=int, default=100,
                        help="palnsize is the grain size for genome planning, nligation=0 here (default=100bp)")
	parser.add_argument('-q', '--mixfreq', dest='mixfreq', default="1",
                        help="generated spikein-in mixing frequency, homozygous=1")
	parser.add_argument('-n', '--nprocs', dest='nprocs', default=1,
                        help="split into multiple processes (default=1)")
	parser.add_argument('-l', '--layout', dest='layout', action='store_true', default=False,
                        help="only dry layout (default=False)")
	args = parser.parse_args()
	main(args)

if __name__ == '__main__':
	starttime=time.time()
	run()
	print "total runtime: ", time.time()-starttime, " seconds"
