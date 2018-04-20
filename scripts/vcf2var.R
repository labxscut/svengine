#!/usr/bin/env Rscript

####Externalities####
ptm=proc.time() #time keeper
library(vcfR)

option_list2 <- list(
  make_option(c("-o", "--outPrefix"), default="none",
              help="none=vcfFile.$sample.meta [default %default]"),
  make_option(c("-s", "--sampleName"), default="none",
              help="none=first sample [default %default]"),
  make_option(c("-p", "--planSize"), default=5000,
              help="plan size to avoid overlap events, default 5000"),
  make_option(c("-f", "--seqFiles"), defaul='{"REF":"human_g1k_v37.fasta","ME":"human_g1k_v37.me.fasta"}',
              help="use hg37 fasta for default [default %default]"),
  make_option(c("-r", "--reference"), default="hg37",
              help="use hg37 for default [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-u", "--isUCSC"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] vcfFile",
                       option_list=option_list2)
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
if(length(cmd$args)!=1 & !cmd$options$debug){ print_help(parser); quit() }
print(cmd$args[1])

vcfFile=cmd$args[1]
seqFiles=cmd$options$seqFiles
outPrefix=cmd$options$outPrefix
sampleName=cmd$options$sampleName
reference=cmd$options$reference
planSize=cmd$options$planSize
isUCSC=cmd$options$isUCSC
mydebug=cmd$options$debug
noQuiet=cmd$options$noQuiet

if(mydebug){ #debugging parameters
  vcfFile = "~/ws/2018_NA12878/NA12878_1KG_SV.vcf"
  seqFiles = '{"REF":"human_g1k_v37.fasta","ME":"human_g1k_v37.me.fasta"}' #a mobile element file
  outPrefix = "none"
  sampleName = "none"
  reference = "hg37"
  planSize = 5000
  isUCSC = FALSE
}

outPrefix = ifelse(outPrefix=="none", vcfFile, outPrefix)
seqSource =  rjson::fromJSON(seqFiles)
vcfEntries = vcfR::read.vcfR(vcfFile)
vcfSamples = colnames(extract.gt(vcfEntries)) #all samples
if(sampleName=="none") sampleName=vcfSamples[1] else sampleName=stringr::str_split(sampleName, ",")[[1]]
if(outPrefix=="none") outPrefix=paste(vcfFile,"meta",sep=".")
vcfMetaData = vcfR::vcfR2tidy( vcfEntries )$meta
vcfGTData = vcfR::vcfR2tidy( vcfEntries )$gt
vcfFixData = vcfR::vcfR2tidy( vcfEntries )$fix
vcfFixData = vcfFixData %>% dplyr::filter(FILTER=="PASS") 
vcfFixData = vcfFixData %>% dplyr::mutate(END=ifelse(is.na(END),POS+1,END)) #external ins has NA END, change to POS+1
#sel = ! (grepl("INS",vcfFixData$ALT) & (is.na(vcfFixData$MSTART) | is.na(vcfFixData$END)))
#vcfFixData = vcfFixData %>% dplyr::filter(sel) 
vcfRanges = GenomicRanges::GRanges(
  seqnames = vcfFixData$CHROM,
  ranges = IRanges::IRanges(start=vcfFixData$POS - planSize,
                            end=vcfFixData$END + planSize),
  id = vcfFixData$ID,
  svtype = vcfFixData$SVTYPE,
  alt = vcfFixData$ALT
)
ov = tibble::as.tibble(GenomicRanges::findOverlaps(vcfRanges, vcfRanges))
total = seq(1,length(vcfRanges))
duped = unlist(unique(ov[duplicated(ov[,1]),1]))
solo = total[!total %in% duped]
# vcfSoloRanges = vcfRanges[solo]
# GenomicRanges::findOverlaps(vcfSoloRanges,vcfSoloRanges) #verified

vcfSoloData = vcfFixData[solo, ] #cleaned up for META, total 20861
#wrote all me seq names
#write.table(na.omit(unique(sapply(stringr::str_split(vcfSoloData$MEINFO,","),"[",1))), "/Users/charlie/setup/svengine/test/na12878/me.all", row.names=F, quote=F, col.names=F)

#TEST
# i=21
# hap=unique(vcfSoloData$ALT)[i]
# entry=vcfSoloData[vcfSoloData$ALT==hap,][i,]
# hap_of_entry(hap, entry, mode)

hap_of_entry = function(hap, entry, mode) { #handling one haplotype alternative allele
  e = entry # for each hap generate a hap_allele
  if(hap == "0") {
    hap_allele = list(varchr=e$CHROM, varpos=e$POS, vardel="False", 
                      vardelsize=0, varins="False", varinsseq="None")
  } else if(grepl("CN", hap)) { # if see a <CNx> pattern
    stopifnot( mode %in% c("DEL", "DUP") ) # must be copy number events
    hap_cn = as.integer(gsub(">", "", gsub("<CN", "", hap))) #extract the copy number
    if(hap_cn == 0) { # if is a deletion event
      hap_allele = list(varchr=e$CHROM, varpos=e$POS, vardel="True", 
                        vardelsize=ifelse(is.na(e$SVLEN), e$END-e$POS+1, e$SVLEN),
                        varins="False", varinsseq="None")
    } else { # if is a duplication or no-change event
      hap_allele = list(varchr=e$CHROM, varpos=e$POS, vardel="True", 
                        vardelsize=ifelse(is.na(e$SVLEN), e$END-e$POS+1, e$SVLEN),
                        varins="True", 
                        varinsseq=paste(seqSource['REF'],",",e$CHROM,":",e$POS+1,"-",e$END,",",hap_cn,",f",sep=""))
    }
  } else if(hap %in% c("A","C","G","T")) {  #deal with A, C, G, T
    stopifnot( mode %in% c("DEL") ) # must be copy number events
    hap_allele = list(varchr=e$CHROM, varpos=e$POS, vardel="True", 
                      vardelsize=ifelse(is.na(e$SVLEN), e$END-e$POS+1, e$SVLEN),
                      varins="False", varinsseq="None")
  } else if(grepl("INS:MT", hap)) { # deal with <INS:MT>
    stopifnot( mode %in% c("INS") ) # must be copy number events
    hap_chrom = ifelse(isUCSC, "ChrM", "MT")
    hap_allele = list(varchr=e$CHROM, varpos=e$POS, vardel="False", vardelsize=0, varins="True", 
                      varinsseq=paste(seqSource['REF'],",",hap_chrom,":",e$MSTART+1,"-",e$MEND,",",1,",f",sep="")) #why there is no strand info
  } else if(grepl("INS", hap)) { # deal with <INS:ME:LINE1>, <INS:ME:ALU>, <INS:ME:SVA>
    stopifnot( mode %in% c("INS") ) # must be copy number events
    hap_type = dplyr::last(stringr::str_split( gsub(">|<", "", hap), ":", simplify = TRUE))
    hap_source = seqSource['ME']
    hap_info = stringr::str_split(e$MEINFO, ",", simplify = T)
    hap_chrom = hap_info[1]
    hap_chrom = ifelse(hap_chrom=="LINE1", "L1", hap_chrom) #Repbase use L1 instead of LINE1
    hap_chrom = ifelse(hap_chrom=="AluUndef", "AluY", hap_chrom) #Repbase does not have AluUndef, use AluY
    hap_start = as.integer(hap_info[2])
    hap_end = as.integer(hap_info[3])
    hap_strand = ifelse(hap_info[4] == "+", "f", "r")
    hap_allele = list(varchr=e$CHROM, varpos=e$POS, vardel="False", vardelsize=0, varins="True", 
                      varinsseq=paste(hap_source,",",hap_chrom,":",hap_start,"-",hap_end,",",1,",f",sep="")) 
  } else if(grepl("INV", hap)) { # deal with <INV>
    stopifnot( mode %in% c("INV") ) # must be copy number events
    hap_allele = list(varchr=e$CHROM, varpos=e$POS, vardel="True", 
                      vardelsize=ifelse(is.na(e$SVLEN), e$END-e$POS+1, e$SVLEN),
                      varins="True", 
                      varinsseq=paste(seqSource['REF'],",",e$CHROM,":",e$POS,"-",e$END,",",1,",r",sep=""))
  } else {
    stop("unhandled hap")
  }
  return(hap_allele)
}
  
vcf_to_meta = function(i, vcf) { # vcf =  tmpSoloData; i = 1
  e = vcf[i,]
  alt = c("0", strsplit(e$ALT, ",")[[1]])
  gt = ifelse(e$gt_GT=="0|0", "1|1", e$gt_GT)
  gt = as.integer(strsplit(gt, "\\|")[[1]]) + 1 # 0|0 is homozygous
  stopifnot(length(gt)==2) # two haplotypes
  stopifnot(max(gt)<=length(alt)) # two haplotypes
  #print(alt)
  #print(gt)
  if(e$SVTYPE %in% c("DEL", "DEL_LINE1", "DEL_SVA", "DEL_HERV", "DEL_ALU")) {
    m = "DEL"
    alt0 = hap_of_entry(alt[gt[1]], e, m)
    alt1 = hap_of_entry(alt[gt[2]], e, m)
  } else if(e$SVTYPE %in% c("INS", "LINE1", "ALU", "SVA")) {
    m = "INS"
    alt0 = hap_of_entry(alt[gt[1]], e, m)
    alt1 = hap_of_entry(alt[gt[2]], e, m)
  } else if(e$SVTYPE %in% c("CNV", "DUP")) {
    m = "DUP"
    alt0 = hap_of_entry(alt[gt[1]], e, m)
    alt1 = hap_of_entry(alt[gt[2]], e, m)
  } else if(e$SVTYPE %in% c("INV")) {
    m = "INV"
    alt0 = hap_of_entry(alt[gt[1]], e, m)
    alt1 = hap_of_entry(alt[gt[2]], e, m)
  } else {
    stop("unrecognized SVTYPE")
  }
  alt0$mid = paste(e$SVTYPE, i, sep="_")
  alt1$mid = alt0$mid
  alt0$varhap = 0 
  alt1$varhap = 1
  alt0$varfreq = 1.0
  alt1$varfreq = 1.0
  res = rbind(alt0, alt1)
  row.names(res) = NULL
  return(res)
}

for(vs in sampleName) { #vs="NA12878"
  tmpGTData = vcfGTData %>% dplyr::filter(Indiv == vs)
  tmpSoloData = vcfSoloData %>% dplyr::left_join(vcfGTData, by=c("ChromKey", "POS"))
  vcfMetaData = lapply(seq.int(nrow(tmpSoloData)), vcf_to_meta, vcf=tmpSoloData)
  vcfMetaTable = as.data.frame(do.call("rbind", vcfMetaData))
  vcfMetaTable$vid = seq.int(nrow(vcfMetaTable))
  colnames(vcfMetaTable) = toupper(colnames(vcfMetaTable))
  vcfMetaTable = vcfMetaTable[, c("VID","MID", "VARFREQ", "VARHAP",	"VARCHR",	"VARPOS",	"VARDEL",	"VARDELSIZE",	"VARINS",	"VARINSSEQ")] #order it correctly
  vcfMetaTable$VID = paste("VAR",vcfMetaTable$VID,sep="_")
  vcfMetaTable = apply(vcfMetaTable,2,as.character)
  write.table(vcfMetaTable, file=paste(outPrefix, vs, "meta", sep="."), col.names=T, quote=F, sep="\t", row.names=F)
}

message("VCF2VARDONE")

# "DUP" "CNV"
# "DEL"       "DEL_ALU"   "INV"       "INS"       "DEL_LINE1"
# "DEL_SVA"   "DEL_HERV"  "LINE1"     "ALU"       "SVA"

# ALT of SVTYPE
# "DUP" =  <CN0>,<CN2>       <CN2> <CN2>,<CN3> 
#                  309        5713           3       
# "DEL" = <CN0> <CN0>,<CN2>           A           C           G           T 
# 33105          23        2262        1690        1745        2150 
# "DEL_ALU" = <CN0>     A     C     G     T 
#2   199   133   118   786
# "INV" = <INV> 
#786 
# "INS" = <INS:MT> 
#168
# "DEL_LINE1" = <CN0>     A     C     G     T 
#15     8     6     1    26
# "DEL_SVA" = <CN0>     C     T 
# 1     2     6 
# "DEL_HERV" = <CN0> 
# 1 
# "LINE1" = <INS:ME:LINE1> 
#150
# "ALU" = <INS:ME:ALU> 
# 285
# "SVA" = <INS:ME:SVA> 
#16
# "CNV" = <CN0>,<CN2> 
# 2636 
# <CN0>,<CN2>,<CN3> 
#   92 
# <CN0>,<CN2>,<CN3>,<CN4> 
#   34 
# <CN0>,<CN2>,<CN3>,<CN4>,<CN5> 
#   23 
# <CN0>,<CN2>,<CN3>,<CN4>,<CN5>,<CN6> 
#   9 
# <CN0>,<CN2>,<CN3>,<CN4>,<CN5>,<CN6>,<CN7> 
#   3 
# <CN2>,<CN3> 
#   98 
# <CN2>,<CN3>,<CN4> 
#   21 
# <CN2>,<CN3>,<CN4>,<CN5> 
#   7 
# <CN2>,<CN3>,<CN4>,<CN5>,<CN6> 
#   1 
# <CN2>,<CN3>,<CN4>,<CN5>,<CN6>,<CN7> 
#   3 
# <CN2>,<CN3>,<CN4>,<CN5>,<CN6>,<CN7>,<CN8> 
#   1 
# <CN2>,<CN3>,<CN4>,<CN5>,<CN6>,<CN7>,<CN8>,<CN9> 
#   1

# DEL, INS, INV, DUP

#Fix, remove overlapping entries
#Fix, more than two allele => two allele
#Fix, DEL/DEL_* ALT = A, C, G, T => SVLEN=LEN(REF)-1; use CHROM, POS, END
#Fix, DUP/CNV
#Fix, 


#what's the plan size, say 5k
#discard anything with overlapping 5k
