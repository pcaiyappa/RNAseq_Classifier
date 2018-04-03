import pysam
#extracts reference base for each position
#extracts reference seq for each window
#extracts gc content for each window
def extractRef(windowFileName,genomeFileName,refBaseFileName,refSeqFileName,\
               gcFileName,windowSize):
    windows = open(windowFileName)
    refBase = open(refBaseFileName,'w')
    refSeq = open(refSeqFileName,'w')
    gcContent = open(gcFileName,'w')
    GC = {"G","g","C","c"}; #symbols for G or C
    with pysam.FastaFile(genomeFileName) as fh:
        for line in windows:
            [chrom,windowStart,windowEnd] = line.split();
            windowStart = int(windowStart); windowEnd = int(windowEnd)
            seq = fh.fetch(chrom,windowStart,windowEnd);
            #write per position ref base
            for i in range(0,windowSize):
                pos = windowStart+i+1; #samfile is 1-based
                refBase.write("{}\t{}\t{}\n".format(chrom,pos,seq[i].upper()))
            #write per window ref seq
            refSeq.write("{}\t{}\t{}\t{}\n".\
                         format(chrom,windowStart,windowEnd,seq.upper()))
            #calculate gc
            numGC = 0.0
            for base in seq:
                if base in GC: numGC += 1
            gcContent.write("{}\t{}\t{}\t{}\n".\
                            format(chrom,windowStart,windowEnd,\
                                   numGC/windowSize))         
    windows.close()
    refBase.close()
    refSeq.close()
    gcContent.close()
