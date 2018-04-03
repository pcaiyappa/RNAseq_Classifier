import pysam

#extract mismatch content for each position in each window
def extractMismatch(refBaseFileName,bamFileName,windowFileName,mismatchFileName,\
                    windowSize):
    windows = open(windowFileName)
    refBase = open(refBaseFileName)
    mismatch = open(mismatchFileName ,'w')
    alignment = pysam.AlignmentFile(bamFileName,"rb")
    for line in windows:
        [chrom,windowStart,windowEnd] = line.split()
        windowStart = int(windowStart); windowEnd = int(windowEnd)
        #pileup reads in the window
        for col in alignment.pileup(chrom, windowStart, windowEnd):
            if col.pos>=windowStart and col.pos<windowEnd:
                #get reference Base
                ref = refBase.readline().split()
                numMismatch = 0.0;
                numRead = 0.0;
                #get base in reads
                for read in col.pileups:
                    numRead += 1
                    if not read.is_del and not read.is_refskip:
                        base = read.alignment.query_sequence[read.query_position]
                        if base.upper()!= ref[-1].upper():
                            numMismatch += 1
                #bam-file is 1-based
                mismatch.write("{}\t{}\t{}\n".\
                               format(chrom,col.pos+1,numMismatch/numRead))
    windows.close()
    refBase.close()
    mismatch.close()
