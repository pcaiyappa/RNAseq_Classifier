from extractCov import extractCov
from extractRef import extractRef
from extractMismatch import extractMismatch
from mfe import Compute_MFE

#input files
##regionFileName = "../window50/negativeSample50.txt"
##bamFileName = "replicate_b/bowtie_alignment_sorted.bam"
##genomeFileName = "hg19.fa"

def extractFeature(regionFileName,bamFileName,genomeFileName):
    
    #intermediate files
    windowFileName = "temp/windows.bed"
    covFileName = "temp/coverage.txt"
    refBaseFileName = "temp/refBase.txt"
    refSeqFileName = "temp/refSeq.txt"
    gcFileName = "temp/gcContent.txt"
    mismatchFileName = "temp/mismatch.txt"
    mfeFileName = "temp/mfe.txt"

    #list of files containing the features for classfication
    fileNameList = [covFileName,refBaseFileName,gcFileName,mismatchFileName,\
                    mfeFileName]
    #parameters
    windowSize = 40
    stepSize = 10

    #extract features
    extractCov(regionFileName,bamFileName,covFileName,windowSize,stepSize,\
              windowFileName)
    extractRef(windowFileName,genomeFileName,refBaseFileName,refSeqFileName,\
               gcFileName,windowSize)
    extractMismatch(refBaseFileName,bamFileName,windowFileName,mismatchFileName,\
                        windowSize)
    Compute_MFE(refSeqFileName,mfeFileName)
    return fileNameList
