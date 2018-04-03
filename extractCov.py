import pysam

#extract coverage and for each position in all regions specified;
#create sliding windows with non-zero coverage at each position for all regions;

def extractCov(regionFileName,bamFileName,covFileName,windowSize,stepSize,\
               windowFileName):
    regions = open(regionFileName)
    alignment = pysam.AlignmentFile(bamFileName,"rb")
    windows = open(windowFileName,'w')
    coverage = open(covFileName,'w')
    for line in regions:
        region = line.split()
        extractCovPerRegion(region,alignment,coverage,windowSize,stepSize,\
                            windows)
    regions.close()
    windows.close()
    coverage.close()

#extract coverage and for each position in one region;
#create sliding windows with non-zero coverage at each position in this region;
def extractCovPerRegion(region,alignment,coverage,windowSize,stepSize,\
               windows):
    [chrom,regionStart,regionEnd] = region
    regionStart = int(regionStart); regionEnd = int(regionEnd)
    cov = {} #stores the coverage and position in the region for each coordinate
    posInPileup = 0; #order of non-zero positions
    for col in alignment.pileup(chrom,regionStart,regionEnd):
        cov[col.pos] = [col.n,posInPileup]
        posInPileup += 1
    # create sliding windows in this region;
    windowStart = regionStart; windowEnd = regionStart+windowSize
    while windowEnd <= regionEnd:
        # skip windows containing missing (zero coverage) positions
        if (windowStart in cov) and (windowEnd in cov) and \
           ((cov[windowEnd][1]-cov[windowStart][1]) == windowSize):
            windows.write("{}\t{}\t{}\n".format(chrom,windowStart,windowEnd))
            for pos in range(windowStart,windowEnd):
                #samFile is 1-based
                coverage.write("{}\t{}\t{}\n".format(chrom,pos+1,cov[pos][0]))          
        windowStart += stepSize; windowEnd = windowStart+windowSize
