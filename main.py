import os
from sys import argv
from extractFeature import extractFeature
import aligner as align
from testing import MODEL_TESTING

global num_threads

def main():

    if len(argv) == 1:
        print ("instructions")
    else:
        command = argv[1]
        if command == "find":
            current_directory = os.getcwd()
            temp_directory = current_directory + '/temp'
            if not os.path.exists(temp_directory):
                os.makedirs(temp_directory)
            bamFileName = argv[2]
            genomeFileName = argv[3]
            regionFileName = argv[4]
            outputFileName = argv[5]

            featFileList = extractFeature(regionFileName,bamFileName,genomeFileName)
            coverage_file=featFileList[0]
            base_file = featFileList[1]
            gc_file = featFileList[2]
            mismatch_file=featFileList[3]
            mfe_file=featFileList[4]
            output_file=outputFileName
            MODEL_TESTING(coverage_file, base_file, gc_file, mismatch_file, mfe_file,output_file)

        elif command == "align":
            trimmed_fastq_filename = align.Cutadapt_Routine(argv[2])
            align.FastQC_Routine(trimmed_fastq_filename)
            reference_index_file = align.Index_Reference_Genome(argv[3])
            alignment_directory,parent_directory = align.Bowtie_aligner(reference_index_file,trimmed_fastq_filename)
            os.chdir(alignment_directory)
            align.Convert_to_bam(parent_directory)
            align.Sort_Index_bam(parent_directory,argv[4])
            print "Alignment Complete. Please use the .bam file in the bowtie2_alignment folder to run the \'find\' method."
        else:
            print ("Please reconsider Arguments!")

if __name__ == '__main__':
    main()