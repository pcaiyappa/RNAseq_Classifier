import subprocess
from sys import argv
import shlex
from re import sub
import os
import pysam

num_threads = 1

def Cutadapt_Routine(fastq_file):
    print ("\n1.Initiating Quality Control Routine with Cutadapt 1.9.1\n")
    outfile = sub('[.].*','',fastq_file) + '_trimmed.fastq'
    command = 'cutadapt -m 30 -q 20,20 -e 0.06 -a adaptors.fa ' + fastq_file + ' -o ' + outfile
    # -m controls minimum read length at 30bp, -q trims bases from both ends with PHRED score <20 , -e controls maximum mismatch rate while aligning the adaptor at 6%
    arguments = shlex.split(command)
    with open('CutAdapt_log.txt', "w+") as file:
        file.write('**************************Cutadapt Standard Output******************************************\n')
        cutadapt_process = subprocess.Popen(arguments,stdout = file)
        cutadapt_process.wait()
        cutadapt_error_code = str(cutadapt_process.returncode)
        file.write('\nCutadapt terminated with exit code ' + cutadapt_error_code)
        print("\nCutadapt terminated with exit code " + cutadapt_error_code + ". Please check the file 'CutAdapt_log.txt' for statistics.")
        print("\nThe trimmed fastq file is in the local folder.")
        print("*******************************************************************************************")
        return (outfile)

def FastQC_Routine(input_file) :
    print ("\n2.Generating FastQC report. Please check the \"FastQC_report\" folder for FastQC Log.\n\n")
    current_directory = os.getcwd()
    fastqc_report_directory = current_directory + '/' + 'FastQC_report'
    if not os.path.exists(fastqc_report_directory):
        os.makedirs(fastqc_report_directory)

    os.chdir(fastqc_report_directory)
    command = current_directory + "/FastQC/fastqc " + current_directory + '/' + input_file
    subprocess.call( command + ' 2> fastqc_log.txt', shell=True)
    print("\nFastQC has finished running. Please check local folder for the html report.")
    print("*******************************************************************************************")
    os.chdir(current_directory)

def Index_Reference_Genome(genome_file):
        print ("\n3.Indexing the Reference Genome as required by the Bowtie2 Aligner.\n\n")

        reference_index_name = sub('^.*/','',genome_file)
        reference_index_name = sub('.fa.*','',reference_index_name)

        current_directory = os.getcwd()
        command = 'bowtie2-build ' + current_directory + '/' + genome_file + ' ' + reference_index_name
        arguments = shlex.split(command)

        reference_index_directory = current_directory + '/' + reference_index_name + '_bowtie2_index'
        if not os.path.exists(reference_index_directory):
                os.makedirs(reference_index_directory)
        os.chdir(reference_index_directory)

        with open('Bowtie2_Build_log.txt','w') as file:
            file.write('**************************Bowtie2-Indexer Standard Output******************************************\n')
            bowtie_indexer_process = subprocess.Popen(arguments,stdout=file)
            bowtie_indexer_process.wait()
            bowtie_indexer_error_code = str(bowtie_indexer_process.returncode)
            file.write('\nBowtie2_Build terminated with exit code ' + bowtie_indexer_error_code)

        os.chdir(current_directory)
        reference_index = reference_index_directory + '/' + reference_index_name
        print "\nReference Genome Indexed in \'foo_bowtie_index\' folder."
        print("*******************************************************************************************")
        return(reference_index)

def Bowtie_aligner(index_file,fastq_file):
        print ("\n4.Aligning the reads using Bowtie2.\n\n")
        num_threads = input("Please enter the number of cpu threads to be used by Bowtie\n")
        print("\n")
        current_directory = os.getcwd()
        threads_argument = '-p ' + str(num_threads)
        command = current_directory + '/bowtie2/' + 'bowtie2 -x ' + index_file + ' ' + threads_argument + ' -U ' + current_directory + '/' + fastq_file + ' -S ' + 'alignment.sam'
        arguments = shlex.split(command)

        alignment_directory = current_directory + '/' + 'bowtie2_alignment'
        if not os.path.exists(alignment_directory):
            os.makedirs(alignment_directory)
        os.chdir(alignment_directory)

        with open('Bowtie2_alignment_log.txt', 'w+') as file :
            file.write('**************************Bowtie2-Aligner Standard Output******************************************\n')
            bowtie_aligner_process = subprocess.Popen(arguments,stdout = file)
            bowtie_aligner_process.wait()
            bowtie_aligner_return_code = str(bowtie_aligner_process.returncode)
            file.write("\nBowtie2 terminated with exit code " + bowtie_aligner_return_code)

        print "\nBowtie finished Alignment. Sam file created in \'Bowtie2_alignment\' folder."
        print("*******************************************************************************************")

        return alignment_directory,current_directory


def Convert_to_bam(parent_directory) :
    print ("\n5.Converting sam file to bam file using samtools\n\n")
    samtools_location = parent_directory + '/samtools/bin/samtools'
    threads_argument = '-@ ' + str(num_threads)
    with open(os.devnull, 'w') as devnull:
        subprocess.call(samtools_location + ' view ' + threads_argument + ' -q 1 -b alignment.sam -o alignment.bam',shell = True,stdout=devnull,stderr=devnull)
    os.remove('alignment.sam')
    print "\n SAM To BAM conversion complete."
    print("*******************************************************************************************")


def Sort_Index_bam(parent_directory,outfile):
    print ("6.Sorting and Indexing the bam file\n\n")
    samtools_location = parent_directory + '/samtools/bin/samtools'
    threads_argument = '-@ ' + str(num_threads)
    with open(os.devnull, 'w') as devnull:
        subprocess.call(samtools_location + ' sort ' + threads_argument + ' alignment.bam ' + outfile, shell=True,stdout=devnull,stderr=devnull)
        subprocess.call(samtools_location + ' index ' + outfile + '.bam', shell=True,stdout=devnull,stderr=devnull)
        os.remove('alignment.bam')
    print "\nFinal bam file is sorted. Index .bai file created. Both can be found in the \'bowtei2_alignment\' folder."
    return(0)
