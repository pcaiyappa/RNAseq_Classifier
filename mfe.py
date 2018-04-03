import subprocess
import os
import shlex
from re import sub

def Compute_MFE(reference_sequence_file,out_file) :
    with open(reference_sequence_file) as file:
        with open('sequences_RNAFold.txt','w+') as sequences_file:
            count = 0
            for line in file:
                columns = line.split('\t')
                fasta_header = '>' + columns[0] + ':' +  columns[1] + ':' + columns[2]
                sequence = columns[3]
                if not 'N' in sequence:
                    sequences_file.write("%s\n%s"%(fasta_header,sequence))
                    count += 1

    command = './RNAfold --noPS -i sequences_RNAFold.txt'
    mfe_process = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
    folds_data = mfe_process.communicate()
    lines = folds_data[0].split('\n')
    # print folds_data
    mfe_values = []
    for line in lines:
        if '-' in line:
            temp = sub('[()]','',line)
            mfe_value_temp = '-' + temp.split('-')[1]
            # print mfe_value_temp
            mfe_values.append(mfe_value_temp)
        elif '0.00' in line:
            mfe_value_temp = '0.00'
            mfe_values.append(mfe_value_temp)

    with open(reference_sequence_file) as original_file:
        with open(out_file,'w+') as out_file:
            index_count = 0
            for line in original_file:
                line = line.strip()
                # print mfe_values
                columns = line.split('\t')
                if not 'N' in columns[3]:
                    # print index_count,len(mfe_values)
                    replacement = columns[0] + '\t' + columns[1] + '\t' + columns[2] + '\t' + mfe_values[index_count]
                    replacement = replacement.strip()
                    index_count += 1
                    out_file.write("%s\n"%(replacement))

#
# sequence_file = 'refSeqNeg.txt'
# Compute_MFE(sequence_file)
# os.remove('sequences_RNAFold.txt')
