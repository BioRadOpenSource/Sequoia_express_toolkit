#!/usr/bin/env python3
# coding: utf-8

"""
Downsamples fastq files
Accepts single R1 or R1/R2 pairs
version launched: February 8, 2021
edited by Evan Bursey to remove the sizename from the output file filename
"""

import os
import sys
import time
import io
import subprocess
import random
import gzip
import shutil

def generatename(inputfilename,sizename):
    if '_R1' in inputfilename:
        delimiter = '_R1'
    elif '_R2' in inputfilename:
        delimiter = '_R2'
        
    splitloc = len(inputfilename.split('/')) - 1
    basefilename = inputfilename.split('/')[splitloc].split('.')[0]
    filename = basefilename.split(delimiter)[0] + delimiter + basefilename.split(delimiter)[1] +'.fastq'
#    filename = basefilename.split(delimiter)[0] + '_' + sizename + delimiter + basefilename.split(delimiter)[1] +'.fastq'

    return filename

if __name__ == "__main__":
    
    version_num = 1.4
    
    start_time = time.time()
    
    print('*'*100)
    print('FASTQ downsampler')
    print('2020 Mark Shulewiz Bio-Rad Laboratories, Inc.')
    print('version ' + str(version_num))
    print('*'*100)
    
# parse arguments
    cli_input_filenames = []
    input_filenames = []
    output_dir = os.getcwd()
    input_number_of_reads = 0
    output_filenames = []

    if len(sys.argv) == 2:
        if sys.argv[1] == '--help' or sys.argv[1] == '--h':
            print('*'*10)
            print('Usage: ' + os.path.basename(__file__) + ' --I [input files] --O [output directory] --size [number of reads]')
            print('Number of reads accepts M for millions and K for thousands (ie: 1M, 100K)')
            print('*'*10)
            exit()
        else:
            print('Error 1 command line: unknown argument {}'.format(sys.argv[1]))
            exit()
    elif len(sys.argv) > 2:
        input_trigger = False
        for a in range(0, len(sys.argv)):
            if input_trigger == True:
                if sys.argv[a][0:2] == '--':
                    input_trigger = False
                else:
                    cli_input_filenames.append(sys.argv[a])
            if input_trigger == False:
                if a < len(sys.argv) - 1:
                    if sys.argv[a] == '--I':
                        input_trigger = True
                    elif sys.argv[a] == '--O':
                        output_dir = sys.argv[a + 1]
                    elif sys.argv[a].lower() == '--size':
                        input_number_of_reads = sys.argv[a + 1]
    else:
        print('Error 2: something is missing in the commandline')
        exit()
        
# process number of reads
    if input_number_of_reads != 0:
        try:
            number_of_reads = int(input_number_of_reads)
        except:
            if input_number_of_reads[len(input_number_of_reads) - 1] == 'M':
                number_of_reads = int(input_number_of_reads[:len(input_number_of_reads) - 1])*1000000
            elif input_number_of_reads[len(input_number_of_reads) - 1] == 'K':
                number_of_reads = int(input_number_of_reads[:len(input_number_of_reads) - 1])*1000
            else:
                print('Error 3: invalid entry for size: {}'.format(input_number_of_reads))
                exit()
    else:
        print('Error 4: invalid entry for size: {}'.format(input_number_of_reads))
        exit()
        
# format output number of reads  
    if str(number_of_reads)[-3:] == '000':
        if str(number_of_reads)[-6:] == '000000':
            readable_reads = str(number_of_reads)[:-6] + 'M'
        else:
            readable_reads = str(number_of_reads)[:-3] + 'K'
    else:
        readable_reads = str(number_of_reads)
       
# check for output directory and create if necessary
    if output_dir != os.getcwd():
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
            
    print('Output directory: {}'.format(output_dir))
    print('*'*10)
    
    print('Number of output reads: {}'.format(readable_reads))
    print('*'*10)
        
# process input filenames
    input_filenames = []
    if len(cli_input_filenames) == 0:
        print('Error 5: Missing input file(s)')
        exit()
    elif len(cli_input_filenames) == 1:
        if 'R1' not in cli_input_filenames[0]:
            print('Error 6: {} is not an R1 file'.format(cli_input_filenames[0]))
        else:
            input_filenames.append([cli_input_filenames[0]])
    elif len(cli_input_filenames) > 1:
        for b in range(0, len(cli_input_filenames)):
            if 'R1' in cli_input_filenames[b]:
                if b < len(cli_input_filenames) - 1:
                    if 'R2' in cli_input_filenames[b + 1]:
                        if cli_input_filenames[b].replace('R1','R2') == cli_input_filenames[b + 1]:
                            input_filenames.append([cli_input_filenames[b],cli_input_filenames[b + 1]])
                    else:
                        input_filenames.append([cli_input_filenames[b]])
                else:
                    input_filenames.append([cli_input_filenames[b]])
            elif 'R2' not in cli_input_filenames[b]:
                print('Error 7: {} is not an R2 file'.format(cli_input_filenames[b]))  
                
# check input file names for .gz extension
    for input_filename in input_filenames:
        if len(input_filename) == 1:
            if '.gz' not in input_filename[0]:
                print('Error 8: Input file {} does not end in .gz'.format(input_filename[0]))
                exit()
        elif len(input_filename) == 2:
            if '.gz' not in input_filename[0]:
                print('Error 9: Input file {} does not end in .gz'.format(input_filename[0]))
                exit()
            if '.gz' not in input_filename[1]:
                print('Error 10: Input file {} does not end in .gz'.format(input_filename[1]))
                exit()
                                             
    print('Processing files')
    print('*'*5)
    
    for input_filename in input_filenames:
        
 # create temporary directory
        tmp_dir = output_dir + '/tmp'
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)
        os.makedirs(tmp_dir)
        
# process single reads
        if len(input_filename) == 1:
            print('Processing R1 file: {}'.format(input_filename[0]))
            
# unzip fastq
            print('unzipping file...')
    
            cmd = 'pigz -d -c {} | paste - - - - > {}'.format(input_filename[0],tmp_dir + '/compressed.tsv')
            process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
            process.wait()
          
            inputreadnum = subprocess.run(['wc','-l',tmp_dir + '/compressed.tsv'], stdout=subprocess.PIPE)
            inputreadnum = int(str(inputreadnum.stdout).split()[0].strip("b'"))

            if inputreadnum < int(number_of_reads):
                print('Input number of reads ({}) less than downsampled number of reads ({}); will not downsample'.format(str(inputreadnum),number_of_reads))

# pass reads through unchanged so downstream still receives an output file
                r1_output_filename = output_dir + '/' + generatename(input_filename[0],readable_reads) + '.gz'
                print('copying input reads to output unchanged...')
                shutil.copyfile(input_filename[0], r1_output_filename)
            else:

# shuffle
                print('shuffling reads...')
                cmd = 'shuf {} > {}'.format(tmp_dir + '/compressed.tsv',tmp_dir + '/shuffled.tsv')
                process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
                process.wait()

# downsample
                print('downsampling reads...')
                cmd = 'head -{} {} > {}'.format(number_of_reads,tmp_dir + '/shuffled.tsv',tmp_dir + '/downsampled.tsv')
                process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
                process.wait()

# generate output filename
                r1_output_filename = output_dir + '/' + generatename(input_filename[0],readable_reads) + '.gz'

# split compressed file and write output
                print('writing output files...')
                cmd = 'cat {} | cut -f 1-4 | tr "\t" "\n" | pigz > {}'.format(tmp_dir + '/downsampled.tsv',r1_output_filename)
                process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
                process.wait()

# process read pairs
        elif len(input_filename) == 2:
            print('Processing R1 file: {}'.format(input_filename[0]))
            print('Processing R2 file: {}'.format(input_filename[1]))
            
# unzip fastq files
            print('unzipping files...')
            cmd = 'paste <(pigz -d -c {} | paste - - - -) <(pigz -d -c {} | paste - - - -) > {}'.format(input_filename[0],input_filename[1],tmp_dir + '/compressed.tsv')
            process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
            process.wait()
            
            inputreadnum = subprocess.run(['wc','-l',tmp_dir + '/compressed.tsv'], stdout=subprocess.PIPE)
            inputreadnum = int(str(inputreadnum.stdout).split()[0].strip("b'"))

            if inputreadnum < int(number_of_reads):
                print('Input number of reads ({}) less than downsampled number of reads ({}); will not downsample'.format(str(inputreadnum),number_of_reads))

# pass reads through unchanged so downstream still receives output files
                r1_output_filename = output_dir + '/' + generatename(input_filename[0],readable_reads) + '.gz'
                r2_output_filename = output_dir + '/' + generatename(input_filename[1],readable_reads) + '.gz'
                print('copying input reads to output unchanged...')
                shutil.copyfile(input_filename[0], r1_output_filename)
                shutil.copyfile(input_filename[1], r2_output_filename)
            else:

# shuffle
                print('shuffling reads...')
                cmd = 'shuf {} > {}'.format(tmp_dir + '/compressed.tsv',tmp_dir + '/shuffled.tsv')
                process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
                process.wait()

# downsample
                print('downsampling reads...')
                cmd = 'head -{} {} > {}'.format(number_of_reads,tmp_dir + '/shuffled.tsv',tmp_dir + '/downsampled.tsv')
                process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
                process.wait()

# generate output filenames
                r1_output_filename = output_dir + '/' + generatename(input_filename[0],readable_reads) + '.gz'
                r2_output_filename = output_dir + '/' + generatename(input_filename[1],readable_reads) + '.gz'

# split compressed file and write output
                print('writing output files...')
                cmd = 'cat {} | tee >(cut -f 1-4 | tr "\t" "\n" | pigz > {}) | cut -f 5-8 | tr "\t" "\n" | pigz > {}'.format(tmp_dir + '/downsampled.tsv',r1_output_filename,r2_output_filename)
                process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
                process.wait()
                                          
# delete temporary directory
        shutil.rmtree(tmp_dir)
    
        print('*'*10)
                
    elapsed_time_sec = time.time() - start_time
    elapsed_time_minutes = (elapsed_time_sec - (elapsed_time_sec % 60))/60
    print('finished, elapsed time: {} min {} sec'.format(elapsed_time_minutes,int(elapsed_time_sec % 60)))
