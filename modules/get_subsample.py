# !/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
https://github.com/bichangwei/PMAT

This is part of the PMAT. Random interception of 
fasta sequences, with a random seed default of 6.
"""

import sys
import random
import time
import os
from log import Log
import gzip
from check_file import mkdir_file_path
# from tqdm import tqdm
from progressbar import ProgressBar, Percentage, Bar

log = Log()

widgets = ['Progress: ', Percentage(), ' ', Bar('#'), ' ']

def is_gzipped_file(filename):
    root, ext = os.path.splitext(filename)
    if ext == '.gz':
        return True
    else:
        return False

def subsample(output, corrected_seq, factor, seed=6):
    mkdir_file_path(f'{output}/subsample')
    if factor < 1:
        log.Info(f'The seed of the random number generator is {seed}')
        log.section_header('Random selection sequence start...')
        assembly_seq = os.path.join(f'{output}/subsample', f'assembly_seq_subset.{factor}.fasta')
        if os.path.exists(assembly_seq):
            log.Error(assembly_seq)
            log.Error('The file already exists and rewrites')
        # build index
        dict_contig={}
        dict_seq = {}
        count=0 # Statistical sequence number
        if is_gzipped_file(corrected_seq):
            with gzip.open(corrected_seq,'rt') as unzip_fa:
                    start_time = time.time()
                    count_seq = ""
                    # for line in tqdm(unzip_fa, desc="Random selected subsets ", ascii=True):
                    for line in unzip_fa:
                        if line.startswith('>'):
                            count_seq = ""
                            count += 1
                            dict_contig[count]=line.strip()
                        else:
                            count_seq = count_seq + line.strip()
                            dict_seq[count] = count_seq
                        elapsed_time = time.time() - start_time
                        print(f">>>>>> Random selected subsets {elapsed_time:.2f}s <<<<<<", end="\r")
        else:
            with open(corrected_seq, 'r') as fp_read:
                start_time = time.time()
                count_seq = ""
                # for line in tqdm(fp_read, desc="Random selected subsets ", ascii=True):
                for line in fp_read:
                    if line.startswith('>'):
                        count_seq = ""
                        count += 1
                        dict_contig[count]=line.strip()
                    else:
                        count_seq = count_seq + line.strip()
                        dict_seq[count] = count_seq
                    elapsed_time = time.time() - start_time
                    print(f">>>>>> Random selected subsets {elapsed_time:.2f}s <<<<<<", end="\r")

        # Random number sequence generation
        length = int(factor*count)
        indexlist_random=[]

        i=1
        while i<= length:
            random.seed(seed+i)
            indexlist_random.append(1+int(random.random() * count))
            i=i+1
        with open(assembly_seq,'w') as fp_result:
            # for indexid in tqdm(indexlist_random, desc="Extracting subsets ", ascii=True, bar_format="{l_bar}{bar}"):
            for indexid in ProgressBar(widgets=widgets)(indexlist_random):
                fp_result.write(dict_contig[indexid]+'\n')
                fp_result.write(dict_seq[indexid]+'\n')
            fp_result.close()

        assembly_seq_path = assembly_seq

        log.section_tail("Random selection sequence end.")
        log.get_path(f'The random sequences path : {assembly_seq_path}')
    
    else:
        assembly_seq = os.path.join(f'{output}/subsample', f'assembly_seq_subset.1.0.fasta')
        if is_gzipped_file(corrected_seq):
            with open(assembly_seq, 'w') as fp_result:
                with gzip.open(corrected_seq,'rt') as unzip_fa:
                        start_time = time.time()
                        # for line in tqdm(unzip_fa, desc="Random selected subsets ", ascii=True):
                        for line in unzip_fa:
                            line = line.strip()
                            fp_result.write(line+'\n')
                            elapsed_time = time.time() - start_time
                            print(f">>>>>> Random selected subsets {elapsed_time:.2f}s <<<<<<", end="\r")
        else:
            with open(assembly_seq, 'w') as fp_result:
                with open(corrected_seq, 'r') as fp_read:
                    start_time = time.time()
                    # for line in tqdm(fp_read, desc="Random selected subsets ", ascii=True):
                    for line in fp_read:
                        line = line.strip()
                        fp_result.write(line+'\n')
                        elapsed_time = time.time() - start_time
                        print(f">>>>>> Random selected subsets {elapsed_time:.2f}s <<<<<<", end="\r")

        assembly_seq_path = assembly_seq

    return assembly_seq_path

if __name__ == '__main__':
    factor = float(sys.argv[3])
    subsample(sys.argv[1], sys.argv[2], factor)