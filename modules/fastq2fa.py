# !/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
https://github.com/bichangwei/PMAT

This is part of PMAT and is used to convert fastq files to fasta.
"""


import sys
import os
import gzip
# from tqdm import tqdm
from log import Log
from check_file import check_file_format
import time

log = Log()
def fq2fa(filename):
    # root, ext = os.path.splitext(filename)
    # if ext == '.gz':
    #     log.section_header('Unzip file...')
    #     faroot, gaext = os.path.splitext(root)
    #     with gzip.open(filename, 'rb') as raw_data:
    #         seq_store = raw_data.readlines()
    #     firstline = seq_store[0].decode().strip()
    #     log.section_header('convert the fastq to fasta...')
    #     if firstline.startswith(">"):
    #         pass
    #     elif firstline.startswith("@"):
    #         with open(faroot + '.fa', 'w') as ungzip_fa:
    #             for index in tqdm(range(0, len(seq_store), 4), desc="fastq2fa ", ascii=True, bar_format="{l_bar}{bar}"):
    #                 ungzip_fa.write(f'>{index//4+1}\n')
    #                 ungzip_fa.write(seq_store[index+1].decode())
    #     log.section_tail('convert end.')
    #     fa_seq = f'{root}+.fa'

    root, ext = os.path.splitext(filename)
    if ext == '.gz':
        log.section_header('Unzip file...')
        faroot, gaext = os.path.splitext(root)
        with gzip.open(filename, 'rt') as raw_data:  # Open file for reading with gzip and text mode
            with open(faroot + '.fa', 'w') as ungzip_fa:
                if check_file_format(filename) == 'fastq':
                    id_num = 1
                    line_num = 0
                    start_time = time.time()
                    # for line in tqdm(raw_data, desc="Convert format ", ascii=True):
                    for line in raw_data:
                        line_num += 1
                        id_line = 1 + (id_num - 1)*4
                        line = line.strip()
                        if line_num == id_line:
                            ungzip_fa.write(f'>{id_num}\n')  # Write FASTA header
                        elif line_num == id_line + 1:
                            id_num += 1
                            ungzip_fa.write(line+'\n')
                        elapsed_time = time.time() - start_time
                        print(f">>>>>> Convert format {elapsed_time:.2f}s <<<<<<", end="\r")

                elif check_file_format(filename) == 'fasta':
                    id_num = 0
                    start_time = time.time()
                    # for line in tqdm(raw_data, desc="Unzip file ", ascii=True):
                    for line in raw_data:
                        line = line.strip()
                        if line.startswith('>'):
                            id_num += 1
                            ungzip_fa.write(f'>{id_num}\n')  # Write FASTA header
                        else:
                            ungzip_fa.write(line+'\n')
                        elapsed_time = time.time() - start_time
                        print(f">>>>>> Convert format {elapsed_time:.2f}s <<<<<<", end="\r")

                else:
                    log.Info('Input data error')

        log.section_tail('convert end.')
        fa_seq = f'{faroot}.fa'

    elif check_file_format(filename) == 'fastq':
        log.section_header('convert the fastq to fasta...')
        with open(filename, 'r') as raw_data:
            with open(root + '.fa', 'w') as ungzip_fa:
                id_num = 1
                line_num = 0
                start_time = time.time()
                # for line in tqdm(raw_data, desc="Convert format ", ascii=True):
                for line in raw_data:
                    line_num += 1
                    id_line = 1 + (id_num - 1)*4
                    line = line.strip()
                    if line_num == id_line:
                        ungzip_fa.write(f'>{id_num}\n')  # Write FASTA header
                    elif line_num == id_line + 1:
                        id_num += 1
                        ungzip_fa.write(line+'\n')
                    elapsed_time = time.time() - start_time
                    print(f">>>>>> Convert format {elapsed_time:.2f}s <<<<<<", end="\r")

        log.section_tail('convert end.')
        fa_seq = f'{root}.fa'

    elif check_file_format(filename) == 'fasta':
        fa_seq = filename
        pass

    else:
        log.Info(log.bold_purple('file error'))
    
    return fa_seq
 
if __name__ == '__main__':
    fq2fa(sys.argv[1])
