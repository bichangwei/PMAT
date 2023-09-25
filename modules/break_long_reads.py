#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
https://github.com/bichangwei/PMAT

This is part of the PMAT. Truncation of reads greater than 
30,000bp, default truncation length is 20,000bp.
"""

import os
from log import Log

log = Log()

def BreakLongReads(corrected_seq, break_length=20000):
    '''
    Newbler can only assemble seqs with a length below 30000bp;
    Cutoff for corrected ONT, CLR or HiFi;
    corrected_seq -> corrected ONT, CLR or HiFi
    break_length -> default:20000bp
    '''
    
    log.section_header("Reads break start ...")

    seq_count = 0
    name = 1

    cut_length = f'cut{int(break_length * 0.001)}K'
    assembly_seq =  os.path.join(os.path.dirname(corrected_seq), f'assembly_seq.{cut_length}.fasta')
    if os.path.exists(assembly_seq):
        log.Error(assembly_seq)
        log.Error('The file already exists and rewrites')
    with open(corrected_seq, 'r') as inseq, open(assembly_seq, 'w') as outseq:
        for line in inseq:
            line = line.rstrip()  # remove trailing newline character
            
            if line.startswith('>'):
                seq_count += 1
            else:
                seq = line.strip()
                seq_len = len(seq)

                if seq_len < 30000:
                    outseq.write(f'>{name}\n{seq}\n')
                    name += 1
                else:
                    num = seq_len // break_length
                    for i in range(num+1):
                        data = seq[i*break_length:(i+1)*break_length]
                        outseq.write(f'>{name}\n{data}\n')
                        name += 1
    
    log.Info(f'The raw sequence number is: {seq_count}')
    log.Info(f'The processed sequence number is: {name-1}')

    log.section_tail("Reads break end.")
    assembly_seq_cut_path = assembly_seq

    return assembly_seq_cut_path