# !/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
https://github.com/bichangwei/PMAT

This is part of PMAT and is used to convert gfa files to fasta.
"""

def convert_seq(seq, flag):
    if flag == "F":
        final_seq = seq
    elif flag == "R":
        seq = seq[::-1]
        replacements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 'T', 't': 'A', 'c': 'G', 'g': 'C'}
        final_seq = ''.join(replacements[base] for base in seq)
    return final_seq

def gfa2fa(gfa, output):
    with open(gfa, 'r') as gr:
        gfa_info = gr.readlines()

    gfa_seq = {}

    for gfa_line in gfa_info:
        gfa_line = gfa_line.strip().split("\t")
        if gfa_line[0] == "S":
            gfa_seq[str(gfa_line[1])] = gfa_line[2]

    
    unrep_seed_set = set()
    unrep_seed_list = []
    unrep_seed_dict = {}
    loop_seq = ''
    for gfa_line in gfa_info:
        gfa_line = gfa_line.strip().split("\t")
        if gfa_line[0] == "L":
            ctgL = gfa_line[1]
            ctgR = gfa_line[3]
            if str(ctgL) not in unrep_seed_set and str(ctgR) not in unrep_seed_set:
                unrep_seed_set.add(str(ctgL))
                unrep_seed_set.add(str(ctgR))
                unrep_seed_list.append(str(ctgL))
                unrep_seed_list.append(str(ctgR))

                if gfa_line[2] == "+":
                    unrep_seed_dict[ctgL] = "F"
                    loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgL)], "F")
                    if gfa_line[2] == gfa_line[4]:
                        unrep_seed_dict[str(ctgR)] = "F"
                        loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgR)], "F")
                    else:
                        unrep_seed_dict[str(ctgR)] = "R"
                        loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgR)], "R")
                else:
                    unrep_seed_dict[ctgL] = "R"
                    loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgL)], "R")
                    if gfa_line[2] == gfa_line[4]:
                        unrep_seed_dict[str(ctgR)] = "R"
                        loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgR)], "R")
                    else:
                        unrep_seed_dict[str(ctgR)] = "F"
                        loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgR)], "F")

            elif str(ctgL) in unrep_seed_set and str(ctgR) not in unrep_seed_set:
                unrep_seed_set.add(str(ctgR))
                unrep_seed_list.append(str(ctgR))

                if gfa_line[2] == gfa_line[4]:
                    loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgR)], unrep_seed_dict[str(ctgL)])
                    unrep_seed_dict[str(ctgR)] = unrep_seed_dict[str(ctgL)]
                elif gfa_line[2] != gfa_line[4]:
                    if unrep_seed_dict[str(ctgL)] == "F":
                        loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgR)], "R")
                        unrep_seed_dict[str(ctgR)] = "R"
                    else:
                        loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgR)], "F")
                        unrep_seed_dict[str(ctgR)] = "F"

            elif str(ctgL) not in unrep_seed_set and str(ctgR) in unrep_seed_set:
                unrep_seed_set.add(str(ctgL))
                unrep_seed_list.append(str(ctgL))

                if gfa_line[2] == gfa_line[4]:
                    loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgL)], unrep_seed_dict[str(ctgR)])
                    unrep_seed_dict[str(ctgL)] = unrep_seed_dict[str(ctgR)]
                elif gfa_line[2] != gfa_line[4]:
                    if unrep_seed_dict[str(ctgR)] == "F":
                        loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgL)], "R")
                        unrep_seed_dict[str(ctgL)] = "R"
                    else:
                        loop_seq = loop_seq + convert_seq(gfa_seq[str(ctgL)], "F")
                        unrep_seed_dict[str(ctgL)] = "F"    

    with open(output, 'w') as ow:
        header = '-'.join(unrep_seed_list)
        ow.write(f'>{header}\n{loop_seq}')
