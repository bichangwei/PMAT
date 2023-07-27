#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
https://github.com/bichangwei/PMAT

This is part of the PMAT. The connection relationships are generated 
into gfa files for further visualisation and analysis.
"""

import re
import time

class AssemblyGraph:

    def __init__(self, initial_connections, file_data_fna_name, id_length, id_depth, simple_pairs, contig_dict=None):
        self.del_connections = initial_connections
        self.initial_connections = initial_connections.copy()
        self.file_data_fna_name = file_data_fna_name
        self.id_length = id_length
        self.id_depth = id_depth
        self.simple_pairs = simple_pairs
        self.contig_dict = contig_dict
        pass

    def delete_single_edge(self):
        '''
        Delete contigs with no connection information at both sides
        '''
        n = 1
        m = 0
        simple_sampleIDs = set() # store the contig id after simplified
        self.del_connections = self.initial_connections
        while n != m:
            set2 = set()
            contig_edge = {} # store id and matching edge
            n = len(self.del_connections)
            for init_connection in self.del_connections:
                for left_cnt, right_cnt in init_connection.items():
                    left_id = left_cnt.split()[0]
                    left_edge = left_cnt.split()[1]
                    right_id = right_cnt.split()[0]
                    right_edge = right_cnt.split()[1]
                    if left_id not in set2:
                        set2.add(left_id)
                        contig_edge[left_id] = set(list(left_edge)[0])

                        if right_id not in set2:
                            set2.add(right_id)
                            contig_edge[right_id] = set(list(right_edge)[0])
                            pass
                        else:
                            contig_edge[right_id].add(list(right_edge)[0])
                            pass
                    else:
                        contig_edge[left_id].add(list(left_edge)[0])
                        if right_id not in set2:
                            set2.add(right_id)
                            contig_edge[right_id] = set(list(right_edge)[0])
                            pass
                        else:
                            contig_edge[right_id].add(list(right_edge)[0])
                            pass
            for cnt_id, cnt_edge in contig_edge.items():
                if len(cnt_edge) == 1:
                    set2.remove(cnt_id)
                pass
            for i, init_connection in enumerate(self.del_connections):
                for left_cnt, right_cnt in init_connection.items():
                    left_id = left_cnt.split()[0]
                    left_edge = left_cnt.split()[1]
                    right_id = right_cnt.split()[0]
                    right_edge = right_cnt.split()[1]
                    if left_id not in set2 or right_id not in set2:
                        del self.del_connections[i]
                    else:
                        simple_sampleIDs.update([left_id, right_id])
                    pass
            m = len(self.del_connections)
            simple_connectios = self.del_connections
            pass
        return simple_connectios, simple_sampleIDs
    
    def Target_seeds(self, connections):
        '''
        Return initial seeds and main seeds according to the provided parameters
        Traget_seeds(self.initial_connections) -> initial_seeds 
        Traget_seeds(self.delete_single_edge()[0]) -> main_seeds
        '''
        if connections == 'init':
            connections = self.initial_connections
        else:
            connections = self.delete_single_edge()[0]

        seeds = set()  # A collection used to store seeds
        for contig_connection in connections:
            for left, right in contig_connection.items():
                seeds.add(self.simple_pairs[str(left.split()[0])])
                seeds.add(self.simple_pairs[str(right.split()[0])])
                
        #seeds = sorted(seeds)
        seeds=list(seeds)
        for i in range(0, len(seeds)):
            seeds[i] = re.sub('contig0*', '', seeds[i])
            seeds[i]=int(seeds[i])
        seeds=sorted(seeds)
        for i in range(0, len(seeds)):
            seeds[i] = str(seeds[i])
        return seeds

    def save_gfa(self, gfa, seeds, init_main_connections):
        '''
        Return initial_gfa and main_gfa according to the provided parameters
        save_gfa(*_initial.gfa, Traget_seeds(initial_seeds)) -> *_initial.gfa
        save_gfa(*_main.gfa, Traget_seeds(simple_seeds)) -> *_main.gfa
        '''

        start_time = time.time()

        for i, ctg in enumerate(seeds):
            ctg_RC = int(self.id_length[str(ctg)])*float(self.id_depth[str(ctg)])
            gfa.write("S\t{}\t{}\tLN:i:{}\tRC:i:{}\n".format(ctg, self.contig_dict[str(ctg)], self.id_length[str(ctg)], ctg_RC))

        if init_main_connections == 'init':
            init_main_connections = self.initial_connections
        else:
            init_main_connections = self.delete_single_edge()[0]

        #Output connections between seeds
        for contig_connection in init_main_connections:
            elapsed_time = time.time() - start_time
            for left_contig, right_contig in contig_connection.items():
                gfa.write('L' + '\t' + left_contig.split()[0] + '\t')
                left_edge = left_contig.split()[1]
                right_edge = right_contig.split()[1]
                if left_contig.split()[1]!=right_contig.split()[1]:
                    if(left_contig.split()[1]=="3'"):
                        left_edge='+'
                        right_edge='+'
                    else:
                        left_edge='-'
                        right_edge='-'
                else:
                    if (left_contig.split()[1] == "3'"):
                        left_edge = '+'
                        right_edge = '-'
                    else:
                        left_edge= '-'
                        right_edge= '+'
                gfa.write(left_edge + '\t')
                gfa.write(right_contig.split()[0] + '\t')
                gfa.write(right_edge + '\t' + '0M' + '\n')
            pass
