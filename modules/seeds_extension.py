#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
https://github.com/bichangwei/PMAT

This is part of the PMAT. Connection relations are extended using 
candidate seeds and filtered according to the depth of the contig 
and the depth of the link.
"""

# from tqdm import tqdm
import os
from progressbar import ProgressBar, ProgressBar, Percentage, Bar


class Extend_seeds:

    def __init__(self, all_connections, id_depth, id_length, link_depth, dynamic_sampleIDs, simple_pairs, proleptic_connections, minPath, genomesize, assemblysize):
        self.all_connections = all_connections
        self.id_depth = id_depth
        self.id_length = id_length
        self.link_depth = link_depth
        self.dynamic_sampleIDs = dynamic_sampleIDs
        self.simple_pairs = simple_pairs
        self.proleptic_connections = proleptic_connections
        self.minPath = minPath
        # self.minDep = minDep
        # self.minLen = minLen
        self.nucl_contig_depth = assemblysize / genomesize

    def _seed_extend(self, cycl):
        '''
        dynamic_sampleIDs: a dynamic list, and the contigID of the UTR connection is added by loop.
        The initial element is the seed ID of the simplified version if seed='contig00002' dynamic_sampleIDs=['2']
        Carry on according to the seed sequence provided
        Contig with seed ID were connected by cycling
        ContigID stores the list
        Join proleptic_connection to a dictionary  proleptic_connection
        '''

        widgets = [f'Extension No.{cycl}: ', Percentage(), ' ', Bar('#'), ' ']
        # progree_connections = tqdm(self.all_connections, desc="Extension No.%d" % cycl, ascii=True, bar_format="{l_bar}{bar}")
        # for connection in progree_connections:
        for connection in ProgressBar(widgets=widgets)(self.all_connections):
            for left_contig, right_contig in connection.items():
                left_contig_id = left_contig.split()[0]
                left_contig_edge = left_contig.split()[1]
                right_contig_id = right_contig.split()[0]
                right_contig_edge = right_contig.split()[1]

                if self.minPath is None:
                    minPath = min(float(self.id_depth[left_contig_id]), float(self.id_depth[right_contig_id])) / 5
                else:
                    minPath = self.minPath

                # At first dynamic_sampleIDs had only the 4004 seed
                for sampleID in self.dynamic_sampleIDs:  
                    if int(left_contig_id) == int(sampleID):
                        # Filter based on the depth and length of the contig
                        # if (float(self.id_depth[left_contig_id]) > self.nucl_contig_depth*2 or int(self.id_length[left_contig_id]) < self.minLen) and (float(self.id_depth[right_contig_id]) > self.nucl_contig_depth*2 or int(self.id_length[right_contig_id]) < self.minLen):
                        if float(self.id_depth[left_contig_id]) > self.nucl_contig_depth*2 and float(self.id_depth[right_contig_id]) > self.nucl_contig_depth*2:
                            # Filter based on the path depth
                            if int(self.link_depth[(left_contig_id, right_contig_id)]) > float(minPath):
                                proleptic_connection = {}
                                source = ""  # source :"contig*" + "edge"
                                source = left_contig_id + " " + left_contig_edge
                                target = ""  # target :"contig*" + "edge"
                                target = right_contig_id + " " + right_contig_edge
                                proleptic_connection[source] = target
                                self.proleptic_connections.append(proleptic_connection)
                                # Dynamic list proleptic_connections, the contig
                                # pairs generated in each cycle are stored in the list (the list composed of the dictionary proleptic_connection
                                if right_contig_id not in self.dynamic_sampleIDs:
                                    self.dynamic_sampleIDs.append(right_contig_id)  # Dynamic list dynamic_sampleIDs, each cycle corresponding to
                                    # the new contigID into the dynamic list, extend down
                                    pass

                            pass
                    elif int(right_contig_id) == int(sampleID):
                        # Filter based on the depth and length of the contig
                        if float(self.id_depth[left_contig_id]) > self.nucl_contig_depth*2 and float(self.id_depth[right_contig_id]) > self.nucl_contig_depth*2:
                            # Filter based on the path depth
                            if int(self.link_depth[(left_contig_id, right_contig_id)]) > float(minPath):
                                proleptic_connection = {}
                                source_full = ""
                                source_full = right_contig_id + " " + right_contig_edge
                                target_full = ""
                                target_full = left_contig_id + " " + left_contig_edge
                                proleptic_connection[source_full] = target_full
                                self.proleptic_connections.append(proleptic_connection)
                                if left_contig_id not in self.dynamic_sampleIDs:
                                    self.dynamic_sampleIDs.append(left_contig_id)
                                    pass
                                pass
                        pass
        # progree_connections.close()

        return self.dynamic_sampleIDs, self.simple_pairs, self.proleptic_connections


    def update_seed_extend(self):
        '''
        Automatically  updates to  seed_extend
        '''
        set1 = set()  #Set a set (non-duplicate) that stores contigs in all loops
        initial_connections= []  # The list stores contig connections
        n = 1
        cycles = 1
        while n != 0:
            self._seed_extend(cycles)
            for nl in self.proleptic_connections:
                for left, right in nl.items():
                    r_temp_line = {}#The opposite of this contig_connection
                    r_temp_line[right] = left
                t = tuple(nl.items())
                t1 = tuple(r_temp_line.items())
                if t not in set1 and t1 not in set1:
                    set1.add(t)
                    initial_connections.append(nl)
                pass
            if n != len(set1):
                n = len(set1)
            else:
                n = 0
            cycles += 1
        return initial_connections