"""
https://github.com/bichangwei/PMAT

This script is a part of PMAT and is used for 
automatically resolving loops in the master results.
"""

# import numpy as np
# from sklearn.mixture import GaussianMixture
import re
import os
import shutil
# from log import Log
from collections import Counter
import subprocess
import itertools
import random
from log import Log 

log = Log()

if os.path.exists(os.path.abspath(os.path.dirname(__file__) + '/../Conserved_PCGs_db/Plant_conserved_mtgene_nt.fa')):
    db_path = os.path.join(os.path.abspath(os.path.dirname(__file__) + '/../Conserved_PCGs_db'), "Plant_conserved_mtgene_nt.fa")
elif shutil.which("PMAT"):
    db_path = os.path.join(os.path.abspath(os.path.dirname(shutil.which("PMAT")) + '/../Conserved_PCGs_db'), "Plant_conserved_mtgene_nt.fa")
else:
    log.Warning("Please check if the Conserved_PCGs_db file is installed correctly!")


class MasterLoops:
    def __init__(self, id_depth, id_length, connections, main_gfa,output):
        # self.contiggraph = contiggraph
        self.id_depth = id_depth
        self.id_length = id_length
        self.connections = connections
        self.main_gfa = main_gfa
        self.output = output
        self.conserved_PCGs = "atp1 atp4 atp6 atp8 atp9 cob cox1 cox2 cox3 mttB matR rpl2 rpl10 rpl16 rps1 rps3 rps4 rps7 rps12 nad1 nad2 nad3 nad4 nad4L nad5 nad6 nad7 nad9 ccmB ccmC ccmFn ccmFc sdh3 sdh4".split(' ')
        self._initialize_parser()

    def _initialize_parser(self):
        # self.simple_pairs, self.all_connections, self.id_depth, self.id_length, self.id_seq, self.link_depth = self._get_all_connections_from_datafile()
        self.median_seed, self.median_depth = self.get_median()
        self.main_length, self.main_idnum, self.main_seed, self.begin_seed, self.main_connections = self._get_all_connections_from_gfa()
        self.slim_connections, self.slim_main_idnum = self._check_bubble()
        self.final_contig_id, self.final_rep_contig, self.final_rep_edge = self.get_final_path("5'")
        self.head_contig_id, self.head_rep_contig, self.head_rep_edge = self.get_final_path("3'")

    def _PCGs_length(self):
        fna = open(db_path, 'r')
        PCGs_line = fna.readlines()
        fna.close()
        
        PCGs_seq = {}
        PCGs_len = {}
        index = 0
        # print(PCGs_line)
        while index < len(PCGs_line):
            if re.match('>', PCGs_line[index]):
                PCGs_seq[PCGs_line[index]] = PCGs_line[index+1]
            index = index+1

        for PCGs in self.conserved_PCGs:
            list_length = []

            for PCG_id, PCG_seq in PCGs_seq.items():
                if PCGs == str(PCG_id.strip().split('_')[2]):
                    list_length.append(len(PCG_seq))

            PCGs_len[PCGs] = int(sum(list_length) / len(list_length))

        return PCGs_len


    def _Run_blastn(self, qryseq):
            Blastn_command = ['blastn', '-db', db_path, '-query', qryseq, '-outfmt', '6', '-num_threads', '30', '-num_alignments', '1', '-max_hsps', '1']
            Blastn_process = subprocess.Popen(Blastn_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # output the result of blastn and error.
            blastn_out, blastn_err = Blastn_process.communicate(timeout=None) 

            # Error output
            if blastn_err:
                print('\nBLASTn encountered an error:\n' + blastn_err.decode())
            return blastn_out
    
    def condidate_seeds(self, qryseq):
        #Find the target contig
        PCGs_len = self._PCGs_length()
        blast_info = [] # [['atp1', 'contig00001', '90.9', '501'], ['apt2', 'contig00002', '891', '1030'] ...]
        redundant_seed = set()
        for line in self._Run_blastn(qryseq).decode().splitlines():
            lines = line.split()
            PCGs = re.sub(r".*_", "", lines[1])
            PCGs = re.sub(r"-.*", "", PCGs)
            # select contigs
            if PCGs in PCGs_len.keys() and int(lines[3]) > float(PCGs_len[PCGs]) * 0.9 and float(lines[2]) > 0.95:
                blast_info.append(str(re.sub('contig0*', '', lines[0])))
        return blast_info
    
    def get_median(self):
        '''
        find the median seed 
        '''
        temp_file = os.path.join(self.output, f'temp_blast.fa')
        with open(temp_file, 'w') as wf:
            with open(self.main_gfa, 'r') as fg:
                main_list = fg.readlines()
                for line in main_list:
                    if re.match('S', line):
                        lines = line.strip().split()
                        wf.write('>{}\n{}\n'.format(lines[1], lines[2]))
        blast_info = self.condidate_seeds(temp_file)
        if len(blast_info) > 0:
            main_depth = {seed: self.id_depth[seed] for seed in blast_info}
            sorted_value = sorted(main_depth.items(), key=lambda x: x[1])
            # sorted_value = sorted(main_dict.items(), key=lambda x: x[1])
            n = len(sorted_value)
            if n % 2 == 0:
                median_index = n // 2
            else:
                median_index = (n - 1) // 2

            median_seed, median_depth = sorted_value[median_index]
        else:
            log.Info("No result")
        os.remove(temp_file)
        return median_seed, median_depth
    
    def _get_all_connections_from_gfa(self):
        '''
        Reading the main gfa file;
        obtaining information on main contig lengths, depths and connection relationships;
        '''
        main_connections = self.connections
        idstore = []
        # main_length = {} #
        main_depth = {}
        rep_seed = []
        unrep_seed = []

        # for line in main_gfa:
        #     #
        #     if re.match('L', line):
        #         connection = {}
        #         lio = line.strip().split('\t')
        #         contig_and_UTR1 = lio[1] + " " + lio[2]
        #         contig_and_UTR2 = lio[3] + " " + lio[4]  # Put the corresponding contig and edge number into the list
        #         connection[contig_and_UTR1] = contig_and_UTR2
        #         # main_connections.append(connection)
        #         idstroe.extend([lio[1], lio[3]])
        for connection in main_connections:
            for lef, rig in connection.items():
                idstore.extend([lef.split(' ')[0], rig.split(' ')[0]])
        main_idnum = dict(Counter(idstore))
        for seed, rep in main_idnum.items():
            if rep > 2:
                rep_seed.append(seed)
            elif rep == 2:
                unrep_seed.append(seed)
            else:
                print(seed,rep)

        main_length = {seed: self.id_length[seed] for seed in idstore}
        main_depth = {seed: self.id_depth[seed] for seed in idstore}
        
        # find the longest seed as main seed
        main_seed = sorted(unrep_seed, key=lambda x: main_length[x], reverse=True)
        for seed in main_seed:
            if float(main_depth[seed]) > float(self.median_depth)*4:
                pass
            else:
                begin_seed = seed
                break

        return main_length, main_idnum, main_seed, begin_seed, main_connections

    
    def _check_bubble(self):
        del_num = []
        slim_connections = self.connections
        for num, connect in enumerate(slim_connections):
            for left_contig, right_contig in connect.items():
                left_contig_id, left_contig_edge = left_contig.split()
                right_contig_id, right_contig_edge = right_contig.split()

                if int(self.main_idnum[left_contig_id]) > 2 and int(self.main_idnum[right_contig_id]) > 2:
                    left_bubble_fact = []
                    right_bubble_fact = []
                    for bubble_connect in slim_connections[:num] + slim_connections[num+1:]:

                        for bubble_left_contig, bubble_right_contig in bubble_connect.items():
                            bubble_left_contig_id, bubble_left_contig_edge = bubble_left_contig.split()
                            bubble_right_contig_id, bubble_right_contig_edge = bubble_right_contig.split()
                            if int(left_contig_id) == int(bubble_left_contig_id) and bubble_left_contig_edge == left_contig_edge:
                                left_bubble_fact.append({bubble_right_contig_id:bubble_right_contig_edge})
                            elif int(left_contig_id) == int(bubble_right_contig_id) and bubble_right_contig_edge == left_contig_edge:
                                left_bubble_fact.append({bubble_left_contig_id:bubble_left_contig_edge})
                            elif int(right_contig_id) == int(bubble_left_contig_id) and bubble_left_contig_edge == right_contig_edge:
                                right_bubble_fact.append({bubble_right_contig_id:bubble_right_contig_edge})
                            elif int(right_contig_id) == int(bubble_right_contig_id) and bubble_right_contig_edge == right_contig_edge:
                                right_bubble_fact.append({bubble_left_contig_id:bubble_left_contig_edge})
                    
                    if len(left_bubble_fact) > 0 and len(right_bubble_fact) > 0:
                        for lef in left_bubble_fact:
                            for left_bubble_id, left_bubble_edge in lef.items():
                                for rig in right_bubble_fact:
                                    for right_bubble_id, right_bubble_edge in rig.items():
                                        if int(left_bubble_id) == int(right_bubble_id) and left_bubble_edge != right_bubble_edge:
                                            del_num.append(num)
        if len(del_num) > 0:
            for del_index in sorted(del_num, reverse=True):
                del slim_connections[del_index]
        
        del_num = []
        for num, connect in enumerate(slim_connections):
            for left_contig, right_contig in connect.items():
                left_contig_id, left_contig_edge = left_contig.split()
                right_contig_id, right_contig_edge = right_contig.split()
                if float(self.id_depth[left_contig_id]) > float(self.median_depth)*3 and float(self.id_depth[right_contig_id]) > float(self.median_depth)*3:

                    del_num.append(num)
        if len(del_num) > 0:
            for del_index in sorted(del_num, reverse=True):
                del slim_connections[del_index]

        idstore = []
        for connection in slim_connections:
            for lef, rig in connection.items():
                idstore.extend([lef.split(' ')[0], rig.split(' ')[0]])
        slim_main_idnum = dict(Counter(idstore))

        return slim_connections, slim_main_idnum
 
    def get_final_path(self, utr):
        cycl = True
        for connection in self.slim_connections:
            for left_contig, right_contig in connection.items():
                left_contig_id = left_contig.strip().split()[0]
                left_contig_edge = left_contig.strip().split()[1]
                right_contig_id = right_contig.strip().split()[0]
                right_contig_edge = right_contig.strip().split()[1]
                if int(self.begin_seed) == int(left_contig_id) and left_contig_edge == utr:
                    if self.slim_main_idnum[right_contig_id] > 2:
                        cycl = False
                        final_contig_id = self.begin_seed
                        final_contig_edge = utr
                        final_rep_contig = right_contig_id
                        final_rep_edge = right_contig_edge
                    else:
                        final_contig_id = right_contig_id
                        final_contig_edge = right_contig_edge
                    break
                elif int(self.begin_seed) == int(right_contig_id) and right_contig_edge == utr:
                    if self.slim_main_idnum[left_contig_id] > 2:
                        cycl = False
                        final_contig_id = self.begin_seed
                        final_contig_edge = utr
                        final_rep_contig = left_contig_id
                        final_rep_edge = left_contig_edge
                    else:
                        final_contig_id = left_contig_id
                        final_contig_edge = left_contig_edge
                    break

        while cycl:
            if self.slim_main_idnum[final_contig_id] > 2:
                cycl = False
            for connection in self.slim_connections:
                for left_contig, right_contig in connection.items():
                    left_contig_id = left_contig.strip().split()[0]
                    left_contig_edge = left_contig.strip().split()[1]
                    right_contig_id = right_contig.strip().split()[0]
                    right_contig_edge = right_contig.strip().split()[1]
                    if int(final_contig_id) == int(left_contig_id) and left_contig_edge != final_contig_edge:
                        if self.slim_main_idnum[right_contig_id] == 2:
                            final_contig_id = right_contig_id
                            final_contig_edge = right_contig_edge
                        elif self.slim_main_idnum[right_contig_id] > 2:
                            final_rep_contig = right_contig_id
                            final_rep_edge = right_contig_edge
                            cycl = False
                    elif int(final_contig_id) == int(right_contig_id) and right_contig_edge != final_contig_edge: 
                        if self.slim_main_idnum[left_contig_id] == 2:
                            final_contig_id = left_contig_id
                            final_contig_edge = left_contig_edge
                        elif self.slim_main_idnum[left_contig_id] > 2:
                            final_rep_contig = left_contig_id
                            final_rep_edge = left_contig_edge
                            cycl = False

        return final_contig_id, final_rep_contig, final_rep_edge

    def rep_path_order(self):
        rep_path = {}
        for seed, num in self.slim_main_idnum.items():
            if num > 2:
                path_list = []
                for connection in self.slim_connections:
                    for left_contig, right_contig in connection.items():
                        left_contig_id = left_contig.strip().split()[0]
                        left_contig_edge = left_contig.strip().split()[1]
                        right_contig_id = right_contig.strip().split()[0]
                        right_contig_edge = right_contig.strip().split()[1]
                        if int(left_contig_id) == int(seed):
                            path_list.append(right_contig+" "+left_contig_edge)
                            # rep_path[seed] = path_list
                        elif int(right_contig_id) == int(seed):
                            path_list.append(left_contig+" "+right_contig_edge)
                            # rep_path[seed] = path_list
                
                if len(path_list) > 1:
                    path_list = sorted(path_list, key=lambda x: float(self.id_depth[x.split()[0]]), reverse=True)
                    rep_path[seed] = path_list
                    # match_list = []
                    # for mat in path_list:
                    #     match_list.append(mat.split()[0])
                    # if final_contig_id in match_list:
                    #     for ctg in path_list:
                    #         if int(ctg.split()[0]) == int(final_contig_id):
                    #             rep_path[seed].remove(ctg)
                    #             rep_path[seed].append(ctg)
        input_dict = rep_path

        keys = list(input_dict.keys())
        values = list(input_dict.values())

        value_permutations = [list(itertools.permutations(v)) for v in values]

        all_combinations = []
        for permuted_values in itertools.product(*value_permutations):
            permuted_dict = {keys[i]: list(permuted_values[i]) for i in range(len(keys))}
            # for rep_seed, match_id in permuted_dict.items():
            #     for ctg in match_id:
            #         if int(ctg.split()[0]) == int(self.final_contig_id):
            #             permuted_dict[rep_seed].remove(ctg)
            #             permuted_dict[rep_seed].append(ctg)
            all_combinations.append(permuted_dict)
        return all_combinations
    
    def get_final_mt(self, rep_path, begin_utr = "5'"):

        def go_path(dy_ctg, dy_utr, dy_path):
            if dy_ctg not in rep_path.keys():
                found  = False
                for connection in self.slim_connections:
                    for left_contig, right_contig in connection.items():
                        left_contig_id = left_contig.strip().split()[0]
                        left_contig_edge = left_contig.strip().split()[1]
                        right_contig_id = right_contig.strip().split()[0]
                        right_contig_edge = right_contig.strip().split()[1]
                        if int(left_contig_id) == int(dy_ctg) and left_contig_edge != dy_utr:
                            dy_path.append(connection)
                            if right_contig_id in rep_path.keys():
                                        rep_path[right_contig_id].remove(left_contig+' '+right_contig_edge)
                                        rep_path[right_contig_id].append(left_contig+' '+right_contig_edge)
                            dy_ctg = right_contig_id
                            dy_utr = right_contig_edge
                            found = True
                            break
                        elif int(right_contig_id) == int(dy_ctg) and right_contig_edge != dy_utr:
                            dy_path.append(connection)
                            if left_contig_id in rep_path.keys():
                                        rep_path[left_contig_id].remove(right_contig+' '+left_contig_edge)
                                        rep_path[left_contig_id].append(right_contig+' '+left_contig_edge)
                            dy_ctg = left_contig_id
                            dy_utr = left_contig_edge
                            found = True
                            break
                    if found:
                        break

            elif dy_ctg in rep_path.keys():
                if len(rep_path[dy_ctg]) > 0:
                    found = False
                    tep_path = rep_path[dy_ctg]
                    for rep_ctg in rep_path[dy_ctg]:
                        match_id, match_edge, rep_edge = rep_ctg.split()
                        match_ctg = match_id + ' ' + match_edge 
                        for connection in self.slim_connections:
                            for left_contig, right_contig in connection.items():
                                left_contig_id = left_contig.strip().split()[0]
                                left_contig_edge = left_contig.strip().split()[1]
                                right_contig_id = right_contig.strip().split()[0]
                                right_contig_edge = right_contig.strip().split()[1]
                                if int(left_contig_id) == int(dy_ctg) and int(right_contig_id) == int(match_id) and (int(right_contig_id) != int(self.head_contig_id) or (int(self.final_contig_id) == int(self.head_contig_id) and right_contig_edge == "5'")) and left_contig_edge != dy_utr:
                                    dy_path.append(connection)
                                    rep_path[dy_ctg].remove(rep_ctg)
                                    rep_path[dy_ctg].append(rep_ctg)
                                    if right_contig_id in rep_path.keys():
                                        rep_path[right_contig_id].remove(left_contig+' '+right_contig_edge)
                                        rep_path[right_contig_id].append(left_contig+' '+right_contig_edge)
                                    dy_ctg = right_contig_id
                                    dy_utr = right_contig_edge
                                    found = True
                                    break
                                elif int(right_contig_id) == int(dy_ctg) and int(left_contig_id) == int(match_id) and (int(right_contig_id) != int(self.head_contig_id) or (int(self.final_contig_id) == int(self.head_contig_id) and left_contig_edge == "5'")) and right_contig_edge != dy_utr:
                                    dy_path.append(connection)
                                    rep_path[dy_ctg].remove(rep_ctg)
                                    rep_path[dy_ctg].append(rep_ctg)
                                    if left_contig_id in rep_path.keys():
                                        rep_path[left_contig_id].remove(right_contig+' '+left_contig_edge)
                                        rep_path[left_contig_id].append(right_contig+' '+left_contig_edge)
                                    dy_ctg = left_contig_id
                                    dy_utr = left_contig_edge
                                    found = True
                                    break
                            if found:
                                break
                        if found:
                            break

            return dy_ctg, dy_utr, dy_path

        dy_path = []
        final_path = []
        cycl = True
        dy_ctg = self.begin_seed
        dy_utr = begin_utr
        # dy_ctg, dy_utr, dy_path = go_path(dy_ctg, dy_utr, dy_path)
        dy_track_path = rep_path[self.final_rep_contig]
        blur_path = []
        for blur_match in dy_track_path:
            blur_ctg, blur_edge, rep_edge = blur_match.strip().split()
            if rep_edge == self.final_rep_edge:
                blur_path.append(blur_ctg)

        while cycl:
            dy_ctg, dy_utr, dy_path = go_path(dy_ctg, dy_utr, dy_path)
            if int(dy_ctg) == int(self.begin_seed) and dy_utr == begin_utr:
                final_path = dy_path
                cycl = False
            elif int(dy_ctg) == int(self.begin_seed) and dy_utr != begin_utr and len(dy_path) > 2:
                final_path = False
                cycl = False
            else:
                final_path = dy_path

        return final_path, blur_path


    def getloop(self):
        BFS_path_set = set()
        BFS_path_list = []
        all_combinations = self.rep_path_order()
        unrep_combinations_set = set()
        unrep_combinations_list = []
        
        for rep_path in all_combinations:
            if str(rep_path) not in unrep_combinations_set:
                unrep_combinations_set.add(str(rep_path))
                unrep_combinations_list.append(rep_path)
            
        random.seed = 666
        if len(unrep_combinations_list) > 100000:
            sim_combinations = random.sample(unrep_combinations_list, 100000)
        else:
            sim_combinations = unrep_combinations_list
        blur_flag = []
        for rep_path in sim_combinations:
            final_path, blur_path = self.get_final_mt(rep_path=rep_path)
            blur_flag.append([blur_path,final_path])
            if final_path and int(blur_path[0]):
                BFS_path_set.add(str(final_path))
                BFS_path_list.append(final_path)

        if len(blur_flag) > 0:
            blur_path_set = set()
            blur_path_list = []
            blur_path_flag = []
            for fg in blur_flag:
                if int(fg[0][-1]) == int(self.final_contig_id) and fg[1]:
                    blur_path_set.add(str(fg[1]))
                    blur_path_list.append(fg[1])
            for blur_path in blur_path_list:
                if str(blur_path) in blur_path_set:

                    blur_idstore = []
                    for blur_match in blur_path:
                        for lef, rig in blur_match.items():
                            blur_idstore.extend([lef.split(' ')[0], rig.split(' ')[0]])
                    blur_rep_dict = dict(Counter(blur_idstore))
                    blur_rep_seed = [str(seed) for seed, num in blur_rep_dict.items() if num > 2]
                    copy_num = []
                    blur_path_renm = []
                    for blur_match in blur_path:
                        for left_contig, right_contig in blur_match.items():
                            left_contig_id, left_contig_edge = left_contig.strip().split()
                            right_contig_id, right_contig_edge = right_contig.strip().split()
                            if str(left_contig_id) in blur_rep_seed and str(right_contig_id) not in blur_rep_seed:
                                copy_num.append(str(left_contig_id))
                                if copy_num.count(str(left_contig_id)) % 2 == 1:
                                    left_contig_id = f'{left_contig_id}_copy_{(copy_num.count(str(left_contig_id))) // 2 + 1}'
                                else:
                                    left_contig_id = f'{left_contig_id}_copy_{(copy_num.count(str(left_contig_id)) - 1) // 2 + 1}'
                                temp_dict = {}
                                temp_dict[f'{left_contig_id} {left_contig_edge}'] = right_contig
                                blur_path_renm.append(temp_dict)
                            elif str(left_contig_id) not in blur_rep_seed and str(right_contig_id) in blur_rep_seed:
                                copy_num.append(str(right_contig_id))
                                if copy_num.count(str(right_contig_id)) % 2 == 1:
                                    right_contig_id = f'{right_contig_id}_copy_{(copy_num.count(str(right_contig_id))) // 2 +1}'
                                else:
                                    right_contig_id = f'{right_contig_id}_copy_{(copy_num.count(str(right_contig_id)) - 1) // 2 + 1}'
                                temp_dict = {}
                                temp_dict[f'{right_contig_id} {right_contig_edge}'] = left_contig
                                blur_path_renm.append(temp_dict)
                            elif str(left_contig_id) in blur_rep_seed and str(right_contig_id) in blur_rep_seed:
                                copy_num.append(str(left_contig_id))
                                if copy_num.count(str(left_contig_id)) % 2 == 1:
                                    left_contig_id = f'{left_contig_id}_copy_{(copy_num.count(str(left_contig_id))) // 2 + 1}'
                                else:
                                    left_contig_id = f'{left_contig_id}_copy_{(copy_num.count(str(left_contig_id)) - 1) // 2 + 1}'
                                
                                copy_num.append(str(right_contig_id))
                                if copy_num.count(str(right_contig_id)) % 2 == 1:
                                    right_contig_id = f'{right_contig_id}_copy_{(copy_num.count(str(right_contig_id))) // 2 + 1}'
                                else:
                                    right_contig_id = f'{right_contig_id}_copy_{(copy_num.count(str(right_contig_id)) - 1) // 2 + 1}'

                                temp_dict = {}
                                temp_dict[f'{left_contig_id} {left_contig_edge}'] = f'{right_contig_id} {right_contig_edge}'
                                blur_path_renm.append(temp_dict)
                            else:
                                blur_path_renm.append(blur_match)

                    blur_path_flag.append(blur_path_renm)
                    blur_path_set.remove(str(blur_path))
            if len(blur_path_flag) > 1:
                blur_path_flag = sorted(blur_path_flag, key=len, reverse=True)

        if len(BFS_path_set) > 0:
            BFS_path_final = []
            for final_path in BFS_path_list:
                if str(final_path) in BFS_path_set:
                    final_idstore = []
                    for final_match in final_path:
                        for lef, rig in final_match.items():
                            final_idstore.extend([lef.split(' ')[0], rig.split(' ')[0]])
                    BFS_rep_dict = dict(Counter(final_idstore))
                    BFS_rep_seed = [str(seed) for seed, num in BFS_rep_dict.items() if num > 2]
                    copy_num = []
                    BFS_path_renm = []
                    for BFS_match in final_path:
                        for left_contig, right_contig in BFS_match.items():
                            left_contig_id, left_contig_edge = left_contig.strip().split()
                            right_contig_id, right_contig_edge = right_contig.strip().split()
                            if str(left_contig_id) in BFS_rep_seed and str(right_contig_id) not in BFS_rep_seed:
                                copy_num.append(str(left_contig_id))
                                if copy_num.count(str(left_contig_id)) % 2 == 1:
                                    left_contig_id = f'{left_contig_id}_copy_{(copy_num.count(str(left_contig_id))) // 2 + 1}'
                                else:
                                    left_contig_id = f'{left_contig_id}_copy_{(copy_num.count(str(left_contig_id)) - 1) // 2 + 1}'
                                temp_dict = {}
                                temp_dict[f'{left_contig_id} {left_contig_edge}'] = right_contig
                                BFS_path_renm.append(temp_dict)
                            elif str(left_contig_id) not in BFS_rep_seed and str(right_contig_id) in BFS_rep_seed:
                                copy_num.append(str(right_contig_id))
                                if copy_num.count(str(right_contig_id)) % 2 == 1:
                                    right_contig_id = f'{right_contig_id}_copy_{(copy_num.count(str(right_contig_id))) // 2 +1}'
                                else:
                                    right_contig_id = f'{right_contig_id}_copy_{(copy_num.count(str(right_contig_id)) - 1) // 2 + 1}'
                                temp_dict = {}
                                temp_dict[f'{right_contig_id} {right_contig_edge}'] = left_contig
                                BFS_path_renm.append(temp_dict)
                            elif str(left_contig_id) in BFS_rep_seed and str(right_contig_id) in BFS_rep_seed:
                                copy_num.append(str(left_contig_id))
                                if copy_num.count(str(left_contig_id)) % 2 == 1:
                                    left_contig_id = f'{left_contig_id}_copy_{(copy_num.count(str(left_contig_id))) // 2 + 1}'
                                else:
                                    left_contig_id = f'{left_contig_id}_copy_{(copy_num.count(str(left_contig_id)) - 1) // 2 + 1}'
                                
                                copy_num.append(str(right_contig_id))
                                if copy_num.count(str(right_contig_id)) % 2 == 1:
                                    right_contig_id = f'{right_contig_id}_copy_{(copy_num.count(str(right_contig_id))) // 2 + 1}'
                                else:
                                    right_contig_id = f'{right_contig_id}_copy_{(copy_num.count(str(right_contig_id)) - 1) // 2 + 1}'

                                temp_dict = {}
                                temp_dict[f'{left_contig_id} {left_contig_edge}'] = f'{right_contig_id} {right_contig_edge}'
                                BFS_path_renm.append(temp_dict)
                            else:
                                BFS_path_renm.append(BFS_match)

                    BFS_path_final.append(BFS_path_renm)
                    BFS_path_set.remove(str(final_path))
            if len(BFS_path_final) > 1:
                BFS_path_final = sorted(BFS_path_final, key=len, reverse=True)

        if len(blur_path_flag) > 0:
            loop_result = blur_path_flag[:10]
        elif len(BFS_path_final) > 0:
            loop_result = BFS_path_final[:10]
            log.Info("Poor results with automatic ring closure")
        else:
            loop_result = False

        return loop_result
