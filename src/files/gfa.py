
import networkx as nx 
import itertools
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
from Bio.Seq import Seq
import re
from tqdm import tqdm 
import copy

# https://www.pnas.org/doi/10.1073/pnas.171285098

reverse_complement = lambda seq : str(Seq(seq).reverse_complement())
get_node_number = lambda node_id : re.search(r'\d+', node_id).group(0)


# def get_unique_nodes(node_ids):
#     node_ids = sorted(node_ids)

#     unique_node_ids = list()
#     node_numbers = list()
#     for node_id in node_ids:
#         if get_node_number(node_id) in node_numbers:
#             continue 
#         else:
#             unique_node_ids.append(node_id)
#             node_numbers.append(get_node_number(node_id))
#     return unique_node_ids



# def merge(seqs):
#     seq = seqs[0]
#     for seq_ in seqs[1:]:
#         max_overlap_length = min(len(seq), len(seq_))
#         overlaps = [seq_[:i + 1] for i in range(max_overlap_length)]
#         overlaps = [overlap for overlap in overlaps if (overlap == seq[-len(overlap):])] # Filter for matching overlaps. 
#         assert len(overlaps) > 0, f'merge: No overlap found between input sequences.'
#         best_overlap = max(overlaps, key=len)
#         seq = seq + seq_.replace(best_overlap, '')
#     return seq 


class GFAFile():
    def __init__(self):
        pass 

    @staticmethod
    def _parse_node(line):
        depth = float(re.search(r'DP:f:([^\s]+)', line).group(1))
        fields = ['node_id', 'seq']
        line = line.strip().split()[1:3]
        
        node = dict(zip(fields, line))
        node['depth'] = depth 
        return node 
    
    @staticmethod
    def _parse_edge(line):
        overlap = int(re.search(r'(\d+)M', line).group(1))
        fields = ['from_node_id', 'from_orientation', 'to_node_id', 'to_orientation']
        line = line.strip().split()[1:5]
        edge = dict(zip(fields, line))
        edge['overlap'] = overlap
        return edge
    
    def _get_nodes_with_self_edges(self):
        node_ids = [node_id for node_id in self.graph.nodes if (self.graph.has_edge(node_id, node_id))]
        print(f'GFAFile._get_nodes_with_self_edges: Found {len(node_ids)} nodes with self-edges.')
        return node_ids
    
    def _remove_cycles(self):
        removed_node_ids = list()
        while not nx.is_directed_acyclic_graph(self.graph):
            cycle = next(nx.simple_cycles(self.graph))
            remove_node_id = min([node_id for node_id in cycle], key=lambda node_id : self.graph.nodes[node_id]['depth'])
            self.graph.remove_node(remove_node_id)
            removed_node_ids.append(remove_node_id)

        print(f'GFAFile._remove_cycles: Removed {len(removed_node_ids)} nodes from the graph: {removed_node_ids}')
    
    def _remove_singleton_nodes_with_self_edges(self):
        n = 0
        for node_id in self._get_nodes_with_self_edges():
            n_self_edges = self.graph.number_of_edges(node_id, node_id)
            if (self.graph.in_degree(node_id) == n_self_edges) and (self.graph.out_degree(node_id) == n_self_edges):
                self.graph.remove_node(node_id)
                n += 1
        print(f'GFAFile._remove_singleton_nodes_with_self_edges: Removed {n} singleton nodes with self-edges.')
        
    @classmethod
    def from_file(cls, path):

        file = cls()
        graph = nx.DiGraph()

        with open(path, 'r') as f:
            for line in f.readlines():

                if line.startswith('S'):
                    node = GFAFile._parse_node(line)
                    graph.add_node(node['node_id'] + '+', orientation='+', seq=node['seq'], depth=node['depth'])
                    graph.add_node(node['node_id'] + '-', orientation='-', seq=reverse_complement(node['seq']), depth=node['depth'])

                flip_orientation = lambda orientation : '+' if (orientation == '-') else '-'
                if line.startswith('L'):
                    edge = GFAFile._parse_edge(line)
                    # if edge['from_node_id'] == edge['to_node_id']:
                    #     continue
                    graph.add_edge(edge['from_node_id'] + edge['from_orientation'], edge['to_node_id'] + edge['to_orientation'], overlap=edge['overlap'])
                    graph.add_edge(edge['to_node_id'] + flip_orientation(edge['to_orientation']), edge['from_node_id'] + flip_orientation(edge['from_orientation']), overlap=edge['overlap'])

        file.graph = graph 
        file._remove_singleton_nodes_with_self_edges()
        file._remove_cycles()
        return file 
    
    def __len__(self):
        return len(self.graph.nodes)
    
    def remove_edge(self, *node_ids):
        pattern = '|'.join([str(node_id) for node_id in node_ids])
        pattern = rf'({pattern})(\+|\-)'
        edges = self.graph.copy().edges
        for from_node_id, to_node_id in edges:
            if re.match(pattern, from_node_id) and re.match(pattern, to_node_id):
                self.graph.remove_edge(from_node_id, to_node_id)

    def remove_node(self, node_id):
        pattern = rf'({node_id})(\+|\-)' 
        node_ids = self.graph.copy().nodes
        for node_id in node_ids:
            if re.match(pattern, node_id):
                self.graph.remove_node(node_id)

    def get_node_data(self, node_id):
        pattern = fr'{node_id}([\+\-])'
        for node_id, data in self.graph.nodes(data=True):
            if re.match(pattern, node_id):
                return data 

    def get_node_ids(self, unique:bool=True):

        node_ids = self.graph.copy().nodes()
        if unique:
            node_ids = sorted(list(set([get_node_number(node_id) for node_id in node_ids])))
        return sorted(node_ids)

    def to_fasta(self, path:str, min_length:int=150, unique:bool=True):   
        
        nodes = dict()
        for node_id in self.get_node_ids(unique=unique):
            nodes[node_id] = self.graph.nodes[f'{node_id}+']['seq']
        nodes = {id_:seq for id_, seq in nodes.items() if (len(seq) > min_length)}
        print(f'GFAFile.to_fasta: {len(nodes)} nodes exceeding length {min_length}.')

        file = open(path, 'w')
        for id_, seq in nodes.items():
            file.write(f'>{id_}\n')
            file.write('\n'.join([seq[j:j + 80] for j in range(0, len(seq), 80)]))
            file.write('\n')
        file.close()
        return path 


    def get_seq_from_path(self, path:list):
        seq = [self.graph.nodes[path[0]]['seq']]
        for i in range(len(path) - 1):
            from_node_id, to_node_id = path[i], path[i + 1]
            overlap = self.graph.edges[from_node_id, to_node_id]['overlap']
            seq.append(self.graph.nodes[to_node_id]['seq'][overlap:])
        return ''.join(seq)
    
    # By construction, the number of source nodes and terminal nodes should be equal. 

    def _get_source_nodes(self):
        return [node_id for node_id in self.graph.nodes if (self.graph.in_degree(node_id) == 0)]

    def get_paths(self, node_id:str):
        assert not not nx.is_directed_acyclic_graph(self.graph), f'GFAFile.get_paths: Graph has cycles, so cannot use depth-first search.'

        if self.graph.out_degree(node_id) == 0:
            return [[node_id]]
        
        paths = list()
        for child_node_id in self.graph.successors(node_id):
            child_paths = self.get_paths(child_node_id) # Returns a list of lists. 
            for child_path in child_paths:
                paths.append([node_id] + child_path)

        return paths 
    
    def prune(self, min_depth=5, min_length:int=300):
        '''Remove low-support leaf nodes from the graph.'''
        n = 0 
        remove = lambda node_id : (self.graph.nodes[node_id]['depth'] < min_depth) and (len(self.graph.nodes[node_id]['seq']) < min_length)
        is_leaf_node = lambda node_id : (self.graph.out_degree(node_id) == 0) or (self.graph.in_degree(node_id) == 0)

        while True:
            leaf_node_ids = [node_id for node_id in self.graph.nodes if is_leaf_node(node_id)]
            leaf_node_ids = [node_id for node_id in leaf_node_ids if remove(node_id)]
            if len(leaf_node_ids) == 0:
                break 
            self.graph.remove_nodes_from(leaf_node_ids)
            n += len(leaf_node_ids)

        print(f'GFAFile.prune: Removed {n} terminal nodes shorter than {min_length} with depth less than {min_depth}')

    def get_all_paths(self, path:str=None, min_length:int=0):

        source_node_ids = self._get_source_nodes() # This includes singleton nodes. 
        print(f'GFAFile: Found {len(source_node_ids)} source nodes.')

        longest_path = 0
        i = 0

        file = open(path, 'w') 

        pbar = tqdm(source_node_ids, desc='GFAFile.get_all_paths')
        for source_node_id in pbar:
            pbar.set_description(f'GFAFile.get_all_paths: Computing all paths starting from {source_node_id}. Longest path is {longest_path} base pairs.')
            
            paths_ = self.get_paths(source_node_id)
            seqs_ = [self.get_seq_from_path(path_) for path_ in paths_]

            for path_, seq in zip(paths_, seqs_):

                if len(seq) < min_length:
                    continue 
                
                longest_path = max(longest_path, len(seq))

                file.write(f'>path_{i + 1} ' + ','.join(path_) + '\n')
                file.write('\n'.join([seq[j:j + 80] for j in range(0, len(seq), 80)]))
                file.write('\n')
                i += 1

        file.close()
        
        print(f'GFAFile.get_all_paths: Wrote {i} unique paths to {path}.')


    def draw(self, path=None):

        fig, ax = plt.subplots()

        graph = self.graph.copy()
        graph.remove_nodes_from(list(nx.isolates(graph)))
        if path is not None:
            graph = graph.subgraph(path)

        layout = nx.spring_layout(graph)
        nx.draw_networkx(graph, pos=layout, node_size=2)
        ax.set_axis_off()
        plt.show()

    def get_n_cycles(self):
        return len(list(nx.simple_cycles(self.graph))) 
    

    
    # @staticmethod
    # def _get_unique_endpoints(endpoints):
    #     get_node_number = lambda node_id : int(re.search(r'\d+', node_id).group(0))
    #     is_equal = lambda e, e_: (get_node_number(e_[0]) == get_node_number(e[1])) and (get_node_number(e_[1]) == get_node_number(e[0]))

    #     unique_endpoints = list()
    #     for e in endpoints:
    #         if not any(is_equal(e, e_) for e_ in unique_endpoints):
    #             unique_endpoints.append(e)
    #     return unique_endpoints
    
    # @staticmethod
    # def _get_unique_nodes(node_ids):
    #     get_node_number = lambda node_id : int(re.search(r'\d+', node_id).group(0))

    #     unique_node_ids = list()
    #     for n in node_ids:
    #         if not any(get_node_number(n) == get_node_number(n_) for n_ in unique_node_ids):
    #             unique_node_ids.append(n)
    #     return unique_node_ids

    
    # def get_paths(self, from_node_id, to_node_id, n_paths:int=np.inf):    
    #     i, seqs, paths = 0, list(), list()

    #     from_node_id, to_node_id = self._has_path(from_node_id=from_node_id, to_node_id=to_node_id, error=True)[0]

    #     for path in nx.algorithms.shortest_simple_paths(self.graph, from_node_id, to_node_id):
    #         seq = None
    #         for next_seq in [self.graph.nodes[node_id]['seq'] for node_id in path]:
    #             seq = next_seq if (seq is None) else GFAFile._merge_seqs(seq, next_seq)

    #         seqs += [seq]
    #         paths += [path]
    #         i += 1
    #         if i >= n_paths:
    #             break

    #     print(f'GFAFile.get_paths: Found {len(paths)} paths between {from_node_id} and {to_node_id}.')
    #     return seqs, paths


    # def get_all_paths(self):
    #     seqs, paths = list(), list()
    #     source_node_ids = self._get_source_nodes()
    #     terminal_node_ids = self._get_terminal_nodes()
    #     singleton_node_ids = self._get_unique_nodes(self._get_singleton_nodes())

    #     seqs += [self.get_node(node_id) for node_id in singleton_node_ids]
    #     paths += [[node_id] for node_id in singleton_node_ids]

    #     endpoints = [(from_node_id, to_node_id) for (from_node_id, to_node_id) in itertools.product(source_node_ids, terminal_node_ids)]
    #     endpoints = [(from_node_id, to_node_id) for from_node_id, to_node_id in endpoints if ((from_node_id != to_node_id) and self._has_path(from_node_id, to_node_id))]
    #     endpoints = GFAFile._get_unique_endpoints(endpoints)

    #     print(f'GFAFile: Found {len(endpoints)} valid sets of termini.')
    #     pbar = tqdm(endpoints, desc='get_all_paths')
    #     for (from_node_id, to_node_id) in pbar:
    #         pbar.set_description(f'get_all_paths: Computing all paths between {from_node_id} and {to_node_id}. Current number of paths is {len(paths)}.')
    #         seqs_, paths_ = self._get_paths(from_node_id, to_node_id)
    #         seqs += seqs_
    #         paths += paths_

    #     print(f'GFAFile: Found {len(paths)} unique paths.')
    #     return seqs, paths


    # def _get_paths(self, from_node_id:str, to_node_id:str, n_paths:int=np.inf):
    #     paths, seqs, i = list(), list(), 0
    #     for path in nx.algorithms.shortest_simple_paths(self.graph, from_node_id, to_node_id):
    #         # seqs += [self.get_seq_from_path(path)]
    #         paths += [path]
    #         i += 1
    #         if i >= n_paths:
    #             break
    #     return seqs, paths

    # def _has_path(self, from_node_id=None, to_node_id=None, error=True):
    #     endpoints = list()
    #     for orientations in [('+', '+'), ('+', '-'), ('-', '-'), ('-', '+')]:
    #         try:
    #             paths = nx.algorithms.shortest_simple_paths(self.graph, f'{from_node_id}{orientations[0]}', f'{to_node_id}{orientations[1]}')
    #             next(paths)
    #         except:
    #             continue
    #         # print(f'GFAFile._has_path: At least one path found between {from_node_id}{orientations[0]} and {to_node_id}{orientations[1]}')
    #         endpoints += [(f'{from_node_id}{orientations[0]}', f'{to_node_id}{orientations[1]}')]

    #     if error:
    #         assert len(endpoints) > 0, f'GFAFile._has_path: No valid path found between nodes {from_node_id} and {to_node_id}'
    #     # assert (len(endpoints) == 2) or (len(endpoints) == 0), f'GFAFile._has_path: Expected 2 sets of endpoints, one for each orientation, but found {len(endpoints)} between nodes {from_node_id} and {to_node_id}.'
    #     return endpoints