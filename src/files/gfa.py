import networkx as nx 
import itertools
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
from Bio.Seq import Seq
import re

# https://www.pnas.org/doi/10.1073/pnas.171285098

reverse_complement = lambda seq : str(Seq(seq).reverse_complement())


def merge(seqs):
    seq = seqs[0]
    for seq_ in seqs[1:]:
        max_overlap_length = min(len(seq), len(seq_))
        overlaps = [seq_[:i + 1] for i in range(max_overlap_length)]
        overlaps = [overlap for overlap in overlaps if (overlap == seq[-len(overlap):])] # Filter for matching overlaps. 
        assert len(overlaps) > 0, f'merge: No overlap found between input sequences.'
        best_overlap = max(overlaps, key=len)
        seq = seq + seq_.replace(best_overlap, '')
    return seq 


class GFAFile():
    def __init__(self):
        pass 

    @staticmethod
    def _parse_node(line):
        fields = ['node_id', 'seq']
        line = line.strip().split()[1:3]
        return dict(zip(fields, line))
    
    @staticmethod
    def _parse_edge(line):
        fields = ['from_node_id', 'from_orientation', 'to_node_id', 'to_orientation']
        line = line.strip().split()[1:5]
        return dict(zip(fields, line))
    
    def _get_nodes_with_self_edges(self):
        node_ids = [node_id for node_id in self.graph.nodes if (self.graph.has_edge(node_id, node_id))]
        print(f'GFAFile._get_nodes_with_self_edges: Found {len(node_ids)} nodes with self-edges.')
        return node_ids
    
    def _remove_disconnected_nodes_with_self_edges(self):
        n = 0
        for node_id in self._get_nodes_with_self_edges():
            n_self_edges = self.graph.number_of_edges(node_id, node_id)
            if (self.graph.in_degree(node_id) == n_self_edges) and (self.graph.out_degree(node_id) == n_self_edges):
                self.graph.remove_node(node_id)
                n += 1
        print(f'GFAFile._remove_disconnected_nodes_with_self_edges: Removed {n} singleton nodes with self-edges.')
        
    @classmethod
    def from_file(cls,  path):

        file = cls()
        graph = nx.DiGraph()

        with open(path, 'r') as f:
            for line in f.readlines():

                if line.startswith('S'):
                    node = GFAFile._parse_node(line)
                    graph.add_node(node['node_id'] + '+', orientation='+', seq=node['seq'])
                    graph.add_node(node['node_id'] + '-', orientation='-', seq=reverse_complement(node['seq']))

                flip_orientation = lambda orientation : '+' if (orientation == '-') else '-'
                if line.startswith('L'):
                    edge = GFAFile._parse_edge(line)
                    # if edge['from_node_id'] == edge['to_node_id']:
                    #     continue
                    graph.add_edge(edge['from_node_id'] + edge['from_orientation'], edge['to_node_id'] + edge['to_orientation'])
                    graph.add_edge(edge['to_node_id'] + flip_orientation(edge['to_orientation']), edge['from_node_id'] + flip_orientation(edge['from_orientation']))

        file.graph = graph 
        file._remove_disconnected_nodes_with_self_edges()
        return file 
    
    def __len__(self):
        return len(self.graph.nodes)
    
    def _has_path(self, from_node_id=None, to_node_id=None, error=True):
        endpoints = list()
        for orientations in [('+', '+'), ('+', '-'), ('-', '-'), ('-', '+')]:
            try:
                paths = nx.algorithms.shortest_simple_paths(self.graph, f'{from_node_id}{orientations[0]}', f'{to_node_id}{orientations[1]}')
                next(paths)
            except:
                continue
            # print(f'GFAFile._has_path: At least one path found between {from_node_id}{orientations[0]} and {to_node_id}{orientations[1]}')
            endpoints += [(f'{from_node_id}{orientations[0]}', f'{to_node_id}{orientations[1]}')]

        if error:
            assert len(endpoints) > 0, f'GFAFile._has_path: No valid path found between nodes {from_node_id} and {to_node_id}'
        # assert (len(endpoints) == 2) or (len(endpoints) == 0), f'GFAFile._has_path: Expected 2 sets of endpoints, one for each orientation, but found {len(endpoints)} between nodes {from_node_id} and {to_node_id}.'
        return endpoints
    
    def get_node(self, node_id):
        return self.graph.nodes[node_id]['seq']
    
    def remove_edge(self, *node_ids):
        pattern = '|'.join([str(node_id) for node_id in node_ids])
        pattern = rf'({pattern})(\+|\-)'
        edges = self.graph.copy().edges
        for from_node_id, to_node_id in edges:
            if re.match(pattern, from_node_id) and re.match(pattern, to_node_id):
                self.graph.remove_edge(from_node_id, to_node_id)


    def _get_seq_from_path(self, path:list):
        return merge([self.graph.nodes[node_id]['seq'] for node_id in path])

    # By construction, the number of source nodes and terminal nodes should be equal. 

    def _get_singleton_nodes(self):
        return [node_id for node_id in self.graph.nodes if (self.graph.in_degree(node_id) == 0) and (self.graph.out_degree(node_id) == 0)]

    def _get_source_nodes(self):
        return [node_id for node_id in self.graph.nodes if (self.graph.in_degree(node_id) == 0)]
    
    def _get_terminal_nodes(self):
        return [node_id for node_id in self.graph.nodes if (self.graph.out_degree(node_id) == 0)]

    def _get_paths(self, from_node_id:str, to_node_id:str, n_paths:int=np.inf):
        paths, seqs, i = list(), list(), 0
        for path in nx.algorithms.shortest_simple_paths(self.graph, from_node_id, to_node_id):
            seqs += [self._get_seq_from_path(path)]
            paths += [path]
            i += 1
            if i >= n_paths:
                break
        return seqs, paths
    
    def _has_path(self, from_node_id, to_node_id):
        try:
            self._get_paths(from_node_id, to_node_id, n_paths=1)
        except:
            return False
        return True
    
    @staticmethod
    def _get_unique_endpoints(endpoints):
        get_node_number = lambda node_id : int(re.search(r'\d+', node_id).group(0))
        is_equal = lambda e, e_: (get_node_number(e_[0]) == get_node_number(e[1])) and (get_node_number(e_[1]) == get_node_number(e[0]))

        unique_endpoints = list()
        for e in endpoints:
            if not any(is_equal(e, e_) for e_ in unique_endpoints):
                unique_endpoints.append(e)
        return unique_endpoints
    
    @staticmethod
    def _get_unique_nodes(node_ids):
        get_node_number = lambda node_id : int(re.search(r'\d+', node_id).group(0))

        unique_node_ids = list()
        for n in node_ids:
            if not any(get_node_number(n) == get_node_number(n_) for n_ in unique_node_ids):
                unique_node_ids.append(n)
        return unique_node_ids

    
    def get_all_paths(self):
        seqs, paths = list(), list()
        source_node_ids = self._get_source_nodes()
        terminal_node_ids = self._get_terminal_nodes()
        singleton_node_ids = self._get_unique_nodes(self._get_singleton_nodes())

        seqs += [self.get_node(node_id) for node_id in singleton_node_ids]
        paths += [[node_id] for node_id in singleton_node_ids]

        endpoints = [(from_node_id, to_node_id) for (from_node_id, to_node_id) in itertools.product(source_node_ids, terminal_node_ids)]
        endpoints = [(from_node_id, to_node_id) for from_node_id, to_node_id in endpoints if ((from_node_id != to_node_id) and self._has_path(from_node_id, to_node_id))]
        endpoints = GFAFile._get_unique_endpoints(endpoints)

        print(f'GFAFile: Found {len(endpoints)} valid sets of termini.')
        for (from_node_id, to_node_id) in endpoints:
            seqs_, paths_ = self._get_paths(from_node_id, to_node_id)
            seqs += seqs_
            paths += paths_

        print(f'GFAFile: Found {len(paths)} unique paths.')
        return seqs, paths

    

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