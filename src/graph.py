import numpy as np 
import os 
import pandas as pd 
import subprocess
from src.files import FASTAFile, BamFile, FLAGS
from src.align import align_filter 
import io
from Bio.Seq import Seq
import networkx as nx 
import re
import matplotlib.pyplot as plt


get_read_pair_id = lambda read_id : read_id.split('.')[0]


class StringGraph():

    # Idea is to use knowledge from the seed contig to enforce direction, which should limit the complexity of the problem. 
    def _add_edge(self, row):
        '''Assumes all edges overlap properly, and that all alignments are in the same direction.'''
        if row.qstart > 1: # If alignment is at the right extremus of the query, then direction is query to target. 
            metadata = {'start':row.qstart, 'stop':row.tend, 'score':row.fident * row.alnlen}
            self.G.add_edge(f'{row.target}', f'{row.query}', **metadata) # I might only need to add the beginning to beginning edge because of enforced directionality. 
        else: # If alignment is at the left extremus of the target, then direction is target to query. 
            metadata = {'start':row.tstart, 'stop':row.qend, 'score':row.fident * row.alnlen}
            self.G.add_edge(f'{row.query}', f'{row.target}', **metadata)

    def _get_read_pair(self, read_pair_id:str):
        get_read_id_pattern = lambda read_number: read_pair_id + '.' + read_number + '.(R|F)'
        get_read_ids = lambda read_number : [r for r in self.G.nodes() if re.match(get_read_id_pattern(read_number), r)]
        return get_read_ids('0'), get_read_ids('1')
    
    def _get_valid_mate_id(self, read_id:str):
        read_pair_id = read_id.split('.')[0]
        read_orientation = 'R' if ('F' in read_id) else 'F'
        read_number = '0' if (read_id.split('.')[1] == '1') else '1'
        return f'{read_pair_id}.{read_number}.{read_orientation}'
    
    def _remove_invalid_pairs(self):
        read_ids = list(self.G.nodes()) # These reads will contain the strand and pair information. 
        read_pair_ids = list(set(get_read_pair_id(read_id) for read_id in read_ids))

        keep_read_ids = list()
        for read_pair_id in read_pair_ids:
            read_1_ids, read_2_ids = self._get_read_pair(read_pair_id)

            degrees = [self.G.degree(mate_read_id) for mate_read_id in read_1_ids + read_2_ids]
            best_read_id = np.array(read_1_ids + read_2_ids)[np.argmax(degrees)]
            best_mate_id = self._get_valid_mate_id(best_read_id)

            assert (best_mate_id in read_ids), f'StringGraph._remove_invalid_pairs: Read {best_mate_id} missing from the graph.'
            keep_read_ids += [best_read_id, best_mate_id]
        self.G = self.G.subgraph(self.G).copy()
    
    def _remove_symmetric_edges(self):
        remove_edges = list()
        for a, b in self.G.edges():
            if (b, a) in self.G.edges():
                remove_edges += [min([(a, b), (b, a)], key=lambda e: self.G[e[0]][e[1]]['score'])]
        print(f'StringGraph._remove_symmetric_edges: Removing {(len(remove_edges))} symmetric edges.')
        self.G.remove_edges_from(remove_edges)

    def _remove_contained_reads(self):
        remove_read_ids = [read_id for read_id in self.G.nodes() if self._is_contained(read_id)]
        print(f'StringGraph._remove_contained_reads: Removing {(len(remove_read_ids))} contained reads.')
        self.G.remove_nodes_from(remove_read_ids)
    
    @staticmethod
    def _get_T(reads_df, align_df):

        mutual_contained_read_ids = list()

        T_info = dict()
        is_query_contained = lambda row : (row.qcov >= 0.99) # & (row.tcov < 0.99)
        is_target_contained = lambda row : (row.tcov >= 0.99) #  & (row.qcov < 0.99)

        for row in align_df.itertuples():
            if not (is_query_contained(row) or is_target_contained(row)):
                continue 
            if (is_query_contained(row) and is_target_contained(row)):
                mutual_contained_read_ids += [(row.target, row.query)]

            child_read_id, parent_read_id = (row.target, row.query) if is_query_contained(row) else (row.query, row.target)
            parent_length = len(reads_df.loc[parent_read_id].seq)

            if (child_read_id in T_info) and (T_info[child_read_id][1] < parent_length):
                T_info[child_read_id] = (parent_read_id, parent_length)
            elif (child_read_id not in T_info):
                T_info[child_read_id] = (parent_read_id, parent_length)
            
        T = nx.DiGraph()
        T.add_nodes_from(reads_df.index.values)
        T.add_edges_from([(parent_read_id, child_read_id, {'parent_length':parent_length}) for child_read_id, (parent_read_id, parent_length) in T_info.items()])
        for cycle in nx.simple_cycles(T):
            T.remove_edge(cycle[1], cycle[0])

        assert nx.is_directed_acyclic_graph(T), 'StringGraph._get_T: Containment graph has cycles.'
        # print(f'StringGraph._get_T: Identified {len([read_id for read_id in T.nodes if T.in_degree(read_id) > 0])} contained reads.')
        # print(f'StringGraph._get_T: Identified {len(mutual_contained_read_ids)} mutually-contained reads.')
        return T
    
    def _get_parent(self, read_id:str):
        while self.T.in_degree(read_id) > 0: # Walk up the containment tree. 
            read_id = next(self.T.predecessors(read_id))
        return read_id

    def _is_contained(self, read_id:str):
        return self.T.in_degree(read_id) > 0

    # Want to initialize with both the alignments and the complete read_ids list. 
    def __init__(self, align_path:str=None, reads_path:str=None):

        align_df = align_filter(align_path)
        reads_df = FASTAFile.from_file(reads_path).to_df(parse_description=False)

        self.T = StringGraph._get_T(reads_df, align_df)
        self.G = nx.DiGraph()

        self.G.add_nodes_from([(row.Index, {'length':len(row.seq), 'seq':row.seq}) for row in reads_df.itertuples()])
        for row in align_df.itertuples():
            if self._get_parent(row.query) == self._get_parent(row.target):
                continue # Omit edges between sequences which are contained by the same read
            self._add_edge(row)

        self._remove_invalid_pairs()
        self._remove_symmetric_edges()
        self._remove_contained_reads()

        self.strongly_connected_components = [comp for comp in nx.strongly_connected_components(self.G) if (len(comp) > 1)]
        self.strongly_connected_components = sorted(self.strongly_connected_components, key=lambda comp : len(comp))

        n = 0
        while not nx.is_directed_acyclic_graph(self.G):
            for component in [comp for comp in nx.strongly_connected_components(self.G) if (len(comp) > 1)]:
                g = self.G.subgraph(component)
                edge = min(g.edges, key=lambda e : self.G[e[0]][e[1]]['score']) # Get the cycle edge with the worst score. 
                self.G.remove_edge(*edge)
                n += 1

        print(n)

            # while not nx.is_directed_acyclic_graph(g):
            #     cycle = next(nx.simple_cycles(g))
            #     edges = [(cycle[i], cycle[(i + 1) % len(cycle)]) for i in range(len(cycle))]
            #     remove_edges += [min(edges, key=lambda e : self.G[e[0]][e[1]]['score'])] # Get the cycle edge with the worst score. 
            #     g.remove_edge(*remove_edges[-1])

    # def _is_redundant(self, read_id:str):
    #     parent_read_id = self._get_parent(read_id)
    #     if parent_read_id == read_id:
    #         return False
        
    #     prev_read_ids, next_read_ids = set(self.G.predecessors(read_id)), set(self.G.successors(read_id))
    #     prev_read_ids.discard(parent_read_id)
    #     next_read_ids.discard(parent_read_id)

    #     # parent_prev_read_ids, parent_next_read_ids = self.G.predecessors(parent_read_id), self.G.successors(parent_read_id)
    #     redundant =  all(self.G.has_edge(a, parent_read_id) for a in prev_read_ids)
    #     redundant = redundant and all(self.G.has_edge(parent_read_id, b) for b in next_read_ids)
    #     return redundant

        # print(f'StringGraph.__init__: Initialized a graph with {len(self.G.edges)} edges and {len(self.G.nodes)} nodes.')
        # if not nx.is_directed_acyclic_graph(self.G):
        #     # Strongly connected components are clusters of nodes where every node can reach every other node in the component by following the edges. 
        #     # These can contain multiple cycles. 
        #     G = self.G.copy()
        #     n = 0
        #     while not nx.is_directed_acyclic_graph(G):
        #         for comp in [comp for comp in nx.strongly_connected_components(G) if (len(comp) > 1)]:
        #             g = G.subgraph(comp)
        #             a, b, _ = min(g.edges(data=True), key=lambda e: e[2]['score'])
        #             G.remove_edge(a, b)
        #             n += 1
        #     print(f'StringGraph.__init__: Removed {n} edges.')


    # def _get_node_groups(self): 
    #     G = nx.DiGraph()
    #     G.add_edges_from((a, b) for (a, b, data) in self.G.edges(data=True) if data['contained'])
    #     root_node_ids = [n for n, d in G.in_degree() if d == 0]
    #     groups = {node_id:root_node_id for root_node_id in root_node_ids for node_id in nx.descendants(G, root_node_id) | {root_node_id}}
    #     return groups

    # def collapse(self):
    #     groups = self._get_node_groups()
    #     nodes = [(rep, self.G.nodes[rep]) for rep in set(groups.values())]
    #     edges = set([(groups[a], groups[b], self.G[groups[a]][groups[b]]) for a, b in self.G.edges() if (groups[a] != groups[b])])
    #     # get_overlap_length = lambda data : min(data['stop_a'], data['stop_b']) - max(data['start_a'], data['start_b'])
    #     G = nx.DiGraph()
    #     G.add_edges_from([(a, b, data) for (a, b), data in edges.items()])
    #     G.add_nodes_from(nodes)
    #     self.G = G

        # edges = dict()
        # for a, b, data in self.G.edges(data=True):
        #     # Need to somehow recompute alignments so they are all relative to the contained sequence, 
        #     a_rep, b_rep, data = groups[a], groups[b], data.copy()# By construction, b_rep contains b and a_rep contains a. 
        #     # By graph construction, inner read should always be downstream of the containing read. 
        #     print(a, a_rep)
        #     print(self.G[a_rep][a])
        #     print(self.G[a][a_rep]) # Error is here.
        #     delta_a = (self.G[a_rep][a]['start_a'], -self.G[a_rep][a]['stop_a']) if (a != a_rep) else (0, 0)
        #     delta_b = (self.G[b_rep][b]['start_a'], -self.G[b_rep][b]['stop_a']) if (b != b_rep) else (0, 0)
        #     data['start_a'], data['stop_a'] = data['start_a'] + delta_a[0], data['stop_a'] + delta_a[1] 
        #     data['start_b'], data['stop_b'] = data['start_a'] + delta_b[0], data['stop_b'] + delta_b[1] 

        #     # Always just want to keep the edge corresponding to the longest alignment between nodes. 
        #     if ((a_rep, b_rep) not in edges) or (get_overlap_length(data) > get_overlap_length(edges[(a_rep, b_rep)])):
        #         edges[(a_rep, b_rep)] = data

    # def traverse(self):


        # G_undirected = self.G.to_undirected()
        # for _, mate_read_ids in self._get_pair_map().items():
        #     assert len(mate_read_ids) <= 2, 'Graph.__init__: There should not be more than two nodes per read in the graph.'
        #     print(mate_read_ids)
        #     degrees = [self.G.degree(mate_read_id) for mate_read_id in mate_read_ids]
        #     if (len(mate_read_ids) < 2) or (degrees[0] == degrees[-1]):
        #         continue 
        #     print(_, degrees)
        #     remove_nodes = nx.node_connected_component(G_undirected, np.array(mate_read_ids)[np.argmin(degrees)])
        #     self.G.remove_nodes_from(remove_nodes)


        # print(f'StringGraph._get_read_groups: Found {len(edges)} edges in the graph corresponding to contained alignments.')
        # # print(f'StringGraph._get_read_groups: Contained alignments produce {len(list(nx.connected_components(G)))} node groups.')
        # print(f'StringGraph._get_read_groups: Created {n_groups} groups from {len(self.G.nodes())} nodes.')


            # assert (len(read_1_ids) <= 2) and (len(read_2_ids) <= 2), 'Graph.__init__: There should not be more than two nodes per read in the graph.'

# # Representing every read_id of length n as a collection of n - k + 1 overlapping k-tuples (continuous short strings of fixed length k)
# class KmerGraph():
#     def _get_sequences(read_ids_df):

#         # Include read_ids in the orientation they were mapped to the contig. If the read_id's mate is mapped in one orientation, even if the read_id
#         # itself is unmapped, enforce the opposite orientation
#         pass 

#     def __init__(self, read_ids_df:pd.DataFrame, k:int=10):
#         self.k = k 

    # def _get_orphan_nodes(self):
    #     return [n for n, d in self.G.degree() if (d == 0)]
    
    # def _get_start_nodes(self):
    #     return [n for n, d in self.G.in_degree() if (d == 0)]
    
    # def _get_end_nodes(self):
    #     return [n for n, d in self.G.out_degree() if (d == 0)]


# def _get_node_groups(self, align_df:pd.DataFrame):
#     align_df = align_df[(align_df['query'].isin(self.G.nodes()))]

#     lengths = {row.query:row.qlen for row in align_df.itertuples()}
#     lengths.update({row.query:row.qlen for row in align_df.itertuples()})
#     G = nx.Graph()
#     G.add_edges_from([(row.target, row.query)] for row in align_df.itertuples() if ((row.tcov == 1) or (row.qcov == 1)))
#     groups = dict()
#     for group in nx.connected_components(G):
#         # Use the longest sequence in each containment group as the representative of the group. 
#         ids = np.array(list(group))
#         rep_id = ids[np.argmax([lengths[id_] for id_ in ids])]
#         groups[rep_id] = ids 
#     return groups 
