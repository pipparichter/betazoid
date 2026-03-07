import numpy as np 
import os 
import pandas as pd 
import subprocess
from src.files import FASTAFile, BamFile, FLAGS
import io
from Bio.Seq import Seq
import networkx as nx 
import re

# What are the conditions I want to use for collapsing a node?
# Key thing is that the path it produces must be identical to the path produced through another node.
# In other words, if there is a bubble, collapse the bubble if all generated paths are identical. 

# I actually might not want to collapse all contained alignments in case one of the pairs leads off on a different path (i.e. to a different strain). 
# Instead, want to collapse sub-paths that have the same start and end point. 

# Maybe best thing to do is to collapse the contained alignments, but keep a record of the paired mates.

# Based on paper, I think we want to collapse contained reads into a single node, while still preserving the pair information. 
# I might try just not doing this at first to see what things look like. 
def get_contained_alignments(align_df:pd.DataFrame):
    lengths = dict()
    contained = list()
    for row in align_df.itertuples():
        if (row.tcov == 1) or (row.qcov == 1):
            contained.append((row.target, row.query))
            lengths[row.target], lengths[row.query] = row.tlen, row.qlen 
    graph = nx.Graph()
    graph.add_edges_from(contained)
    groups = dict()
    for group in nx.connected_components(graph):
        # Use the longest sequence in each containment group as the representative of the group. 
        ids = np.array(list(group))
        rep_id = ids[np.argmax([lengths[id_] for id_ in ids])]
        groups[rep_id] = ids 
    return groups 



# Now need to think about how to contruct a graph from the alignments. When traversing the graph, want to build a sequence which does not include the overlap regions,
# so will accrue like (read_a_not_aligned) + (aligned_region) + (read_a_not_aligned)
# Probably want to store a tuple like (qaln/taln, )

# I think my approach will be to build a normal K-mer graph, requiring overlap between reads in order to form a valid path. Also should disallow any paths 
# which connect paired reads in the wrong orientation

# I think the best way to do this is just going to be an mmseqs pairwise alignment between all recruited reads, and then reading in the file. 




class StringGraph():

    # Idea is to use knowledge from the seed contig to enforce direction, which should limit the complexity of the problem. 
    def _add_edge(self, row):
        '''Assumes all edges overlap properly, and that all alignments are in the same direction. '''
        if row.qstart > 1: # If alignment is at the right extremus of the query, then direction is query to target. 
            metadata = {'start_a':row.qstart, 'start_b':row.tstart, 'end_a':row.qend, 'end_b':row.tend}
            self.G.add_edge(f'{row.target}', f'{row.query}', **metadata) # I might only need to add the beginning to beginning edge because of enforced directionality. 
            # graph.add_edge(f'{row.query}.R', f'{row.target}.R', **metadata)
        else: # If alignment is at the left extremus of the target, then direction is target to query. 
            metadata = {'start_a':row.tstart, 'start_b':row.qstart, 'end_a':row.tend, 'end_b':row.qend}
            self.G.add_edge(f'{row.query}', f'{row.target}', **metadata)
            # graph.add_edge(f'{row.target}.R', f'{row.query}.R', **metadata)

    def _get_pair_map(self, reads):

        get_mate_pattern = lambda read : read.split('.')[0] + '.' + ('1' if (read.split('.')[1] == '0') else '0') + '.(R|F)'
        return {read:[r for r in reads if re.match(get_mate_pattern(read), r)] for read in reads}

    def _get_start_nodes(self):
        return [n for n, d in self.G.in_degree() if d == 0]
    
    def _get_end_nodes(self):
        return [n for n, d in self.G.out_degree() if d == 0]

    # Want to initialize with both the alignments and the complete reads list. 
    def __init__(self, reads:list, align_df:pd.DataFrame, overlap_length:int=10):

        self.G = nx.DiGraph()

        # self.G.add_node(row.)
        for row in align_df.itertuples():
            self._add_edge(row)

        print(f'StringGraph.__init__: Initialized a graph with {len(self.G.edges)} edges.')
        print(f'StringGraph.__init__: {len(self._get_start_nodes())} start nodes.')
        print(f'StringGraph.__init__: {len(self._get_end_nodes())} end nodes.')

    # def traverse(self):





# # Representing every read of length n as a collection of n - k + 1 overlapping k-tuples (continuous short strings of fixed length k)
# class KmerGraph():
#     def _get_sequences(reads_df):

#         # Include reads in the orientation they were mapped to the contig. If the read's mate is mapped in one orientation, even if the read
#         # itself is unmapped, enforce the opposite orientation
#         pass 

#     def __init__(self, reads_df:pd.DataFrame, k:int=10):
#         self.k = k 
        

# Rules for drawing an edge ar



