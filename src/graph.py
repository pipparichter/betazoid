import numpy as np 
import os 
import pandas as pd 
import subprocess
from src.files import FASTAFile, BamFile, FLAGS
import io
from Bio.Seq import Seq
import networkx as nx 


MMSEQS_FIELDS = ['query', 'target', 'alnlen', 'qcov', 'tcov', 'qstart', 'qend', 'tstart', 'tend', 'fident', 'qseq', 'tseq', 'qaln', 'taln', 'qlen', 'tlen']
MMSEQS_FIELDS = '"' + ','.join(MMSEQS_FIELDS) + '"'

MMSEQS_ALIGN_PARAMS = dict()
MMSEQS_ALIGN_PARAMS['--alignment-mode'] = 2 # Semi-global alignment. 
MMSEQS_ALIGN_PARAMS['--min-seq-id'] = 0.95
MMSEQS_ALIGN_PARAMS['--min-aln-len'] = 20 
MMSEQS_ALIGN_PARAMS = ['-a'] + [f'{param} {value}' for param, value in MMSEQS_ALIGN_PARAMS.items()]

MMSEQS_PREFILTER_PARAMS = dict()
MMSEQS_PREFILTER_PARAMS['-k'] = 7 # K-mer size to use for prefiltering 
MMSEQS_PREFILTER_PARAMS['--max-seqs'] = 20 # Controls the maximum number of prefiltering results per query sequence.
MMSEQS_PREFILTER_PARAMS['--mask'] = 0 # Turn off low-complexity matching.
MMSEQS_PREFILTER_PARAMS['--min-ungapped-score'] = 15 # The min. score of an ungapped seed alignment that must exist before a candidate pair is passed to the full alignment stage; roughly equivalent to required number of exactly-matching base pairs.
MMSEQS_PREFILTER_PARAMS = [f'{param} {value}'  for param, value in MMSEQS_PREFILTER_PARAMS.items()]

MMSEQS_CONVERTALIS_PARAMS = ['--format-output', MMSEQS_FIELDS, '--search-type 3']


# https://academic.oup.com/bioinformatics/article/32/9/1323/1744460

def align_reads(path, output_dir:str=None):

    # tmp_dir = os.path.join(output_dir, 'tmp')
    # os.makedirs(tmp_dir, exist_ok=True)

    database_path = os.path.join(output_dir, 'readsDB')
    prefilter_database_path = os.path.join(output_dir, 'readsDB_prefilter')
    aligned_database_path = os.path.join(output_dir, 'readsDB_aligned')
    alignment_path = os.path.join(output_dir, 'alignments.tsv')

    cmds = [' '.join(['mmseqs', 'createdb', path, database_path, '--dbtype', '2'])]
    cmds += [' '.join(['mmseqs', 'prefilter', database_path, database_path, prefilter_database_path] + MMSEQS_PREFILTER_PARAMS)]
    cmds += [' '.join(['mmseqs', 'align', database_path, database_path, prefilter_database_path, aligned_database_path] + MMSEQS_ALIGN_PARAMS)]
    cmds += [' '.join(['mmseqs', 'convertalis', database_path, database_path, aligned_database_path, alignment_path] + MMSEQS_CONVERTALIS_PARAMS)]
    for cmd in cmds:
        print(cmd)
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


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



# Idea is to use knowledge from the seed contig to enforce direction, which should limit the complexity of the problem. 
def add_edge(row, graph:nx.Graph):
    '''Assumes all edges overlap properly, and that all alignments are in the same direction. '''
    if row.qstart > 1: # If alignment is at the right extremus of the query, then direction is query to target. 
        metadata = {'start_a':row.qstart, 'start_b':row.tstart, 'end_a':row.qend, 'end_b':row.tend}
        graph.add_edge(f'{row.target}', f'{row.query}', **metadata) # I might only need to add the beginning to beginning edge because of enforced directionality. 
        # graph.add_edge(f'{row.query}.R', f'{row.target}.R', **metadata)
    else: # If alignment is at the left extremus of the target, then direction is target to query. 
        metadata = {'start_a':row.tstart, 'start_b':row.qstart, 'end_a':row.tend, 'end_b':row.qend}
        graph.add_edge(f'{row.query}', f'{row.target}', **metadata)
        # graph.add_edge(f'{row.target}.R', f'{row.query}.R', **metadata)



# class StringGraph():

#     def __init__(self, df:pd.DataFrame, overlap_length:int=10):

#         reads = get_reads(df)

#         for alignment in align_df.itertuples():




# Representing every read of length n as a collection of n - k + 1 overlapping k-tuples (continuous short strings of fixed length k)
class KmerGraph():
    def _get_sequences(reads_df):

        # Include reads in the orientation they were mapped to the contig. If the read's mate is mapped in one orientation, even if the read
        # itself is unmapped, enforce the opposite orientation
        pass 

    def __init__(self, reads_df:pd.DataFrame, k:int=10):
        self.k = k 
        

# Rules for drawing an edge ar



