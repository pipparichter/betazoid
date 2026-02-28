import numpy as np 
import os 
import pandas as pd 
import subprocess
from src.files import FASTAFile, BamFile, FLAGS
import io
from Bio.Seq import Seq
import networkx as nx 

# https://academic.oup.com/bioinformatics/article/21/suppl_2/ii79/227189?login=true
# https://www.pnas.org/doi/10.1073/pnas.171285098


get_reverse_complement = lambda seq : str(Seq(seq).reverse_complement())

# When writing ambiguous alignments, each alignment becomes a separate SAM record
BBMAP_PARAMS = dict()
BBMAP_PARAMS['minid'] = 0.95 # Minimum percent identity for an alignment to be considered mapped, applied only to the aligned portion of the read. 
BBMAP_PARAMS['idfilter'] = 0.97 # Filters the entire read for identity after alignment, applied over the entire read.
BBMAP_PARAMS['ambiguous'] = 'all' # Ambiguous means all alignments are equally good. 
BBMAP_PARAMS['mappedonly'] = 't' # If true, treats out like outm.
BBMAP_PARAMS['minratio'] = 0.5
BBMAP_PARAMS['editfilter'] = 5 # Consider reads with fewer than 5 indels or mismatches as mapped.
BBMAP_PARAMS['local'] = 'f' # Allow local alignments. Actually, don't.
# BBMAP_PARAMS['unmapped'] = 't' # Include reads that fail the thresholds. 
BBMAP_PARAMS['pigz'] = 't'
BBMAP_PARAMS['unpigz'] = 't'
BBMAP_PARAMS['threads'] = 64
BBMAP_PARAMS['k'] = 13 # K-mer size (13 is default for short-read data)
BBMAP_PARAMS['minhits'] = 1 # Minimum number of K-mer matches to consider a read (1 is default)
BBMAP_PARAMS['nodisk'] = 't' # No temporary files written.
BBMAP_PARAMS['notags'] = 't' # Turn off optional tags.
# BBMAP_PARAMS['secondary'] = 'f' # Keep secondary alignments, which are worse than the primary alignments. 
# BBMAP_PARAMS['pairlenient'] = 't' # Allow unmapped mates to remain. This is not a real option.  
# outm: Write only mapped reads to this file. Includes unmapped paired reads with a mapped mate.

def run_bbmap(ref_path, reads_path_1:str=None, reads_path_2:str=None, output_path:str=None):
    ''''''
    cmd = ['bbmap.sh']
    cmd += [f'in1={reads_path_1}']
    cmd += [f'in2={reads_path_2}']
    cmd += [f'ref={ref_path}']
    cmd += [f'out={output_path}']

    cmd += [f'{param}={str(value)}' for param, value in BBMAP_PARAMS.items()]
    cmd = ' '.join(cmd)

    # cmd += f'| shrinksam | sambam > {output_path}'
    print('run_bbmap:', cmd)
    # Need to pipe output into samtools view so as to not store basically the entire library (which is the typical behavior)?

    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)


# def check_reads(df):
#     '''Confirm that the recruited reads match what I am assuming about them.'''
#     checks = dict()
#     checks['one_mate_mapped'] = np.all((~df.mate_unmapped) | (~df.unmapped))
#     checks['read_strand_assigned'] = np.any(df.strand.isnull())
#     for name, value in checks.items():
#         assert value, f'check_reads: Failed check {name}.'


def recruit_reads(job_name, ref_path:str, n_iters:int=5, output_dir:str='.', reads_path_1:str=None, reads_path_2:str=None):
    output_paths = list()
    # For first iterations, don't care about non-primary alignments. 
    output_paths.append(os.path.join(output_dir, f'{job_name}.reads.{0}.bam'))
    run_bbmap(ref_path, reads_path_1=reads_path_1, reads_path_2=reads_path_2, output_path=output_paths[0])

    for i in range(1, n_iters):
        ref_path_i, n = BamFile.from_file(output_paths[-1]).to_fasta(include_flags=FLAGS['unmapped'] + FLAGS['read_paired'], exclude_flags=FLAGS['mate_unmapped'])
        if (n == 0): # n is the number of sequences written to the FASTA file. 
            print(f'recruit_reads: No additional unmapped reads recruited. Exiting at iteration {i}.')
            break 
        print(f'recruit_reads: Recruiting using {n} unmapped reads at iteration {i}.')
        output_path_i = os.path.join(output_dir, f'{job_name}.reads.{i}.bam')
        if not os.path.exists(output_path_i):
            run_bbmap(ref_path_i, output_path=output_path_i, reads_path_1=reads_path_1, reads_path_2=reads_path_2)
        output_paths.append(output_path_i)

    df = pd.concat([BamFile.from_file(path).to_df() for path in output_paths])
    df['read_number'] = np.select([(df.read_paired & ~df.read_1), (df.read_paired & df.read_2), (~df.read_paired)], ['1', '2', ''], default='')
    df = df.drop_duplicates(['read_id', 'strand', 'read_number'])
    return df


TMP_DIR = '/home/philippar/tmp'

# A handful of possible scenarios... 
# (1) Both mates map in opposite orientations: Include both reads in the specified orientation. 
# (2) One mate maps in one orientation: Include the mapped read and the unmapped read in the opposite orientation. 
# (3) One mate maps in both orientations: Include both members of the pair in both orientations. 
# (4) Both mates map in the same orientation: Include both members of the pair in both orientations. 
def get_reads(df, path:str=os.path.join(TMP_DIR, 'reads.fasta')):
    '''Takes the DataFrame read from a BAM file as input. Outputs a list of sequences to use for building the graph, with IDS encoding the direction and read ID.'''
    # Include reads in the orientation they were mapped to the contig. If the read's mate is mapped in one orientation, even if the read
    # itself is unmapped, enforce the opposite orientation
    id_map = {read_id:i for i, read_id in enumerate(np.sort(df.read_id.unique()))} # Map read IDS to integers to make life easier.
    ids, seqs = list(), list()
    for row in df.itertuples():
        n = row.read_number
        i = id_map[row.read_id]
        assert row.orientation != 'XX', 'get_reads: There should not be any reads where both members of the pair are unmapped.'
        if (row.orientation == 'RF') or (row.orientation == 'FR'):
            ids += [f'{id_map[row.read_id]}.{row.orientation[0]}']
            seqs += [get_reverse_complement(row.seq) if (row.orientation[0] == 'R') else row.seq]
        elif (row.orientation == 'RR') or (row.orientation == 'FF'):
            ids += [f'{id_map[row.read_id]}.F', f'{id_map[row.read_id]}.R']
            seqs += [row.seq, get_reverse_complement(row.seq)]
        elif (row.orientation == 'XF') or (row.orientation == 'XR'):
            orientation = 'F' if (row.orientation[-1] == 'R') else 'F'
            ids += [f'{id_map[row.read_id]}.{orientation}']
            seqs += [row.seq, get_reverse_complement(row.seq)]
        else:
            continue
    print(f'get_reads: Writing {len(ids)} sequences to {path} from {len(id_map)} pairs.')
    fasta_file = FASTAFile()
    fasta_file.ids, fasta_file.seqs = ids, seqs
    fasta_file.write(path)
    return path


MMSEQS_FIELDS = ['query', 'target', 'alnlen', 'qcov', 'tcov', 'qstart', 'qend', 'tstart', 'tend', 'fident', 'qseq', 'tseq', 'qaln', 'taln', 'qlen', 'tlen']
MMSEQS_FIELDS = '"' + ','.join(MMSEQS_FIELDS) + '"'

MMSEQS_ALIGN_PARAMS = dict()
MMSEQS_ALIGN_PARAMS['alignment-mode'] = 2 # Semi-global alignment. 
MMSEQS_ALIGN_PARAMS['min-seq-id'] = 0.95
MMSEQS_ALIGN_PARAMS['min-aln-len'] = 20 
# MMSEQS_ALIGN_PARAMS['a'] = '' # Add backtrace string (convert to alignments with mmseqs convertalis module). Basically allows you to request qaln and taln later with convertalis.
# MMSEQS_ALIGN_PARAMS = ['-a'] + list(np.ravel([[f'--{param}', str(value)] for param, value in MMSEQS_ALIGN_PARAMS.items()]))
MMSEQS_ALIGN_PARAMS = ['-a', '--forward-strand-only'] + list(np.ravel([[f'--{param}', str(value)] for param, value in MMSEQS_ALIGN_PARAMS.items()]))

MMSEQS_PREFILTER_PARAMS = dict()
MMSEQS_PREFILTER_PARAMS['k'] = 7 # K-mer size to use for prefiltering 
MMSEQS_PREFILTER_PARAMS['max-seqs'] = 20 # Controls the maximum number of prefiltering results per query sequence.
MMSEQS_PREFILTER_PARAMS['mask'] = 0 # Turn off low-complexity matching.
MMSEQS_PREFILTER_PARAMS['min-ungapped-score'] = 15 # The min. score of an ungapped seed alignment that must exist before a candidate pair is passed to the full alignment stage; roughly equivalent to required number of exactly-matching base pairs.
MMSEQS_PREFILTER_PARAMS = list(np.ravel([[f'--{param}', str(value)] for param, value in MMSEQS_PREFILTER_PARAMS.items()]))


# https://academic.oup.com/bioinformatics/article/32/9/1323/1744460

def align_reads(path):

    tmp_dir = os.path.join(TMP_DIR, 'tmp')
    database_path = os.path.join(TMP_DIR, 'readsDB')
    prefilter_database_path = os.path.join(TMP_DIR, 'readsDB_prefilter')
    aligned_database_path = os.path.join(TMP_DIR, 'readsDB_aligned')
    alignment_path = os.path.join(TMP_DIR, 'alignments.tsv')

    subprocess.run(['mmseqs', 'createdb', path, database_path, '--dbtype', '2'], check=True)
    subprocess.run(['mmseqs', 'prefilter', database_path, database_path, prefilter_database_path, tmp_dir] + MMSEQS_ALIGN_PARAMS, check=True)
    subprocess.run(['mmseqs', 'align', database_path, database_path, prefilter_database_path, aligned_database_path, tmp_dir] + MMSEQS_PREFILTER_PARAMS, check=True)
    subprocess.run(['mmseqs', 'convertalis', database_path, database_path, aligned_database_path, tmp_dir, alignment_path] + ['--format-output', MMSEQS_FIELDS], check=True)


# Based on paper, I think we want to collapse contained reads into a single node, while still preserving the pair information. 

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



