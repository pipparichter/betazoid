import os 
import subprocess
import pandas as pd 
import re
import numpy as np 


# About semi-global alignment: https://academic.oup.com/bioinformatics/article/32/9/1323/1744460

MMSEQS_FIELDS = ['query', 'target', 'alnlen', 'qcov', 'tcov', 'qstart', 'qend', 'tstart', 'tend', 'fident', 'qseq', 'tseq', 'qaln', 'taln', 'qlen', 'tlen']

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

MMSEQS_CONVERTALIS_PARAMS = ['--format-output', '"' + ','.join(MMSEQS_FIELDS) + '"', '--search-type 3']



def align(path, output_dir:str=None):

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


id_pattern = f'(\d+.\d).(R|F|X)'
get_id_no_orientation = lambda id_ : re.match(id_pattern, id_).group(1) # Just removing the orientation from the ID.

def align_load(path:str):
    '''Load the output of MMseqs alignment.'''
    df = pd.read_csv(path, header=None, names=MMSEQS_FIELDS, sep='\t')
    # IDs are of the form {read_pair_integer}.{read_number}.{orientation}
    df['query_id_no_orientation'] = df.query.apply(get_id_no_orientation)
    df['target_id_no_orientation'] = df.target.apply(get_id_no_orientation)
    return df


# is_contained_alignment = lambda df : 

def align_filter(path, max_distance_from_extrema:int=3):

    df = align_load(df)
    print(f'align_filter: Loaded {len(df)} alignments from {path}')
    # Based on how the mapped reads were extracted from the BAM file (both forward and reverse versions were extracted independently
    # depending on how they mapped to the scaffold), we only want alignments where both are on the forward strand.

    filters = dict()
    filters['opposite_strand_alignments'] = (df.qstart > df.qend) | (df.tstart > df.tend)
    filters['self_alignments'] = df.target_id_no_orientation == df.query_id_no_orientation
    filters['non_extreme_endpoints'] = ((df.qstart <= max_distance_from_extrema) + ((df.tlen - df.tend) <= max_distance_from_extrema) + ((df.qlen - df.qend) <= max_distance_from_extrema) + (df.start <= max_distance_from_extrema)) < 2

    masks = np.ones(len(df))
    for name, mask in filter.items():
        print(f'align_filter: {mask.sum()} alignments failing filter {name}')
        masks = masks & ~mask 
    
    df = df[masks].copy()
    print(f'align_filter: {len(df)} alignments remaining after filtering')


