import sys 
sys.path.append('/home/prichter/Documents/banfield/betazoid/src/files/')
sys.path.append('/home/prichter/Documents/banfield/betazoid/src/')

from fasta import FASTAFile
from gfa import * 
from bam import BamFile

import os 
import re 
import pandas as pd 
import glob
import numpy as np 
import seaborn as sns 
from scipy.stats import gmean 
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LinearSegmentedColormap, to_hex
from Bio.Seq import Seq

PROJECT_IDS = pd.read_csv('project_ids.csv', index_col=0).project_id.to_dict()
FORWARD_READS_PATHS = pd.read_csv('reads_paths.csv', index_col=0).forward_reads_path.to_dict()
REVERSE_READS_PATHS = pd.read_csv('reads_paths.csv', index_col=0).reverse_reads_path.to_dict()

START_CODONS = ['ATG', 'GTG', 'TTG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']

reverse_complement = lambda seq : str(Seq(seq).reverse_complement())
# get_gene_id = lambda path : os.path.basename(path).replace('_pred.txt', '').replace('genes_', '') # Remove the prefix and file extension. 
get_gene_id = lambda string : re.search(r'orfm.bz_\d+\.\d+_\d+', string).group(0)
get_genome_id = lambda string : re.search(r'bz_\d+', string).group(0)

def apply_filters(filters:dict, df:pd.DataFrame, return_idxs=False):
    keep = None
    for name, filter_ in filters.items(): # mask should correspond to things failing the test. 
        keep = ~filter_ if (keep is None) else (keep & ~filter_)
        print(f'apply_filters: {filter_.sum()} entries removed by {name}.')

    if return_idxs:
        return df[keep].copy(), df[~keep].index.values  
    else:
        return df[keep].copy()
    

REDUCED_ALPHABET = {'A':'A', 'V':'A', 'L':'A', 'I':'A', 'M':'A', 'F':'R', 'W':'R','Y':'R','K':'+','R':'+','H':'+','D':'-','E':'-','S':'P','T':'P','N':'P','Q':'P','G':'G','P':'P','C':'C', '.':'.'}

# ! muscle -align ../data/genes/cluster_1.fa -output ../data/genes/cluster_1.afa
class MSA():
    gap_symbol = '.'
    def __init__(self, path:str='../data/genes/cluster_1.afa'):
        self.df = FASTAFile.from_file(path).to_df()
        self.df['seq'] = self.df.seq.str.replace('-', MSA.gap_symbol)
        self.arr = np.array([list(seq) for seq in self.df.seq])
        self.seqs = self.df.seq.values

        self.n_cols = len(self.seqs[0])
    
    def to_arr(self, alphabet:dict=None):
        '''Convert the MSA loaded from the FASTA file into a two-dimensional numpy array, where each entry is a single residue.'''
        return np.vectorize(alphabet.get)(self.arr.copy()) if (alphabet is not None) else self.arr.copy()
    
    def to_df(self, alphabet:dict=None):
        df = self.df.copy()
        df['seq'] = df.seq.str.replace(alphabet) if (alphabet is not None) else df['seq']
        return df        
    
    def map_idx_from_msa(self, idx, gene_id:str):
        '''Convert the index of a residue in the MSA to the index of a residue in one of the aligned sequences.'''
        seq = self.df.loc[gene_id].seq
        assert seq[idx] != MSA.gap_symbol, f'get_idx: The input index corresponds to a gap in the aligned {gene_id}.'
        n_gaps = seq[:idx].count('-') # Get the number of gaps which occur before the requested index.
        return idx - n_gaps
    
    def map_idx_to_msa(self, idx, gene_id:str):
        '''Convert the index of a residue in one of the sequences to an index in the MSA.'''
        n, n_gaps = 0, 0
        for aa in self.df.loc[gene_id].seq:
            if (n == idx) and (aa != MSA.gap_symbol):
                break # So that it doesn't exit if idx = 0 and the first symbol is a gap. 
            n += int(aa != MSA.gap_symbol)
            n_gaps += int(aa == MSA.gap_symbol)
        return n + n_gaps
    
    def show(self, start:int=None, stop:int=None):
        string = ''
        start = 0 if (start is None) else start
        stop = self.n_cols if (stop is None) else stop
        for row in self.df.sort_index().itertuples():
            seq = row.seq[start:stop]
            string += f'{start}\t{seq}\t{stop}\t{row.Index}\n'
        print(string)
    

# def download_ncbi(ids, db='protein', output_dir='.'):
#     for id_ in tqdm(ids, desc='download_ncbi: Downloading entries from NCBI.'):
        
#         output_path = os.path.join(output_dir, f'{id_}.gbk')
#         if os.path.exists(output_path):
#             continue 
#         try:
#             with entrez.efetch(id=id_, db=db, rettype='gb', retmode='text') as result:
#                 content = result.read()
#         except:
#             print(f'download_ncbi: Failed to obtain data for {id_}.')
#             continue 

#         with open(output_path, 'w') as f:
#             f.write(content)
