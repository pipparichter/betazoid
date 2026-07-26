import sys 
sys.path.append('/home/prichter/Documents/banfield/betazoid/src/files/')
sys.path.append('/home/prichter/Documents/banfield/betazoid/src/')

from fasta import FASTAFile
from tmhmm import TMHMMFile
from gfa import GFAFile 
from bam import BamFile
from dssp import DSSPFile
from alphafold import AlphaFoldInputFile, AlphaFoldOutput, AlphaFoldServerOutput
from blast import BLASTFile
from files.pdb import PDBFile

import os 
import re 
import pandas as pd 
import io 
from Bio import SeqIO
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
from Bio.Align import PairwiseAligner 
import itertools 
import subprocess
import json

INTERPROSCAN_FIELDS = ['gene_id', 'checksum', 'length', 'analysis', 'accession', 'description', 'start','stop', 'e_value', 'status', 'date', 'interpro_accession', 'interpro_description', 'go_terms', 'pathways']

PROJECT_IDS = pd.read_csv('project_ids.csv', index_col=0).project_id.to_dict()
FORWARD_READS_PATHS = pd.read_csv('reads_paths.csv', index_col=0).forward_reads_path.to_dict()
REVERSE_READS_PATHS = pd.read_csv('reads_paths.csv', index_col=0).reverse_reads_path.to_dict()

colors = ['#5E5247','#6C6154','#7B7061','#8A8170','#84886C','#737F69','#65776B','#5F7472','#63727D','#676F88','#6B6B93']
colors = ['#6A5A4B', '#7A6956', '#8B7A62', '#9A8B68','#8B9164','#768C66','#64836C','#5D7F76','#607B83','#66778C','#6B6F91']
colors = ['#8A5E3B', '#A0704F', '#B38A5A', '#9AA45A', '#7FA65E', '#62A46A', '#4F9D78', '#4F9A8C', '#5A93A0', '#6A89B0', '#7A7FB8']
GENOME_ID_COLORS = dict(zip(['bz_0', 'bz_1', 'bz_2', 'bz_3', 'bz_4', 'bz_5', 'bz_7', 'bz_8', 'bz_9', 'bz_10', 'bz_11'], colors))

START_CODONS = ['ATG', 'GTG', 'TTG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']

reverse_complement = lambda seq : str(Seq(seq).reverse_complement())
# get_gene_id = lambda path : os.path.basename(path).replace('_pred.txt', '').replace('genes_', '') # Remove the prefix and file extension. 
get_gene_id = lambda string : re.search(r'orfm.bz_\d+\.\d+_\d+', string).group(0)
get_genome_id = lambda string : re.search(r'bz_\d+', string).group(0)

def load_colabfold_plddts(path:str='../data/genes/colabfold/scores.json', gene_ids=None):
    with open(path, 'r') as f:
        plddts = json.load(f)
    if gene_ids is not None:
        plddts = {gene_id:plddts_ for gene_id, plddts_ in plddts.items() if (gene_id in gene_ids)}     
    return plddts

COLABFOLD_FILE_NAME_PATTERN = r'({gene_ids})_relaxed_rank_001_alphafold2_ptm_model_\d_seed_000.pdb'

def copy_structures(gene_ids:list, source_dir:str=None, file_name_pattern=COLABFOLD_FILE_NAME_PATTERN, output_dir:str=None):
    '''Copy structures from the specified source directory to the output directory, and create the output directory if it
    does not already exist. The structure files in the output directory will be renamed according to their gene ID.
    
    : param gene_ids: The list of gene IDs to copy the structures for. 
    : param source_dir: The directory containing the structures. 
    : param file_name_pattern: The regex pattern matching the structure files to copy over. There should be one capturing group 
        matching the gene ID in the file name. 
    : param output_dir: The destination directory for the copied structure files.
    '''
    os.makedirs(output_dir, exist_ok=True)
    file_name_pattern = file_name_pattern.format(gene_ids='|'.join(gene_ids))

    for path in glob.glob(os.path.join(source_dir, '**', '*')):
        match_ = re.search(file_name_pattern, path)
        if match_ is None:
            continue 
        else:
            ext = os.path.splitext(path)[-1] # Get the source file extension. 
            output_path = os.path.join(output_dir, match_.group(1) + f'{ext}')
            subprocess.run(f'cp {path} {output_path}', shell=True, check=True)



def apply_filters(filters:dict, df:pd.DataFrame, return_idxs=False):
    keep = None
    for name, filter_ in filters.items(): # mask should correspond to things failing the test. 
        keep = ~filter_ if (keep is None) else (keep & ~filter_)
        print(f'apply_filters: {filter_.sum()} entries removed by {name}.')

    if return_idxs:
        return df[keep].copy(), df[~keep].index.values  
    else:
        return df[keep].copy()


def get_pairwise_alignments(seqs, mode:str='global', metric='identities', normalize:bool=True):
    '''Returns an n-by-n DataFrame containing the pairwise identity of each sequence in the input. 
    
    : param seqs: pd.Series containing the sequences to compare.
    : param mode: The alignment mode, either local or global. The mode is global by default. 
    '''
    n = len(seqs)
    df = pd.DataFrame(np.zeros((n, n)), columns=seqs.index, index=seqs.index)
    aligner =  PairwiseAligner()
    aligner.mode = mode 

    for ((id_1, seq_1), (id_2, seq_2)) in itertools.combinations(seqs.to_dict().items(), 2):
        alignment = aligner.align(seq_1, seq_2)[0]

        n = len(alignment[0]) if normalize else 1
        df.loc[id_1, id_2] = getattr(alignment.counts(), metric) / n
        df.loc[id_2, id_1] = getattr(alignment.counts(), metric) / n
    
    # np.fill_diagonal(df.values, 1) # Because the self-alignments are not computed. 
    # for id_ in seqs.index:
    #     df.loc[id_, id_] = np.nan
    return df



REDUCED_ALPHABET = {'A':'A', 'V':'A', 'L':'A', 'I':'A', 'M':'A', 'F':'R', 'W':'R','Y':'R','K':'+','R':'+','H':'+','D':'-','E':'-','S':'P','T':'P','N':'P','Q':'P','G':'G','P':'P','C':'C', '.':'.'}

# ! muscle -align ../data/genes/cluster_1.fa -output ../data/genes/cluster_1.afa
class MSA():
    gap_symbol = '.'
    def __init__(self, arr, ids):

        self.ids = ids 
        self.arr = arr 
        self.seqs = np.array([''.join(row) for row in arr])
        self.n_cols = arr.shape[-1]

    @classmethod
    def from_file(cls, path:str='../data/genes/cluster_1.afa'):

        df = FASTAFile.from_file(path).to_df()
        # print(df)
        ids, arr = df.index.values, np.array([list(seq.replace('-', MSA.gap_symbol)) for seq in df.seq])
        return MSA.from_array(arr, ids) 

    @classmethod
    def from_string(cls, content:str, af3:bool=True):
        f = io.StringIO(content)
        records = [record for record in SeqIO.parse(f, 'fasta')]
        ids = [record.id for record in records]

        seqs = [str(record.seq).replace('-', MSA.gap_symbol) for record in records]
        if af3: # AF3 format is different, contains lowercase characters to indicate insertions relative to the query.
            seqs = [re.sub(r'[a-z]', '', seq) for seq in seqs]
        n_cols, n_rows = len(seqs[0]), len(seqs)
        arr = np.array([list(seq) for seq in seqs])
        return cls(arr, ids)

    def __len__(self):
        return len(self.ids)
    
    def __getitem__(self, id_):
        assert id_ in self.ids, f'MSA.__getitem__: ID {id_} is missing in the MSA.'
        return self.seqs[self.ids == id_][0]
    
    def to_array(self, alphabet:dict=None):
        '''Convert the MSA loaded from the FASTA file into a two-dimensional numpy array, where each entry is a single residue.'''
        return np.vectorize(alphabet.get)(self.arr.copy()) if (alphabet is not None) else self.arr.copy()
    
    def to_df(self, alphabet:dict=None):
        df = pd.DataFrame(index=self.ids)
        df['seq'] = self.seqs
        df['seq'] = df.seq.str.replace(alphabet) if (alphabet is not None) else df['seq']
        return df        
    
    def map_idx_from_msa(self, idx, gene_id:str):
        '''Convert the index of a residue in the MSA to the index of a residue in one of the aligned sequences.'''
        seq = self[gene_id]
        assert seq[idx] != MSA.gap_symbol, f'get_idx: The input index corresponds to a gap in the aligned {gene_id}.'
        n_gaps = seq[:idx].count(MSA.gap_symbol) # Get the number of gaps which occur before the requested index.
        return idx - n_gaps
    
    def map_idx_to_msa(self, idx, gene_id:str):
        '''Convert the index of a residue in one of the sequences to an index in the MSA.'''
        n, n_gaps = 0, 0
        for aa in self[gene_id]:
            if (n == idx) and (aa != MSA.gap_symbol):
                break # So that it doesn't exit if idx = 0 and the first symbol is a gap. 
            n += int(aa != MSA.gap_symbol)
            n_gaps += int(aa == MSA.gap_symbol)
        return n + n_gaps
    
    @classmethod
    def from_array(cls, arr, ids=None):
        ids = np.arange(len(arr)) if (ids is None) else ids
        arr = np.array(arr)
        ids = np.array(ids)
        return MSA(arr, ids)

    def get_mean_gap_fraction(self):
        return np.mean([np.mean(col == MSA.gap_symbol) for col in self.arr.T])


    def get_consensus(self):
        consensus = list()
        for col in self.arr.T:
            symbols, counts = np.unique(col, return_counts=True)
            consensus.append(symbols[np.argsort(counts)][-1])
        return np.array(consensus)

    def show(self, start:int=None, stop:int=None):
        start = 0 if (start is None) else start
        stop = self.n_cols if (stop is None) else stop

        for id_, seq in zip(self.ids, self.seqs):
            print(f'{start}\t{seq[start:stop]}\t{stop}\t{id_}')

    def get_map_idxs(self):
        map_idxs = dict()
        for id_, row in zip(self.ids, self.arr.copy()):
            row_map_idxs  = np.zeros(len(row), dtype=int)
            row_map_idxs[np.where(row != MSA.gap_symbol)[0]] = np.arange((row != MSA.gap_symbol).sum()) # Fill in the values with the ungapped index positions. 
            row_map_idxs[np.where(row == MSA.gap_symbol)[0]] = -1
            map_idxs[id_] = row_map_idxs
        return map_idxs



def plot_foldseek_search_coverage(msa:MSA, foldseek_search_df:pd.DataFrame):

    fig, ax = plt.subplots(figsize=(5, 4))

    figure_df = foldseek_search_df.copy()
    # Get the indices relative to the MSA. 
    figure_df['msa_start'] = figure_df.apply(lambda row : msa.map_idx_to_msa(row.qstart, getattr(row, 'query')), axis=1) 
    figure_df['msa_stop'] = figure_df.apply(lambda row : msa.map_idx_to_msa(row.qend, getattr(row, 'query')), axis=1) 

    positions = np.arange(msa.n_cols)
    heights = np.zeros(len(positions))
    for row in figure_df.itertuples():
        heights[np.arange(row.msa_start, row.msa_stop)] += 1 

    ax.bar(positions, heights, color='lightgray', lw=0, edgecolor='black')
    ax.set_ylabel('num. Foldseek hits covering residue')
    ax.set_xlabel('position')
    plt.show()


# Want to map original indices to an alignment. Return an array which is the same shape as the alignment, but with each position
# indicating the corresponding sequence position. 
# We want to be able to do something like X[i, idxs] = new_chars to map other symbols (e.g. structural symbols) onto the alignment.




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
