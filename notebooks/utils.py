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
from files.pdb import PDBFile, cif_to_pdb, sph_to_pdb, ATOMS
import orjson

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

import matplotlib.pyplot as plt
from cycler import cycler
from datetime import date

today = date.today().strftime("%m%d%Y")


def get_alphafold_metadata(paths:str, output_path:str=None) -> pd.DataFrame:
    '''
    '''
    if (output_path is not None) and os.path.exists(output_path):
        return pd.read_csv(output_path)
    
    df = list() 
    for path in paths:

        try:
            output = AlphaFoldOutput(path)
        except Exception as err:
            print(f'get_alphafold_metadata: Problem loading output at {path}, {err}')
            continue 

        assert output.get_num_proteins() == 1, f'get_alphafold_metadata: Expected 1 unique protein per structure, but got {output.get_num_proteins()} in {path}'

        row = dict()
        row['path'] = os.path.abspath(path)
        row['name'] = output.name 
        row['msa_unpaired_num_seqs'] = output.get_msas()[0]['unpaired'].count('>') - 1
        row['msa_paired_num_seqs'] = output.get_msas()[0]['unpaired'].count('>') - 1
        row['iptms'] = list(output.get_iptms(mean_pool=False).values())
        row['ptms'] = list(output.get_ptms(mean_pool=False).values())
        row['best_model'] = output.best_model
        row['iptm_best_model'] = list(output.get_iptms(mean_pool=False, models=[output.best_model]).values())[0]
        row['ptm_best_model'] = list(output.get_ptms(mean_pool=False, models=[output.best_model]).values())[0]
        row['name'] = output.name 
        row['num_seeds'] = len(output.data['modelSeeds'])
        row['num_protein_chains'] = output.get_num_protein_chains()
        row['num_proteins'] = output.get_num_protein_chains()
        df.append(row)

    df = pd.DataFrame(df)

    if (output_path is not None):
        df.to_csv(output_path, index=False)

    return df

# Functions for JSON parsing. 
# -------------------------------------------------------------------------------------------------------------------------------------------------

def default(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    raise TypeError


def load_json(path):
    '''Read a JSON file from the input path using the faster orjson parser. '''
    try:
        with open(path, 'rb') as f:
            data = f.read()
            data = orjson.loads(data)
        return data
    except Exception as err:
        print(f'load_json: Could not decode {path}, {err.msg}')
        print(data)
        return None
# -------------------------------------------------------------------------------------------------------------------------------------------------


INTERPROSCAN_FIELDS = ['gene_id', 'checksum', 'length', 'analysis', 'accession', 'description', 'start','stop', 'e_value', 'status', 'date', 'interpro_accession', 'interpro_description', 'go_terms', 'pathways']

PROJECT_IDS = pd.read_csv('project_ids.csv', index_col=0).project_id.to_dict()
FORWARD_READS_PATHS = pd.read_csv('reads_paths.csv', index_col=0).forward_reads_path.to_dict()
REVERSE_READS_PATHS = pd.read_csv('reads_paths.csv', index_col=0).reverse_reads_path.to_dict()

colors = ['#5E5247','#6C6154','#7B7061','#8A8170','#84886C','#737F69','#65776B','#5F7472','#63727D','#676F88','#6B6B93']
colors = ['#6A5A4B', '#7A6956', '#8B7A62', '#9A8B68','#8B9164','#768C66','#64836C','#5D7F76','#607B83','#66778C','#6B6F91']
colors = ['#8A5E3B', '#A0704F', '#B38A5A', '#9AA45A', '#7FA65E', '#62A46A', '#4F9D78', '#4F9A8C', '#5A93A0', '#6A89B0', '#7A7FB8']
GENOME_ID_COLORS = dict(zip(['bz_0', 'bz_1', 'bz_2', 'bz_3', 'bz_4', 'bz_5', 'bz_7', 'bz_8', 'bz_9', 'bz_10', 'bz_11'], colors))

plt.rcParams['axes.prop_cycle'] = cycler(color=colors)


START_CODONS = ['ATG', 'GTG', 'TTG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']

reverse_complement = lambda seq : str(Seq(seq).reverse_complement())
# get_gene_id = lambda path : os.path.basename(path).replace('_pred.txt', '').replace('genes_', '') # Remove the prefix and file extension. 
get_gene_id = lambda string : re.search(r'orfm.bz_\d+\.\d+_\d+', string).group(0) if (re.search(r'orfm.bz_\d+\.\d+_\d+', string) is not None) else None
get_genome_id = lambda string : re.search(r'bz_\d+', string).group(0)

def load_colabfold_plddts(path:str='../data/genes/colabfold/scores.json', gene_ids=None):
    with open(path, 'r') as f:
        plddts = json.load(f)
    if gene_ids is not None:
        plddts = {gene_id:plddts_ for gene_id, plddts_ in plddts.items() if (gene_id in gene_ids)}     
    return plddts


# Functions for analyzing interfaces in multi-chain AlphaFold structures. 
# -------------------------------------------------------------------------------------------------------------------------------------------------

def get_interfaces(contact_probs_df:np.ndarray, min_contact_prob=0.5):
    '''Use the contact_probs in the AlphaFold output to locate potential inter-chain residue contacts.
    
    :param contact_probs_df: A square DataFrame containing the contact_probs for a given structure. 
    :param min_contact_prob: The minimum contact_prob to say whether or not two residues are likely to be in contact. 
    :returns: A two-tuple containing (1) the number of predicted contacts and (2) a square DataFrame containing boolean
        values indicating the interface contacts. 
    '''
    token_chain_ids = contact_probs_df.index.to_numpy()
    mask = (contact_probs_df.values > min_contact_prob) # Require a minimum contact probability. 
    mask = mask & (np.expand_dims(token_chain_ids, axis=1) != token_chain_ids)  # Don't include intra-chain contacts. 
    # print(f'get_interfaces: {mask.sum().sum()} inter-chain residues predicted to be in contact.')
    return mask.ravel().sum(), pd.DataFrame(mask, index=token_chain_ids, columns=token_chain_ids)


def get_interface(contact_probs_df, chain_ids=['A', 'B'], min_contact_prob:float=0.5):
    '''Get the interface specifically between the two specified chains.
    
    :param contact_probs_df: A square DataFrame containing the contact_probs for a given structure. 
    :param chain_ids: The IDs for the chains to find the interface between. No more than two chains are expected. 
    :param min_contact_prob: The minimum contact_prob to say whether or not two residues are likely to be in contact. 
    :returns: A square DataFrame containing boolean values indicating the interface contacts. 
    ''' 
    assert len(chain_ids) == 2, f'get_interface: Expected two chains, but got {len(chain_ids)}.'
    df = get_interfaces(contact_probs_df, min_contact_prob=min_contact_prob)[1] # Get all interface contacts. 
    mask = (df.index.values == chain_ids[0]).reshape(-1, 1) & (df.columns.values == chain_ids[1])
    mask = mask & (df.values) # Also make sure the residues are in contact. 
    n = mask.sum(axis=None)
    # print(f'get_interface: Found {n} residues at the interface of chains {chain_ids[0]} and {chain_ids[1]}.')
    return mask 


def has_interface(contact_probs_df, chain_ids=None, min_contact_prob:float=0.5):
    '''

    :param contact_probs_df: A square DataFrame containing the contact_probs for a given structure. 
    :param chain_ids: The IDs for the chains to find the interface between. No more than two chains are expected. 
    :param min_contact_prob: The minimum contact_prob to say whether or not two residues are likely to be in contact. 
    :returns: A two-tuple with the first element being the number of predicted contacts at the specified threshold, and the second element
        being a boolean indicating whether or not the two chains are in contact. 
    '''
    mask = get_interface(contact_probs_df, chain_ids=chain_ids, min_contact_prob=min_contact_prob)
    n = mask.sum(axis=None)
    return n, n > 0


def get_interface_idxs(contact_probs_df, chain_ids=['A', 'B'], min_contact_prob:float=0.5, shift:bool=True):
    '''Get the indices of the residues in the two specified chains participating in an interface.
    
    :param contact_probs_df: A square DataFrame containing the contact_probs for a given structure. 
    :param chain_ids: The IDs for the chains to find the interface between. No more than two chains are expected. 
    :param min_contact_prob: The minimum contact_prob to say whether or not two residues are likely to be in contact. 
    '''
    token_chain_ids = contact_probs_df.index.values # Get the token chain IDs from the interface DataFrame. 

    # def _shift_idxs(idxs:np.ndarray, chain_id:str='A'):
    #     delta = np.where(token_chain_ids == chain_id)[0][0] # Get the first occurrence of the chain ID. 
    #     return idxs - delta # Subtract the shift from the indices. 

    def _get_shift(chain_id:str):
        return np.where(token_chain_ids == chain_id)[0][0] # Get the first occurrence of the chain ID. 

    mask = get_interface(contact_probs_df, chain_ids=chain_ids, min_contact_prob=min_contact_prob)
    idxs =  dict(zip(chain_ids, np.where(mask)))
    shifts = {chain_id:_get_shift(chain_id) for chain_id in chain_ids}
    
    if shift:
        return idxs, {chain_id:idxs_ - shifts[chain_id] for chain_id, idxs_ in idxs.items()}
    else:
        return idxs
# -------------------------------------------------------------------------------------------------------------------------------------------------


COLABFOLD_FILE_NAME_PATTERN = r'({gene_ids})_relaxed_rank_001_alphafold2_ptm_model_\d_seed_000.pdb'

def copy_structures(gene_ids:list, source_dir:str=None, file_name_pattern=COLABFOLD_FILE_NAME_PATTERN, output_dir:str=None):
    '''Copy structures from the specified source directory to the output directory, and create the output directory if it
    does not already exist. The structure files in the output directory will be renamed according to their gene ID.
    
    :param gene_ids: The list of gene IDs to copy the structures for. 
    :param source_dir: The directory containing the structures. 
    :param file_name_pattern: The regex pattern matching the structure files to copy over. There should be one capturing group 
        matching the gene ID in the file name. 
    :param output_dir: The destination directory for the copied structure files.
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
    
    :param seqs: pd.Series containing the sequences to compare.
    :param mode: The alignment mode, either local or global. The mode is global by default. 
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

# def load_bin_analysis(path:str, gene_name='rps3', drop=['SR-VP_05_06_2024_coassembly']):
#     '''Convert bin analysis output from wide-form to long-form data. Note that an rps3 sequence observed in multiple samples is identified according to a 
#     99 percent identity grouping.
#     column 0 : Protein identifier in the form {gene_id} | {bin_id}. The gene ID specified here is repeated under the corresponding sample column, so this can be excluded.
#     column 1 : Taxonomy, which is not very informative, as well as the confidence of the assignment.
#     column 2+ : Column names are the ggKbase sample ID, and entries are of the form "{gene_id} in {bin_id} ({bin_coverage}x) with contig info - feature cnt: {num_genes} | size: {contig_length} | cov: {contig_coverage}x | gc: {contig_gc_percent}%
#     '''
#     df = pd.read_csv(path, sep='\t').drop(columns=drop)
#     print(f'load_bin_analysis: Loaded {len(df)} entries from {path}.')
#     df = df.iloc[:, 1:] # Ignore the first column.
#     df = df.rename(columns={'phylogeny winner':'taxonomy'})
#     sample_ids = [col for col in df.columns if (col != 'taxonomy')]
    
#     contig_info_pattern = r'feature cnt: (?P<num_genes>[\d]+) \| size: (?P<contig_length>[\d]+) \| cov: (?P<coverage>[^\s]+)x \| gc: (?P<gc_percent>[^\s]+)%'
#     pattern = fr'(?P<ggkbase_gene_id>[^\s]+) in (?P<ggkbase_bin_id>[^\s]+) \((?P<bin_coverage>[^x]+)x\) with contig info - {contig_info_pattern}'
    
#     df_ = list()    
#     for i, row in df.iterrows():
#         for sample_id in sample_ids:
#             row_ = {'gene_id':f'{gene_name}_{i}', 'sample_id':sample_id, 'taxonomy':row['taxonomy']}
#             if type(row[sample_id]) == float: # Check if the value is empty. 
#                 row_.update({'coverage':0, 'contig_info':None, 'ggkbase_gene_id':None}) 
#             else:
#                 row_.update(re.search(pattern, row[sample_id]).groupdict())
#             df_.append(row_)
#     df_ = pd.DataFrame(df_)
#     df_['coverage'] = df_['coverage'].astype(float)
#     return df_



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
