import sys 
sys.path.append('/home/prichter/Documents/banfield/betazoid/src/files/')
sys.path.append('/home/prichter/Documents/banfield/betazoid/src/')
sys.path.append('/home/prichter/Documents/banfield/betazoid/scripts/')

from fasta import FASTAFile
from tmhmm import TMHMMFile
from gfa import GFAFile
from msa import MSAFile
from bam import BamFile
from dssp import DSSPFile
from alphafold import AlphaFoldInputFile, AlphaFoldOutput, AlphaFoldServerOutput
from colabfold import ColabFoldOutput
from blast import BLASTFile
from files.pdb import PDBFile, cif_to_pdb, sph_to_pdb, ATOMS
import orjson
import ast

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
from matplotlib.colors import LinearSegmentedColormap, to_hex
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner 
import itertools 
import subprocess
import json

import matplotlib.pyplot as plt
from cycler import cycler
from datetime import date

FOLDSEEK_FIELD_MAP = dict()
FOLDSEEK_FIELD_MAP['query'] = 'query_id'
FOLDSEEK_FIELD_MAP['target'] = 'target_id'
FOLDSEEK_FIELD_MAP['evalue'] = 'e_value'
FOLDSEEK_FIELD_MAP['gapopen'] = 'num_gaps'
FOLDSEEK_FIELD_MAP['pident'] = 'percent_identity'
FOLDSEEK_FIELD_MAP['fident'] = 'fraction_identical'
FOLDSEEK_FIELD_MAP['nident'] = 'num_identical'
FOLDSEEK_FIELD_MAP['qstart'] = 'query_start'
FOLDSEEK_FIELD_MAP['qend'] = 'query_end'
FOLDSEEK_FIELD_MAP['qlen'] = 'query_length'
FOLDSEEK_FIELD_MAP['tstart'] = 'target_start'
FOLDSEEK_FIELD_MAP['tend'] = 'target_end'
FOLDSEEK_FIELD_MAP['tlen'] = 'target_length'
FOLDSEEK_FIELD_MAP['alnlen'] = 'alignment_length'
FOLDSEEK_FIELD_MAP['bits'] = 'bit_score'
FOLDSEEK_FIELD_MAP['cigar'] = 'cigar'
FOLDSEEK_FIELD_MAP['qseq'] = 'query_seq'
FOLDSEEK_FIELD_MAP['tseq'] = 'target_seq'
FOLDSEEK_FIELD_MAP['qheader'] = 'query_header'
FOLDSEEK_FIELD_MAP['theader'] = 'target_header'
FOLDSEEK_FIELD_MAP['qaln'] = 'query_alignment'
FOLDSEEK_FIELD_MAP['taln'] = 'target_alignment'
FOLDSEEK_FIELD_MAP['mismatch'] = 'num_mismatches'
FOLDSEEK_FIELD_MAP['qcov'] = 'query_coverage'
FOLDSEEK_FIELD_MAP['tcov'] = 'target_coverage'
FOLDSEEK_FIELD_MAP['taxid'] = 'taxonomy_id'
FOLDSEEK_FIELD_MAP['taxname'] = 'taxonomy'
FOLDSEEK_FIELD_MAP['taxlineage'] = 'lineage'
FOLDSEEK_FIELD_MAP['lddt'] = 'lddt'
FOLDSEEK_FIELD_MAP['lddtfull'] = 'lddt_full'
FOLDSEEK_FIELD_MAP['qtmscore'] = 'query_tm_score'
FOLDSEEK_FIELD_MAP['ttmscore'] = 'target_tm_score'
FOLDSEEK_FIELD_MAP['alntmscore'] = 'alignment_tm_score'
FOLDSEEK_FIELD_MAP['rmsd'] = 'rmsd'
FOLDSEEK_FIELD_MAP['prob'] = 'probability'

FOLDSEEK_FIELDS = 'query target evalue gapopen pident fident nident qstart qend qlen tstart tend tlen alnlen bits cigar qseq tseq qheader theader qaln taln mismatch qcov tcov taxid taxname taxlineage lddt lddtfull qtmscore ttmscore alntmscore rmsd prob'
FOLDSEEK_FIELDS = [FOLDSEEK_FIELD_MAP.get(field) for field in FOLDSEEK_FIELDS.split(' ')]

today = date.today().strftime("%m%d%Y")

INTERPROSCAN_FIELDS = ['gene_id', 'checksum', 'length', 'analysis', 'accession', 'description', 'start','stop', 'e_value', 'status', 'date', 'interpro_accession', 'interpro_description', 'go_terms', 'pathways']

PROJECT_IDS = pd.read_csv('project_ids.csv', index_col=0).project_id.to_dict()
FORWARD_READS_PATHS = pd.read_csv('reads_paths.csv', index_col=0).forward_reads_path.to_dict()
REVERSE_READS_PATHS = pd.read_csv('reads_paths.csv', index_col=0).reverse_reads_path.to_dict()


colors = ['#8A5E3B', '#A0704F', '#B38A5A', '#9AA45A', '#7FA65E', '#62A46A', '#4F9D78', '#4F9A8C', '#5A93A0', '#6A89B0', '#7A7FB8']
GENOME_ID_COLORS = dict(zip(['bz_0', 'bz_1', 'bz_2', 'bz_3', 'bz_4', 'bz_5', 'bz_7', 'bz_8', 'bz_9', 'bz_10', 'bz_11'], colors))

plt.rcParams['axes.prop_cycle'] = cycler(color=colors)

START_CODONS = ['ATG', 'GTG', 'TTG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']

reverse_complement = lambda seq : str(Seq(seq).reverse_complement())
# get_gene_id = lambda path : os.path.basename(path).replace('_pred.txt', '').replace('genes_', '') # Remove the prefix and file extension. 
get_gene_id = lambda string : re.search(r'orfm.bz_\d+\.\d+_\d+', string).group(0) if (re.search(r'orfm.bz_\d+\.\d+_\d+', string) is not None) else None
get_genome_id = lambda string : re.search(r'bz_\d+', string).group(0)

REDUCED_ALPHABET = {'A':'A', 'V':'A', 'L':'A', 'I':'A', 'M':'A', 'F':'R', 'W':'R','Y':'R','K':'+','R':'+','H':'+','D':'-','E':'-','S':'P','T':'P','N':'P','Q':'P','G':'G','P':'P','C':'C', '.':'.'}


HYDROPHOBICITY_SCALE = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5}
CHARGE_SCALE = {'D': -1, 'E': -1, 'K': 1, 'R': 1, 'H': 0, 'A': 0, 'C': 0, 'F': 0, 'G': 0, 'I': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}


def plot_hydrophobicity(msa:MSAFile, window_size:int=5, step_size:int=1, x_min:int=0, x_max:int=100, ax:plt.Axes=None, legend:bool=False, palette=dict()):
    '''
    
    '''
    is_gap = lambda seq : np.all(np.array(list(seq)) == msa.gap_symbol)
    get_hydrophobicity = lambda seq : None if is_gap(seq) else np.mean([HYDROPHOBICITY_SCALE[aa] for aa in seq if (aa != msa.gap_symbol)])
    # get_charge = lambda seq : None if is_gap(seq) else np.mean([CHARGE_SCALE[aa] for aa in seq if (aa != msa.gap_symbol)])

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))

    x = list(range(0, msa.n_cols, step_size))
    y = list()

    for id_, seq in zip(msa.ids, msa.seqs):
        seq = seq.replace('X', msa.gap_symbol)
        windows = [seq[i:i + window_size] for i in x]
        y.append(np.array([get_hydrophobicity(window) for window in windows]))
        ax.plot(x, y[-1], label=id_, color=palette.get(id_, 'lightgray'))

    get_mean = lambda arr : None if np.all(np.isnan(arr)) else np.mean(arr[~np.isnan(arr)])

    y = np.array(y)
    y_mean = [get_mean(y_.astype(float)) for y_ in y.T]

    ax.scatter(x, y_mean, color='black', zorder=100)

    ax.set_ylabel('hydrophobicity')
    ax.set_xlabel('position')
    ax.set_xlim(xmin=x_min, xmax=x_max)

    if legend:
        ax.legend()



def get_msas_unpaired(gene_ids:list, dir_path:str='../data/genes/mmseqs/'):
    '''Load the unpaired MSAs for the specified gene IDs.
    
    :param gene_ids: The IDs of the genes to construct paired MSAs for. These should all be from the same genome.   
    :param dir_path: The path to the directory where the MSAs are stored. This function assumes the a3m files have names matching the gene ID of the query 
        sequence. 
    '''
    paths = {gene_id:os.path.join(dir_path, f'{gene_id}.a3m') for gene_id in gene_ids}
    return {gene_id:str(FASTAFile.from_file(path)) for gene_id, path in paths.items()}


 
def map_dssp_to_msa(msa, dssp_output_dir:str='../data/genes/dssp/'):
    '''Read in the DSSP results and map the structure codes to the MSAFile object.
    
    :param msa: An MSAFile object loaded from an afa file.
    :param dddp_output_dir: The directory where the DSSP output is written.
    :returns: A new MSAFile object with the amino acids replaced by DSSP structure codes. 
    '''

    ids, arr = list(), list()
    for gene_id, map_idxs in msa.get_map_idxs().items():
        path = os.path.join(dssp_output_dir, f'{gene_id}.dssp') # Expects the files to be named according to the gene ID. 
        row = np.array(list(str(DSSPFile.from_file(path))) + ['.'])
        arr.append(row[map_idxs])
        ids.append(gene_id)
        
    return MSAFile.from_array(arr, ids)


def get_fold_metadata(paths:str, output_path:str=None, parser=AlphaFoldOutput) -> pd.DataFrame:
    '''Collect metadata about the AlphaFold or ColabFoldjobs stored at the given paths in a pandas DataFrame. 

    :param paths: A list of paths specifying the output location. For AlphaFoldOutputs, this is a directory name. For
        ColabFold outputs, which are not organized into directories, this is a file prefix.
    :param output_path:
    :param parser:
    '''
    if (output_path is not None) and os.path.exists(output_path):
        df = pd.read_csv(output_path)
        df['iptms'] = df.iptms.apply(ast.literal_eval) # The list of ipTMS gets stored as a string, so need to convert to a list. 
        df['ptms'] = df.ptms.apply(ast.literal_eval) # The list of pTMS gets stored as a string, so need to convert to a list. 
        return df 
    
    df = list() 
    for path in tqdm(paths, desc='get_fold_metadata'):

        output = parser(path)
        assert len(output.get_chain_to_chain_id_map(chain_type='protein')) == 1, f'get_fold_metadata: Expected 1 unique protein per structure, but got {output.get_num_proteins()} in {path}'

        row = dict()
        row['path'] = os.path.abspath(path)
        row['name'] = output.name 
        row['msa_num_seqs'] = output.get_msas()[0].get('unpaired', '').count('>')
        row['iptms'] = list(output.get_iptms(mean_pool=False).values())
        row['ptms'] = list(output.get_ptms(mean_pool=False).values())
        row['best_model'] = output.best_model
        row['iptm_best_model'] = output.get_iptms(best_model=True)
        row['ptm_best_model'] = output.get_ptms(best_model=True)
        row['num_seeds'] = output.get_num_seeds()
        row['num_protein_chains'] = len(output.get_chain_ids(chain_type='protein'))
        row['num_proteins'] = len(output.get_chain_to_chain_id_map(chain_type='protein'))

        if isinstance(output, ColabFoldOutput):
            row.update(output.get_msa_metadata()[0])
            row['plddts_best_model'] = output.get_plddts(best_model=True)

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




# # Functions for analyzing interfaces in multi-chain AlphaFold structures. 
# # -------------------------------------------------------------------------------------------------------------------------------------------------

# def get_interfaces(contact_probs_df:np.ndarray, min_contact_prob=0.5):
#     '''Use the contact_probs in the AlphaFold output to locate potential inter-chain residue contacts.
    
#     :param contact_probs_df: A square DataFrame containing the contact_probs for a given structure. 
#     :param min_contact_prob: The minimum contact_prob to say whether or not two residues are likely to be in contact. 
#     :returns: A two-tuple containing (1) the number of predicted contacts and (2) a square DataFrame containing boolean
#         values indicating the interface contacts. 
#     '''
#     token_chain_ids = contact_probs_df.index.to_numpy()
#     mask = (contact_probs_df.values > min_contact_prob) # Require a minimum contact probability. 
#     mask = mask & (np.expand_dims(token_chain_ids, axis=1) != token_chain_ids)  # Don't include intra-chain contacts. 
#     # print(f'get_interfaces: {mask.sum().sum()} inter-chain residues predicted to be in contact.')
#     return mask.ravel().sum(), pd.DataFrame(mask, index=token_chain_ids, columns=token_chain_ids)


# def get_interface(contact_probs_df, chain_ids=['A', 'B'], min_contact_prob:float=0.5):
#     '''Get the interface specifically between the two specified chains.
    
#     :param contact_probs_df: A square DataFrame containing the contact_probs for a given structure. 
#     :param chain_ids: The IDs for the chains to find the interface between. No more than two chains are expected. 
#     :param min_contact_prob: The minimum contact_prob to say whether or not two residues are likely to be in contact. 
#     :returns: A square DataFrame containing boolean values indicating the interface contacts. 
#     ''' 
#     assert len(chain_ids) == 2, f'get_interface: Expected two chains, but got {len(chain_ids)}.'
#     df = get_interfaces(contact_probs_df, min_contact_prob=min_contact_prob)[1] # Get all interface contacts. 
#     mask = (df.index.values == chain_ids[0]).reshape(-1, 1) & (df.columns.values == chain_ids[1])
#     mask = mask & (df.values) # Also make sure the residues are in contact. 
#     n = mask.sum(axis=None)
#     # print(f'get_interface: Found {n} residues at the interface of chains {chain_ids[0]} and {chain_ids[1]}.')
#     return mask 


# def has_interface(contact_probs_df, chain_ids=None, min_contact_prob:float=0.5):
#     '''

#     :param contact_probs_df: A square DataFrame containing the contact_probs for a given structure. 
#     :param chain_ids: The IDs for the chains to find the interface between. No more than two chains are expected. 
#     :param min_contact_prob: The minimum contact_prob to say whether or not two residues are likely to be in contact. 
#     :returns: A two-tuple with the first element being the number of predicted contacts at the specified threshold, and the second element
#         being a boolean indicating whether or not the two chains are in contact. 
#     '''
#     mask = get_interface(contact_probs_df, chain_ids=chain_ids, min_contact_prob=min_contact_prob)
#     n = mask.sum(axis=None)
#     return n, n > 0


# def get_interface_idxs(contact_probs_df, chain_ids=['A', 'B'], min_contact_prob:float=0.5, shift:bool=True):
#     '''Get the indices of the residues in the two specified chains participating in an interface.
    
#     :param contact_probs_df: A square DataFrame containing the contact_probs for a given structure. 
#     :param chain_ids: The IDs for the chains to find the interface between. No more than two chains are expected. 
#     :param min_contact_prob: The minimum contact_prob to say whether or not two residues are likely to be in contact. 
#     '''
#     token_chain_ids = contact_probs_df.index.values # Get the token chain IDs from the interface DataFrame. 

#     # def _shift_idxs(idxs:np.ndarray, chain_id:str='A'):
#     #     delta = np.where(token_chain_ids == chain_id)[0][0] # Get the first occurrence of the chain ID. 
#     #     return idxs - delta # Subtract the shift from the indices. 

#     def _get_shift(chain_id:str):
#         return np.where(token_chain_ids == chain_id)[0][0] # Get the first occurrence of the chain ID. 

#     mask = get_interface(contact_probs_df, chain_ids=chain_ids, min_contact_prob=min_contact_prob)
#     idxs =  dict(zip(chain_ids, np.where(mask)))
#     shifts = {chain_id:_get_shift(chain_id) for chain_id in chain_ids}
    
#     if shift:
#         return idxs, {chain_id:idxs_ - shifts[chain_id] for chain_id, idxs_ in idxs.items()}
#     else:
#         return idxs
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
