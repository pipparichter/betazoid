import sys 
sys.path.append('/home/prichter/Documents/banfield/betazoid/src/files/')

from src.files.fasta import * 
from src.files.gfa import * 
from src.files.bam import BamFile

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
from src.files import FASTAFile
from src.files.fasta import * 
import subprocess
from src.files.kofamscan import KofamscanFile
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LinearSegmentedColormap, to_hex
from Bio.Seq import Seq


PROJECT_IDS = pd.read_csv('project_ids.csv', index_col=0).project_id.to_dict()
FORWARD_READS_PATHS = pd.read_csv('reads_paths.csv', index_col=0).forward_reads_path.to_dict()
REVERSE_READS_PATHS = pd.read_csv('reads_paths.csv', index_col=0).reverse_reads_path.to_dict()

reverse_complement = lambda seq : str(Seq(seq).reverse_complement())


def load_seqs(ids:list, path:str=None):
    '''Load the specified sequences from the FASTA file at path. Output is a dictionary mapping the IDs to the sequences.'''
    ids = set(ids) # Remove duplicate contig IDs.
    fasta_file = FASTAFile.from_file(path, filter_=lambda record : record.id in ids)
    seqs = dict(zip(fasta_file.ids, fasta_file.seqs))
    assert len(seqs) == len(ids), f'load_seqs: {len(ids)} sequence IDs provided, but only {len(ids)} loaded from {path}'
    return seqs 


def get_bbmap_command(ref_path, output_dir:str=None, forward_reads_path:str=None, reverse_reads_path:str=None, verbose:bool=True, **_params):
    '''Generate a command for mapping reads files to a reference contig. Returns the bbmap command, as well as the output path
    for the generated BAM file. 
    
    :param ref_path : 
    :param output_dir : 
    :param forward_reads_path : 
    :param reverse_reads_path
    '''
    ref_id = os.path.basename(ref_path).replace('.fasta', '')
    output_path = os.path.join(output_dir, f'{ref_id}.bam')

    # Define the default parameters. 
    params = dict()
    params['local'] = 'f'
    params['pairedonly'] = 'f'
    params['ambiguous'] = 'random'
    params['minid'] = 0.9
    params['idfilter'] = 0.95
    params['editfilter'] = -1
    params.update(_params)

    if verbose:
        print('get_bbmap_command: Using the following mapping settings:')
        for param, value in params.items():
            print(f'\t{param}: {value}')

    params['out'] = 'stdout.sam'
    params['threads'] = 64
    params['pigz'] = 't'
    params['unpigz'] = 't'

    params = ' '.join([f'{param}={value}' for param, value in params.items()])
    cmd = f'bbmap.sh {params} in1={forward_reads_path} in2={reverse_reads_path} ref={ref_path} nodisk | shrinksam | sambam > {output_path}'
    return cmd, output_path



# # TODO: Speed this up with https://pypi.org/project/pyfaidx/.
# def get_sequences(results_df:pd.DataFrame, path:str='../data/databases/nantong_groundwater/contigs.fa'):
#     seqs = dict()
#     pattern = re.compile(rf"({'|'.join(results_df.index.tolist())})")

#     seq, contig_id, read = [], None, False
#     pbar = tqdm(total=len(results_df), desc='get_sequences')

#     f = open(path, 'r')
#     for line in f:
#         if (not read):
#             match_ = re.search(pattern, line)
#             if match_:
#                 read = True 
#                 contig_id = match_.group(1)
#         elif read:
#             if line.startswith('>'):
#                 seqs[contig_id] = ''.join(seq) 
#                 seq = []
#                 read = False
#                 pbar.update(1)
#             else:
#                 seq.append(line.strip())
#     results_df['seq'] = results_df.index.map(seqs)
#     return results_df 