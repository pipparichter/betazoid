import re 
import pandas as pd 
import os 
import glob
import numpy as np 
import seaborn as sns 
from scipy.stats import gmean 
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from src.files import FASTAFile
import requests 
import scipy 
import warnings 
import json
import itertools
import subprocess
from src.files.kofamscan import KofamscanFile
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LinearSegmentedColormap, to_hex
from Bio.Seq import Seq 

reverse_complement = lambda seq : str(Seq(seq).reverse_complement())

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