import os 
import re 
import pandas as pd 
import io 
from Bio import SeqIO
import glob
import numpy as np 
from fasta import FASTAFile

# ! muscle -align ../data/genes/cluster_1.fa -output ../data/genes/cluster_1.afa
class MSAFile():
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
        ids, arr = df.index.values, np.array([list(seq.replace('-', MSAFile.gap_symbol)) for seq in df.seq])
        return MSAFile.from_array(arr, ids) 

    @classmethod
    def from_string(cls, content:str, af3:bool=True):
        f = io.StringIO(content)
        records = [record for record in SeqIO.parse(f, 'fasta')]
        ids = [record.id for record in records]

        seqs = [str(record.seq).replace('-', MSAFile.gap_symbol) for record in records]
        if af3: # AF3 format is different, contains lowercase characters to indicate insertions relative to the query.
            seqs = [re.sub(r'[a-z]', '', seq) for seq in seqs]
        n_cols, n_rows = len(seqs[0]), len(seqs)
        arr = np.array([list(seq) for seq in seqs])
        return cls(arr, ids)

    def __len__(self):
        return len(self.ids)
    
    def __getitem__(self, id_):
        assert id_ in self.ids, f'MSAFile.__getitem__: ID {id_} is missing in the MSAFile.'
        return self.seqs[self.ids == id_][0]
    
    def to_array(self, alphabet:dict=None):
        '''Convert the MSAFile loaded from the FASTA file into a two-dimensional numpy array, where each entry is a single residue.'''
        return np.vectorize(alphabet.get)(self.arr.copy()) if (alphabet is not None) else self.arr.copy()
    
    def to_df(self, alphabet:dict=None):
        df = pd.DataFrame(index=self.ids)
        df['seq'] = self.seqs
        df['seq'] = df.seq.str.replace(alphabet) if (alphabet is not None) else df['seq']
        return df        
    
    def map_idx_from_msa(self, idx, gene_id:str):
        '''Convert the index of a residue in the MSAFile to the index of a residue in one of the aligned sequences.'''
        seq = self[gene_id]
        assert seq[idx] != MSAFile.gap_symbol, f'get_idx: The input index corresponds to a gap in the aligned {gene_id}.'
        n_gaps = seq[:idx].count(MSAFile.gap_symbol) # Get the number of gaps which occur before the requested index.
        return idx - n_gaps
    
    def map_idx_to_msa(self, idx, gene_id:str):
        '''Convert the index of a residue in one of the sequences to an index in the MSAFile.'''
        n, n_gaps = 0, 0
        for aa in self[gene_id]:
            if (n == idx) and (aa != MSAFile.gap_symbol):
                break # So that it doesn't exit if idx = 0 and the first symbol is a gap. 
            n += int(aa != MSAFile.gap_symbol)
            n_gaps += int(aa == MSAFile.gap_symbol)
        return n + n_gaps
    
    @classmethod
    def from_array(cls, arr, ids=None):
        ids = np.arange(len(arr)) if (ids is None) else ids
        arr = np.array(arr)
        ids = np.array(ids)
        return MSAFile(arr, ids)

    def get_mean_gap_fraction(self):
        return np.mean([np.mean(col == MSAFile.gap_symbol) for col in self.arr.T])


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
            row_map_idxs[np.where(row != MSAFile.gap_symbol)[0]] = np.arange((row != MSAFile.gap_symbol).sum()) # Fill in the values with the ungapped index positions. 
            row_map_idxs[np.where(row == MSAFile.gap_symbol)[0]] = -1
            map_idxs[id_] = row_map_idxs
        return map_idxs


