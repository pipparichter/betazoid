import os 
import re 
import pandas as pd 
import io 
from Bio import SeqIO
import glob
import numpy as np 
from files.fasta import FASTAFile


def load_clipkit_log(path:str):
    '''Load the clipkit log file stored at the specified path.'''
    fields = ['idx', 'status', 'informative', 'gap_fraction']
    df = pd.read_csv(path, sep=r' ', names=fields)
    df['idx'] = df['idx'] - 1 # Make sure this is zero-indexed. 
    return df

    # If I want to map the MSA index back to the original sequence, I need to account for the dropped columns. 
    # For example, if I am looking for the location of residue i in a sequence in the trimmed MSA, but n columns prior to 
    # i have been removed, then the mapped index will be n positions too far to the right. 

    # Actually, there is not enough information in the log file to map from the sequence to the MSA, because we don't know if the removed 
    # column had a gap in the target sequence or was filled.  

    # I will need to load the original file, get the index in that file, and then determine if the column was retained. 



# ! muscle -align ../data/genes/cluster_1.fa -output ../data/genes/cluster_1.afa
class MSAFile():
    gap_symbol = '.'
    def __init__(self, arr, ids):

        self.ids = ids 
        self.arr = arr 
        self.seqs = np.array([''.join(row) for row in arr])
        self.n_cols = arr.shape[-1]

    @classmethod
    def from_file(cls, path:str='../data/genes/cluster_1.afa', clipkit_log_path:str=None):

        df = FASTAFile.from_file(path).to_df()
        # print(df)
        ids, arr = df.index.values, np.array([list(seq.replace('-', MSAFile.gap_symbol)) for seq in df.seq])
        obj = MSAFile.from_array(arr, ids) 

        if clipkit_log_path is not None:
            obj.clipkit_log_df = load_clipkit_log(clipkit_log_path)

        return obj

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

    
    def map_idx_to_msa(self, idx, gene_id: str):
        '''Convert the index of a residue in one of the sequences to an index in the MSAFile.'''
        seq_idx = 0

        for msa_idx, aa in enumerate(self[gene_id]):
            if aa != MSAFile.gap_symbol:
                if seq_idx == idx:
                    return msa_idx
                seq_idx += 1 # Increment the sequence index only if there is not a gap symbol in the MSA.
        return msa_idx # Return if the index is the last residue of the gene_id sequence. 


    def map_idx_to_trimmed_msa(self, idx, gene_id):
        ''''''
        idx = self.map_idx_to_msa(idx, gene_id)
        original_msa_idxs = self.clipkit_log_df[self.clipkit_log_df.status == 'keep'].idx.values 
        idx_map = dict(zip(original_msa_idxs, np.arange(len(original_msa_idxs))))
        return idx_map.get(idx, None)

    def map_idx_from_trimmed_msa(self, idx, gene_id):
        ''''''
        original_msa_idxs = self.clipkit_log_df[self.clipkit_log_df.status == 'keep'].idx.values 
        idx_map = dict(zip(np.arange(len(original_msa_idxs)), original_msa_idxs))
        idx = idx_map[idx]
        return self.map_idx_from_msa(idx, gene_id)

    
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


    def get_entropy(self, alphabet:dict=None, max_gap_fraction:float=0.2):

        '''
        :param max_gap_fraction: If more than this fraction of columns has a gap, then return None instead of an entropy value. 
        '''

        alphabet_size = 20 if (alphabet is None) else len(np.unique(list(alphabet.values())))
        get_gap_fraction = lambda col : np.mean(col == self.gap_symbol)

        entropy = list()
        for col in self.to_array(alphabet=alphabet).T:
            if get_gap_fraction(col) > max_gap_fraction:
                entropy.append(np.nan)
                continue 
            _, counts = np.unique(col[col != self.gap_symbol], return_counts=True)
            frequencies = counts / counts.sum()
            entropy.append(sum(-(frequencies * np.log(frequencies) / np.log(alphabet_size))))
        return np.array(entropy)