
import pandas as pd 
import io 
import numpy as np 
import seaborn as sns 
import matplotlib.pyplot as plt
import os
import re

float_pattern = r'[-+]?\d+\.\d+'

class RNAFoldOutput():

    def _parse_output(self):

        with open(self.paths['output'], 'r') as f:
            text = f.read()

        nt_seq = text.split('\n')[1].replace('U', 'T')

        output = dict()
        output['mfe_frequency_in_ensemble'] = float(re.search(r'frequency of mfe structure in ensemble ([-+]?\d+\.\d+[e\-\d]*)', text).group(1))
        output['ensemble_diversity'] = float(re.search(f'ensemble diversity ({float_pattern})', text).group(1))
        output['mfe_free_energy'] = float(re.search(r'\(\s*([-+]?\d+\.\d+)\)', text).group(1))
        output['ensemble_free_energy'] = float(re.search(r'\[\s*([-+]?\d+\.\d+)\]', text).group(1))
        output['centroid_free_energy'] = float(re.search(r'\{\s*([-+]?\d+\.\d+)', text).group(1))
        output['centroid_distance_to_mfe'] = float(re.search(r'd=([-+]?\d+\.\d+)', text).group(1))
        return output, nt_seq

    def __init__(self, path:str):
        self.name = os.path.basename(path)
        self.dir_path = os.path.dirname(path)

        self.paths = dict()
        self.paths['dp'] = path + '_dp.ps' if os.path.exists(path + '_dp.ps') else None
        self.paths['ss'] = path + '_ss.ps' if os.path.exists(path + '_ss.ps') else None
        self.paths['output'] = path + '.out' if os.path.exists(path + '.out') else None

        self.output, self.nt_seq = self._parse_output()

    def get_data(self, nt_seq:bool=False):
        ''''''
        data = self.output.copy()
        data['name'] = self.name
        if nt_seq:
            data['nt_seq'] = self.nt_seq
        return data 
    
    def _parse_dotplot(self, as_array:bool=True, min_probability=0.8) -> pd.DataFrame:
        '''Read in the RNAFold output dot plot, which contains the base-pairing probabilities. This reflects the frequency that these bases are paired in the ensemble.
        
        :param path:
        '''
        path = self.paths['dp']
        with open(path, 'r') as f:
            lines = f.readlines()
        
        start_idx = np.where(np.array(lines) == '%start of base pair probability data\n')[0][0]
        end_idx = np.where(np.array(lines) == 'showpage\n')[0][0]
        text = ''.join(lines[start_idx + 1:end_idx])
        df = pd.read_csv(io.StringIO(text), sep=r'\s+', names=['base_idx_1', 'base_idx_2', 'probability_transformed', 'box'])
        df['probability'] = df.probability_transformed ** 2 # Undo the square-root transform.
        df = df[df.box == 'ubox'].copy()

        n = max(df.base_idx_1.max(), df.base_idx_2.max()) # Get the sequence length. 
        df = df[df.probability >= min_probability].copy()

        if as_array:
            arr = np.zeros((n, n))
            for row in df.itertuples():
                # Note that the bases are one-indexed, so need to shift by 1.
                arr[row.base_idx_1 - 1, row.base_idx_2 - 1] = row.probability
                arr[row.base_idx_2 - 1, row.base_idx_1 - 1] = row.probability
            return arr 

        return df

    # rnafold_df = pd.DataFrame([rnafold_from_file(path) for path in glob.glob(os.path.join(RNAFOLD_DIR, '*'))]).set_index('gene_id')


    def dotplot(self, start:int=None, stop:int=None, min_probability:float=0.8, ax:plt.Axes=None):

        arr = self._parse_dotplot(min_probability=min_probability, as_array=True)

        if ax is None:
            fig, ax = plt.subplots(figsize=(5, 5))

        start = 0 if (start is None) else start
        stop = -1 if (start is None) else stop

        arr = arr[start:stop, start:stop]
        sns.heatmap(arr, cbar=False, cmap='Grays', ax=ax)

        for _, spine in ax.spines.items():
            spine.set_visible(True)
        