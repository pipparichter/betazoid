import pandas as pd 
import numpy as np 
import re 


class TMHMMFile():

    patterns = dict()
    patterns['num_tmhs'] = r'# (?P<gene_id>[^\s]+) Number of predicted TMHs:\s+(?P<value>\d+)'
    patterns['prob_n_inside'] = r'# (?P<gene_id>[^\s]+) Total prob of N-in:\s+(?P<value>[\d\.]+)'
    # "Total prob of N-in" is the posterior probability that the N-terminus of the protein is on the cytoplasmic side of the membrane.
    patterns['num_aa_in_tmhs'] = r'# (?P<gene_id>[^\s]+) Exp number of AAs in TMHs:\s+(?P<value>[\d\.]+)'
    # If the first TMH could also be interpreted as a signal peptide, there is a POSSIBLE N-term signal sequence flag. 

    fields = ['gene_id', 'version', 'location', 'start', 'end']

    def __init__(self, path:str=None):
        self.path = path  
        with open(path, 'r') as f:
            self.content = f.read() 
    
    @classmethod
    def from_file(cls, path):
        obj = cls(path=path)
        return obj 
    
    def get_num_tmhs(self):
        pattern = TMHMMFile.patterns['num_tmhs']
        num_tmhs = {match_.group(1):int(match_.group(2)) for match_ in re.finditer(pattern, self.content)}
        return num_tmhs
    
    def to_df(self):
        df = list()
        
        for line in self.content.split('\n'):
            for field, pattern in TMHMMFile.patterns.items():
                match = re.search(pattern, line)
                if match is not None:
                    info = match.groupdict()
                    info['field'] = field
                    df.append(info)
                    
        df = pd.DataFrame(df)
        df = df.pivot(index='gene_id', columns='field', values='value')
        df.columns.name = ''
        assert np.all(df.index.value_counts() == 1), f'load_tmhmm: There should be only one metadata entry per gene.'

        df_ = pd.read_csv(self.path, comment='#', sep=r'\s+', names=TMHMMFile.fields).set_index('gene_id')
        df = df_.merge(df, left_index=True, right_index=True, how='left')
        df = df.astype({'num_tmhs':int})

        return df.reset_index(names='gene_id')