import pandas as pd 
import numpy as np 
import re 
import io 
import os

class DeepTMHMMFile():
    fields = ['id', 'annotation', 'start', 'end']

    def __init__(self):
        pass 

    @staticmethod
    def _parse_3line(path:str) -> pd.DataFrame:
        '''Parse a 3line output file. These files are in a FASTA-like format.'''
        with open(path, 'r') as f:
            lines = f.readlines()
            lines = [line for line in lines if (not line.startswith('#'))] # Remove comments. 
        
        assert (len(lines) % 3) == 0, f'DeepTMHMMFile._parse_3line: The number of lines in the 3line file should be a multiple of three {lines}.'
        headers = [lines[i].replace('>', '') for i in range(0, len(lines), 3)]

        ids = [header.split('|')[0].strip() for header in headers]
        df = pd.DataFrame(index=pd.Series(ids, name='id'))
        df['topology_type'] = [header.split('|')[1].strip() for header in headers] 
        df['topology_string_full'] = [lines[i + 2].replace('>', '') for i in range(0, len(lines), 3)] # Topologies are on the second line. 
        return df

    @staticmethod
    def _parse_gff3(path:str) -> pd.DataFrame:

        with open(path, 'r') as f:
            lines = f.readlines()

        keep = lambda line : not (line.startswith('//') or line.startswith('#'))
        text = '\n'.join([line for line in lines if keep(line)])

        df = pd.read_csv(io.StringIO(text), sep=r'\s+', names=DeepTMHMMFile.fields).set_index('id')

        topology_codes = {'signal':'S', 'TMhelix':'M', 'outside':'O', 'inside':'I'}
        get_topology_string = lambda df : ''.join([topology_codes.get(annotation, 'X') for annotation in df.annotation])

        df['topology_string'] = df.index.map(df.groupby(df.index).apply(get_topology_string))

        return df



    @classmethod
    def from_file(cls, path):
        
        df = DeepTMHMMFile._parse_gff3(path)

        if os.path.exists(path.replace('.gff3', '.3line')):
            # print(f'DeepTMHMMFile: Found 3line file at {path.replace('.gff3', '.3line')}')
            df = df.merge(DeepTMHMMFile._parse_3line(path.replace('.gff3', '.3line')), left_index=True, right_index=True)
        
        obj = cls()
        obj.df = df
        return obj 

    def to_df(self):
        return self.df.copy()



class TMHMMFile():

    patterns = dict()
    patterns['num_tmhs'] = r'# (?P<gene_id>[^\s]+) Number of predicted TMHs:\s+(?P<value>\d+)'
    patterns['prob_n_inside'] = r'# (?P<gene_id>[^\s]+) Total prob of N-in:\s+(?P<value>[\d\.]+)'
    # "Total prob of N-in" is the posterior probability that the N-terminus of the protein is on the cytoplasmic side of the membrane.
    patterns['num_aa_in_tmhs'] = r'# (?P<gene_id>[^\s]+) Exp number of AAs in TMHs:\s+(?P<value>[\d\.]+)'
    # If the first TMH could also be interpreted as a signal peptide, there is a POSSIBLE N-term signal sequence flag. 

    fields = ['gene_id', 'version', 'annotation', 'start', 'end']

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