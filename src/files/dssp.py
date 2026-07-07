
import pandas as pd 
import numpy as np 

class DSSPFile():
    # RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
    # fields = ['residue', 'aa', 'structure', 'bp1', 'bp2', 'acc', 'n_h_o_1', 'o_h_n_1', 'n_h_o_2', 'o_h_n_2', 'tco', 'kappa', 'alpha', 'phi', 'psi', 'x_ca', 'y_ca', 'z_ca']
    code_map = dict()
    code_map['H'] = 'alpha helix'
    code_map['G'] = '3-10 helix' # Tighter than an alpha helix, sometimes occurs at helix boundaries. 
    code_map['I'] = 'pi helix' # Looser than an alpha helix. 
    code_map['E'] = 'beta strand'
    code_map['B'] = 'beta bridge'
    code_map[' '] = 'none'
    code_map['T'] = 'turn'
    code_map['S'] = 'bend'
    
    def __init__(self):
        pass 
    
    @classmethod
    def from_file(cls, path):
        '''DSSP files are fixed-width, meaning that splitting them along certain characters will not work. '''

        df = list()
        with open(path, 'r') as f:

            line = f.readline()

            is_in_header = lambda line : (line.strip()[0] != '#')
            while is_in_header(line):
                line = f.readline()
            
            while line:
                line = f.readline()
                try:
                    row = dict()
                    row['residue'] = int(line[5:11])
                    row['chain'] = line[11]
                    row['aa'] = line[13]
                    row['code'] = line[16]
                    df.append(row)
                except Exception as err:
                    continue
                    # print(f'DSSPFile.from_file: Failed to parse line, {err}')
        df = pd.DataFrame(df)

        obj = cls()
        obj.df = df 
        return obj
    
    def __str__(self):
        codes = self.df.code.tolist()
        return ''.join(codes)

    def to_df(self):
        df = self.df.copy()
        df['structure'] = df.code.map(DSSPFile.code_map)
        return df
    
 