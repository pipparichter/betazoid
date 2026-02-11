import pandas as pd 
import re 
import numpy as np 
import os 

def parse_enzyme_commission_code(definition:str):
    codes = re.search(r'\[EC:(.*)\]', definition, flags=re.DOTALL)
    if codes is not None:
        codes = codes.group(0).split()
    else:
        return {f'ec_{i + 1}':['-'] for i in range(4)}
    
    parsed_codes = {'ec_1':[], 'ec_2':[], 'ec_3':[], 'ec_4':[]}
    for code in codes: # Some enzymes have multiple assignments. 
        pattern = r'(\d+|-).(\d+|-).(\d+|-).(\d+|-)'
        match = re.search(pattern, code)
        if match is None:
            print(code)
        for i in range(1, 5):
            parsed_codes[f'ec_{i}'].append(match.group(i))
    return parsed_codes 

class KofamscanFile():


    def __init__(self):
        pass 

    @classmethod
    def from_file(cls, path:str):
        cols = ['gene_id', 'ko', 'threshold', 'score', 'e_value'] # , 'definition']
        # pattern = r'\s+(\d+_\d+)\s+(K\d{5})\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+"(.+)"'
        definition_pattern = r'"(.+)"'

        df = list()
        with open(path, 'r') as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue 
                line = line.replace('*', '')
                line = line.split()
                row = dict(zip(cols, line[:len(cols)]))
                row['definition'] = ' '.join(line[len(cols):])
                df += [row]

        df = pd.DataFrame(df)
        df['e_value'] = pd.to_numeric(df.e_value)
        df['score'] = pd.to_numeric(df.score)
        df['threshold'] = pd.to_numeric(df.threshold)

        obj = cls()
        obj.df = df  
        return obj

    def to_df(self, parse_ecs:bool=True):
        df = self.df.copy()
        if parse_ecs:
            # ec_df = pd.DataFrame([parsed_code for definition in df.definition for parsed_code in parse_enzyme_commission_code(definition)])
            ec_df = pd.DataFrame([parse_enzyme_commission_code(definition) for definition in df.definition])
            df = pd.concat([df, ec_df], axis=1)
            df = df.explode(['ec_1', 'ec_2', 'ec_3', 'ec_4'])
            df['definition'] = df.definition.str.replace(r'\[.*\]', '', regex=True)
        return df

