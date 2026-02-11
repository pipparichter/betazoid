import pandas as pd 
import re 
import numpy as np 
import os 


class InterProScanFileTSV():

    fields = ['id', 'md5', 'length', 'signature_analysis', 'signature_accession', 'signature_description', 'start', 'stop', 'e_value', 'status', 'date', 'interpro_accession', 'interpro_description']

    def __init__(self):
        pass 

    @classmethod
    def from_file(cls, path:str):

        df = pd.read_csv(path, usecols=list(range(len(InterProScanFileTSV.fields))), names=InterProScanFileTSV.fields, header=None, sep='\t')
        # df = df.set_index('id')
        
        obj = cls()
        obj.df = df 
        return obj 

    def to_df(self):
        df = self.df.copy()
        df['e_value'] = df.e_value.replace('-', np.nan)
        df['e_value'] = pd.to_numeric(df.e_value)
        return df

