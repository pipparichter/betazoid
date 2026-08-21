import pandas as pd 
import numpy as np 
import re


def chimerx_write_attrs(name, attrs_df:pd.DataFrame, path:str=None):
    assert 'residue_num' in attrs_df.columns, 'chimerx_write_attrs: attributes_df is missing the residue_num column.'
    assert 'chain_id' in attrs_df.columns, 'chimerx_write_attrs: attributes_df is missing the chain_id column.'
    assert name in attrs_df.columns, f'chimerx_write_attrs: attributes_df is missing the {name} column.'

    content = [f'attribute: {name}\nrecipient: residues\nmatch mode: 1-to-1']
    for row in attrs_df.itertuples():
        value = getattr(row, name)
        content += [f'\t/{row.chain_id}:{row.residue_num}\t{value}']

    with open(path, 'w') as f:
        f.write('\n'.join(content))