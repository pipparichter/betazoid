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


# chimerx_commands = ['color gray']
# chimerx_commands += ['sel :lys; color sel blue']
# chimerx_commands += ['sel :arg; color sel red']
# chimerx_commands += ['sel :tyr; color sel green']
# chimerx_commands += ['sel :trp; color sel yellow']
# chimerx_commands += ['sel :phe; color sel orange']

# print('; '.join(chimerx_commands))

# idxs = [139]
# chains = 'FC'
# atom = 'NZ'
# chimerx_commands = list()
# for idx in idxs:
#     chimerx_commands += [f'select /{chains[0]}:{idx}@{atom} /{chains[1]}:{idx}@{atom}; distance sel']
# chimerx_commands += ['color white pseudobonds; color white label']
# print('; '.join(chimerx_commands))
# # select /B:224@NZ /E:224@NZ
# # distance sel 

# # color white pseudobonds 
# # color white label

# # Oddly, the residues which are conserved do not appear to be the ones lining the channel (based on examination of the structures). 
# # It seems possible that the structures are off, or possibly adopt different conformations. Might be worth folding with dsDNA and ssDNA and then look at 
# # which residues are closest to the nucleic acid through the pore. 

# # Some of the pentameric configurations (e.g. bz_3 bound to Mg and ATP) have a pore size that seems too small to be biologically-relevant
# # (like 6 angstroms). This would be on the small side, even if the substrate were ssDNA. 