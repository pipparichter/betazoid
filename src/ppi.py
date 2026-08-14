import numpy as np 
import pandas as pd 
from tqdm import tqdm 
import os 

def get_residue(idx, token_chain_ids:np.ndarray):
    '''

    :param idx: The zero-indexed contact location relative to the token_chain_ids array. 
    :param token_chain_ids: The token_chain_ids array from the AlphaFold output. 
    :returns: A two-tuple containing (1) the chain ID and (2) zero-indexed residue number relative to the start of the chain. 
    '''
    chain_start_idxs = dict(zip(*np.unique(token_chain_ids, return_index=True))) # This is a dictionary mapping the chain IDs to their start indices.
    chain_id = token_chain_ids[idx] # Get the chain that the specified index belongs to. 
    return (chain_id, idx - chain_start_idxs[chain_id])


def get_contact_idxs(contact_probs:np.ndarray, token_chain_ids:np.ndarray, min_contact_prob=0.3):
    '''Use the contact_probs in the AlphaFold output to locate potential inter-chain residue contacts.
    
    :param contact_probs: A square numpy array the contact_probs for a given structure. 
    :param idx: The zero-indexed contact location relative to the token_chain_ids array. 
    :param min_contact_prob: The minimum contact_prob to say whether or not two residues are likely to be in contact. 
    :returns: A list of two-tuples containing the contact_probs_df indices where contacts occur.
    '''
    assert isinstance(contact_probs, np.ndarray), f'get_contact_idxs: Expected contact_probs to be a numpy array, but got {type(contact_probs)}'
    mask = (contact_probs > min_contact_prob) # Require a minimum contact probability. 
    mask = mask & (np.expand_dims(token_chain_ids, axis=1) != np.expand_dims(token_chain_ids, axis=0))  # Don't include intra-chain contacts. 
    mask[np.tril_indices_from(mask, k=-1)] = False # Don't include duplicate contacts, as the matrix is symmetric.
    return list(zip(*np.where(mask))) # Gets the contact indices as a list of two-tuples (i1, j1), (12, j2), ...


# NOTE: Everything here should be zero-indexed. 
def get_contacts(output, min_contact_prob:float=0.5):

    contact_probs = output.get_contact_probs(mean_pool=False, best_model=False) # Get the contact_probs as a dictionary mapping the model name to the contact_probs_df. 
    paes = output.get_paes(mean_pool=False, best_model=False)
    token_chain_ids = output.get_token_chain_ids()
    chains = output.get_chains()

    df = list()
    for model in contact_probs.keys():
        contact_idxs = get_contact_idxs(contact_probs[model], token_chain_ids=token_chain_ids, min_contact_prob=min_contact_prob)

        for idxs in contact_idxs:
            row = dict()
            for i, idx in enumerate(idxs):
                row[f'idx_{i}'] = idx 
                row[f'chain_id_{i}'], residue_num = get_residue(idx, token_chain_ids)

                seq = chains.get(row[f'chain_id_{i}'], None) # Only the protein chains are in the chains dictionary, but there could be contacts defined with ligands.
                if seq is not None:
                    row[f'residue_num_{i}'] = residue_num
                    row[f'residue_{i}'] = seq[residue_num]
                    row[f'seq_{i}'] = seq # Need this to match the chain with the gene ID (not a better way to do this without chain descriptions). 

            row['pae'] = paes[model][idxs[0], idxs[1]] # PAE is not symmetric. 
            row['contact_prob'] = contact_probs[model][idxs[0], idxs[1]] # This should be symmetric. 
            row['name'] = output.name 
            row['model'] = model

            df.append(row)
    return pd.DataFrame(df)