import numpy as np 
import pandas as pd 
from tqdm import tqdm
import re  
import os 
import matplotlib.pyplot as plt 
import networkx as nx
import seaborn as sns

get_group_id = lambda node_id : node_id.split(':')[0]
get_position = lambda node_id : int(node_id.split(':')[1])
get_edge_type = lambda contact_id : '-'.join(sorted([get_group_id(node_id) for node_id in contact_id.split('-')]))


# TODO: Should move this over to the AlphaFoldOutput object. 
def _get_residue(idx, token_chain_ids:np.ndarray):
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
    '''
    
    '''    
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
                row[f'chain_id_{i}'], residue_num = _get_residue(idx, token_chain_ids)

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
    if len(df) == 0:
        print(f'get_contacts: No contacts found for {output.name}')
    return pd.DataFrame(df)




def _get_pos(graph):
    '''Obtain positions for nodes in a contact graph. Positions are determined by the node ID format. Node IDs are of the format {group_id}:{n}, 
    where the group_1 determines the horizontal position of a node, and the n determines the vertical position.
    
    :param graph: A graph which is the output of get_contact_graph.
    :returns: A dictionary mapping each node ID in the input graph to an (x, y) position. 
    '''


    groups = dict()
    for node_id in list(graph.nodes):
        group_id = get_group_id(node_id)
        if group_id not in groups:
            groups[group_id] = list()
        groups[group_id].append(node_id)

    # If there is only one group, then use a circular layout. 
    if len(groups) == 1:
        return nx.circular_layout(graph)

    else:

        x_pos, y_pos = 0, 0
        pos = dict()
        for group_id, node_ids in groups.items():
            for node_id in sorted(node_ids, key=get_position): # Order according to the residue number
                pos[node_id] = (x_pos, y_pos)
                y_pos += 1
            x_pos += 1
            y_pos = 0 # Reset the y position for the next column. 

        return pos 


def get_contact_graph(contacts_df, edge_types:list=None, min_contact_prob:float=0.5):
    '''Build an undirected graph encoding the contacts provided in the contacts_df. 
    
    '''
    graph = nx.Graph()
    graph_df = contacts_df[contacts_df.contact_prob > min_contact_prob].copy()

    graph_df['edge_type'] = graph_df.contact_id.apply(get_edge_type)

    # for edge_type, df in graph_df.groupby('edge_type'):
        # print(f'get_contact_graph: Number of {edge_type} edges:', len(df))

    if edge_types is not None: # If edge types are specified, filter for those before constructing the graph. 
        graph_df = graph_df[graph_df.edge_type.isin(edge_types)].copy()

    node_ids = np.unique([node_id for node_ids in graph_df.contact_id for node_id in node_ids.split('-')])
    assert np.all([re.match(r'.+:\d+', node_id) for node_id in list(graph.nodes)]), f'get_contact_graph: Some of the node IDs are not in the correct format.'

    graph.add_nodes_from(node_ids)

    for contact_id in graph_df.contact_id:
        graph.add_edge(*contact_id.split('-'))
    return graph 


def get_contact_profile(contacts_df, group_id:str='cluster_3', edge_types=None, min_contact_prob:float=0, metric:str='num_contacts', max_position:int=None):
    
    assert metric in ['num_contacts', 'has_contact'], f'get_contact_profile: Specified metric {metric} is invalid.'
    graph = get_contact_graph(contacts_df, edge_types=edge_types, min_contact_prob=min_contact_prob)

    df = pd.DataFrame(graph.degree, columns=['node_id', 'num_contacts']) # Now each entry in the DataFrame is (1) a node and (2) the number of contacts the node participates in.
    df['has_contact'] = np.where(df.num_contacts > 0, 1, 0)
    df['position'] = df.node_id.apply(get_position) # n is the position encoded
    df = df[df.node_id.apply(get_group_id) == group_id].copy()

    assert len(df) > 0, 'get_contact_profile: No contacts to plot!'

    max_position = df.position.max() if (max_position is None) else max_position # Standardize profile lengths, if specified.
    df = pd.concat([df, pd.DataFrame({'position':[i for i in range(max_position) if i not in df.position.values]})]) # Get the positions with missing values. 
    df = df.fillna(0) # Fill in the missing degree values with zeroes. 
    df = df.sort_values('position')
    return df[metric].values


def plot_contact_profile(contacts_df, group_id:str='cluster_3', edge_types=None, ax:plt.Axes=None, x_min:int=0, x_max:int=100, min_contact_prob:float=0, metric:str='num_contacts', **kwargs):
    '''Construct a networkx Graph object using the contacts in the contacts_df.

    :param contacts_df:
    :param metric
    :returns None:
    '''

    figure_df = pd.DataFrame({metric:get_contact_profile(contacts_df, group_id=group_id, edge_types=edge_types, metric=metric, min_contact_prob=min_contact_prob)})
    figure_df = figure_df.reset_index(names='position', drop=False)

    if ax is None:
        fig, ax = plt.subplots(figsize=(0.1 * len(figure_df), 4))

    color = kwargs.get('color', 'gray')
    alpha = kwargs.get('alpha', 1)
    label = kwargs.get('label', None)
    legend = kwargs.get('legend', False)

    # sns.lineplot(figure_df, x='position', y='degree', color=color, ax=ax, alpha=alpha, edgecolor='black', linewidth=0.5)
    sns.lineplot(figure_df, x='position', y=metric, color=color, ax=ax, alpha=alpha, linewidth=1, label=label, legend=legend)
    ax.set_xticks(np.arange(x_min, x_max, 10), labels=np.arange(x_min, x_max, 10))
    ax.set_xlim(xmin=x_min, xmax=x_max)



def plot_contact_graph(contacts_df:pd.DataFrame, min_contact_prob:float=0.5, palette={'cluster_1':'lightgray', 'cluster_3':'gray'}):
    '''Construct a networkx Graph object using the contacts in the contacts_df, and plot the graph with a custom
    layout, where nodes are positioned in columns corresponding to groups (see get_pos)
    
    :param contacts_df:
    :param palette:
    :returns None:
    '''

    graph = get_contact_graph(contacts_df, min_contact_prob=min_contact_prob)
    node_colors = [palette.get(get_group_id(node_id), 'black') for node_id in graph.nodes] 
    node_ids = list(graph.nodes)

    fig, ax = plt.subplots(figsize=(10, 10))

    pos = _get_pos(graph)
    nx.draw_networkx_nodes(graph, pos=pos, node_color=node_colors, label=node_ids)
    nx.draw_networkx_edges(graph, pos=pos, edge_color='black')
    nx.draw_networkx_labels(graph, pos, labels=dict(zip(node_ids, node_ids)))

    for _, spine in ax.spines.items():
        spine.set_visible(False)

    plt.show()