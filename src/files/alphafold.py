import pandas as pd 
import numpy as np 
import os 
import glob 
import json
import re
from fasta import FASTAFile
import itertools 


CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# https://github.com/google-deepmind/alphafold3/blob/main/docs/input
class AlphaFoldInputFile():
    '''Generate a JSON input file for the AlphaFold server. This is based on the job_request.json file which is included in the data downloaded
    from the online AlphaFold server.'''

    # def __init__(self, name:str, dialect='alphafoldserver', version:int=1, num_seeds:int=5):
    def __init__(self, name:str, dialect='alphafold3', version:int=2, num_seeds:int=5):
        self.name = name 
        
        self.info = {'dialect':dialect, 'version':version, 'modelSeeds':[seed for seed in range(1, num_seeds + 1)]} 
        self.info['name'] = name
        self.info['sequences'] = list()

        self.dialect = dialect

        self.n = 0 # Stores the total number of entities in the input. 

        self.seq_types = dict()
        self.seq_types['protein'] = 'proteinChain' if (dialect == 'alphafoldserver') else 'protein'
        self.seq_types['dna'] = 'dnaSequence' if (dialect == 'alphafoldserver') else 'dna'

        self.ligand_types = dict()
        self.ligand_types['atp'] = ('ligand', 'CCD_ATP')
        self.ligand_types['adp'] = ('ligand', 'CCD_ADP')
        self.ligand_types['mg'] = ('ion', 'MG') if (dialect == 'alphafoldserver') else ('ligand', 'MG')

        self.ligands = {ligand_type:0 for ligand_type in self.ligand_types.keys()}

    def add_seq(self, seq:str, n:int=1, type_:str='protein', paired_msa:str=None, unpaired_msa:str=None):
        ''' '''
        assert type_ in ['protein', 'dna'], f'AlphaFoldInputFile.add_seq: Sequence type must be either dna or protein.'
        assert not ((paired_msa is None) ^ (unpaired_msa is None)), f'AlphaFoldInputFile.add_seq: Both unpairedMsa and pairedMsa have to either be both set (i.e. non-null), or both unset (i.e. both null, explicitly or implicitly).'

        info = {'sequence':seq}
        if self.dialect == 'alphafoldserver':
            info['count'] = n 
        elif self.dialect == 'alphafold3':
            info['id'] = self._get_chain_ids(n)

        if (paired_msa is not None) and (unpaired_msa is not None):
            info['pairedMsa'] = paired_msa
            info['unpairedMsa'] = unpaired_msa

            query_seqs = np.array([paired_msa.split('\n')[1], unpaired_msa.split('\n')[1]])
            assert np.all(seq == query_seqs), 'AlphaFoldInputFile.add_seq: The sequence in the first line of each MSA must match the query.'

        self.info['sequences'] += [{self.seq_types[type_]:info}]
        
    
    def _get_chain_ids(self, n):
        assert (self.n + n) <= len(CHAIN_IDS), 'AlphaFoldInput._get_chain_ids: Not enough available chain IDs left!'
        chain_ids = [CHAIN_IDS[i] for i in range(self.n, self.n + n)]
        self.n += n 
        return chain_ids 

    def _add_ligand(self, type_:str=None, n:int=1):
        category, symbol = self.ligand_types[type_]
        self.ligands[type_] += n

        if self.dialect == 'alphafoldserver':
            self.info['sequences'] += [{category:{category:symbol, 'count':n}}]
        elif self.dialect == 'alphafold3':
            self.info['sequences'] += [{category:{'ccdCodes':[symbol], 'id':self._get_chain_ids(n)}}]
        else:
            raise Exception(f'AlphaFoldInput._add_ligand: Dialect {self.dialect} is not recognized.')
        
    def add_atp(self, n:int=1):
        self._add_ligand(type_='atp', n=n)

    def add_adp(self, n:int=1):
        self._add_ligand(type_='adp', n=n)

    def add_mg(self, n:int=1):
        self._add_ligand(type_='mg', n=n)
        
    def write(self, path:str):
        with open(path, 'w') as f:
            json.dump(self.info, f)

    def get_info(self):
        return [self.info]



class AlphaFoldOutput():
    patterns = dict()
    patterns['summary'] = r'(seed-\d+_sample-\d+)\/summary_confidences'
    patterns['full'] = r'(seed-\d+_sample-\d+)\/confidences'


    @staticmethod
    def _get_confidences(dir_path:str):
        '''Recursively search the AlphaFold output directory for the data paths containing the full and summary model confidence data.
        These output files are then organized into models according to the sub-model they correspond to. 
        The JSON file containing the full confidence data includes:
        (1)
        The JSON file containing the summary confidence data includes:
        (1)

        : param dir_path: The root directory for the AlphaFold prediction. 
        '''
        confidences = dict()

        for path in glob.glob(os.path.join(dir_path, '**', '*')):

            for file_type, pattern in AlphaFoldOutput.patterns.items():
                if re.search(pattern, path) is None:
                    continue 
                model = re.search(pattern, path).group(1)

                if model in confidences:
                    confidences[model][file_type] = path 
                else:
                    confidences[model] = {file_type:path}

        return confidences 


    def __init__(self, dir_path=None):
        # Paths to the confidence output for each model. 
        self.confidences = AlphaFoldOutput._get_confidences(dir_path)
        self.name = os.path.basename(dir_path)

        with open(os.path.join(dir_path, f'{self.name}_data.json'), 'r') as f:
            self.data = json.load(f)

        get_chain_id = lambda entry : [info['id'] for info in entry.values()][0]
        self.chain_ids = [get_chain_id(entry) for entry in self.data['sequences']]
        self.num_chains = len(self.chain_ids)


    def get_paes(self):
        paes = dict()
        for model, paths in self.confidences.items():
            with open(paths['full'], 'r') as f:
                data = json.load(f)
                paes[model] = pd.DataFrame(data['pae'], index=data['token_chain_ids'], columns=data['token_chain_ids'])
        return paes
    
    
    def get_iptms(self):
        iptms = dict()

        for model, paths in self.confidences.items():
            with open(paths['summary'], 'r') as f:
                try:
                    iptms[model] = json.load(f)['iptm']
                except:
                    print(f'AlphaFoldOutput.get_iptms: Could not decode {f.read()}')
        return iptms
    
    def get_ptms(self):
        ptms = dict()
        for model, paths in self.confidences.items():
            with open(paths['summary'], 'r') as f:
                try:
                    ptms[model] = json.load(f)['ptm']
                except:
                    print(f'AlphaFoldOutput.get_ptms: Could not decode {f.read()}')
        return ptms

    def get_protein_chain_ids(self):
        '''Obtain the chains corresponding to actual protein sequences (not ligands or DNA) using the data.json file.'''
        protein_chain_ids = list()
        for entry in self.data['sequences']:
            if 'protein' in entry:
                protein_chain_ids += entry['protein']['id']
        return protein_chain_ids

    def get_msas(self):

        msas = list()

        for entry in self.data['sequences']:
            if 'protein' in entry:
                msa = dict()
                msa['paired'] = entry['protein'].get('pairedMsa', '')
                msa['unpaired'] = entry['protein'].get('unpairedMsa', '')
                msas.append(msa)

        return msas






class AlphaFoldServerOutput(AlphaFoldOutput):

    patterns = dict()
    patterns['summary'] = r'.*summary_confidences_(\d+)'
    patterns['full'] = r'full_data_(\d+)'
    
    
    def __init__(self, dir_path=None):
        # Paths to the confidence output for each model. 
        self.confidences = AlphaFoldServerOutput._get_confidences(dir_path)
        self.name = os.path.basename(dir_path)
        self.dir_path = dir_path

        with open(os.path.join(dir_path, f'{self.name}_job_request.json'), 'r') as f:
            self.data = json.load(f)[0]

        get_num_chains = lambda entry : [info['count'] for info in entry.values()][0]
        self.num_chains = sum([get_num_chains(entry) for entry in self.data['sequences']])
        self.chain_ids = [CHAIN_IDS[i] for i in range(self.num_chains)]

    def get_protein_chain_ids(self):
        '''Obtain the chains corresponding to actual protein sequences (not ligands or DNA) using the data.json file.'''
        
        protein_chain_idxs = list()
        i = 0

        for entry in self.data['sequences']:
            type_ = list(entry.keys())[0]
            n = entry[type_]['count']

            if type_ == 'protein':
                protein_chain_idxs += list(range(i, i + n))
            i += n
        return [CHAIN_IDS[i] for i in protein_chain_idxs]

    def get_msas(self):
        # MSA file names look like: fold_a7xxr1_bp742_1mer_mg_atp_paired_msa_chains_a_b_c.a3m
        msa_pattern = r'(\w+)_(unpaired|paired)_msa_chains_(\w+)\.a3m'
        msa_dir_path = os.path.join(self.dir_path, 'msas')

        # get_chain_ids = lambda path : [chain_id.upper() for chain_id in re.search(msa_pattern, path).group(3).split('_')]
        # chain_ids = set([chain_id for path in glob.glob(os.path.join(msa_dir_path, '*')) for chain_id in get_chain_ids(path)])
        # print(f'AlphaFoldServerOutput.get_msas: Found MSA files for {len(chain_ids)} protein chains, {' '.join(chain_ids)}')

        msas = dict()

        for path in glob.glob(os.path.join(msa_dir_path, '*')):
            match_ = re.search(msa_pattern, path)

            if match_.group(3) not in msas:
                msas[match_.group(3)] = dict()

            with open(path, 'r') as f:
                msa = f.read()
            
            msas[match_.group(3)][match_.group(2)] = msa 

        return list(msas.values()) # Make sure paired and unpaired for the same chains are associated here. 






    


def get_interfaces(self, min_contact_prob=0.7, pairwise=True, chain_ids=None):

    with open(path, 'r') as f:
        data = json.load(f)
        contact_probs = np.array(data['contact_probs'])
        token_chain_ids = np.array(data['token_chain_ids'])

    chain_ids = list(np.unique(token_chain_ids)) if (chain_ids is None) else chain_ids

    mask = (contact_probs > min_contact_prob) # Require a minimum contact probability. 
    mask = mask & (token_chain_ids.reshape(-1, 1) != token_chain_ids)  # Don't include intra-chain contacts. 
    mask[np.tril_indices_from(mask, k=-1)] = False # Don't double-count pairwise interactions (remove everything below the diagonal)

    if pairwise:
        masks = dict()
        for chain_id_pair in itertools.product(chain_ids, chain_ids):
            masks[chain_id_pair] = mask & (np.isin(token_chain_ids, chain_id_pair).reshape(-1, 1) & np.isin(token_chain_ids, chain_id_pair))
        return masks
    else:
        return mask

#     # Want to not consider intra-chain contacts. 

# def get_num_contacts(path, chain_ids:list=None, min_contact_prob=0.5):
#     df = pd.DataFrame()
#     for chain_id_pair, mask in get_interfaces(path, min_contact_prob=min_contact_prob, pairwise=True, chain_ids=chain_ids).items():
#         df.loc[chain_id_pair[0], chain_id_pair[1]] = mask.sum()
#         df.loc[chain_id_pair[1], chain_id_pair[0]] = mask.sum()
#     return df.fillna(0)


# def get_interface_pae(path, chain_ids:list=None, pairwise:bool=True, min_contact_prob:float=0.5, func=np.min):
#     '''Locates the interfaces between subunits using the contact_probs field and returns the PAE for the residues along the interface.
    
#     :param path 
#     :param chain_ids
#     :param pairwise
#     :param min_contact_prob
#     :param func
    
#     '''
#     with open(path, 'r') as f:
#         paes = np.array(json.load(f)['pae'])

#     if pairwise:
#         df = pd.DataFrame()
#         for chain_id_pair, mask in get_interfaces(path, min_contact_prob=min_contact_prob, pairwise=True, chain_ids=chain_ids).items():
#             if mask.sum() == 0:
#                 print(f'get_interface_pae: No interface found between chains {chain_id_pair[0]} and {chain_id_pair[1]}.')
#                 continue 
#             df.loc[chain_id_pair[0], chain_id_pair[1]] = func(paes[mask])
#             df.loc[chain_id_pair[1], chain_id_pair[0]] = func(paes[mask])
#         return df.fillna(0)

#     else:
#         mask = get_interfaces(path, min_contact_prob=min_contact_prob, pairwise=False, chain_ids=chain_ids)
#         return paes[mask].mean()
    

