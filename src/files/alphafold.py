import pandas as pd 
import numpy as np 
import os 
import glob 
import json
import orjson
import re
from fasta import FASTAFile
import itertools 

CHAIN_IDS = list()
for n in range(1, 3):
    for chain_id in itertools.product('ABCDEFGHIJKLMNOPQRSTUVWXYZ', repeat=n):
        CHAIN_IDS.append(''.join(chain_id[::-1]))


def load_json(path):
    '''Read a JSON file from the input path using the faster orjson parser. '''
    try:
        with open(path, 'rb') as f:
            data = f.read()
            data = orjson.loads(data)
        return data
    except Exception as err:
        print(f'load_json: Could not decode {path}, {err.msg}')
        print(data)
        return None


# https://github.com/google-deepmind/alphafold3/blob/main/docs/input
class AlphaFoldInputFile():
    '''Generate a JSON input file for the AlphaFold server. This is based on the job_request.json file which is included in the data downloaded
    from the online AlphaFold server.'''

    # def __init__(self, name:str, dialect='alphafoldserver', version:int=1, num_seeds:int=5):
    def __init__(self, name:str, dialect='alphafold3', version:int=4, num_seeds:int=2):
        self.name = name 
        self.version = version

        self.info = {'dialect':dialect, 'version':version, 'modelSeeds':[seed for seed in range(1, num_seeds + 1)]} 
        self.info['name'] = name
        self.info['sequences'] = list()

        self.dialect = dialect

        self.n = 0 # Stores the total number of entities in the input. 

        self.seq_types = dict()
        self.seq_types['protein'] = 'proteinChain' if (dialect == 'alphafoldserver') else 'protein'
        self.seq_types['dna'] = 'dnaSequence' if (dialect == 'alphafoldserver') else 'dna'

        self.ligand_types = dict()
        self.ligand_types['atp'] = ('ligand', 'CCD_ATP') if (dialect == 'alphafoldserver') else ('ligand', 'ATP')
        self.ligand_types['adp'] = ('ligand', 'CCD_ADP') if (dialect == 'alphafoldserver') else ('ligand', 'ADP')
        self.ligand_types['mg'] = ('ion', 'MG') if (dialect == 'alphafoldserver') else ('ligand', 'MG')

        self.ligands = {ligand_type:0 for ligand_type in self.ligand_types.keys()}

    def add_seq(self, seq:str, n:int=1, type_:str='protein', paired_msa:str=None, unpaired_msa:str=None, description:str=None):
        ''' '''
        assert type_ in ['protein', 'dna'], f'AlphaFoldInputFile.add_seq: Sequence type must be either dna or protein.'
        # assert not ((paired_msa is None) ^ (unpaired_msa is None)), f'AlphaFoldInputFile.add_seq: Both unpairedMsa and pairedMsa have to either be both set (i.e. non-null), or both unset (i.e. both null, explicitly or implicitly).'

        info = {'sequence':seq}
        if self.dialect == 'alphafoldserver':
            info['count'] = n 
        elif self.dialect == 'alphafold3':
            info['id'] = self._get_chain_ids(n)


        if (unpaired_msa is not None):
            info['unpairedMsa'] = unpaired_msa
            info['pairedMsa'] = ''
        if (paired_msa is not None):
            info['pairedMsa'] = paired_msa
        if (description is not None) and (self.version == 4):
            info['description'] = description

            # query_seqs = np.array([paired_msa.split('\n')[1], unpaired_msa.split('\n')[1]])
            # assert np.all(seq == query_seqs), 'AlphaFoldInputFile.add_seq: The sequence in the first line of each MSA must match the query.'

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
        info = [self.info] if (self.dialect == 'alphafoldserver') else self.info
        with open(path, 'w') as f:
            json.dump(info, f)

    def get_info(self):
        return self.info.copy()



class AlphaFoldOutput():
    patterns = dict()
    patterns['summary'] = r'(seed-\d+_sample-\d+)\/summary_confidences'
    patterns['full'] = r'(seed-\d+_sample-\d+)\/confidences'

    def _get_best_model(self):
        '''Get the name of the model with the highest ranking score from the ranking_scores.csv file in the AlphaFold3 output root directory.
        
        :returns: The name of the model with the highest ranking score, which is of the form seed-x-sample-y. 
        '''
        df = pd.read_csv(os.path.join(self.dir_path, 'ranking_scores.csv')) # This file contains the model rankings. 
        df = df.sort_values('ranking_score', ascending=False)
        seed, sample = df.iloc[0]['seed'], df.iloc[0]['sample']
        return f'seed-{int(seed)}_sample-{int(sample)}' # The name of the model corresponding to the highest score. 

    def _get_confidences(self):
        '''Recursively search the AlphaFold output directory for the data paths containing the full and summary model confidence data.
        These output files are then organized into models according to the sub-model they correspond to. 
        The JSON file containing the full confidence data includes:
        (1)
        The JSON file containing the summary confidence data includes:
        (1)
        '''
        # print(f'AlphaFoldOutput._get_confidences: Using file patterns:')
        # for file_type, pattern in self.patterns.items():
        #     print(f'\t{file_type}\t{pattern}')

        confidences = dict()
        for path in glob.glob(os.path.join(self.dir_path, '**', '*'), recursive=True):

            for file_type, pattern in self.patterns.items():
                if re.search(pattern, path) is None:
                    continue 
                model = re.search(pattern, path).group(1)

                if model in confidences:
                    confidences[model][file_type] = os.path.abspath(path)
                else:
                    confidences[model] = {file_type:os.path.abspath(path)}
        assert len(confidences) > 0, f'AlphaFoldOutput._get_confidences: Was unable to find any confidence JSON files in {self.dir_path}.'
        return confidences 


    def __init__(self, dir_path=None):
        # Paths to the confidence output for each model. 
        self.dir_path = os.path.abspath(dir_path)
        self.confidences = self._get_confidences()
        self.name = os.path.basename(dir_path)

        self.data = load_json(os.path.join(dir_path, f'{self.name}_data.json'))

        get_chain_id = lambda entry : [info['id'] for info in entry.values()][0]
        self.chain_ids = [get_chain_id(entry) for entry in self.data['sequences']]
        self.num_chains = len(self.chain_ids)
        self.best_model = self._get_best_model()
        self.best_model_path = os.path.join(self.dir_path, self.best_model, 'model.cif')
        assert os.path.exists(self.best_model_path), f'AlphaFoldOutput.__init__: Best model path {self.best_model_path} does not exist.'



    def get_token_chain_ids(self):
        '''Get the token chain IDs from the full data output. This is a list of strings indicating which chain each position in a matrix (e.g. 
        the PAE matrix) belongs to.'''
        # Load the full data for each model; these SHOULD be the same across models, but load them all just to be safe.
        token_chain_ids = np.array([load_json(paths['full'])['token_chain_ids'] for paths in self.confidences.values()])
        assert np.all(np.expand_dims(token_chain_ids[0], axis=0) == token_chain_ids), 'AlphaFoldOutput.get_token_chain_ids: Token chain IDs are inconsistent.'
        return token_chain_ids[0]
    
    def _get_summary_data(self, field:str='iptm', mean_pool:bool=False, models:list=None):
        '''
        
        '''
        models = list(self.confidences.keys()) if (models is None) else models
        data = {model:load_json(paths['summary'])[field] for model, paths in self.confidences.items() if (model in models)}
        if mean_pool:
            data = np.mean(list(data.values()))
        return data

    def get_ptms(self, mean_pool:bool=False, models:list=None):
        return self._get_summary_data('ptm', mean_pool=mean_pool, models=models)

    def get_iptms(self, mean_pool:bool=False, models:list=None):
        return self._get_summary_data('iptm', mean_pool=mean_pool, models=models)


    def _get_full_data(self, field:str='contact_probs', mean_pool:bool=False, models:list=None):
        '''
        
        '''
        assert field in ['contact_probs', 'pae'], f'AlphaFoldOutput._get_full_data: Input field {field} is not recognized.'
        token_chain_ids = self.get_token_chain_ids()
        # Get the confidence data for the specified models (or all models if none are specified).
        data = {model:load_json(paths['full'])[field] for model, paths in self.confidences.items() if ((model is None) or (model in models))}

        if mean_pool:
            data = np.mean(list(data.values()), axis=0)
            return pd.DataFrame(data, index=token_chain_ids, columns=token_chain_ids)
        else: # Convert each individual matrix to DataFrames at the end if not mean pooling. 
            data = {model:pd.DataFrame(data_, index=token_chain_ids, columns=token_chain_ids) for model, data_ in data.items()}
        return data


    def get_contact_probs(self, mean_pool:bool=True, models:list=None):
        return self._get_full_data(field='contact_probs', mean_pool=mean_pool, models=models)   
    
    def get_paes(self, mean_pool:bool=True, models:list=None):
        return self._get_full_data(field='pae', mean_pool=mean_pool, models=models)


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

    def _get_best_model(self):
        '''Get the name of the model with the highest ranking score from the ranking_scores.csv file in the AlphaFold3 output root directory.
        
        :returns: The name of the model with the highest ranking score, which is of the form seed-x-sample-y. 
        '''
        models, scores = list(), list()
        for model, paths in self.confidences.items():
            score = load_json(paths['summary'])['ranking_score']
            scores.append(score)
            models.append(model)

        best_model = models[np.argmax(scores)]
        return best_model

    @staticmethod
    def _get_inputs(path:str):
        ''''''
        type_map = {'proteinChain':'protein', 'dnaSequence':'dna'}
        data = load_json(path)[0] # The entire thing is wrapped in a list, so grab the first index.

        df, i = list(), 0
        for group, info in enumerate(data['sequences']):
            type_ = list(info.keys())[0]
            count, seq = info[type_]['count'], info[type_].get('sequence', None)
            df += [{'type':type_, 'chain_id':CHAIN_IDS[j], 'group':group, 'seq':seq} for j in range(i, i + count)]
            i += count 

        df = pd.DataFrame(df)
        df['type'] = df['type'].apply(lambda type_: type_map.get(type_, type_))
        return df, data['dialect']


    
    def __init__(self, dir_path=None):
        # Paths to the confidence output for each model. 
        self.dir_path = dir_path 
        self.confidences = self._get_confidences()
        self.name = os.path.basename(dir_path)
        self.dir_path = dir_path

        self.inputs, dialect = AlphaFoldServerOutput._get_inputs(os.path.join(dir_path, f'{self.name}_job_request.json'))
        assert dialect == 'alphafoldserver', f'AlphaFoldServerOutput: Expected dialect alphafoldserver, but got {dialect}.'

        self.num_chains = len(self.inputs)
        self.chain_ids = self.inputs.chain_id.unique()
        self.best_model = self._get_best_model()


    def get_chain_ids(self, type_:str='protein', group:bool=False):

        inputs = self.inputs[self.inputs['type'] == type_].copy()
        if group:
            return [df.chain_id.unique().tolist() for _, df in inputs.groupby('group')]
        else:
            return inputs.chain_id.unique().tolist()
    

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
    

