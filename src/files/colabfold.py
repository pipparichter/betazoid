import pandas as pd 
import numpy as np 
import os 
import glob 
import json
import orjson
import re
import itertools 

# Note that ColabFold does not support multimers and/or ligands, which vastly simplifies parsing.

def load_json(path, errors:str='raise'):
    '''Read a JSON file from the input path using the faster orjson parser. '''
    try:
        with open(path, 'rb') as f:
            data = f.read()
            data = orjson.loads(data)
        return data
    except Exception as err:
        message = f'load_json: Could not decode {path}, {err}'
        if (errors == 'raise'):
            Exception(message)
        elif (errors == 'ignore'):
            print(message)
            return None



class ColabFoldOutput():

    patterns = dict()
    patterns['scores'] = r'{name}_scores_rank_(\d+)_alphafold2_ptm_model_(\d+)_seed_(\d+)\.json'
    patterns['model_relaxed'] = r'{name}_relaxed_rank_(\d+)_alphafold2_ptm_model_(\d+)_seed_(\d+)\.pdb'
    patterns['model_unrelaxed'] = r'{name}_relaxed_rank_(\d+)_alphafold2_ptm_model_(\d+)_seed_(\d+)\.pdb'
    # patterns['pae'] = r'{name}_predicted_aligned_error_v1.json'
    patterns['msa'] = r'{name}\.a3m'


    def _get_scores(self, paths:list):
        '''Map each model to the path to the JSON file with the corresponding scores. 
        '''
        scores = dict()
        for path in paths:
            for model in self.models:
                if model in path:
                    scores[model] = path 
    
        assert len(scores) > 0, f'ColabFoldOutput._get_scores: Was unable to find any confidence JSON files in {self.dir_path}.'
        assert len(scores) == self.num_models, f'ColabFoldOutput._get_scores: There is not one scores file per model. Found {len(scores)}, but expected {self.num_models}.'
        return scores 


    def __init__(self, path:str, use_relaxed:bool=False):
        '''Because ColabFold does not automatically organize the outputs into directories, this all needs to be done
        by pattern matching using the gene ID.
        '''
        # Get a subset of all paths in the directory to reduce search time. 
        self.path = path 
        self.dir_path = os.path.dirname(path)
        self.name = os.path.basename(path)
        self.patterns = {file_type:pattern.format(name=self.name) for file_type, pattern in ColabFoldOutput.patterns.items()}
        # print(f'ColabFoldOutput__init__: Loading structure data for {self.name} from {self.dir_path}')

        paths = {file_type:list() for file_type in self.patterns}
        # Organize the paths which contain the specified name into groups according to file type. 
        for path in glob.glob(self.path + '*'):
            for file_type, pattern in self.patterns.items():
                if re.search(pattern, path) is not None:
                    paths[file_type].append(os.path.abspath(path))

        self.has_relaxed = len(paths['model_relaxed']) > 0
        self.model_paths = paths['model_relaxed'] if (self.has_relaxed and use_relaxed) else paths['model_unrelaxed']
        self.models = [re.search(r'model_\d+_seed_\d+', path).group(0) for path in self.model_paths]
        self.model_ranks = [int(re.search(r'rank_(\d+)', path).group(1)) for path in self.model_paths]

        self.best_model = self.models[np.argmin(self.model_ranks)]
        self.best_model_path = self.model_paths[np.argmin(self.model_ranks)]

        self.num_models = len(self.models)
        self.scores = self._get_scores(paths['scores'])

        assert len(paths['msa']) == 1, f'ColabFoldOutput__init__: There should be exactly one MSA file.'
        self.msa_path = paths['msa'][0]

    def get_num_seeds(self):
        '''Get the number of seeds used for a single ColabFold prediction.'''
        seeds = [model.split('_')[-1] for model in self.models]
        seeds = set(seeds)
        for seed in seeds:
            assert re.match(r'\d+', seed) is not None, f'ColabFoldOutput: Expected all seeds to be made up of integers, but got {seed}.'
        return len(seeds)




    def get_chain_id_to_chain_map(self, **kwargs):
        '''
        :returns: A dictionary mapping the protein chain ID to the sequence. 
        '''
        return {'A':self._get_chain()}

    def get_chain_to_chain_id_map(self, **kwargs) -> dict:
        '''
        
        :returns: A dictionary mapping each protein sequence to the corresponding chain IDs.'''
        return {self._get_chain():['A']}


    def get_chain_ids(self, **kwargs):
        return ['A']
    
    def _get_score_data(self, field:str='pae', mean_pool:bool=False, models:list=None, best_model:bool=False):
        '''Obtain data from the scores file for each model. 

        :param field: The field to extract from the scores.json file. 
        :param mean_pool: If True, then mean pool the extracted data across all models or a subset of models. 
        
        '''
        if best_model:
            data = load_json(self.scores[self.best_model])[field]
        else:
            models = list(self.scores.keys()) if (models is None) else models
            data = {model:load_json(path)[field] for model, path in self.scores.items() if (model in models)}
            data = np.mean(list(data.values())) if mean_pool else data
        return data

    def get_ptms(self, mean_pool:bool=False, models:list=None, best_model:bool=False):
        return self._get_score_data('ptm', mean_pool=mean_pool, models=models, best_model=best_model)

    def get_iptms(self, mean_pool:bool=False, models:list=None, best_model:bool=False):
        if best_model or mean_pool:
            return None
        models = list(self.scores.keys()) if (models is None) else models
        return {model:None for model in models}

    def get_paes(self, mean_pool:bool=False, models:list=None, best_model:bool=False):
        return self._get_score_data('pae', mean_pool=mean_pool, models=models, best_model=best_model)
    
    def get_plddts(self, mean_pool:bool=False, models:list=None, best_model:bool=False):
        return self._get_score_data('plddt', mean_pool=mean_pool, models=models, best_model=best_model)
    
    def get_msa(self) -> str:
        with open(self.msa_path, 'r') as f:
            msa = f.read()
        return msa

    def get_msas(self) -> list:
        '''Load the MSAs used to generate the structures for all models. ColabFold only supports monomers, so only 
        unpaired MSAs are returned in the output.
        
        :returns: A list containing a single dictionary, with one entry containing the unpaired MSA. There will always only 
            be one MSA for the ColabFold structures, but this approach allows interface consistency with the AlphaFoldOutput class.
        '''
        msa = self.get_msa()
        return [{'unpaired':msa}]

    def _get_chain(self):
        '''The only way to get the sequence out of the ColabFold ouput is to read it from the MSA.'''
        msa = self.get_msa()
        seqs = re.split(r'^>.*$', msa, flags=re.MULTILINE)
        return seqs[0].replace('\n', '')

    def get_msa_metadata(self):
        '''Read and parse the headers of a ColabFold MSA file. These MSAs were constructed using MMseqs, and each header line contains 
        information about the alignment quality. 
        
        :param path: The path to the MSA file. 
        :returns: A dictionary containing (1) the full MSA text, (2) the number of sequences in the MSA, (3) the mean bit score, (4) the minimum bit score, 
            and (5) the maximum bit score. 
        '''
        msa = self.get_msa()

        # It seems as though ColabFold embeds statistics from the MMseqs search in the FASTA headers. 
        fields = ['id', 'bit_score', 'identity', 'e_value', 'query_start', 'query_end', 'query_length', 'target_start', 'target_end', 'target_length']
        headers = [line for line in msa.split('\n') if line.startswith('>')]

        df = pd.DataFrame([dict(zip(fields, header.split())) for header in headers])
        if len(df.columns) == 1:
            return {'msa':msa, 'msa_num_seqs':0}, None
        df = df[~df.bit_score.isnull()].copy()

        metadata = dict()
        metadata['msa_num_seqs'] = len(df)
        for metric in ['bit_score', 'identity', 'e_value']:            
            metadata[f'msa_mean_{metric}'] = df[metric].astype(float).mean()
            metadata[f'msa_min_{metric}'] = df[metric].astype(float).min()
            metadata[f'msa_max_{metric}'] = df[metric].astype(float).max()
        return metadata, df
