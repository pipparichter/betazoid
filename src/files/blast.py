import re 
import json 
import pandas as pd 
import numpy as np 
import os 



class BLASTFile():

    field_map = dict()
    field_map['accession'] = 'subject_id'
    field_map['query_title'] = 'id'
    field_map['title'] = 'subject_description'
    field_map['sciname'] = 'subject_taxon'
    field_map['taxid'] = 'subject_taxonomy_id'
    field_map['bit_score'] = 'bit_score'
    field_map['evalue'] = 'e_value'
    field_map['identity'] = 'identity'
    field_map['positive'] = 'positive'
    field_map['hit_from'] = 'subject_alignment_start'
    field_map['hit_to'] = 'subject_alignment_stop'
    field_map['query_from'] = 'query_alignment_start'
    field_map['query_to'] = 'query_alignment_stop'
    field_map['gaps'] = 'n_gaps'
    field_map['align_len'] = 'alignment_length'
    field_map['qseq'] = 'query_seq'
    field_map['hseq'] = 'subject_seq'
    field_map['len'] = 'subject_length'
    field_map['query_len'] = 'query_length'
    field_map['midline'] = 'alignment'

    def __init__(self):
        pass 

    def _parse_json_file(path:str) -> pd.DataFrame:

        with open(path, 'r') as f:
            content = json.load(f)

        df = []
        for query in content['BlastOutput2']:
            results = query['report']['results']['search']
            query_info = {field:value for field, value in results.items() if (field != 'hits')}

            if (len(results['hits']) == 0):
                continue
            
            # Probably fine to just get the first high-scoring pair. There is generally just one anyway. 
            for hit in results['hits']:
                row = hit['description'][0].copy()
                row['len'] = hit['len']
                row.update(query_info)
                row.update(hit['hsps'][0]) # Get the best high-scoring pair. 
                df.append(row)

        return pd.DataFrame(df)
        # df = df[list(BLASTFileJSON.field_map.keys())].rename(columns=BLASTFileJSON.field_map) # Rename columns according to the field map. 

        

    @staticmethod
    def _from_ggkbase_text_file(path:str) -> pd.DataFrame:
        '''For parsing BLAST output files from ggKbase.
        Example BLAST target entry looks like this. The first number after the > is the feature ID, which can be located :

        >128941808 Mayberry_contig_91_7252045_length_5790_multi_2_in_0_out_0_2
        Length=750

        Score = 285 bits (315),  Expect = 7e-72
        Identities = 513/742 (69%), Gaps = 16/742 (2%)
        Strand=Plus/Plus

        The query entries look like:

        Query= Final_SR-
        VP_05_06_2024_coassembly_19kb_linear_ECE_26_1334_complete_38

        Length=297
        '''

        patterns = [r'>(?P<feature_number>\d+) (?P<feature_id>[^\n]+)']
        patterns += [r'Length=(?P<subject_length>\d+)']
        patterns += [r'Score = (?P<bit_score>[\d\.]+) bits']
        patterns += [r'Expect = (?P<e_value>[e\-\.\d]+)']
        patterns += [r'Identities = (?P<n_identical>\d+)/(?P<alignment_length>\d+)']
        patterns += [r'Gaps = (?P<n_gaps>\d+)/(?P<alignment_length>\d+)']

        with open(path, 'r') as f:
            content = f.read()
            
        # print(f"BLASTFile._from_ggkbase_text_file: Found {content.count('Query=')} query entries in the file.")

        query_info_pattern = r'Query=\s*(?P<query_id>.+?)\n\nLength=(?P<query_length>\d+)' # Make sure the first match is non-greedy.

        queries = list(re.finditer(query_info_pattern, content, flags=re.DOTALL|re.MULTILINE))
        print(f'BLASTFile._from_ggkbase_text_file: Found {len(queries)} query entries in the file.')

        df = list()
        for i, query in enumerate(queries):

            start, end = query.end(), queries[i + 1].start() if ((i + 1) < len(queries)) else -1
            content_ = content[start:end]
            content_ = re.sub(r'Lambda.*', '', content_, flags=re.DOTALL) # remove the extra data at the end of the chunk. 

            query_info = query.groupdict()
            query_info['query_id'] = re.sub(r'\(.*\)', '', query_info['query_id'])  # Strip any (id=...) stuff tagged on to the query ID. 
            query_info['query_id'] = query_info['query_id'].replace('\n', '').replace(' ', '')         
  

            for target in re.findall(r'(?=>)(.*?)(?=^>)', content_, flags=re.MULTILINE|re.DOTALL):  # Make sure to use .+? so the match is lazy, not greedy. 
                if len(target) == 0: # Gets empty strings for some reason.
                    continue

                row = query_info.copy()
                for pattern in patterns:
                    match_ = re.search(pattern, target, flags=re.DOTALL)
                    if match_ is not None:
                        row.update(match_.groupdict())
                df.append(row)


        df = pd.DataFrame(df)
        df['percent_identity'] = df.n_identical.astype(int) / df.alignment_length.astype(int)

        cols = [col for col in df.columns if ('length' in col)]
        for col in ['bit_score', 'n_gaps', 'n_identical', 'e_value'] + cols:
            df[col] = df[col].apply(pd.to_numeric)

        return df 


    @classmethod
    def from_file(cls, path:str, format:str=None):

        format = os.path.splitext(os.path.basename(path))[-1] if (format is None) else format
        format = format.replace('.', '').lower() # Remove any dots or capital letters in the file extension. 
        assert format in ['txt', 'json'], 'BLASTFile.from_file: File format {format} is not supported. Must be one of ' + str(['txt', 'json'])

        parser = BLASTFile._from_ggkbase_text_file if (format == 'txt') else BLASTFile._from_json_file
        df = parser(path)

        obj = cls()
        obj.df = df 
        obj.path = os.path.abspath(path)
        return obj

    def to_df(self):
        return self.df.copy()

