import os 
import re
import pandas as pd 
import numpy as np 


class GenBankFile():

    def __int__(self):
        pass 
    
    @staticmethod
    def _parse_metadata(metadata:str):
        metadata = re.split(r'^(?=\w)', metadata.strip(), flags=re.MULTILINE)
        metadata = [re.sub(r'\n|\s{2,}', ' ', entry.strip()) for entry in metadata]

        info = dict()
        for entry in metadata: # Split metadata into lines that don't start with an indent. 
            if entry.startswith('KEYWORDS'):
                info['keywords'] = entry.replace('KEYWORDS ', '')
            if entry.startswith('DBSOURCE'):
                info['dbsource'] = entry.replace('DBSOURCE ', '')
            if entry.startswith('SOURCE'):
                entry = entry.replace('SOURCE', '').strip()
                if 'ORGANISM' in entry:
                    info['source'], info['organism'] = entry.split(' ORGANISM ')
                else:
                    info['source'] = entry
        return info 

    @staticmethod
    def _parse_features(features:str):
        df = list()

        features = re.split(r'^\s{5}(?=\w)', features.strip(), flags=re.MULTILINE) # Split into individual indented blocks. 
        for feature in features:
            info = dict()

            feature = feature.strip()
            feature = [re.sub(r'\n|\s{2,}', ' ', entry.strip()) for entry in re.split(r'\/', feature)]

            try:
                # First line is the feature type and coordinates. 
                info['feature_type'], info['coordinate'] = feature[0].split()
            except:
                # print(f'GenBankFile._parse_features: Failed to obtain coordinates from feature {feature[0]}')
                continue 

            for entry in feature[1:]:
                match = re.match(r'([\w\W_]+)="(.+)"', entry)
                if match is None:
                    continue
                info[match.group(1)] = match.group(2)
            df.append(info)

        return pd.DataFrame(df)


#      Site            order(4726,4730,4780,4782..4783,4786..4787,4789,4798..4806,4814,4816..4817,4819..4820,4823..4824,4827)
                    #  /site_type="other"
                    #  /note="putative actin binding site [polypeptide binding]"
                    #  /db_xref="CDD:409049"

    @classmethod
    def from_file(cls, path:str, load_origins_only:bool=False):

        with open(path, 'r') as f:
            content = f.read()

        # re.split returns everything before the first match as the first element, which we want to discard with [1:] slicing. This should be an empty string.x
        loci = re.split(r'(?=^LOCUS\b)', content, flags=re.MULTILINE)[1:]
        locus_ids = [locus.strip().split()[1] for locus in loci]

        obj = cls()
        origins = dict()
        df = list()

        for locus, locus_id in zip(loci, locus_ids):
            metadata, features, origin = re.split(r'^ORIGIN|FEATURES.+', locus, flags=re.MULTILINE)
            origins[locus_id] = re.sub(r'[^a-zA-Z]', '', origin)
            if not load_origins_only:
                metadata = GenBankFile._parse_metadata(metadata)
                metadata.update({'locus_id':locus_id})
                df.append(GenBankFile._parse_features(features).assign(**metadata))
        df = pd.concat(df) if (len(df) > 0) else None

        obj.df = df
        obj.origins = origins
        obj.locus_ids = locus_ids 
        return obj


    def to_df(self):
        df = self.df.copy()
        df = df[df.feature_type != 'source'].copy() # Source is not really important.
        for field in df.columns:
            if np.all(df[field].isnull()):
                df = df.drop(columns=[field]) 
        # df['strand'] = ['-' if ('comp' in coordinate) else '+' for coordinate in  df.coordinate]
        # df['stop'] = [int(re.search(r'\.\.[<>]*(\d+)', coordinate).group(1)) for coordinate in df.coordinate]
        # df['start'] = [int(re.search(r'(\d+)[<>]*\.\.', coordinate).group(1)) for coordinate in df.coordinate]
        return df


    
    # def to_fasta(self, path:str, mode:str='w'):
    #     fasta_file = FASTAFile.from_df(self.df.set_index('gene_id'))
    #     fasta_file.write(path, mode=mode)
    



    # def to_fasta(path:str):


    # def to_gff(self, path:str, model='Prodigal', contig_id:str=None, mode='w'):
    #     if mode == 'w':
    #         with open(path, 'w') as f:
    #             f.write('##gff-version  3\n')
        
    #     cols = ['contig_id', 'model', 'feature_type', 'start', 'stop', 'score', 'strand', 'frame', 'description']
    #     gff_df = self.to_df()
    #     gff_df['model'] = model
    #     gff_df['feature_type'] = 'CDS'
    #     gff_df['start'] = gff_df.start + 1
    #     gff_df['frame'] = 0
    #     gff_df['score'] = '-'
    #     if contig_id is not None:
    #         gff_df['contig_id'] = gff_df.contig_id.fillna()
    #     gff_df['description'] = [f'ID={gene_id}' for gene_id in self.df.gene_id]
    #     gff_df = gff_df[cols].copy()
    #     gff_df.to_csv(path, header=None, sep='\t', index=False, mode='a') # Append to the file which already has the header. 


        