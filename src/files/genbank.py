import os 
import re
import pandas as pd 
from src.files.fasta import FASTAFile

feature_pattern = r'^\s{0,5}(\S+)\s+(\S+)'
seq_pattern = r'/translation="([^"]+)"'
product_pattern = r'/product="([^"]+)"'
locus_tag_pattern = r'/locus_tag="([^"]+)"'
# locus_tag_pattern = r'/locus_tag="([^"]+)"'


# def genbank_read_features(path):
#     content= ''
#     with open(path, 'r') as f:
#         read = False
#         for line in f.readlines():
#             if 'FEATURES' in line:
#                 read = True 
#                 continue 
#             if 'ORIGIN' in line:
#                 read = False 
#                 break 
#             if read:
#                 content += line
#     return content

def _clean_product(product:str):
    product = re.sub(r'[a-z]{3}:[^\s]+' , '', product)
    product = re.sub(r'\(.*\)' , '', product)
    product = re.sub(r'\b\w*_\w*\b', '', product) # The \b character represents a word boundary.
    # product = re.sub(r'\b\w*=\w*\b', '', product) # The \b character represents a word boundary.
    product = re.sub(r'\s*\S*=\S*.*$', '', product) # Remove any word with an equal sign and all words after it. 
    if 'Tax' in product:
        product = 'hypothetical protein'
    product = product.strip()
    if len(product) < 2:
        return 'none'
    # product = product[0].lower() + product[1:]
    return product.strip()



def genbank_read_features(path):
    with open(path, 'r') as f:
        content = f.read()
    # Contig label should be the first non-whitespace thing after the split. 
    contig_ids = [re.search(r'([^\s]+)', contig, flags=re.MULTILINE).group(0) for contig in content.split('LOCUS')[1:]]
    features = dict()

    for contig_id, contig in zip(contig_ids, content.split('LOCUS')[1:]):
        read = False
        s = ''
        for line in contig.split('\n'):
            if 'FEATURES' in line:
                read = True 
                continue 
            if 'ORIGIN' in line:
                read = False 
                break 
            if read:
                s += line + '\n'
        features[contig_id] = s
    return features

def genbank_parse_features(features:str):
    # feature_pattern = r'(?=^\s{0,5}(\S+)\s+(\S+))'
    features = re.split(feature_pattern, features, flags=re.MULTILINE)
    features = [feature for feature in features if len(feature) > 0]
    assert (len(features) % 3) == 0, 'GenBankFile: Expected the number of features to be divisible by three.'
    features = [(features[i], features[i + 1], features[i + 2]) for i in range(0, len(features), 3)]
    return features 
    

def genbank_parse_feature(feature):
    # feature = re.sub(r'[\s\n]{2,}', '', feature) # Remove any whitespace bigger than one. 
    info = dict()
    info['seq'] = re.search(seq_pattern, feature, flags=re.MULTILINE).group(1).replace('\n', '').replace(' ', '').replace('*', '')
    info['product'] = re.search(product_pattern, feature, flags=re.MULTILINE).group(1).replace('\n', '')
    info['product'] = re.sub(r'\s{2,}', ' ', info['product']) # Remove all the extra whitespace.
    # info['locus_tag'] = re.search(locus_tag_pattern, feature, flags=re.MULTILINE).group(1).replace('\n', '')
    return info 


class GenBankFile():

    aas = 'ACDEFGHIKLMNPQRSTVWYX*'
    seq_pattern = r'/transl="([^"]+)"'

    def __int__(self):
        pass 

    @classmethod
    def from_file(cls, path:str):

        df = list()

        contig_index = 1
        for contig_id, features in genbank_read_features(path).items():
            contig_df = list()
            for feature_type, coordinate, feature in genbank_parse_features(features):
                if feature_type != 'CDS':
                    continue 
                row = {'feature_type':feature_type, 'coordinate':coordinate}
                row.update(genbank_parse_feature(feature))
                contig_df.append(row)
            contig_df = pd.DataFrame(contig_df)
            contig_df['contig_index'] = contig_index
            contig_df['contig_id'] = contig_id
            contig_df['gene_id'] = [f'{contig_index}_{i + 1}' for i in range(len(contig_df))]
            contig_df = contig_df.rename(columns={'translation':'seq'})
            # contig_df['locus_tag'] = [f'{contig_id}_{i + 1}' for i in range(len(contig_df))]
            df.append(contig_df)
            contig_index += 1 

        df = pd.concat(df)
        
        obj = cls()
        obj.df = df.copy()
        obj.path = path 
        obj.file_name = os.path.basename(path)

        return obj 
    
    def to_fasta(self, path:str, mode:str='w'):
        fasta_file = FASTAFile.from_df(self.df.set_index('gene_id'))
        fasta_file.write(path, mode=mode)
    
    def to_gff(self, path:str, model='Prodigal', contig_id:str=None, mode='w'):
        if mode == 'w':
            with open(path, 'w') as f:
                f.write('##gff-version  3\n')
        
        cols = ['contig_id', 'model', 'feature_type', 'start', 'stop', 'score', 'strand', 'frame', 'description']
        gff_df = self.to_df()
        gff_df['model'] = model
        gff_df['feature_type'] = 'CDS'
        gff_df['start'] = gff_df.start + 1
        gff_df['frame'] = 0
        gff_df['score'] = '-'
        if contig_id is not None:
            gff_df['contig_id'] = gff_df.contig_id.fillna()
        gff_df['description'] = [f'ID={gene_id}' for gene_id in self.df.gene_id]
        gff_df = gff_df[cols].copy()
        gff_df.to_csv(path, header=None, sep='\t', index=False, mode='a') # Append to the file which already has the header. 


    def to_df(self, clean_product:bool=False):
        df = self.df.copy()
        df = df.set_index('gene_id')
        df['strand'] = ['-' if ('comp' in coordinate) else '+' for coordinate in  df.coordinate]
        df['stop'] = [int(re.search(r'\.\.[<>]*(\d+)', coordinate).group(1)) for coordinate in df.coordinate]
        df['start'] = [int(re.search(r'(\d+)[<>]*\.\.', coordinate).group(1)) for coordinate in df.coordinate]
        if clean_product:
            df['product'] = [_clean_product(product) for product in df['product']]

        return df

    # def to_fasta(path:str):

        