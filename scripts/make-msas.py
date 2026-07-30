MMSEQS_DIR = '../data/genes/mmseqs'
MMSEQS_DATABASE_DIR = '../data/genes/mmseqs/db'
MMSEQS_TMP_DIR = '../data/tmp'

@dataclass
class Databases:
    input: str
    output_search: str
    output_align: str
    output_msa: str

name = 'level_1_conserved_clusters'
databases = Databases(f'{MMSEQS_DATABASE_DIR}/{name}', f'{MMSEQS_DATABASE_DIR}/{name}.out', f'{MMSEQS_DATABASE_DIR}/{name}.aln', f'{MMSEQS_DATABASE_DIR}/{name}.msa')


def make_msas(path:str, databases=databases, output_dir_path:str=MMSEQS_DIR, sensitivity=9.5):
    '''Use MMseqs align utilities to construct a3m-format alignment files for each sequence in the input FASTA file. The pipeline
    is as follows:
        (1) Use the input FASTA file to construct an MMseqs database. 
        (2) Generate a search result database by querying the input database against itself. 
        (3) Construct full alignments using the results of the search. 
        (4) Convert the alignments to MSAs. 
        (5) Unpack the alignment output info a3m files compatible with AlphaFold; there will be one file per query.
        
    :param path: The path to the FASTA file containing all sequences to create MSAs for. 
    :param databases: The Databases dataclass containing the paths and names for all requisite intermediate MMseqs databases.    
    :param output_dir_path: The directory where all final a3m files will be deposited. 
    :param sensitivity: The search sensitivity for the initial MMseqs search step. 
    '''
    subprocess.run(f'mmseqs createdb {path} {databases.input}', shell=True, check=True)
    subprocess.run(f'mmseqs search {databases.input} {databases.input} {databases.output_search} {MMSEQS_TMP_DIR} -s {sensitivity}', shell=True, check=True)
    subprocess.run(f'mmseqs align {databases.input} {databases.input} {databases.output_search} {databases.output_align}', shell=True, check=True)
    subprocess.run(f'mmseqs result2msa {databases.input} {databases.input} {databases.output_align} {databases.output_msa}', shell=True, check=True)
    subprocess.run(f'mmseqs unpackdb {databases.output_msa} {output_dir_path} --unpack-suffix .a3m', shell=True, check=True)

    # The unpacked database files just have numerical IDs, so want to rename according to their actual gene IDs. 
    # This is done by mapping the numerical ID to the original gene ID using the lookup table associated with the input database. 
    num_to_gene_id_map = pd.read_csv(databases.input + '.lookup', sep=r'\s+', names=['num', 'gene_id'], index_col=0, usecols=[0, 1]).gene_id.to_dict()
    for file_name in os.listdir(output_dir_path):
        if not re.match(r'\d+\.a3m', file_name):
            continue 
        n = int(re.match(r'(\d+)\.a3m', file_name).group(1))
        old_path = os.path.join(output_dir_path, file_name)
        new_path = os.path.join(output_dir_path, f'{num_to_gene_id_map[n]}.a3m')
        subprocess.run(f'mv {old_path} {new_path}', shell=True, check=True)