import numpy as np 
import os 
import pandas as pd 
import subprocess
from src.files import FASTAFile, BamFile, FLAGS
import io
from Bio.Seq import Seq

get_reverse_complement = lambda seq : str(Seq(seq).reverse_complement())

BBMAP_PARAMS = dict()
BBMAP_PARAMS['minid'] = 0.95 # Minimum percent identity for an alignment to be considered mapped, applied only to the aligned portion of the read. 
BBMAP_PARAMS['idfilter'] = 0.97 # Filters the entire read for identity after alignment, applied over the entire read.
BBMAP_PARAMS['ambiguous'] = 'all'
BBMAP_PARAMS['minratio'] = 0.5
BBMAP_PARAMS['editfilter'] = 5 # Consider reads with fewer than 5 indels or mismatches as mapped.
BBMAP_PARAMS['local'] = 't' # Allow local alignments. 
# BBMAP_PARAMS['unmapped'] = 't' # Include reads that fail the thresholds. 
BBMAP_PARAMS['pigz'] = 't'
BBMAP_PARAMS['unpigz'] = 't'
BBMAP_PARAMS['threads'] = 64
BBMAP_PARAMS['k'] = 13 # K-mer size (13 is default for short-read data)
BBMAP_PARAMS['minhits'] = 2 # Minimum number of K-mer matches to consider a read (2 is default)

# bbmap.sh pigz=t unpigz=t ambiguous=random minid=0.96 idfilter=0.97 threads=64 out=stdout.sam editfilter=5 in1=/groups/banfield/sequences/2014/16ft_4/raw.d/16ft_4_CZBZ.6237.3.40316_trim_clean.PE.1.fastq.gz in2=/groups/banfield/sequences/2014/16ft_4/raw.d/16ft_4_CZBZ.6237.3.40316_trim_clean.PE.2.fastq.gz ref=/home/philippar/scaffold_2135.fasta nodisk | shrinksam | sambam > /home/philippar/scaffold_2135.bam
def run_bbmap(ref_path, reads_path_1:str=None, reads_path_2:str=None, output_path:str=None):
    ''''''
    cmd = ['bbmap.sh']
    cmd += [f'in1={reads_path_1}']
    cmd += [f'in2={reads_path_2}']
    cmd += [f'ref={ref_path}']
    cmd += [f'out={output_path}']

    cmd += [f'{param}={str(value)}' for param, value in BBMAP_PARAMS.items()]
    cmd = ' '.join(cmd)
    print('run_bbmap:', cmd)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)


def check_reads(df):
    '''Confirm that the recruited reads match what I am assuming about them.'''
    checks = dict()
    checks['one_mate_mapped'] = np.all((~df.mate_unmapped) | (~df.unmapped))
    checks['read_strand_assigned'] = np.any(df.strand.isnull())
    for name, value in checks.items():
        assert value, f'check_reads: Failed check {name}.'


def recruit_reads(job_name, ref_path:str, n_iters:int=5, output_dir:str='.', reads_path_1:str=None, reads_path_2:str=None):
    output_paths = list()
    # For first iterations, don't care about non-primary alignments. 
    output_paths.append(os.path.join(output_dir, f'{job_name}.reads.{0}.bam'))
    run_bbmap(ref_path, reads_path_1=reads_path_1, reads_path_2=reads_path_2, output_path=output_paths[0])

    for i in range(1, n_iters):
        ref_path_i, n = BamFile.from_file(output_paths[-1]).to_fasta(include_flags=FLAGS['unmapped'] + FLAGS['read_paired'], exclude_flags=FLAGS['mate_unmapped'])
        if (n == 0): # n is the number of sequences written to the FASTA file. 
            print(f'recruit_reads: No additional unmapped reads recruited. Exiting at iteration {i}.')
            break 
        print(f'recruit_reads: Recruiting using {n} unmapped reads at iteration {i}.')
        output_path_i = os.path.join(output_dir, f'{job_name}.reads.{i}.bam')
        run_bbmap(ref_path_i, output_path=output_path_i, reads_path_1=reads_path_1, reads_path_2=reads_path_2)
        output_paths.append(output_path_i)

    df = pd.concat([BamFile.from_file(path).to_df() for path in output_paths])
    df['read_number'] = np.select([(df.read_paired & ~df.read_1), (df.read_paired & df.read_2), (~df.read_paired)], ['1', '2', ''], default='')
    df = df.drop_duplicates(['read_id', 'strand', 'read_number'])
    return df



# Need to think a bit about the best way to do this... 
# Two approaches are k-mer or overlap graph. Apparently k-mer graphs are better for short reads, but overlap graphs more closely mimic the manual
# curation process and might allow easier integration of pair information. 

# https://academic.oup.com/bioinformatics/article/21/suppl_2/ii79/227189?login=true
# https://www.pnas.org/doi/10.1073/pnas.171285098

TMP_DIR = '/home/philippar/tmp'

def get_reads(df):
    df = df[~df.strand.isnull()].copy() # Remove all cases where the strand is unknown, which occurs if both reads are mapped in the same direction. 
    # Include reads in the orientation they were mapped to the contig. If the read's mate is mapped in one orientation, even if the read
    # itself is unmapped, enforce the opposite orientation
    ids, seqs = list(), list()
    for row in df.itertuples():
        ids += [f'{row.read_id} {row.read_number} {row.strand}']
        seqs += [get_reverse_complement(row.seq) if (row.strand == '-') else row.seq]
    return ids, seqs

MMSEQS_FIELDS = ['query', 'target', 'alnlen', 'qcov', 'qstart', 'qend', 'tstart', 'tend', 'fident', 'qseq', 'tseq', 'qaln', 'taln']

def align_reads(reads):
    tmp_dir = os.path.join(TMP_DIR, 'tmp')
    fasta_file = FASTAFile()
    fasta_file.ids, fasta_file.seqs = reads
    fasta_file.write(os.path.join(TMP_DIR, 'reads.fasta'))

    database_path = os.path.join(TMP_DIR, 'readsDB')
    prefilter_database_path = os.path.join(TMP_DIR, 'readsDB_prefilter')
    aligned_database_path = os.path.join(TMP_DIR, 'readsDB_aligned')
    alignment_path = os.path.join(TMP_DIR, 'alignments.tsv')

    prefilter_options = ['--k', '7', '--kmers-per-seq', '300', '--min-ungapped-score', '15']
    fields = ','.join(MMSEQS_FIELDS)
    align_options = ['-a', '--mask', 'no', '--alignment-mode', '2', '--max-seqs', '10', '--min-seq-id', '95', '--format-output', f'"{fields}"']

    subprocess.run(['mmseqs', 'createdb', os.path.join(TMP_DIR, 'reads.fasta'), database_path, '--dbtype', '2'], check=True)
    subprocess.run(['mmseqs', 'prefilter', database_path, database_path, prefilter_database_path, tmp_dir] + prefilter_options, check=True)
    subprocess.run(['mmseqs', 'align', database_path, database_path, prefilter_database_path, aligned_database_path, tmp_dir] + align_options, check=True)
    subprocess.run(['mmseqs', 'convertalis', database_path, database_path, aligned_database_path, tmp_dir, alignment_path] + align_options, check=True)




# I think the best way to do this is just going to be an mmseqs pairwise alignment between all recruited reads, and then reading in the file. 

class StringGraph():

    def __init__(self, df:pd.DataFrame, overlap_length:int=10):

        reads = get_reads(df)


# Representing every read of length n as a collection of n - k + 1 overlapping k-tuples (continuous short strings of fixed length k)
class KmerGraph():
    def _get_sequences(reads_df):

        # Include reads in the orientation they were mapped to the contig. If the read's mate is mapped in one orientation, even if the read
        # itself is unmapped, enforce the opposite orientation
        pass 

    def __init__(self, reads_df:pd.DataFrame, k:int=10):
        self.k = k 
        

# Rules for drawing an edge ar



