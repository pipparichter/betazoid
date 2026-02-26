import numpy as np 
import os 
import pandas as pd 
import subprocess
from src.files import FASTAFile, BamFile, FLAGS
import io



BBMAP_PARAMS = dict()
BBMAP_PARAMS['minid'] = 0.95 # Minimum percent identity for an alignment to be considered mapped, applied only to the aligned portion of the read. 
BBMAP_PARAMS['idfilter'] = 0.97 # Filters the entire read for identity after alignment, applied over the entire read.
BBMAP_PARAMS['ambiguous'] = 'all'
BBMAP_PARAMS['minratio'] = 0.95 
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
    cmd += [f'in1={reads_path_2}']
    cmd += [f'ref={ref_path}']
    cmd += [f'out={output_path}']

    cmd += [f'{param}={str(value)}' for param, value in BBMAP_PARAMS.items()]
    cmd = ' '.join(cmd)
    print('run_bbmap:', cmd)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)


def recruit_reads(job_name, ref_path:str, n_iters:int=5, output_dir:str='.', reads_path_1:str=None, reads_path_2:str=None):
    output_paths = list()
    # For first iterations, don't care about non-primary alignments. 
    output_paths.append(os.path.join(output_dir, f'{job_name}.reads.{0}.bam'))
    run_bbmap(ref_path, reads_path_1=reads_path_1, reads_path_2=reads_path_2, output_path=output_paths[0])

    for i in range(1, n_iters):
        ref_path_i, n = BamFile.from_file(output_paths[-1]).to_fasta(include_flags=FLAGS['unmapped'])
        if (n == 0): # n is the number of sequences written to the FASTA file. 
            print(f'recruit_reads: No additional unmapped reads recruited. Exiting at iteration {i}.')
            break 
        print(f'recruit_reads: Recruited {n} additional reads at iteration {i}.')
        output_path_i = os.path.join(output_dir, f'{job_name}.reads.{i}.bam')
        run_bbmap(ref_path_i, output_path=output_path_i, reads_path_1=reads_path_1, reads_path_2=reads_path_2)
        output_paths.append(output_path_i)

    df = pd.concat([BamFile.from_file(path).to_df() for path in output_paths])
    df['strand'] = np.select([(df.reverse_strand & ~df.unmapped), (~df.reverse_strand & ~df.unmapped), (df.unmapped)], ['-', '+', 'none'], default=None)
    df['read_number'] = np.select([(df.read_paired & ~df.read_1), (df.read_paired & df.read_2), (~df.read_paired)], ['(1)', '(2)', ''], default='none')
    df['read_id'] = [read_id + read_number for read_id, read_number in zip(df.read_id, df.read_number)]
    df = df.drop_duplicates(['read_id', 'strand'])

    return df




class AssemblyGraph():
    pass 

# Rules for drawing an edge are



