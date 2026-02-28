import numpy as np 
import os 
import pandas as pd 
import subprocess
from src.files import FASTAFile, BamFile, FLAGS
import io
from Bio.Seq import Seq

# https://academic.oup.com/bioinformatics/article/21/suppl_2/ii79/227189?login=true
# https://www.pnas.org/doi/10.1073/pnas.171285098


get_reverse_complement = lambda seq : str(Seq(seq).reverse_complement())

BBMAP_PARAMS = dict()
BBMAP_PARAMS['minid'] = 0.95 # Minimum percent identity for an alignment to be considered mapped, applied only to the aligned portion of the read. 
BBMAP_PARAMS['idfilter'] = 0.97 # Filters the entire read for identity after alignment, applied over the entire read.
BBMAP_PARAMS['ambiguous'] = 'all' # Ambiguous means all alignments are equally good. When writing ambiguous alignments, each alignment becomes a separate SAM record.
# BBMAP_PARAMS['mappedonly'] = 't' # If true, treats out like outm.
BBMAP_PARAMS['minratio'] = 0.5
BBMAP_PARAMS['editfilter'] = 5 # Consider reads with fewer than 5 indels or mismatches as mapped.
BBMAP_PARAMS['local'] = 'f' # Allow local alignments. Actually, don't.
# BBMAP_PARAMS['unmapped'] = 't' # Include reads that fail the thresholds. 
BBMAP_PARAMS['pigz'] = 't'
BBMAP_PARAMS['unpigz'] = 't'
BBMAP_PARAMS['threads'] = 64
BBMAP_PARAMS['k'] = 13 # K-mer size (13 is default for short-read data)
BBMAP_PARAMS['minhits'] = 1 # Minimum number of K-mer matches to consider a read (1 is default)
BBMAP_PARAMS['nodisk'] = 't' # No temporary files written.
BBMAP_PARAMS['notags'] = 't' # Turn off optional tags. 
# outm: Write only mapped reads to this file. Includes unmapped paired reads with a mapped mate.

def _run_bbmap(ref_path, reads_path_1:str=None, reads_path_2:str=None, output_path:str=None):
    ''''''
    if os.path.exists(output_path):
        print(f'_run_bbmap: Using existing output stored at {output_path}')
        return 
    
    cmd = ['bbmap.sh']
    cmd += [f'in1={reads_path_1}']
    cmd += [f'in2={reads_path_2}']
    cmd += [f'ref={ref_path}']
    cmd += [f'out=stdout.sam']
    # cmd += [f'out={output_path}']

    cmd += [f'{param}={str(value)}' for param, value in BBMAP_PARAMS.items()]
    cmd = ' '.join(cmd)

    print('_run_bbmap:', cmd)
    # Getting different results if I don't do this. 
    cmd += f'| shrinksam | sambam > {output_path}'

    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)


# A handful of possible scenarios... 
# (1) Both mates map in opposite orientations: Include both reads in the specified orientation. 
# (2) One mate maps in one orientation: Include the mapped read and the unmapped read in the opposite orientation. 
# (3) One mate maps in both orientations: Include both members of the pair in both orientations. 
# (4) Both mates map in the same orientation: Include both members of the pair in both orientations. 
# If the reads were recruited in iterations after the first iteration, their orientation relative to the reference is no longer important. 
def _write_reads_to_fasta(df, output_path:str=None):
    '''Takes the DataFrame read from a BAM file as input. Outputs a list of sequences to use for building the graph, with IDS encoding the direction and read ID.'''
    # Include reads in the orientation they were mapped to the contig. If the read's mate is mapped in one orientation, even if the read
    # itself is unmapped, enforce the opposite orientation
    id_map = {read_id:i for i, read_id in enumerate(np.sort(df.read_id.unique()))} # Map paired read IDS to integers to make life easier.
    ids, seqs = list(), list()
    for row in df.itertuples():
        
        n = row.read_number
        i = id_map[row.read_id]
        assert row.orientation != 'XX', 'get_reads: There should not be any reads where both members of the pair are unmapped.'

        if (row.iteration > 0): # We don't know anything about the true orientation relative to the seed contig for reads recruited after the first iteration. 
            ids += [f'{i}.{n}.F', f'{i}.{n}.R']
            seqs += [row.seq, get_reverse_complement(row.seq)]

        elif (row.orientation == 'RF') or (row.orientation == 'FR'):
            ids += [f'{i}.{n}.{row.orientation[0]}']
            seqs += [get_reverse_complement(row.seq) if (row.orientation[0] == 'R') else row.seq]
        elif (row.orientation == 'RR') or (row.orientation == 'FF'):
            ids += [f'{i}.{n}.F', f'{i}.{n}.R']
            seqs += [row.seq, get_reverse_complement(row.seq)]
        elif (row.orientation == 'XF') or (row.orientation == 'XR'):
            orientation = 'F' if (row.orientation[-1] == 'R') else 'F'
            ids += [f'{i}.{n}.{orientation}']
            seqs += [row.seq, get_reverse_complement(row.seq)]
        else:
            continue
    print(f'get_reads: Writing {len(ids)} sequences to {output_path} from {len(id_map)} pairs.')
    fasta_file = FASTAFile()
    fasta_file.ids, fasta_file.seqs, fasta_file.descriptions = ids, seqs, []
    fasta_file.write(output_path)


def recruit(ref_path:str, n_iters:int=5, output_dir:str='.', reads_path_1:str=None, reads_path_2:str=None):
    output_paths = list()
    # For first iterations, don't care about non-primary alignments. 
    output_paths.append(os.path.join(output_dir, f'reads.{0}.bam'))
    _run_bbmap(ref_path, reads_path_1=reads_path_1, reads_path_2=reads_path_2, output_path=output_paths[0])

    for i in range(1, n_iters):
        ref_path_i, n = BamFile.from_file(output_paths[-1]).to_fasta(include_flags=FLAGS['unmapped'] + FLAGS['read_paired'], exclude_flags=FLAGS['mate_unmapped'])
        if (n == 0): # n is the number of sequences written to the FASTA file. 
            print(f'recruit: No additional unmapped reads recruited. Exiting at iteration {i}.')
            break 
        print(f'recruit: Recruiting using {n} unmapped reads at iteration {i}.')
        output_path_i = os.path.join(output_dir, f'reads.{i}.bam')
        _run_bbmap(ref_path_i, output_path=output_path_i, reads_path_1=reads_path_1, reads_path_2=reads_path_2)
        output_paths.append(output_path_i)

    # Need to store the iteration to properly length-normalize reads for the graph.
    df = pd.concat([BamFile.from_file(path).to_df(include_flags=FLAGS['read_paired']).assign(iteration=i) for i, path in enumerate(output_paths)])
    _write_reads_to_fasta(df, output_path=os.path.join(output_dir, 'reads.fasta'))
    return df

