import pandas as pd 
from Bio.Seq import Seq
import numpy as np 


START_CODONS = ['ATG', 'GTG', 'TTG']
STOP_CODONS = ['TAA', 'TAG', 'TGA']

get_codons = lambda nt_seq : np.array([nt_seq[i:i + 3] for i in range(0, len(nt_seq), 3)])
get_reverse_complement = lambda seq : str(Seq(seq).reverse_complement())
get_seq = lambda nt_seq : str(Seq(nt_seq).translate(to_stop=False))


def orf_check(nt_seq:str):
    ''''''
    codons = get_codons(nt_seq)
    assert len(nt_seq) % 3 == 0, 'orf_check: The ORF should have a length which is a multiple of three.'
    assert codons[0] in START_CODONS, f'orf_check: The ORF does not have a valid start codon: {codons[0]}.'
    assert codons[-1] in STOP_CODONS, f'orf_check: The ORF does not have a valid stop codon: {codons[-1]}.'



def _get_orfs(nt_seq:str, min_length:int=30, start_codons:list=START_CODONS):
    ''''''
    codons = get_codons(nt_seq)     
    start_codon_idxs = np.where(np.isin(codons, start_codons))[0].tolist()

    rows = list()
    for i in start_codon_idxs:
        row = {'start_codon':codons[i]}
        row['nt_seq'] = ''.join(codons[i:])
        row['seq'] = get_seq(row['nt_seq'])
        if len(row['seq']) < min_length:
            continue 
        rows.append(row)

    return rows


def orf_adjust_coordinates():
    pass


def orf_get_alternates(start:int, stop:int, strand:int, extension_length:int=30, contig:str=None, frame_shift:int=0, min_length:int=30, start_codons:list=START_CODONS):
    '''Extend the gene between start and stop at the N-terminus by the specified number amino acids, ignoring any in-frame stop codons.
    Then obtain the corresponding extended sequences, as well as all in-frame subsequences within the extended region starting from a valid start codon.

    :param start: The one-indexed (inclusive) starting coordinate of the ORF. 
    :param stop: The one-indexed (inclusive) ending coordinate of the ORF.
    :param strand: The strand of the ORF, either - or +.
    :param contig: The contig on which the gene is found.
    :param extension_length: The number of amino acids to extend the gene by. 
    :param frame_shift: If not zero, this specifies how much to shift the reading frame when translating the extended N-terminus. Specifically, it 
        shifts the N-terminal extension stop coordinate to the left by frame_shift.
    :param min_length: The minumum length of the final sequence, in units of amino acids. 
    :param start_codons: The list of start codons to consider when looking for alternate start sites. 
    :returns:
    '''
    assert ((stop - start) % 3) == 2, 'orf_get_alternates: Gene coordinates are not one-indexed inclusive.'

    start = start - 1

    if strand == '-':
        old_start, old_stop = start, stop
        start = (len(contig) - old_stop) 
        stop = (len(contig) - old_start)
        contig = get_reverse_complement(contig)

    extension_stop = start - frame_shift # Extended region will end at the start position. Apply a frame shift, if specified.
    extension_start = extension_stop - (3 * extension_length)

    assert extension_start > 0,  f'orf_get_alternates: Extension start coordinate cannot be negative, so {extension_start} is invalid.'

    nt_seq = contig[extension_start:extension_stop] + contig[start:stop]
    
    assert (len(contig[start:stop]) % 3) == 0, f'orf_get_alternates: Expected the original ORF to have a length divisible by 3.'
    assert contig[start:start + 3] in START_CODONS, f'orf_get_alternates: Expected the region to start with a start codon, but got {contig[start:start + 3]}.'
     
    df = list()
    df += _get_orfs(nt_seq, min_length=min_length, start_codons=start_codons)
    df = pd.DataFrame(df)
    df['frame_shift'] = frame_shift
    df['original_length'] = (stop - start) // 3
    df['length'] = df.seq.str.replace(r'\*', '').apply(len)
    return df


# It seems as though many of the longer ORFs have alternate start sites which would truncate the ORF to around 80-90 residues.
# However, most of the longer ones do not have an ATG start codon that would shift the length to ~60 (with the 
# possible exception of SR-VP_9_9_2021_59_4A_0_85m_scaffold_34976_2). This suggests that if there was an error in ORF boundary
# prediction, it is more likely that the 60 residue genes were N-terminally truncated.

# def orf_get_n_terminal_truncation(start:int, stop:int, strand:int, min_length:int=30, contig:str=None, start_codons:list=START_CODONS):
#     '''Look for alternate start codons for the predicted ORF.

#     :param start: The one-indexed (inclusive) starting coordinate of the ORF. 
#     :param stop: The one-indexed (inclusive) ending coordinate of the ORF.
#     :param strand: The strand of the ORF, either - or +.
#     :param contig: The contig on which the gene is found.
#     :param min_length: The minumum length of the final sequence, in units of amino acids. 
#     :param start_codons: The list of start codons to consider when looking for alternate start sites. 
#     :return: 
#     '''

#     if strand == '-': # Reverse complement if the gene is on the opposite strand.
#         contig = get_reverse_complement(contig[start - 1:stop])

#     nt_seq = contig[start:stop]
#     orf_check(nt_seq)

#     df = _get_orfs(nt_seq, min_length=min_length, start_codons=start_codons)
#     df = pd.DataFrame(df)
#     return df 




# for row in ggkbase_df.sort_values('contig_id').itertuples():
#     if row.contig_id not in contigs:
#         continue 

#     ext, codons = extend(row.start, row.stop, row.strand, contig=contigs[row.contig_id], frame_shift=0)
#     has_start_codon = np.any(np.isin(['ATG'], codons))

#     # Get the potential start codons which are not followed by a stop codon (so possible valid extensions)
#     stop_codon_idxs = np.where(np.isin(codons, STOP_CODONS))[0]
#     start_codon_idxs = np.where(np.isin(codons, START_CODONS))[0]
#     start_codon_idxs = start_codon_idxs[start_codon_idxs > stop_codon_idxs.max()]

#     if len(start_codon_idxs) > 0:
#         seq_extended = ext[start_codon_idxs.min():] + row.seq
#         ggkbase_df.loc[row.Index, 'seq_extended'] = seq_extended
#         print(f'Extended {row.Index} ({len(row.seq)} aa to {len(seq_extended)} aa)\t{seq_extended}')
#     else:
#         ggkbase_df.loc[row.Index, 'seq_extended'] = row.seq
#         print(f'No extension for {row.Index} ({len(row.seq)} aa)')

#     # print(row.Index, f'({len(row.seq)} aa)\tHas start codon? {has_start_codon}')
#     # print(ext, row.seq, end='\n\n')

# # It seems unlikely that the 61-residue ones cannot be extended (I tried all three frames). However, the proteins in the 88-residue
# # range can be reasonably extended. If these is indeed a signal peptide, an N-terminal extension would add a run of positive residues that
# # might make the signal peptide prediction more robust. 

# # See https://pmc.ncbi.nlm.nih.gov/articles/PMC2323981/ for information on archaeal signal peptides.
# # Betazoid signal peptides should resemble those of archaea (as they would be using host secretion machinery.)
    