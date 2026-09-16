def get_random_seqs_gc_content(nt_seq, n:int=100, seed:int=42):
    '''Generates n random nucleotide sequences of the same length and GC content as the input sequence. 
    
    :param nt_seq: The nucleotide sequence to use as reference.
    :param n: The number of sequences to generate. 
    '''
    if seed is not None:
        np.random.seed(seed)

    gc_content = get_gc_content(nt_seq)
    length = len(nt_seq)
    gc_count = int(length * gc_content)

    nt_seqs = list()
    for _ in range(n):
        # Randomly select a new set of bases with each iteration. 
        bases = np.random.choice(['G', 'C'], size=gc_count, replace=True).tolist()
        bases += np.random.choice(['A', 'T'], size=length - gc_count, replace=True).tolist()
        np.random.shuffle(bases) # Modifies inplace.
        nt_seqs.append(''.join(bases))

    return nt_seqs


# AMINO_ACID_TO_CODONS_MAP = {'A':['GCU','GCC','GCA','GCG'],'R':['CGU','CGC','CGA','CGG','AGA','AGG'],'N':['AAU','AAC'],'D':['GAU','GAC'],'C':['UGU','UGC'],'E':['GAA','GAG'],'Q':['CAA','CAG'],'G':['GGU','GGC','GGA','GGG'],'H':['CAU','CAC'],'I':['AUU','AUC','AUA'],'L':['UUA','UUG','CUU','CUC','CUA','CUG'],'K':['AAA','AAG'],'M':['AUG'],'F':['UUU','UUC'],'P':['CCU','CCC','CCA','CCG'],'S':['UCU','UCC','UCA','UCG','AGU','AGC'],'T':['ACU','ACC','ACA','ACG'],'W':['UGG'],'Y':['UAU','UAC'],'V':['GUU','GUC','GUA','GUG'],'*':['UAA','UAG','UGA']}
AMINO_ACID_TO_CODONS_MAP = {'A':['GCT','GCC','GCA','GCG'],'R':['CGT','CGC','CGA','CGG','AGA','AGG'],'N':['AAT','AAC'],'D':['GAT','GAC'],'C':['TGT','TGC'],'E':['GAA','GAG'],'Q':['CAA','CAG'],'G':['GGT','GGC','GGA','GGG'],'H':['CAT','CAC'],'I':['ATT','ATC','ATA'],'L':['TTA','TTG','CTT','CTC','CTA','CTG'],'K':['AAA','AAG'],'M':['ATG'],'F':['TTT','TTC'],'P':['CCT','CCC','CCA','CCG'],'S':['TCT','TCC','TCA','TCG','AGT','AGC'],'T':['ACT','ACC','ACA','ACG'],'W':['TGG'],'Y':['TAT','TAC'],'V':['GTT','GTC','GTA','GTG'],'*':['TAA','TAG','TGA']}
CODON_TO_AMINO_ACID_MAP = {codon:aa for aa, codons in AMINO_ACID_TO_CODONS_MAP.items() for codon in codons}



def get_codon_usage(nt_seqs):
    '''Get the codon usage frequencies for a particular set of genes'''
    counts = {codon:0 for codon in CODON_TO_AMINO_ACID_MAP.keys()}
    for nt_seq in nt_seqs:
        for codon, n in zip(*np.unique(get_codons(nt_seq), return_counts=True)):
            counts[str(codon)] += n

    df = pd.DataFrame.from_dict(counts, orient='index', columns=['count']).reset_index(names='codon')
    df['amino_acid'] = df.codon.map(CODON_TO_AMINO_ACID_MAP)
    df = df[df.amino_acid != '*'].copy()
    # df['frequency'] = df.groupby('amino_acid').transform(lambda df :print(df))
    df['total'] = df.groupby('amino_acid')['count'].transform('sum')
    df['frequency'] = df['count'] / df['total']
    df = df.fillna(0)
    return df.set_index('codon').frequency.to_dict()



def get_random_seqs_amino_acids(seq, n:int=1000, codon_usage:dict=dict(), seed:int=42):
    '''Generates n random nucleotide sequences constrained by the given amino acid sequence. If a codon_usage
    dictionary is provided, then codons are selected with a probability concordant with the background.
    '''
    if seed is not None:
        np.random.seed(seed)

    seqs = list()
    for aa in seq: # Iterate over each amino acid and store possible codons in a (len(seq), n)-dimensional array.
        codons = AMINO_ACID_TO_CODONS_MAP[aa] # Get all possible codons for the amino acid.
        frequencies = [codon_usage.get(codon, 1/len(codons)) for codon in codons]
        seqs.append(np.random.choice(codons, p=frequencies, size=n, replace=True))
    seqs = np.array(seqs).T
    seqs = [''.join(seq) for seq in seqs]
    return seqs
