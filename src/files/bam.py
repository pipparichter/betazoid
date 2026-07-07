import io 
import pandas as pd 
import numpy as np 
import subprocess
import re 


# BAM files contain the following information. 
# 1. RNAME: The name of the reference contig the read is aligned to. 
# 2. POS: Leftmost mapping position of the read on the reference. 
# 3. MAPQ: Mapping quality score, which is confidence in the placement. 
# 4. CIGAR: Compact representation of how the read aligns to the reference (matches, insertions, deletions)
# 5. FLAGS: Binary flags that encode information such as strand, paired-end status, whether the read is mapped, etc.


# -f option includes reads with flag, -F excludes reads with flag. 
FLAGS = dict()
FLAGS['read_paired'] = 1
FLAGS['proper_pair'] = 2
FLAGS['unmapped'] = 4
FLAGS['mate_unmapped'] = 8 # Read is mapped but pair is not. 
FLAGS['reverse_strand'] = 16
FLAGS['mate_reverse_strand'] = 32
FLAGS['read_1'] = 64
FLAGS['read_2'] = 128
FLAGS['secondary'] = 256 # Each read should have exactly one primary alignment in the BAM file. If this flag is set, alignment is secondary. Can't have a read with only a chimeric alignment. 
FLAGS['qc_fail'] = 512
FLAGS['duplicate'] = 1024 # PCR duplicates, safe to exclude. 
FLAGS['supplementary'] = 2048 # Flag is set when a portion of the read aligns to a different location (split read). Primary alignments will not have this set. This can occur if reads span structural variants (SVs) like deletions, inversions, or translocations

CIGAR_CODES = ['I', 'D', 'X', '=', 'S', 'H']
CIGAR_PATTERN = r'(\d+)(I|D|X|=|S|H)'

def parse_cigar(cigar:str):
    '''
    
    (1) H: Hard clips, ... Clips only occur at read ends. 
    (2) S: Soft clips, ... Clips only occur at read ends. 
    (3) M/=: Matches
    (4) I
    (5) D
    (6) X 
    '''    
    info = {code:0 for code in CIGAR_CODES}

    left_clip_length, right_clip_length = 0, 0

    alignment = list()
    for i, (n, code) in enumerate(re.findall(CIGAR_PATTERN, cigar)):
        info[code] += int(n)

        if (code == '=') or (code == 'X'):
            alignment += [1] * int(n)
        if (code == 'I'):
            continue
        if (code == 'D'):
            alignment += [0] * int(n) 
        if (code == 'H') or (code == 'S'):
            left_clip_length = int(n) if (i == 0) else left_clip_length
            right_clip_length = int(n) if (i > 0) else right_clip_length 

    info['n_errors'] = info['X'] + info['D'] + info['I'] 
    info['n_clips'] = info['S'] + info['H'] 
    info['ref_alignment_length'] = info['='] + info['X'] + info['D'] # Reference-consuming operations. 
    info['read_alignment_length'] = info['='] + info['X'] + info['I'] + info['S'] # Read-consuming operations. 
    info['left_clip_length'] = left_clip_length
    info['right_clip_length'] = right_clip_length
    return info, np.array(alignment)


class BamFile():
    # SAM file will always have 11 columns, but can have a variable number of tags after the first 11 (which I ignore here).
    fields = ['read_id', 'flag', 'ref_id', 'position', 'mapping_quality', 'cigar', 'mate_ref_id', 'mate_position', 'template_length', 'seq', 'quality_string']

    def __init__(self, path:str=None):
        self.path = path  

    @classmethod
    def from_file(cls, path:str):
        return cls(path=path)
    
    @staticmethod
    def _parse_flags(flags):

        metadata = dict()
        for label, code in FLAGS.items():
            metadata[label] = [bool(code & flag) for flag in flags]
        return metadata
    
    def _read_header(self):
        '''Read the BAM file header, which contains the reference length information.'''
        cmd = ['samtools', 'view', '-H', self.path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        
        df = list()
        for line in result.split('\n'):
            if line.startswith('@SQ'):
                info = dict()
                info['ref_id'] = re.search(r'SN:([^\s]+)', line).group(1)
                info['ref_length'] = re.search(r'LN:([\d]+)', line).group(1)
                df.append(info)
        df = pd.DataFrame(df).astype({'ref_id':'string', 'ref_length':'int'})
        return df 
    

    def _read(self, include_flags:int=None, exclude_flags:int=None):
        cmd = ['samtools', 'view']
        if include_flags is not None:
            cmd += ['-f', str(include_flags)]
        if exclude_flags is not None:
            cmd += ['-F', str(exclude_flags)]
        cmd += [self.path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        df = pd.read_csv(io.StringIO(result), sep='\t', header=None, names=BamFile.fields, usecols=range(len(BamFile.fields)))
        df = df.astype({'ref_id':'string'}) 
        df = df.merge(self._read_header(), on='ref_id', how='left') # Add header information. 
        return df 
    

    @staticmethod
    def _get_pair_orientation(df:pd.DataFrame):
        strands = np.select([(~df.unmapped & df.reverse_strand), (~df.unmapped & ~df.reverse_strand)], ['R', 'F'], default='X')
        mate_strands = np.select([(~df.mate_unmapped & df.mate_reverse_strand), (~df.mate_unmapped & ~df.mate_reverse_strand)], ['R', 'F'], default='X')
        return  [f'{a}{b}' for a, b in zip(strands, mate_strands)]


    def to_alignments(self, max_n_errors:int=5):
        
        df = self._read(exclude_flags=4) # Exclude unmapped reads. 
        alignments = {row.ref_id:np.zeros(int(row.ref_length)) for row in self._read_header().itertuples()}
        
        is_aligned_to_left = lambda info : (info['position'] == 1)
        is_aligned_to_right = lambda info : (info['position'] + info['ref_alignment_length']) > info['ref_length']
        has_internal_clips = lambda info : ((info['left_clip_length'] > 0) and not is_aligned_to_left(info)) or (info['right_clip_length'] > 0) and not is_aligned_to_right(info)
        is_invalid_alignment = lambda info : (info['n_errors'] > max_n_errors) or has_internal_clips(info)

        n = 0
        for ref_id, df_ in df.groupby('ref_id'):    
            for row in df_.to_dict(orient='records'): # Iterate rows as dictionaries. 
                info, alignment = parse_cigar(row['cigar']) # Get the alignment array and info. 
                info.update(row)

                if is_invalid_alignment(info):
                    n += 1
                    continue

                start_position = info['position'] - 1
                stop_position = info['position'] - 1 + len(alignment)
                alignments[ref_id][start_position:stop_position] += alignment

        print(f'BamFile.to_alignments: Discarded {n} alignments with internal clips or max_n_errors > {max_n_errors}.')
        return alignments 


    def to_df(self, include_flags:int=None, exclude_flags:int=None):

        df = self._read(include_flags=include_flags, exclude_flags=exclude_flags)

        # Parse the flags and add stored metadata to the DataFrame. 
        for col, data in BamFile._parse_flags(df.flag).items():
            df[col] = data 
        df['orientation'] = BamFile._get_pair_orientation(df) # Assign orientation to each pair.
        df['read_number'] = np.select([(df.read_paired & ~df.read_1), (df.read_paired & df.read_2), (~df.read_paired)], [1, 2, 0], default=0)

        df = pd.concat([df, pd.DataFrame([parse_cigar(cigar)[0] for cigar in df.cigar])], axis=1)
        return df 
    

    def to_fasta(self, path:str=None, include_flags:int=None, exclude_flags:int=None):
        # fasta_path = os.path(OUTPUT_DIR, f'{JOB_NAME}.reads.{i}.fasta')
        path = self.path.replace('.bam', '.fasta') if (path is None) else path

        cmd = self._read('fasta', include_flags=include_flags, exclude_flags=exclude_flags)
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        n = result.stdout.count('>') # Number of sequences in the FASTA file. 
        with open(path, 'w') as f:
            f.write(result.stdout)
        return path, n
        

# def bam_to_fastq(bam_path:str):
#     fastq_path_paired  = bam_path.str.replace('.bam', '.fastq.paired') # For orphan reads. 
#     fastq_path_orphan  = bam_path.str.replace('.bam', '.fastq.orphan') 

#     cmd = ['samtools', 'fastq', '-i'] # -i option means to interleave paired reads.
#     cmd += ['-0', fastq_path_orphan]
#     cmd += ['-o', fastq_path_paired]
#     subprocess.run(cmd,check=True)
#     return fastq_path_paired, fastq_path_orphan

    # def _get_strands(df):
    #     '''Note that if both mate pairs are mapped, but the orientation is the same, the strand will be set to None 
    #     to flag the ambiguity in the mapping.'''
    #     conditions = [(df.unmapped & df.mate_reverse_strand) | (~df.unmapped & (~df.reverse_strand))]
    #     conditions += [(df.unmapped & ~df.mate_reverse_strand) | (~df.unmapped & (df.reverse_strand))]
    #     return np.select(conditions, ['+', '-'], default=None)
    