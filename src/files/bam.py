import io 
import pandas as pd 
import numpy as np 
import subprocess
from src.files import FASTAFile


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
    
    def _get_samtools_cmd(self, subcommand:str, include_flags:int=None, exclude_flags:int=None, options:list=[]):
        cmd = ['samtools', subcommand] + options
        if include_flags is not None:
            cmd += ['-f', str(include_flags)]
        if exclude_flags is not None:
            cmd += ['-F', str(exclude_flags)]
        cmd += [self.path]
        return cmd 
    
    @staticmethod
    def _get_pair_orientation(df:pd.DataFrame):
        strands = np.select([(~df.unmapped & df.reverse_strand), (~df.unmapped & ~df.reverse_strand)], ['R', 'F'], default='X')
        mate_strands = np.select([(~df.mate_unmapped & df.mate_reverse_strand), (~df.mate_unmapped & ~df.mate_reverse_strand)], ['R', 'F'], default='X')
        return  [f'{a}{b}' for a, b in zip(strands, mate_strands)]
        
    def to_df(self, include_flags:int=None, exclude_flags:int=None):

        cmd = self._get_samtools_cmd('view', include_flags=include_flags, exclude_flags=exclude_flags)
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        df = pd.read_csv(io.StringIO(result.stdout), sep='\t', header=None, names=BamFile.fields, usecols=range(len(BamFile.fields)))
        # Parse the flags and add stored metadata to the DataFrame. 
        for col, data in BamFile._parse_flags(df.flag).items():
            df[col] = data 
        df['orientation'] = BamFile._get_pair_orientation(df) # Assign orientation to each pair.

        return df 
 

    def to_fasta(self, path:str=None, include_flags:int=None, exclude_flags:int=None):
        # fasta_path = os.path(OUTPUT_DIR, f'{JOB_NAME}.reads.{i}.fasta')
        path = self.path.replace('.bam', '.fasta') if (path is None) else path

        cmd = self._get_samtools_cmd('fasta', include_flags=include_flags, exclude_flags=exclude_flags)
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
    