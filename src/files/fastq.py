# class FastqFile():

#     def __init__(self, paired:bool=False, seqs:list=[], ids:list=[]):
#         self.seqs = seqs 
#         self.ids = ids 
#         self.paired = paired 
    
#     @classmethod
#     def from_file(cls, path:str, paired:bool=True):

#         seqs, ids = list(), list()
#         if paired:
#             iterator = SeqIO.parse(path, 'fastq')
#             while iterator:
#                 record_1, record_2 = next(iterator), next(iterator)
#                 seqs.append((str(record_1.seq), str(record_2.seq)))
#                 ids.append((record_1.id, record_2.id))
#         else:
#             for record in SeqIO.parse(path, 'fastq'):
#                 seqs.append(str(record.seq))
#                 ids.append(record.id)
#         return cls(paired=paired, seqs=seqs, ids=ids)
    

#     def to_df(self):

#         df = list()
#         if self.paired:
#             pass 
#         return df