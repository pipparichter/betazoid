import Bio.PDB as PDB
from Bio.PDB import PDBParser, MMCIFParser, Structure, Model, Chain, Select
from scipy.spatial.distance import cdist
import pandas as pd 
import numpy as np 

class Select(PDB.Select):

    def __init__(self, start:int, stop:int, chain_id:str='A'):
        self.chain_id = chain_id
        self.start = start
        self.stop = stop 

    def accept_chain(self, chain):
        return self.chain_id == chain.id 

    def accept_residue(self, residue):
        accept = (self.chain_id == residue.get_parent().id)

        i = residue.id[1] # residue.id is a 3-tuple, and the middle value is the residue number. 
        accept = accept and ((i > self.start) and (i <= self.stop)) # Account for the fact that the residue number is 1-indexed.
        return accept 


is_residue_atp = lambda residue : residue.id[0] == 'H_ATP'
is_residue_mg = lambda residue : residue.id[0] == 'H_MG'

is_chain_atp = lambda chain : (len(list(chain.get_residues())) == 1) and (is_residue_atp(list(chain.get_residues())[0]))
is_chain_mg = lambda chain : (len(list(chain.get_residues())) == 1) and (is_residue_mg(list(chain.get_residues())[0]))

class PDBFile():

    def __init__(self):
        pass 

    @classmethod
    def from_file(cls, path, format:str=None):
        
        format = path.split('.')[-1] if (format is None) else format
        assert format in ['cif', 'pdb'], 'PDBFile.from_file: Format must be either pdb or cif.'
        
        parser = MMCIFParser(QUIET=True) if (format == 'cif') else PDBParser(QUIET=True)
        structure = parser.get_structure('protein', path)[0]

        obj = cls()
        obj.path = path  
        obj.structure = structure
        return obj 
    
    def get_chain_ids(self):
        return [chain.get_id() for chain in self.structure.get_chains() if (not is_chain_atp(chain)) and (not is_chain_mg(chain))]
        # else:
        #     return [chain.get_id() for chain in self.structure.get_chains()]
    
    def to_df(self, atoms:list=['NZ', 'CZ', 'CA', 'CB'], residues=['LYS', 'ARG', 'PHE', 'TYR']):

        df = list()
        for chain in self.structure.get_chains():
            for residue in chain.get_residues():
                if residue.get_resname() not in residues:
                    continue
                row = {'chain_id':chain.get_id(), 'residue':residue.get_resname()}
                row['residue_num'] = residue.get_id()[1]
                row.update({atom:(residue[atom].coord if (atom in residue) else None) for atom in atoms})
                df.append(row)
        df = pd.DataFrame(df)
        return df

    def write(self, path:str):
        io = PDB.PDBIO()
        io.set_structure(self.structure)
        io.save(path)

    def subset(self, output_path:str, start:int=None, stop:int=None, chain_id:str='A'):
        assert start is not None, f'PDBFile.subset: The start and stop indices need to be numbers, got {type(start)}.'
        assert stop is not None, f'PDBFile.subset: The start and stop indices need to be numbers, got {type(start)}.'
        io = PDB.PDBIO()
        io.set_structure(self.structure)
        io.save(output_path,  Select(start, stop, chain_id=chain_id))

