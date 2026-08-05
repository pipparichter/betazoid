import Bio.PDB as PDB
from Bio.PDB import PDBParser, MMCIFParser, Select 
from Bio.PDB import Structure, Model, Chain, Residue, Atom
from Bio.PDB.PDBExceptions import PDBConstructionException
from scipy.spatial.distance import cdist
import pandas as pd 
import numpy as np
import os
from tqdm import tqdm
import re

AMINO_ACID_CODES = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
ATOMS = ['N', 'CA', 'C', 'O', 'OXT', 'CB', 'CG', 'CG1', 'CG2', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CZ', 'CZ2', 'CZ3', 'CH2', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'SD', 'SG']


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

# NOTE: A traditional PDB file cannot represent multi-character chain IDs; this requires MMCIF. 

class SPHParser():
    '''Parser for HOLE output files, which have the .sph extension.'''

    # PDB files have a fixed-width column format, so these indices always match the particular field. 
    fields = dict()
    fields['record'] = (0, 6)
    fields['residue_num'] = (22, 26) # This will be the position along the pore. 
    fields['x'] = (30, 38)
    fields['y'] = (38, 46)
    fields['z'] = (46, 54)
    # fields['occupancy'] = (54, 60) # In HOLE output, this stores the radius. 
    fields['b_factor'] = (60, 66) # In HOLE output, this stores the radius. 

    dtypes = dict()
    dtypes['x'] = float 
    dtypes['y'] = float 
    dtypes['z'] = float 
    dtypes['residue_num'] = int 
    dtypes['b_factor'] = float 

    def __init__(self, name:str='none', QUIET:bool=True):
        self.structure = Structure.Structure(name)
        self.model =  Model.Model(0)
        self.chain = Chain.Chain('S') # HOLE output calls the chain S. 
        self.residues = dict()
        self.build = True

        self.n_atoms = 0 # For assigning serial numbers to the atoms, which is needed because all atoms have a serial number of 1 in the HOLE output.

    @staticmethod
    def _parse_line(line:str):
        '''Parse a line in the SPH file, ensuring the output types are correct.'''
        data = dict()
        for field, (start, stop) in SPHParser.fields.items():
            value = line[start:stop].strip()
            if len(value) == 0: # This happens if all whitespace is removed by strip. 
                continue 
            dtype = SPHParser.dtypes.get(field, str)
            data[field] = dtype(value)
        assert isinstance(data['residue_num'], int), f'SPHParser._parse_line: residue_num should be an integer, but got ' + str(type(data['residue_num']))
        return data 

    def _get_atom(self, line:str):
        '''Takes a line of the SPH file as input, and returns an Atom object, as well as the residue number.'''
        line = SPHParser._parse_line(line)

        kwargs = {'name':'QSS', 'fullname':' QSS'}
        kwargs['bfactor'] = line['b_factor']
        kwargs['occupancy'] = 1.0
        kwargs['serial_number'] = self.n_atoms
        kwargs['altloc'] = ' '
        kwargs['coord'] = np.array([line['x'], line['y'], line['z']])
        kwargs['element'] = 'C' # This is just a placeholder element. 

        self.n_atoms += 1 # Increment the serial number. 
        atom = Atom.Atom(**kwargs)
        return line['residue_num'], atom

    def _get_residue(self, residue_num:int):
        if residue_num in self.residues:
            return self.residues[residue_num]
        else:
            # The residue ID is typically (het_flag, residue_num, i_code). There is no het_flag because these are pseudo-atoms. 
            residue = Residue.Residue((' ', residue_num, ' '), 'SPH', '')
            self.residues[residue_num] = residue 
            return residue

    def _parse(self, path:str):

        with open(path, 'r') as f:
            lines = f.readlines()

        # for line in tqdm(lines, desc='SPHParser._parse'):
        for line in lines:
            if line.startswith('LAST-REC-END'):
                continue

            residue_num, atom = self._get_atom(line)
            if residue_num == -888: # Not sure what these are for, but they don't seem to be real parts of the channel. 
                continue 

            residue = self._get_residue(residue_num)
            try:
                residue.add(atom)
            except PDBConstructionException:
                continue 
                # print(f'SPHParser._parse: Skipping residue {residue_num}, which has already been assigned an Atom object.')

    def _build(self):
        ''''''
        self.build = False

        self.structure.add(self.model)
        self.model.add(self.chain)
        for residue in self.residues.values():
            self.chain.add(residue)


    def get_structure(self, type_:str, path:str):

        if self.build:
            self._parse(path)
            self._build()

        return self.structure


class PDBFile():

    parsers  = dict()
    parsers['cif'] = MMCIFParser
    parsers['pdb'] = PDBParser
    parsers['sph'] = SPHParser

    amino_acid_codes = AMINO_ACID_CODES

    def __init__(self):
        self.format = None

    @classmethod
    def from_file(cls, path, format:str=None):
        
        format = path.split('.')[-1] if (format is None) else format.lower()
        assert format in ['cif', 'pdb', 'sph'], 'PDBFile.from_file: Format must be either pdb, cif, or sph.'
        
        parser = PDBFile.parsers[format](QUIET=True)
        structure = parser.get_structure('protein', path)[0]

        obj = cls()
        obj.format = format
        obj.path = path  
        obj.structure = structure
        return obj 
    
    def get_chain_ids(self):
        return [chain.get_id() for chain in self.structure.get_chains() if (not is_chain_atp(chain)) and (not is_chain_mg(chain))]
        # else:
        #     return [chain.get_id() for chain in self.structure.get_chains()]
    
    def to_df(self, atoms:list=ATOMS, residues=AMINO_ACID_CODES):
        '''Converts the data in the PDB file to a DataFrame with one entry per residue (not per atom). The DataFrame has the following columns:
            (1) chain_id
            (2) residue_num
            (3) residue: The residue code, e.g. LYS, SPH, ARG, etc. 
            (4) b_factor
            (5+) {atom}: Stores the x, y, z coordinate of an atom associated with the residue. 

        :param atoms: 
        :param residues: 
        '''
        df = list()
        for chain in self.structure.get_chains():
            for residue in chain.get_residues():
                if residue.get_resname() not in residues:
                    continue
                row = {'chain_id':chain.get_id(), 'residue':residue.get_resname()}
                row['residue_num'] = residue.get_id()[1]
                row.update({atom:(residue[atom].coord if (atom in residue) else None) for atom in atoms})

                b_factor = [residue[atom].bfactor for atom in atoms if (atom in residue)]
                # assert len(b_factors) == 1, f'PDBFile.to_df: Expected one b factor per residue, but got {b_factors}.' 
                row['b_factor'] = np.mean(b_factor)

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



def cif_to_pdb(path:str, output_path:str=None, output_dir:str=None, overwrite:bool=False):
    '''Convert a CIF file to a PDB file.
    
    '''
    file_name, ext = os.path.splitext(os.path.basename(path))
    assert ext == '.cif', f'cif_to_pdb: Expected the input file to have the .cif extension, but got {ext}.'

    if output_path is None:
        output_dir = os.path.dirname(path) if (output_dir is None) else output_dir
        output_path = os.path.join(output_dir, file_name + '.pdb')

    if (not os.path.exists(output_path)) or overwrite:
        file = PDBFile.from_file(path, format='cif')
        file.write(output_path)
    return output_path


def sph_to_pdb(path:str, output_dir:str=None, overwrite:bool=False):
    '''Convert a SPH file to a PDB file.
    
    '''
    file_name, ext = os.path.splitext(os.path.basename(path))
    assert ext == '.sph', f'sph_to_pdb: Expected the input file to have the .sph extension, but got {ext}.'
    output_dir = os.path.dirname(path) if (output_dir is None) else output_dir
    output_path = os.path.join(output_dir, file_name + '.pdb')

    if (not os.path.exists(output_path)) or overwrite:
        file = PDBFile.from_file(path, format='sph')
        file.write(output_path)
    return output_path


