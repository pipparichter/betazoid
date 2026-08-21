import sys
import numpy as np 
import pandas as pd

HOLE_RAD_PATH = '../data/genes/hole/rad/simple.rad'


# SAMPLE specifies the distance between the planes used in the HOLE procedure. The default value is 0.25 A, this should be reasonable for most purposes; if you wish to visualize a very tight 
# constriction then specify a smaller value. This essentially increases the resolution along the pore. 

HOLE_INPUT = '''coord {pdb_path}
radius {rad_path}
Conserved motifs
endrad 40
sample 0.1
sphpdb {output_path} 
cpoint {cpoint}
cvect {cvect} 
'''

# ! Everything preceded by a "!" is a comment and will be ignored by HOLE
# ! https://www.holeprogram.org/doc/index.html 
# ! https://www.holeprogram.org/doc/old/hole_d03.html Contains more documentation details. 

# ! The arguments below are required. 
# coord {path} ! Path to the PDB file. 
# radius ~/hole2/rad/simple.rad ! File that contains information about the Van Der Waals radii for different residues. 

# ! The arguments below are optional. 
# sphpdb {output_path} ! PDB output file. 
# endrad 5. ! Specifies the radius above which the program regards a result as an indicating that the end of the pore has been reached. The default value is 15.0 angstroms.
# cpoint {x} {y} {z} ! The coordinates for the pore start site. 
# cvect {x} {y} {z} !  vector which lies in the direction of the channel/pore


def hole_clean_pdb(path:str):
    '''HOLE does not support non-amino acid residues in the PDB file (e.g. Mg or ATP), so these need
    to be removed prior to running the program.'''

    assert 'alphafold' not in path, f'hole_clean_pdb: Make sure you are not overwriting an original PDB file dummy, {path}'

    with open(path, 'r') as f:
        lines = f.readlines()
        lines = [line for line in lines if not (line.startswith('HETATM'))] # Remove all non-amino acid entries.
    with open(path, 'w') as f:
        f.write('\n'.join(lines))
    return path


def hole_get_axis(df:pd.DataFrame, start_residue_num:int=1, end_residue_num:int=220, atom:str='CA', chain_ids:list=None):
    '''Find the start position (CPOINT) and direction (CVECT) of the pore through a cluster_3 ATPase. 
    This code relies on the multimer being symmetric and a homo-oligomer.
    
    :param df: Path to the PDB converted to a DataFrame. 
    :param start_residue_num:
    :param end_residue_num:
    :param chain_ids
    '''
    

    if chain_ids is not None:
        df = df[df.chain_id.isin(chain_ids)].copy()
        print('get_hole_axis: Computing the pore axis using chains', ', '.join(chain_ids))

    coords = dict()
    coords['start'] = np.array(df[df.residue_num == start_residue_num][atom].tolist())
    coords['end'] = np.array(df[df.residue_num == end_residue_num][atom].tolist())

    # Taking the average over corresponding residues in all chains should be the centroid, or where the pore lies. 
    cpoints = (coords['start'].mean(axis=0), coords['end'].mean(axis=0))
    cvect = cpoints[1] - cpoints[0]

    return cpoints[0], cvect


def hole_check_input(input, max_line_length:int=90):
    '''HOLE is old Fortran code, and some versions have fixed-length character buffers (often 80 or 255 characters). Long filenames 
    or paths can cause silent failures.'''
    lines = input.split('\n')
    for line in lines:
        assert len(line) <= max_line_length, f'check_hole_input: Found a line that is {len(line)} characters long.\n{line}'


def hole_write_input(input_path, rad_path=HOLE_RAD_PATH, **kwargs):
    '''Create a HOLE input file using the specified keyword arguments.'''
    input = HOLE_INPUT.format(rad_path=rad_path, **kwargs)

    hole_check_input(input)
    hole_clean_pdb(kwargs.get('pdb_path'))

    with open(input_path, 'w') as f:
        f.write(input)
