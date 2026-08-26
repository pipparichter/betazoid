import sys
import numpy as np 
import pandas as pd
from files.structure import StructureFile, cif_to_pdb, sph_to_pdb
import os 
import subprocess
from scipy.spatial.distance import pdist, cdist

HOLE_RAD_PATH = '../data/genes/hole/rad/simple.rad'


# SAMPLE specifies the distance between the planes used in the HOLE procedure. The default value is 0.25 A, this should be reasonable for most purposes; if you wish to visualize a very tight 
# constriction then specify a smaller value. This essentially increases the resolution along the pore. 

HOLE_INPUT = '''coord {pdb_path}
radius {rad_path}
Conserved motifs
endrad 40
sample 0.2
sphpdb {output_path} 
cpoint {cpoint}
cvect {cvect} 
'''
# Something seems to have gone wrong with the bz_3 cluster_3 multimer HOLE prediction when I did an initial run, as the channel seems to have stopped partway through the pore. 
# I am unsure why, as the pore did not seem especially narrow, although perhaps the issue was that the ENDRAD threshold was met). I decreased the
# SAMPLE parameter from the default (0.25) to 0.1 and increased ENDRAD parameter from 25 to 40. This produced much better results; HOLE captured the complete
# channel, and also was more consistent across runs. I then decided to re-run HOLE on all PDB structures using a larger ENRAD and decreased SAMPLE. 


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
    This code relies on all chains specified by chain_ids being identical.
    
    :param df: Path to the PDB converted to a DataFrame. 
    :param start_residue_num: The residue to use to compute the coordinates of the hole entrance; this id done by taking the average of the 
        coordinates of residue residue_num in each chain specified by chain_ids.
    :param end_residue_num: The residue to use to compute the coordinates of the hole exit; this id done by taking the average of the 
        coordinates of residue residue_num in each chain specified by chain_ids.
    :param chain_ids:
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


def hole_main(paths, output_dir:str='../data/genes/hole/output', overwrite:bool=False, chain_ids:list=None, start_residue_num:int=1, end_residue_num:int=220,):
    '''
    
    :param paths: A list of PDB files to run HOLE on. 
    :param chain_ids: The IDs for the chains to use to compute the HOLE starting point and direction. 
    '''
    output_paths = list()
    for path in paths: 
        name = os.path.basename(path).replace('.pdb', '')
        output_path = os.path.join(output_dir, f'{name}.sph')
        input_path = os.path.join(os.path.dirname(path), f'{name}.inp')

        # assert os.path.exists(input_path), f'hole_main: File {input_path} does not exist.'
        assert os.path.exists(path), f'hole_main: File {path} does not exist.'
        
        if os.path.exists(output_path) and (not overwrite):
            continue 
        else:
            cpoint, cvect = hole_get_axis(StructureFile.from_file(path).to_df(), chain_ids=chain_ids, start_residue_num=start_residue_num, end_residue_num=end_residue_num)
            cpoint, cvect = ' '.join([str(n) for n in cpoint]), ' '.join([str(n) for n in cvect])
            hole_write_input(input_path, pdb_path=path, cvect=cvect, cpoint=cpoint, output_path=output_path)

            cmd = f'hole < {input_path} > /dev/null'
            print('hole_main:', cmd)
            subprocess.run(cmd, shell=True, check=True)
            
        output_paths.append(sph_to_pdb(output_path, overwrite=overwrite))
    return output_paths # Return the paths to the PDB outputs. 



def hole_load_sph_center_coords(name:str, output_dir:str='../data/genes/hole/output'):

    
    path = os.path.join(output_dir, f'{name}.pdb')

    pdb_df = StructureFile.from_file(path, format='pdb').to_df(atoms=['QSS'], residues=['SPH'])
    pdb_df = pdb_df.sort_values('residue_num') # Residue number is the position along the pore. 
    # This should be in order of position through the channel (which is the pseudo-"residue number" in the PDB file). 
    # Should be a numpy array with shape (n_atoms, 3). 
    sph_center_coords = np.array(pdb_df['QSS'].tolist()) 
    return sph_center_coords


def hole_load_coords(name:str, input_dir:str='../data/genes/hole/output'):    
    path = os.path.join(input_dir, f'{name}.pdb')
    pdb_df = StructureFile.from_file(path, format='pdb').to_df().drop(columns='b_factor')
    pdb_df = pdb_df.melt(value_name='coord', var_name='atom', id_vars=['residue_num', 'residue', 'chain_id'])
    pdb_df = pdb_df.dropna()

    coords = np.array(pdb_df['coord'].tolist())
    return pdb_df, coords



def hole_process_output(name:str, input_dir:str='../data/genes/hole/input', output_dir:str='../data/genes/hole/output') -> pd.DataFrame:
    '''Locates the atom in each chain closest to the pore center for each sphere in the HOLE output. This is accomplished by computing the pairwise distance between
    each sphere and all other atoms in the PDB file, and then returning the atom from each chain which is closest to the sphere center. 
    For each sphere, there will be {num_chains} nearest atoms.
    
    :param name:
    :returns: 
    '''
    hole_df = list()

    sph_center_coords = hole_load_sph_center_coords(name, output_dir=output_dir)
    pdb_df, coords = hole_load_coords(name, input_dir=input_dir)

    for i, coord in enumerate(sph_center_coords):
        coord = np.expand_dims(coord, axis=0)
        dists = cdist(coord, coords).ravel()

        df = pdb_df.copy().assign(distance=dists, sph_idx=i)
        df = df.sort_values('distance', ascending=True).drop_duplicates('chain_id')
        hole_df.append(df)

    hole_df = pd.concat(hole_df)
    hole_df['name'] = name 
    hole_df['hole_length'] =  max(pdist(sph_center_coords))
    return hole_df


    