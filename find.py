import os 
import glob 
import pandas as pd 
import numpy as np 
import argparse 
import re 



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('pattern')
    parser.add_argument('--root', default='/groups/banfield/sequences/')
    parser.add_argument('--root', default='/groups/banfield/scratch/projects/environmental/RES/int/LDS/sandpiper_not_in_Logan/')
    parser.add_argument('--levels', default=2)
    parser.add_argument('--match-filename', action='store_true')
    args = parser.parse_args()

    if not args.match_filename:
        dirs = list()
        dir_pattern = [args.root] + ['*'] * args.levels
        dir_pattern = os.path.join(*dir_pattern)
        # print(list(glob.glob(dir_pattern, recursive=True)))
        for path in glob.glob(dir_pattern, recursive=True):
            dir_name = os.path.basename(os.path.normpath(path)) 
            if os.path.isdir(path) and (re.search(args.pattern, dir_name) is not None):
                dirs.append(path)

        paths = dict()
        for dir_path in dirs:
            reads_path_pattern = os.path.join(dir_path, '**', '*.gz')
            paths[dir_path] = sorted(list(glob.glob(reads_path_pattern, recursive=True)))

        for dir_path, reads_paths in paths.items():
            dir_name = os.path.basename(os.path.normpath(dir_path)) 
            print(f'{dir_name}\t{reads_paths[0]}\t{reads_paths[-1]}')
    
    else:
        reads_path_pattern = [args.root] + ['*'] * args.levels + ['*.gz']
        reads_path_pattern = os.path.join(*reads_path_pattern)
        reads_paths = sorted(list(glob.glob(reads_path_pattern, recursive=True)))
        assert len(reads_path_pattern) == 2, f'Expected two matching reads paths, got: {reads_paths}'
        print(f'{reads_paths[0]}\t{reads_paths[-1]}')