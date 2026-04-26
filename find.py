import os 
import glob 
import pandas as pd 
import numpy as np 
import argparse 
import re 



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('pattern', dtype=str)
    parser.add_argument('--root', default='/groups/banfield/sequences/')
    # parser.add_argument('--levels', default=2, dtype=int)
    args = parser.parse_args()

    dirs = list()
    for path in glob.glob(os.path.join(args.root, '*'), recursive=True):
        name = os.path.basename(path)
        if os.path.isdir(path) and (re.search(args.pattern, name) is not None):
            dirs.append(path)

    reads_paths = dict()
    for dir_path in dirs:
        reads_path_pattern = os.path.join(dir_path, '**', '*.gz')
        reads_paths[dir_path] = sorted(list(glob.glob(reads_path_pattern, recursive=True)))

    print(reads_paths)