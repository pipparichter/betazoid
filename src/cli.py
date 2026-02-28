import argparse
# import src.graph as graph 
import src.align 
import src.recruit 
import os 
import pandas as pd 

NAME = 'test'
OUTPUT_DIR = f'/home/philippar/out/{NAME}'
if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)


READS_PATH_1 = '/groups/banfield/sequences/2014/16ft_4/raw.d/16ft_4_CZBZ.6237.3.40316_trim_clean.PE.1.fastq.gz'
READS_PATH_2 = '/groups/banfield/sequences/2014/16ft_4/raw.d/16ft_4_CZBZ.6237.3.40316_trim_clean.PE.2.fastq.gz'


def recruit():
    # recruit(ref_path:str, n_iters:int=5, reads_path_1:str=None, reads_path_2:str=None)
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref-path', type=str)
    parser.add_argument('--output-dir', type=str, default=OUTPUT_DIR)
    parser.add_argument('--n-iters', type=int, default=2)
    parser.add_argument('--reads-path-1', type=str, default=READS_PATH_1)
    parser.add_argument('--reads-path-2', type=str, default=READS_PATH_2)
    args = parser.parse_args()

    output_path = os.path.join(OUTPUT_DIR, 'reads.csv')
    df = src.recruit(args.ref_path, n_iters=args.n_iters, reads_path_1=args.reads_path_1, reads_path_2=args.reads_path_2, output_dir=args.output_dir)
    df.to_csv(output_path)
    print('recruit: Recruited read information written to', output_path)


def align():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reads-path', type=str, default=os.path.join(OUTPUT_DIR, 'reads.fasta')) # This is the output of recruit reads. 
    args = parser.parse_args()

    src.align.align(os.path.join(OUTPUT_DIR, 'reads.fasta'), output_dir=OUTPUT_DIR)