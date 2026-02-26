import argparse
import src.graph as graph 



def recruit_reads():
    # recruit_reads(ref_path:str, n_iters:int=5, reads_path_1:str=None, reads_path_2:str=None)
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref-path', type=str)
    parser.add_argument('--n-iters', type=int, default=5)
    parser.add_argument('--reads-path-1', type='str')
    parser.add_argument('--reads-path-2', type='str')

    bam_df
