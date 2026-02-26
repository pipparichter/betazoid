import argparse
import src.graph as graph 


OUTPUT_DIR = '/home/philippar/out/'
JOB_NAME = 'test'

READS_PATH_1 = '/groups/banfield/sequences/2014/16ft_4/raw.d/16ft_4_CZBZ.6237.3.40316_trim_clean.PE.1.fastq.gz'
READS_PATH_2 = '/groups/banfield/sequences/2014/16ft_4/raw.d/16ft_4_CZBZ.6237.3.40316_trim_clean.PE.2.fastq.gz'


def recruit_reads():
    # recruit_reads(ref_path:str, n_iters:int=5, reads_path_1:str=None, reads_path_2:str=None)
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref-path', type=str)
    parser.add_argument('--job-name', type=str, default=JOB_NAME)
    parser.add_argument('--n-iters', type=int, default=2)
    parser.add_argument('--reads-path-1', type='str', default=READS_PATH_1)
    parser.add_argument('--reads-path-2', type='str', default=READS_PATH_2)
    args = parser.parse_args()

    df = graph.recruit_reads(args.job_name, args.ref_path, n_iters=args.n_iters, reads_path_1=args.reads_path_1, reads_path_2=args.reads_path_2, output_dir=OUTPUT_DIR)
    df.to_csv(f'{args.job_name}.csv')