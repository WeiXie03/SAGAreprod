from _utils import *
from _pipeline import *
from _reproducibility import *
from _interpret import *
from _visualize import *
from _concat import *
from _chromhmm import *
from get_data import *
from matchlabels import *
from run import *
import pandas as pd
import glob, os 
import mp as mp
from functools import partial
import matplotlib
import argparse
from pathlib import Path

'''
This script is to be run after steps 1 - 2 in `run.py`,
especially creating genomedatas for Segway.

Segway:
    3 - run segway (separate and concatenated)

ChromHMM:
    2 - Binarize data
    3 - run chromhmm
    
4 - parse and bin results
4.5  - biological label interpretation
5 - match labels (for separate runs)
6 - cluster labels
7 - reproducibility analysis (at different levels of abstraction)

'''
def RunParse_segway_replicates_single_ct_paramd(celltype_dir: str, stw: float, res: int, output_dir: str, random_seed: int):
    '''
    Same as `RunParse_segway_replicates` but with user-specified
    segment transition weight, [0, 1], and resolution, an int > 0

    `celltype_dir` must containing the genomedatas for the cell type
    This should have no trailing '/'.

    `output_dir` is the path for the directory for the Segway output
    with the given hyperparameters pair, with a trailing '/'.
    '''
    # path "leaf" of `celltype_dir`, which is the cell type's name
    name_sig = celltype_dir.split('/')[-1]
    num_tracks = len(
        [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)])

    params_dict_1 = {
        "random_seed":random_seed, "track_weight":0.01,
        "stws":stw, "ruler_scale":100, "prior_strength":1, "resolution":res, 
        "mini_batch_fraction":0.03, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
        "name_sig":output_dir+name_sig+'_rep1', 
        "genomedata_file":celltype_dir+'/rep1.genomedata', 
        "traindir":output_dir+name_sig+'rep1'+'_train', 
        "posteriordir":output_dir+name_sig+'rep1'+'_posterior'
    }
    
    if os.path.exists(params_dict_1['name_sig']) == False:
        print('Running Segway celltype {} Rep1'.format(celltype_dir))
        run_segway_and_post_process(params_dict_1)

    else:
        print(params_dict_1['name_sig'], "already exists")

    params_dict_2 = {
        "random_seed":random_seed, "track_weight":0.01,
        "stws":stw, "ruler_scale":100, "prior_strength":1, "resolution":res, 
        "mini_batch_fraction":0.03, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
        "name_sig":output_dir+name_sig+'_rep2', 
        "genomedata_file":celltype_dir+'/rep2.genomedata', 
        "traindir":output_dir+name_sig+'rep2'+'_train', 
        "posteriordir":output_dir+name_sig+'rep2'+'_posterior'
    }
    
    if os.path.exists(params_dict_2['name_sig']) == False:
        print('Running Segway celltype {} Rep2'.format(celltype_dir))
        run_segway_and_post_process(params_dict_2)

    else:
        print(params_dict_2['name_sig'], "already exists")

def RunParse_segway_replicates_paramd(stw: float, res: int, master_out_dir: str, celltype_dirs: list[str], rand_seed: int, mp_pool=None):
    '''
    Runs and parses replicates for all specified cell types with
    given segment transition weight and resolution parameters.

    Cell types are specified by `celltype_dirs`, a list of
    paths of the directories containing their genomedata's.
    These should have no trailing '/'.

    `master_out_dir` should include a trailing '/'.
    '''
    out_dir_path = f"{master_out_dir}stw{str(stw)}_res{str(res)}"
    if not os.path.isdir(out_dir_path):
        os.mkdir(out_dir_path)
    else:
        print(out_dir_path, "already exists, adding to it")

    if type(mp_pool) == mp.pool.Pool:
        # Run segway replicates     MP
        partial_runs_i = partial(RunParse_segway_replicates_single_ct_paramd,
                                 stw=stw, res=res, output_dir=out_dir_path, random_seed=rand_seed)
        mp_pool.map(partial_runs_i, celltype_dirs)

def hyperparams2Dexpmt_RunParse_segway_replicates(master_out_dir: str, celltype_dirs: str, stw_range_params: tuple[float], resolution_range_params: tuple[int], n_procs: int, resolution_range_base=4, random_seed=73):
    '''
    Runs and parses replicates with Segway with a 2D grid of settings for
    the segment transition weight and resolution hyperparameters.

    Output
    ===
    Stores the output directories of Segway runs in the following hierarchy:
    - `master_out_dir/`
        - `stw0.1_res100/`
            - `hepatocyte/`
            - `HeLA-S3/`
            .
            .
            .
        .
        .
        .

    Arguments
    ===
        stw_range_params: a tuple of the arguments for range, (start, stop, range)
        resolution_range_params: a tuple of the arguments for the resolution exponent range, (start, stop, range)
    '''
    stws_array = list(range(*stw_range_params))
    res_array = [resolution_range_base**res for res in range(*resolution_range_params)]
    hyperparam_pairs = list(itertools.product(stws_array, res_array))
    # calculate number of hyperparameters configurations/settings/pairs
    n_configs = len(hyperparam_pairs)

    with mp.Manager() as manager:
        # create a process pool, local to the main thread/process
        with mp.Pool(n_procs) as pool:
            # issue all top-level tasks
            _partial_f = functools.partial(RunParse_segway_replicates_paramd,
                                           rand_seed=random_seed,
                                           mp_pool=pool,
                                           master_out_dir=master_out_dir,
                                           celltype_dirs=celltype_dirs)
            pool.starmap(_partial_f, hyperparam_pairs)

# def RunParse_segway_psdreps(celltype_dir, output_dir, random_seed=73):
#     name_sig = celltype_dir.split('/')[-1]
#     num_tracks = len(
#         [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)])

#     concat_param_dict = {
#         "random_seed":random_seed, "track_weight":0.01,"stws":1, "ruler_scale":100, 
#         "prior_strength":1, "resolution":100, "mini_batch_fraction":0.03,
#         "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
#         "name_sig_concat":output_dir+name_sig+'_rep1psd_concat', 
#         "name_sig_rep1":output_dir+name_sig+'_concat_rep1psd1', 
#         "name_sig_rep2":output_dir+name_sig+'_concat_rep1psd2',

#         "genomedata_file_concat":celltype_dir+'/rep1psd_concatenated.genomedata',  
#         "genomedata_file_rep1":celltype_dir+'/rep1psd1_concat.genomedata', 
#         "genomedata_file_rep2":celltype_dir+'/rep1psd2_concat.genomedata',

#         "traindir":output_dir+name_sig+'_rep1psd_concat_train', 
#         "posteriordir_rep1":output_dir+name_sig+'_rep1psd_concat_posterior_1',
#         "posteriordir_rep2":output_dir+name_sig+'_rep1psd_concat_posterior_2'}

#     if os.path.exists(concat_param_dict['name_sig_rep1']) == False or \
#          os.path.exists(concat_param_dict['name_sig_rep2']) == False:
#         print('Running Segway concatenated celltype {}...'.format(celltype_dir))
#         concat_segwayrun_and_postprocess(concat_param_dict)

# def RunParse_segway_param_init(celltype_dir, replicate_number, random_seeds, output_dir):
#     # replicate number should be in format "repN" -> i.e. rep1, rep2
#     name_sig = celltype_dir.split('/')[-1]
#     num_tracks = len(
#         [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)])

#     params_dict_1 = {
#         "random_seed":random_seeds[0], "track_weight":0.01,
#         "stws":1, "ruler_scale":100, "prior_strength":1, "resolution":100, 
#         "mini_batch_fraction":0.03, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
#         "name_sig":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[0]), 
#         "genomedata_file":celltype_dir+'/{}.genomedata'.format(replicate_number), 
#         "traindir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[0])+'_train', 
#         "posteriordir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[0])+'_posterior'
#     }

#     params_dict_2 = {
#         "random_seed":random_seeds[1], "track_weight":0.01,
#         "stws":1, "ruler_scale":100, "prior_strength":1, "resolution":100, 
#         "mini_batch_fraction":0.03, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
#         "name_sig":output_dir+name_sig+'_{}_rs{}'.format(replicate_number,random_seeds[1]),
#         "genomedata_file":celltype_dir+'/{}.genomedata'.format(replicate_number), 
#         "traindir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[1])+'_train', 
#         "posteriordir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[1])+'_posterior'
#     }

#     if os.path.exists(params_dict_1['name_sig']) == False:
#         print('Running Segway parameter initialization test on celltype {}, {}, with random seed {}'.format(
#             celltype_dir, replicate_number, random_seeds[0]))
#         run_segway_and_post_process(params_dict_1)

#     else:
#         print(params_dict_1['name_sig'], "already exists")

#     if os.path.exists(params_dict_2['name_sig']) == False:
#         print('Running Segway parameter initialization test on celltype {}, {}, with random seed {}'.format(
#             celltype_dir, replicate_number, random_seeds[1]))
#         run_segway_and_post_process(params_dict_2)

#     else:
#         print(params_dict_2['name_sig'], "already exists")

# def intersect_parsed_posteriors(parsed_df_dir_1, parsed_df_dir_2):
#     is_chmm_concat = bool(
#         ("chromhmm_runs" in parsed_df_dir_1) and ("chromhmm_runs" in parsed_df_dir_2) and 
#         ("concat" in parsed_df_dir_1) and ("concat" in parsed_df_dir_2))

#     if is_chmm_concat:
#         parsed_df_dir_1 = parsed_df_dir_1.replace("parsed_posterior.csv", "parsed_posterior_rep1.csv")
#         parsed_df_dir_2 = parsed_df_dir_2.replace("parsed_posterior.csv", "parsed_posterior_rep2.csv")
    
#     df1 = pd.read_csv(parsed_df_dir_1, on_bad_lines="skip", encoding_errors="ignore").drop("Unnamed: 0", axis=1)
#     df1.iloc[:, 3:] = df1.iloc[:, 3:].astype("float16")
#     df2 = pd.read_csv(parsed_df_dir_2, on_bad_lines="skip", encoding_errors="ignore").drop("Unnamed: 0", axis=1)
#     df2.iloc[:, 3:] = df2.iloc[:, 3:].astype("float16")

#     if df1.columns[3] == "posterior0" and df2.columns[3] == "posterior1":
#         df2.columns =  list(df2.columns[:3]) + ["posterior"+str(i)for i in range(len(df2.columns)-3)]
#     elif df1.columns[3] == "posterior1" and df2.columns[3] == "posterior0":
#         df1.columns =  list(df1.columns[:3]) + ["posterior"+str(i)for i in range(len(df1.columns)-3)]

        
#     # to handle concat indexing
#     if "_1" in df1.iloc[0, 0] or "_2" in df1.iloc[0, 0]:
#         chrdf1 = list(df1.chr)
#         for i in range(len(chrdf1)):
#             if "_1" in chrdf1[i]:
#                 chrdf1[i] = chrdf1[i].replace("_1", "")
#             elif "_2" in chrdf1[i]:
#                 chrdf1[i] = chrdf1[i].replace("_2", "")
#         df1.chr = np.array(chrdf1)

#     if "_1" in df2.iloc[0, 0] or "_2" in df2.iloc[0, 0]:
#         chrdf2 = list(df2.chr)
#         for i in range(len(chrdf2)):
#             if "_2" in chrdf2[i]:
#                 chrdf2[i] = chrdf2[i].replace("_2", "")
#             elif "_1" in chrdf2[i]:
#                 chrdf2[i] = chrdf2[i].replace("_1", "")
#         df2.chr = np.array(chrdf2)

#     intersect = pd.merge(
#         df1, 
#         df2, 
#         how='inner', on=['chr', 'start', 'end'])

#     df1 = [intersect.chr, intersect.start, intersect.end]
#     df2 = [intersect.chr, intersect.start, intersect.end]

#     for c in intersect.columns:
#         if c[-1] == 'x':
#             df1.append(intersect[c])
#         elif c[-1] == 'y':
#             df2.append(intersect[c])

#     df1 = pd.concat(df1, axis=1)
#     df1.columns = [c.replace("_x", "") for c in df1.columns]
    
#     df2 = pd.concat(df2, axis=1)
#     df2.columns = [c.replace("_y", "") for c in df2.columns]

#     return df1, df2

def init_argparser(parser: argparse.ArgumentParser) -> argparse.Namespace:
    """
    - downloads dir
    - dir for segway runs
    - dir for final reports
    """
    parser.add_argument("download_dir", type=Path, help="Path to directory of data (contains subdirectories for cell types)")
    parser.add_argument("segway_runs_dir", type=Path, help="Path to directory in which to put output of Segway, including posteriors")
    parser.add_argument("results_dir", type=Path, help="Path to directory in which to put results of reproducibility analysis")
    return parser.parse_args()

"""when running the whole script from start to end to generate (and reproduce) results
remember to put label interpretation in try blocks (skippable) to prevent any kind of
dependency issue of segtools to cause issues in reproducibility of results"""

if __name__=="__main__":
    CellType_list = np.array(
        ['K562', 'MCF-7', 'GM12878', 'HeLa-S3', 'CD14-positive monocyte'])

    # download_dir = 'files/'
    # segway_dir = 'segway_runs/'
    # res_dir = 'reprod_results/'
    parser = argparse.ArgumentParser()
    args = init_argparser(parser)
    download_dir = str(args.download_dir)+'/'
    segway_dir = str(args.segway_runs_dir)+'/'
    res_dir = str(args.results_dir)+'/'
    
    if os.path.exists(res_dir) == False:
        os.mkdir(res_dir)

    print('list of target celltypes', CellType_list)
    CellType_list = [ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+ct)]
    # potential space characters in directory names should
    # already have been cleaned up to prevent later issues

    if os.path.exists(segway_dir) == False:
        os.mkdir(segway_dir)

    print("Starting Segway runs")
    
    # ======= Hyperparameter Ranges =======
    STW_RANGE_PARAMS = (0, 1, 0.1)
    RES_BASE = 4
    RES_EXP_RANGE_PARAMS = (0, 12)

    hyperparams2Dexpmt_RunParse_segway_replicates(
        segway_dir, CellType_list, STW_RANGE_PARAMS, RES_EXP_RANGE_PARAMS,
        mp.cpu_count(), RES_BASE)

    exit()

    # parse_posteriors 
    print('Checking for unparsed posteriors...')
    list_of_seg_runs = [
        d for d in os.listdir(segway_dir) if os.path.isdir(segway_dir+'/'+d)]
    print(list_of_seg_runs)
    for d in list_of_seg_runs:
        print('-Checking for {}  ...'.format(segway_dir+'/'+d+'/parsed_posterior.csv'))

        if os.path.exists(segway_dir+'/'+d+'/parsed_posterior.csv') == False:
            parse_posterior_results(segway_dir+'/'+d, 100, mp=False)

        else:
            print('-Exists!')

    print('All parsed!')


def genomedata_from_dict(input_dict):
    tracklist = ''