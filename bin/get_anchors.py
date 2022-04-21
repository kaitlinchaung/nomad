#!/usr/bin/env python3

import pandas as pd
import gzip
import os
import re
import sys
import argparse
import numpy as np
from Bio import SeqIO
import math
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from datetime import datetime
import logging
import time
import utils
import Dicts


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n_iterations",
        type=int
    )
    parser.add_argument(
        "--max_reads",
        type=int
    )
    parser.add_argument(
        "--kmer_size",
        type=int
    )
    parser.add_argument(
        "--lookahead",
        type=int
    )
    parser.add_argument(
        "--samplesheet",
        type=str
    )
    parser.add_argument(
        "--target_counts_threshold",
        type=int
    )
    parser.add_argument(
        "--anchor_counts_threshold",
        type=int
    )
    parser.add_argument(
        "--anchor_freeze_threshold",
        type=int
    )
    parser.add_argument(
        "--read_freeze_threshold",
        type=int
    )
    parser.add_argument(
        "--anchor_mode",
        type=str
    )
    parser.add_argument(
        "--window_slide",
        type=int
    )
    parser.add_argument(
        "--num_keep_anchors",
        help='up or down',
        type=int
    )
    parser.add_argument(
        "--anchor_score_threshold",
        help='up or down',
        type=int
    )
    parser.add_argument(
        "--c_type",
        type=str
    )
    parser.add_argument(
        "--outfile",
        help='up or down',
        type=str
    )
    parser.add_argument(
        "--use_std",
        type=str
    )
    parser.add_argument(
        "--compute_target_distance",
        type=str
    )
    parser.add_argument(
        "--max_ignorelist",
        type=int
    )
    parser.add_argument(
        "--distance_type",
        type=str
    )
    parser.add_argument(
        "--score_type",
        type=str
    )

    args = parser.parse_args()
    return args


def main():
    args = get_args()

    # read in samples from samplesheet, formatted as [fastq_file, optional group_id]
    sample_list = pd.read_csv(
        args.samplesheet,
        header=None
    )

    # redefine bool
    if args.use_std == "true":
        use_std = True
    elif args.use_std == "false":
        use_std = False
    # if samplesheet only has 1 column, force use_std
    if sample_list.shape[1] == 1:
        use_std = True

    # set up logging
    logging.basicConfig(
        filename = f'get_anchors.log',
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        filemode='w',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # # logging: print input parametrs
    # utils.log_params(args, use_std)

    # get list of samples from fastq_files (if fastq_file = "file1.fastq.gz", sample = "file1")
    samples = (
        sample_list
        .iloc[:,0]
        .apply(
            lambda x:
            os.path.basename(x).split('.')[0]
        )
        .tolist()
    )

    # get sample-related dicts
    group_ids_dict, sample_index_dict = utils.get_sample_dict(
        sample_list,
        samples,
        use_std
    )

    # initialise objects
    anchor_abundance = Dicts.AnchorAbundance(args.anchor_freeze_threshold, len(samples))
    anchor_targets_samples = Dicts.AnchorTargetsSamples(sample_index_dict)
    anchor_targets_scores = Dicts.AnchorTargetsScores()
    anchor_targets_distances = Dicts.AnchorTargetsDistances()
    ignorelist = Dicts.Ignorelist(
        anchor_abundance,
        anchor_targets_samples,
        anchor_targets_scores,
        anchor_targets_distances,
        args.max_ignorelist
        )

    # begin iterations
    for iteration in range(1, args.n_iterations+1):

        num_reads = iteration * args.max_reads

        # get reads from this iteration's fastq files
        read_chunk = utils.get_read_chunk(
            iteration,
            samples,
            args.n_iterations
        )

        # accumulate anchors for this iteration
        utils.get_iteration_summary_scores(
            read_chunk,
            iteration,
            anchor_abundance,
            anchor_targets_samples,
            anchor_targets_scores,
            anchor_targets_distances,
            ignorelist,
            group_ids_dict,
            args.kmer_size,
            args.lookahead,
            args.max_reads,
            args.target_counts_threshold,
            args.anchor_counts_threshold,
            args.anchor_freeze_threshold,
            args.anchor_mode,
            args.window_slide,
            args.distance_type,
            args.score_type
        )

        # abundance threshold: ignorelist anchors that don't meet an abundance requirement
        utils.ignorelist_abundance(
            anchor_abundance,
            ignorelist,
            num_reads,
            args.c_type,
            samples
        )

        # get scores for each anchor
        for anchor in anchor_targets_scores.dict:
            scores, _, _ = utils.compute_phase_1_scores(
                anchor_targets_samples.get_anchor_counts_df(anchor),
                group_ids_dict,
                args.kmer_size,
                args.distance_type,
                "force"
            )

            # updates scores for each anchor
            anchor_targets_scores.initialise(anchor, scores)

    # after all iterations are done accumulating, calculate the summary score
    keep_anchors = anchor_targets_scores.get_final_anchors(args.num_keep_anchors, use_std)

    # filter for these anchors and drop any duplicates
    final_anchors = (
        pd.DataFrame(keep_anchors, columns = ['anchor'])
        .drop_duplicates()
    )

    ## output this final anchors list
    final_anchors.to_csv(args.outfile, index=False, sep='\t')



main()
