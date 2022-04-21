#!/usr/bin/env python3

import pandas as pd
import gzip
import os
import re
import sys
import argparse
import logging
import numpy as np
from Bio import SeqIO
import math
import nltk
import time
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from datetime import datetime


def get_sample_dict(sample_list, samples, use_std):
    """
    Return group_ids_dict and sample_index_dict
    """
    # if not using standard deviation score, make dict of {sample : group_id}
    group_ids_dict = {}
    if not use_std or sample_list.shape[1] != 1:
        group_ids = sample_list.iloc[:,1].tolist()
        for i in range(0, len(samples)):
            group_ids_dict[samples[i]] = group_ids[i]
    else:
        for i in range(0, len(samples)):
            group_ids_dict[samples[i]] = 1

    # create sample index dict, of {sample : sample_id}
    sample_index_dict = {}
    for i in range(len(samples)):
        sample_index_dict[samples[i]] = i

    return group_ids_dict, sample_index_dict


def get_read_chunk(iteration, samples, n_iterations):
    """
    Returns reads from all files for this iteration
    """
    read_chunk = []
    for sample in samples:

        try:
            n_digits = len(str(n_iterations))
            fastq_file = f'{sample}_{str(iteration).zfill(n_digits)}.fastq.gz'

            with gzip.open(fastq_file, 'rt') as fastq:
                for record in SeqIO.parse(fastq, 'fastq'):
                    read = str(record.seq)
                    read_chunk.append((read, sample))
        except:
            pass

    return read_chunk


def get_anchor_target_list(read, anchor_mode, window_slide, lookahead, kmer_size):
    """
    Returns a list of valid (anchor, target) from a read
    """
    # initialise
    anchor_target_list = []

    # define, chunk will be disjoint while tile will be overlapping
    if anchor_mode == 'chunk':
        step_size = kmer_size
    elif anchor_mode == 'tile':
        step_size = window_slide

    # only parse read up to where it is possible to have a valid target
    last_base = len(read) - (lookahead + 2 * kmer_size)

    for i in range(0, last_base, step_size):
        # get anchor
        anchor = read[0+i : kmer_size+i]

        # get target start and end positions, as a function of anchor end
        target_start = (kmer_size+i) + lookahead
        target_end = target_start + kmer_size

        # get target
        target = read[target_start:target_end]

        if "N" not in anchor and "N" not in target:
            anchor_target_list.append((anchor, target))

    return anchor_target_list


def is_valid_anchor(anchor, anchor_abundance, anchor_freeze_threshold, ignorelist):
    """
    Return True if anchor and targets are both valid for further computation
    """

    if anchor_abundance.is_full():
        if anchor_abundance.contains(anchor):
            return True
        else:
            False
    else:
        if anchor_abundance.contains(anchor):
            return True
        else:
            if not ignorelist.contains(anchor):
                return True
            else:
                return False


    ## if we are no longer accumulating new anchors, only continue with already active anchors
    if anchor_abundance.is_full():
        if anchor.active:
            is_valid = True
        else:
            is_valid = False
    # if we are accumulating new anchors, continue with either:
    #   1. already active anchors
    #   2. new anchors that are not ignorelisted
    else:
        if anchor.active or (not anchor.active and not anchor.ignorelisted):
            is_valid = True
        else:
            is_valid = False

    return is_valid


def get_distance(seq1, seq2, distance_type):
    """
    Returns the distance between 2 strings
    """
    if distance_type.lower() == 'hamming':
        distance = 0
        for x,y in zip(seq1, seq2):
            distance = distance + (x != y)

    elif distance_type.lower() == 'jaccard':

        # hardcode kmer size for now
        k = 7
        set1 = set([seq1[i:i+k] for i in range(len(seq1)-k+1)])
        set2 = set([seq2[i:i+k] for i in range(len(seq2)-k+1)])

        distance = len(set1.intersection(set2)) / len(set1.union(set2))

    return distance


def compute_phase_1_scores(anchor_counts, group_ids_dict, kmer_size, distance_type, score_type):
    # intialise a df for targets x distances
    distance_df = pd.DataFrame(index=anchor_counts['target'], columns=['distance'])

    # iterate over targets sorted by decreasing abundance
    # get min hamming distance of each target to its preceeding targets
    for index, row in anchor_counts.iterrows():
        target = row['target']
        candidates = (
            anchor_counts['target']
            .loc[0:index-1]
            .to_numpy()
        )

        # if it's the first aka most abundant target, initialise the distance to be 0
        if index == 0:
            min_dist = 0

        # else, get the min distance to all preceeding targets
        else:
            min_dist = min([get_distance(target, x, distance_type) for x in candidates])

        # update target distances
        distance_df.loc[distance_df.index==target, 'distance'] = min_dist
        anchor_counts.loc[anchor_counts['target']==target, 'distance'] = min_dist

    # intialise
    p = pd.DataFrame()
    p['i'] = range(0, kmer_size+1)
    p['counts'] = 0

    # for each distance, get the the total number of occurrences of that distance
    for i, df in anchor_counts.groupby('distance'):
        count = (
            df
            .drop(['target', 'distance'], axis=1)
            .values
            .sum()
        )
        p.loc[p['i']==i, 'counts'] = count

    # get proportion of counts with that distance
    p['p_hat'] = p['counts'] / p['counts'].sum()

    # get u
    mu = (p['i'] * p['p_hat']).sum()

    if score_type == 'slow' or score_type == "force":
        # center the distances with u
        anchor_counts['distance_centered'] = anchor_counts['distance'] - mu

        # prepare df to get terms
        anchor_counts = (
            anchor_counts
            .drop(['distance'], axis=1)
            .set_index('target')
            )

        m = anchor_counts.drop('distance_centered', axis=1)
        distance_centered = anchor_counts['distance_centered']

        numerator = (
            m
            .mul(pd.Series(group_ids_dict))
            .mul(distance_centered, axis=0)
            .sum(axis=0)
        )

        denominator = m.sum(axis=0)

        phase_1_scores = (
            numerator
            .divide(denominator, fill_value=0, axis=0)
        )

        distances = distance_df['distance']

    elif score_type == "fast":
        phase_1_scores = pd.DataFrame()
        distances = pd.DataFrame()

    return phase_1_scores, distances, mu


def ignorelist_abundance(anchor_abundance, ignorelist, num_reads, c_type, samples):
    """
    Update ignorelist with an anchor if its coutns do not meet an abundance requirement
    """

    # define
    k = math.floor(num_reads / 100000)
    if c_type == "num_samples":
        c = len(samples)
    elif c_type == "constant":
        c = 1

    # define min number of counts required per anchor to prevent ignorelisting
    anchor_min_count = k * c

    # get the anchors that do not meet the abundance requirement of anchor_min_count
    ignorelist_anchors_abundance = anchor_abundance.get_low_abundance_anchors(anchor_min_count)

    # add those anchors to the ignorelist
    for anchor in ignorelist_anchors_abundance:
        ignorelist.update(anchor)

    return

def get_iteration_summary_scores(
        read_chunk,
        iteration,
        anchor_abundance,
        anchor_targets_samples,
        anchor_targets_scores,
        anchor_targets_distances,
        ignorelist,
        group_ids_dict,
        kmer_size,
        lookahead,
        max_reads,
        target_counts_threshold,
        anchor_counts_threshold,
        anchor_freeze_threshold,
        anchor_mode,
        window_slide,
        distance_type,
        score_type
    ):
    """
    Return the summary scores for the reads of one iteration
    """
    # logging: intialise
    start_time = time.time()
    invalid_anchor = 0
    valid_anchor = 0
    phase_0 = 0
    phase_1 = 0
    phase_1_compute_score = 0
    phase_1_ignore_score = 0

    # set mu_threshold
    if distance_type == 'hamming':
        mu_threshold = 2
    elif distance_type == 'jaccard':
        mu_threshold = 0.2

    # get reads for this iteration
    for read_tuple in read_chunk:

        # define read_tuple
        read, sample = read_tuple

        # get list of all anchors from each
        anchor_target_list = get_anchor_target_list(read, anchor_mode, window_slide, lookahead, kmer_size)

        # loop over each anchor in the list of all anchors from each read
        for anchor, target in anchor_target_list:

            # if we are still accumulating new anchors and this is not an ignorelisted anchor
            if not anchor_abundance.is_full() and not ignorelist.contains(anchor):
                valid_anchor += 1

                # intialise anchor
                anchor_abundance.initialise(anchor, iteration)
                anchor_targets_samples.initialise(anchor, target, sample)

            # if we are done accumulating anchors and this is an anchor we have seen before
            elif anchor_abundance.is_full() and anchor_abundance.contains(anchor):
                valid_anchor += 1

                # update anchor counts
                anchor_abundance.update(anchor, iteration)

                # check if we are doing phase 2 computations
                if not anchor_targets_scores.contains(anchor):

                    # update targets since we are still accumulating targets at this point
                    anchor_targets_samples.update(anchor, target, sample)

                    # check if we have reached phase 1
                    phase_1_conditions = [
                        anchor_targets_samples.num_targets(anchor) == target_counts_threshold,
                        anchor_abundance.count(anchor) >= anchor_counts_threshold
                    ]

                    if any(phase_1_conditions):
                        phase_1 += 1
                        # compute phase_1 score
                        scores, topTargets_distances, mu = compute_phase_1_scores(
                            anchor_targets_samples.get_anchor_counts_df(anchor),
                            group_ids_dict,
                            kmer_size,
                            distance_type,
                            score_type
                        )

                        if mu < mu_threshold:
                            ignorelist.update(anchor)
                            phase_1_ignore_score += 1

                        else:
                            anchor_targets_scores.initialise(anchor, scores)
                            phase_1_compute_score += 1
                    else:
                        phase_0 += 1

                # do phase_2 stuff
                else:
                    pass

            else:
                invalid_anchor += 1

    # logging: get total run time
    run_time = time.time() - start_time

    # logging: get anchor percentages
    try:
        valid_anchor_percent = round((valid_anchor / (valid_anchor + invalid_anchor)) * 100, 2)
        invalid_anchor_percent = round((invalid_anchor / (valid_anchor + invalid_anchor)) * 100, 2)
        phase_0_percent = round((phase_0 / (phase_0+phase_1)) * 100, 2)
        phase_1_percent = round((phase_1 / (phase_0+phase_1)) * 100, 2)

    except:
        valid_anchor_percent = 0
        invalid_anchor_percent = 0
        phase_0_perecent = 0
        phase_1_percent = 0

    # logging: output for this iteration
    logging.info("")
    logging.info("-----------------------------------------------------------------------------------------------------------------")
    logging.info("")
    logging.info(f'i = {iteration}')
    logging.info("")
    logging.info(f'\tIteration run time = {round(run_time, 2)} seconds')
    logging.info("")
    logging.info(f'\tinvalid anchors = {invalid_anchor} ({invalid_anchor_percent}%)')
    logging.info("")
    logging.info(f'\tvalid anchors = {valid_anchor} ({valid_anchor_percent}%)')
    logging.info(f'\t\tphase_0 anchors = {phase_0} ({phase_0_percent}%)')
    logging.info(f'\t\tphase_1 anchors = {phase_1} ({phase_1_percent}%)')
    logging.info(f'\t\t\tscore passed diversity condition = {phase_1_compute_score}')
    logging.info(f'\t\t\tscore failed diversity condition = {phase_1_ignore_score}')
    logging.info("")
    logging.info(f'\tnumber of anchors with candidate scores = {len(anchor_targets_scores.dict)}')
    logging.info(f'\t\tsize of anchor_counts dict = {len(anchor_abundance.dict)}')
    logging.info(f'\t\tsize of anchor_targets_samples dict = {len(anchor_targets_samples.dict)}')
    logging.info(f'\t\tsize of anchor_topTargets_scores dict = {len(anchor_targets_scores.dict)}')
    logging.info(f'\t\tsize of anchor_topTargets_distances dict = {len(anchor_targets_distances.dict)}')
    logging.info("")
    logging.info(f'\tignorelist size = {len(ignorelist.dict)}')
    logging.info(f'\t\tanchors that have failed the diversity requirement = {phase_1_ignore_score}')




