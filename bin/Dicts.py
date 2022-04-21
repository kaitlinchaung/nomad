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
import nltk
import math
import logging
import utils
import Anchors


class AnchorAbundance():
    def __init__(self, anchor_freeze_threshold, num_samples):
        self.num_anchors = 0
        self.dict = {}
        self.anchor_freeze_threshold = anchor_freeze_threshold
        self.num_samples = num_samples

    def initialise(self, anchor, iteration):
        # increment total count
        self.num_anchors += 1

        # if this is a new anchor, increment relative to sample number and iteration
        L = math.floor(iteration / 100000)
        increment = self.num_samples * L
        self.dict[anchor] = 1 + increment

    def update(self, anchor, iteration):
        self.dict[anchor] += self.num_samples

    def is_full(self):
        if self.num_anchors >= self.anchor_freeze_threshold:
            return True
        else:
            return False

    def count(self, anchor):
        return self.dict[anchor]

    def decrement_num_anchors(self):
        self.num_anchors -= 1

    def get_low_abundance_anchors(self, anchor_min_count):
        """
        Return anchors that do not have at least anchor_min_count number of counts
        """
        return [anchor for anchor, count in self.dict.items() if count < anchor_min_count]

    def contains(self, anchor):
        if anchor in self.dict:
            return True
        else:
            return False


class AnchorTargetsSamples():
    def __init__(self, sample_index_dict):
        self.dict = {}
        self.sample_index_dict = sample_index_dict
        self.name_dict = {v:k for k,v in self.sample_index_dict.items()}

    def initialise(self, anchor, target, sample):
        """
        Initialise counts for (anchor, target) for a sample
        """
        # get sample index
        sample_index = self.sample_index_dict[sample]

        # add anchor to dict
        self.dict[anchor] = {target : [0] * len(self.sample_index_dict)}

        # update
        self.dict[anchor][target][sample_index] = 1

    def update(self, anchor, target, sample):
        """
        Update counts for (anchor, target) for a sample
        """
        # get sample index
        sample_index = self.sample_index_dict[sample]

        # add target to dict
        if target not in self.dict[anchor]:
            self.dict[anchor].update({target : [0] * len(self.sample_index_dict)})

        # update
        self.dict[anchor][target][sample_index] += 1

    def num_targets(self, anchor):
        """
        Returns the number of targets an anchor currently has
        """
        return len(self.dict[anchor])

    def get_anchor_counts_df(self, anchor):
        """
        Return a table of anchors x samples counts for a specific anchor, sorted by abundance
        """
        df = pd.DataFrame(self.dict[anchor]).T
        df.columns = [self.name_dict[c] for c in df.columns]
        df['count'] = df.sum(axis=1)
        df = (
            df
            .sort_values('count', ascending=False)
            .drop('count', axis=1)
            .reset_index()
            .rename(columns={'index':'target'})
        )
        return df


class Ignorelist():

    def __init__(self, anchor_abundance, anchor_targets_samples, anchor_targets_scores, anchor_targets_distances, max_ignorelist):
        self.dict = {}
        self.anchor_abundance = anchor_abundance
        self.anchor_targets_samples = anchor_targets_samples
        self.anchor_targets_distances = anchor_targets_distances
        self.anchor_targets_scores = anchor_targets_scores
        self.max_ignorelist = max_ignorelist

    def contains(self, anchor):
        if anchor in self.dict:
            return True
        else:
            return False

    def update(self, anchor):
        """
        Updates ignorelist and removes anchor from all of our dictionaries
        """
        # only add anchor to ignorlist if not in read_counter_freeze or if ignorelist is not larger than 4M
        if len(self.dict) < self.max_ignorelist:
            self.dict[anchor] = True

        # decrement total count
        self.anchor_abundance.decrement_num_anchors()

        # pop anchor from all dicts
        self.anchor_abundance.dict.pop(anchor, None)
        self.anchor_targets_samples.dict.pop(anchor, None)
        self.anchor_targets_samples.dict.pop(anchor, None)
        self.anchor_targets_scores.dict.pop(anchor, None)
        self.anchor_targets_distances.mu.pop(anchor, None)


class AnchorTargetsScores():
    def __init__(self):
        self.dict = {}

    def initialise(self, anchor, scores):
        # increment total count
        self.dict[anchor] = scores

    def contains(self, anchor):
        if anchor in self.dict:
            return True
        else:
            return False

    def get_final_anchors(self, num_keep_anchors, use_std):
        """
        Get final anchors for parse_anchors
        """
        if len(self.dict) < num_keep_anchors:
            anchors = list(self.dict.keys())
        else:
            df = pd.DataFrame(self.dict).T
            if use_std:
                df['z'] = (
                    df
                    .std(axis=1)
                    .abs()
                )
            else:
                df['z'] = (
                    df
                    .sum(axis=1)
                    .abs()
                )

            anchors = (
                df['z']
                .sort_values(ascending=False)
                .head(num_keep_anchors)
                .index
                .to_list()
            )
        return anchors


class AnchorTargetsDistances(dict):
    """
    Dictionary object to keep track of anchors, targets, and their distances
    """
    def __init__(self):
        self.dict = {}
        self.mu = {}
