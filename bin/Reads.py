#!/usr/bin/env python3

from Bio import SeqIO

class Read():
    """
    Object for a fastq read
    """
    def __init__(self, read):
        self.read = read

    def get_anchor_target_list(self, anchor_mode, window_slide, lookahead, kmer_size):
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
        last_base = len(self.read) - (lookahead + 2 * kmer_size)

        for i in range(0, last_base, step_size):
            # get anchor
            anchor = self.read[0+i : kmer_size+i]

            # get target start and end positions, as a function of anchor end
            target_start = (kmer_size+i) + lookahead
            target_end = target_start + kmer_size

            # get target
            target = self.read[target_start:target_end]

            if "N" not in anchor and "N" not in target:
                anchor_target_list.append((anchor, target))
